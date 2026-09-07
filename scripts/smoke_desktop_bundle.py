#!/usr/bin/env python3
"""Smoke-test a frozen fp-tools GUI executable."""

from __future__ import annotations

import argparse
import os
import socket
import signal
import subprocess
import tempfile
import time
import urllib.request
from pathlib import Path


def write_signature_fixture(root: Path) -> dict[str, Path]:
    """Create a tiny three-cell-type input set for frozen plotting smoke tests."""

    import anndata as ad
    import numpy as np
    import pandas as pd

    root.mkdir(parents=True, exist_ok=True)
    barcodes = ["B1-1", "B2-1", "M1-1", "M2-1", "T1-1", "T2-1"]
    cell_types = ["B_cell", "B_cell", "Monocyte", "Monocyte", "T_NK_cell", "T_NK_cell"]
    annotations = pd.DataFrame(
        {
            "barcode": barcodes,
            "cell_type": cell_types,
            "snap_cell_type": cell_types,
            "umap_1": [-3.0, -2.8, 0.0, 0.2, 3.0, 3.2],
            "umap_2": [1.0, 1.2, -1.0, -0.8, 1.0, 1.2],
        }
    )
    annotations_path = root / "annotations.tsv"
    annotations.to_csv(annotations_path, sep="\t", index=False)

    h5ad_path = root / "cells.h5ad"
    adata = ad.AnnData(
        X=np.asarray(
            [
                [8, 2],
                [7, 2],
                [4, 5],
                [4, 6],
                [2, 8],
                [2, 7],
            ],
            dtype=np.float32,
        ),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(
            {"selected": [True, True]},
            index=["chr1:0-500", "chr1:500-1000"],
        ),
    )
    adata.write_h5ad(h5ad_path)

    marker_centers = {"STAT6": 100, "CEBPA": 300, "ZNF683": 700}
    marker_groups = {
        "STAT6": "B_cell",
        "CEBPA": "Monocyte",
        "ZNF683": "T_NK_cell",
    }
    site_dir = root / "motif_sites"
    site_dir.mkdir()
    fragment_lines = []
    for marker, center in marker_centers.items():
        (site_dir / f"{marker}.motif_hits.bed").write_text(
            f"chr1\t{center - 2}\t{center + 2}\t{marker}\n",
            encoding="utf-8",
        )
        for barcode, cell_type in zip(barcodes, cell_types, strict=True):
            if cell_type == marker_groups[marker]:
                start, end = center - 13, center + 14
            else:
                start, end = center - 2, center + 3
            fragment_lines.append(f"chr1\t{start}\t{end}\t{barcode}\t1")
    fragments_path = root / "fragments.tsv"
    fragments_path.write_text("\n".join(fragment_lines) + "\n", encoding="utf-8")

    score_rows = []
    score_patterns = {
        "STAT6_M1": ("STAT6", "B_cell", [2.0, 1.8, -0.8, -0.7, -0.6, -0.5]),
        "CEBPA_M2": ("CEBPA", "Monocyte", [-0.7, -0.6, 2.0, 1.8, -0.5, -0.4]),
        "ZNF683_M3": ("ZNF683", "T_NK_cell", [-0.6, -0.5, -0.7, -0.6, 2.0, 1.8]),
    }
    for motif_id, (tf_name, dominant, values) in score_patterns.items():
        means = {
            group: float(np.mean([value for value, cell in zip(values, cell_types, strict=True) if cell == group]))
            for group in ("B_cell", "Monocyte", "T_NK_cell")
        }
        score_rows.append(
            {
                "motif_id": motif_id,
                "tf_name": tf_name,
                "dominant_cell_type": dominant,
                "dynamic_range": float(max(values) - min(values)),
                "B_cell_mean_z": means["B_cell"],
                "Monocyte_mean_z": means["Monocyte"],
                "T_NK_cell_mean_z": means["T_NK_cell"],
                **dict(zip(barcodes, values, strict=True)),
            }
        )
    all_scores_path = root / "all_motif_scores.tsv"
    pd.DataFrame(score_rows).to_csv(all_scores_path, sep="\t", index=False)
    return {
        "annotations": annotations_path,
        "fragments": fragments_path,
        "h5ad": h5ad_path,
        "site_dir": site_dir,
        "all_scores": all_scores_path,
    }


def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as handle:
        handle.bind(("127.0.0.1", 0))
        return int(handle.getsockname()[1])


def stop_process_tree(process: subprocess.Popen[str]) -> str:
    """Stop the GUI and any server child while preserving captured output."""

    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True,
            text=True,
            check=False,
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return ""
    try:
        return process.communicate(timeout=20)[0] or ""
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        return process.communicate(timeout=20)[0] or ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("--timeout", type=float, default=90.0)
    args = parser.parse_args()

    executable = Path(args.executable).resolve()
    help_run = subprocess.run(
        [str(executable), "--fp-tools-internal-smoke-command-helps"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if help_run.returncode != 0 or "Validated " not in help_run.stdout:
        raise SystemExit(
            f"Internal command dispatch audit failed:\n{help_run.stdout}\n{help_run.stderr}"
        )

    raw_read_run = subprocess.run(
        [str(executable), "--fp-tools-internal-command", "prepare-atac", "--help"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if raw_read_run.returncode == 0:
        raise SystemExit("Desktop bundle unexpectedly exposes prepare-atac")

    examples_run = subprocess.run(
        [str(executable), "--fp-tools-internal-list-gui-examples"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if examples_run.returncode != 0 or "bulk_footprinting_bam.yml" not in examples_run.stdout:
        raise SystemExit(
            f"Desktop bundle is missing GUI examples:\n{examples_run.stdout}\n{examples_run.stderr}"
        )

    with tempfile.TemporaryDirectory(prefix="fp-tools-desktop-smoke-") as run_dir:
        run_path = Path(run_dir)
        fixture_root = Path(__file__).resolve().parents[1] / "test_data"
        call_output = run_path / "call_footprints"
        call_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-command",
                "call-footprints",
                "--signals",
                str(fixture_root / "Bcell_corrected.bw"),
                "--regions",
                str(fixture_root / "plot_regions.bed"),
                "--outdir",
                str(call_output),
                "--score",
                "footprint",
                "--cores",
                "1",
            ],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            cwd=run_dir,
        )
        if call_run.returncode != 0 or not any(call_output.glob("*.bw")):
            raise SystemExit(
                f"Frozen call-footprints smoke failed:\n{call_run.stdout}\n{call_run.stderr}"
            )

        normalize_output = run_path / "normalize_bigwig"
        normalize_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-command",
                "normalize-bigwig",
                "--bigwigs",
                str(fixture_root / "Bcell_corrected.bw"),
                str(fixture_root / "Tcell_corrected.bw"),
                "--background",
                str(fixture_root / "plot_regions.bed"),
                "--outdir",
                str(normalize_output),
                "--method",
                "background-scale",
                "--stat",
                "q95",
                "--target",
                "median",
                "--workers",
                "2",
            ],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            cwd=run_dir,
        )
        if normalize_run.returncode != 0 or len(list(normalize_output.glob("*.bw"))) != 2:
            raise SystemExit(
                "Frozen multiprocessing normalize-bigwig smoke failed:\n"
                f"{normalize_run.stdout}\n{normalize_run.stderr}"
            )

        signature_fixture = write_signature_fixture(run_path / "signature_fixture")
        signature_output = run_path / "find_signature_fp"
        signature_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-command",
                "find-signature-fp",
                "--annotations",
                str(signature_fixture["annotations"]),
                "--fragments",
                str(signature_fixture["fragments"]),
                "--h5ad",
                str(signature_fixture["h5ad"]),
                "--tf-site-dir",
                str(signature_fixture["site_dir"]),
                "--all-motif-score-table",
                str(signature_fixture["all_scores"]),
                "--markers",
                "STAT6,CEBPA,ZNF683",
                "--marker-groups",
                "STAT6:B_cell,CEBPA:Monocyte,ZNF683:T_NK_cell",
                "--knn",
                "1",
                "--flank",
                "20",
                "--center-half-width",
                "2",
                "--flank-inner",
                "5",
                "--flank-outer",
                "15",
                "--no-create-fragment-index",
                "--top-motif-signatures-per-cell-type",
                "1",
                "--outdir",
                str(signature_output),
            ],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            cwd=run_dir,
        )
        signature_svgs = [
            path for path in signature_output.glob("*.svg") if path.stat().st_size > 0
        ]
        signature_pdfs = [
            path for path in signature_output.glob("*.pdf") if path.stat().st_size > 0
        ]
        if signature_run.returncode != 0 or not signature_svgs or not signature_pdfs:
            raise SystemExit(
                "Frozen find-signature-fp SVG/PDF smoke failed:\n"
                f"{signature_run.stdout}\n{signature_run.stderr}\n"
                f"SVG files: {len(signature_svgs)}; PDF files: {len(signature_pdfs)}"
            )

        port = free_port()
        process = subprocess.Popen(
            [
                str(executable),
                "--fp-tools-internal-gui-server",
                "--port",
                str(port),
                "--run-dir",
                run_dir,
                "--no-browser",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=run_dir,
            env={**os.environ, "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false"},
            start_new_session=os.name != "nt",
        )
        try:
            deadline = time.monotonic() + args.timeout
            last_error: Exception | None = None
            ready = False
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    output = process.stdout.read() if process.stdout else ""
                    raise SystemExit(f"GUI exited before becoming ready ({process.returncode}):\n{output}")
                try:
                    with urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2) as response:
                        if response.status == 200:
                            print(f"Desktop GUI ready on port {port}")
                            ready = True
                            break
                except Exception as exc:  # startup polling
                    last_error = exc
                    time.sleep(0.5)
            if not ready:
                raise SystemExit(f"GUI did not become ready within {args.timeout:g} seconds: {last_error}")
        finally:
            output = stop_process_tree(process)
            if process.returncode not in (0, -15, 1):
                print(output)

        native_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-native-window-smoke",
                "--run-dir",
                str(run_path / "native-window"),
            ],
            capture_output=True,
            text=True,
            timeout=max(args.timeout, 150.0),
            cwd=run_dir,
            env={
                **os.environ,
                "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false",
                "QTWEBENGINE_CHROMIUM_FLAGS": "--disable-gpu",
            },
        )
        if native_run.returncode != 0:
            raise SystemExit(
                "Native desktop-window smoke failed:\n"
                f"{native_run.stdout}\n{native_run.stderr}"
            )
        print("Native fp-tools window rendered successfully")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
