"""Build the static multi-comparison browser used by fp-tools reports."""

from __future__ import annotations

import base64
import csv
from datetime import date
import gzip
import hashlib
import html
from importlib.resources import files
import json
from pathlib import Path
import re
import shutil


PROFILE_SHARDS = 16


def _safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.")
    if not name:
        raise ValueError(f"Could not create a safe filename from {value!r}")
    return name


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_gzip_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    path.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))


def _read_gzip_json(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def write_embedded_static_browser(review_payload: dict, output: str | Path) -> Path:
    """Package the shared browser and a review payload into one portable HTML file."""
    comparisons = review_payload.get("comparisons") or []
    if not comparisons:
        raise ValueError("No comparison payloads were supplied")
    for record in comparisons:
        payload = record.get("payload") or {}
        if len(payload.get("conditions") or []) != 2 or not payload.get("points"):
            raise ValueError("Each input must be a two-condition diff-footprints report")

    raw = json.dumps(
        review_payload,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")
    payload_b64 = base64.b64encode(
        gzip.compress(raw, compresslevel=9, mtime=0)
    ).decode("ascii")
    has_aggregates = any(
        ((record.get("payload") or {}).get("aggregate") or {}).get("motifs")
        for record in comparisons
    )
    template_root = files("fp_tools.resources.static_browser")
    document = template_root.joinpath("index.html").read_text(encoding="utf-8")
    stylesheet = template_root.joinpath("styles.css").read_text(encoding="utf-8")
    application = template_root.joinpath("app.js").read_text(encoding="utf-8")
    plot_controls = template_root.joinpath("plot_controls.js").read_text(
        encoding="utf-8"
    )
    title = html.escape(
        str(
            review_payload.get("title")
            or "Review multiple differential footprint comparisons"
        )
    )
    bootstrap = (
        f'<script>const reportPayloadB64="{payload_b64}",'
        f'hasAggregateProfiles={str(has_aggregates).lower()};'
        "window.fpToolsBrowserBootstrap={mode:\"embedded\","
        "payloadB64:reportPayloadB64};</script>"
    )
    document = document.replace(
        '<link rel="stylesheet" href="styles.css" />',
        f"<style>\n{stylesheet}\n</style>",
    )
    document = document.replace(
        '<script src="plot_controls.js" defer></script>',
        f"<script>\n{plot_controls}\n</script>",
    )
    document = document.replace(
        '<script src="app.js" defer></script>',
        f"{bootstrap}\n<script>\n{application}\n</script>",
    )
    document = document.replace(
        "<title>Differential footprint report</title>",
        f"<title>{title}</title>",
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")
    return output


def _profile_shard(prefix: str) -> int:
    return int(hashlib.sha256(prefix.encode("utf-8")).hexdigest()[:2], 16) % PROFILE_SHARDS


def split_browser_payload(payload: dict) -> tuple[dict, list[list[dict]]]:
    """Separate fast motif summaries from lazily loaded profile arrays."""
    aggregate = payload.get("aggregate") or {}
    motifs = aggregate.get("motifs") or []
    if not motifs:
        raise ValueError("A comparison payload has no aggregate motif profiles")
    core = {key: value for key, value in payload.items() if key not in {"aggregate", "motif_matrices"}}
    core["logos"] = {}
    summaries = []
    shards: list[list[dict]] = [[] for _ in range(PROFILE_SHARDS)]
    for motif in motifs:
        prefix = str(motif.get("prefix") or "")
        if not prefix:
            continue
        shard = _profile_shard(prefix)
        conditions = []
        for condition in motif.get("conditions") or []:
            conditions.append(
                {
                    **{key: value for key, value in condition.items() if key != "profile"},
                    "samples": [
                        {key: value for key, value in sample.items() if key != "profile"}
                        for sample in condition.get("samples") or []
                    ],
                }
            )
        summaries.append(
            {
                **{key: value for key, value in motif.items() if key != "conditions"},
                "conditions": conditions,
                "profile_shard": shard,
            }
        )
        shards[shard].append(motif)
    motif_matrices = payload.get("motif_matrices") or {}
    core["motif_matrices"] = {
        prefix: motif_matrices[prefix]
        for prefix in (motif["prefix"] for motif in summaries)
        if prefix in motif_matrices
    }
    core["aggregate"] = {**{key: value for key, value in aggregate.items() if key != "motifs"}, "motifs": summaries}
    return core, shards


def _logo_pngs(payloads: list[dict]) -> dict[str, bytes]:
    logos: dict[str, bytes] = {}
    for payload in payloads:
        for prefix, record in (payload.get("logos") or {}).items():
            safe_prefix = _safe_name(prefix)
            if safe_prefix != prefix or prefix in logos:
                continue
            uri = str((record or {}).get("png") or "")
            if not uri.startswith("data:image/png;base64,"):
                continue
            image = base64.b64decode(uri.split(",", 1)[1])
            if image.startswith(b"\x89PNG\r\n\x1a\n"):
                logos[prefix] = image
    return logos


def _condition_records(payloads: list[dict]) -> list[dict]:
    names: list[str] = []
    samples: dict[str, list[str]] = {}
    for payload in payloads:
        for condition in payload.get("conditions") or []:
            condition = str(condition)
            if condition not in names:
                names.append(condition)
            samples.setdefault(condition, [])
        motifs = (payload.get("aggregate") or {}).get("motifs") or []
        if not motifs:
            continue
        for condition in motifs[0].get("conditions") or []:
            name = str(condition.get("name") or "")
            if not name:
                continue
            for sample in condition.get("samples") or []:
                sample_name = str(sample.get("name") or "")
                if sample_name and sample_name not in samples.setdefault(name, []):
                    samples[name].append(sample_name)
    return [{"name": name, "samples": samples.get(name, [])} for name in names]


def _write_all_results(payloads: list[dict], comparisons: list[str], output: Path) -> None:
    columns = [
        "comparison", "condition1", "condition2", "prefix", "name", "motif_id",
        "group", "n_profile_sites", "n_motif_regions_condition1",
        "n_motif_regions_condition2", "effect",
        "ci_lower", "ci_upper", "pvalue", "qvalue", "statistical_method",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for comparison, payload in zip(comparisons, payloads):
            conditions = payload.get("conditions") or ["condition1", "condition2"]
            aggregate = {
                item.get("prefix"): item
                for item in (payload.get("aggregate") or {}).get("motifs") or []
            }
            for point in payload.get("points") or []:
                writer.writerow(
                    {
                        "comparison": comparison,
                        "condition1": conditions[0],
                        "condition2": conditions[1],
                        "prefix": point.get("prefix", ""),
                        "name": point.get("name", ""),
                        "motif_id": point.get("motif_id", ""),
                        "group": point.get("group", ""),
                        "n_profile_sites": aggregate.get(point.get("prefix"), {}).get(
                            "n_sites", ""
                        ),
                        "n_motif_regions_condition1": point.get("n_motif_regions_set_1", ""),
                        "n_motif_regions_condition2": point.get("n_motif_regions_set_2", ""),
                        "effect": point.get("change", ""),
                        "ci_lower": point.get("ci_lower", ""),
                        "ci_upper": point.get("ci_upper", ""),
                        "pvalue": point.get("pvalue", ""),
                        "qvalue": point.get("fdr", ""),
                        "statistical_method": point.get("statistical_method", ""),
                    }
                )


def _resolve_default_motifs(payload: dict, selectors: list[str] | None) -> list[str]:
    motifs = (payload.get("aggregate") or {}).get("motifs") or []
    if not selectors:
        configured = (payload.get("aggregate") or {}).get("default_motifs") or []
        return [str(prefix) for prefix in configured if prefix]
    selected = []
    for selector in selectors:
        token = str(selector).strip().casefold()
        matches = {
            str(motif.get("prefix"))
            for motif in motifs
            if token in {
                str(motif.get("prefix") or "").strip().casefold(),
                str(motif.get("motif_id") or "").strip().casefold(),
                str(motif.get("name") or "").strip().casefold(),
            }
        }
        matches.discard("")
        if not matches:
            raise ValueError(f"Unknown default aggregate motif: {selector}")
        if len(matches) > 1:
            raise ValueError(
                f"Ambiguous default aggregate motif {selector!r}; use an exact motif ID or prefix"
            )
        prefix = next(iter(matches))
        if prefix not in selected:
            selected.append(prefix)
    return selected


def build_static_browser(
    payloads: list[dict],
    output_dir: str | Path,
    title: str,
    default_comparison: tuple[str, str] | list[str] | None = None,
    default_motifs: list[str] | None = None,
    default_aggregate_plots: int | None = None,
    documentation_url: str | None = None,
) -> Path:
    """Write an ENCODE-demo-compatible static browser bundle."""
    if not payloads:
        raise ValueError("No comparison payloads were supplied")
    for payload in payloads:
        if len(payload.get("conditions") or []) != 2 or not payload.get("points"):
            raise ValueError("Each input must be a two-condition diff-footprints report")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    if data_dir.exists():
        shutil.rmtree(data_dir)
    reports_dir = data_dir / "reports"
    profiles_dir = data_dir / "profiles"
    logos_dir = data_dir / "logos"
    reports_dir.mkdir(parents=True)
    profiles_dir.mkdir(parents=True)
    logos_dir.mkdir(parents=True)

    comparison_ids: list[str] = []
    metadata_records = []
    seen_pairs: set[frozenset[str]] = set()
    for payload in payloads:
        condition1, condition2 = [str(value) for value in payload["conditions"]]
        pair = frozenset((condition1, condition2))
        if condition1 == condition2 or pair in seen_pairs:
            raise ValueError(f"Duplicate or invalid comparison: {condition1} vs {condition2}")
        seen_pairs.add(pair)
        comparison = f"{_safe_name(condition1)}_vs_{_safe_name(condition2)}"
        comparison_ids.append(comparison)
        compact, shards = split_browser_payload(payload)
        report_path = reports_dir / f"{comparison}.json.gz"
        _write_gzip_json(compact, report_path)
        shard_records = []
        for shard_id, motifs in enumerate(shards):
            shard_path = profiles_dir / comparison / f"{shard_id:02x}.json.gz"
            _write_gzip_json({"motifs": motifs}, shard_path)
            shard_records.append(
                {
                    "id": shard_id,
                    "file": f"data/profiles/{comparison}/{shard_id:02x}.json.gz",
                    "sha256": _sha256(shard_path),
                    "motifs": len(motifs),
                }
            )
        metadata_records.append(
            {
                "comparison": comparison,
                "condition1": condition1,
                "condition2": condition2,
                "file": f"data/reports/{comparison}.json.gz",
                "payload_sha256": _sha256(report_path),
                "profile_shards": shard_records,
                "motifs": len(payload.get("points") or []),
                "aggregate_motifs": len((payload.get("aggregate") or {}).get("motifs") or []),
                "aggregate_site_set": (payload.get("aggregate") or {}).get("site_set", ""),
            }
        )

    _write_all_results(payloads, comparison_ids, data_dir / "all_pairwise_results.tsv.gz")
    logos = _logo_pngs(payloads)
    for prefix, image in logos.items():
        (logos_dir / f"{prefix}.png").write_bytes(image)

    default = metadata_records[0]
    if default_comparison:
        first, second = map(str, default_comparison)
        match = next(
            (
                record
                for record in metadata_records
                if {record["condition1"], record["condition2"]} == {first, second}
            ),
            None,
        )
        if match is None:
            raise ValueError(f"Unknown default comparison: {first} vs {second}")
        default = match
    default_payload = payloads[metadata_records.index(default)]
    resolved_default_motifs = _resolve_default_motifs(default_payload, default_motifs)
    if default_aggregate_plots is None:
        default_aggregate_plots = int(
            (default_payload.get("aggregate") or {}).get("default_plot_count") or 4
        )
    if not 1 <= int(default_aggregate_plots) <= 12:
        raise ValueError("default aggregate plots must be between 1 and 12")
    metadata = {
        "schema": "fp-tools.static-comparison-browser.v1",
        "release_date": date.today().isoformat(),
        "method": "fp-tools differential footprinting",
        "title": title,
        "conditions": _condition_records(payloads),
        "comparisons": metadata_records,
        "default_comparison": {
            "condition1": str(default_comparison[0]) if default_comparison else default["condition1"],
            "condition2": str(default_comparison[1]) if default_comparison else default["condition2"],
        },
        "default_aggregate_motifs": resolved_default_motifs,
        "default_aggregate_plots": min(
            int(default_aggregate_plots), max(1, len(resolved_default_motifs) or 12)
        ),
        "documentation_url": documentation_url or "",
        "downloads": {"all_results": "data/all_pairwise_results.tsv.gz"},
        "logos": {"base": "data/logos", "format": "png", "count": len(logos)},
    }
    (data_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    template_root = files("fp_tools.resources.static_browser")
    for name in ("index.html", "plot_controls.js", "app.js", "styles.css"):
        (output_dir / name).write_bytes(template_root.joinpath(name).read_bytes())
    return output_dir / "index.html"


def read_static_browser_review(bundle_root: str | Path) -> dict:
    """Reconstruct a review payload from a browser bundle without duplicate data."""
    bundle_root = Path(bundle_root)
    metadata = json.loads((bundle_root / "data" / "metadata.json").read_text(encoding="utf-8"))
    comparisons = []
    for record in metadata.get("comparisons") or []:
        payload = _read_gzip_json(bundle_root / record["file"])
        profiles = {}
        for shard in record.get("profile_shards") or []:
            shard_path = bundle_root / shard["file"]
            if shard.get("sha256") and _sha256(shard_path) != shard["sha256"]:
                raise ValueError(f"Profile shard checksum mismatch: {shard_path}")
            profiles.update(
                {motif["prefix"]: motif for motif in _read_gzip_json(shard_path).get("motifs") or []}
            )
        summaries = (payload.get("aggregate") or {}).get("motifs") or []
        payload["aggregate"]["motifs"] = [
            profiles[summary["prefix"]]
            for summary in summaries
            if summary.get("prefix") in profiles
        ]
        comparisons.append(
            {
                "label": f"{record['condition1']} vs {record['condition2']}",
                "payload": payload,
            }
        )
    return {
        "schema": "fp-tools.review-multi-comparisons.v1",
        "title": metadata.get("title") or "Review multiple differential footprint comparisons",
        "comparisons": comparisons,
    }
