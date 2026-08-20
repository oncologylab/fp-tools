"""Review multiple diff-footprints HTML reports in one interactive page."""

from __future__ import annotations

import argparse
import base64
import copy
import gzip
import json
import re
from pathlib import Path
from fp_tools.tools.static_comparison_browser import (
    build_static_browser,
    write_embedded_static_browser,
)
from fp_tools.utils.project_layout import comparisons_dir, is_project_layout, project_root, review_output_path


DEFAULT_INPUT_GLOB = "diff_footprints_*.html"


def _compressed_json_b64(payload: dict) -> str:
    text = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    return base64.b64encode(gzip.compress(text.encode("utf-8"), compresslevel=9)).decode("ascii")


def _decode_payload_b64(payload_b64: str) -> dict:
    return json.loads(gzip.decompress(base64.b64decode(payload_b64)).decode("utf-8"))


def read_diff_html_payload(path: str | Path) -> dict:
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    match = re.search(r'const\s+reportPayloadB64\s*=\s*"([^"]+)"', text)
    if not match:
        match = re.search(r'reportPayloadB64\s*=\s*"([^"]+)"', text)
    if not match:
        raise ValueError(f"Could not find reportPayloadB64 in {path}")
    payload = _decode_payload_b64(match.group(1))
    if "points" not in payload:
        raise ValueError(f"{path} does not look like a diff-footprints HTML payload")
    return payload


def discover_input_htmls(inputs: list[str | Path]) -> list[Path]:
    paths: list[Path] = []
    for value in inputs:
        path = Path(value)
        if path.is_dir():
            paths.extend(sorted(path.rglob(DEFAULT_INPUT_GLOB)))
        else:
            paths.append(path)
    seen = set()
    unique = []
    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(path)
    if not unique:
        raise ValueError("No diff-footprints HTML inputs were found")
    return unique


def _comparison_label(path: Path, payload: dict, override: str | None = None) -> str:
    if override:
        return override
    report_label = str(payload.get("report_label") or "").strip()
    if report_label:
        method = report_label.split(";")[0].replace("Method:", "").replace("Normalization:", "").strip()
        if method:
            return method[:80]
    title = str(payload.get("title") or "").strip()
    title_match = re.match(r"^Differential footprint report\s*\((.+)\)$", title)
    if title_match:
        return title_match.group(1).strip()[:80]
    if title and title != "Differential footprint report":
        return title[:80]
    return path.parent.name or path.stem


def build_review_payload(paths: list[str | Path], labels: list[str] | None = None, title: str = "Review multiple differential footprint comparisons") -> dict:
    html_paths = discover_input_htmls(paths)
    if labels and len(labels) != len(html_paths):
        raise ValueError("--labels must have the same length as resolved comparison HTML inputs")
    comparisons = []
    for idx, path in enumerate(html_paths):
        payload = read_diff_html_payload(path)
        comparisons.append({
            "label": _comparison_label(Path(path), payload, labels[idx] if labels else None),
            "path": str(path),
            "payload": payload,
        })
    return {"schema": "fp-tools.review-multi-comparisons.v1", "title": title, "comparisons": comparisons}


def _aggregate_prefixes(payload: dict) -> set[str]:
    return {str(motif.get("prefix")) for motif in (payload.get("aggregate") or {}).get("motifs") or [] if motif.get("prefix")}


def count_missing_aggregate_profiles(review_payload: dict) -> tuple[int, int]:
    missing = 0
    total = 0
    for item in review_payload.get("comparisons") or []:
        payload = item.get("payload") or {}
        aggregate_prefixes = _aggregate_prefixes(payload)
        for point in payload.get("points") or []:
            prefix = point.get("prefix")
            if not prefix:
                continue
            total += 1
            if str(prefix) not in aggregate_prefixes:
                missing += 1
    return missing, total


def _infer_aggregate_flank(review_payload: dict, requested: str | int | None = "auto") -> int:
    if requested not in (None, "auto"):
        flank = int(requested)
        if flank < 1:
            raise ValueError("--aggregate-flank must be at least 1")
        return flank
    for item in review_payload.get("comparisons") or []:
        x_values = ((item.get("payload") or {}).get("aggregate") or {}).get("x") or []
        numeric = []
        for value in x_values:
            try:
                numeric.append(float(value))
            except (TypeError, ValueError):
                pass
        if numeric:
            return max(1, int(max(abs(min(numeric)), abs(max(numeric)) + 1)))
    return 100


def _coerce_profile_to_axis(profile: list, source_x: list, target_x: list) -> list[float]:
    target_len = len(target_x)
    values = []
    for value in profile or []:
        try:
            values.append(float(value))
        except (TypeError, ValueError):
            values.append(float("nan"))
    if not target_x:
        return values
    if len(values) == target_len:
        return values
    if source_x and len(source_x) == len(values):
        try:
            import numpy as np

            src = np.asarray(source_x, dtype=float)
            val = np.asarray(values, dtype=float)
            tgt = np.asarray(target_x, dtype=float)
            finite = np.isfinite(src) & np.isfinite(val)
            if int(finite.sum()) >= 2:
                return [round(float(x), 6) for x in np.interp(tgt, src[finite], val[finite])]
        except Exception:
            pass
    if not values:
        return [0.0] * target_len
    if len(values) > target_len:
        return values[:target_len]
    return values + [values[-1]] * (target_len - len(values))


def _coerce_aggregate_to_axis(aggregate: dict, target_x: list) -> dict:
    source_x = aggregate.get("x") or target_x
    for condition in aggregate.get("conditions") or []:
        if "profile" in condition:
            condition["profile"] = _coerce_profile_to_axis(condition.get("profile") or [], source_x, target_x)
        for sample in condition.get("samples") or []:
            sample["profile"] = _coerce_profile_to_axis(sample.get("profile") or [], source_x, target_x)
    aggregate["x"] = list(target_x)
    return aggregate


def fill_missing_aggregate_profiles(
    review_payload: dict,
    project: str | Path | None = None,
    fill_missing: bool = True,
    recompute_missing: bool = False,
    aggregate_flank: str | int | None = "auto",
    cores: int | None = None,
) -> dict:
    before_missing, total = count_missing_aggregate_profiles(review_payload)
    if before_missing == 0:
        return {"before_missing": 0, "after_missing": 0, "filled": 0, "total": total}
    if recompute_missing and project is None:
        raise ValueError("--recompute-missing-aggregate-profiles requires --outdir in project layout")

    from fp_tools.tools.motif_aggregate_grid import ordered_comparisons, prepare_aggregate_maps

    flank = _infer_aggregate_flank(review_payload, aggregate_flank)
    aggregate_maps = prepare_aggregate_maps(
        review_payload,
        project=project,
        fill_missing=fill_missing or recompute_missing,
        recompute_missing=recompute_missing,
        flank=flank,
        cores=cores,
    )
    filled = 0
    for comparison in ordered_comparisons(review_payload):
        payload = comparison.payload
        aggregate_block = payload.setdefault("aggregate", {})
        target_x = aggregate_block.get("x") or list(range(-flank, flank))
        aggregate_block["x"] = list(target_x)
        motifs = aggregate_block.setdefault("motifs", [])
        existing = _aggregate_prefixes(payload)
        for point in payload.get("points") or []:
            prefix = str(point.get("prefix") or "")
            if not prefix or prefix in existing:
                continue
            aggregate, source = aggregate_maps.get((comparison.index, prefix), (None, "missing"))
            if not aggregate:
                continue
            aggregate = _coerce_aggregate_to_axis(copy.deepcopy(aggregate), target_x)
            aggregate["profile_source"] = source
            motifs.append(aggregate)
            existing.add(prefix)
            filled += 1
    after_missing, _ = count_missing_aggregate_profiles(review_payload)
    return {"before_missing": before_missing, "after_missing": after_missing, "filled": filled, "total": total}


def write_review_html(
    review_payload: dict,
    output: str | Path,
    display_panels: int = 4,
    aggregate_legends: str = "show",
) -> None:
    if display_panels < 4 or display_panels > 8:
        raise ValueError("display_panels must be between 4 and 8")
    if aggregate_legends not in {"show", "hide"}:
        raise ValueError("aggregate_legends must be 'show' or 'hide'")
    write_embedded_static_browser(review_payload, output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="review-multi-comparisons", description="Combine diff-footprints reports into a static browser bundle or one self-contained HTML file.")
    parser.add_argument("--inputs", nargs="+", help="diff-footprints HTML files or directories containing diff_footprints_*.html files; directories are searched recursively.")
    parser.add_argument("--labels", nargs="*", help="Optional labels, one per resolved input HTML.")
    outputs = parser.add_mutually_exclusive_group()
    outputs.add_argument("--output-dir", help="Output directory for index.html, JavaScript, CSS, and compact static data files.")
    outputs.add_argument("--output-html", help="One self-contained HTML report; aggregate profiles are optional.")
    outputs.add_argument("--output", dest="output_dir", help=argparse.SUPPRESS)
    parser.add_argument("--outdir", help="Project directory used with --layout project.")
    parser.add_argument("--layout", choices=["custom", "project"], default="project", help="Use fp-tools standard project output layout under --outdir (default: project when only --outdir is provided).")
    parser.add_argument("--display-panels", type=int, default=4, help=argparse.SUPPRESS)
    parser.add_argument("--default-comparison", nargs=2, metavar=("<group1>", "<group2>"), help="Region or condition pair initially shown in the static browser")
    parser.add_argument("--default-aggregate-motifs", nargs="+", metavar="<motif>", help="Ordered motif IDs, names, or output prefixes initially shown")
    parser.add_argument("--default-aggregate-plots", type=int, metavar="<int>", help="Number of aggregate profiles initially shown (default: 4; maximum: 12)")
    parser.add_argument("--documentation-url", metavar="<url>", help="Optional link back to the documentation site")
    parser.add_argument("--aggregate-legends", choices=["show", "hide"], default="show", help=argparse.SUPPRESS)
    parser.add_argument("--fill-missing-aggregate-profiles", action="store_true", help="Fill missing motif aggregate panels from profiles embedded elsewhere in the combined review payload.")
    parser.add_argument("--recompute-missing-aggregate-profiles", action="store_true", help="Recompute still-missing motif aggregate panels from project sample bigWigs and match-motifs BEDs.")
    parser.add_argument("--aggregate-flank", default="auto", help="Flank used when recomputing missing aggregate profiles, or 'auto' to match the existing report axis (default: auto).")
    parser.add_argument("--cores", type=int, default=None, help="Worker processes for --recompute-missing-aggregate-profiles (default: all available cores).")
    parser.add_argument("--title", default="Review multiple differential footprint comparisons")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project = None
    if is_project_layout(args.layout) and args.outdir:
        project = project_root(args.outdir)
        if not args.inputs:
            args.inputs = [str(comparisons_dir(project))]
            if not args.output_dir and not args.output_html:
                legacy_path = review_output_path(project)
                args.output_dir = str(legacy_path.with_suffix(""))
    if not args.inputs:
        parser.error("provide --inputs or use --layout project with --outdir")
    if not args.output_dir and not args.output_html:
        parser.error("provide --output-dir, --output-html, or use --layout project with --outdir")
    if args.recompute_missing_aggregate_profiles and project is None:
        parser.error("--recompute-missing-aggregate-profiles requires --outdir in project layout")
    try:
        payload = build_review_payload(args.inputs, labels=args.labels, title=args.title)
        if args.fill_missing_aggregate_profiles or args.recompute_missing_aggregate_profiles:
            fill_stats = fill_missing_aggregate_profiles(
                payload,
                project=project,
                fill_missing=args.fill_missing_aggregate_profiles or args.recompute_missing_aggregate_profiles,
                recompute_missing=args.recompute_missing_aggregate_profiles,
                aggregate_flank=args.aggregate_flank,
                cores=args.cores,
            )
        else:
            if args.output_html:
                fill_stats = None
            else:
                before_missing, total = count_missing_aggregate_profiles(payload)
                fill_stats = {
                    "before_missing": before_missing,
                    "after_missing": before_missing,
                    "filled": 0,
                    "total": total,
                }
        if args.output_html:
            write_review_html(
                payload,
                args.output_html,
                display_panels=args.display_panels,
                aggregate_legends=args.aggregate_legends,
            )
            index_path = Path(args.output_html)
        else:
            index_path = build_static_browser(
                [item["payload"] for item in payload["comparisons"]],
                args.output_dir,
                title=args.title,
                default_comparison=args.default_comparison,
                default_motifs=args.default_aggregate_motifs,
                default_aggregate_plots=(
                    args.default_aggregate_plots
                    if args.default_aggregate_plots is not None
                    else args.display_panels
                ),
                documentation_url=args.documentation_url,
            )
    except ValueError as exc:
        parser.error(str(exc))
    print(f"Wrote {index_path}")
    if fill_stats and fill_stats["before_missing"]:
        print(f"Aggregate profiles: filled {fill_stats['filled']} missing panels; {fill_stats['after_missing']} remain missing of {fill_stats['total']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
