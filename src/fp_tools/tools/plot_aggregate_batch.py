"""Interactive batch aggregate HTML reports for match-motifs outputs."""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path

import numpy as np
import pyBigWig


def _read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [dict(row) for row in reader]
    required = {"sample", "signal", "match_dir"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"Manifest is missing required column(s): {', '.join(sorted(missing))}")
    return rows


def _read_bed_centers(path: Path) -> list[tuple[str, int]]:
    centers = []
    if not path.exists():
        return centers
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                centers.append((fields[0], (int(fields[1]) + int(fields[2])) // 2))
            except ValueError:
                continue
    return centers


def _discover_motifs(match_dir: str | Path) -> list[dict[str, str | int | float]]:
    root = Path(match_dir)
    result_files = sorted(root.glob("*_results.txt"))
    motifs = []
    if result_files:
        with result_files[0].open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                prefix = row.get("output_prefix") or row.get("name") or row.get("motif_id")
                if not prefix:
                    continue
                total = float(row.get("total_tfbs") or 0)
                motifs.append({"prefix": prefix, "name": row.get("name") or prefix, "score": total, "sites": int(total)})
    if not motifs:
        for bed in sorted(root.glob("*/beds/*_all.bed")):
            prefix = bed.parent.parent.name
            centers = _read_bed_centers(bed)
            motifs.append({"prefix": prefix, "name": prefix, "score": len(centers), "sites": len(centers)})
    motifs.sort(key=lambda row: (-float(row["score"]), str(row["name"])))
    return motifs


def _mean_profile(signal: str | Path, centers: list[tuple[str, int]], flank: int) -> list[float]:
    profiles = []
    with pyBigWig.open(str(signal)) as bw:
        chroms = bw.chroms()
        for chrom, center in centers:
            if chrom not in chroms:
                continue
            start = center - flank
            end = center + flank
            if start < 0 or end > chroms[chrom] or end <= start:
                continue
            values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            if values.size != flank * 2:
                continue
            profiles.append(np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0))
    if not profiles:
        return [0.0] * (flank * 2)
    return [round(float(v), 6) for v in np.nanmean(np.vstack(profiles), axis=0)]


def build_payload(rows: list[dict[str, str]], flank: int, top_n: int) -> dict:
    sample_payload = []
    x = list(range(-flank, 0)) + list(range(1, flank + 1))
    for row in rows:
        motifs = _discover_motifs(row["match_dir"])[:top_n]
        motif_payload = []
        for motif in motifs:
            prefix = str(motif["prefix"])
            bed = Path(row["match_dir"]) / prefix / "beds" / f"{prefix}_all.bed"
            centers = _read_bed_centers(bed)
            motif_payload.append({
                "prefix": prefix,
                "name": motif["name"],
                "score": motif["score"],
                "sites": len(centers),
                "profile": _mean_profile(row["signal"], centers, flank),
            })
        sample_payload.append({
            "sample": row["sample"],
            "label": row.get("label") or row["sample"],
            "condition": row.get("condition") or row["sample"],
            "motifs": motif_payload,
        })
    return {"x": x, "samples": sample_payload}


def write_html(payload: dict, output: str | Path, title: str) -> None:
    escaped_title = html.escape(title)
    data = json.dumps(payload)
    document = f"""<!doctype html>
<html lang=\"en\">
<head>
<meta charset=\"utf-8\">
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<title>{escaped_title}</title>
<style>
body {{ margin: 0; font-family: Arial, Helvetica, sans-serif; background: #f4f7fb; color: #152133; font-weight: 700; }}
.wrap {{ max-width: 1180px; margin: 24px auto; padding: 0 18px; }}
.panel {{ background: #fff; border: 1px solid #d9e2ec; border-radius: 10px; box-shadow: 0 8px 28px rgba(21,33,51,.08); overflow: hidden; }}
.head {{ padding: 18px 22px; border-bottom: 1px solid #e7edf4; }}
h1 {{ margin: 0; font-size: 24px; }}
.controls {{ display: flex; flex-wrap: wrap; gap: 12px; padding: 14px 22px; border-bottom: 1px solid #e7edf4; background: #fbfdff; }}
label {{ font-size: 12px; color: #52606d; text-transform: uppercase; letter-spacing: .06em; }}
select {{ display: block; min-width: 190px; padding: 7px; border: 1px solid #cbd5e1; border-radius: 7px; font-weight: 700; background: white; }}
.body {{ display: grid; grid-template-columns: 360px 1fr; gap: 0; }}
.left {{ padding: 18px; border-right: 1px solid #e7edf4; }}
.right {{ padding: 18px; }}
svg {{ width: 100%; height: auto; display: block; }}
.bar {{ cursor: pointer; fill: #1f77b4; opacity: .75; }}
.bar:hover {{ opacity: 1; }}
.axis {{ stroke: #64748b; }}
@media (max-width: 900px) {{ .body {{ grid-template-columns: 1fr; }} .left {{ border-right: 0; border-bottom: 1px solid #e7edf4; }} }}
</style>
</head>
<body>
<div class=\"wrap\"><div class=\"panel\">
<div class=\"head\"><h1>{escaped_title}</h1></div>
<div class=\"controls\"><div><label>Sample<select id=\"sample\"></select></label></div><div><label>Motif<select id=\"motif\"></select></label></div><div><label>Layout<select id=\"layout\"><option value=\"single\">Single panel</option><option value=\"wide\">Wide panel</option></select></label></div></div>
<div class=\"body\"><div class=\"left\"><svg id=\"bars\" viewBox=\"0 0 340 520\"></svg></div><div class=\"right\"><svg id=\"profile\" viewBox=\"0 0 720 520\"></svg></div></div>
</div></div>
<script>
const data = {data};
const sampleSel = document.getElementById('sample');
const motifSel = document.getElementById('motif');
const bars = document.getElementById('bars');
const profile = document.getElementById('profile');
sampleSel.innerHTML = data.samples.map((s, i) => `<option value="${{i}}">${{s.label}}</option>`).join('');
function esc(s) {{ return String(s).replace(/[&<>]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;'}}[c])); }}
function currentSample() {{ return data.samples[Number(sampleSel.value || 0)]; }}
function drawBars() {{
  const s = currentSample();
  motifSel.innerHTML = s.motifs.map((m, i) => `<option value="${{i}}">${{m.name}}</option>`).join('');
  const maxScore = Math.max(...s.motifs.map(m => Number(m.score) || 0), 1);
  let parts = [`<rect width="340" height="520" fill="#fff"/>`, `<text x="12" y="24" font-size="14" font-weight="700">Top motifs</text>`];
  s.motifs.forEach((m, i) => {{
    const y = 46 + i * 22;
    const w = 210 * ((Number(m.score) || 0) / maxScore);
    parts.push(`<rect class="bar" x="118" y="${{y-12}}" width="${{w.toFixed(1)}}" height="15" data-i="${{i}}"/>`);
    parts.push(`<text x="112" y="${{y}}" text-anchor="end" font-size="11" font-weight="700">${{esc(m.name).slice(0,18)}}</text>`);
    parts.push(`<text x="${{124+w}}" y="${{y}}" font-size="10" fill="#52606d">${{m.sites}}</text>`);
  }});
  bars.innerHTML = parts.join('');
  bars.querySelectorAll('.bar').forEach(el => el.addEventListener('click', () => {{ motifSel.value = el.dataset.i; drawProfile(); }}));
  drawProfile();
}}
function drawProfile() {{
  const s = currentSample();
  const m = s.motifs[Number(motifSel.value || 0)] || s.motifs[0];
  if (!m) {{ profile.innerHTML = '<text x="40" y="40">No motifs found</text>'; return; }}
  const x = data.x;
  const y = m.profile;
  const ymin = Math.min(...y, 0), ymax = Math.max(...y, 1e-9);
  const pad = (ymax - ymin || 1) * .12;
  const y0 = ymin - pad, y1 = ymax + pad;
  const sx = v => 70 + ((v - x[0]) / (x[x.length-1] - x[0] || 1)) * 590;
  const sy = v => 420 - ((v - y0) / (y1 - y0 || 1)) * 320;
  const d = y.map((v, i) => `${{i ? 'L' : 'M'}}${{sx(x[i]).toFixed(1)}},${{sy(v).toFixed(1)}}`).join(' ');
  profile.innerHTML = `<rect width="720" height="520" fill="#fff"/><text x="360" y="36" text-anchor="middle" font-size="18" font-weight="700">${{esc(s.label)}} - ${{esc(m.name)}} (${{m.sites}} sites)</text><line x1="70" y1="420" x2="660" y2="420" class="axis"/><line x1="70" y1="100" x2="70" y2="420" class="axis"/><line x1="${{sx(0).toFixed(1)}}" y1="100" x2="${{sx(0).toFixed(1)}}" y2="420" stroke="#94a3b8" stroke-dasharray="4 4"/><path d="${{d}}" fill="none" stroke="#d62728" stroke-width="2.2"/><text x="360" y="470" text-anchor="middle" font-size="13" font-weight="700">bp from motif center</text>`;
}}
sampleSel.addEventListener('change', drawBars);
motifSel.addEventListener('change', drawProfile);
drawBars();
</script>
</body>
</html>"""
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(document, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create an interactive aggregate HTML report from match-motifs outputs.")
    parser.add_argument("--manifest", required=True, help="TSV with sample, signal, and match_dir columns.")
    parser.add_argument("--output", required=True, help="Output self-contained HTML file.")
    parser.add_argument("--flank", type=int, default=100, help="Flank around motif centers for aggregate profiles (default: 100).")
    parser.add_argument("--top-n", type=int, default=30, help="Number of motifs to preload per sample (default: 30).")
    parser.add_argument("--title", default="Aggregate motif footprint browser")
    args = parser.parse_args(argv)
    rows = _read_manifest(args.manifest)
    payload = build_payload(rows, flank=max(1, args.flank), top_n=max(1, args.top_n))
    write_html(payload, args.output, args.title)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
