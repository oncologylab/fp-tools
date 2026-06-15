"""Interactive batch aggregate HTML reports for match-motifs and diff-footprints outputs."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyBigWig


DEFAULT_COLORS = ["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#f97316", "#0891b2", "#7c3aed", "#64748b"]


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
                motifs.append({
                    "prefix": prefix,
                    "name": row.get("name") or prefix,
                    "motif_id": row.get("motif_id") or "",
                    "score": total,
                    "sites": int(total),
                })
    if not motifs:
        for bed in sorted(root.glob("*/beds/*_all.bed")):
            prefix = bed.parent.parent.name
            centers = _read_bed_centers(bed)
            motifs.append({"prefix": prefix, "name": prefix, "motif_id": "", "score": len(centers), "sites": len(centers)})
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


def _x_values(flank: int) -> list[int]:
    return list(range(-flank, 0)) + list(range(1, flank + 1))


def _condition_colors(conditions: list[str]) -> dict[str, str]:
    return {cond: DEFAULT_COLORS[idx % len(DEFAULT_COLORS)] for idx, cond in enumerate(conditions)}


def _normalize_profile(profile: list[float], mode: str, baseline: float | None = None, scale: float | None = None) -> list[float]:
    arr = np.asarray(profile, dtype=float)
    if mode == "none" or arr.size == 0:
        return [round(float(v), 6) for v in arr]
    center = float(np.nanmedian(arr)) if baseline is None else baseline
    spread = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10)) if scale is None else scale
    if not np.isfinite(spread) or abs(spread) < 1e-12:
        spread = 1.0
    return [round(float(v), 6) for v in ((arr - center) / spread)]


def _profile_mean(profiles: list[list[float]]) -> list[float]:
    if not profiles:
        return []
    return [round(float(v), 6) for v in np.nanmean(np.asarray(profiles, dtype=float), axis=0)]


def build_payload(rows: list[dict[str, str]], flank: int, top_n: int, normalization: str = "none") -> dict:
    """Build the compact batch aggregate payload from manifest rows."""

    x = _x_values(flank)
    normalization = (normalization or "none").replace("_", "-")
    sample_rows = []
    motif_meta: dict[str, dict[str, str | int | float]] = {}
    motif_scores: dict[str, float] = defaultdict(float)

    for idx, row in enumerate(rows):
        sample = row["sample"]
        condition = row.get("condition") or sample
        label = row.get("label") or sample
        sample_rows.append({"sample": sample, "condition": condition, "label": label, "row": row, "idx": idx})
        for motif in _discover_motifs(row["match_dir"]):
            prefix = str(motif["prefix"])
            motif_meta.setdefault(prefix, motif)
            motif_scores[prefix] = max(motif_scores[prefix], float(motif.get("score") or 0.0))

    selected_prefixes = sorted(motif_scores, key=lambda pfx: (-motif_scores[pfx], str(motif_meta[pfx].get("name", pfx))))[:top_n]
    raw_profiles: dict[tuple[str, str], list[float]] = {}
    profile_values_for_norm = []
    for sample_info in sample_rows:
        row = sample_info["row"]
        sample = sample_info["sample"]
        for prefix in selected_prefixes:
            bed = Path(row["match_dir"]) / prefix / "beds" / f"{prefix}_all.bed"
            centers = _read_bed_centers(bed)
            profile = _mean_profile(row["signal"], centers, flank)
            raw_profiles[(sample, prefix)] = profile
            profile_values_for_norm.extend(profile)

    global_baseline = global_scale = None
    if normalization == "sample-quantile" and profile_values_for_norm:
        arr = np.asarray(profile_values_for_norm, dtype=float)
        global_baseline = float(np.nanmedian(arr))
        global_scale = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10))

    conditions = []
    for sample_info in sample_rows:
        if sample_info["condition"] not in conditions:
            conditions.append(sample_info["condition"])
    colors = _condition_colors(conditions)

    motifs = []
    for prefix in selected_prefixes:
        meta = motif_meta[prefix]
        series = []
        condition_profiles: dict[str, list[list[float]]] = defaultdict(list)
        for sample_info in sample_rows:
            sample = sample_info["sample"]
            condition = sample_info["condition"]
            profile = raw_profiles.get((sample, prefix), [0.0] * len(x))
            if normalization == "condition-quantile":
                cond_values = []
                for other in sample_rows:
                    if other["condition"] == condition:
                        cond_values.extend(raw_profiles.get((other["sample"], prefix), []))
                arr = np.asarray(cond_values, dtype=float) if cond_values else np.asarray(profile, dtype=float)
                baseline = float(np.nanmedian(arr))
                scale = float(np.nanpercentile(arr, 90) - np.nanpercentile(arr, 10))
                profile = _normalize_profile(profile, normalization, baseline, scale)
            else:
                profile = _normalize_profile(profile, normalization, global_baseline, global_scale)
            condition_profiles[condition].append(profile)
            series.append({"id": f"sample::{sample}", "label": sample_info["label"], "kind": "sample", "condition": condition, "profile": profile})
        for condition, profiles in condition_profiles.items():
            series.append({"id": f"condition::{condition}", "label": f"{condition} mean", "kind": "condition", "condition": condition, "profile": _profile_mean(profiles)})
        motifs.append({"prefix": prefix, "name": str(meta.get("name") or prefix), "motif_id": str(meta.get("motif_id") or ""), "score": round(float(motif_scores[prefix]), 6), "sites": int(meta.get("sites") or 0), "series": series})

    return {"schema": "fp-tools.aggregate.batch.v2", "x": x, "motifs": motifs, "conditions": conditions, "colors": colors, "normalization": normalization, "x_label": "Distance from motif center (bp)", "y_label": "Corrected cut-site signal (a.u.)" if normalization == "none" else "Normalized corrected cut-site signal (a.u.)"}


def _compressed_json_b64(payload: dict) -> str:
    text = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    return base64.b64encode(gzip.compress(text.encode("utf-8"), compresslevel=9)).decode("ascii")


def _decode_payload_b64(payload_b64: str) -> dict:
    return json.loads(gzip.decompress(base64.b64decode(payload_b64)).decode("utf-8"))


def read_embedded_payload(path: str | Path) -> dict:
    text = Path(path).read_text(encoding="utf-8")
    match = re.search(r'reportPayloadB64="([^"]+)"', text)
    if not match:
        match = re.search(r"const\s+reportPayloadB64\s*=\s*['\"]([^'\"]+)['\"]", text)
    if not match:
        raise ValueError(f"Could not find reportPayloadB64 in {path}")
    return _decode_payload_b64(match.group(1))


def _series_from_diff_payload(payload: dict, source_label: str) -> dict:
    aggregate = payload.get("aggregate") or {}
    motifs = []
    conditions_seen = []
    colors = dict(payload.get("colors") or {})
    for motif in aggregate.get("motifs") or []:
        series = []
        for cond in motif.get("conditions") or []:
            cond_name = str(cond.get("name") or "condition")
            if cond_name not in conditions_seen:
                conditions_seen.append(cond_name)
            if cond_name not in colors:
                colors[cond_name] = DEFAULT_COLORS[len(colors) % len(DEFAULT_COLORS)]
            for sample in cond.get("samples") or []:
                label = str(sample.get("name") or f"{source_label} {cond_name}")
                series.append({"id": f"{source_label}::sample::{label}", "label": label, "kind": "sample", "condition": cond_name, "profile": sample.get("profile") or []})
            if cond.get("profile") is not None:
                series.append({"id": f"{source_label}::condition::{cond_name}", "label": f"{source_label} {cond_name} mean" if source_label else f"{cond_name} mean", "kind": "condition", "condition": cond_name, "profile": cond.get("profile") or []})
        motifs.append({"prefix": str(motif.get("prefix") or motif.get("name") or "motif"), "name": str(motif.get("name") or motif.get("prefix") or "motif"), "motif_id": str(motif.get("motif_id") or ""), "score": abs(float(motif.get("change") or 0.0)), "sites": int(motif.get("n_sites") or motif.get("sites") or 0), "series": series})
    return {"schema": "fp-tools.aggregate.batch.v2", "x": aggregate.get("x") or [], "motifs": motifs, "conditions": conditions_seen, "colors": {cond: colors.get(cond) or DEFAULT_COLORS[idx % len(DEFAULT_COLORS)] for idx, cond in enumerate(conditions_seen)}, "normalization": aggregate.get("normalization") or payload.get("normalization") or "none", "x_label": aggregate.get("x_label") or "Distance from motif center (bp)", "y_label": aggregate.get("y_label") or "Corrected cut-site signal (a.u.)"}


def _ensure_batch_payload(payload: dict, source_label: str = "") -> dict:
    if payload.get("schema") == "fp-tools.aggregate.batch.v2":
        return payload
    if "aggregate" in payload:
        return _series_from_diff_payload(payload, source_label)
    if "samples" in payload:
        motifs_by_prefix: dict[str, dict] = {}
        conditions = []
        for sample in payload.get("samples") or []:
            condition = str(sample.get("condition") or sample.get("sample") or "sample")
            if condition not in conditions:
                conditions.append(condition)
            for motif in sample.get("motifs") or []:
                prefix = str(motif.get("prefix") or motif.get("name") or "motif")
                entry = motifs_by_prefix.setdefault(prefix, {"prefix": prefix, "name": str(motif.get("name") or prefix), "motif_id": str(motif.get("motif_id") or ""), "score": float(motif.get("score") or 0.0), "sites": int(motif.get("sites") or 0), "series": []})
                entry["series"].append({"id": f"sample::{sample.get('sample')}", "label": str(sample.get("label") or sample.get("sample") or condition), "kind": "sample", "condition": condition, "profile": motif.get("profile") or []})
        return {"schema": "fp-tools.aggregate.batch.v2", "x": payload.get("x") or [], "motifs": list(motifs_by_prefix.values()), "conditions": conditions, "colors": _condition_colors(conditions), "normalization": payload.get("normalization") or "none", "x_label": payload.get("x_label") or "Distance from motif center (bp)", "y_label": payload.get("y_label") or "Corrected cut-site signal (a.u.)"}
    raise ValueError("Unsupported aggregate HTML payload schema")


def merge_payloads(payloads: list[dict]) -> dict:
    if not payloads:
        raise ValueError("No aggregate payloads were provided")
    normalized = [_ensure_batch_payload(payload, f"report{idx + 1}") for idx, payload in enumerate(payloads)]
    merged = {"schema": "fp-tools.aggregate.batch.v2", "x": normalized[0].get("x") or [], "motifs": [], "conditions": [], "colors": {}, "normalization": ", ".join(sorted({str(p.get("normalization") or "none") for p in normalized})), "x_label": normalized[0].get("x_label") or "Distance from motif center (bp)", "y_label": normalized[0].get("y_label") or "Corrected cut-site signal (a.u.)"}
    motifs_by_prefix: dict[str, dict] = {}
    for payload in normalized:
        for cond in payload.get("conditions") or []:
            if cond not in merged["conditions"]:
                merged["conditions"].append(cond)
        merged["colors"].update(payload.get("colors") or {})
        if not merged["x"] and payload.get("x"):
            merged["x"] = payload["x"]
        for motif in payload.get("motifs") or []:
            prefix = str(motif.get("prefix") or motif.get("name") or "motif")
            entry = motifs_by_prefix.setdefault(prefix, {"prefix": prefix, "name": str(motif.get("name") or prefix), "motif_id": str(motif.get("motif_id") or ""), "score": float(motif.get("score") or 0.0), "sites": int(motif.get("sites") or 0), "series": []})
            entry["score"] = max(float(entry.get("score") or 0.0), float(motif.get("score") or 0.0))
            entry["sites"] = max(int(entry.get("sites") or 0), int(motif.get("sites") or 0))
            entry["series"].extend(motif.get("series") or [])
    if not merged["colors"]:
        merged["colors"] = _condition_colors(merged["conditions"])
    else:
        for idx, cond in enumerate(merged["conditions"]):
            merged["colors"].setdefault(cond, DEFAULT_COLORS[idx % len(DEFAULT_COLORS)])
    merged["motifs"] = sorted(motifs_by_prefix.values(), key=lambda m: (str(m.get("name") or ""), str(m.get("prefix") or "")))
    return merged


def write_html(payload: dict, output: str | Path, title: str, default_layout: str = "2x2") -> None:
    payload = _ensure_batch_payload(payload)
    payload["default_layout"] = default_layout
    escaped_title = html.escape(title)
    payload_b64 = _compressed_json_b64(payload)
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>{escaped_title}</title><style>
:root{{--ink:#152133;--muted:#596579;--line:#d9e2ec;--grid:#e8eef5;--panel:#fff;--bg:#eef3f8;--accent:#173b73;--soft:#f7fafc}}*{{box-sizing:border-box}}body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--ink);font-weight:700}}.wrap{{max-width:min(1840px,calc(100vw - 20px));margin:6px auto;padding:0 6px}}.panel{{background:#fff;border:1px solid var(--line);border-radius:8px;box-shadow:0 14px 34px rgba(21,33,51,.10);overflow:hidden}}.head{{padding:8px 16px 6px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#fff 0%,#f7fafc 100%)}}h1{{margin:0;font-size:20px;line-height:1.1;font-weight:900}}.sub{{margin:2px 0 0;color:var(--muted);font-size:11px}}.top-row{{display:grid;grid-template-columns:280px minmax(300px,1fr) 330px 210px;gap:8px;padding:7px 10px;border-bottom:1px solid var(--line);background:#fbfdff}}.card{{border:1px solid var(--line);border-radius:7px;background:#fff;padding:6px;min-height:70px}}.section-title{{font-size:10px;line-height:1.05;text-transform:uppercase;letter-spacing:.08em;color:#728197;margin:0 0 4px;font-weight:900}}.controls{{display:flex;flex-wrap:wrap;gap:5px}}label{{font-size:10px;color:#52606d;text-transform:uppercase;letter-spacing:.06em;font-weight:900}}select,input{{border:1px solid #cbd5e1;border-radius:6px;background:white;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:800;padding:5px 7px}}button{{border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:5px 7px;font-family:Arial,Helvetica,sans-serif;font-size:10px;font-weight:900;cursor:pointer}}button:hover{{background:#f2f6fb}}.color-row{{display:flex;align-items:center;justify-content:space-between;gap:6px;font-size:11px;color:#334e68;font-weight:800;border:1px solid #e6edf5;border-radius:999px;padding:3px 6px;background:#fbfdff}}.color-row input{{width:26px;height:18px;padding:0}}.layout-controls{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:5px;padding:7px 10px;border-bottom:1px solid var(--line);background:#fff}}.slot-control{{display:grid;gap:3px}}.grid{{display:grid;gap:8px;padding:8px;background:#f8fbff}}.grid.g1x1{{grid-template-columns:1fr}}.grid.g1x2{{grid-template-columns:1fr 1fr}}.grid.g2x2{{grid-template-columns:1fr 1fr}}.grid.g2x3{{grid-template-columns:1fr 1fr 1fr}}.plot-card{{background:#fff;border:1px solid var(--line);border-radius:7px;padding:5px;min-width:0}}svg{{width:100%;height:auto;display:block;background:#fff}}.combo{{position:relative}}.combo input{{width:100%}}.combo-list{{display:none;position:absolute;z-index:20;left:0;right:0;top:29px;max-height:260px;overflow:auto;border:1px solid var(--line);border-radius:6px;background:#fff;box-shadow:0 10px 24px rgba(21,33,51,.16)}}.combo-list.open{{display:block}}.combo-option{{padding:7px 10px;font-size:13px;font-weight:800;cursor:pointer;border-bottom:1px solid #eef3f8}}.combo-option small{{display:block;color:#728197;font-size:11px;font-weight:700}}.combo-option:hover,.combo-option.active{{background:#edf4ff}}.axis{{stroke:#3b4552;stroke-width:1.35}}.grid-line{{stroke:var(--grid);stroke-width:1}}.zero{{stroke:#677386;stroke-width:1.4;stroke-dasharray:4 4}}.tick{{font-family:Arial,Helvetica,sans-serif;font-size:11px;fill:var(--muted);font-weight:800}}.axis-label{{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:var(--ink);font-weight:900}}.plot-title{{font-family:Arial,Helvetica,sans-serif;font-size:14px;fill:var(--ink);font-weight:900}}@media(max-width:1100px){{.top-row{{grid-template-columns:1fr 1fr}}.layout-controls{{grid-template-columns:1fr 1fr}}.grid.g1x2,.grid.g2x2,.grid.g2x3{{grid-template-columns:1fr}}}}@media(max-width:760px){{.top-row,.layout-controls{{grid-template-columns:1fr}}}}
</style></head><body><div class="wrap"><div class="panel"><div class="head"><h1>{escaped_title}</h1><p class="sub">Standalone multi-sample, multi-TF aggregate report</p></div><div class="top-row"><div><p class="section-title">Groups</p><div class="card controls" id="color-controls"></div></div><div><p class="section-title">Selected motif</p><div class="card"><div class="combo"><input id="motif-search" type="text" autocomplete="off" placeholder="Search motif"><div id="motif-options" class="combo-list"></div></div><p id="motif-detail" class="sub">Loading report</p></div></div><div><p class="section-title">Layout</p><div class="card controls"><label>Grid<select id="layout"><option value="1x1">1x1</option><option value="1x2">1x2</option><option value="2x2">2x2</option><option value="2x3">2x3</option></select></label></div></div><div><p class="section-title">Export editable SVG</p><div class="card controls"><button id="download-grid">Download grid SVG</button><button id="download-panel">Download selected panel SVG</button></div></div></div><div id="slot-controls" class="layout-controls"></div><div id="plot-grid" class="grid g2x2"></div></div></div><script>
const DEFAULT_COLORS={json.dumps(DEFAULT_COLORS)};const reportPayloadB64="{payload_b64}";let payload=null,selectedPrefix=null,activeOptionIndex=0,slotChoices=[],selectedPanel=0;const motifSearch=document.getElementById('motif-search'),motifOptions=document.getElementById('motif-options'),motifDetail=document.getElementById('motif-detail'),colorControls=document.getElementById('color-controls'),layoutSel=document.getElementById('layout'),slotControls=document.getElementById('slot-controls'),plotGrid=document.getElementById('plot-grid');
function escText(v){{return String(v??'').replace(/[&<>\"']/g,c=>({{'&':'&amp;','<':'&lt;','>':'&gt;','\"':'&quot;',"'":'&#39;'}}[c]))}}function b64ToBytes(b64){{return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}}async function decodePayload(){{if(!('DecompressionStream'in window))throw new Error('This report needs a modern browser with gzip DecompressionStream support.');const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}}function motifLabel(m){{return m?(m.motif_id?`${{m.name}} (${{m.motif_id}})`:m.name):''}}function currentMotif(){{return payload.motifs.find(m=>m.prefix===selectedPrefix)||payload.motifs[0]}}function currentColors(){{const out={{...payload.colors}};document.querySelectorAll('[data-color-condition]').forEach(inp=>out[inp.dataset.colorCondition]=inp.value);return out}}function niceTicks(min,max,n){{const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}}function fmt(v){{return Math.abs(v)>=1?v.toFixed(1).replace('-0.0','0.0'):v.toFixed(2).replace('-0.00','0.00')}}
function renderColors(){{colorControls.innerHTML=(payload.conditions||[]).map(c=>`<label class="color-row"><span>${{escText(c)}}</span><input type="color" data-color-condition="${{escText(c)}}" value="${{payload.colors[c]||'#64748b'}}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',renderAll))}}function filteredMotifs(){{const q=motifSearch.value.trim().toLowerCase();return [...payload.motifs].sort((a,b)=>(a.name||'').localeCompare(b.name||'')||(a.prefix||'').localeCompare(b.prefix||'')).filter(m=>!q||(`${{m.name}} ${{m.motif_id||''}} ${{m.prefix}}`).toLowerCase().includes(q)).slice(0,100)}}function renderMotifOptions(){{const list=filteredMotifs();activeOptionIndex=Math.min(activeOptionIndex,Math.max(0,list.length-1));motifOptions.innerHTML=list.map((m,i)=>`<div class="combo-option${{i===activeOptionIndex?' active':''}}" data-prefix="${{escText(m.prefix)}}"><span>${{escText(motifLabel(m))}}</span><small>${{escText(m.prefix)}} · ${{m.sites||0}} sites</small></div>`).join('');motifOptions.classList.add('open');motifOptions.querySelectorAll('.combo-option').forEach(el=>el.addEventListener('mousedown',ev=>{{ev.preventDefault();setMotif(el.dataset.prefix);motifOptions.classList.remove('open')}}))}}function setMotif(prefix){{selectedPrefix=prefix;const m=currentMotif();motifSearch.value=motifLabel(m);motifDetail.textContent=`${{motifLabel(m)}} · ${{m.sites||0}} sites · ${{m.series.length}} available series`;renderAll()}}
function availableChoices(motif){{const choices=[{{id:'all_conditions',label:'All condition means'}},{{id:'all_samples',label:'All samples'}}];const conds=[...new Set(motif.series.map(s=>s.condition))];conds.forEach(c=>choices.push({{id:`condition_mean::${{c}}`,label:`${{c}} mean`}},{{id:`condition_samples::${{c}}`,label:`${{c}} samples`}}));motif.series.filter(s=>s.kind==='sample').forEach(s=>choices.push({{id:s.id,label:s.label}}));return choices}}function slotCount(layout){{return layout==='1x1'?1:layout==='1x2'?2:layout==='2x2'?4:6}}function ensureSlots(){{const n=slotCount(layoutSel.value);while(slotChoices.length<n)slotChoices.push(slotChoices.length===0?'all_conditions':'all_samples');slotChoices=slotChoices.slice(0,n)}}function renderSlotControls(){{ensureSlots();const choices=availableChoices(currentMotif());slotControls.innerHTML=slotChoices.map((choice,idx)=>`<div class="slot-control"><label>Panel ${{idx+1}}<select data-slot="${{idx}}">${{choices.map(c=>`<option value="${{escText(c.id)}}" ${{c.id===choice?'selected':''}}>${{escText(c.label)}}</option>`).join('')}}</select></label></div>`).join('');slotControls.querySelectorAll('select').forEach(sel=>sel.addEventListener('change',()=>{{slotChoices[Number(sel.dataset.slot)]=sel.value;renderPlots()}}))}}function seriesForChoice(motif,choice){{if(choice==='all_conditions')return motif.series.filter(s=>s.kind==='condition');if(choice==='all_samples')return motif.series.filter(s=>s.kind==='sample');if(choice.startsWith('condition_mean::')){{const c=choice.split('::')[1];return motif.series.filter(s=>s.kind==='condition'&&s.condition===c)}}if(choice.startsWith('condition_samples::')){{const c=choice.split('::')[1];return motif.series.filter(s=>s.kind==='sample'&&s.condition===c)}}return motif.series.filter(s=>s.id===choice)}}function pathD(profile,x,sx,sy){{return profile.map((y,i)=>`${{i?'L':'M'}}${{sx(x[i]).toFixed(2)}},${{sy(y).toFixed(2)}}`).join(' ')}}function drawPanel(motif,choice,idx){{const series=seriesForChoice(motif,choice),x=payload.x,width=760,height=320,margin={{top:38,right:22,bottom:50,left:78}},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,colors=currentColors();const allY=series.flatMap(s=>s.profile).filter(Number.isFinite);let ymin=Math.min(...allY,0),ymax=Math.max(...allY,1e-9);const pad=Math.max((ymax-ymin||1)*.18,1e-6);ymin-=pad;ymax+=pad;const sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH;const xTicks=[x[0],Math.round(x[0]/2),0,Math.round(x[x.length-1]/2),x[x.length-1]],yTicks=niceTicks(ymin,ymax,5);let parts=[`<svg class="aggregate-panel" data-panel="${{idx}}" viewBox="0 0 ${{width}} ${{height}}"><rect width="${{width}}" height="${{height}}" fill="#fff"/><text x="${{width/2}}" y="20" text-anchor="middle" class="plot-title">${{escText(motifLabel(motif))}}</text><text x="${{width/2}}" y="35" text-anchor="middle" class="tick">${{escText(choice.replace('::',' '))}}</text>`];yTicks.forEach(v=>parts.push(`<line x1="${{margin.left}}" y1="${{sy(v)}}" x2="${{margin.left+innerW}}" y2="${{sy(v)}}" class="grid-line"/><text x="${{margin.left-9}}" y="${{sy(v)+4}}" class="tick" text-anchor="end">${{fmt(v)}}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${{sx(v)}}" y1="${{margin.top}}" x2="${{sx(v)}}" y2="${{margin.top+innerH}}" class="grid-line"/><text x="${{sx(v)}}" y="${{margin.top+innerH+22}}" class="tick" text-anchor="middle">${{v}}</text>`));parts.push(`<line x1="${{sx(0)}}" y1="${{margin.top}}" x2="${{sx(0)}}" y2="${{margin.top+innerH}}" class="zero"/><line x1="${{margin.left}}" y1="${{margin.top+innerH}}" x2="${{margin.left+innerW}}" y2="${{margin.top+innerH}}" class="axis"/><line x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top+innerH}}" class="axis"/>`);series.forEach((s,i)=>{{const color=colors[s.condition]||DEFAULT_COLORS[i%DEFAULT_COLORS.length],wide=s.kind==='condition';parts.push(`<path d="${{pathD(s.profile,x,sx,sy)}}" fill="none" stroke="${{color}}" stroke-width="${{wide?2.2:1.0}}" stroke-opacity="${{wide?0.95:0.35}}"><title>${{escText(s.label)}}</title></path>`);if(i<8)parts.push(`<text x="${{margin.left+8}}" y="${{margin.top+14+i*13}}" font-family="Arial,Helvetica,sans-serif" font-size="10" font-weight="900" fill="${{color}}">${{escText(s.label).slice(0,34)}}</text>`)}});parts.push(`<text x="${{margin.left+innerW/2}}" y="${{height-13}}" class="axis-label" text-anchor="middle">${{escText(payload.x_label)}}</text><text x="20" y="${{margin.top+innerH/2}}" class="axis-label" text-anchor="middle" transform="rotate(-90 20 ${{margin.top+innerH/2}})">${{escText(payload.y_label)}}</text></svg>`);return parts.join('')}}function renderPlots(){{const motif=currentMotif();plotGrid.className=`grid g${{layoutSel.value}}`;plotGrid.innerHTML=slotChoices.map((choice,idx)=>`<div class="plot-card" data-card="${{idx}}">${{drawPanel(motif,choice,idx)}}</div>`).join('');plotGrid.querySelectorAll('.plot-card').forEach(card=>card.addEventListener('click',()=>{{selectedPanel=Number(card.dataset.card)}}))}}function renderAll(){{renderSlotControls();renderPlots()}}function svgBlob(svgNode){{const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');return new Blob([new XMLSerializer().serializeToString(clone)],{{type:'image/svg+xml;charset=utf-8'}})}}function downloadBlob(blob,filename){{const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}}function downloadGrid(){{const svgs=[...document.querySelectorAll('.aggregate-panel')];if(!svgs.length)return;const w=760,h=320,n=slotCount(layoutSel.value),cols=layoutSel.value==='2x3'?3:(layoutSel.value==='1x1'?1:2),rows=Math.ceil(n/cols);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${{cols*w}} ${{rows*h}}">`];svgs.forEach((svg,i)=>{{parts.push(`<g transform="translate(${{(i%cols)*w}},${{Math.floor(i/cols)*h}})">${{svg.innerHTML}}</g>`)}});parts.push('</svg>');downloadBlob(new Blob(parts,{{type:'image/svg+xml;charset=utf-8'}}),'plot_aggregate_batch_grid.svg')}}document.getElementById('download-grid').addEventListener('click',downloadGrid);document.getElementById('download-panel').addEventListener('click',()=>{{const svg=document.querySelector(`.aggregate-panel[data-panel="${{selectedPanel}}"]`)||document.querySelector('.aggregate-panel');if(svg)downloadBlob(svgBlob(svg),'plot_aggregate_batch_panel.svg')}});layoutSel.addEventListener('change',()=>{{ensureSlots();renderAll()}});motifSearch.addEventListener('focus',()=>{{activeOptionIndex=0;renderMotifOptions()}});motifSearch.addEventListener('input',()=>{{activeOptionIndex=0;renderMotifOptions()}});motifSearch.addEventListener('keydown',ev=>{{const list=filteredMotifs();if(ev.key==='ArrowDown'){{ev.preventDefault();activeOptionIndex=Math.min(activeOptionIndex+1,Math.max(0,list.length-1));renderMotifOptions()}}else if(ev.key==='ArrowUp'){{ev.preventDefault();activeOptionIndex=Math.max(activeOptionIndex-1,0);renderMotifOptions()}}else if(ev.key==='Enter'){{ev.preventDefault();if(list[activeOptionIndex]){{setMotif(list[activeOptionIndex].prefix);motifOptions.classList.remove('open')}}}}else if(ev.key==='Escape')motifOptions.classList.remove('open')}});document.addEventListener('click',ev=>{{if(!motifSearch.parentElement.contains(ev.target))motifOptions.classList.remove('open')}});decodePayload().then(data=>{{payload=data;layoutSel.value=payload.default_layout||'2x2';renderColors();if(payload.motifs&&payload.motifs.length)setMotif(payload.motifs[0].prefix);else motifDetail.textContent='No motifs found'}}).catch(err=>{{motifDetail.textContent=`Could not open report payload: ${{err.message}}`}});
</script></body></html>"""
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(document, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create an interactive aggregate HTML report from match-motifs or embedded diff-footprints outputs.")
    parser.add_argument("--manifest", help="TSV with sample, signal, and match_dir columns.")
    parser.add_argument("--input-html", nargs="*", default=[], help="Existing aggregate/diff-footprints HTML report(s) with embedded reportPayloadB64 payloads.")
    parser.add_argument("--output", required=True, help="Output self-contained HTML file.")
    parser.add_argument("--flank", type=int, default=100, help="Flank around motif centers for aggregate profiles (default: 100).")
    parser.add_argument("--top-n", type=int, default=30, help="Number of motifs to preload from manifest mode (default: 30).")
    parser.add_argument("--normalization", choices=["none", "sample-quantile", "condition-quantile"], default="none", help="Profile scaling for manifest mode (default: none).")
    parser.add_argument("--default-layout", choices=["1x1", "1x2", "2x2", "2x3"], default="2x2", help="Initial panel grid layout (default: 2x2).")
    parser.add_argument("--title", default="Aggregate motif footprint browser")
    args = parser.parse_args(argv)
    payloads = []
    if args.manifest:
        payloads.append(build_payload(_read_manifest(args.manifest), flank=max(1, args.flank), top_n=max(1, args.top_n), normalization=args.normalization))
    for path in args.input_html:
        payloads.append(read_embedded_payload(path))
    if not payloads:
        parser.error("provide --manifest and/or --input-html")
    payload = merge_payloads(payloads)
    write_html(payload, args.output, args.title, default_layout=args.default_layout)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
