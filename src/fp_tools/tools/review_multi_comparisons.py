"""Review multiple diff-footprints HTML reports in one interactive page."""

from __future__ import annotations

import argparse
import base64
import copy
import gzip
import html
import json
import re
from pathlib import Path
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
    payload_b64 = _compressed_json_b64(review_payload)
    escaped_title = html.escape(str(review_payload.get("title") or "Review multiple differential footprint comparisons"))
    template = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>__TITLE__</title><style>
:root{--bg:#f5f8fc;--panel:#fff;--ink:#172033;--muted:#64748b;--border:#d8e2ef;--grid:#e7edf5}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:700}.wrap{max-width:2600px;margin:0 auto;padding:5px}.sub{margin:2px 0 0;color:var(--muted);font-weight:800}.sr-only{position:absolute;left:-10000px;top:auto;width:1px;height:1px;overflow:hidden}.board{display:grid;grid-template-columns:390px minmax(0,1fr);grid-template-rows:auto auto;gap:5px}.side{grid-column:1;grid-row:1/3;display:grid;align-content:start;gap:5px}.plots{grid-column:2;grid-row:1}.aggregate-card{grid-column:2;grid-row:2;overflow-x:auto}.card,.comparison-card,.aggregate-card{background:var(--panel);border:1px solid var(--border);border-radius:8px;padding:5px}.section-title{margin:0 0 3px;font-size:10px;text-transform:uppercase;letter-spacing:.04em;color:#46566c;font-weight:900}.export-stack{display:grid;gap:4px}.rows-control{display:grid;grid-template-columns:42px minmax(0,1fr) 48px;gap:4px;align-items:center;margin-top:4px}.rows-control input[type=number],.rows-control select{width:48px}.check-control{display:flex;align-items:center;gap:6px;margin-top:4px;font-size:10px;font-weight:900;color:#26364d}.check-control input{width:14px;height:14px}button,select,input{height:22px;border:1px solid #b9c8da;border-radius:5px;background:#fff;color:#26364d;font-family:Arial,Helvetica,sans-serif;font-size:10px;font-weight:900}button{padding:0 7px;cursor:pointer}.motif-select{width:100%;margin-bottom:4px}.motif-logo{height:132px;display:flex;align-items:center;justify-content:center;overflow:hidden;border:1px solid #e1e8f0;border-radius:7px;background:#fff}.motif-logo svg,.motif-logo img{max-width:100%;max-height:126px}.detail-grid p{margin:2px 0}.sample-style-panel{display:grid;gap:5px;max-height:470px;overflow:auto}.sample-style-row{display:grid;grid-template-columns:18px minmax(138px,1fr) 34px 42px 42px 54px;gap:4px;align-items:center;min-height:22px}.sample-style-name{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.sample-style-row input[type=color]{width:34px;padding:1px}.sample-style-row input[type=number],.sample-style-row select{width:100%;min-width:0}.comparison-grid{display:grid;grid-template-columns:repeat(var(--comparison-cols,2),minmax(0,1fr));gap:5px}.comparison-card{min-width:0}.comparison-head{display:flex;align-items:center;gap:5px;margin-bottom:2px}.comparison-head select{flex:1;min-width:0}.pair{display:grid;grid-template-columns:minmax(250px,.9fr) minmax(300px,1.1fr);gap:4px;align-items:start}.comparison-grid.compact-panels .pair{grid-template-columns:1fr}.comparison-grid.compact-panels .plot-box svg{max-height:255px}.plot-box{min-width:0;overflow:hidden}.plot-box svg,.aggregate-plot svg{width:100%;height:auto;display:block;background:#fff}.aggregate-head{display:flex;align-items:center;justify-content:space-between;gap:5px;margin-bottom:2px}.aggregate-grid{display:grid;grid-template-columns:repeat(var(--aggregate-cols,4),max-content);gap:8px;justify-content:start}.aggregate-tile{border:1px solid #e1e8f0;border-radius:7px;background:#fff;min-width:0;display:grid;grid-template-columns:180px max-content;gap:4px 7px;align-items:start;padding:3px;width:max-content;max-width:100%}.aggregate-grid.hide-legends .aggregate-tile{grid-template-columns:180px}.aggregate-grid.hide-legends .aggregate-legend-mini{display:none}.aggregate-tile-label{grid-column:1/-1;font-size:9px;font-weight:900;color:#172033;line-height:1.1;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}.aggregate-plot{width:180px;aspect-ratio:1/1}.aggregate-plot svg{height:100%;aspect-ratio:1/1}.aggregate-legend-mini{display:grid;gap:3px;align-content:start;padding:3px 2px 0 0;background:#fff;overflow:visible}.agg-legend-row{display:grid;grid-template-columns:24px max-content;gap:4px;align-items:center;font-family:Arial,Helvetica,sans-serif;font-size:8.5px;font-weight:900;color:#000;line-height:1.15}.agg-legend-row span{white-space:nowrap;color:#000}.agg-legend-line{width:22px;height:0;border-top-style:solid}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:var(--grid);stroke-width:1}.tick{font-size:10px;fill:#526176;font-weight:800}.axis-label{font-size:11px;fill:#243247;font-weight:900}.plot-title{font-size:13px;font-weight:900;fill:#172033}.pt,.rank-bar,.rank-name{cursor:pointer}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.30))}.rank-bar.selected{stroke:#111827;stroke-width:1.5}@media(max-width:1450px){.board{grid-template-columns:340px minmax(0,1fr)}.comparison-grid{grid-template-columns:1fr}.pair,.comparison-grid.compact-panels .pair{grid-template-columns:minmax(250px,.9fr) minmax(300px,1.1fr)}.comparison-grid.compact-panels .plot-box svg{max-height:none}}@media(max-width:900px){.board{display:block}.side,.plots,.aggregate-card{margin-bottom:8px}.comparison-grid,.pair,.comparison-grid.compact-panels .pair{display:block}.comparison-card,.aggregate-tile{margin-bottom:8px}}@media(max-width:520px){.aggregate-plot{width:180px}}
</style></head><body><div class="wrap"><p class="sr-only" id="report-detail">Loading report</p><main class="board"><aside class="side"><section class="card"><p class="section-title">Sample line styles</p><div id="sample-style-panel" class="sample-style-panel"></div></section><section class="card"><p class="section-title">Export editable SVG</p><div class="export-stack"><button id="download-logo">Download motif logo panel</button><button id="download-rank">Download bar plot panel</button><button id="download-volcano">Download volcano plot panel</button><button id="download-aggregate">Download motif aggregate panel</button><button id="download-panel">Download combined panel</button></div><label class="rows-control">Top rows <input id="rank-rows-slider" type="range" min="2" max="200" step="2" value="20"><input id="rank-rows" type="number" min="2" max="200" step="2" value="20"></label><label class="rows-control">Panels <span></span><select id="panel-count" aria-label="Number of comparison panels"><option value="4">4</option><option value="5">5</option><option value="6">6</option><option value="7">7</option><option value="8">8</option></select></label><label class="check-control"><input id="aggregate-legends" type="checkbox">Aggregate legends</label></section><section class="card"><p class="section-title">Selected motif</p><select id="motif-select" class="motif-select"></select><div id="motif-logo" class="motif-logo"></div></section></aside><section id="comparison-grid" class="plots comparison-grid"></section><section class="aggregate-card"><div class="aggregate-head"><p class="section-title">Motif aggregate review</p><span class="sub">Group autoscale</span></div><div id="aggregate-grid" class="aggregate-grid"></div></section></main></div><script>
const reportPayloadB64="__PAYLOAD__",aggregateDisplayBp=60,initialDisplayPanels=__DISPLAY_PANELS__,initialAggregateLegends="__AGGREGATE_LEGENDS__",plotSvgStyle='svg,text{font-family:Arial,Helvetica,sans-serif}.axis{stroke:#000;stroke-width:1.2}.zero{stroke:#555;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e7edf5;stroke-width:1}.tick{font-size:10px;fill:#000;font-weight:900}.axis-label{font-size:11px;fill:#000;font-weight:900}.plot-title{font-size:13px;font-weight:900;fill:#000}';let review=null,slotComparisons=[],activePrefix=null,sampleLineStyles={},aggregateDomain=null,showAggregateLegends=initialAggregateLegends!=="hide";const comparisonGrid=document.getElementById('comparison-grid'),aggregateGrid=document.getElementById('aggregate-grid'),sampleStylePanel=document.getElementById('sample-style-panel'),motifSelect=document.getElementById('motif-select'),motifLogo=document.getElementById('motif-logo'),rankRowsSel=document.getElementById('rank-rows'),rankRowsSlider=document.getElementById('rank-rows-slider'),panelCountSel=document.getElementById('panel-count'),aggregateLegendsToggle=document.getElementById('aggregate-legends'),reportDetail=document.getElementById('report-detail');
function escText(v){return String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}function b64ToBytes(b64){return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}async function decodePayload(){const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}function compPayload(idx){return review.comparisons[idx]?.payload||{points:[],aggregate:{motifs:[]},conditions:[],colors:{},groups:[]}}function compLabel(idx){return review.comparisons[idx]?.label||`Comparison ${idx+1}`}function motifLabel(item){if(!item)return'';const id=item.motif_id||item.id||'';return id?`${item.name} (${id})`:item.name}function sampleDisplayName(sample,condition){const name=String(sample?.name??sample??''),cond=String(condition??sample?.condition??'').trim();return cond&&name&&!name.startsWith(cond+'_')?`${cond}_${name}`:name}function allMotifs(){const map=new Map();review.comparisons.forEach(c=>[...((c.payload.aggregate||{}).motifs||[]),...(c.payload.points||[])].forEach(m=>{if(m&&m.prefix&&!map.has(m.prefix))map.set(m.prefix,m)}));return [...map.values()].sort((a,b)=>motifLabel(a).localeCompare(motifLabel(b),undefined,{sensitivity:'base'}))}function pointByPrefix(payload,prefix){return (payload.points||[]).find(p=>p.prefix===prefix)}function aggregateByPrefix(payload,prefix){return ((payload.aggregate||{}).motifs||[]).find(m=>m.prefix===prefix)}function allSampleRows(){const rows=[],seen=new Set();review.comparisons.forEach(c=>((c.payload.aggregate||{}).motifs||[]).forEach(m=>(m.conditions||[]).forEach(cond=>(cond.samples||[]).forEach(s=>{if(!seen.has(s.name)){seen.add(s.name);rows.push({name:s.name,condition:cond.name})}}))));return rows}function conditionColor(payload,condition){const colors=payload.colors||{};return colors[condition+'_up']||colors[condition]||'#64748b'}function groupColor(payload,group){return (payload.colors||{})[group]||'#8a94a6'}function sampleStyleKey(compIdx,name){return `${compIdx}::${name}`}function sampleStyle(key,defaults={}){const s=sampleLineStyles[key]||{};return{visible:s.visible!==false,color:s.color||defaults.color||'#64748b',alpha:s.alpha??.9,width:s.width||.7,type:s.type||'solid'}}function lineDash(t){return t==='dash'?'6 4':(t==='dot'?'1.2 3':'')}function borderStyle(t){return t==='dash'?'dashed':(t==='dot'?'dotted':'solid')}function dashAttr(t){const d=lineDash(t);return d?` stroke-dasharray="${d}"`:''}function lineWidth(v,f){const n=Number(v);return Number.isFinite(n)&&n>0?Math.max(.1,Math.min(8,n)):f}function alpha(v,f){const n=Number(v);return Number.isFinite(n)?Math.max(.05,Math.min(1,n)):f}function niceStep(raw){if(!Number.isFinite(raw)||raw<=0)return 1;const pow=Math.pow(10,Math.floor(Math.log10(raw))),f=raw/pow;return(f<=1?1:f<=1.5?1.5:f<=2.5?2.5:f<=5?5:10)*pow}function niceLimit(v){if(!Number.isFinite(v)||v<=0)return 1;const step=niceStep(Math.abs(v)/5);return Math.max(step,Math.ceil(Math.abs(v)/step)*step)}function niceTicks(min,max,n){const step=niceStep((max-min)/Math.max(1,n-1)),start=Math.ceil(min/step)*step,end=Math.floor(max/step)*step,out=[];for(let v=start;v<=end+step/2;v+=step)out.push(Number(v.toPrecision(12)));return out.length?out:[0]}function fmt(v){const a=Math.abs(v);if(!Number.isFinite(v))return'';if(a===0)return'0';if(a>=1)return v.toFixed(1);if(a>=.01)return v.toFixed(2);if(a>=.001)return v.toFixed(3);return v.toExponential(1)}
function boundedPanelCount(value){return Math.max(4,Math.min(8,Math.floor(Number(value)||4)))}function panelColumnCount(n){return n<=4?Math.min(2,n):Math.ceil(n/2)}function setPanelGridShape(){const n=Math.max(1,slotComparisons.length),cols=panelColumnCount(n);comparisonGrid.style.setProperty('--comparison-cols',cols);comparisonGrid.classList.toggle('compact-panels',cols>2);aggregateGrid.style.setProperty('--aggregate-cols',n);aggregateGrid.classList.toggle('hide-legends',!showAggregateLegends)}function setPanelCount(value){const available=review.comparisons.length,target=Math.min(available,boundedPanelCount(value)),previous=slotComparisons.slice(),next=[];for(let i=0;i<target;i++){let idx=previous[i];if(!Number.isInteger(idx)||idx<0||idx>=available)idx=i;next.push(idx)}slotComparisons=next;setPanelGridShape()}function initState(){panelCountSel.value=String(boundedPanelCount(initialDisplayPanels));aggregateLegendsToggle.checked=showAggregateLegends;setPanelCount(panelCountSel.value);const ranked=[...(review.comparisons[0]?.payload.points||[])].sort((a,b)=>Math.abs(b.change)-Math.abs(a.change));activePrefix=ranked[0]?.prefix||allMotifs()[0]?.prefix||null}
function comparisonSampleRows(compIdx){const rows=[],seen=new Set(),payload=compPayload(compIdx);((payload.aggregate||{}).motifs||[]).forEach(m=>(m.conditions||[]).forEach(cond=>(cond.samples||[]).forEach(s=>{const key=sampleStyleKey(compIdx,s.name);if(!seen.has(key)){seen.add(key);rows.push({key,name:s.name,label:sampleDisplayName(s,cond.name),condition:cond.name,color:conditionColor(payload,cond.name)})}})));return rows}function renderSampleStyles(){sampleStylePanel.innerHTML=slotComparisons.map((compIdx,slot)=>{const rows=comparisonSampleRows(compIdx);return `<div class="sample-style-block"><p class="section-title">Comparison ${slot+1}: ${escText(compLabel(compIdx))}</p>${rows.map(row=>{const st=sampleStyle(row.key,{color:row.color});return `<label class="sample-style-row"><input data-visible="${escText(row.key)}" type="checkbox" ${st.visible?'checked':''}><span class="sample-style-name" title="${escText(row.label)}">${escText(row.label)}</span><input data-color="${escText(row.key)}" type="color" value="${st.color}"><input data-alpha="${escText(row.key)}" type="number" min="0.05" max="1" step="0.05" value="${st.alpha}"><input data-width="${escText(row.key)}" type="number" min="0.2" max="5" step="0.1" value="${st.width}"><select data-type="${escText(row.key)}"><option value="solid"${st.type==='solid'?' selected':''}>Solid</option><option value="dash"${st.type==='dash'?' selected':''}>Dash</option><option value="dot"${st.type==='dot'?' selected':''}>Dot</option></select></label>`}).join('')}</div>`}).join('');sampleStylePanel.querySelectorAll('[data-visible]').forEach(el=>el.addEventListener('change',()=>{sampleLineStyles[el.dataset.visible]={...(sampleLineStyles[el.dataset.visible]||{}),visible:el.checked};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-color]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.color]={...(sampleLineStyles[el.dataset.color]||{}),color:el.value};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-alpha]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.alpha]={...(sampleLineStyles[el.dataset.alpha]||{}),alpha:alpha(el.value,.9)};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-width]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.width]={...(sampleLineStyles[el.dataset.width]||{}),width:lineWidth(el.value,.7)};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-type]').forEach(el=>el.addEventListener('change',()=>{sampleLineStyles[el.dataset.type]={...(sampleLineStyles[el.dataset.type]||{}),type:el.value};renderAll(false)}))}
function drawRank(payload,slot){const points=payload.points||[],limit=Math.max(2,Math.floor(Number(rankRowsSel.value||20))),perDir=Math.max(1,Math.floor(limit/2)),positive=points.filter(p=>p.change>0).sort((a,b)=>b.change-a.change||a.pvalue-b.pvalue).slice(0,perDir),negative=points.filter(p=>p.change<0).sort((a,b)=>a.change-b.change||a.pvalue-b.pvalue).slice(0,perDir),shown=[...negative,...positive],w=340,rowH=9,gap=2,margin={top:34,bottom:30,left:118,right:8},h=Math.max(260,margin.top+shown.length*(rowH+gap)+margin.bottom+4),mid=205,xW=96,maxAbs=niceLimit(Math.max(...shown.map(p=>Math.abs(p.change)),.01)),sx=v=>mid+(v/maxAbs)*xW,colors=payload.colors||{};let y=margin.top,parts=[`<svg class="rank-svg" data-slot="${slot}" viewBox="0 0 ${w} ${h}"><style>${plotSvgStyle}</style><rect width="${w}" height="${h}" fill="#fff"/><text x="${w/2}" y="13" class="plot-title" text-anchor="middle">Top differential motifs</text><text x="${mid-6}" y="28" text-anchor="end" font-size="11" font-weight="900" fill="${colors[payload.conditions?.[1]+'_up']||'#2563eb'}">${escText(payload.conditions?.[1]||'left')}_up</text><text x="${mid+6}" y="28" text-anchor="start" font-size="11" font-weight="900" fill="${colors[payload.conditions?.[0]+'_up']||'#dc2626'}">${escText(payload.conditions?.[0]||'right')}_up</text><line x1="${mid}" y1="22" x2="${mid}" y2="${h-margin.bottom}" stroke="#172033" stroke-width="2"/>`];function row(p){const bw=Math.abs(p.change)/maxAbs*xW,x=p.change>=0?mid:mid-bw,color=groupColor(payload,p.group),active=p.prefix===activePrefix;parts.push(`<text class="rank-name" data-prefix="${escText(p.prefix)}" x="5" y="${y+rowH-1}" font-size="8.5" font-weight="${active?'900':'700'}" fill="${active?color:'#526176'}">${escText(motifLabel(p)).slice(0,20)}</text><rect class="rank-bar${active?' selected':''}" data-prefix="${escText(p.prefix)}" x="${x}" y="${y}" width="${bw}" height="${rowH}" fill="${color}" fill-opacity="${active?.95:.72}"/><text x="${p.change>=0?x-3:x+bw+3}" y="${y+rowH-1}" class="tick" text-anchor="${p.change>=0?'end':'start'}">${fmt(p.change)}</text>`);y+=rowH+gap}negative.forEach(row);y+=4;positive.forEach(row);niceTicks(-maxAbs,maxAbs,5).forEach(t=>parts.push(`<text x="${sx(t)}" y="${h-19}" class="tick" text-anchor="middle">${fmt(t)}</text>`));parts.push(`<line x1="${sx(-maxAbs)}" y1="${h-27}" x2="${sx(maxAbs)}" y2="${h-27}" class="axis"/><text x="${mid}" y="${h-4}" class="axis-label" text-anchor="middle">${escText(payload.change_label||'Differential footprint score')}</text></svg>`);return parts.join('')}
function drawVolcano(payload,slot){const points=payload.points||[],w=430,h=270,margin={top:13,right:14,bottom:31,left:50},innerW=w-margin.left-margin.right,innerH=h-margin.top-margin.bottom,maxAbs=niceLimit(Math.max(1e-9,...points.map(p=>Math.abs(p.change)))*1.05),maxY=niceLimit(Math.max(1,...points.map(p=>p.neglog10p||0))*1.03),sx=v=>margin.left+((v+maxAbs)/(2*maxAbs))*innerW,sy=v=>margin.top+innerH-(v/maxY)*innerH,colors=payload.colors||{};let parts=[`<svg class="volcano-svg" data-slot="${slot}" viewBox="0 0 ${w} ${h}"><style>${plotSvgStyle}</style><rect width="${w}" height="${h}" fill="#fff"/>`];niceTicks(0,maxY,5).forEach(t=>parts.push(`<line x1="${margin.left}" y1="${sy(t)}" x2="${margin.left+innerW}" y2="${sy(t)}" class="grid"/><text x="${margin.left-6}" y="${sy(t)+3}" class="tick" text-anchor="end">${fmt(t)}</text>`));niceTicks(-maxAbs,maxAbs,5).forEach(t=>parts.push(`<line x1="${sx(t)}" y1="${margin.top}" x2="${sx(t)}" y2="${margin.top+innerH}" class="grid"/><text x="${sx(t)}" y="${margin.top+innerH+13}" class="tick" text-anchor="middle">${fmt(t)}</text>`));parts.push(`<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top+innerH}" class="zero"/><line x1="${margin.left}" y1="${margin.top+innerH}" x2="${margin.left+innerW}" y2="${margin.top+innerH}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+innerH}" class="axis"/>`);points.forEach(p=>{const active=p.prefix===activePrefix,color=groupColor(payload,p.group),r=active?5:1.8,opacity=p.group==='n.s.'?.45:.78;parts.push(`<circle class="pt${active?' selected':''}" data-prefix="${escText(p.prefix)}" cx="${sx(p.change)}" cy="${sy(p.neglog10p||0)}" r="${r}" fill="${color}" fill-opacity="${opacity}" stroke="${active?'#111827':'none'}"><title>${escText(motifLabel(p))}</title></circle>`)});parts.push(`<text x="${margin.left+4}" y="${margin.top+innerH-7}" font-size="12" font-weight="900" fill="${colors[payload.conditions?.[1]+'_up']||'#2563eb'}">${escText(payload.conditions?.[1]||'left')}_up</text><text x="${margin.left+innerW-4}" y="${margin.top+innerH-7}" text-anchor="end" font-size="12" font-weight="900" fill="${colors[payload.conditions?.[0]+'_up']||'#dc2626'}">${escText(payload.conditions?.[0]||'right')}_up</text><text x="${margin.left+innerW/2}" y="${h-4}" class="axis-label" text-anchor="middle">${escText(payload.change_label||'Differential footprint score')}</text><text x="13" y="${margin.top+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 13 ${margin.top+innerH/2})">-log10(p-value)</text></svg>`);return parts.join('')}
function renderComparisons(){comparisonGrid.innerHTML=slotComparisons.map((compIdx,slot)=>{const payload=compPayload(compIdx);const options=review.comparisons.map((c,i)=>`<option value="${i}" ${i===compIdx?'selected':''}>${escText(c.label)}</option>`).join('');return `<section class="comparison-card"><div class="comparison-head"><span class="section-title">Comparison ${slot+1}</span><select data-comparison-slot="${slot}">${options}</select></div><div class="pair"><div class="plot-box">${drawRank(payload,slot)}</div><div class="plot-box">${drawVolcano(payload,slot)}</div></div></section>`}).join('');comparisonGrid.querySelectorAll('[data-comparison-slot]').forEach(sel=>sel.addEventListener('change',()=>{slotComparisons[Number(sel.dataset.comparisonSlot)]=Number(sel.value);renderAll(true)}));comparisonGrid.querySelectorAll('[data-prefix]').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}
	function setSelectedMotif(prefix){if(!prefix)return;activePrefix=prefix;renderAll(false)}function motifLogoSvgFromCounts(counts,attrs=''){if(!Array.isArray(counts)||counts.length!==4||!Array.isArray(counts[0]))return'';const n=counts[0].length;if(!n)return'';const w=320,h=132,left=38,right=8,top=9,bottom=34,plotW=w-left-right,plotH=h-top-bottom,bases=['A','C','G','T'],colors={A:'#198754',C:'#0d6efd',G:'#f59f00',T:'#dc3545'},colW=plotW/Math.max(1,n);let parts=[`<svg ${attrs} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}"><rect width="100%" height="100%" fill="#fff"/><line x1="${left}" y1="${top+plotH}" x2="${left+plotW}" y2="${top+plotH}" stroke="#26364d" stroke-width="1.3"/><line x1="${left}" y1="${top}" x2="${left}" y2="${top+plotH}" stroke="#26364d" stroke-width="1.3"/><text x="13" y="${top+plotH/2}" transform="rotate(-90 13 ${top+plotH/2})" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="900" fill="#152133">bits</text><text x="${left+plotW/2}" y="${h-5}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="900" fill="#152133">position</text>`];[0,1,2].forEach(t=>{const y=top+plotH-(t/2)*plotH;parts.push(`<line x1="${left-4}" y1="${y.toFixed(2)}" x2="${left}" y2="${y.toFixed(2)}" stroke="#26364d" stroke-width="1"/>`,`<text x="${left-7}" y="${(y+3.5).toFixed(2)}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="900" fill="#334155">${t}</text>`)});for(let pos=0;pos<n;pos++){const col=counts.map(r=>Number(r[pos])||0),sum=col.reduce((a,b)=>a+b,0)||1,p=col.map(v=>v/sum),entropy=-p.reduce((a,v)=>a+(v>0?v*Math.log2(Math.max(v,1e-12)):0),0),bits=p.map(v=>v*Math.max(0,2-entropy)),order=[0,1,2,3].sort((a,b)=>bits[a]-bits[b]);let y=top+plotH,x=left+pos*colW+colW/2;if(n<=18||pos===0||pos===n-1||(pos+1)%5===0)parts.push(`<text x="${x.toFixed(2)}" y="${top+plotH+12}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8.5" font-weight="900" fill="#334155">${pos+1}</text>`);order.forEach(idx=>{const val=bits[idx];if(val<=.015)return;const lh=Math.max(2.5,val/2*plotH);y-=lh;const base=bases[idx],fs=Math.max(7,Math.min(30,lh*1.25));parts.push(`<text x="${x.toFixed(2)}" y="${(y+lh*.88).toFixed(2)}" text-anchor="middle" font-family="Arial Black,Arial,Helvetica,sans-serif" font-size="${fs.toFixed(2)}" font-weight="900" fill="${colors[base]}">${base}</text>`)})}parts.push('</svg>');return parts.join('')}function motifLogoHtml(prefix){for(const c of review.comparisons){const p=c.payload,counts=(p.motif_matrices||{})[prefix],logo=(p.logos||{})[prefix]||{};if(Array.isArray(counts)){const svg=motifLogoSvgFromCounts(counts);if(svg)return svg}if(logo.png)return `<img src="${logo.png}" alt="Motif logo">`;if(logo.svg)return logo.svg}return '<span class="tick">Motif logo unavailable</span>'}
function renderSelected(){const motifs=allMotifs();motifSelect.innerHTML=motifs.map(m=>`<option value="${escText(m.prefix)}" ${m.prefix===activePrefix?'selected':''}>${escText(motifLabel(m))}</option>`).join('');motifLogo.innerHTML=motifLogoHtml(activePrefix);motifSelect.onchange=()=>setSelectedMotif(motifSelect.value)}
function aggregateSamples(motif,compIdx,payload){const out=[];(motif.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>{const key=sampleStyleKey(compIdx,s.name),st=sampleStyle(key,{color:conditionColor(payload,c.name)});if(st.visible)out.push({...s,condition:c.name,style:st})}));return out}function pathD(profile,x,sx,sy){return profile.map((y,i)=>`${i?'L':'M'}${sx(x[i]).toFixed(2)},${sy(y).toFixed(2)}`).join(' ')}
function computeAggregateDomain(prefix){let vals=[];slotComparisons.forEach(compIdx=>{const payload=compPayload(compIdx),motif=aggregateByPrefix(payload,prefix);if(motif)aggregateSamples(motif,compIdx,payload).forEach(s=>vals=vals.concat((s.profile||[]).filter(Number.isFinite)))});let rawMin=Math.min(...vals,0),rawMax=Math.max(...vals,1e-9),pad=Math.max((rawMax-rawMin||1)*.06,1e-6),step=niceStep((rawMax-rawMin+2*pad)/4),ymin=Math.floor((rawMin-pad)/step)*step,ymax=Math.ceil((rawMax+pad)/step)*step;return[ymin,ymax]}
function drawAggregate(payload,prefix,compIdx){const motif=aggregateByPrefix(payload,prefix),w=180,h=180;if(!motif)return `<svg class="aggregate-panel" viewBox="0 0 ${w} ${h}"><style>${plotSvgStyle}</style><rect width="${w}" height="${h}" fill="#fff"/><text x="${w/2}" y="${h/2}" class="axis-label" text-anchor="middle">No aggregate profile</text></svg>`;const rawX=(payload.aggregate||{}).x||[],keep=rawX.map((v,i)=>({v,i})).filter(p=>p.v>=-aggregateDisplayBp&&p.v<=aggregateDisplayBp),x=keep.length?keep.map(p=>p.v):rawX,samples=aggregateSamples(motif,compIdx,payload).map(s=>({...s,profile:keep.length?keep.map(p=>s.profile[p.i]):s.profile})),domain=aggregateDomain||computeAggregateDomain(prefix);let ymin=domain[0],ymax=domain[1];const margin={top:4,right:10,bottom:18,left:30},innerW=w-margin.left-margin.right,innerH=h-margin.top-margin.bottom,sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH;let parts=[`<svg class="aggregate-panel" viewBox="0 0 ${w} ${h}"><style>${plotSvgStyle}</style><rect width="${w}" height="${h}" fill="#fff"/>`];niceTicks(ymin,ymax,4).forEach(t=>parts.push(`<text x="${margin.left-4}" y="${sy(t)+3}" class="tick" text-anchor="end">${fmt(t)}</text>`));[-aggregateDisplayBp,0,aggregateDisplayBp].forEach(t=>parts.push(`<line x1="${sx(t)}" y1="${margin.top}" x2="${sx(t)}" y2="${margin.top+innerH}" class="grid"/><text x="${sx(t)}" y="${margin.top+innerH+12}" class="tick" text-anchor="${t===-aggregateDisplayBp?'start':(t===aggregateDisplayBp?'end':'middle')}">${t}</text>`));parts.push(`<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top+innerH}" class="zero"/><line x1="${margin.left}" y1="${margin.top+innerH}" x2="${margin.left+innerW}" y2="${margin.top+innerH}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+innerH}" class="axis"/><text x="${margin.left+5}" y="${margin.top+innerH-5}" class="tick">${motif.n_sites||0}</text>`);samples.sort((a,b)=>(Number(b.fp_score||0)-Number(a.fp_score||0))).forEach(s=>{const st=s.style,label=sampleDisplayName(s,s.condition);parts.push(`<path d="${pathD(s.profile,x,sx,sy)}" fill="none" stroke="${st.color}" stroke-width="${lineWidth(st.width,.7)}"${dashAttr(st.type)} stroke-opacity="${alpha(st.alpha,.9)}"><title>${escText(label)}</title></path>`)});parts.push('</svg>');return parts.join('')}
function aggregateLegendHtml(payload,prefix,compIdx){const motif=aggregateByPrefix(payload,prefix);if(!motif)return'<div class="aggregate-legend-mini"></div>';const samples=aggregateSamples(motif,compIdx,payload);return `<div class="aggregate-legend-mini">${samples.map(s=>{const st=s.style,label=sampleDisplayName(s,s.condition);return `<div class="agg-legend-row"><i class="agg-legend-line" style="border-top-color:${st.color};border-top-width:${lineWidth(st.width,.7)}px;border-top-style:${borderStyle(st.type)};opacity:${alpha(st.alpha,.9)}"></i><span title="${escText(label)}">${escText(label)}</span></div>`}).join('')}</div>`}
function renderAggregates(){aggregateDomain=computeAggregateDomain(activePrefix);aggregateGrid.classList.toggle('hide-legends',!showAggregateLegends);aggregateGrid.innerHTML=slotComparisons.map((compIdx,slot)=>{const payload=compPayload(compIdx),label=`Comparison ${slot+1}`;return `<div class="aggregate-tile"><div class="aggregate-tile-label" title="${escText(label+': '+compLabel(compIdx))}">${escText(label)}</div><div class="aggregate-plot">${drawAggregate(payload,activePrefix,compIdx)}</div>${showAggregateLegends?aggregateLegendHtml(payload,activePrefix,compIdx):''}</div>`}).join('')}
	function downloadBlob(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}function downloadSvgList(selector,name){const svgs=[...document.querySelectorAll(selector)];if(!svgs.length)return;const w=320,h=320,cols=Math.min(4,svgs.length),rows=Math.ceil(svgs.length/cols);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${cols*w} ${rows*h}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`];svgs.forEach((s,i)=>parts.push(`<g transform="translate(${(i%cols)*w},${Math.floor(i/cols)*h}) scale(${w/Number(s.viewBox.baseVal.width||w)})">${s.innerHTML}</g>`));parts.push('</svg>');downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),name)}function aggregateTileMarkup(tile,x,y){const svg=tile.querySelector('.aggregate-panel');if(!svg)return{markup:'',width:0,height:0};const title=tile.querySelector('.aggregate-tile-label')?.textContent||'',titleH=title?14:0,plotW=180,plotH=180,rows=[...tile.querySelectorAll('.agg-legend-row')].map(row=>{const line=row.querySelector('.agg-legend-line'),span=row.querySelector('span');return{label:span?.textContent||'',color:line?.style.borderTopColor||'#000',width:lineWidth(parseFloat(line?.style.borderTopWidth),.7),type:line?.style.borderTopStyle||'solid',alpha:alpha(line?.style.opacity,.9)}}),legendW=rows.length?Math.max(134,...rows.map(r=>r.label.length*5.2+40)):0,width=plotW+(rows.length?12+legendW:6),height=Math.max(plotH+6+titleH,12+titleH+rows.length*13);let parts=[`<g transform="translate(${x},${y})"><rect width="${width}" height="${height}" rx="7" fill="#fff" stroke="#e1e8f0"/>`];if(title)parts.push(`<text x="5" y="10" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="900" fill="#172033">${escText(title)}</text>`);parts.push(`<g transform="translate(3,${3+titleH})">${svg.innerHTML}</g>`);rows.forEach((r,i)=>{const yy=titleH+15+i*13,dash=r.type==='dashed'?' stroke-dasharray="6 4"':(r.type==='dotted'?' stroke-dasharray="1.2 3"':'');parts.push(`<line x1="${plotW+14}" y1="${yy}" x2="${plotW+38}" y2="${yy}" stroke="${r.color}" stroke-width="${r.width}" stroke-opacity="${r.alpha}"${dash}/><text x="${plotW+44}" y="${yy+3}" font-family="Arial,Helvetica,sans-serif" font-size="8.5" font-weight="900" fill="#000">${escText(r.label)}</text>`) });parts.push('</g>');return{markup:parts.join(''),width,height}}function aggregateTilesData(){const tiles=[...document.querySelectorAll('.aggregate-tile')];if(!tiles.length)return{markup:'',width:0,height:0};const tileData=tiles.map(t=>aggregateTileMarkup(t,0,0)),cols=tileData.length,gap=8,colW=Math.max(...tileData.map(t=>t.width),1),rowH=Math.max(...tileData.map(t=>t.height),1),width=cols*colW+(cols-1)*gap,height=rowH;let parts=[];tiles.forEach((tile,i)=>parts.push(aggregateTileMarkup(tile,i*(colW+gap),0).markup));return{markup:parts.join(''),width,height}}function aggregateTilesSvg(){const data=aggregateTilesData();if(!data.markup)return'';return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${data.width} ${data.height}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>${data.markup}</svg>`}function downloadAggregateTiles(){const svg=aggregateTilesSvg();if(svg)downloadBlob(new Blob([svg],{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_aggregate.svg')}function downloadCombinedPanel(){const svgs=[...document.querySelectorAll('.rank-svg,.volcano-svg')],w=320,h=320,cols=Math.min(4,Math.max(1,svgs.length)),rows=Math.ceil(svgs.length/cols),agg=aggregateTilesData(),topW=cols*w,topH=rows*h,gap=18,width=Math.max(topW,agg.width),height=topH+(agg.markup?gap+agg.height:0);let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`];svgs.forEach((s,i)=>parts.push(`<g transform="translate(${(i%cols)*w},${Math.floor(i/cols)*h}) scale(${w/Number(s.viewBox.baseVal.width||w)})">${s.innerHTML}</g>`));if(agg.markup)parts.push(`<g transform="translate(0,${topH+gap})">${agg.markup}</g>`);parts.push('</svg>');downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_panel.svg')}function logoMarkup(prefix,x,y,w,h){for(const c of review.comparisons){const p=c.payload,counts=(p.motif_matrices||{})[prefix],logo=(p.logos||{})[prefix]||{};if(Array.isArray(counts)){const svg=motifLogoSvgFromCounts(counts,`x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet"`);if(svg)return svg}if(logo.png)return `<image x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet" href="${logo.png}"/>`;if(logo.svg)return logo.svg.replace(/<\\?xml[^>]*>/g,'').replace(/<!DOCTYPE[^>]*>/g,'').replace(/<svg\\b/i,`<svg x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet"`)}return `<text x="${x+w/2}" y="${y+h/2}" class="tick" text-anchor="middle">Motif logo unavailable</text>`}function downloadLogoPanel(){const motif=allMotifs().find(m=>m.prefix===activePrefix)||{prefix:activePrefix,name:activePrefix};let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 310 190" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/><rect x="1" y="1" width="308" height="188" rx="7" fill="#fff" stroke="#d8e2ef"/><text x="10" y="20" class="plot-title">${escText(motifLabel(motif)).slice(0,38)}</text>${logoMarkup(activePrefix,16,34,278,144)}</svg>`];downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_motif_logo_panel.svg')}
function syncRankRows(source){const max=Math.max(2,Number(rankRowsSel.max)||200),value=Math.max(2,Math.min(max,Math.floor(Number(source.value)||20)));rankRowsSel.value=value;rankRowsSlider.value=value;renderAll(false)}
function renderAll(refreshStyles=true){setPanelGridShape();if(refreshStyles)renderSampleStyles();renderComparisons();renderSelected();renderAggregates();reportDetail.textContent=`${review.comparisons.length} comparisons - ${allMotifs().length} motifs - ${allSampleRows().length} samples - ${slotComparisons.length} displayed panels - 1 selected motif`}document.getElementById('download-rank').addEventListener('click',()=>downloadSvgList('.rank-svg','review_multi_comparisons_barplots.svg'));document.getElementById('download-volcano').addEventListener('click',()=>downloadSvgList('.volcano-svg','review_multi_comparisons_volcano.svg'));document.getElementById('download-aggregate').addEventListener('click',downloadAggregateTiles);document.getElementById('download-panel').addEventListener('click',downloadCombinedPanel);document.getElementById('download-logo').addEventListener('click',downloadLogoPanel);rankRowsSel.addEventListener('input',()=>syncRankRows(rankRowsSel));rankRowsSlider.addEventListener('input',()=>syncRankRows(rankRowsSlider));panelCountSel.addEventListener('change',()=>{setPanelCount(panelCountSel.value);renderAll(true)});aggregateLegendsToggle.addEventListener('change',()=>{showAggregateLegends=aggregateLegendsToggle.checked;renderAll(false)});decodePayload().then(data=>{review=data;initState();const maxRows=Math.max(20,...review.comparisons.map(c=>(c.payload.points||[]).length));rankRowsSel.max=maxRows;rankRowsSlider.max=maxRows;renderAll()}).catch(err=>{reportDetail.textContent=`Could not open review payload: ${err.message}`});
</script></body></html>'''
    document = template.replace("__TITLE__", escaped_title).replace("__PAYLOAD__", payload_b64).replace("__DISPLAY_PANELS__", str(display_panels)).replace("__AGGREGATE_LEGENDS__", aggregate_legends)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="review-multi-comparisons", description="Review multiple diff-footprints HTML reports in one interactive HTML file.")
    parser.add_argument("--inputs", nargs="+", help="diff-footprints HTML files or directories containing diff_footprints_*.html files; directories are searched recursively.")
    parser.add_argument("--labels", nargs="*", help="Optional labels, one per resolved input HTML.")
    parser.add_argument("--output", help="Output standalone review HTML.")
    parser.add_argument("--outdir", help="Project directory used with --layout project.")
    parser.add_argument("--layout", choices=["custom", "project"], default="project", help="Use fp-tools standard project output layout under --outdir (default: project when only --outdir is provided).")
    parser.add_argument("--display-panels", type=int, default=4, help="Initial number of comparison panels to display in the HTML report, from 4 to 8 (default: 4).")
    parser.add_argument("--aggregate-legends", choices=["show", "hide"], default="show", help="Initial visibility for legends beside motif aggregate subplots (default: show). Use hide to fit 4-8 aggregate panels in one row.")
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
            if not args.output:
                args.output = str(review_output_path(project))
    if not args.inputs:
        parser.error("provide --inputs or use --layout project with --outdir")
    if not args.output:
        parser.error("provide --output or use --layout project with --outdir")
    if args.display_panels < 4 or args.display_panels > 8:
        parser.error("--display-panels must be between 4 and 8")
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
            before_missing, total = count_missing_aggregate_profiles(payload)
            fill_stats = {"before_missing": before_missing, "after_missing": before_missing, "filled": 0, "total": total}
        write_review_html(payload, args.output, display_panels=args.display_panels, aggregate_legends=args.aggregate_legends)
    except ValueError as exc:
        parser.error(str(exc))
    print(f"Wrote {args.output}")
    if fill_stats["before_missing"]:
        print(f"Aggregate profiles: filled {fill_stats['filled']} missing panels; {fill_stats['after_missing']} remain missing of {fill_stats['total']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
