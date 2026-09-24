"use strict";

// Classic presentation reconstructed from the verified pre-control report.
// Scientific payloads are supplied by the shared loader and are never edited.
function createClassicView(host) {
  const byId = id => document.getElementById('classic-' + id);
  host.innerHTML = "<div class=\"wrap\"><p class=\"sr-only\" id=\"classic-report-detail\">Loading report</p><main class=\"board\"><aside class=\"side\"><select id=\"classic-layout\" class=\"classic-view-choice\" aria-label=\"Report view\"><option value=\"single\">Single comparison</option><option value=\"side\">Side by side</option><option value=\"classic\" selected>Classic</option></select><section class=\"card\"><p class=\"section-title\">Sample line styles</p><div id=\"classic-sample-style-panel\" class=\"sample-style-panel\"></div></section><section class=\"card\"><p class=\"section-title\">Export editable SVG</p><div class=\"export-stack\"><button id=\"classic-download-logo\">Download motif logo panel</button><button id=\"classic-download-rank\">Download bar plot panel</button><button id=\"classic-download-volcano\">Download volcano plot panel</button><button id=\"classic-download-aggregate\">Download motif aggregate panel</button><button id=\"classic-download-panel\">Download combined panel</button></div><label class=\"rows-control\">Top motifs <input id=\"classic-rank-rows-slider\" type=\"range\" min=\"1\" max=\"200\" step=\"1\" value=\"20\"><input id=\"classic-rank-rows\" type=\"number\" min=\"1\" max=\"200\" step=\"1\" value=\"20\"></label><label class=\"rows-control\">Panels <span></span><select id=\"classic-panel-count\" aria-label=\"Number of comparison panels\"><option value=\"4\">4</option><option value=\"5\">5</option><option value=\"6\">6</option><option value=\"7\">7</option><option value=\"8\">8</option></select></label></section><section class=\"card\"><p class=\"section-title\">Selected motif</p><select id=\"classic-motif-select\" class=\"motif-select\"></select><div id=\"classic-motif-logo\" class=\"motif-logo\"></div></section><section id=\"classic-controls\" class=\"card classic-controls\"><p class=\"section-title\">Plot controls</p><div id=\"classic-filter\"></div><details open><summary>Plot ranges</summary><label>Apply ranges to <select id=\"classic-edit-panel\" aria-label=\"Apply ranges to\"></select></label><p class=\"classic-range-help\">All includes comparisons opened later. Editing All replaces individual limits for that axis.</p><div id=\"classic-ranges\"></div></details><p id=\"classic-status\" class=\"classic-status\" role=\"status\"></p></section></aside><section id=\"classic-comparison-grid\" class=\"plots comparison-grid\"></section><section class=\"aggregate-card\"><div class=\"aggregate-head\"><p class=\"section-title\">Motif aggregate review</p><span class=\"sub\">Group autoscale</span></div><div id=\"classic-aggregate-grid\" class=\"aggregate-grid\"></div></section></main></div>";
  const aggregateDisplayBp=60, initialDisplayPanels=8;
  const plotSvgStyle="svg,text{font-family:Arial,Helvetica,sans-serif}.axis{stroke:#000;stroke-width:1.2}.zero{stroke:#555;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e7edf5;stroke-width:1}.tick{font-size:10px;fill:#000;font-weight:900}.axis-label{font-size:11px;fill:#000;font-weight:900}.plot-title{font-size:13px;font-weight:900;fill:#000}";
  let review=null,slotComparisons=[],activePrefix=null,sampleLineStyles={},aggregateDomain=null;
  let generation=0, editing='all', clipSerial=0;
  const profiles=new Map(), failures=new Map(), ranges=new Map(), profileLoads=new Map();
  const comparisonGrid=byId('comparison-grid'),aggregateGrid=byId('aggregate-grid'),
    sampleStylePanel=byId('sample-style-panel'),motifSelect=byId('motif-select'),motifLogo=byId('motif-logo'),
    rankRowsSel=byId('rank-rows'),rankRowsSlider=byId('rank-rows-slider'),panelCountSel=byId('panel-count'),
    reportDetail=byId('report-detail');
function escText(v){return String(v??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}

function compPayload(idx){return review.comparisons[idx]?.payload||{points:[],aggregate:{motifs:[]},conditions:[],colors:{},groups:[]}}

function compLabel(idx){return review.comparisons[idx]?.label||`Comparison ${idx+1}`}

function motifLabel(item){if(!item)return'';const id=item.motif_id||item.id||'';return id?`${item.name} (${id})`:item.name}

function sampleDisplayName(sample,condition){const name=String(sample?.name??sample??''),cond=String(condition??sample?.condition??'').trim();return cond&&name&&!name.startsWith(cond+'_')?`${cond}_${name}`:name}

function allMotifs(){const map=new Map();review.comparisons.forEach(c=>[...((c.payload.aggregate||{}).motifs||[]),...(c.payload.points||[])].forEach(m=>{if(m&&m.prefix&&!map.has(m.prefix))map.set(m.prefix,m)}));return [...map.values()].sort((a,b)=>motifLabel(a).localeCompare(motifLabel(b),undefined,{sensitivity:'base'}))}

function pointByPrefix(payload,prefix){return (payload.points||[]).find(p=>p.prefix===prefix)}

function allSampleRows(){const rows=[],seen=new Set();review.comparisons.forEach(c=>((c.payload.aggregate||{}).motifs||[]).forEach(m=>(m.conditions||[]).forEach(cond=>(cond.samples||[]).forEach(s=>{if(!seen.has(s.name)){seen.add(s.name);rows.push({name:s.name,condition:cond.name})}}))));return rows}

function conditionColor(payload,condition){const colors=payload.colors||{};return colors[condition+'_up']||colors[condition]||'#64748b'}

function sampleStyleKey(compIdx,name){return `${compIdx}::${name}`}

function sampleStyle(key,defaults={}){const s=sampleLineStyles[key]||{};return{visible:s.visible!==false,color:s.color||defaults.color||'#64748b',alpha:s.alpha??.9,width:s.width||.7,type:s.type||'solid'}}

function lineDash(t){return t==='dash'?'6 4':(t==='dot'?'1.2 3':'')}

function borderStyle(t){return t==='dash'?'dashed':(t==='dot'?'dotted':'solid')}

function dashAttr(t){const d=lineDash(t);return d?` stroke-dasharray="${d}"`:''}

function lineWidth(v,f){const n=Number(v);return Number.isFinite(n)&&n>0?Math.max(.1,Math.min(8,n)):f}

function alpha(v,f){const n=Number(v);return Number.isFinite(n)?Math.max(.05,Math.min(1,n)):f}

function niceStep(raw){if(!Number.isFinite(raw)||raw<=0)return 1;const pow=Math.pow(10,Math.floor(Math.log10(raw))),f=raw/pow;return(f<=1?1:f<=1.5?1.5:f<=2.5?2.5:f<=5?5:10)*pow}

function niceLimit(v){if(!Number.isFinite(v)||v<=0)return 1;const step=niceStep(Math.abs(v)/5);return Math.max(step,Math.ceil(Math.abs(v)/step)*step)}

function niceTicks(min,max,n){const step=niceStep((max-min)/Math.max(1,n-1)),start=Math.ceil(min/step)*step,end=Math.floor(max/step)*step,out=[];for(let v=start;v<=end+step/2;v+=step)out.push(Number(v.toPrecision(12)));return out.length?out:[0]}

function fmt(v){const a=Math.abs(v);if(!Number.isFinite(v))return'';if(a===0)return'0';if(a>=1)return v.toFixed(1);if(a>=.01)return v.toFixed(2);if(a>=.001)return v.toFixed(3);return v.toExponential(1)}

function panelColumnCount(n){return n<=4?Math.min(2,n):Math.ceil(n/2)}

function setPanelGridShape(){const n=Math.max(1,slotComparisons.length),cols=panelColumnCount(n);comparisonGrid.style.setProperty('--comparison-cols',cols);comparisonGrid.classList.toggle('compact-panels',cols>2);aggregateGrid.style.setProperty('--aggregate-cols',cols)}

function setPanelCount(value){const available=review.comparisons.length,target=Math.min(available,boundedPanelCount(value)),previous=slotComparisons.slice(),next=[];for(let i=0;i<target;i++){let idx=previous[i];if(!Number.isInteger(idx)||idx<0||idx>=available)idx=i;next.push(idx)}slotComparisons=next;setPanelGridShape()}

function initState(){panelCountSel.value=String(boundedPanelCount(initialDisplayPanels));setPanelCount(panelCountSel.value);const ranked=[...(review.comparisons[0]?.payload.points||[])].sort((a,b)=>Math.abs(b.change)-Math.abs(a.change));activePrefix=ranked[0]?.prefix||allMotifs()[0]?.prefix||null}

function comparisonSampleRows(compIdx){const rows=[],seen=new Set(),payload=compPayload(compIdx);((payload.aggregate||{}).motifs||[]).forEach(m=>(m.conditions||[]).forEach(cond=>(cond.samples||[]).forEach(s=>{const key=sampleStyleKey(compIdx,s.name);if(!seen.has(key)){seen.add(key);rows.push({key,name:s.name,label:sampleDisplayName(s,cond.name),condition:cond.name,color:conditionColor(payload,cond.name)})}})));return rows}

function renderSampleStyles(){sampleStylePanel.innerHTML=slotComparisons.map((compIdx,slot)=>{const rows=comparisonSampleRows(compIdx);return `<div class="sample-style-block"><p class="section-title">Comparison ${slot+1}: ${escText(compLabel(compIdx))}</p>${rows.map(row=>{const st=sampleStyle(row.key,{color:row.color});return `<label class="sample-style-row"><input data-visible="${escText(row.key)}" type="checkbox" ${st.visible?'checked':''}><span class="sample-style-name" title="${escText(row.label)}">${escText(row.label)}</span><input data-color="${escText(row.key)}" type="color" value="${st.color}"><input data-alpha="${escText(row.key)}" type="number" min="0.05" max="1" step="0.05" value="${st.alpha}"><input data-width="${escText(row.key)}" type="number" min="0.2" max="5" step="0.1" value="${st.width}"><select data-type="${escText(row.key)}"><option value="solid"${st.type==='solid'?' selected':''}>Solid</option><option value="dash"${st.type==='dash'?' selected':''}>Dash</option><option value="dot"${st.type==='dot'?' selected':''}>Dot</option></select></label>`}).join('')}</div>`}).join('');sampleStylePanel.querySelectorAll('[data-visible]').forEach(el=>el.addEventListener('change',()=>{sampleLineStyles[el.dataset.visible]={...(sampleLineStyles[el.dataset.visible]||{}),visible:el.checked};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-color]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.color]={...(sampleLineStyles[el.dataset.color]||{}),color:el.value};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-alpha]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.alpha]={...(sampleLineStyles[el.dataset.alpha]||{}),alpha:alpha(el.value,.9)};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-width]').forEach(el=>el.addEventListener('input',()=>{sampleLineStyles[el.dataset.width]={...(sampleLineStyles[el.dataset.width]||{}),width:lineWidth(el.value,.7)};renderAll(false)}));sampleStylePanel.querySelectorAll('[data-type]').forEach(el=>el.addEventListener('change',()=>{sampleLineStyles[el.dataset.type]={...(sampleLineStyles[el.dataset.type]||{}),type:el.value};renderAll(false)}))}

function renderComparisons(){comparisonGrid.innerHTML=slotComparisons.map((compIdx,slot)=>{const payload=compPayload(compIdx);const options=review.comparisons.map((c,i)=>`<option value="${i}" ${i===compIdx?'selected':''}>${escText(c.label)}</option>`).join('');return `<section class="comparison-card"><div class="comparison-head"><span class="section-title">Comparison ${slot+1}</span><select data-comparison-slot="${slot}">${options}</select></div><div class="pair"><div class="plot-box">${drawRank(payload,slot)}</div><div class="plot-box">${drawVolcano(payload,slot)}</div></div></section>`}).join('');comparisonGrid.querySelectorAll('[data-comparison-slot]').forEach(sel=>sel.addEventListener('change',()=>{slotComparisons[Number(sel.dataset.comparisonSlot)]=Number(sel.value);renderAll(true)}));comparisonGrid.querySelectorAll('[data-prefix]').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}

function setSelectedMotif(prefix){if(!prefix)return;activePrefix=prefix;renderAll(false)}

function motifLogoSvgFromCounts(counts,attrs=''){if(!validMatrix(counts))return '';if(!Array.isArray(counts)||counts.length!==4||!Array.isArray(counts[0]))return'';const n=counts[0].length;if(!n)return'';const w=320,h=132,left=38,right=8,top=9,bottom=34,plotW=w-left-right,plotH=h-top-bottom,bases=['A','C','G','T'],colors={A:'#198754',C:'#0d6efd',G:'#f59f00',T:'#dc3545'},colW=plotW/Math.max(1,n);let parts=[`<svg ${attrs} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}"><rect width="100%" height="100%" fill="#fff"/><line x1="${left}" y1="${top+plotH}" x2="${left+plotW}" y2="${top+plotH}" stroke="#26364d" stroke-width="1.3"/><line x1="${left}" y1="${top}" x2="${left}" y2="${top+plotH}" stroke="#26364d" stroke-width="1.3"/><text x="13" y="${top+plotH/2}" transform="rotate(-90 13 ${top+plotH/2})" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="900" fill="#152133">bits</text><text x="${left+plotW/2}" y="${h-5}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="900" fill="#152133">position</text>`];[0,1,2].forEach(t=>{const y=top+plotH-(t/2)*plotH;parts.push(`<line x1="${left-4}" y1="${y.toFixed(2)}" x2="${left}" y2="${y.toFixed(2)}" stroke="#26364d" stroke-width="1"/>`,`<text x="${left-7}" y="${(y+3.5).toFixed(2)}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="900" fill="#334155">${t}</text>`)});for(let pos=0;pos<n;pos++){const col=counts.map(r=>Number(r[pos])||0),sum=col.reduce((a,b)=>a+b,0)||1,p=col.map(v=>v/sum),entropy=-p.reduce((a,v)=>a+(v>0?v*Math.log2(Math.max(v,1e-12)):0),0),bits=p.map(v=>v*Math.max(0,2-entropy)),order=[0,1,2,3].sort((a,b)=>bits[a]-bits[b]);let y=top+plotH,x=left+pos*colW+colW/2;if(n<=18||pos===0||pos===n-1||(pos+1)%5===0)parts.push(`<text x="${x.toFixed(2)}" y="${top+plotH+12}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="8.5" font-weight="900" fill="#334155">${pos+1}</text>`);order.forEach(idx=>{const val=bits[idx];if(val<=.015)return;const lh=Math.max(2.5,val/2*plotH);y-=lh;const base=bases[idx],fs=Math.max(7,Math.min(30,lh*1.25));parts.push(`<text x="${x.toFixed(2)}" y="${(y+lh*.88).toFixed(2)}" text-anchor="middle" font-family="Arial Black,Arial,Helvetica,sans-serif" font-size="${fs.toFixed(2)}" font-weight="900" fill="${colors[base]}">${base}</text>`)})}parts.push('</svg>');return parts.join('')}

function motifLogoHtml(prefix){for(const c of review.comparisons){const p=c.payload,counts=(p.motif_matrices||{})[prefix],logo=(p.logos||{})[prefix]||{};if(validMatrix(counts)){const svg=motifLogoSvgFromCounts(counts);if(svg)return svg}if(logo.png)return `<img src="${logo.png}" alt="Motif logo">`;if(logo.svg)return logo.svg}return '<span class="tick">Motif logo unavailable</span>'}

function renderSelected(){const motifs=allMotifs();motifSelect.innerHTML=motifs.map(m=>`<option value="${escText(m.prefix)}" ${m.prefix===activePrefix?'selected':''}>${escText(motifLabel(m))}</option>`).join('');motifLogo.innerHTML=motifLogoHtml(activePrefix);motifSelect.onchange=()=>setSelectedMotif(motifSelect.value)}

function aggregateSamples(motif,compIdx,payload){const out=[];(motif.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>{const key=sampleStyleKey(compIdx,s.name),st=sampleStyle(key,{color:conditionColor(payload,c.name)});if(st.visible)out.push({...s,condition:c.name,style:st})}));return out}

function aggregateLegendHtml(payload,prefix,compIdx){const motif=aggregateByPrefix(payload,prefix);if(!motif)return'<div class="aggregate-legend-mini"></div>';const samples=aggregateSamples(motif,compIdx,payload);return `<div class="aggregate-legend-mini">${samples.map(s=>{const st=s.style,label=sampleDisplayName(s,s.condition);return `<div class="agg-legend-row"><i class="agg-legend-line" style="border-top-color:${st.color};border-top-width:${lineWidth(st.width,.7)}px;border-top-style:${borderStyle(st.type)};opacity:${alpha(st.alpha,.9)}"></i><span title="${escText(label)}">${escText(label)}</span></div>`}).join('')}</div>`}

function downloadBlob(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}

function aggregateTileMarkup(tile,x,y){const svg=tile.querySelector('.aggregate-panel');if(!svg)return{markup:'',width:0,height:0};const title=tile.querySelector('.aggregate-tile-label')?.textContent||'',titleH=title?14:0,plotW=180,plotH=180,rows=[...tile.querySelectorAll('.agg-legend-row')].map(row=>{const line=row.querySelector('.agg-legend-line'),span=row.querySelector('span');return{label:span?.textContent||'',color:line?.style.borderTopColor||'#000',width:lineWidth(parseFloat(line?.style.borderTopWidth),.7),type:line?.style.borderTopStyle||'solid',alpha:alpha(line?.style.opacity,.9)}}),legendW=Math.max(134,...rows.map(r=>r.label.length*5.2+40)),width=plotW+12+legendW,height=Math.max(plotH+6+titleH,12+titleH+rows.length*13);let parts=[`<g transform="translate(${x},${y})"><rect width="${width}" height="${height}" rx="7" fill="#fff" stroke="#e1e8f0"/>`];if(title)parts.push(`<text x="5" y="10" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="900" fill="#172033">${escText(title)}</text>`);parts.push(`<g transform="translate(3,${3+titleH})">${svg.innerHTML}</g>`);rows.forEach((r,i)=>{const yy=titleH+15+i*13,dash=r.type==='dashed'?' stroke-dasharray="6 4"':(r.type==='dotted'?' stroke-dasharray="1.2 3"':'');parts.push(`<line x1="${plotW+14}" y1="${yy}" x2="${plotW+38}" y2="${yy}" stroke="${r.color}" stroke-width="${r.width}" stroke-opacity="${r.alpha}"${dash}/><text x="${plotW+44}" y="${yy+3}" font-family="Arial,Helvetica,sans-serif" font-size="8.5" font-weight="900" fill="#000">${escText(r.label)}</text>`) });parts.push('</g>');return{markup:parts.join(''),width,height}}

function aggregateTilesData(){const tiles=[...host.querySelectorAll('.aggregate-tile')];if(!tiles.length)return{markup:'',width:0,height:0};const tileData=tiles.map(t=>aggregateTileMarkup(t,0,0)),cols=Math.min(4,tileData.length),gap=8,colW=Math.max(...tileData.map(t=>t.width),1),rowH=Math.max(...tileData.map(t=>t.height),1),rows=Math.ceil(tileData.length/cols),width=cols*colW+(cols-1)*gap,height=rows*rowH+(rows-1)*gap;let parts=[];tiles.forEach((tile,i)=>parts.push(aggregateTileMarkup(tile,(i%cols)*(colW+gap),Math.floor(i/cols)*(rowH+gap)).markup));return{markup:parts.join(''),width,height}}

function aggregateTilesSvg(){const data=aggregateTilesData();if(!data.markup)return'';return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${data.width} ${data.height}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>${data.markup}</svg>`}

function downloadAggregateTiles(){const svg=aggregateTilesSvg();if(svg)downloadBlob(new Blob([svg],{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_aggregate.svg')}

function logoMarkup(prefix,x,y,w,h){for(const c of review.comparisons){const p=c.payload,counts=(p.motif_matrices||{})[prefix],logo=(p.logos||{})[prefix]||{};if(validMatrix(counts)){const svg=motifLogoSvgFromCounts(counts,`x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet"`);if(svg)return svg}if(logo.png)return `<image x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet" href="${logo.png}"/>`;if(logo.svg)return logo.svg.replace(/<\?xml[^>]*>/g,'').replace(/<!DOCTYPE[^>]*>/g,'').replace(/<svg\b/i,`<svg x="${x}" y="${y}" width="${w}" height="${h}" preserveAspectRatio="xMidYMid meet"`)}return `<text x="${x+w/2}" y="${y+h/2}" class="tick" text-anchor="middle">Motif logo unavailable</text>`}

function downloadLogoPanel(){const motif=allMotifs().find(m=>m.prefix===activePrefix)||{prefix:activePrefix,name:activePrefix};let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 310 190" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/><rect x="1" y="1" width="308" height="188" rx="7" fill="#fff" stroke="#d8e2ef"/><text x="10" y="20" class="plot-title">${escText(motifLabel(motif)).slice(0,38)}</text>${logoMarkup(activePrefix,16,34,278,144)}</svg>`];downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_motif_logo_panel.svg')}

  function validMatrix(counts) {
    return Array.isArray(counts) && counts.length===4 && Array.isArray(counts[0]) && counts[0].length>0 &&
      counts.every(row=>Array.isArray(row)&&row.length===counts[0].length&&row.every(v=>plotControls.numeric(v)&&Number(v)>=0));
  }
  function boundedPanelCount(value){return Math.max(1,Math.min(8,Math.floor(Number(value)||8)));}
  function aggregateByPrefix(payload,prefix){
    const index=review.comparisons.findIndex(c=>c.payload===payload);
    return profiles.get(index+':'+prefix)||((payload.aggregate||{}).motifs||[]).find(m=>m.prefix===prefix);
  }
  function pass(point){return plotControls.passesDisplayFilter(point,view.filter);}
  function groupColor(payload,group){return (payload.colors||{})[group]||'#8a94a6';}
  function pointColor(payload,p){
    return pass(p)?conditionColor(payload,payload.conditions[Number(p.change)<0?1:0]):'#8a94a6';
  }
  function caption(){
    return (view.filter.metric==='fdr'?'FDR':'p-value')+' ≤ '+view.filter.alpha+'; |ΔFP| ≥ '+view.filter.delta;
  }
  function ranked(payload){
    const r=plotControls.rankMotifs((payload.points||[]).filter(pass),'effect',Number(rankRowsSel.value));
    return [...r.negative,...r.positive];
  }
  function observations(payload,kind){
    return (kind==='rankX'?ranked(payload):(payload.points||[])).map(p=>kind==='volcanoY'?plotControls.negLog10P(p):p.change).filter(plotControls.numeric).map(Number);
  }
  function autoRange(kind,target='all'){
    if(kind==='aggregateY')return computeAggregateDomain(activePrefix,target);
    const payloads=target==='all'?review.comparisons.map(c=>c.payload):[compPayload(Number(target))];
    const values=payloads.flatMap(p=>observations(p,kind));
    if(kind==='volcanoY')return [0,niceLimit(values.reduce((m,v)=>Math.max(m,v),1)*1.03)];
    const maximum=values.reduce((m,v)=>Math.max(m,Math.abs(v)),.01);
    const bound=niceLimit(maximum*(kind==='volcanoX'?1.05:1));
    return [-bound,bound];
  }
  function rangeKey(kind,target){
    return kind+':'+target;
  }
  function targetRange(kind,target){
    const key=rangeKey(kind,target);
    // An explicit null is an individual autoscale, not an inherited All range.
    if(ranges.has(key))return ranges.get(key)||autoRange(kind,target);
    return ranges.get(rangeKey(kind,'all'))||autoRange(kind);
  }
  function setRange(kind,target,next){
    if(target==='all'){
      for(const key of ranges.keys())if(key.startsWith(kind+':'))ranges.delete(key);
    }
    ranges.set(rangeKey(kind,target),next);
  }
  function rangeFor(kind,slot){
    return targetRange(kind,String(slotComparisons[slot]));
  }
  function computeAggregateDomain(prefix,target='all'){
    const values=[];
    review.comparisons.forEach((comparison,index)=>{
      if(target!=='all'&&index!==Number(target))return;
      const payload=compPayload(index),motif=aggregateByPrefix(payload,prefix);
      if(motif)aggregateSamples(motif,index,payload).forEach(s=>(s.profile||[]).forEach((v,i)=>{
        const x=payload.aggregate?.x?.[i];
        if(x>=-60&&x<=60&&plotControls.numeric(v))values.push(Number(v));
      }));
    });
    let low=0,high=1e-9;
    values.forEach(v=>{low=Math.min(low,v);high=Math.max(high,v);});
    const pad=Math.max((high-low||1)*.06,1e-6),step=niceStep((high-low+2*pad)/4);
    return [Math.floor((low-pad)/step)*step,Math.ceil((high+pad)/step)*step];
  }
  function frame(w,h,kind,slot,xRange,yRange){
    return '<svg class="'+kind+'" data-slot="'+slot+'"'+
      (xRange?' data-x-range="'+escText(JSON.stringify(xRange))+'"':'')+
      (yRange?' data-y-range="'+escText(JSON.stringify(yRange))+'"':'')+
      ' viewBox="0 0 '+w+' '+h+'"><rect width="'+w+'" height="'+h+'" fill="#fff"/>';
  }
  function clipping(x,y,width,height){
    const id='classic-clip-'+(++clipSerial);
    return {id,markup:'<defs><clipPath id="'+id+'"><rect x="'+x+'" y="'+y+'" width="'+width+'" height="'+height+'"/></clipPath></defs>'};
  }
  function drawRank(payload,slot){
    const shown=ranked(payload),[low,high]=rangeFor('rankX',slot),w=340,rowH=9,gap=2,left=118,right=330;
    const h=Math.max(260,48+shown.length*(rowH+gap)+38),bottom=h-34,sx=v=>left+(v-low)/(high-low)*(right-left);
    const clip=clipping(left,32,right-left,bottom-32),parts=[frame(w,h,'rank-svg',slot,[low,high]),clip.markup];
    parts.push('<text x="170" y="13" class="plot-title" text-anchor="middle">Top differential motifs</text>');
    const clipped=shown.some(p=>p.change<low||p.change>high);
    parts.push('<text x="170" y="25" font-size="7" text-anchor="middle">'+escText(caption()+(clipped?' · Range clips data':''))+'</text>');
    const mid=sx(0);
    if(low<=0&&high>=0)parts.push('<text x="'+(mid-6)+'" y="37" text-anchor="end" font-size="11" fill="'+conditionColor(payload,payload.conditions[1])+'">'+escText(payload.conditions[1])+'_up</text><text x="'+(mid+6)+'" y="37" font-size="11" fill="'+conditionColor(payload,payload.conditions[0])+'">'+escText(payload.conditions[0])+'_up</text>');
    let y=48;
    for(const p of shown){
      const value=Number(p.change),active=p.prefix===activePrefix,color=pointColor(payload,p);
      parts.push('<text class="rank-name" data-prefix="'+escText(p.prefix)+'" x="5" y="'+(y+8)+'" font-size="8.5" fill="'+(active?color:'#526176')+'">'+escText(motifLabel(p).slice(0,20))+'</text>');
      parts.push('<g clip-path="url(#'+clip.id+')"><rect class="rank-bar'+(active?' selected':'')+'" data-prefix="'+escText(p.prefix)+'" x="'+Math.min(sx(0),sx(value))+'" y="'+y+'" width="'+Math.abs(sx(value)-sx(0))+'" height="9" fill="'+color+'" fill-opacity="'+(active?.95:.72)+'"/></g>');
      if(low<=0&&high>=0)parts.push('<text x="'+(sx(0)+(value>=0?-3:3))+'" y="'+(y+8)+'" font-size="8" text-anchor="'+(value>=0?'end':'start')+'">'+fmt(value)+'</text>');
      y+=rowH+gap;
    }
    if(low<=0&&high>=0)parts.push('<line x1="'+sx(0)+'" x2="'+sx(0)+'" y1="32" y2="'+bottom+'" class="axis"/>');
    parts.push('<line x1="'+left+'" x2="'+right+'" y1="'+bottom+'" y2="'+bottom+'" class="axis"/>');
    niceTicks(low,high,5).forEach(v=>parts.push('<line x1="'+sx(v)+'" x2="'+sx(v)+'" y1="'+bottom+'" y2="'+(bottom+3)+'" class="axis"/><text x="'+sx(v)+'" y="'+(bottom+13)+'" class="tick" text-anchor="middle">'+fmt(v)+'</text>'));
    parts.push('<text x="'+((left+right)/2)+'" y="'+(h-5)+'" class="axis-label" text-anchor="middle">'+escText(payload.change_label||'Differential footprint score')+'</text>');
    if(!shown.length)parts.push('<text x="170" y="130" class="tick" text-anchor="middle">No motifs meet the display filter</text>');
    return parts.join('')+'</svg>';
  }
  function drawVolcano(payload,slot){
    const [xmin,xmax]=rangeFor('volcanoX',slot),[ymin,ymax]=rangeFor('volcanoY',slot);
    const w=430,h=270,left=50,top=15,iw=366,ih=224;
    const sx=v=>left+(v-xmin)/(xmax-xmin)*iw,sy=v=>top+ih-(v-ymin)/(ymax-ymin)*ih;
    const clip=clipping(left,top,iw,ih),parts=[frame(w,h,'volcano-svg',slot,[xmin,xmax],[ymin,ymax]),clip.markup];
    niceTicks(ymin,ymax,5).forEach(t=>parts.push('<line x1="'+left+'" y1="'+sy(t)+'" x2="'+(left+iw)+'" y2="'+sy(t)+'" class="grid"/><text x="'+(left-6)+'" y="'+(sy(t)+3)+'" class="tick" text-anchor="end">'+fmt(t)+'</text>'));
    niceTicks(xmin,xmax,5).forEach(t=>parts.push('<line x1="'+sx(t)+'" y1="'+top+'" x2="'+sx(t)+'" y2="'+(top+ih)+'" class="grid"/><text x="'+sx(t)+'" y="'+(top+ih+13)+'" class="tick" text-anchor="middle">'+fmt(t)+'</text>'));
    if(xmin<=0&&xmax>=0)parts.push('<line x1="'+sx(0)+'" x2="'+sx(0)+'" y1="'+top+'" y2="'+(top+ih)+'" class="zero"/>');
    parts.push('<line x1="'+left+'" x2="'+(left+iw)+'" y1="'+(top+ih)+'" y2="'+(top+ih)+'" class="axis"/><line x1="'+left+'" x2="'+left+'" y1="'+top+'" y2="'+(top+ih)+'" class="axis"/>');
    let clipped=false;
    parts.push('<g clip-path="url(#'+clip.id+')">');
    for(const p of payload.points||[]){
      const yy=plotControls.negLog10P(p);
      if(!plotControls.numeric(p.change)||!Number.isFinite(yy)||(view.hide&&!pass(p)))continue;
      const active=p.prefix===activePrefix&&pass(p),color=pointColor(payload,p);
      clipped ||= p.change<xmin||p.change>xmax||yy<ymin||yy>ymax;
      parts.push('<circle class="pt'+(active?' selected':'')+'" data-prefix="'+escText(p.prefix)+'" cx="'+sx(Number(p.change))+'" cy="'+sy(yy)+'" r="'+(active?5:1.8)+'" fill="'+color+'" fill-opacity="'+(pass(p)?.78:.45)+'" stroke="'+(active?'#111827':'none')+'"><title>'+escText(motifLabel(p))+'</title></circle>');
    }
    parts.push('</g><text x="215" y="9" font-size="7" text-anchor="middle">'+escText(caption())+(clipped?' · Range clips data':'')+'</text>');
    parts.push('<text x="'+left+'" y="'+(top+ih-5)+'" font-size="11" fill="'+conditionColor(payload,payload.conditions[1])+'">'+escText(payload.conditions[1])+'_up</text><text x="'+(left+iw)+'" y="'+(top+ih-5)+'" text-anchor="end" font-size="11" fill="'+conditionColor(payload,payload.conditions[0])+'">'+escText(payload.conditions[0])+'_up</text>');
    parts.push('<text x="'+(left+iw/2)+'" y="268" class="axis-label" text-anchor="middle">'+escText(payload.change_label||'Differential footprint score')+'</text><text x="12" y="127" transform="rotate(-90 12 127)" class="axis-label" text-anchor="middle">−log10(p-value)</text>');
    return parts.join('')+'</svg>';
  }
  function pathD(profile,x,sx,sy){
    let open=false;
    return profile.map((v,i)=>{
      if(!plotControls.numeric(v)||!Number.isFinite(x[i])){open=false;return '';}
      const cmd=open?'L':'M';open=true;
      return cmd+sx(x[i]).toFixed(2)+','+sy(Number(v)).toFixed(2);
    }).join(' ');
  }
  function drawAggregate(payload,prefix,index,slot){
    const motif=aggregateByPrefix(payload,prefix),domain=rangeFor('aggregateY',slot);
    if(!motif||failures.has(index+':'+prefix))return frame(180,180,'aggregate-panel',slot,null,domain)+'<text x="90" y="90" font-size="9" text-anchor="middle">No aggregate profile</text></svg>';
    const axis=payload.aggregate?.x||[],keep=axis.map((v,i)=>({v,i})).filter(p=>p.v>=-60&&p.v<=60);
    const x=keep.map(p=>p.v),samples=aggregateSamples(motif,index,payload),[ymin,ymax]=domain;
    const left=30,top=4,iw=140,ih=158,sx=v=>left+(v+60)/120*iw,sy=v=>top+ih-(v-ymin)/(ymax-ymin)*ih;
    const clip=clipping(left,top,iw,ih),parts=[frame(180,180,'aggregate-panel',slot,null,domain),clip.markup];
    niceTicks(ymin,ymax,4).forEach(v=>parts.push('<text x="26" y="'+(sy(v)+3)+'" class="tick" text-anchor="end">'+fmt(v)+'</text>'));
    [-60,0,60].forEach(v=>parts.push('<line x1="'+sx(v)+'" x2="'+sx(v)+'" y1="'+top+'" y2="'+(top+ih)+'" class="grid"/><text x="'+sx(v)+'" y="174" class="tick" text-anchor="'+(v<0?'start':v>0?'end':'middle')+'">'+v+'</text>'));
    parts.push('<line x1="'+sx(0)+'" x2="'+sx(0)+'" y1="4" y2="162" class="zero"/><path d="M30,4 V162 H170" fill="none" class="axis"/>');
    parts.push('<text x="35" y="157" class="tick">'+(plotControls.numeric(motif.n_sites)?Number(motif.n_sites):'NA')+'</text><g clip-path="url(#'+clip.id+')">');
    samples.sort((a,b)=>Number(b.fp_score||0)-Number(a.fp_score||0)).forEach(s=>{
      const st=s.style,profile=keep.map(p=>s.profile?.[p.i]);
      parts.push('<path data-sample="'+escText(s.name)+'" d="'+pathD(profile,x,sx,sy)+'" fill="none" stroke="'+st.color+'" stroke-width="'+lineWidth(st.width,.7)+'"'+dashAttr(st.type)+' stroke-opacity="'+alpha(st.alpha,.9)+'"><title>'+escText(sampleDisplayName(s,s.condition))+'</title></path>');
    });
    return parts.join('')+'</g></svg>';
  }
  function renderAggregates(){
    const overridden=review.comparisons.some((_,i)=>ranges.has(rangeKey('aggregateY',String(i))));
    host.querySelector('.aggregate-head .sub').textContent=overridden?'Individual Y limits applied':'All comparisons: common Y scale';
    aggregateGrid.innerHTML=slotComparisons.map((index,slot)=>{
      const p=compPayload(index),passText=pass(pointByPrefix(p,activePrefix)||{})?'':' · Does not meet display filter';
      const motif=aggregateByPrefix(p,activePrefix),domain=rangeFor('aggregateY',slot);
      const clipped=motif&&aggregateSamples(motif,index,p).some(s=>(s.profile||[]).some((v,i)=>p.aggregate?.x?.[i]>=-60&&p.aggregate.x[i]<=60&&plotControls.numeric(v)&&(v<domain[0]||v>domain[1])));
      const note=passText+(clipped?' · Range clips data':'');
      return '<div class="aggregate-tile"><div class="aggregate-tile-label" title="'+escText(compLabel(index)+note)+'">Comparison '+(slot+1)+escText(note)+'</div><div class="aggregate-plot">'+drawAggregate(p,activePrefix,index,slot)+'</div>'+aggregateLegendHtml(p,activePrefix,index)+'</div>';
    }).join('');
  }
  async function ensureProfiles(){
    const prefix=activePrefix;
    async function load(index){
      const key=index+':'+prefix,p=compPayload(index),m=aggregateByPrefix(p,prefix);
      if(!m||profiles.has(key)||failures.has(key))return;
      if(m.conditions?.some(c=>c.samples?.some(s=>Array.isArray(s.profile))))return;
      if(!profileLoads.has(key))profileLoads.set(key,(async()=>{
        try{
          const context=await comparisonContext(index,0),record=await withView(context,()=>profileRecord(prefix));
          profiles.set(key,{...m,conditions:m.conditions.map(c=>({...c,samples:c.samples.map(s=>({...s,profile:record.samples[s.name]}))}))});
        }catch(error){failures.set(key,error.message);}
      })());
      await profileLoads.get(key);
    }
    // All includes hidden comparisons. Cache profiles and bound concurrent I/O.
    let next=0;
    await Promise.all(Array.from({length:Math.min(4,review.comparisons.length)},async()=>{
      while(next<review.comparisons.length)await load(next++);
    }));
  }
  function refreshRangeControls(){
    const select=byId('edit-panel');
    select.innerHTML='<option value="all">All</option>'+review.comparisons.map((_,index)=>'<option value="'+index+'">'+escText(compLabel(index))+'</option>').join('');
    select.value=editing;
    select.onchange=()=>{editing=select.value;refreshRangeControls();};
    const target=byId('ranges');target.replaceChildren();
    for(const [kind,label] of Object.entries({volcanoX:'Volcano X',volcanoY:'Volcano Y',rankX:'Waterfall X',aggregateY:'Aggregate Y'})){
      const selected=editing,key=rangeKey(kind,selected);
      const custom=!!(ranges.has(key)?ranges.get(key):ranges.get(rangeKey(kind,'all')));
      const editor=plotControls.axisEditor(label,targetRange(kind,selected),custom,next=>{
        setRange(kind,selected,next);
        renderAll(false,false);
        if(!next)refreshRangeControls();
      },kind==='volcanoY');
      editor.querySelector('button').textContent='Auto scale '+label;
      target.append(editor);
    }
  }
  async function renderAll(refreshStyles=true,refreshRanges=true){
    const request=++generation;
    host.dataset.ready='false';
    host.querySelectorAll('.export-stack button').forEach(b=>b.disabled=true);
    await ensureProfiles();
    if(request!==generation)return;
    setPanelGridShape();
    if(refreshStyles)renderSampleStyles();
    renderComparisons();renderSelected();renderAggregates();
    if(refreshRanges)refreshRangeControls();
    const missing=review.comparisons.filter((_,i)=>!aggregateByPrefix(compPayload(i),activePrefix)||failures.has(i+':'+activePrefix)).length;
    const status='Display only: '+caption()+'. '+review.comparisons.length+' comparisons; '+slotComparisons.length+' displayed.'+(missing?' '+missing+' comparisons have no aggregate profile for this motif; excluded from aggregate autoscaling.':'')+' Original statistics are unchanged.';
    byId('status').textContent=status;reportDetail.textContent=status;
    host.querySelectorAll('.export-stack button').forEach(b=>b.disabled=false);
    host.dataset.ready='true';
  }
  function syncRankRows(source){
    const value=Math.max(1,Math.min(Number(rankRowsSel.max)||200,Math.floor(Number(source.value)||20)));
    rankRowsSel.value=value;rankRowsSlider.value=value;renderAll(false);
  }
  function packedSvg(nodes){
    const cols=Math.min(panelColumnCount(slotComparisons.length),nodes.length)||1,width=340;
    const heights=nodes.map(s=>width*s.viewBox.baseVal.height/s.viewBox.baseVal.width+18),rowHeights=[];
    heights.forEach((h,i)=>rowHeights[Math.floor(i/cols)]=Math.max(rowHeights[Math.floor(i/cols)]||0,h));
    let markup='';
    nodes.forEach((s,i)=>{
      const row=Math.floor(i/cols),y=rowHeights.slice(0,row).reduce((a,b)=>a+b,0),slot=Number(s.dataset.slot);
      markup+='<g transform="translate('+(i%cols*width)+','+y+')"><text x="5" y="10" font-family="Arial" font-size="8">'+escText('Comparison '+(slot+1)+': '+compLabel(slotComparisons[slot]))+'</text><g transform="translate(0,18) scale('+(width/s.viewBox.baseVal.width)+')">'+s.innerHTML+'</g></g>';
    });
    return {markup,width:cols*width,height:rowHeights.reduce((a,b)=>a+b,0)};
  }
  function exportDocument(data){
    return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 '+data.width+' '+data.height+'"><style>'+plotSvgStyle+'</style><rect width="100%" height="100%" fill="white"/>'+data.markup+'</svg>';
  }
  function downloadSvgList(selector,name){
    downloadBlob(new Blob([exportDocument(packedSvg([...host.querySelectorAll(selector)]))],{type:'image/svg+xml;charset=utf-8'}),name);
  }
  function downloadCombinedPanel(){
    const top=packedSvg([...host.querySelectorAll('.rank-svg,.volcano-svg')]),agg=aggregateTilesData();
    downloadBlob(new Blob([exportDocument({width:Math.max(top.width,agg.width),height:top.height+18+agg.height,markup:top.markup+'<g transform="translate(0,'+(top.height+18)+')">'+agg.markup+'</g>'})],{type:'image/svg+xml;charset=utf-8'}),'review_multi_comparisons_classic_panel.svg');
  }
  function setupControls(){
    byId('filter').innerHTML='<details open><summary>Display filters</summary><label>Significance <select id="classic-filter-metric"><option value="fdr">False discovery rate (FDR)</option><option value="pvalue">Raw p-value</option></select></label><label>Maximum <input id="classic-alpha" type="number" min="0" max="1" step="any" value=".05"></label><label>Minimum |ΔFP| <input id="classic-delta" type="number" min="0" step="any" value=".1"></label><button id="classic-apply">Apply filters</button><label><input id="classic-hide" type="checkbox">Hide failing volcano points</label></details>';
    byId('filter-metric').value=view.filter.metric;byId('alpha').value=view.filter.alpha;byId('delta').value=view.filter.delta;byId('hide').checked=view.hide;
    byId('apply').onclick=()=>{
      const alpha=byId('alpha').value,delta=byId('delta').value;
      if(!plotControls.numeric(alpha)||+alpha<0||+alpha>1||!plotControls.numeric(delta)||+delta<0){
        byId('status').textContent='Enter a probability from 0 to 1 and a nonnegative |ΔFP| threshold.';return;
      }
      view.filter={metric:byId('filter-metric').value,alpha:+alpha,delta:+delta};
      $('filter-metric').value=view.filter.metric;$('filter-alpha').value=alpha;$('filter-delta').value=delta;
      renderAll(false);
    };
    byId('hide').onchange=()=>{view.hide=byId('hide').checked;$('hide-failing').checked=view.hide;renderAll(false);};
    byId('layout').onchange=()=>{$('report-layout').value=byId('layout').value;$('report-layout').dispatchEvent(new Event('change'));};
    byId('download-rank').onclick=()=>downloadSvgList('.rank-svg','review_multi_comparisons_classic_barplots.svg');
    byId('download-volcano').onclick=()=>downloadSvgList('.volcano-svg','review_multi_comparisons_classic_volcano.svg');
    byId('download-aggregate').onclick=downloadAggregateTiles;
    byId('download-panel').onclick=downloadCombinedPanel;
    byId('download-logo').onclick=downloadLogoPanel;
    rankRowsSel.oninput=()=>syncRankRows(rankRowsSel);rankRowsSlider.oninput=()=>syncRankRows(rankRowsSlider);
    panelCountSel.onchange=()=>{setPanelCount(panelCountSel.value);renderAll(true);};
  }
  return {async show(data){
    if(!review){
      review=data;initState();
      if(review.comparisons.length<4)panelCountSel.innerHTML=Array.from({length:review.comparisons.length},(_,i)=>'<option value="'+(i+1)+'">'+(i+1)+'</option>').join('');
      panelCountSel.value=slotComparisons.length;
      rankRowsSel.max=rankRowsSlider.max=Math.max(1,...review.comparisons.map(c=>c.payload.points.length));
      setupControls();
    }
    byId('layout').value='classic';
    byId('filter-metric').value=view.filter.metric;byId('alpha').value=view.filter.alpha;byId('delta').value=view.filter.delta;byId('hide').checked=view.hide;
    await renderAll();
  }};
}

window.fpToolsClassic={
  instance:null,request:0,
  async show(){
    const request=++this.request;
    let host=document.getElementById('classic-view');
    if(!host){host=document.createElement('div');host.id='classic-view';document.body.append(host);}
    document.body.classList.add('classic-mode');host.hidden=false;
    let payloads;
    try {payloads=await allComparisonData();}
    catch(error){
      let notice=host.querySelector('[role="alert"]');
      if(!notice){notice=document.createElement('p');notice.setAttribute('role','alert');host.prepend(notice);}
      notice.textContent='Classic view could not load comparison data: '+error.message;
      throw error;
    }
    if(request!==this.request||view.layout!=='classic')return;
    if(!this.instance)this.instance=createClassicView(host);
    await this.instance.show({comparisons:state.metadata.comparisons.map((entry,i)=>({label:entry.label||entry.condition1+' vs '+entry.condition2,payload:payloads[i]}))});
  },
  hide(){++this.request;document.body.classList.remove('classic-mode');const host=document.getElementById('classic-view');if(host)host.hidden=true;}
};
