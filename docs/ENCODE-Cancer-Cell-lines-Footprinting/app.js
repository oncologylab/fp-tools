"use strict";

const $ = id => document.getElementById(id);
const state = {metadata:null,raw:null,motifs:[],selected:null,motifIndex:null,profileShards:new Map(),request:0};
const firstColor="#b83b3b",secondColor="#2468a2",neutralColor="#8a97a8";
const svgStyle="text{font-family:Arial,Helvetica,sans-serif}.axis{stroke:#27364a;stroke-width:1}.grid{stroke:#e7edf4;stroke-width:1}.zero{stroke:#7a8798;stroke-width:1;stroke-dasharray:4 4}.tick{fill:#536277;font-size:10px}.axis-label{fill:#26364a;font-size:11px;font-weight:700}.plot-title{fill:#172033;font-size:12px;font-weight:700}";

function esc(value){return String(value??"").replace(/[&<>"']/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));}
function finite(value,fallback=0){const n=Number(value);return Number.isFinite(n)?n:fallback;}
function fmt(value,digits=3){const n=Number(value);if(!Number.isFinite(n))return"NA";const a=Math.abs(n);if(a===0)return"0";if(a<0.001||a>=1000)return n.toExponential(2);return n.toFixed(digits).replace(/\.0+$|(?<=\.[0-9]*[1-9])0+$/g,"");}
function logq(value){return-Math.log10(Math.max(1e-300,finite(value,1)));}
function colorFor(motif){return motif.significant?(motif.effect>=0?firstColor:secondColor):neutralColor;}
function downloadBlob(blob,name){const url=URL.createObjectURL(blob),a=document.createElement("a");a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
async function fetchJson(path){const response=await fetch(path);if(!response.ok)throw new Error(`${path}: HTTP ${response.status}`);return response.json();}
function linePath(values,xScale,yScale,axis){return values.map((value,index)=>`${index?"L":"M"}${xScale(axis[index]).toFixed(2)},${yScale(value).toFixed(2)}`).join(" ");}
function ticks(min,max,count=5){if(!Number.isFinite(min)||!Number.isFinite(max)||min===max)return[min||0];const step=(max-min)/(count-1);return Array.from({length:count},(_,i)=>min+i*step);}

function comparisonEntry(first,second){
  const direct=state.metadata.comparisons.find(item=>item.condition1===first&&item.condition2===second);
  if(direct)return{entry:direct,reversed:false};
  const reverse=state.metadata.comparisons.find(item=>item.condition1===second&&item.condition2===first);
  if(reverse)return{entry:reverse,reversed:true};
  throw new Error(`No comparison for ${first} and ${second}`);
}

function orientedMotif(motif,reversed){
  if(!reversed)return{...motif};
  return{...motif,mean1:motif.mean2,sd1:motif.sd2,mean2:motif.mean1,sd2:motif.sd1,effect:-motif.effect,ci_lower:-motif.ci_upper,ci_upper:-motif.ci_lower,moderated_t:-motif.moderated_t};
}

async function loadComparison(){
  const first=$("condition-1").value,second=$("condition-2").value;
  if(first===second)return;
  const token=++state.request;
  $("status").textContent=`Loading ${first} vs ${second}`;
  const {entry,reversed}=comparisonEntry(first,second);
  const raw=await fetchJson(entry.file);
  if(token!==state.request)return;
  state.raw=raw;
  state.motifs=raw.motifs.map(motif=>orientedMotif(motif,reversed));
  state.selected=state.motifs.slice().sort((a,b)=>a.qvalue-b.qvalue||Math.abs(b.effect)-Math.abs(a.effect))[0]?.prefix||null;
  renderAll();
  $("status").textContent=`${state.motifs.length} motifs; ${first} minus ${second}`;
}

function filteredMotifs(){
  const query=$("search").value.trim().toLowerCase(),sig=$("significant-only").checked;
  return state.motifs.filter(motif=>(!sig||motif.significant)&&(!query||`${motif.name} ${motif.motif_id} ${motif.prefix} ${motif.cluster}`.toLowerCase().includes(query)));
}

function renderTable(){
  const motifs=filteredMotifs().sort((a,b)=>a.qvalue-b.qvalue||Math.abs(b.effect)-Math.abs(a.effect));
  $("motif-count").textContent=`${motifs.length} of ${state.motifs.length}`;
  $("motif-table").innerHTML=motifs.map(motif=>`<button class="motif-row${motif.prefix===state.selected?" selected":""}" data-prefix="${esc(motif.prefix)}" role="option" aria-selected="${motif.prefix===state.selected}"><span class="motif-name" title="${esc(motif.name)} (${esc(motif.motif_id)})">${esc(motif.name)} (${esc(motif.motif_id)})</span><span class="${motif.effect>=0?"positive":"negative"}">${fmt(motif.effect)}</span><span>${fmt(motif.qvalue,2)}</span></button>`).join("");
  $("motif-table").querySelectorAll("[data-prefix]").forEach(node=>node.addEventListener("click",()=>selectMotif(node.dataset.prefix)));
}

function volcanoSvg(){
  const w=640,h=350,m={l:54,r:18,t:12,b:43},iw=w-m.l-m.r,ih=h-m.t-m.b;
  const xmax=Math.max(.01,...state.motifs.map(m=>Math.abs(finite(m.effect))));
  const ymax=Math.max(1,...state.motifs.map(m=>logq(m.qvalue)));
  const sx=x=>m.l+(x+xmax)/(2*xmax)*iw,sy=y=>m.t+ih-y/ymax*ih;
  let parts=[`<svg data-export="volcano" viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg"><style>${svgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`];
  ticks(-xmax,xmax).forEach(t=>parts.push(`<line x1="${sx(t)}" y1="${m.t}" x2="${sx(t)}" y2="${m.t+ih}" class="grid"/><text x="${sx(t)}" y="${h-25}" text-anchor="middle" class="tick">${fmt(t,2)}</text>`));
  ticks(0,ymax).forEach(t=>parts.push(`<line x1="${m.l}" y1="${sy(t)}" x2="${m.l+iw}" y2="${sy(t)}" class="grid"/><text x="${m.l-6}" y="${sy(t)+3}" text-anchor="end" class="tick">${fmt(t,1)}</text>`));
  parts.push(`<line x1="${sx(0)}" y1="${m.t}" x2="${sx(0)}" y2="${m.t+ih}" class="zero"/><line x1="${m.l}" y1="${m.t+ih}" x2="${m.l+iw}" y2="${m.t+ih}" class="axis"/><line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${m.t+ih}" class="axis"/><text x="${m.l+iw/2}" y="${h-5}" text-anchor="middle" class="axis-label">Differential footprint score</text><text x="14" y="${m.t+ih/2}" transform="rotate(-90 14 ${m.t+ih/2})" text-anchor="middle" class="axis-label">-log10(q)</text>`);
  state.motifs.forEach(motif=>parts.push(`<circle data-prefix="${esc(motif.prefix)}" cx="${sx(motif.effect)}" cy="${sy(logq(motif.qvalue))}" r="${motif.prefix===state.selected?4.6:2.7}" fill="${colorFor(motif)}" fill-opacity="${motif.significant?.85:.48}" stroke="${motif.prefix===state.selected?"#111827":"none"}" stroke-width="1.3"><title>${esc(motif.name)} (${esc(motif.motif_id)}): effect ${fmt(motif.effect)}, q ${fmt(motif.qvalue,2)}</title></circle>`));
  parts.push("</svg>");return parts.join("");
}

function rankSvg(){
  const selected=state.motifs.slice().sort((a,b)=>Math.abs(b.effect)-Math.abs(a.effect)).slice(0,20).sort((a,b)=>a.effect-b.effect);
  const w=640,h=350,m={l:150,r:24,t:10,b:38},iw=w-m.l-m.r,ih=h-m.t-m.b,row=ih/Math.max(1,selected.length);
  const xmax=Math.max(.01,...selected.map(motif=>Math.abs(motif.effect))),sx=x=>m.l+(x+xmax)/(2*xmax)*iw;
  let parts=[`<svg data-export="rank" viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg"><style>${svgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`];
  ticks(-xmax,xmax).forEach(t=>parts.push(`<line x1="${sx(t)}" y1="${m.t}" x2="${sx(t)}" y2="${m.t+ih}" class="grid"/><text x="${sx(t)}" y="${h-20}" text-anchor="middle" class="tick">${fmt(t,2)}</text>`));
  parts.push(`<line x1="${sx(0)}" y1="${m.t}" x2="${sx(0)}" y2="${m.t+ih}" class="zero"/><text x="${m.l+iw/2}" y="${h-4}" text-anchor="middle" class="axis-label">Differential footprint score</text>`);
  selected.forEach((motif,index)=>{const y=m.t+index*row+row*.16,height=Math.max(3,row*.68),x0=sx(Math.min(0,motif.effect)),x1=sx(Math.max(0,motif.effect));parts.push(`<text x="${m.l-6}" y="${y+height*.75}" text-anchor="end" class="tick">${esc(motif.name.slice(0,18))}</text><rect data-prefix="${esc(motif.prefix)}" x="${x0}" y="${y}" width="${Math.max(1,x1-x0)}" height="${height}" fill="${colorFor(motif)}" fill-opacity="${motif.significant?.88:.55}" stroke="${motif.prefix===state.selected?"#111827":"none"}"><title>${esc(motif.name)} (${esc(motif.motif_id)}): ${fmt(motif.effect)}</title></rect>`);});
  parts.push("</svg>");return parts.join("");
}

function bindPlotClicks(){document.querySelectorAll("svg [data-prefix]").forEach(node=>node.addEventListener("click",()=>selectMotif(node.dataset.prefix)));}

async function profileRecord(prefix){
  if(!state.motifIndex)state.motifIndex=await fetchJson("data/motif_index.json");
  const shard=state.motifIndex[prefix];if(shard===undefined)throw new Error(`No profile shard for ${prefix}`);
  if(!state.profileShards.has(shard))state.profileShards.set(shard,await fetchJson(`data/profiles/${shard}.json`));
  return state.profileShards.get(shard).motifs[prefix];
}

function meanSd(profiles){
  const rows=profiles.map(values=>values.map(finite)),length=rows[0]?.length||0,mean=[],sd=[];
  for(let i=0;i<length;i++){const values=rows.map(row=>row[i]),avg=values.reduce((a,b)=>a+b,0)/values.length;mean.push(avg);sd.push(values.length>1?Math.sqrt(values.reduce((a,b)=>a+(b-avg)**2,0)/(values.length-1)):0);}
  return{mean,sd};
}

function profileSvg(record,motif){
  const axis=state.metadata.profile_axis,first=$("condition-1").value,second=$("condition-2").value;
  const groups=[{name:first,color:firstColor,samples:state.metadata.conditions.find(c=>c.name===first).samples},{name:second,color:secondColor,samples:state.metadata.conditions.find(c=>c.name===second).samples}];
  const rows=groups.map(group=>{const profiles=group.samples.map(sample=>record.samples[sample]),stats=meanSd(profiles);return{...group,profiles,stats};});
  const all=rows.flatMap(row=>[...row.profiles.flat(),...row.stats.mean.map((value,index)=>value-row.stats.sd[index]),...row.stats.mean.map((value,index)=>value+row.stats.sd[index])]).map(finite),ymin=Math.min(0,...all),ymax=Math.max(.001,...all),pad=(ymax-ymin)*.06||.1;
  const w=640,h=350,m={l:58,r:20,t:20,b:42},iw=w-m.l-m.r,ih=h-m.t-m.b,sx=x=>m.l+(x-axis[0])/(axis[axis.length-1]-axis[0])*iw,sy=y=>m.t+ih-(y-(ymin-pad))/(ymax-ymin+2*pad)*ih;
  let parts=[`<svg data-export="profile" viewBox="0 0 ${w} ${h}" xmlns="http://www.w3.org/2000/svg"><style>${svgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`];
  ticks(ymin,ymax).forEach(t=>parts.push(`<line x1="${m.l}" y1="${sy(t)}" x2="${m.l+iw}" y2="${sy(t)}" class="grid"/><text x="${m.l-6}" y="${sy(t)+3}" text-anchor="end" class="tick">${fmt(t,2)}</text>`));
  [-100,-50,0,50,99].filter(t=>t>=axis[0]&&t<=axis[axis.length-1]).forEach(t=>parts.push(`<line x1="${sx(t)}" y1="${m.t}" x2="${sx(t)}" y2="${m.t+ih}" class="${t===0?"zero":"grid"}"/><text x="${sx(t)}" y="${h-24}" text-anchor="middle" class="tick">${t}</text>`));
  parts.push(`<line x1="${m.l}" y1="${m.t+ih}" x2="${m.l+iw}" y2="${m.t+ih}" class="axis"/><line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${m.t+ih}" class="axis"/><text x="${m.l+iw/2}" y="${h-5}" text-anchor="middle" class="axis-label">Position relative to motif center (bp)</text><text x="14" y="${m.t+ih/2}" transform="rotate(-90 14 ${m.t+ih/2})" text-anchor="middle" class="axis-label">Normalized corrected cut-site signal</text>`);
  rows.forEach((group,groupIndex)=>{const stats=group.stats,upper=stats.mean.map((v,i)=>v+stats.sd[i]),lower=stats.mean.map((v,i)=>v-stats.sd[i]),polygon=[...axis.map((x,i)=>`${sx(x)},${sy(upper[i])}`),...axis.slice().reverse().map((x,j)=>{const i=axis.length-1-j;return`${sx(x)},${sy(lower[i])}`;})].join(" ");parts.push(`<polygon points="${polygon}" fill="${group.color}" fill-opacity=".12"/>`);group.profiles.forEach((profile,index)=>parts.push(`<path d="${linePath(profile,sx,sy,axis)}" fill="none" stroke="${group.color}" stroke-width=".75" stroke-opacity=".48"><title>${esc(group.samples[index])}</title></path>`));parts.push(`<path d="${linePath(stats.mean,sx,sy,axis)}" fill="none" stroke="${group.color}" stroke-width="2.2"/><line x1="${m.l+12+groupIndex*170}" y1="12" x2="${m.l+34+groupIndex*170}" y2="12" stroke="${group.color}" stroke-width="2.2"/><text x="${m.l+40+groupIndex*170}" y="15" class="tick">${esc(group.name)} mean +/- SD</text>`);});
  parts.push("</svg>");return parts.join("");
}

function renderDetails(motif){
  const first=$("condition-1").value,second=$("condition-2").value;
  $("selected-label").textContent=`${motif.name} (${motif.motif_id})`;
  const items=[[`${first} mean`,motif.mean1],[`${second} mean`,motif.mean2],["Effect",motif.effect],["95% CI",`${fmt(motif.ci_lower)} to ${fmt(motif.ci_upper)}`],["Moderated t",motif.moderated_t],["Reference distribution",motif.normal_limit?"Normal limit":`t, df ${fmt(motif.moderated_df)}`],["p-value",motif.pvalue],["BH q-value",motif.qvalue]];
  $("statistics").innerHTML=items.map(([label,value])=>`<div><dt>${esc(label)}</dt><dd>${typeof value==="number"?fmt(value):esc(value)}</dd></div>`).join("");
}

async function renderProfile(){
  const motif=state.motifs.find(item=>item.prefix===state.selected);if(!motif)return;
  const prefix=motif.prefix;$("profile").innerHTML=`<svg viewBox="0 0 640 350"><text x="320" y="175" text-anchor="middle" class="axis-label">Loading profile...</text></svg>`;
  try{const record=await profileRecord(prefix);if(state.selected!==prefix)return;$("profile").innerHTML=profileSvg(record,motif);const used=record.n_profile_sites||record.n_sites,siteText=used<record.n_sites?`${used.toLocaleString()} of ${record.n_sites.toLocaleString()} motif sites`:`${record.n_sites.toLocaleString()} motif sites`;$("selected-detail").textContent=`${siteText}; thin lines are replicates`;}
  catch(error){$("profile").innerHTML=`<svg viewBox="0 0 640 350"><text x="320" y="175" text-anchor="middle" class="axis-label">${esc(error.message)}</text></svg>`;}
}

function renderAll(){
  const first=$("condition-1").value,second=$("condition-2").value,sig=state.motifs.filter(m=>m.significant);
  $("comparison-summary").textContent=`${first} minus ${second}; ${sig.length} motifs with q < 0.05`;
  renderTable();$("volcano").innerHTML=volcanoSvg();$("rank").innerHTML=rankSvg();bindPlotClicks();
  const motif=state.motifs.find(item=>item.prefix===state.selected);if(motif)renderDetails(motif);renderProfile();
}

function selectMotif(prefix){state.selected=prefix;renderAll();}

function comparisonTsv(){
  const first=$("condition-1").value,second=$("condition-2").value;
  const columns=["condition1","condition2","prefix","name","motif_id","cluster","n_sites","mean1","sd1","mean2","sd2","effect","ci_lower","ci_upper","moderated_t","moderated_df","normal_limit","pvalue","qvalue","significant"];
  const rows=[columns.join("\t")];state.motifs.forEach(motif=>rows.push([first,second,...columns.slice(2).map(key=>motif[key])].join("\t")));return rows.join("\n")+"\n";
}

function downloadTsv(){downloadBlob(new Blob([comparisonTsv()],{type:"text/tab-separated-values;charset=utf-8"}),`${$("condition-1").value}_vs_${$("condition-2").value}_fp_tools.tsv`);}

function combinedSvg(){
  const volcano=$("volcano").querySelector("svg"),rank=$("rank").querySelector("svg"),profile=$("profile").querySelector("svg");if(!volcano||!rank||!profile)return;
  const motif=state.motifs.find(item=>item.prefix===state.selected),first=$("condition-1").value,second=$("condition-2").value;
  return{svg:`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1280 720" width="1280" height="720" font-family="Arial,Helvetica,sans-serif"><style>${svgStyle}</style><rect width="100%" height="100%" fill="#fff"/><g transform="translate(0,0)">${volcano.innerHTML}</g><g transform="translate(640,0)">${rank.innerHTML}</g><g transform="translate(0,360)">${profile.innerHTML}</g><g transform="translate(660,390)"><text x="0" y="0" class="plot-title">${esc(motif.name)} (${esc(motif.motif_id)})</text><text x="0" y="28" class="axis-label">${esc(first)} minus ${esc(second)}</text><text x="0" y="54" class="tick">Effect ${fmt(motif.effect)}; 95% CI ${fmt(motif.ci_lower)} to ${fmt(motif.ci_upper)}</text><text x="0" y="76" class="tick">Moderated t ${fmt(motif.moderated_t)}; BH q ${fmt(motif.qvalue,2)}; ${motif.normal_limit?"normal limit":`t df ${fmt(motif.moderated_df)}`}</text></g></svg>`,name:`${first}_vs_${second}_${motif.prefix}`};
}

function saveFigure(){
  const figure=combinedSvg();if(!figure)return;
  const format=$("figure-format").value;
  if(format==="svg"){downloadBlob(new Blob([figure.svg],{type:"image/svg+xml;charset=utf-8"}),`${figure.name}.svg`);return;}
  if(format==="png"){
    const image=new Image(),url=URL.createObjectURL(new Blob([figure.svg],{type:"image/svg+xml;charset=utf-8"}));
    image.onload=()=>{const canvas=document.createElement("canvas");canvas.width=2560;canvas.height=1440;const context=canvas.getContext("2d");context.fillStyle="#fff";context.fillRect(0,0,canvas.width,canvas.height);context.drawImage(image,0,0,canvas.width,canvas.height);URL.revokeObjectURL(url);canvas.toBlob(blob=>{if(blob)downloadBlob(blob,`${figure.name}.png`);},"image/png");};
    image.onerror=()=>{URL.revokeObjectURL(url);$("status").textContent="PNG export could not be rendered";};image.src=url;return;
  }
  const page=window.open("","_blank");if(!page){$("status").textContent="Allow pop-ups to use Print / PDF";return;}
  page.document.write(`<!doctype html><html><head><title>${esc(figure.name)}</title><style>@page{size:landscape;margin:0}html,body{margin:0;width:100%;height:100%;font-family:Arial,Helvetica,sans-serif}svg{display:block;width:100%;height:100%}</style></head><body>${figure.svg}<script>window.onload=()=>{window.print();}<\/script></body></html>`);page.document.close();
}

function keepDistinct(changed,other){if($(changed).value!==$(other).value)return;const replacement=state.metadata.conditions.map(c=>c.name).find(name=>name!==$(changed).value);$(other).value=replacement;}

async function init(){
  try{
    state.metadata=await fetchJson("data/metadata.json");
    const names=state.metadata.conditions.map(item=>item.name),options=names.map(name=>`<option value="${esc(name)}">${esc(name)}</option>`).join("");
    $("condition-1").innerHTML=options;$("condition-2").innerHTML=options;$("condition-1").value="A549";$("condition-2").value="HCT116";
    $("release").textContent=`${state.metadata.sample_count} replicates | ${state.metadata.comparison_count} comparisons | ${state.metadata.motif_count} motifs | ${state.metadata.genome} | ${state.metadata.release_date}`;
    $("download-all").href=state.metadata.downloads.all_results;
    $("condition-1").addEventListener("change",()=>{keepDistinct("condition-1","condition-2");loadComparison();});
    $("condition-2").addEventListener("change",()=>{keepDistinct("condition-2","condition-1");loadComparison();});
    $("swap").addEventListener("click",()=>{const value=$("condition-1").value;$("condition-1").value=$("condition-2").value;$("condition-2").value=value;loadComparison();});
    $("search").addEventListener("input",renderTable);$("significant-only").addEventListener("change",renderTable);$("download-tsv").addEventListener("click",downloadTsv);$("save-figure").addEventListener("click",saveFigure);
    await loadComparison();
  }catch(error){$("status").textContent=`Could not load resource: ${error.message}`;throw error;}
}

init();
