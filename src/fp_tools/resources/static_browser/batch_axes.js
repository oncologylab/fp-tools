"use strict";
// Adapter for aggregate.batch.v2: retain its sample controls and output schema.
const batchAxes = {ranges:new Map(), shared:false, sharedRange:null, serial:0};
const originalBatchDomain = yDomainForPanel;
function batchDomainKey(index) { return `${index}:${slotPrefixes[index]}`; }
function batchSharedDomain() {
  return fpToolsPlotControls.autoDomain(slotPrefixes.flatMap((prefix,index) =>
    samplesForMotif(motifByPrefix(prefix),index).flatMap(sample=>sample.profile)));
}
yDomainForPanel = function(index) {
  return batchAxes.ranges.get(batchDomainKey(index)) ||
    (batchAxes.shared ? batchAxes.sharedRange || batchSharedDomain() : originalBatchDomain(index));
};
const originalBatchDraw = drawAggregate;
drawAggregate = function(motif,index) {
  const holder=document.createElement('div'); holder.innerHTML=originalBatchDraw(motif,index);
  const svg=holder.firstChild, ns='http://www.w3.org/2000/svg', id=`batch-clip-${++batchAxes.serial}`;
  const defs=document.createElementNS(ns,'defs'), clip=document.createElementNS(ns,'clipPath'), rect=document.createElementNS(ns,'rect');
  clip.id=id; for (const [k,v] of Object.entries({x:50,y:28,width:276,height:274})) rect.setAttribute(k,v);
  clip.append(rect);defs.append(clip);svg.prepend(defs);
  svg.querySelectorAll('path').forEach(path=>path.setAttribute('clip-path',`url(#${id})`));
  const range=yDomainForPanel(index);svg.dataset.yRange=JSON.stringify(range);
  if (samplesForMotif(motif,index).some(s=>s.profile.some(v=>v<range[0]||v>range[1]))) {
    svg.insertAdjacentHTML('beforeend','<text x="170" y="27" text-anchor="middle" font-size="9">Range clips data</text>');
  }
  return svg.outerHTML;
};
function redrawBatchAxes() {
  document.querySelectorAll('.aggregate-panel').forEach((svg,index)=>{
    const holder=document.createElement('div');holder.innerHTML=drawAggregate(motifByPrefix(slotPrefixes[index]),index);svg.replaceWith(holder.firstChild);
  });
}
function installBatchAxes() {
  let controls=document.getElementById('batch-axis-settings');
  if (!controls) {
    controls=document.createElement('section');controls.id='batch-axis-settings'; controls.style.padding='10px';
    const label=document.createElement('label'), input=document.createElement('input');input.type='checkbox';input.id='batch-shared-y';
    label.append(input,document.createTextNode('Shared aggregate Y: displayed plots'));controls.append(label);
    const range=document.createElement('div');range.id='batch-shared-range';controls.append(range);
    document.querySelector('.top-controls').after(controls);
    input.addEventListener('change',()=>{batchAxes.shared=input.checked;renderAll(false);});
    const style=document.createElement('style');style.textContent='.axis-editor{padding:6px;font-size:12px;background:#f5f8fc}.axis-editor label{display:flex;align-items:center;gap:4px;margin:4px 0}.axis-editor input[type=number]{width:100px}.axis-editor input[type=range]{width:100px}.axis-editor [aria-invalid=true]{outline:2px solid #b91c1c}.axis-editor [role=status]{color:#b91c1c}.axis-editor summary{cursor:pointer}';document.head.append(style);
  }
  const shared=document.getElementById('batch-shared-range');shared.replaceChildren();shared.hidden=!batchAxes.shared;
  if(batchAxes.shared) shared.append(fpToolsPlotControls.axisEditor('Shared aggregate Y',batchAxes.sharedRange||batchSharedDomain(),!!batchAxes.sharedRange,next=>{
    batchAxes.sharedRange=next;redrawBatchAxes();if(!next)installBatchAxes();
  }));
  document.querySelectorAll('.aggregate-panel').forEach((svg,index)=>{
    const parent=svg.parentElement;parent.querySelector('.batch-axis-control')?.remove();
    const box=document.createElement('div');box.className='batch-axis-control';
    const key=batchDomainKey(index);
    box.append(fpToolsPlotControls.axisEditor('Aggregate Y',yDomainForPanel(index),batchAxes.ranges.has(key),next=>{
      if(next)batchAxes.ranges.set(key,next);else batchAxes.ranges.delete(key);
      redrawBatchAxes();if(!next)installBatchAxes();
    }));parent.append(box);
  });
}
const originalBatchRender=renderAll;
renderAll=function(...args){originalBatchRender(...args);installBatchAxes();};
