"use strict";

// Presentation state is deliberately separate from the embedded scientific payload.
const view = {
  filter: {metric: 'fdr', alpha: .05, delta: .1}, hide: false,
  layout: 'single', slots: [0, 1, 2, 3], count: 4,
  shared: {volcanoX: false, volcanoY: false, rankX: false, aggregateY: false},
  manual: new Map(), groupData: null, sharedAggregate: null,
  rendered: new Map(), styles: new Map(), generation: 0,
};

function displayPass(point) { return plotControls.passesDisplayFilter(point, view.filter); }
function siteCountLabel(value) { return plotControls.numeric(value) ? Number(value).toLocaleString() : 'NA'; }
function setFigureBusy(busy) {
  for (const id of ['download-rank','download-volcano','download-aggregate','download-panel']) $(id).disabled = busy;
}
function filterCaption() {
  return `${view.filter.metric === 'fdr' ? 'FDR' : 'p-value'} ≤ ${view.filter.alpha}; |ΔFP| ≥ ${view.filter.delta}`;
}
function comparisonKey() {
  return `${state.entry?.comparison || state.comparisonIndex}:${state.first}:${state.second}`;
}
function rangeKey(kind, prefix = '') { return `${state.plotInstance || 'main'}:${comparisonKey()}:${kind}:${prefix}`; }
function plotRange(kind, automatic, prefix = '') {
  const manual = view.manual.get(rangeKey(kind, prefix));
  if (manual) return manual;
  if (kind === 'aggregateY' && view.shared.aggregateY)
    return view.manual.get('aggregateY:shared') || view.sharedAggregate || automatic;
  if (!view.shared[kind] || !view.groupData) return automatic;
  const values = [];
  view.groupData.forEach(payload => {
    let points = payload.points || [];
    if (kind === 'rankX') {
      const ranked = plotControls.rankMotifs(points.filter(displayPass), rankMode(), Number($('rank-rows').value));
      points = [...ranked.negative, ...ranked.positive];
    }
    points.forEach(p => values.push(kind === 'volcanoY' ? plotControls.negLog10P(p) :
      kind === 'rankX' ? plotControls.rankMetric(p, rankMode()) : p.change));
  });
  return plotControls.autoDomain(values, kind === 'volcanoY' ? 'positive' : 'symmetric');
}

function withView(context, callback) {
  const original = {...state};
  Object.assign(state, context);
  try { return callback(); } finally { Object.assign(state, original); }
}

function filterNote(svg, points, clipped = false) {
  svg.dataset.filter = filterCaption();
  const note = document.createElementNS('http://www.w3.org/2000/svg', 'text');
  note.setAttribute('x', String(svg.viewBox.baseVal.width / 2));
  note.setAttribute('y', '46'); note.setAttribute('text-anchor', 'middle');
  note.setAttribute('font-size', '10'); note.setAttribute('fill', '#475569');
  note.textContent = `${points.filter(displayPass).length}/${points.length} pass · ${filterCaption()}${clipped ? ' · Range clips data' : ''}`;
  svg.append(note);
}

let plotClipSerial = 0;
function clipPlot(svg, selector, box) {
  const ns = 'http://www.w3.org/2000/svg', id = `fp-clip-${++plotClipSerial}`;
  const defs = document.createElementNS(ns, 'defs'), clip = document.createElementNS(ns, 'clipPath');
  clip.id = id;
  const rect = document.createElementNS(ns, 'rect');
  ['x', 'y', 'width', 'height'].forEach((key, i) => rect.setAttribute(key, box[i]));
  clip.append(rect); defs.append(clip); svg.prepend(defs);
  svg.querySelectorAll(selector).forEach(node => node.setAttribute('clip-path', `url(#${id})`));
}

function aggregateAutoRange(record) {
  const values = [];
  [state.first, state.second].forEach(condition => conditionSamples(condition).forEach((sample, i) => {
    if (!sampleStyle(sample, condition, i).visible) return;
    (record.samples[sample] || []).forEach((v, j) => {
      const x = state.profileAxis[j];
      if (x >= -60 && x <= 60 && plotControls.numeric(v)) values.push(Number(v));
    });
  }));
  return plotControls.autoDomain(values);
}

function registerPlot(key, item) { view.rendered.set(key, item); }
function redrawRanges() {
  for (const item of view.rendered.values()) {
    if (!item.node.isConnected) continue;
    withView(item.context, () => {
      if (item.kind === 'volcano') renderVolcano(item.node);
      else if (item.kind === 'rank') drawRank(item.node);
      else {
        const holder = document.createElement('div');
        holder.innerHTML = profileSvg(item.record, item.motif, item.index);
        item.node.replaceWith(holder.firstElementChild); item.node = holder.firstElementChild || item.node;
        // replaceWith moves the child out of holder.
        item.node = item.parent.querySelector('.aggregate-panel');
      }
    });
  }
}

function axisControls(parent, kind, automatic, prefix = '') {
  const key = rangeKey(kind, prefix),
    label = {volcanoX:'Volcano X', volcanoY:'Volcano Y', rankX:'Waterfall X', aggregateY:'Aggregate Y'}[kind],
    domain = plotRange(kind, automatic, prefix), context = {...state};
  parent.append(plotControls.axisEditor(label, domain, view.manual.has(key), next => {
    if (next) view.manual.set(key, next); else view.manual.delete(key);
    redrawRanges();
    if (!next) {
      const replacement = withView(context, () => {
        const holder = document.createElement('div'); axisControls(holder, kind, automatic, prefix); return holder.firstChild;
      });
      [...parent.children].find(child => child.contains(document.activeElement))?.replaceWith(replacement);
    }
  }, kind === 'volcanoY'));
}

async function allComparisonData() {
  if (view.groupData) return view.groupData;
  const results = [];
  // Bound concurrent requests and never substitute an incomplete shared domain.
  const records = state.metadata.comparisons;
  for (let start = 0; start < records.length; start += 4) {
    const batch = await Promise.all(records.slice(start, start + 4).map(entry =>
      state.mode === 'embedded' ? entry.payload : fetchGzipJsonCached(entry.file)));
    results.push(...batch);
  }
  view.groupData = results;
  return results;
}

function setupReportOptions() {
  const box = document.createElement('details'); box.id = 'plot-settings';
  box.className = 'plot-settings'; box.open = true;
  box.innerHTML = `<summary>Plot ranges and display filters</summary><div class="settings-row">
    <label>View <select id="report-layout"><option value="single">Single comparison</option><option value="side">Side by side</option><option value="classic">Classic</option></select></label>
    <label id="comparison-count-label">Comparisons <input id="comparison-count" type="number" min="2" max="8" value="4"></label>
    <label>Significance <select id="filter-metric"><option value="fdr">False discovery rate (FDR)</option><option value="pvalue">Raw p-value</option></select></label>
    <label>Maximum <input id="filter-alpha" type="number" min="0" max="1" step="any" value="0.05"></label>
    <label>Minimum |ΔFP| <input id="filter-delta" type="number" min="0" step="any" value="0.1"></label>
    <button id="apply-filters" type="button">Apply filters</button>
    <label><input id="hide-failing" type="checkbox"> Hide failing volcano points</label>
    </div><div class="settings-row" id="shared-axes"></div>
    <p id="filter-status" role="status"></p><div id="main-axis-controls" class="settings-row"></div>
    <div id="shared-aggregate-controls"></div>`;
  $('dashboard').before(box);
  const side = document.createElement('section'); side.id = 'comparison-grid'; $('dashboard').after(side);
  $('report-layout').value = ['side','classic'].includes(bootstrap.defaultView) ? bootstrap.defaultView : 'single';
  view.layout = $('report-layout').value;
  $('report-layout').disabled = state.metadata.comparisons.length < 2;
  $('comparison-count').max = Math.min(8, state.metadata.comparisons.length);
  view.count = Math.min(4, state.metadata.comparisons.length);
  $('comparison-count').value = view.count;
  $('report-layout').addEventListener('change', () => {view.layout = $('report-layout').value; renderAll(false);});
  $('comparison-count').addEventListener('change', () => {
    view.count = Math.max(2, Math.min(8, state.metadata.comparisons.length, Math.floor(Number($('comparison-count').value) || 4)));
    $('comparison-count').value = view.count; renderAll(false);
  });
  $('apply-filters').addEventListener('click', () => {
    const alpha = $('filter-alpha').value, delta = $('filter-delta').value;
    if (!plotControls.numeric(alpha) || !plotControls.numeric(delta) || +alpha < 0 || +alpha > 1 || +delta < 0) {
      $('filter-status').textContent = 'Enter a probability from 0 to 1 and a nonnegative |ΔFP| threshold.'; return;
    }
    view.filter = {metric: $('filter-metric').value, alpha:+alpha, delta:+delta}; renderAll(false);
  });
  $('hide-failing').addEventListener('change', () => {view.hide = $('hide-failing').checked; renderAll(false);});
  Object.entries({volcanoX:'Shared volcano X: all comparisons', volcanoY:'Shared volcano Y: all comparisons', rankX:'Shared waterfall X: all comparisons', aggregateY:'Shared aggregate Y: displayed plots'}).forEach(([kind, text]) => {
    const label = document.createElement('label'), input = document.createElement('input');
    input.type = 'checkbox'; input.id = `shared-${kind}`; label.append(input, document.createTextNode(text)); $('shared-axes').append(label);
    input.addEventListener('change', async () => {
      if (input.checked && kind !== 'aggregateY') {
        input.disabled = true; $('filter-status').textContent = 'Loading all comparison statistics for shared scaling…';
        try { await allComparisonData(); } catch (error) {
          input.checked = false; $('filter-status').textContent = `Shared scale unavailable: ${error.message}`; input.disabled = false; return;
        }
        input.disabled = false;
      }
      view.shared[kind] = input.checked; renderAll(false);
    });
  });
  const rowLabel = $('rank-rows').closest('label');
  rowLabel.firstChild.textContent = 'Top motifs (total) ';
  for (const id of ['rank-rows', 'rank-rows-slider']) {
    $(id).min = 1; $(id).max = Math.max(...state.metadata.comparisons.map(c => c.payload?.points?.length || c.motifs || 20));
  }
}

function refreshViewControls() {
  if (!$('plot-settings')) return;
  const side = view.layout === 'side';
  $('dashboard').hidden = side; $('comparison-grid').hidden = !side;
  $('comparison-count-label').hidden = !side;
  $('selected-grid').previousElementSibling.querySelector('.section-title').textContent = side
    ? 'Selected motifs · click one to compare across panels' : 'Selected motifs';
  $('main-axis-controls').hidden = side;
  $('filter-status').textContent = `Display only: ${filterCaption()}. Original statistics and result-table downloads are unchanged.${side ? ' The active selected motif is shown across comparisons.' : ''}`;
  if (side) { renderComparisonGrid().catch(showError); return; }
  view.generation++;
  $('comparison-grid').replaceChildren();
  const controls = $('main-axis-controls'); controls.replaceChildren();
  for (const [kind, values, mode] of [
    ['volcanoX', state.motifs.map(p => p.effect), 'symmetric'],
    ['volcanoY', state.motifs.map(p => p.neglog10p), 'positive'],
    ['rankX', rankedVisible().map(p => plotControls.rankMetric(p, rankMode())), 'symmetric'],
  ]) axisControls(controls, kind, plotControls.autoDomain(values, mode));
  registerPlot('main-volcano', {kind:'volcano', node:$('chart'), context:{...state}});
  registerPlot('main-rank', {kind:'rank', node:$('rank-chart'), context:{...state}});
}

function rankedVisible() {
  const rows = plotControls.rankMotifs(state.motifs.filter(displayPass), rankMode(), Number($('rank-rows').value));
  return [...rows.negative, ...rows.positive];
}

function prepareAggregateRanges(items) {
  const values = items.flatMap(item => withView(item.context, () => aggregateAutoRange(item.record)));
  view.sharedAggregate = values.length ? [Math.min(...values), Math.max(...values)] : [0, 1];
  const holder = $('shared-aggregate-controls');
  if (!holder) return;
  holder.replaceChildren(); holder.hidden = !view.shared.aggregateY;
  if (view.shared.aggregateY) holder.append(plotControls.axisEditor('Shared aggregate Y',
    view.manual.get('aggregateY:shared') || view.sharedAggregate, view.manual.has('aggregateY:shared'), next => {
      if (next) view.manual.set('aggregateY:shared', next); else view.manual.delete('aggregateY:shared');
      redrawRanges(); if (!next) prepareAggregateRanges(items);
    }));
}

function attachAggregateControls(item) {
  withView(item.context, () => {
    item.parent.querySelector('.aggregate-axis-controls')?.remove();
    const controls = document.createElement('div'); controls.className = 'aggregate-axis-controls';
    axisControls(controls, 'aggregateY', aggregateAutoRange(item.record), `${item.motif.prefix}:${item.index}`);
    item.parent.append(controls);
    const pass = displayPass(item.motif), note = document.createElement('p'); note.className = 'aggregate-filter-note';
    note.textContent = pass ? 'Meets display filter' : 'Does not meet display filter; shown for inspection';
    controls.append(note);
    registerPlot(`aggregate:${comparisonKey()}:${item.index}`, {...item, kind:'aggregate'});
  });
}

async function comparisonContext(index, slot) {
  const entry = state.metadata.comparisons[index];
  const payload = state.mode === 'embedded' ? entry.payload : await fetchGzipJsonCached(entry.file);
  const aggregate = new Map((payload.aggregate?.motifs || []).map(m => [m.prefix, m]));
  const [first, second] = payload.conditions;
  const key = entry.comparison;
  if (!view.styles.has(key)) view.styles.set(key, new Map());
  const context = {...state, entry, payload, aggregate, first, second, comparisonIndex:index, plotInstance:`side-${slot}`,
    profileAxis:payload.aggregate?.x || [], sampleStyles:view.styles.get(key),
    colors:{...state.colors}};
  context.motifs = withView(context, () => payload.points.map(p => orientedMotif(p, false)));
  return context;
}

async function renderComparisonGrid() {
  const token = ++view.generation, grid = $('comparison-grid');
  setFigureBusy(true);
  for (let i = view.slots.length; i < view.count; i++) view.slots.push(i);
  const contexts = await Promise.all(view.slots.slice(0, view.count).map(comparisonContext));
  if (token !== view.generation) return;
  grid.replaceChildren();
  const activePrefix = state.selected[state.active] || state.selected[0], items = [];
  for (let slot = 0; slot < contexts.length; slot++) {
    const context = contexts[slot], card = document.createElement('section'); card.className = 'comparison-card';
    const select = document.createElement('select'); select.setAttribute('aria-label', `Comparison panel ${slot + 1}`);
    state.metadata.comparisons.forEach((entry, i) => {const o = new Option(entry.label || `${entry.condition1} vs ${entry.condition2}`, i); select.add(o);});
    select.value = view.slots[slot]; select.addEventListener('change', () => {view.slots[slot] = +select.value; renderAll(false);});
    card.append(select); grid.append(card);
    withView(context, () => {
      const styleDetails = document.createElement('details');
      styleDetails.innerHTML = '<summary>Sample line styles</summary>';
      const stylePanel = document.createElement('div'); stylePanel.className = 'sample-style-panel';
      styleDetails.append(stylePanel); card.append(styleDetails); renderSampleStyles(stylePanel);
      const plots = document.createElement('div'); plots.className = 'comparison-pair'; card.append(plots);
      for (const kind of ['rank', 'volcano']) {
        const box = document.createElement('div'), svg = document.createElementNS('http://www.w3.org/2000/svg','svg');
        svg.setAttribute('viewBox', kind === 'rank' ? '0 0 380 640' : '0 0 760 760');
        svg.classList.add(`side-${kind}`); box.append(svg); plots.append(box);
        if (kind === 'rank') {
          drawRank(svg); axisControls(box, 'rankX', plotControls.autoDomain(rankedVisible().map(p => plotControls.rankMetric(p, rankMode())), 'symmetric'));
        } else {
          renderVolcano(svg); axisControls(box,'volcanoX',plotControls.autoDomain(state.motifs.map(p=>p.effect),'symmetric'));
          axisControls(box,'volcanoY',plotControls.autoDomain(state.motifs.map(p=>p.neglog10p),'positive'));
        }
        registerPlot(`side-${slot}-${kind}`, {kind,node:svg,context});
      }
    });
    const aggregateBox = document.createElement('div'); aggregateBox.className = 'side-aggregate'; card.append(aggregateBox);
    if (!context.aggregate.has(activePrefix)) {
      aggregateBox.innerHTML = `<svg class="aggregate-panel" viewBox="0 0 300 300"><text x="150" y="130" text-anchor="middle" font-size="12">${esc(activePrefix || 'Selected motif')}</text><text x="150" y="155" text-anchor="middle" font-size="12">Aggregate profile unavailable</text></svg>`; continue;
    }
    let record;
    try { record = await withView(context, () => profileRecord(activePrefix)); }
    catch (error) { aggregateBox.textContent = `Aggregate profile unavailable: ${error.message}`; continue; }
    if (token !== view.generation) return;
    const motif = context.motifs.find(p => p.prefix === activePrefix);
    items.push({context,record,motif,index:slot,parent:aggregateBox});
  }
  if (token !== view.generation) return;
  prepareAggregateRanges(items);
  for (const item of items) withView(item.context, () => {
    item.parent.innerHTML = profileSvg(item.record, item.motif, item.index);
    item.node = item.parent.querySelector('svg');
    const legend = document.createElement('div'); legend.className = 'side-legend';
    legend.innerHTML = legendGroups().map(group => `<div><strong>${esc(group.condition)}</strong> · mean shown as a thick line</div>` + group.rows.map(row =>
      `<span style="display:inline-flex;align-items:center;margin-right:12px"><i style="width:22px;border-top:3px ${row.style.type === 'solid' ? 'solid' : row.style.type === 'dot' ? 'dotted' : 'dashed'} ${row.style.color};opacity:${row.style.alpha};margin-right:5px"></i>${esc(sampleDisplayName(row.sample, row.condition))}</span>`).join('')).join('');
    item.parent.prepend(legend); attachAggregateControls(item);
  });
  setFigureBusy(false);
}

function sidePanelSvg(kind) {
  const selector = kind === 'rank' ? '.side-rank' : kind === 'volcano' ? '.side-volcano' : kind === 'aggregate' ? '.aggregate-panel' : 'svg';
  const nodes = [...$('comparison-grid').querySelectorAll(selector)];
  const cellW = 760, cellH = 820, columns = Math.min(2, nodes.length) || 1;
  const rows = Math.ceil(nodes.length / columns) || 1;
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${columns*cellW} ${rows*cellH}">` + nodes.map((node,i)=>{
    const clone = styledClone(node);
    const item = [...view.rendered.values()].find(item => item.node === node && item.kind === 'aggregate');
    if (item) withView(item.context, () => {
      let y = 315;
      legendGroups().forEach(group => {
        clone.insertAdjacentHTML('beforeend', `<text x="10" y="${y}" font-size="10" font-weight="bold">${esc(group.condition)} · thick line: mean</text>`); y += 14;
        group.rows.forEach(row => {
          clone.insertAdjacentHTML('beforeend', `<line x1="10" x2="35" y1="${y-3}" y2="${y-3}" stroke="${row.style.color}" stroke-width="2"${dashAttribute(row.style.type)} stroke-opacity="${row.style.alpha}"/><text x="40" y="${y}" font-size="9">${esc(sampleDisplayName(row.sample,row.condition))}</text>`); y+=13;
        });
      });
      clone.setAttribute('viewBox', `0 0 300 ${y+5}`);
    });
    const box=clone.viewBox.baseVal;
    const scale=Math.min(cellW/box.width,(cellH-30)/box.height);
    const title=node.closest('.comparison-card')?.querySelector('select')?.selectedOptions[0]?.text || '';
    return `<g transform="translate(${i%columns*cellW},${Math.floor(i/columns)*cellH})"><text x="10" y="18" font-family="Arial" font-size="12">${esc(title)}</text><g transform="translate(0,30) scale(${scale})">${clone.innerHTML}</g></g>`;
  }).join('')+'</svg>';
}
