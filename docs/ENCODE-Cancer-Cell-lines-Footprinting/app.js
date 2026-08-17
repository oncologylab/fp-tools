"use strict";

const $ = (id) => document.getElementById(id);
const state = {
  metadata: null,
  entry: null,
  payload: null,
  motifs: [],
  aggregate: new Map(),
  profileAxis: [],
  first: "K562",
  second: "HepG2",
  selected: [],
  active: 0,
  sampleStyles: new Map(),
  payloadCache: new Map(),
  logoDataCache: new Map(),
  colors: { first: "#dc2626", second: "#2563eb", neutral: "#8a94a6" },
  request: 0,
  renderRequest: 0,
};
const plotSvgStyle =
  "svg,text{font-family:Helvetica,Arial,sans-serif}.plot-title{font-size:15px;font-weight:900;fill:#172033}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e3eaf3;stroke-width:1}.tick{font-size:11px;fill:#526176;font-weight:700}.axis-label{font-size:12px;fill:#243247;font-weight:900}.rank-bar.active{stroke:#111827;stroke-width:1.5}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.28))}";

function esc(value) {
  return String(value ?? "").replace(
    /[&<>"']/g,
    (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        c
      ],
  );
}
function finite(value, fallback = 0) {
  const number = Number(value);
  return Number.isFinite(number) ? number : fallback;
}
function fmt(value, digits = 3) {
  const number = Number(value);
  if (!Number.isFinite(number)) return "NA";
  const magnitude = Math.abs(number);
  if (number === 0) return "0";
  if (magnitude < 0.001 || magnitude >= 1000) return number.toExponential(2);
  return number.toFixed(digits).replace(/\.0+$|(?<=\.[0-9]*[1-9])0+$/g, "");
}
function fmtSci(value) {
  const number = Number(value);
  return Number.isFinite(number)
    ? number.toExponential(1).replace("e-0", "e-").replace("e+0", "e+")
    : "NA";
}
function motifLabel(motif) {
  return motif
    ? `${motif.name}${motif.motif_id ? ` (${motif.motif_id})` : ""}`
    : "";
}
function logp(value) {
  return -Math.log10(Math.max(1e-300, finite(value, 1)));
}
function niceStep(range, count = 5) {
  const raw = Math.max(range / count, 1e-12),
    power = 10 ** Math.floor(Math.log10(raw)),
    scaled = raw / power;
  return (scaled <= 1 ? 1 : scaled <= 2 ? 2 : scaled <= 5 ? 5 : 10) * power;
}
function niceLimit(value) {
  const number = Math.max(Math.abs(finite(value)), 1e-9),
    power = 10 ** Math.floor(Math.log10(number)),
    scaled = number / power;
  return (scaled <= 1 ? 1 : scaled <= 2 ? 2 : scaled <= 5 ? 5 : 10) * power;
}
function niceTicks(min, max, count = 5) {
  if (!Number.isFinite(min) || !Number.isFinite(max) || min === max)
    return [min || 0];
  const step = niceStep(max - min, count - 1),
    start = Math.ceil(min / step) * step,
    end = Math.floor(max / step) * step,
    out = [];
  for (let value = start; value <= end + step * 0.25; value += step)
    out.push(Math.abs(value) < step / 1e6 ? 0 : value);
  return out.length ? out : [min, max];
}
function colorFor(motif) {
  if (!motif.significant) return state.colors.neutral;
  return motif.effect >= 0 ? state.colors.first : state.colors.second;
}
function groupFor(motif) {
  if (!motif.significant) return "n.s.";
  return motif.effect >= 0 ? `${state.first}_up` : `${state.second}_up`;
}
function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob),
    anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = name;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
async function fetchJson(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`${path}: HTTP ${response.status}`);
  return response.json();
}
async function fetchGzipJson(path) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`${path}: HTTP ${response.status}`);
  if (!("DecompressionStream" in window))
    throw new Error(
      "This report needs a modern browser with gzip DecompressionStream support.",
    );
  const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
  return JSON.parse(await new Response(stream).text());
}

function fetchGzipJsonCached(path) {
  if (state.payloadCache.has(path)) return state.payloadCache.get(path);
  const request = fetchGzipJson(path).catch((error) => {
    state.payloadCache.delete(path);
    throw error;
  });
  state.payloadCache.set(path, request);
  return request;
}

function comparisonEntry(first, second) {
  const direct = state.metadata.comparisons.find(
    (item) => item.condition1 === first && item.condition2 === second,
  );
  if (direct) return { entry: direct, reversed: false };
  const reverse = state.metadata.comparisons.find(
    (item) => item.condition1 === second && item.condition2 === first,
  );
  if (reverse) return { entry: reverse, reversed: true };
  throw new Error(`No comparison for ${first} and ${second}`);
}

function orientedMotif(point, reversed) {
  const aggregate = state.aggregate.get(point.prefix),
    effect = (reversed ? -1 : 1) * finite(point.change),
    significant = point.group !== "n.s.",
    group = significant
      ? `${effect >= 0 ? state.first : state.second}_up`
      : "n.s.";
  return {
    ...point,
    effect,
    pvalue: finite(point.pvalue, 1),
    neglog10p: finite(point.neglog10p, logp(point.pvalue)),
    qvalue: finite(point.fdr, 1),
    significant,
    group,
    n_sites: aggregate?.n_sites ?? "",
  };
}

function conditionSamples(condition) {
  const motif = state.aggregate.values().next().value,
    record = motif?.conditions?.find((item) => item.name === condition);
  return (
    record?.samples?.map((item) => item.name) ||
    state.metadata.conditions.find((item) => item.name === condition)
      ?.samples ||
    []
  );
}
function plotCount() {
  return Math.max(1, Math.min(12, Number($("plot-count").value) || 4));
}
function sortedMotifs() {
  return state.motifs.slice().sort((a, b) =>
    motifLabel(a).localeCompare(motifLabel(b), undefined, {
      sensitivity: "base",
    }),
  );
}

function defaultSelected(target) {
  const withProfiles = state.motifs.filter((item) =>
      state.aggregate.has(item.prefix),
    ),
    positive = withProfiles
      .filter((item) => item.effect > 0)
      .sort((a, b) => b.effect - a.effect || a.pvalue - b.pvalue);
  const negative = withProfiles
    .filter((item) => item.effect < 0)
    .sort((a, b) => a.effect - b.effect || a.pvalue - b.pvalue);
  const output = [],
    negativeCount = Math.floor(target / 2),
    positiveCount = target - negativeCount;
  [
    ...positive.slice(0, positiveCount),
    ...negative.slice(0, negativeCount),
  ].forEach((item) => {
    if (!output.includes(item.prefix)) output.push(item.prefix);
  });
  withProfiles
    .slice()
    .sort(
      (a, b) => Math.abs(b.effect) - Math.abs(a.effect) || a.pvalue - b.pvalue,
    )
    .forEach((item) => {
      if (output.length < target && !output.includes(item.prefix))
        output.push(item.prefix);
    });
  return output.slice(0, target);
}

function ensureSelected(reset = false) {
  const target = plotCount(),
    valid = new Set(state.motifs.map((item) => item.prefix));
  if (reset) state.selected = [];
  state.selected = state.selected
    .filter((prefix) => valid.has(prefix))
    .slice(0, target);
  defaultSelected(target).forEach((prefix) => {
    if (state.selected.length < target && !state.selected.includes(prefix))
      state.selected.push(prefix);
  });
  state.active = Math.max(0, Math.min(state.active, target - 1));
}

function updateHeader() {
  $("title-cond1").textContent = state.first;
  $("title-cond2").textContent = state.second;
  $("title-cond1").style.color = state.colors.first;
  $("title-cond2").style.color = state.colors.second;
  document.title = `Differential footprint report (${state.first} vs ${state.second})`;
  $("report-method").textContent = state.payload?.report_label || "";
}

function renderColorControls() {
  const rows = [
    { key: "first", label: `${state.first}_up` },
    { key: "second", label: `${state.second}_up` },
    { key: "neutral", label: "n.s." },
  ];
  $("color-controls").innerHTML = rows
    .map(
      (row) =>
        `<label class="color-row"><span title="${esc(row.label)}">${esc(row.label)}</span><input type="color" data-group-color="${row.key}" value="${state.colors[row.key]}"></label>`,
    )
    .join("");
  $("color-controls")
    .querySelectorAll("[data-group-color]")
    .forEach((input) =>
      input.addEventListener("input", () => {
        state.colors[input.dataset.groupColor] = input.value;
        renderAll(false);
      }),
    );
}

function defaultSampleStyle(sample, condition, index) {
  const conditionColor =
    condition === state.first ? state.colors.first : state.colors.second;
  return {
    visible: true,
    color: conditionColor,
    alpha: 0.9,
    width: 2,
    type: "solid",
  };
}

function sampleStyle(sample, condition, index) {
  if (!state.sampleStyles.has(sample))
    state.sampleStyles.set(
      sample,
      defaultSampleStyle(sample, condition, index),
    );
  return state.sampleStyles.get(sample);
}

function renderSampleStyles() {
  const conditions = [state.first, state.second];
  $("sample-style-panel").innerHTML = conditions
    .map((condition) => {
      const samples = conditionSamples(condition),
        conditionColor =
          condition === state.first ? state.colors.first : state.colors.second;
      const rows = samples
        .map((sample, index) => {
          const style = sampleStyle(sample, condition, index);
          return `<div class="sample-style-row"><input type="checkbox" data-sample-visible="${esc(sample)}" ${style.visible ? "checked" : ""} aria-label="Show ${esc(sample)}"><span class="sample-style-name" title="${esc(sample)}">${esc(sample)}</span><input type="color" data-sample-color="${esc(sample)}" value="${style.color}" aria-label="Color for ${esc(sample)}"><input type="number" data-sample-alpha="${esc(sample)}" min="0.1" max="1" step="0.1" value="${style.alpha}" aria-label="Opacity for ${esc(sample)}"><input type="number" data-sample-width="${esc(sample)}" min="0.3" max="4" step="0.1" value="${style.width}" aria-label="Width for ${esc(sample)}"><select data-sample-type="${esc(sample)}" aria-label="Line type for ${esc(sample)}"><option value="solid" ${style.type === "solid" ? "selected" : ""}>Solid</option><option value="dash" ${style.type === "dash" ? "selected" : ""}>Dash</option><option value="dot" ${style.type === "dot" ? "selected" : ""}>Dot</option></select></div>`;
        })
        .join("");
      return `<div class="sample-style-group"><div class="sample-style-group-title"><i class="sample-style-dot" style="background:${conditionColor}"></i>${esc(condition)}</div><div class="sample-style-row sample-style-head"><span>Show</span><span>Sample</span><span>Color</span><span>Alpha</span><span>Width</span><span>Type</span></div>${rows}</div>`;
    })
    .join("");
  const panel = $("sample-style-panel");
  panel.querySelectorAll("[data-sample-visible]").forEach((input) =>
    input.addEventListener("change", () => {
      state.sampleStyles.get(input.dataset.sampleVisible).visible =
        input.checked;
      renderAll(false);
    }),
  );
  panel.querySelectorAll("[data-sample-color]").forEach((input) =>
    input.addEventListener("input", () => {
      state.sampleStyles.get(input.dataset.sampleColor).color = input.value;
      renderAll(false);
    }),
  );
  panel.querySelectorAll("[data-sample-alpha]").forEach((input) =>
    input.addEventListener("input", () => {
      state.sampleStyles.get(input.dataset.sampleAlpha).alpha = Math.max(
        0.1,
        Math.min(1, finite(input.value, 0.7)),
      );
      renderAll(false);
    }),
  );
  panel.querySelectorAll("[data-sample-width]").forEach((input) =>
    input.addEventListener("input", () => {
      state.sampleStyles.get(input.dataset.sampleWidth).width = Math.max(
        0.3,
        Math.min(4, finite(input.value, 2)),
      );
      renderAll(false);
    }),
  );
  panel.querySelectorAll("[data-sample-type]").forEach((input) =>
    input.addEventListener("change", () => {
      state.sampleStyles.get(input.dataset.sampleType).type = input.value;
      renderAll(false);
    }),
  );
}

function logoPath(prefix) {
  const base = state.metadata?.logos?.base || "data/logos";
  return `${base}/${encodeURIComponent(prefix)}.png`;
}

function motifLogoHtml(prefix) {
  return `<img alt="${esc(prefix)} motif logo" src="${esc(logoPath(prefix))}" width="1000" height="250" style="display:block;max-width:100%;max-height:60px;width:auto;height:auto;object-fit:contain">`;
}

function logoDataUri(prefix) {
  if (state.logoDataCache.has(prefix)) return state.logoDataCache.get(prefix);
  const request = fetch(logoPath(prefix))
    .then((response) => {
      if (!response.ok)
        throw new Error(`Motif logo ${prefix}: HTTP ${response.status}`);
      return response.blob();
    })
    .then(
      (blob) =>
        new Promise((resolve, reject) => {
          const reader = new FileReader();
          reader.onload = () => resolve(reader.result);
          reader.onerror = () => reject(reader.error);
          reader.readAsDataURL(blob);
        }),
    )
    .catch((error) => {
      state.logoDataCache.delete(prefix);
      throw error;
    });
  state.logoDataCache.set(prefix, request);
  return request;
}

function renderSelectedCards() {
  const motifs = sortedMotifs(),
    byPrefix = new Map(state.motifs.map((item) => [item.prefix, item]));
  $("selected-grid").innerHTML = state.selected
    .slice(0, plotCount())
    .map((prefix, index) => {
      const motif = byPrefix.get(prefix),
        options = motifs
          .map(
            (item) =>
              `<option value="${esc(item.prefix)}" ${item.prefix === prefix ? "selected" : ""}>${esc(motifLabel(item))}</option>`,
          )
          .join(""),
        group = groupFor(motif);
      return `<article class="selected-motif${index === state.active ? " active" : ""}" data-selected-panel="${index}"><div class="selected-head"><select class="panel-tf" data-panel-tf="${index}" aria-label="Motif for aggregate plot ${index + 1}">${options}</select></div><div class="motif-logo">${motifLogoHtml(prefix)}</div><div class="detail-grid"><p class="motif-group" style="color:${colorFor(motif)}">${esc(group)}</p><p class="metric-line">ΔFP = ${fmt(motif.effect, 4)}</p><p class="metric-line">FDR = ${fmtSci(motif.qvalue)}</p></div></article>`;
    })
    .join("");
  $("selected-grid")
    .querySelectorAll("[data-selected-panel]")
    .forEach((card) =>
      card.addEventListener("click", (event) => {
        if (event.target.closest("select")) return;
        state.active = Number(card.dataset.selectedPanel);
        renderAll(false);
      }),
    );
  $("selected-grid")
    .querySelectorAll("[data-panel-tf]")
    .forEach((select) =>
      select.addEventListener("change", () => {
        state.selected[Number(select.dataset.panelTf)] = select.value;
        state.active = Number(select.dataset.panelTf);
        renderAll(false);
      }),
    );
}

function visibleSelected() {
  return new Set(state.selected.slice(0, plotCount()));
}

function renderVolcano() {
  const width = 760,
    height = 760,
    margin = { top: 34, right: 48, bottom: 60, left: 84 },
    innerWidth = width - margin.left - margin.right,
    innerHeight = height - margin.top - margin.bottom,
    xValues = state.motifs.map((item) => item.effect),
    yValues = state.motifs.map((item) => item.neglog10p),
    xLimit = niceLimit(
      Math.max(...xValues.map((value) => Math.abs(value)), 1e-9) * 1.05,
    ),
    rawYMax = Math.max(...yValues, 0),
    yMax = rawYMax > 0 ? rawYMax * 1.05 : 1,
    sx = (value) =>
      margin.left + ((value + xLimit) / (2 * xLimit)) * innerWidth,
    sy = (value) => margin.top + innerHeight - (value / yMax) * innerHeight,
    tickStyle =
      "font-size:15px;font-weight:900;font-family:Helvetica,Arial,sans-serif",
    axisStyle =
      "font-size:17px;font-weight:900;font-family:Helvetica,Arial,sans-serif",
    selected = visibleSelected(),
    parts = [
      `<style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><rect x="${margin.left}" y="${margin.top}" width="${innerWidth}" height="${innerHeight}" fill="#fbfdff" stroke="#d9e2ec"/>`,
    ];
  niceTicks(0, yMax, 7).forEach((value) =>
    parts.push(
      `<line x1="${margin.left}" y1="${sy(value)}" x2="${margin.left + innerWidth}" y2="${sy(value)}" class="grid"/><text x="${margin.left - 12}" y="${sy(value) + 5}" class="tick" style="${tickStyle}" text-anchor="end">${fmt(value, 1)}</text>`,
    ),
  );
  niceTicks(-xLimit, xLimit, 7).forEach((value) =>
    parts.push(
      `<line x1="${sx(value)}" y1="${margin.top}" x2="${sx(value)}" y2="${margin.top + innerHeight}" class="grid"/><text x="${sx(value)}" y="${margin.top + innerHeight + 24}" class="tick" style="${tickStyle}" text-anchor="middle">${fmt(value, 3)}</text>`,
    ),
  );
  parts.push(
    `<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top + innerHeight}" class="zero"/><line x1="${margin.left}" y1="${margin.top + innerHeight}" x2="${margin.left + innerWidth}" y2="${margin.top + innerHeight}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + innerHeight}" class="axis"/><text x="${margin.left + innerWidth / 2}" y="${height - 10}" class="axis-label" style="${axisStyle}" text-anchor="middle">${esc(state.payload?.change_label || "Differential footprint score")}</text><text x="16" y="${margin.top + innerHeight / 2}" class="axis-label" style="${axisStyle}" text-anchor="middle" transform="rotate(-90 16 ${margin.top + innerHeight / 2})">-log10(p-value)</text><text x="${margin.left + 18}" y="${margin.top + innerHeight - 14}" font-size="24" font-weight="900" fill="${state.colors.second}">${esc(state.second)}_up</text><text x="${margin.left + innerWidth - 18}" y="${margin.top + innerHeight - 14}" text-anchor="end" font-size="24" font-weight="900" fill="${state.colors.first}">${esc(state.first)}_up</text>`,
  );
  state.motifs
    .map((item) => ({ item, selected: selected.has(item.prefix) }))
    .sort((a, b) => Number(a.selected) - Number(b.selected))
    .forEach(({ item, selected: isSelected }) =>
      parts.push(
        `<circle class="pt${isSelected ? " selected" : ""}" data-prefix="${esc(item.prefix)}" cx="${sx(item.effect).toFixed(2)}" cy="${sy(item.neglog10p).toFixed(2)}" r="${isSelected ? 7.2 : 4.2}" fill="${colorFor(item)}" fill-opacity="${isSelected ? 0.98 : 0.76}" stroke="${isSelected ? "#111827" : "#fff"}" stroke-width="${isSelected ? 2.7 : 0.9}"><title>${esc(motifLabel(item))}: ΔFP ${fmt(item.effect, 4)}, FDR ${fmtSci(item.qvalue)}</title></circle>`,
      ),
    );
  $("chart").innerHTML = parts.join("");
  $("chart")
    .querySelectorAll("[data-prefix]")
    .forEach((node) =>
      node.addEventListener("click", () =>
        setSelectedMotif(node.dataset.prefix),
      ),
    );
}

function drawRank() {
  const limit = Math.max(
      2,
      Math.min(200, Math.floor(Number($("rank-rows").value) || 20)),
    ),
    perDirection = Math.max(1, Math.floor(limit / 2)),
    positive = state.motifs
      .filter((item) => item.effect > 0)
      .sort((a, b) => b.effect - a.effect || a.pvalue - b.pvalue)
      .slice(0, perDirection),
    negative = state.motifs
      .filter((item) => item.effect < 0)
      .sort((a, b) => a.effect - b.effect || a.pvalue - b.pvalue)
      .slice(0, perDirection),
    shown = [...negative, ...positive],
    width = 380,
    rowHeight = 14,
    rowGap = 3,
    sectionGap = 8,
    margin = { top: 64, bottom: 68, left: 128, right: 14 },
    height = Math.max(
      430,
      margin.top +
        shown.length * (rowHeight + rowGap) +
        sectionGap +
        margin.bottom,
    ),
    xMiddle = 246,
    xWidth = 112,
    maxAbs = niceLimit(
      Math.max(...shown.map((item) => Math.abs(item.effect)), 1e-9),
    ),
    sx = (value) => xMiddle + (value / maxAbs) * xWidth,
    axisY = height - 60,
    selected = visibleSelected(),
    parts = [
      `<style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><text x="${width / 2}" y="18" class="plot-title" text-anchor="middle">Top differential motifs</text><line x1="${xMiddle}" y1="${margin.top - 36}" x2="${xMiddle}" y2="${axisY}" stroke="#172033" stroke-width="2.2"/><text x="${xMiddle - 6}" y="${margin.top - 22}" text-anchor="end" font-size="14" font-weight="900" fill="${state.colors.second}">${esc(state.second)}_up</text><text x="${xMiddle + 6}" y="${margin.top - 22}" text-anchor="start" font-size="14" font-weight="900" fill="${state.colors.first}">${esc(state.first)}_up</text>`,
    ];
  niceTicks(-maxAbs, maxAbs, 5).forEach((value) =>
    parts.push(
      `<line x1="${sx(value)}" y1="${axisY - 4}" x2="${sx(value)}" y2="${axisY + 4}" class="axis"/><text x="${sx(value)}" y="${axisY + 17}" class="tick" text-anchor="middle">${fmt(value, 3)}</text>`,
    ),
  );
  parts.push(
    `<line x1="${sx(-maxAbs)}" y1="${axisY}" x2="${sx(maxAbs)}" y2="${axisY}" class="axis"/><text x="${xMiddle}" y="${height - 8}" class="axis-label" text-anchor="middle">${esc(state.payload?.change_label || "Differential footprint score")}</text>`,
  );
  let y = margin.top;
  const drawRows = (rows) =>
    rows.forEach((item) => {
      const barWidth = (Math.abs(item.effect) / maxAbs) * xWidth,
        x = item.effect >= 0 ? xMiddle : xMiddle - barWidth,
        isSelected = selected.has(item.prefix),
        name = motifLabel(item).slice(0, 20),
        labelY = y + rowHeight - 2;
      parts.push(
        `<text class="rank-name${isSelected ? " active" : ""}" data-prefix="${esc(item.prefix)}" x="6" y="${labelY}" font-size="10" font-weight="${isSelected ? 900 : 700}" fill="${isSelected ? colorFor(item) : "#526176"}">${esc(name)}</text><rect class="rank-bar${isSelected ? " active" : ""}" data-prefix="${esc(item.prefix)}" x="${x}" y="${y}" width="${Math.max(1, barWidth)}" height="${rowHeight}" fill="${colorFor(item)}" fill-opacity="${isSelected ? 0.95 : 0.72}"><title>${esc(motifLabel(item))}: ${fmt(item.effect, 4)}</title></rect><text x="${item.effect >= 0 ? x - 3 : x + barWidth + 3}" y="${labelY}" class="tick" text-anchor="${item.effect >= 0 ? "end" : "start"}">${fmt(item.effect, 3)}</text>`,
      );
      y += rowHeight + rowGap;
    });
  drawRows(negative);
  y += sectionGap;
  drawRows(positive);
  $("rank-chart").setAttribute("viewBox", `0 0 ${width} ${height}`);
  $("rank-chart").innerHTML = parts.join("");
  $("rank-chart")
    .querySelectorAll("[data-prefix]")
    .forEach((node) =>
      node.addEventListener("click", () =>
        setSelectedMotif(node.dataset.prefix),
      ),
    );
}

async function profileRecord(prefix) {
  let motif = state.aggregate.get(prefix);
  if (!motif) throw new Error(`No aggregate profile for ${prefix}`);
  const hasProfiles = motif.conditions?.some((condition) =>
    condition.samples?.some((sample) => Array.isArray(sample.profile)),
  );
  if (!hasProfiles) {
    const shard = state.entry?.profile_shards?.find(
      (item) => Number(item.id) === Number(motif.profile_shard),
    );
    if (!shard) throw new Error(`No profile shard for ${prefix}`);
    const shardPayload = await fetchGzipJsonCached(shard.file);
    motif = shardPayload.motifs.find((item) => item.prefix === prefix);
    if (!motif) throw new Error(`Profile shard does not contain ${prefix}`);
  }
  const samples = {},
    sampleMeta = {};
  motif.conditions.forEach((condition) =>
    condition.samples.forEach((sample) => {
      samples[sample.name] = sample.profile;
      sampleMeta[sample.name] = sample;
    }),
  );
  return { samples, sampleMeta, n_profile_sites: motif.n_sites || 0 };
}

function dashAttribute(type) {
  return type === "dash"
    ? ' stroke-dasharray="7 4"'
    : type === "dot"
      ? ' stroke-dasharray="2 3"'
      : "";
}
function linePath(profile, axis, sx, sy) {
  return profile
    .map(
      (value, index) =>
        `${index ? "L" : "M"}${sx(axis[index]).toFixed(2)},${sy(value).toFixed(2)}`,
    )
    .join(" ");
}

function profileSvg(record, motif, index) {
  const rawAxis = state.profileAxis,
    keep = rawAxis
      .map((value, i) => ({ value, i }))
      .filter((item) => item.value >= -60 && item.value <= 60),
    axis = keep.map((item) => item.value),
    series = [];
  [
    [state.first, conditionSamples(state.first)],
    [state.second, conditionSamples(state.second)],
  ].forEach(([condition, samples]) =>
    samples.forEach((sample, sampleIndex) => {
      const style = sampleStyle(sample, condition, sampleIndex),
        raw = record.samples[sample] || [];
      if (style.visible)
        series.push({
          sample,
          condition,
          style,
          fpScore: finite(record.sampleMeta?.[sample]?.fp_score),
          profile: keep.map((item) => finite(raw[item.i])),
        });
    }),
  );
  const values = series.flatMap((item) => item.profile),
    rawMin = Math.min(...values, 0),
    rawMax = Math.max(...values, 1e-9),
    padding = Math.max((rawMax - rawMin || 1) * 0.18, 1e-6),
    yMin = rawMin - padding,
    yMax = rawMax + padding,
    width = 300,
    height = 300,
    margin = { top: 30, right: 8, bottom: 34, left: 36 },
    innerWidth = width - margin.left - margin.right,
    innerHeight = height - margin.top - margin.bottom,
    sx = (value) =>
      margin.left +
      ((value - axis[0]) / (axis[axis.length - 1] - axis[0] || 1)) * innerWidth,
    sy = (value) =>
      margin.top +
      innerHeight -
      ((value - yMin) / (yMax - yMin || 1)) * innerHeight,
    parts = [
      `<svg class="aggregate-panel" data-panel="${index}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg"><style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><text x="${width / 2}" y="18" class="plot-title" text-anchor="middle">${esc(motifLabel(motif))}</text>`,
    ];
  niceTicks(yMin, yMax, 4).forEach((value) =>
    parts.push(
      `<line x1="${margin.left}" y1="${sy(value)}" x2="${margin.left + innerWidth}" y2="${sy(value)}" class="grid"/><text x="${margin.left - 6}" y="${sy(value) + 3}" class="tick" text-anchor="end">${fmt(value, 3)}</text>`,
    ),
  );
  [-60, 0, 60].forEach((value) =>
    parts.push(
      `<line x1="${sx(value)}" y1="${margin.top}" x2="${sx(value)}" y2="${margin.top + innerHeight}" class="${value === 0 ? "zero" : "grid"}"/><text x="${sx(value)}" y="${margin.top + innerHeight + 17}" class="tick" text-anchor="middle">${value}</text>`,
    ),
  );
  parts.push(
    `<line x1="${margin.left}" y1="${margin.top + innerHeight}" x2="${margin.left + innerWidth}" y2="${margin.top + innerHeight}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + innerHeight}" class="axis"/><text x="${margin.left + 7}" y="${margin.top + innerHeight - 8}" class="tick" fill="#94a3b8">${record.n_profile_sites.toLocaleString()}</text>`,
  );
  series
    .sort((a, b) => b.fpScore - a.fpScore)
    .forEach((item) =>
      parts.push(
        `<path d="${linePath(item.profile, axis, sx, sy)}" fill="none" stroke="${item.style.color}" stroke-width="${item.style.width}"${dashAttribute(item.style.type)} stroke-opacity="${item.style.alpha}"><title>${esc(item.sample)}</title></path>`,
      ),
    );
  parts.push(
    `<text x="${margin.left + innerWidth / 2}" y="${height - 6}" class="axis-label" text-anchor="middle">${esc(state.payload?.aggregate?.x_label || "Distance from motif center (bp)")}</text><text x="10" y="${margin.top + innerHeight / 2}" class="axis-label" text-anchor="middle" transform="rotate(-90 10 ${margin.top + innerHeight / 2})">${esc(state.payload?.aggregate?.y_label || "Corrected cut-site signal")}</text></svg>`,
  );
  return parts.join("");
}

function aggregateShape(count) {
  if (count <= 1) return { columns: 1, rows: 1 };
  if (count <= 2) return { columns: 2, rows: 1 };
  if (count <= 4) return { columns: 2, rows: 2 };
  if (count <= 6) return { columns: 3, rows: 2 };
  if (count <= 9) return { columns: 3, rows: 3 };
  return { columns: 4, rows: 3 };
}

function renderLegend() {
  const rows = [];
  [
    [state.first, conditionSamples(state.first)],
    [state.second, conditionSamples(state.second)],
  ].forEach(([condition, samples]) =>
    samples.forEach((sample, index) => {
      const style = sampleStyle(sample, condition, index);
      if (style.visible) rows.push({ sample, style });
    }),
  );
  $("aggregate-legend").innerHTML = rows
    .map(
      (row) =>
        `<div class="legend-row"><i class="legend-line" style="border-top-color:${row.style.color};border-top-width:${Math.max(2,row.style.width)}px;border-top-style:${row.style.type === "dash" ? "dashed" : row.style.type === "dot" ? "dotted" : "solid"};opacity:${row.style.alpha}"></i><span title="${esc(row.sample)}">${esc(row.sample)}</span></div>`,
    )
    .join("");
}

async function renderAggregateGrid() {
  const token = ++state.renderRequest,
    prefixes = state.selected.slice(0, plotCount()),
    shape = aggregateShape(prefixes.length),
    grid = $("aggregate-grid");
  grid.style.setProperty("--aggregate-cols", shape.columns);
  grid.style.setProperty("--aggregate-rows", shape.rows);
  grid.innerHTML = prefixes
    .map(
      (prefix, index) =>
        `<div class="aggregate-tile${index === state.active ? " active" : ""}" data-tile="${index}"><svg viewBox="0 0 300 300"><text x="150" y="150" text-anchor="middle" class="axis-label">Loading profile…</text></svg></div>`,
    )
    .join("");
  try {
    const records = await Promise.all(prefixes.map(profileRecord));
    if (token !== state.renderRequest) return;
    const byPrefix = new Map(state.motifs.map((item) => [item.prefix, item]));
    grid.innerHTML = prefixes
      .map(
        (prefix, index) =>
          `<div class="aggregate-tile${index === state.active ? " active" : ""}" data-tile="${index}">${profileSvg(records[index], byPrefix.get(prefix), index)}</div>`,
      )
      .join("");
    grid.querySelectorAll("[data-tile]").forEach((tile) =>
      tile.addEventListener("click", () => {
        state.active = Number(tile.dataset.tile);
        renderAll(false);
      }),
    );
  } catch (error) {
    if (token !== state.renderRequest) return;
    grid.innerHTML = `<div class="aggregate-tile"><svg viewBox="0 0 300 300"><text x="150" y="145" text-anchor="middle" class="axis-label">Profile unavailable</text><text x="150" y="165" text-anchor="middle" class="tick">${esc(error.message)}</text></svg></div>`;
  }
}

function setSelectedMotif(prefix) {
  state.selected[state.active] = prefix;
  renderAll(false);
}

function renderAll(refreshControls = true) {
  ensureSelected();
  updateHeader();
  if (refreshControls) {
    renderColorControls();
    renderSampleStyles();
  }
  renderSelectedCards();
  drawRank();
  renderVolcano();
  renderLegend();
  renderAggregateGrid();
}

async function loadComparison(reset = true) {
  const first = $("condition-1").value,
    second = $("condition-2").value;
  if (first === second) return;
  state.first = first;
  state.second = second;
  const token = ++state.request;
  $("status").textContent = `Loading ${first} vs ${second}…`;
  const { entry } = comparisonEntry(first, second),
    payload = await fetchGzipJsonCached(entry.file);
  if (token !== state.request) return;
  state.entry = entry;
  state.payload = payload;
  state.aggregate = new Map(
    (payload.aggregate?.motifs || []).map((item) => [item.prefix, item]),
  );
  state.profileAxis = payload.aggregate?.x || [];
  const reversed = payload.conditions[0] !== first;
  if (
    new Set(payload.conditions).size !== 2 ||
    !payload.conditions.includes(first) ||
    !payload.conditions.includes(second)
  )
    throw new Error(`Payload conditions do not match ${first} and ${second}`);
  state.motifs = payload.points.map((item) => orientedMotif(item, reversed));
  state.colors = {
    first: payload.colors[`${first}_up`] || "#dc2626",
    second: payload.colors[`${second}_up`] || "#2563eb",
    neutral: payload.colors["n.s."] || "#8a94a6",
  };
  state.sampleStyles = new Map();
  state.active = 0;
  ensureSelected(reset);
  renderAll(true);
  const significant = state.motifs.filter((item) => item.significant).length;
  $("status").textContent =
    `${state.motifs.length.toLocaleString()} motifs | ${significant.toLocaleString()} significant | ${first} minus ${second}`;
}

function handleConditionChange(changed) {
  const oldFirst = state.first,
    oldSecond = state.second;
  let first = $("condition-1").value,
    second = $("condition-2").value;
  if (first === second) {
    if (changed === "first") second = oldFirst;
    else first = oldSecond;
    $("condition-1").value = first;
    $("condition-2").value = second;
  }
  loadComparison(true).catch(showError);
}

function comparisonTsv() {
  const columns = [
      "condition1",
      "condition2",
      "prefix",
      "name",
      "motif_id",
      "group",
      "n_sites",
      "effect",
      "pvalue",
      "qvalue",
      "significant",
    ],
    rows = [columns.join("\t")];
  state.motifs.forEach((motif) =>
    rows.push(
      [
        state.first,
        state.second,
        ...columns.slice(2).map((key) => motif[key] ?? ""),
      ].join("\t"),
    ),
  );
  return rows.join("\n") + "\n";
}

function styledClone(node) {
  const clone = node.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("font-family", "Helvetica,Arial,sans-serif");
  if (!clone.querySelector("style"))
    clone.insertAdjacentHTML("afterbegin", `<style>${plotSvgStyle}</style>`);
  return clone;
}
function serializeSvg(node) {
  return new XMLSerializer().serializeToString(styledClone(node));
}

async function logoPanelSvg() {
  const cards = [...$("selected-grid").querySelectorAll(".selected-motif")],
    cardWidth = 240,
    cardHeight = 220,
    gap = 10,
    columns = Math.max(1, Math.min(4, cards.length)),
    rows = Math.ceil(cards.length / columns),
    byPrefix = new Map(state.motifs.map((item) => [item.prefix, item])),
    prefixes = state.selected.slice(0, plotCount()),
    images = await Promise.all(prefixes.map(logoDataUri)),
    parts = [
      `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${columns * cardWidth + (columns - 1) * gap} ${rows * cardHeight + (rows - 1) * gap}" font-family="Helvetica,Arial,sans-serif"><style>${plotSvgStyle}.motif-card-title{font-size:13px;font-weight:900;fill:#172033}.motif-card-metric{font-size:12px;font-weight:800;fill:#334155}</style><rect width="100%" height="100%" fill="#fff"/>`,
    ];
  prefixes.forEach((prefix, index) => {
    const motif = byPrefix.get(prefix),
      x = (index % columns) * (cardWidth + gap),
      y = Math.floor(index / columns) * (cardHeight + gap);
    parts.push(
      `<g transform="translate(${x},${y})"><rect width="${cardWidth}" height="${cardHeight}" rx="7" fill="#fff" stroke="${index === state.active ? "#93c5fd" : "#d8e2ef"}" stroke-width="${index === state.active ? 3 : 1}"/><text x="10" y="20" class="motif-card-title">${esc(motifLabel(motif)).slice(0, 32)}</text><image x="24" y="34" width="${cardWidth - 48}" height="96" preserveAspectRatio="xMidYMid meet" href="${images[index]}"/><text x="10" y="150" font-size="13" font-weight="900" fill="${colorFor(motif)}">${esc(groupFor(motif))}</text><text x="10" y="171" class="motif-card-metric">ΔFP = ${fmt(motif.effect, 4)}</text><text x="10" y="192" class="motif-card-metric">FDR = ${fmtSci(motif.qvalue)}</text></g>`,
    );
  });
  parts.push("</svg>");
  return parts.join("");
}

function aggregateGridSvg() {
  const svgs = [...document.querySelectorAll(".aggregate-panel")],
    shape = aggregateShape(svgs.length),
    plotWidth = 300,
    plotHeight = 300,
    gridWidth = shape.columns * plotWidth,
    gridHeight = shape.rows * plotHeight,
    legendRows = [...$("aggregate-legend").querySelectorAll(".legend-row")],
    legendWidth = legendRows.length ? 170 : 0,
    totalWidth = gridWidth + (legendWidth ? 12 + legendWidth : 0),
    parts = [
      `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalWidth} ${gridHeight}" font-family="Helvetica,Arial,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`,
    ];
  svgs.forEach((svg, index) => {
    const clone = styledClone(svg);
    clone.querySelector("style")?.remove();
    parts.push(
      `<g transform="translate(${(index % shape.columns) * plotWidth},${Math.floor(index / shape.columns) * plotHeight})">${clone.innerHTML}</g>`,
    );
  });
  if (legendRows.length) {
    const x = gridWidth + 12,
      height = legendRows.length * 16 + 14;
    parts.push(
      `<g transform="translate(${x},8)"><rect width="${legendWidth}" height="${height}" rx="5" fill="#fff" stroke="#d8e2ef"/>`,
    );
    legendRows.forEach((row, index) => {
      const line = row.querySelector(".legend-line"),
        label = row.querySelector("span").textContent,
        y = 14 + index * 16,
        style = line.style,
        dash =
          style.borderTopStyle === "dashed"
            ? ' stroke-dasharray="7 4"'
            : style.borderTopStyle === "dotted"
              ? ' stroke-dasharray="2 3"'
              : "";
      parts.push(
        `<line x1="8" y1="${y - 3}" x2="38" y2="${y - 3}" stroke="${style.borderTopColor}" stroke-width="${Math.max(2, parseFloat(style.borderTopWidth) || 2)}"${dash} stroke-opacity="${style.opacity || 1}"/><text x="44" y="${y}" class="tick">${esc(label)}</text>`,
      );
    });
    parts.push("</g>");
  }
  parts.push("</svg>");
  return parts.join("");
}

function combinedPanelSvg() {
  const rank = styledClone($("rank-chart")),
    volcano = styledClone($("chart")),
    aggregate = aggregateGridSvg(),
    aggregateDocument = new DOMParser().parseFromString(
      aggregate,
      "image/svg+xml",
    ).documentElement,
    rankBox = rank.viewBox.baseVal,
    rankWidth = rankBox.width || 380,
    rankHeight = rankBox.height || 600,
    volcanoSize = 760,
    panelHeight = Math.max(760, rankHeight),
    rankScale = panelHeight / rankHeight,
    rankDisplayWidth = rankWidth * rankScale,
    aggregateBox = aggregateDocument.viewBox.baseVal,
    aggregateScale = panelHeight / (aggregateBox.height || 600),
    aggregateWidth = (aggregateBox.width || 600) * aggregateScale,
    gap = 20,
    totalWidth = rankDisplayWidth + volcanoSize + aggregateWidth + gap * 2;
  rank.querySelector("style")?.remove();
  volcano.querySelector("style")?.remove();
  aggregateDocument.querySelector("style")?.remove();
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalWidth} ${panelHeight}" font-family="Helvetica,Arial,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/><g transform="scale(${rankScale})">${rank.innerHTML}</g><g transform="translate(${rankDisplayWidth + gap},0)">${volcano.innerHTML}</g><g transform="translate(${rankDisplayWidth + volcanoSize + gap * 2},0) scale(${aggregateScale})">${aggregateDocument.innerHTML}</g></svg>`;
}

function svgDimensions(svg) {
  const doc = new DOMParser().parseFromString(svg, "image/svg+xml"),
    viewBox = doc.documentElement
      .getAttribute("viewBox")
      ?.split(/\s+/)
      .map(Number) || [0, 0, 1200, 800];
  return { width: viewBox[2] || 1200, height: viewBox[3] || 800 };
}
function exportSvg(svg, name) {
  const format = $("figure-format").value;
  if (format === "svg") {
    downloadBlob(
      new Blob([svg], { type: "image/svg+xml;charset=utf-8" }),
      `${name}.svg`,
    );
    return;
  }
  if (format === "png") {
    const dimensions = svgDimensions(svg),
      image = new Image(),
      url = URL.createObjectURL(
        new Blob([svg], { type: "image/svg+xml;charset=utf-8" }),
      );
    image.onload = () => {
      const scale = Math.min(3, Math.max(1, 2400 / dimensions.width)),
        canvas = document.createElement("canvas");
      canvas.width = Math.round(dimensions.width * scale);
      canvas.height = Math.round(dimensions.height * scale);
      const context = canvas.getContext("2d");
      context.fillStyle = "#fff";
      context.fillRect(0, 0, canvas.width, canvas.height);
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob((blob) => {
        if (blob) downloadBlob(blob, `${name}.png`);
      }, "image/png");
    };
    image.onerror = () => {
      URL.revokeObjectURL(url);
      showError(new Error("PNG export could not be rendered"));
    };
    image.src = url;
    return;
  }
  const page = window.open("", "_blank");
  if (!page) {
    showError(new Error("Allow pop-ups to use Print / PDF"));
    return;
  }
  page.document.write(
    `<!doctype html><html><head><title>${esc(name)}</title><style>@page{size:landscape;margin:0}html,body{margin:0;width:100%;height:100%;font-family:Helvetica,Arial,sans-serif}svg{display:block;width:100%;height:100%}</style></head><body>${svg}<script>window.onload=()=>window.print()<\/script></body></html>`,
  );
  page.document.close();
}

function exportName(suffix) {
  return `${state.first}_vs_${state.second}_${suffix}`.replace(
    /[^A-Za-z0-9_.-]+/g,
    "_",
  );
}
function bindExports() {
  $("download-logo").addEventListener("click", async () => {
    try {
      exportSvg(await logoPanelSvg(), exportName("motif_logos"));
    } catch (error) {
      showError(error);
    }
  });
  $("download-rank").addEventListener("click", () =>
    exportSvg(serializeSvg($("rank-chart")), exportName("barplot")),
  );
  $("download-volcano").addEventListener("click", () =>
    exportSvg(serializeSvg($("chart")), exportName("volcano")),
  );
  $("download-aggregate").addEventListener("click", () =>
    exportSvg(aggregateGridSvg(), exportName("aggregate")),
  );
  $("download-panel").addEventListener("click", () =>
    exportSvg(combinedPanelSvg(), exportName("combined")),
  );
  $("download-tsv").addEventListener("click", () =>
    downloadBlob(
      new Blob([comparisonTsv()], {
        type: "text/tab-separated-values;charset=utf-8",
      }),
      `${state.first}_vs_${state.second}_fp_tools.tsv`,
    ),
  );
}

function syncRows(source) {
  const value = Math.max(
    2,
    Math.min(200, Math.floor(Number(source.value) || 20)),
  );
  $("rank-rows").value = value;
  $("rank-rows-slider").value = value;
  drawRank();
}
function showError(error) {
  $("status").textContent = `Could not load resource: ${error.message}`;
  console.error(error);
}

async function init() {
  try {
    const metadata = await fetchJson("data/metadata.json");
    state.metadata = metadata;
    const names = metadata.conditions.map((item) => item.name),
      options = names
        .map((name) => `<option value="${esc(name)}">${esc(name)}</option>`)
        .join("");
    $("condition-1").innerHTML = options;
    $("condition-2").innerHTML = options;
    $("condition-1").value = state.first;
    $("condition-2").value = state.second;
    $("download-all").href = metadata.downloads.all_results;
    $("condition-1").addEventListener("change", () =>
      handleConditionChange("first"),
    );
    $("condition-2").addEventListener("change", () =>
      handleConditionChange("second"),
    );
    $("plot-count").addEventListener("change", () => {
      ensureSelected();
      renderAll(false);
    });
    $("rank-rows").addEventListener("input", (event) => syncRows(event.target));
    $("rank-rows-slider").addEventListener("input", (event) =>
      syncRows(event.target),
    );
    bindExports();
    if (window.innerWidth >= 1100 && window.innerHeight < 900)
      $("options").removeAttribute("open");
    await loadComparison(true);
  } catch (error) {
    showError(error);
  }
}

init();
