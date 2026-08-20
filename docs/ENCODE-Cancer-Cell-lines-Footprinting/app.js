"use strict";

const $ = (id) => document.getElementById(id);
const bootstrap = window.fpToolsBrowserBootstrap || { mode: "bundle" };
const plotControls = window.fpToolsPlotControls;
const state = {
  mode: bootstrap.mode === "embedded" ? "embedded" : "bundle",
  review: null,
  metadata: null,
  entry: null,
  comparisonIndex: 0,
  payload: null,
  motifs: [],
  aggregate: new Map(),
  profileAxis: [],
  first: "",
  second: "",
  selected: [],
  active: 0,
  sampleStyles: new Map(),
  payloadCache: new Map(),
  logoDataCache: new Map(),
  colors: { first: "#dc2626", second: "#2563eb", neutral: "#8a94a6" },
  request: 0,
  renderRequest: 0,
  hasAggregates: true,
};
const plotSvgStyle =
  "svg,text{font-family:Helvetica,Arial,sans-serif}.plot-title{font-size:15px;font-weight:900;fill:#172033}.summary-label{font-size:10px;font-weight:700;fill:#64748b}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e3eaf3;stroke-width:1}.tick{font-size:11px;fill:#526176;font-weight:700}.axis-label{font-size:12px;fill:#243247;font-weight:900}.rank-bar.active{stroke:#111827;stroke-width:1.5}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.28))}.volcano-user-label{font-size:12px;font-weight:900;fill:#111827}.volcano-label-line{stroke:#475569;stroke-width:1}";
const aggregateLegendLineWidth = 3;

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
function sampleDisplayName(sample, condition = "") {
  if (sample && typeof sample === "object" && sample.display_name)
    return String(sample.display_name);
  const name = String(sample?.name ?? sample ?? "");
  if (name.includes("::")) return name.split("::").slice(1).join("::");
  const prefix = condition ? `${condition}_` : "";
  return prefix && name.startsWith(prefix) ? name.slice(prefix.length) : name;
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

async function decodeEmbeddedPayload(payloadB64) {
  if (!("DecompressionStream" in window))
    throw new Error(
      "This report needs a modern browser with gzip DecompressionStream support.",
    );
  const bytes = Uint8Array.from(atob(payloadB64), (character) =>
      character.charCodeAt(0),
    ),
    stream = new Blob([bytes])
      .stream()
      .pipeThrough(new DecompressionStream("gzip"));
  return JSON.parse(await new Response(stream).text());
}

function payloadHasAggregates(payload) {
  return Boolean(payload?.aggregate?.motifs?.length);
}

function embeddedMetadata(review) {
  const comparisons = (review.comparisons || []).map((record, index) => {
      const payload = record.payload || {},
        conditions = payload.conditions || ["condition1", "condition2"];
      return {
        comparison: `comparison-${index + 1}`,
        ordinal: index,
        label: String(record.label || `Comparison ${index + 1}`),
        condition1: String(conditions[0] || "condition1"),
        condition2: String(conditions[1] || "condition2"),
        payload,
        aggregate_motifs: payload.aggregate?.motifs?.length || 0,
      };
    }),
    conditions = [];
  comparisons.forEach((record) => {
    [record.condition1, record.condition2].forEach((name) => {
      if (!conditions.some((item) => item.name === name))
        conditions.push({ name, samples: [] });
    });
  });
  return {
    schema: "fp-tools.static-comparison-browser.v1",
    selector_mode: "comparison",
    title: review.title || "Review multiple differential footprint comparisons",
    comparisons,
    conditions,
    default_aggregate_plots: comparisons.some(
      (record) => record.aggregate_motifs,
    )
      ? 4
      : 1,
    default_aggregate_motifs: [],
    documentation_url: "",
    downloads: {},
    logos: {},
  };
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
    rawCiLower = Number(point.ci_lower),
    rawCiUpper = Number(point.ci_upper),
    ciLower = reversed && Number.isFinite(rawCiUpper) ? -rawCiUpper : point.ci_lower,
    ciUpper = reversed && Number.isFinite(rawCiLower) ? -rawCiLower : point.ci_upper,
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
    ci_lower: ciLower,
    ci_upper: ciUpper,
    n_motif_regions_set_1: reversed
      ? point.n_motif_regions_set_2
      : point.n_motif_regions_set_1,
    n_motif_regions_set_2: reversed
      ? point.n_motif_regions_set_1
      : point.n_motif_regions_set_2,
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
  if (!state.hasAggregates) return 1;
  return Math.max(1, Math.min(12, Number($("plot-count").value) || 4));
}
function selectableMotifs() {
  return state.hasAggregates
    ? state.motifs.filter((item) => state.aggregate.has(item.prefix))
    : state.motifs.slice();
}
function sortedMotifs() {
  return selectableMotifs()
    .sort((a, b) =>
      motifLabel(a).localeCompare(motifLabel(b), undefined, {
        sensitivity: "base",
      }),
    );
}

function defaultSelected(target) {
  const candidates = selectableMotifs(),
    positive = candidates
      .filter((item) => item.effect > 0)
      .sort((a, b) => b.effect - a.effect || a.pvalue - b.pvalue);
  const negative = candidates
    .filter((item) => item.effect < 0)
    .sort((a, b) => a.effect - b.effect || a.pvalue - b.pvalue);
  const output = [],
    negativeCount = Math.floor(target / 2),
    positiveCount = target - negativeCount;
  (state.metadata.default_aggregate_motifs || []).forEach((prefix) => {
    if (
      output.length < target &&
      candidates.some((item) => item.prefix === prefix) &&
      !output.includes(prefix)
    )
      output.push(prefix);
  });
  [
    ...positive.slice(0, positiveCount),
    ...negative.slice(0, negativeCount),
  ].forEach((item) => {
    if (!output.includes(item.prefix)) output.push(item.prefix);
  });
  candidates
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
  const collectionTitle = state.metadata?.title || "Differential footprint report",
    title = state.mode === "embedded"
      ? state.entry?.label || collectionTitle
      : collectionTitle;
  $("report-title").textContent = title;
  $("title-cond1").textContent = state.first;
  $("title-cond2").textContent = state.second;
  $("title-cond1").style.color = state.colors.first;
  $("title-cond2").style.color = state.colors.second;
  document.title = `${title} (${state.first} vs ${state.second})`;
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
    alpha: 0.38,
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
          const label = sampleDisplayName(sample, condition);
          return `<div class="sample-style-row"><input type="checkbox" data-sample-visible="${esc(sample)}" ${style.visible ? "checked" : ""} aria-label="Show ${esc(label)}"><span class="sample-style-name" title="${esc(label)}">${esc(label)}</span><input type="color" data-sample-color="${esc(sample)}" value="${style.color}" aria-label="Color for ${esc(label)}"><input type="number" data-sample-alpha="${esc(sample)}" min="0.1" max="1" step="0.1" value="${style.alpha}" aria-label="Opacity for ${esc(label)}"><input type="number" data-sample-width="${esc(sample)}" min="0.3" max="4" step="0.1" value="${style.width}" aria-label="Width for ${esc(label)}"><select data-sample-type="${esc(sample)}" aria-label="Line type for ${esc(label)}"><option value="solid" ${style.type === "solid" ? "selected" : ""}>Solid</option><option value="dash" ${style.type === "dash" ? "selected" : ""}>Dash</option><option value="dot" ${style.type === "dot" ? "selected" : ""}>Dot</option></select></div>`;
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

function embeddedLogoUri(prefix) {
  if (state.mode !== "embedded") return "";
  const record = state.payload?.logos?.[prefix] || {};
  return [record.svg, record.png, record.uri, record.data_uri].find(
    (value) => typeof value === "string" && value.startsWith("data:image/"),
  ) || "";
}

function motifLogoSvg(prefix, attributes = "") {
  const counts = state.payload?.motif_matrices?.[prefix];
  if (
    !Array.isArray(counts) ||
    counts.length !== 4 ||
    !Array.isArray(counts[0]) ||
    !counts[0].length
  )
    return "";
  const width = 420,
    height = 150,
    bases = ["A", "C", "G", "T"],
    colors = { A: "#198754", C: "#0d6efd", G: "#f59f00", T: "#dc3545" },
    left = 46,
    right = 14,
    top = 16,
    bottom = 32,
    plotWidth = width - left - right,
    plotHeight = height - top - bottom,
    positions = counts[0].length,
    columnWidth = plotWidth / Math.max(1, positions),
    bits = [[], [], [], []],
    attributeText = attributes ? ` ${attributes}` : "",
    parts = [
      `<svg${attributeText} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${esc(prefix)} motif logo"><rect width="100%" height="100%" fill="#fff"/><line x1="${left}" y1="${top + plotHeight}" x2="${left + plotWidth}" y2="${top + plotHeight}" stroke="#3b4552" stroke-width="1.2"/><line x1="${left}" y1="${top}" x2="${left}" y2="${top + plotHeight}" stroke="#3b4552" stroke-width="1.2"/><text x="18" y="${top + plotHeight / 2}" transform="rotate(-90 18 ${top + plotHeight / 2})" text-anchor="middle" font-size="12" font-weight="700" fill="#152133">bits</text><text x="${left + plotWidth / 2}" y="${height - 7}" text-anchor="middle" font-size="12" font-weight="700" fill="#152133">position</text>`,
    ];
  for (let position = 0; position < positions; position += 1) {
    const total =
        counts.reduce(
          (sum, row) => sum + (Number(row[position]) || 0),
          0,
        ) || 1,
      frequencies = counts.map(
        (row) => (Number(row[position]) || 0) / total,
      ),
      entropy = -frequencies.reduce(
        (sum, value) =>
          sum +
          (value > 0 ? value * Math.log2(Math.max(value, 1e-12)) : 0),
        0,
      ),
      information = Math.max(0, 2 - entropy);
    frequencies.forEach(
      (frequency, base) => (bits[base][position] = frequency * information),
    );
  }
  [0, 1, 2].forEach((tick) => {
    const y = top + plotHeight - (tick / 2) * plotHeight;
    parts.push(
      `<line x1="${left - 4}" y1="${y.toFixed(2)}" x2="${left}" y2="${y.toFixed(2)}" stroke="#3b4552"/><text x="${left - 8}" y="${(y + 4).toFixed(2)}" text-anchor="end" font-size="11" font-weight="700" fill="#56616f">${tick}</text>`,
    );
  });
  for (let position = 0; position < positions; position += 1) {
    let y = top + plotHeight;
    const order = [0, 1, 2, 3].sort(
        (first, second) => bits[first][position] - bits[second][position],
      ),
      x = left + position * columnWidth + columnWidth / 2;
    if (
      positions <= 18 ||
      position === 0 ||
      position === positions - 1 ||
      (position + 1) % 5 === 0
    )
      parts.push(
        `<text x="${x.toFixed(2)}" y="${top + plotHeight + 13}" text-anchor="middle" font-size="9" font-weight="700" fill="#56616f">${position + 1}</text>`,
      );
    order.forEach((baseIndex) => {
      const value = bits[baseIndex][position];
      if (value <= 0.015) return;
      const letterHeight = Math.max(3, (value / 2) * plotHeight),
        base = bases[baseIndex],
        fontSize = Math.max(8, Math.min(40, letterHeight * 1.25));
      y -= letterHeight;
      parts.push(
        `<text x="${x.toFixed(2)}" y="${(y + letterHeight * 0.88).toFixed(2)}" text-anchor="middle" font-family="Arial Black,Helvetica,Arial,sans-serif" font-size="${fontSize.toFixed(2)}" font-weight="900" fill="${colors[base]}">${base}</text>`,
      );
    });
  }
  parts.push("</svg>");
  return parts.join("");
}

function motifLogoHtml(prefix) {
  const svg = motifLogoSvg(prefix);
  if (svg) return svg;
  const embedded = embeddedLogoUri(prefix);
  if (embedded)
    return `<img alt="${esc(prefix)} motif logo" src="${esc(embedded)}" width="1000" height="250" style="display:block;max-width:100%;max-height:60px;width:auto;height:auto;object-fit:contain">`;
  if (state.mode === "embedded")
    return '<span class="logo-empty">Logo unavailable</span>';
  return `<img alt="${esc(prefix)} motif logo" src="${esc(logoPath(prefix))}" width="1000" height="250" style="display:block;max-width:100%;max-height:60px;width:auto;height:auto;object-fit:contain">`;
}

function logoDataUri(prefix) {
  if (state.logoDataCache.has(prefix)) return state.logoDataCache.get(prefix);
  const svg = motifLogoSvg(prefix);
  if (svg) {
    const uri = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
    const request = Promise.resolve(uri);
    state.logoDataCache.set(prefix, request);
    return request;
  }
  const embedded = embeddedLogoUri(prefix);
  if (embedded) {
    const request = Promise.resolve(embedded);
    state.logoDataCache.set(prefix, request);
    return request;
  }
  if (state.mode === "embedded") {
    const label = motifLabel(
        state.motifs.find((item) => item.prefix === prefix) || { name: prefix },
      ),
      placeholder = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 420 150"><rect width="420" height="150" fill="#fff"/><rect x="1" y="1" width="418" height="148" rx="7" fill="none" stroke="#d8e2ef"/><text x="210" y="78" text-anchor="middle" font-family="Helvetica,Arial,sans-serif" font-size="16" font-weight="700" fill="#64748b">${esc(label)}</text></svg>`,
      request = Promise.resolve(
        `data:image/svg+xml;charset=utf-8,${encodeURIComponent(placeholder)}`,
      );
    state.logoDataCache.set(prefix, request);
    return request;
  }
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
      return `<article class="selected-motif${index === state.active ? " active" : ""}" data-selected-panel="${index}"><div class="selected-head"><select class="panel-tf" data-panel-tf="${index}" aria-label="Selected motif ${index + 1}">${options}</select></div><div class="motif-logo">${motifLogoHtml(prefix)}</div><div class="detail-grid"><p class="motif-group" style="color:${colorFor(motif)}">${esc(group)}</p><p class="metric-line">ΔFP = ${fmt(motif.effect, 4)}</p><p class="metric-line">FDR = ${fmtSci(motif.qvalue)}</p></div></article>`;
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

function rankMode() {
  return $("rank-sort-toggle").checked ? "significance" : "effect";
}

function comparisonTitle() {
  return String(
    state.entry?.label || `${state.first || "condition1"} vs ${state.second || "condition2"}`,
  );
}

function volcanoLabelLayout(items, sx, sy, bounds) {
  const minimumGap = 15,
    middle = (bounds.left + bounds.right) / 2,
    groups = { left: [], right: [] };
  items.forEach((item) => {
    const pointX = sx(item.effect),
      pointY = sy(item.neglog10p),
      side = pointX > middle ? "left" : "right";
    groups[side].push({ item, pointX, pointY, labelY: pointY, side });
  });
  Object.values(groups).forEach((rows) => {
    rows.sort((a, b) => a.pointY - b.pointY);
    rows.forEach((row, index) => {
      row.labelY = Math.max(
        bounds.top,
        row.pointY,
        index ? rows[index - 1].labelY + minimumGap : bounds.top,
      );
    });
    for (let index = rows.length - 1; index >= 0; index -= 1) {
      const maximum = index === rows.length - 1
        ? bounds.bottom
        : rows[index + 1].labelY - minimumGap;
      rows[index].labelY = Math.min(rows[index].labelY, maximum);
    }
  });
  return [...groups.left, ...groups.right];
}

function renderVolcano() {
  const width = 760,
    height = 760,
    margin = { top: 54, right: 48, bottom: 60, left: 84 },
    innerWidth = width - margin.left - margin.right,
    innerHeight = innerWidth,
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
    selected = $("volcano-highlight").value === "none"
      ? new Set()
      : visibleSelected(),
    parts = [
      `<style>${plotSvgStyle}</style><text x="${width / 2}" y="24" class="plot-title" text-anchor="middle">Comparison: ${esc(comparisonTitle())}</text><rect x="${margin.left}" y="${margin.top}" width="${innerWidth}" height="${innerHeight}" fill="none" stroke="#d9e2ec"/>`,
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
  const labelItems = plotControls.matchingMotifs(
      state.motifs,
      $("volcano-labels").value,
    ),
    labelLayout = volcanoLabelLayout(labelItems, sx, sy, {
      left: margin.left + 8,
      right: margin.left + innerWidth - 8,
      top: margin.top + 10,
      bottom: margin.top + innerHeight - 10,
    });
  labelLayout.forEach(({ item, pointX, pointY, labelY, side }) => {
    const direction = side === "right" ? 1 : -1,
      labelX = pointX + direction * 12,
      anchor = side === "right" ? "start" : "end";
    parts.push(
      `<line class="volcano-label-line" x1="${pointX.toFixed(2)}" y1="${pointY.toFixed(2)}" x2="${labelX.toFixed(2)}" y2="${labelY.toFixed(2)}" stroke="#475569" stroke-width="1"/><text class="volcano-user-label" x="${(labelX + direction * 2).toFixed(2)}" y="${(labelY + 4).toFixed(2)}" text-anchor="${anchor}" font-family="Helvetica,Arial,sans-serif" font-size="12" font-weight="900" fill="#111827" stroke="none">${esc(item.name || motifLabel(item))}</text>`,
    );
  });
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
    mode = rankMode(),
    ranked = plotControls.rankMotifs(state.motifs, mode, limit),
    positive = ranked.positive,
    negative = ranked.negative,
    shown = [...negative, ...positive],
    width = 380,
    rowHeight = 14,
    rowGap = 3,
    sectionGap = 8,
    margin = { top: 110, bottom: 68, left: 128, right: 14 },
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
      Math.max(
        ...shown.map((item) => Math.abs(plotControls.rankMetric(item, mode))),
        1e-9,
      ),
    ),
    sx = (value) => xMiddle + (value / maxAbs) * xWidth,
    axisY = height - 60,
    selected = visibleSelected(),
    effectColorMax = Math.max(
      ...shown.map((item) => Math.abs(item.effect)),
      1e-9,
    ),
    significanceColorMax = Math.max(
      ...shown.map((item) => plotControls.negLog10P(item)),
      1e-9,
    ),
    colorDomain = mode === "significance"
      ? { maxAbs: effectColorMax }
      : { max: significanceColorMax },
    colorOptions = {
      first: state.colors.first,
      second: state.colors.second,
      neutralCenter: "#f8fafc",
    },
    axisLabel = mode === "significance"
      ? "Signed −log10(p-value)"
      : state.payload?.change_label || "Differential footprint score",
    legendLabel = mode === "significance"
      ? state.payload?.change_label || "Differential footprint score"
      : "−log10(p-value)",
    legendLow = mode === "significance"
      ? -effectColorMax
      : significanceColorMax,
    legendCenter = 0,
    legendHigh = mode === "significance"
      ? effectColorMax
      : significanceColorMax,
    gradientStops = `<stop offset="0%" stop-color="${state.colors.second}"/><stop offset="50%" stop-color="#f8fafc"/><stop offset="100%" stop-color="${state.colors.first}"/>`,
    parts = [
      `<style>${plotSvgStyle}</style><defs><linearGradient id="rank-color-gradient" x1="0" x2="1">${gradientStops}</linearGradient></defs><text x="${width / 2}" y="16" class="plot-title" text-anchor="middle">Top differential motifs</text><text x="${width / 2}" y="34" class="summary-label" text-anchor="middle">Comparison: ${esc(comparisonTitle())}</text><text x="8" y="49" class="summary-label">Color: ${esc(legendLabel)}</text><rect x="8" y="54" width="104" height="7" rx="2" fill="url(#rank-color-gradient)"/><text x="8" y="72" class="tick">${fmt(legendLow, 2)}</text><text x="60" y="72" class="tick" text-anchor="middle">${fmt(legendCenter, 2)}</text><text x="112" y="72" class="tick" text-anchor="end">${fmt(legendHigh, 2)}</text><line x1="${xMiddle}" y1="${margin.top - 20}" x2="${xMiddle}" y2="${axisY}" stroke="#172033" stroke-width="2.2"/><text x="${xMiddle - 6}" y="${margin.top - 27}" text-anchor="end" font-size="14" font-weight="900" fill="${state.colors.second}">${esc(state.second)}_up</text><text x="${xMiddle + 6}" y="${margin.top - 27}" text-anchor="start" font-size="14" font-weight="900" fill="${state.colors.first}">${esc(state.first)}_up</text>`,
    ];
  niceTicks(-maxAbs, maxAbs, 5).forEach((value) =>
    parts.push(
      `<line x1="${sx(value)}" y1="${axisY - 4}" x2="${sx(value)}" y2="${axisY + 4}" class="axis"/><text x="${sx(value)}" y="${axisY + 17}" class="tick" text-anchor="middle">${fmt(value, 3)}</text>`,
    ),
  );
  parts.push(
    `<line x1="${sx(-maxAbs)}" y1="${axisY}" x2="${sx(maxAbs)}" y2="${axisY}" class="axis"/><text x="${xMiddle}" y="${height - 8}" class="axis-label" text-anchor="middle">${esc(axisLabel)}</text>`,
  );
  let y = margin.top;
  const drawRows = (rows) =>
    rows.forEach((item) => {
      const metric = plotControls.rankMetric(item, mode),
        opposite = plotControls.oppositeMetric(item, mode),
        barWidth = (Math.abs(metric) / maxAbs) * xWidth,
        x = metric >= 0 ? xMiddle : xMiddle - barWidth,
        isSelected = selected.has(item.prefix),
        name = motifLabel(item).slice(0, 20),
        labelY = y + rowHeight - 2,
        fill = plotControls.rankColor(item, mode, colorDomain, colorOptions);
      parts.push(
        `<text class="rank-name${isSelected ? " active" : ""}" data-prefix="${esc(item.prefix)}" x="6" y="${labelY}" font-size="10" font-weight="${isSelected ? 900 : 700}" fill="${isSelected ? fill : "#526176"}">${esc(name)}</text><rect class="rank-bar${isSelected ? " active" : ""}" data-prefix="${esc(item.prefix)}" x="${x}" y="${y}" width="${Math.max(1, barWidth)}" height="${rowHeight}" fill="${fill}" fill-opacity="${isSelected ? 1 : 0.82}"><title>${esc(motifLabel(item))}: ΔFP ${fmt(item.effect, 4)}; −log10(p-value) ${fmt(plotControls.negLog10P(item), 3)}</title></rect><text x="${metric >= 0 ? x - 3 : x + barWidth + 3}" y="${labelY}" class="tick" text-anchor="${metric >= 0 ? "end" : "start"}">${fmt(opposite, mode === "significance" ? 3 : 2)}</text>`,
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
    sampleMeta = {},
    conditionCounts = {};
  motif.conditions.forEach((condition) =>
    {
      conditionCounts[condition.name] = Number(condition.n_sites || 0);
      condition.samples.forEach((sample) => {
        samples[sample.name] = sample.profile;
        sampleMeta[sample.name] = sample;
      });
    },
  );
  return {
    samples,
    sampleMeta,
    conditionCounts,
    n_profile_sites: motif.n_sites || 0,
  };
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
  const means = [state.first, state.second]
      .map((condition) => {
        const rows = series.filter((item) => item.condition === condition);
        if (!rows.length) return null;
        return {
          condition,
          profile: rows[0].profile.map(
            (_value, point) =>
              rows.reduce((sum, row) => sum + row.profile[point], 0) /
              rows.length,
          ),
        };
      })
      .filter(Boolean),
    values = series.flatMap((item) => item.profile),
    rawMin = Math.min(...values, 0),
    rawMax = Math.max(...values, 1e-9),
    padding = Math.max((rawMax - rawMin || 1) * 0.18, 1e-6),
    yMin = rawMin - padding,
    yMax = rawMax + padding,
    width = 300,
    height = 300,
    margin = { top: 42, right: 8, bottom: 34, left: 36 },
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
      `<svg class="aggregate-panel" data-panel="${index}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg"><style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><text x="${width / 2}" y="17" class="plot-title" text-anchor="middle">${esc(motifLabel(motif))}</text><text x="${width / 2}" y="32" class="summary-label" text-anchor="middle">regions: ${Number(record.conditionCounts[state.first] || 0).toLocaleString()} / ${Number(record.conditionCounts[state.second] || 0).toLocaleString()}</text>`,
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
    `<line x1="${margin.left}" y1="${margin.top + innerHeight}" x2="${margin.left + innerWidth}" y2="${margin.top + innerHeight}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + innerHeight}" class="axis"/>`,
  );
  series
    .sort((a, b) => b.fpScore - a.fpScore)
    .forEach((item) =>
      parts.push(
        `<path d="${linePath(item.profile, axis, sx, sy)}" fill="none" stroke="${item.style.color}" stroke-width="${item.style.width}"${dashAttribute(item.style.type)} stroke-opacity="${item.style.alpha}"><title>${esc(sampleDisplayName(record.sampleMeta[item.sample], item.condition))}</title></path>`,
      ),
    );
  means.forEach((item) =>
    parts.push(
      `<path d="${linePath(item.profile, axis, sx, sy)}" fill="none" stroke="${item.condition === state.first ? state.colors.first : state.colors.second}" stroke-width="3" stroke-opacity="1"><title>${esc(item.condition)} mean</title></path>`,
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
  if (count <= 8) return { columns: 4, rows: 2 };
  if (count <= 9) return { columns: 3, rows: 3 };
  return { columns: 4, rows: 3 };
}

function legendGroups() {
  return [state.first, state.second]
    .map((condition) => ({
      condition,
      rows: conditionSamples(condition)
        .map((sample, index) => ({
          sample,
          condition,
          style: sampleStyle(sample, condition, index),
        }))
        .filter((row) => row.style.visible),
    }))
    .filter((group) => group.rows.length);
}

function renderLegend() {
  const groups = legendGroups();
  $("aggregate-legend").innerHTML = groups
    .map(
      (group) =>
        `<div class="legend-group"><div class="legend-group-title" title="${esc(group.condition)}">${esc(group.condition)}</div><div class="legend-row"><i class="legend-line" style="border-top-color:${group.condition === state.first ? state.colors.first : state.colors.second};border-top-width:3px"></i><span>Mean</span></div>${group.rows
          .map(
            (row) =>
              `<div class="legend-row"><i class="legend-line" style="border-top-color:${row.style.color};border-top-width:${aggregateLegendLineWidth}px;border-top-style:${row.style.type === "dash" ? "dashed" : row.style.type === "dot" ? "dotted" : "solid"};opacity:${row.style.alpha}"></i><span title="${esc(sampleDisplayName(row.sample, row.condition))}">${esc(sampleDisplayName(row.sample, row.condition))}</span></div>`,
          )
          .join("")}</div>`,
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
  if (state.hasAggregates && !state.aggregate.has(prefix)) {
    const motif = state.motifs.find((item) => item.prefix === prefix);
    $("status").textContent = `${motifLabel(motif || { name: prefix })} has a statistical result, but no embedded aggregate profile.`;
    return;
  }
  state.selected[state.active] = prefix;
  renderAll(false);
}

function renderAll(refreshControls = true) {
  ensureSelected();
  updateHeader();
  if (refreshControls) {
    renderColorControls();
    if (state.hasAggregates) renderSampleStyles();
  }
  renderSelectedCards();
  drawRank();
  renderVolcano();
  if (state.hasAggregates) {
    renderLegend();
    renderAggregateGrid();
  }
}

async function loadComparison(reset = true) {
  let entry,
    payload,
    first,
    second,
    reversed = false;
  if (state.mode === "embedded") {
    state.comparisonIndex = Math.max(
      0,
      Math.min(
        state.metadata.comparisons.length - 1,
        Number($("comparison-selector").value) || 0,
      ),
    );
    entry = state.metadata.comparisons[state.comparisonIndex];
    payload = entry?.payload;
    [first, second] = payload?.conditions || [];
  } else {
    first = $("condition-1").value;
    second = $("condition-2").value;
    if (first === second) return;
    ({ entry } = comparisonEntry(first, second));
  }
  if (!entry || (!payload && state.mode === "embedded"))
    throw new Error("The selected comparison payload is unavailable");
  state.first = first;
  state.second = second;
  const token = ++state.request;
  $("status").textContent = `Loading ${entry.label || `${first} vs ${second}`}…`;
  if (state.mode === "bundle") payload = await fetchGzipJsonCached(entry.file);
  if (token !== state.request) return;
  state.entry = entry;
  state.payload = payload;
  state.hasAggregates = payloadHasAggregates(payload);
  document.body.classList.toggle("no-aggregate", !state.hasAggregates);
  state.aggregate = new Map(
    (payload.aggregate?.motifs || []).map((item) => [item.prefix, item]),
  );
  state.profileAxis = payload.aggregate?.x || [];
  if (state.mode === "bundle") reversed = payload.conditions[0] !== first;
  if (
    new Set(payload.conditions).size !== 2 ||
    !payload.conditions.includes(first) ||
    !payload.conditions.includes(second)
  )
    throw new Error(`Payload conditions do not match ${first} and ${second}`);
  state.motifs = payload.points.map((item) => orientedMotif(item, reversed));
  const colors = payload.colors || {};
  state.colors = {
    first: colors[`${first}_up`] || "#dc2626",
    second: colors[`${second}_up`] || "#2563eb",
    neutral: colors["n.s."] || "#8a94a6",
  };
  state.sampleStyles = new Map();
  state.logoDataCache = new Map();
  state.active = 0;
  ensureSelected(reset);
  renderAll(true);
  const significant = state.motifs.filter((item) => item.significant).length;
  $("status").textContent =
    `${entry.label || `${first} vs ${second}`} | ${state.motifs.length.toLocaleString()} motifs | ${significant.toLocaleString()} significant | ${first} minus ${second}`;
}

function handleConditionChange(changed) {
  if (changed === "first") {
    const first = $("condition-1").value,
      partners = availablePartners(first),
      preferred = partners.includes(state.second) ? state.second : partners[0];
    $("condition-2").innerHTML = optionMarkup(partners);
    $("condition-2").value = preferred || "";
  }
  loadComparison(true).catch(showError);
}

function optionMarkup(names) {
  return names.map((name) => `<option value="${esc(name)}">${esc(name)}</option>`).join("");
}

function availablePartners(condition) {
  const partners = [];
  state.metadata.comparisons.forEach((item) => {
    const partner = item.condition1 === condition ? item.condition2 : item.condition2 === condition ? item.condition1 : null;
    if (partner && !partners.includes(partner)) partners.push(partner);
  });
  return partners;
}

function comparisonTsv() {
  const columns = [
      "condition1",
      "condition2",
      "prefix",
      "name",
      "motif_id",
      "group",
      "n_profile_sites",
      "n_motif_regions_condition1",
      "n_motif_regions_condition2",
      "effect",
      "ci_lower",
      "ci_upper",
      "pvalue",
      "qvalue",
      "significant",
      "statistical_method",
    ],
    rows = [columns.join("\t")];
  state.motifs.forEach((motif) =>
    rows.push(
      [
        state.first,
        state.second,
        motif.prefix ?? "",
        motif.name ?? "",
        motif.motif_id ?? "",
        motif.group ?? "",
        motif.n_sites ?? "",
        motif.n_motif_regions_set_1 ?? "",
        motif.n_motif_regions_set_2 ?? "",
        motif.effect ?? "",
        motif.ci_lower ?? "",
        motif.ci_upper ?? "",
        motif.pvalue ?? "",
        motif.qvalue ?? "",
        motif.significant ?? "",
        motif.statistical_method ?? "",
      ].join("\t"),
    ),
  );
  return rows.join("\n") + "\n";
}

function allComparisonsTsv() {
  const columns = [
      "comparison",
      "condition1",
      "condition2",
      "prefix",
      "name",
      "motif_id",
      "group",
      "effect",
      "ci_lower",
      "ci_upper",
      "pvalue",
      "qvalue",
      "significant",
      "statistical_method",
    ],
    rows = [columns.join("\t")];
  (state.review?.comparisons || []).forEach((record, index) => {
    const payload = record.payload || {},
      conditions = payload.conditions || ["condition1", "condition2"],
      label = record.label || `Comparison ${index + 1}`;
    (payload.points || []).forEach((point) =>
      rows.push(
        [
          label,
          conditions[0] ?? "",
          conditions[1] ?? "",
          point.prefix ?? "",
          point.name ?? "",
          point.motif_id ?? "",
          point.group ?? "",
          point.change ?? "",
          point.ci_lower ?? "",
          point.ci_upper ?? "",
          point.pvalue ?? "",
          point.fdr ?? "",
          point.group !== "n.s.",
          point.statistical_method ?? "",
        ].join("\t"),
      ),
    );
  });
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
    groups = legendGroups(),
    legendRows = Math.max(0, ...groups.map((group) => group.rows.length + 1)),
    legendHeight = groups.length ? 25 + legendRows * 16 : 0,
    totalWidth = gridWidth,
    totalHeight = gridHeight + legendHeight,
    parts = [
      `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalWidth} ${totalHeight}" font-family="Helvetica,Arial,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/>`,
    ];
  svgs.forEach((svg, index) => {
    const clone = styledClone(svg);
    clone.querySelector("style")?.remove();
    parts.push(
      `<g transform="translate(${(index % shape.columns) * plotWidth},${legendHeight + Math.floor(index / shape.columns) * plotHeight})">${clone.innerHTML}</g>`,
    );
  });
  if (groups.length) {
    const groupWidth = gridWidth / groups.length;
    parts.push(
      `<g class="aggregate-export-legend"><rect x="1" y="1" width="${gridWidth - 2}" height="${legendHeight - 4}" rx="5" fill="#fff" stroke="#d8e2ef"/>`,
    );
    groups.forEach((group, groupIndex) => {
      const x = groupIndex * groupWidth + 10;
      parts.push(
        `<text x="${x}" y="14" class="tick" font-weight="900">${esc(group.condition)}</text>`,
      );
      const meanColor =
        group.condition === state.first ? state.colors.first : state.colors.second;
      parts.push(
        `<line x1="${x}" y1="25" x2="${x + 30}" y2="25" stroke="${meanColor}" stroke-width="3"/><text x="${x + 36}" y="28" class="tick">Mean</text>`,
      );
      group.rows.forEach((row, index) => {
        const y = 44 + index * 16,
          dash =
          row.style.type === "dash"
            ? ' stroke-dasharray="7 4"'
            : row.style.type === "dot"
              ? ' stroke-dasharray="2 3"'
              : "";
        parts.push(
          `<line x1="${x}" y1="${y - 3}" x2="${x + 30}" y2="${y - 3}" stroke="${row.style.color}" stroke-width="${aggregateLegendLineWidth}"${dash} stroke-opacity="${row.style.alpha}"/><text x="${x + 36}" y="${y}" class="tick">${esc(sampleDisplayName(row.sample, row.condition))}</text>`,
        );
      });
    });
    parts.push("</g>");
  }
  parts.push("</svg>");
  return parts.join("");
}

function combinedPanelSvg() {
  const rank = styledClone($("rank-chart")),
    volcano = styledClone($("chart")),
    rankBox = rank.viewBox.baseVal,
    rankWidth = rankBox.width || 380,
    rankHeight = rankBox.height || 600,
    volcanoSize = 760,
    panelHeight = Math.max(760, rankHeight),
    rankScale = panelHeight / rankHeight,
    rankDisplayWidth = rankWidth * rankScale,
    gap = 20;
  rank.querySelector("style")?.remove();
  volcano.querySelector("style")?.remove();
  if (!state.hasAggregates) {
    const totalWidth = rankDisplayWidth + volcanoSize + gap;
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalWidth} ${panelHeight}" font-family="Helvetica,Arial,sans-serif"><style>${plotSvgStyle}</style><rect width="100%" height="100%" fill="#fff"/><g transform="scale(${rankScale})">${rank.innerHTML}</g><g transform="translate(${rankDisplayWidth + gap},0)">${volcano.innerHTML}</g></svg>`;
  }
  const aggregate = aggregateGridSvg(),
    aggregateDocument = new DOMParser().parseFromString(
      aggregate,
      "image/svg+xml",
    ).documentElement,
    aggregateBox = aggregateDocument.viewBox.baseVal,
    aggregateScale = panelHeight / (aggregateBox.height || 600),
    aggregateWidth = (aggregateBox.width || 600) * aggregateScale,
    totalWidth = rankDisplayWidth + volcanoSize + aggregateWidth + gap * 2;
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
  const base = state.mode === "embedded" && state.entry?.label
    ? state.entry.label
    : `${state.first}_vs_${state.second}`;
  return `${base}_${suffix}`.replace(
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
      `${exportName("fp_tools")}.tsv`,
    ),
  );
  if (state.mode === "embedded")
    $("download-all").addEventListener("click", (event) => {
      event.preventDefault();
      downloadBlob(
        new Blob([allComparisonsTsv()], {
          type: "text/tab-separated-values;charset=utf-8",
        }),
        "fp_tools_all_comparisons.tsv",
      );
    });
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
    let metadata;
    if (state.mode === "embedded") {
      state.review = await decodeEmbeddedPayload(bootstrap.payloadB64 || "");
      if (state.review.schema !== "fp-tools.review-multi-comparisons.v1")
        throw new Error("Unsupported embedded review payload");
      if (!state.review.comparisons?.length)
        throw new Error("The embedded review contains no comparisons");
      metadata = embeddedMetadata(state.review);
    } else {
      metadata = await fetchJson("data/metadata.json");
    }
    state.metadata = metadata;
    if (metadata.documentation_url) {
      $("documentation-return").href = metadata.documentation_url;
      $("documentation-return").hidden = false;
    }
    const initialPlotCount = Math.max(
      1,
      Math.min(12, Number(metadata.default_aggregate_plots) || 4),
    );
    $("plot-count").value = String(initialPlotCount);
    if (state.mode === "embedded") {
      $("comparison-selector-control").hidden = false;
      $("condition-selector-controls").hidden = true;
      $("comparison-selector").innerHTML = metadata.comparisons
        .map(
          (record, index) =>
            `<option value="${index}">${esc(record.label)}</option>`,
        )
        .join("");
      $("comparison-selector").addEventListener("change", () =>
        loadComparison(true).catch(showError),
      );
      $("download-all").removeAttribute("href");
    } else {
      const names = metadata.conditions.map((item) => item.name);
      $("condition-1").innerHTML = optionMarkup(names);
      const preferred = metadata.default_comparison ||
        metadata.comparisons.find(
          (record) => record.condition1 === "K562" && record.condition2 === "HepG2"
        ) || metadata.comparisons[0];
      state.first = preferred?.condition1 || names[0] || "";
      const partners = availablePartners(state.first);
      state.second = preferred?.condition2 || partners[0] || "";
      $("condition-2").innerHTML = optionMarkup(partners);
      $("condition-1").value = state.first;
      $("condition-2").value = state.second;
      $("download-all").href = metadata.downloads.all_results;
      $("condition-1").addEventListener("change", () =>
        handleConditionChange("first"),
      );
      $("condition-2").addEventListener("change", () =>
        handleConditionChange("second"),
      );
    }
    $("plot-count").addEventListener("change", () => {
      ensureSelected();
      renderAll(false);
    });
    $("rank-rows").addEventListener("input", (event) => syncRows(event.target));
    $("rank-rows-slider").addEventListener("input", (event) =>
      syncRows(event.target),
    );
    $("rank-sort-toggle").addEventListener("change", drawRank);
    $("volcano-highlight").addEventListener("change", renderVolcano);
    $("volcano-labels").addEventListener("input", renderVolcano);
    bindExports();
    if (window.innerWidth >= 1100 && window.innerHeight < 900)
      $("options").removeAttribute("open");
    await loadComparison(true);
  } catch (error) {
    showError(error);
  }
}

init();
