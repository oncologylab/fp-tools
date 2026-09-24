(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  root.fpToolsPlotControls = api;
})(typeof globalThis === "object" ? globalThis : this, function () {
  "use strict";

  function number(value, fallback = 0) {
    const parsed = Number(value);
    return Number.isFinite(parsed) ? parsed : fallback;
  }

  function negLog10P(motif) {
    const embedded = Number(motif?.neglog10p);
    if (Number.isFinite(embedded)) return Math.max(0, embedded);
    const pvalue = Math.max(1e-300, number(motif?.pvalue, 1));
    return -Math.log10(pvalue);
  }

  function rankMetric(motif, mode) {
    const effect = number(motif?.effect ?? motif?.change);
    if (mode !== "significance") return effect;
    return (effect < 0 ? -1 : 1) * negLog10P(motif);
  }

  function oppositeMetric(motif, mode) {
    return mode === "significance"
      ? number(motif?.effect ?? motif?.change)
      : negLog10P(motif);
  }

  function rankMotifs(motifs, mode, limit) {
    const total = Math.max(1, Math.floor(number(limit, 20))),
      negativeCount = Math.floor(total / 2),
      positiveCount = total - negativeCount,
      significanceOrder = (a, b) =>
        negLog10P(b) - negLog10P(a) ||
        Math.abs(number(b.effect ?? b.change)) -
          Math.abs(number(a.effect ?? a.change)) ||
        String(a.prefix || "").localeCompare(String(b.prefix || "")),
      positive = motifs.filter(
        (item) => number(item.effect ?? item.change) > 0,
      ),
      negative = motifs.filter(
        (item) => number(item.effect ?? item.change) < 0,
      );
    if (mode === "significance") {
      positive.sort(significanceOrder);
      negative.sort(significanceOrder);
    } else {
      positive.sort(
        (a, b) =>
          number(b.effect ?? b.change) - number(a.effect ?? a.change) ||
          number(a.pvalue, 1) - number(b.pvalue, 1),
      );
      negative.sort(
        (a, b) =>
          number(a.effect ?? a.change) - number(b.effect ?? b.change) ||
          number(a.pvalue, 1) - number(b.pvalue, 1),
      );
    }
    const n = Math.min(negative.length, negativeCount + Math.max(0, positiveCount - positive.length));
    const p = Math.min(positive.length, positiveCount + Math.max(0, negativeCount - negative.length));
    return { negative: negative.slice(0, n), positive: positive.slice(0, p) };
  }

  function numeric(value) {
    return value !== null && value !== undefined && value !== '' &&
      typeof value !== 'boolean' && Number.isFinite(Number(value));
  }

  function passesDisplayFilter(point, filter) {
    const effect = point.effect ?? point.change;
    const probability = point[filter.metric];
    return numeric(effect) && numeric(probability) &&
      Number(probability) >= 0 && Number(probability) <= 1 &&
      Number(probability) <= filter.alpha && Math.abs(Number(effect)) >= filter.delta;
  }

  function validRange(range) {
    return Array.isArray(range) && range.length === 2 && range.every(numeric) &&
      Number(range[0]) < Number(range[1]);
  }

  function autoDomain(values, mode = 'signed') {
    let lo = 0, hi = 0;
    for (const value of values) if (numeric(value)) {
      lo = Math.min(lo, Number(value)); hi = Math.max(hi, Number(value));
    }
    if (mode === 'symmetric') {
      const bound = Math.max(Math.abs(lo), hi) * 1.05 || 1;
      return [-bound, bound];
    }
    if (mode === 'positive') return [0, hi * 1.05 || 1];
    const pad = Math.max((hi - lo) * .18, 1e-6);
    return [lo - pad, hi + pad];
  }

  // A reusable DOM control; callback receives only valid ranges or null (reset).
  function axisEditor(label, range, custom, onChange, nonnegative = false) {
    const box = document.createElement('details');
    box.className = 'axis-editor';
    const summary = document.createElement('summary');
    summary.textContent = `${label}: ${range.map(v => Number(v.toPrecision(5))).join(' to ')}${custom ? ' · Custom range' : ' · Auto'}`;
    box.append(summary);
    const span = range[1] - range[0], bounds = [range[0] - span, range[1] + span];
    if (nonnegative) bounds[0] = Math.max(0, bounds[0]);
    const numbers = [], sliders = [], error = document.createElement('span');
    error.setAttribute('role', 'status');
    for (let i = 0; i < 2; i++) {
      const row = document.createElement('label');
      row.textContent = i ? 'Maximum ' : 'Minimum ';
      const slider = document.createElement('input'), input = document.createElement('input');
      slider.type = 'range'; slider.min = bounds[0]; slider.max = bounds[1];
      slider.step = span / 1000; slider.value = range[i];
      input.type = 'number'; input.step = 'any'; input.value = Number(range[i].toPrecision(8));
      if (nonnegative) input.min = 0;
      slider.setAttribute('aria-label', `${label} ${i ? 'maximum' : 'minimum'} slider`);
      input.setAttribute('aria-label', `${label} ${i ? 'maximum' : 'minimum'}`);
      row.append(slider, input); box.append(row); numbers.push(input); sliders.push(slider);
      const commit = () => {
        const next = numbers.map(n => n.value === '' ? '' : Number(n.value));
        const valid = validRange(next) && (!nonnegative || next[0] >= 0);
        numbers.forEach(n => n.setAttribute('aria-invalid', String(!valid)));
        error.textContent = valid ? '' : 'Enter finite limits with minimum below maximum.';
        if (!valid) return;
        sliders.forEach((s, j) => {
          s.min = Math.min(Number(s.min), next[0]); s.max = Math.max(Number(s.max), next[1]); s.value = next[j];
        });
        summary.textContent = `${label}: ${next.map(v => Number(v.toPrecision(5))).join(' to ')} · Custom range`;
        onChange(next);
      };
      slider.addEventListener('input', () => {input.value = slider.value; commit();});
      input.addEventListener('change', commit);
    }
    const reset = document.createElement('button'); reset.type = 'button';
    reset.textContent = `Reset ${label}`; reset.addEventListener('click', () => onChange(null));
    box.append(reset, error);
    return box;
  }

  function parseInterestTerms(value) {
    return [...new Set(
      String(value || "")
        .split(/[,;\n]+/)
        .map((term) => term.trim().toLocaleLowerCase())
        .filter(Boolean),
    )];
  }

  function matchingMotifs(motifs, value) {
    const terms = parseInterestTerms(value);
    if (!terms.length) return [];
    const matched = new Set();
    terms.forEach((term) => {
      const records = motifs.map((motif) => {
          const fields = [
            motif?.name,
            motif?.motif_id,
            motif?.prefix,
            motif?.name && motif?.motif_id
              ? `${motif.name} (${motif.motif_id})`
              : "",
          ].map((field) => String(field || "").toLocaleLowerCase());
          return { motif, fields };
        }),
        exact = records.filter((record) => record.fields.includes(term)),
        selected = exact.length
          ? exact
          : records.filter((record) =>
            record.fields.some((field) => field.includes(term)),
          );
      selected.forEach((record) => matched.add(record.motif));
    });
    return motifs.filter((motif) => matched.has(motif));
  }

  function hexToRgb(hex) {
    const value = String(hex || "").replace("#", "");
    if (!/^[0-9a-f]{6}$/i.test(value)) return [128, 128, 128];
    return [0, 2, 4].map((offset) => parseInt(value.slice(offset, offset + 2), 16));
  }

  function rgbToHex(rgb) {
    return `#${rgb
      .map((channel) =>
        Math.max(0, Math.min(255, Math.round(channel)))
          .toString(16)
          .padStart(2, "0"),
      )
      .join("")}`;
  }

  function interpolateColor(from, to, fraction) {
    const start = hexToRgb(from),
      end = hexToRgb(to),
      t = Math.max(0, Math.min(1, number(fraction)));
    return rgbToHex(start.map((value, index) => value + (end[index] - value) * t));
  }

  function rankColor(motif, mode, domain, colors = {}) {
    const effect = number(motif?.effect ?? motif?.change),
      maximum = mode === "significance"
        ? Math.max(1e-12, Math.abs(number(domain?.maxAbs, 1)))
        : Math.max(1e-12, number(domain?.max, 1)),
      fraction = mode === "significance"
        ? Math.min(1, Math.abs(effect) / maximum)
        : Math.min(1, negLog10P(motif) / maximum),
      neutral = colors.neutralCenter || "#f8fafc";
    return interpolateColor(
      neutral,
      effect < 0 ? colors.second || "#2563eb" : colors.first || "#dc2626",
      fraction,
    );
  }

  return {
    numeric,
    passesDisplayFilter,
    validRange,
    autoDomain,
    axisEditor,
    matchingMotifs,
    negLog10P,
    oppositeMetric,
    parseInterestTerms,
    rankColor,
    rankMetric,
    rankMotifs,
  };
});
