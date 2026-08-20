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
    const total = Math.max(2, Math.floor(number(limit, 20))),
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
    return {
      negative: negative.slice(0, negativeCount),
      positive: positive.slice(0, positiveCount),
    };
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
    matchingMotifs,
    negLog10P,
    oppositeMetric,
    parseInterestTerms,
    rankColor,
    rankMetric,
    rankMotifs,
  };
});
