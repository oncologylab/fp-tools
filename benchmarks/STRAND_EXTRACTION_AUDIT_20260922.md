# Research strand-extraction correction

A numerical-parity review found two defects in the research-only strand
profile builder. Querying disjoint motif windows can fetch the same long
alignment more than once; the old implementation assigned its cut to all
matching windows on every fetch. Conversely, a soft-clipped cut can lie in a
window while its aligned span does not overlap that query, so it is missed.
Neither defect depends on occupancy labels or a particular fitted model.

The extractor now visits each requested chromosome once, filters alignments,
computes the established cut coordinate and assigns it once to every containing
motif window. Overlapping sites still each receive the cut. A full chromosome
scan is needed to include arbitrary soft clipping without assuming an unverified
maximum clip length. Its workflow resource cost remains to be measured.
Nonconsecutive DataFrame indices are normalized without changing site order.
Production source and main were not modified.

## Frozen training-only impact audit

The audit used only chromosome 1 coordinates from the existing K562/HepG2
training inputs, with replicate-1 10M seed-2026 BAMs. It read no ChIP-label
columns and performed no detector scoring, model fitting or selection. The
old extraction definitions came from immutable commit `4df283a`. The corrected
extractor was compared with an independent oracle that scans chromosome 1,
sorts strand-specific cut coordinates, and counts separately for each site.

Input BAM/index hashes, source and label-free subset tables, legacy/current
code, runtime, shift +4/-5, mapping-quality threshold 30 and duplicate exclusion
were frozen before execution. The initial v1 freeze was superseded, without
running it, after a fixture exposed the independent oracle's empty-strand dtype
error. Its files remain preserved. The final v2 configuration SHA-256 is
`a500e8191706f17f468fcb986fd484750a928e9b8e17b67c2c4718e21feb14b4`.

| Cell | Sites | Changed sites | Legacy cuts | Corrected cuts |
| --- | ---: | ---: | ---: | ---: |
| HepG2 | 2,381 | 35 | 31,715 | 31,542 |
| K562 | 1,943 | 15 | 40,304 | 40,259 |

These totals count motif-window assignments, not unique library fragments.
The corrected plus/minus matrices and validity masks match the independent
oracle exactly at every site. All observed differences in this slice remove
duplicate assignments; the soft-clipping defect was reproduced by synthetic
fixtures but did not add counts in this slice.

For CTCF, 6/251 HepG2 sites change (4,103 to 4,067 assignments), as do 3/249
K562 sites (6,630 to 6,616). Other affected factors include MEF2A, MEF2D, MYC,
TBP, TCF7L2 and ZNF384. This audit is confined to one training chromosome;
it does not quantify effects on the complete stability study, independent
negative controls, discrimination metrics or protection shapes.

## Evidence disposition and next work

Existing artifacts are checksum-valid historical computations, but checksum
verification does not prove extraction correctness. Preserve them and their
reported numerical results; do not use them alone to certify a new candidate.
The 96-run scores, exact-operator safety estimates and low-depth CTCF diagnoses
require qualification after the extraction correction. Their numerical change
is unknown until independently frozen reassessment. No promotion was valid
before this discovery, and none is claimed now.

The next corrected experiment must use a new output root, frozen coordinates,
model coefficients, dispersions, geometry and thresholds, with no retuning on
examined test labels. It must identify which artifacts used the affected
extractor and verify same-method parity where unaffected. Candidate development
continues only on training and tuning chromosomes; SK-N-SH/GM23338 labels and
final naked-DNA replicate 3 remain locked. A new runtime benchmark must measure
the corrected reader, not the flawed faster extraction path.

Artifacts: `benchmarks/results/strand_extraction_audit_20260922_v2/`.
The result manifest SHA-256 is
`5fd55dcdd330a5f51e88b81855291a6c1ac1e5d3ddceb680a935fb542aa0cf78`.
It records the per-factor table and old/corrected strand arrays. Forty focused
extraction, strand-profile, factorial and stability-extension tests passed in
26.13 seconds. Before the fix, the boundary fixture had three failing tests
and one passing overlapping-site control.
