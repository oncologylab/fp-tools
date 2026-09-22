# All-site DWM control construction exposes region-layout dependence

Constructed the expanded replicate-2 control reference, but it is **not
qualified for detector rescreening**. Observed cuts match the corrected
extraction at every site; expected-profile proportions differ at three of
269 previously positive windows. Preserve both references and investigate
the differences before changing any safety result.

## Frozen construction

Used all 4,200 previously recorded, label-free control sites, represented by
4,192 unique motif intervals. Reused the recorded naked-DNA replicate-2 BAM,
index and hg38 reference, and the existing 25-position dinucleotide weight
matrix (DWM) model adjacent to the original CTCF-only tracks. Inputs, source
and native calculation modules, command and runtime were frozen before
execution. The model's SHA-256 is
`1d8f503cb08144fff26db27d9697a51527dbf48da53fa525f925f1e6354e194d`.

The command used four cores, output-region extension 100, correction window
100, sequence flank 12, read shifts +4/−5 and the existing model via
`--bias-pkl`. It retained the historical peak file while expanding only the
output-region list. Library normalization was disabled; the exact safety
scorer rescales expected profiles to each site's observed total. Overlap
qualification therefore compares expected proportions and observed cut
profiles, rather than requiring the old global scaling factor of
133.94411851375605. The command and settings are recorded in the freeze.

No bias model was refitted. All strand bias and background count arrays in
the produced model exactly match the original arrays. The reused model's
association with the old track is tested through overlap agreement, not
assumed solely from its neighboring filename.

## Coverage and overlap results

- All 4,200 extracted observed profiles exactly match the corrected BAM
  reference, and all sites retain valid coordinates and fresh extraction.
- Positive expected coverage increases from 269 to 2,892 windows. Every one
  of the 21 tasks now has 124–156 positive-expectation windows out of 200.
  Zero expectation still requires an explicit support exclusion.
- Expected proportions match at 266/269 old positive windows under the
  prespecified relative tolerance 1e-5 and absolute tolerance 1e-7.
- Three windows fail. No threshold or detector was changed or evaluated.

| Window | Maximum expected-proportion difference | Old output region | Expanded output region |
| --- | ---: | --- | --- |
| HepG2 CCCTC-binding factor (CTCF), chr10:28991925–28991940 | 0.028374 | 28991825–28992040 | 28991500–28992040 |
| K562 CTCF, chr9:98416493–98416508 | 0.006280 | 98416393–98416608 | 98416281–98416608 |
| K562 FOXA1, chr9:98416381–98416389 | 0.056100 | 98416393–98416608 | 98416281–98416608 |

Reconstructed region processing from the frozen inputs: extend by 162,
merge, then trim 62 bases from each boundary, matching the native output-
region preparation. The two merged regions change start phase modulo ten
from 5 to 0 and from 3 to 1. Native correction fits windows at ten-base
steps anchored to the region. This is consistent with region-layout-induced
expectation drift; a historical-region-only rerun is still required to isolate
that cause. The FOXA1 window was only partly inside the old output domain.

The failed overlap condition is retained in `qualified_for_rescreen: false`.
Expanded coverage alone does not justify replacing the old safety reference,
dropping the changed windows, or revising the earlier CTCF safety failures.

## Verification and provenance

- Eight focused tests passed in 1.91 seconds: label-free region union,
  coordinate validation, absent/zero/positive track distinctions, data-root
  handling and acceptance of global scaling while rejecting shape drift.
- Durable construction exited zero in 41.368155 seconds. Model loading and
  application are diagnostic; this is not the complete fit-and-apply resource
  benchmark or a peak-memory measurement.
- Separate input/code/runtime/output verification passed.
- Both new tracks were extracted and motif-oriented at all recorded sites;
  fresh observed counts agree exactly and all overlap failures are retained.
- Independent interval-based extraction checks all 8,400 observed/expected
  profiles and the 21 task summaries separately from the vector track reader.

Ignored directories under `benchmarks/results/`:
`all_site_control_dwm_20260922_v1/`,
`all_site_control_dwm_review_20260922_v1/`, and
`all_site_control_dwm_jobs_20260922/`.

SHA-256:

```text
configuration.freeze.json 3155703a29b2e1ba036425dc583e3f073c0b58197c83b2df1b3f9fda50e2279f
tracks_manifest.json      0a33db3419a5d8f6d1c78168fe8737edc22f99507bc83da79d4dff28754a45e3
review_manifest.json      c737ac5b13ac703de7fe1409b000cb39c3d93ddb14122dac983f9e8ecce9181b
```

Next freeze a historical-region-only construction with the same current
implementation, recorded BAM, bias model and normalization choice. Compare
that against both old and expanded tracks to isolate output-region effects.
Keep all models, thresholds and existing results unchanged. Main, released
packages, final controls and untouched promotion labels remain untouched.
