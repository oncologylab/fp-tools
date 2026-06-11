# Future Plans

This file records fp-tools modules that exist as prototypes or development-branch
utilities but are not part of the first paper/release scope. The code is kept in
the repository for continuity, testing, and future validation, but these features
should not be presented as first-version supported workflows until their public
benchmarks, examples, and documentation are complete.

## Deferred First-Version Features

### Supervised TFBS Model

Prototype commands:

- `fp-tools-build-tfbs-features`
- `fp-tools-train-tfbs-model`
- `fp-tools-predict-tfbs`

Planned work:

- Finalize public training and held-out evaluation datasets.
- Define stable feature schemas and model-card style reporting.
- Compare directly against motif-only, accessibility-only, footprint-only, and external supervised/deep-learning baselines on the same loci.
- Add calibration and transferability analyses before recommending routine use.

### Motif-Free Candidate Generation

Prototype command:

- `fp-tools-generate-candidates`

Planned work:

- Validate signal-only candidate recovery on motif-removal and ChIP/CUT&RUN labels.
- Tune defaults for score thresholds, local maxima, and per-region limits.
- Document when motif-free candidates are biologically useful versus exploratory.

### Motif-Relaxed Reranking

Prototype command:

- `fp-tools-rerank-candidates`

Planned work:

- Benchmark relaxed-PWM candidates against strict motif baselines.
- Define stable input columns for motif-family, candidate, and optional model scores.
- Add examples that do not depend on unpublished or large local benchmark files.

### Variant Scoring Scaffold

Prototype command:

- `fp-tools-score-variants`

Planned work:

- Validate against allele-specific accessibility, caQTL, or reporter datasets.
- Separate deterministic sequence/motif annotations from predictive claims.
- Add variant-level examples with small public fixtures and clear limitations.

### Competition-Aware Footprint Decomposition

Prototype command:

- `fp-tools-decompose-competition`

Planned work:

- Validate TF-scale versus nucleosome-scale components against matched labels.
- Establish robust defaults for scale bands and competition summaries.
- Add figure examples and interpretation guidance before presenting this as a supported biological analysis.

## Documentation Policy

Until the work above is complete, the README, manuscript, and feature-comparison
tables should describe these items as future plans or development prototypes, not
as first-version supported capabilities.
