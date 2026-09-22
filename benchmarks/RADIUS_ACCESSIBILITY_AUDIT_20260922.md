# Wider-background gains: accessibility and uniform-sequence controls

The radius-500 result needs a broader-accessibility sensitivity analysis.
Matching motif score and central 201-position counts did not balance counts
101–500 base pairs away from the motif. The parametric model nevertheless
outperforms a uniform-sequence spatial control on the three promising
zinc-finger tasks. These findings justify further training-only development,
not promotion or a claim that confounding has been eliminated.

## Frozen scope

All 3,949 matched pairs and all 21 tasks from the radius comparison were kept.
There was no rematching, window selection or refitting. The audit evaluated
five predetermined count features on training chromosomes 1–15: the central
201 positions, absolute offsets 50–100, 50–200, 50–500, and 101–500. The last
feature does not overlap the original central-count matching covariate.

Balance and score associations use log1p counts. Both increasing-count and
decreasing-count label discrimination are reported without choosing the
better direction. The uniform-sequence comparator sets both log-bias tracks
to zero, retaining the same radius 500, half-count prior, fixed deviance
dispersion and conventional detector geometry. It uses the same sites and
zero-read score-tie rule as the parametric arm.

## Broad-count imbalance

Standardized mean difference (SMD) and area under the receiver operating
characteristic curve (AUROC) using **decreasing** counts at offsets 101–500:

| Task | Broad-count SMD | Decreasing-count AUROC | Positive/negative median counts |
| --- | ---: | ---: | ---: |
| HepG2 CTCF | −0.745840 | 0.713521 | 11 / 30 |
| K562 CTCF | −0.522079 | 0.658061 | 8 / 20 |
| HepG2 ZNF362 | −0.069244 | 0.528175 | 22 / 24 |
| HepG2 ZNF384 | −0.210239 | 0.572221 | 15 / 22 |
| K562 ZNF384 | −0.224559 | 0.572847 | 18 / 25 |

CCCTC-binding factor (CTCF) and both ZNF384 tasks exceed the prespecified
absolute SMD balance threshold 0.1 on this new covariate. Positive sites have
fewer distant counts in these tasks despite matched central counts. This can
reflect differences in the broader accessibility profile. It is not proof
that the observed footprint-score gains are entirely confounded, and it does
not establish bias-model overcorrection.

## Uniform-sequence control

| Task | Uniform-background AUROC | Parametric radius-500 AUROC | Increment from sequence correction |
| --- | ---: | ---: | ---: |
| HepG2 ZNF362 | 0.569305 | 0.608994 | +0.039689 |
| HepG2 ZNF384 | 0.581866 | 0.638541 | +0.056675 |
| K562 ZNF384 | 0.568456 | 0.615442 | +0.046986 |
| HepG2 CTCF | 0.580581 | 0.600334 | +0.019753 |
| K562 CTCF | 0.647347 | 0.652347 | +0.005000 |

The uniform comparator fails the numerical canonical-protection sign check
for HepG2 ZNF384 and both CTCF tasks, despite some discrimination. This again
separates occupancy ranking from a footprint interpretation. Uniform aggregate
curves are saved, but this audit does not claim a new visual-review or
uncertainty qualification. The earlier six-arm aggregate review remains valid.

All 21 tasks, including underpowered cases, remain in the output. HepG2 HMGA1
has just two matched pairs and must not be interpreted from its apparently
extreme count association. The full output includes 105 feature associations,
735 feature/score correlations and 21 uniform-comparator metric rows.

## Verification and provenance

Seven relevant tests passed in 2.18 seconds, covering exact feature boundaries,
nonoverlap with central counts, both discrimination directions, constant
features, and the previously tested complete comparison and radius operators.
The supervised audit exited zero in 19.04 seconds. Separate artifact
verification passed. Independent recomputation from saved features and scores
matched all 105 SMD/AUROC associations, all 21 AUROC/area under the
precision-recall curve (AUPRC) metric rows and all 735 correlations, using
exact joined coordinate/label identities and round-trip score parsing.
This is diagnostic execution time, not a full-workflow resource measurement.

Outputs are under `results/radius_accessibility_audit_20260922_v1/`;
launch/exit evidence is under
`results/radius_accessibility_jobs_20260922/execute_attempt1/`.

SHA-256:

- Frozen configuration: `4a9a2dc67c053f87455adb2d183f1ae26be94f0aad042ac63d9dbc2212403b16`
- Audit manifest: `08f3f0b2f2fdbe51a11460aeebacddfbf218768b650fa7e270d1fd02ba3987b6`

Next, freeze a secondary analysis that also balances distant counts while
preserving the primary result and reporting excluded pairs and lost power.
Do not select subsets using footprint scores or aggregate outcomes. Only after
that development check should a candidate be considered for frozen tuning and
exact-operator safety evaluation. New transfer labels and the final negative
control remain unopened. Main and production defaults remain unchanged.
