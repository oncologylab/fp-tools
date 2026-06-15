# Code Availability Statement

All source code for `fp-tools` is openly available at
https://github.com/oncologylab/fp-tools under the MIT license, and the package is
installable from source as `fp-tools-bio`. The repository includes:

- the core command-line tools (`atac-correct`, `call-footprints`,
  `match-motifs`, `diff-footprints`, `plot-aggregate`) and their legacy aliases;
- the optional YAML runner (`run-workflow`) and GUI extra;
- the opt-in scientific modules described in this manuscript, including
  replicate-aware `diff-footprints`, motif-discovery preparation,
  and pseudobulk-fragment utilities;
- the unit-test suite and the continuous-integration workflow used for the
  deterministic golden CLI regressions;
- the benchmark scripts (`benchmarks/scripts/`) and paper figure generators
  (`manuscript/scripts/`);
- `environment.yml`, `Dockerfile`, `Makefile`, `LICENSE`, `CITATION.cff`, and
  `.zenodo.json` for reproducible review and release preparation.

The exact version used for this manuscript is recorded by the repository git tag
and commit hash. Paper reproduction instructions are maintained in `docs/reproduce-paper.md`. No proprietary or closed-source components are required to run the
described workflows.
