# Code Availability Statement

All source code for `fp-tools` is openly available at
https://github.com/oncologylab/fp-tools under the MIT license, and the package is
installable from source as `fp-tools-bio`. The repository includes:

- the core command-line tools (`atac-correct`, `score-footprints`,
  `detect-tf-binding`, `plot-aggregate`) and their legacy aliases;
- the optional YAML runner (`fp-tools-run`) and Streamlit GUI (`fp-tools-gui`);
- the opt-in scientific modules described in this manuscript, including
  replicate-aware `detect-tf-binding`, BINDetect report compatibility support,
  and `fp-tools-decompose-competition`;
- the unit-test suite and the continuous-integration workflow used for the
  deterministic golden CLI regressions;
- the benchmark scripts (`benchmarks/scripts/`) and paper figure generators
  (`paper/scripts/`).

The exact version used for this manuscript is recorded by the repository git tag
and commit hash. No proprietary or closed-source components are required to run the
described workflows.
