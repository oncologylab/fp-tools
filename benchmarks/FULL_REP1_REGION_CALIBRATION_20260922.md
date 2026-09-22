# Higher-count replicate-1 calibration

Status: frozen and launched; scientific results and independent verification
are pending. Do not infer completion from this report or a profile file.

The study uses the audited full filtered replicate-1 alignment at exactly the
same 4,200 label-free chromosome 1–15 sites as the sparse-count calibration.
All 21 tasks remain present. Sequence propensities, strand orientation,
geometry (7, 12, 40), exclusion rules, bootstrap weights and paired likelihood
definition are unchanged. No model is fitted, no threshold is changed, and no
new holdout labels or final naked-DNA control are opened.

This is examined development evidence. Full-input deduplication, mitochondrial
removal and read-level blacklist filtering differ from the first-ten-million
input-pair preparation. It is not a controlled depth-only comparison and does
not qualify independent safety. See `REP1_DEPTH_INPUT_AUDIT_20260922.md`.

## Frozen execution

Script: `benchmarks/scripts/calibrate_full_rep1_regions.py`.
It reuses the verified strand-aware extraction, regional analysis and
independent probability/convolution verifier without changing their source.
Verification rejects absent extraction; partial results are preserved.

Ignored result directory:
`benchmarks/results/full_rep1_region_calibration_20260922_v1/`.

- Study configuration SHA-256:
  `33d716bbaf973027be456da09732decd73672e642995d4dc8123220616b2f039`.
- Job specification SHA-256:
  `4ce505c140baf5bc95703b4b7ed0515331f6762062912208405ff9bd32544617`.
- Supervisor launch: process 958690; child 958691, start ticks 120806942,
  boot `664b7c80-18d6-4147-bdf4-2eae53f766e6`.
- Launch receipt and logs:
  `benchmarks/results/full_rep1_region_calibration_jobs_20260922/execute_attempt1/`.

The job freezes input/code hashes before extraction and freezes the derived
calibration manifest before calculating summaries. It then executes the
independent verifier. Successful process exit must be accompanied by verified
profile and calibration manifests and scientific review.

## Validation so far

Ten focused tests passed in 3.87 seconds: new orchestration checks plus existing
regional calibration, conditional protection and score-ceiling tests.
The orchestration checks cover existing-output preservation, rejection of
changed policy, and refusal to extract missing profiles in verify mode.
The previous checksum-valid sparse-count result was verified without rerunning
its experiment. No full package regression or resource qualification is claimed.

## Resume

Inspect the recorded live process identity and terminal receipt before any
action. Do not restart because a status observation timed out. When terminal,
inspect `stdout.log`, `exit.json`, `profiles_manifest.json` and
`calibration/calibration_manifest.json`; run the script with `--stage verify`
against the same output directory. Report all tasks, calibration discrepancies,
strand/shoulder effects and independent-null assumptions. No promotion follows
from this development calibration alone.
