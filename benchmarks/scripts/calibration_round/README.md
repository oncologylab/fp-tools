# Calibration follow-up (October 2026)

Research-only scripts. Plan: `benchmarks/CALIBRATION_ROUND_PLAN_20261019.md`.

```bash
PY=../fp-tools/.venv/bin/python
OUT=benchmarks/results/calibration_round_20261019_v1
$PY benchmarks/scripts/calibration_round/inputs.py --outdir $OUT   # freeze
$PY benchmarks/scripts/calibration_round/run.py --outdir $OUT      # execute once
```

`inputs.py` writes `inputs.json` (permitted inputs with SHA-256, protected
identities) and `PLAN.freeze.json`. `run.py` refuses to start if the policy or
scientific code changed after the freeze, and never overwrites outputs.
