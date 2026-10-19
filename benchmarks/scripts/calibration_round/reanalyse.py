"""Re-run the frozen calibration-round calculation from a packaged copy of its inputs.

Uses the frozen run.py stage functions unchanged, reading `data/` inside an
extracted evidence package (no raw alignments), writing only to a new
directory, and comparing every regenerated table with the packaged result.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run  # noqa: E402
from protocol import ABLATION, COMPARATOR, DEPTHS, PRIMARY, REFERENCE  # noqa: E402


class PackageAccess:
    def __init__(self, package):
        self.package = Path(package)
        self.registry = json.loads((self.package / 'inputs.json').read_text())['permitted']

    def path(self, relative):
        path = self.package / 'data' / relative
        h = sha256(path.read_bytes()).hexdigest()
        if h != self.registry[relative]['sha256']:
            raise ValueError('packaged input differs from frozen identity: ' + relative)
        return path


def execute(package, out, workers):
    access = PackageAccess(package)
    out.mkdir(parents=True, exist_ok=False)
    for sub in ['tables', 'models']:
        (out / sub).mkdir()
    ledger = run.read(access.path('tables/task_ledger.tsv'))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        controls, coefficients = run.control_stage(access, out, pool)
        scores, _ = run.chromatin_stage(access, out, pool, coefficients)
    folds = run.fold_metrics(scores, run.read(access.path('tables/covariate_baseline_all_libraries.tsv.gz')), ledger)
    run.write(folds, out / 'tables/fold_metrics.tsv')
    run.write(run.arm_means(folds, ledger), out / 'tables/arm_task_metrics.tsv')
    deltas, summaries = [], []
    for depth in DEPTHS:
        for method, reference in [(PRIMARY, COMPARATOR), (ABLATION, COMPARATOR), (REFERENCE, COMPARATOR),
                                  (PRIMARY, REFERENCE), (ABLATION, REFERENCE), (PRIMARY, ABLATION),
                                  ('uniform_count_tail', COMPARATOR), (PRIMARY, 'covariates_only')]:
            rows, summary = run.screen_pairs(folds, ledger, depth, method, reference)
            deltas.extend(rows); summaries.append(summary)
    deltas = pd.DataFrame(deltas)
    run.write(deltas, out / 'tables/task_deltas.tsv')
    run.write(pd.DataFrame(summaries), out / 'tables/screen_summary.tsv')
    shape = run.shape_stage(scores)
    run.write(shape, out / 'tables/regional_shapes.tsv')
    run.write(run.information_stage(scores), out / 'tables/information.tsv')
    run.write(run.power_stage(access, controls, coefficients), out / 'tables/synthetic_power.tsv')
    means = run.read(out / 'tables/control_task_means.tsv'); calls = run.read(out / 'tables/control_calls.tsv')
    decision = run.advancement(out, means, calls, summaries, deltas, shape, ledger)
    mismatches = []
    for path in sorted((out / 'tables').glob('*')):
        packaged = run.read(Path(package) / 'results/tables' / path.name)
        mine = run.read(path)
        try:
            pd.testing.assert_frame_equal(mine, packaged, check_exact=False, rtol=1e-12, atol=1e-12)
        except AssertionError as error:
            mismatches.append(f'{path.name}: {str(error)[:200]}')
    print(json.dumps(dict(decision=decision, tables_compared=len(list((out / 'tables').glob('*'))),
                          mismatches=mismatches), indent=2))
    return 1 if mismatches else 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('package', type=Path)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    a = parser.parse_args()
    sys.exit(execute(a.package.resolve(), a.outdir.resolve(), a.workers))
