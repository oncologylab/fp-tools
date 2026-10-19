"""Independent read-only verifier for the calibration-round evidence package.

Recomputes the headline numbers from packaged derived data with separately
written code: recalibration by generic optimisation, tail probabilities by
explicit convolution, AUROC by ranks and average precision by sorting. It
never writes inside the package and never repairs or regenerates outputs.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import binom, rankdata

THRESHOLD = 1.3010299957
CHROMS = [f'chr{i}' for i in range(1, 16)]
FAILURES = []


def check(condition, message):
    if not condition:
        FAILURES.append(message)
        print('FAIL', message, flush=True)


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', low_memory=False)


def checksums(root):
    listed = {}
    for line in (root / 'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ', 1)
        listed[name] = digest
    present = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name != 'SHA256SUMS'}
    check(set(listed) == present, f'manifest mismatch: missing {sorted(set(listed) - present)[:5]} extra {sorted(present - set(listed))[:5]}')
    bad = [n for n, d in listed.items() if (root / n).exists() and sha256((root / n).read_bytes()).hexdigest() != d]
    check(not bad, f'checksum mismatch {bad[:5]}')
    registry = json.loads((root / 'inputs.json').read_text())['permitted']
    for rel, entry in registry.items():
        check(sha256((root / 'data' / rel).read_bytes()).hexdigest() == entry['sha256'], 'input identity ' + rel)
    plan = json.loads((root / 'PLAN.freeze.json').read_text())
    for rel, digest in plan['code'].items():
        check(sha256((root / 'code' / Path(rel).name).read_bytes()).hexdigest() == digest, 'frozen code identity ' + rel)
    result = json.loads((root / 'results/result.json').read_text())
    check(result['plan_sha256'] == sha256((root / 'PLAN.freeze.json').read_bytes()).hexdigest(), 'plan hash')
    return len(listed)


def nll(theta, c, n, x, slope):
    a, b = (theta[0], theta[1]) if slope else (theta[0], 1.)
    eta = a + b * x
    return -np.sum(c * eta - n * np.logaddexp(0, eta))


def recalibration(root, controls):
    fits = json.loads((root / 'results/models/recalibration.json').read_text())['records']
    train = controls[controls.eligible & controls.LOG21_plus_valid & controls.LOG21_minus_valid
                     & controls.deduplicated_fit_representative]
    coefficients = {}
    for rec in fits:
        rows = train if rec['fold'].startswith('deployable') else train[train.TFBS_chr != rec['fold']]
        c = np.r_[rows.plus_center, rows.minus_center].astype(float)
        n = np.r_[rows.plus_total, rows.minus_total].astype(float)
        x = logit(np.clip(np.r_[rows.LOG21_plus_p, rows.LOG21_minus_p], 1e-6, 1 - 1e-6))
        keep = n > 0
        slope = rec['arm'] == 'LOG21NC_count_tail'
        opt = minimize(nll, [0., 1.] if slope else [0.], args=(c[keep], n[keep], x[keep], slope), method='BFGS',
                       options=dict(gtol=1e-9))
        a, b = (opt.x[0], opt.x[1]) if slope else (opt.x[0], 1.)
        check(abs(a - rec['a']) < 1e-5 and abs(b - rec['b']) < 1e-5, f"recalibration {rec['arm']} {rec['fold']}")
        # Downstream checks use the recorded coefficients once the independent optimum agrees.
        coefficients[(rec['arm'], rec['fold'])] = (rec['a'], rec['b'])
    return coefficients, len(fits)


def probabilities(frame, arm, coefficients):
    if arm == 'uniform_count_tail':
        return np.full(len(frame), 15 / 73), np.full(len(frame), 15 / 73)
    prefix = 'DWM' if arm == 'DWM_count_tail' else 'LOG21'
    plus = frame[prefix + '_plus_p'].to_numpy(float); minus = frame[prefix + '_minus_p'].to_numpy(float)
    if arm in ('LOG21NC_count_tail', 'LOG21NC1_count_tail'):
        ab = np.array([coefficients[(arm, c)] for c in frame.TFBS_chr])
        f = lambda p: expit(ab[:, 0] + ab[:, 1] * logit(np.clip(p, 1e-6, 1 - 1e-6)))  # noqa: E731
        return f(plus), f(minus)
    return plus, minus


def tail(n_plus, n_minus, central, q_plus, q_minus):
    pmf = np.convolve(binom.pmf(np.arange(n_plus + 1), n_plus, q_plus), binom.pmf(np.arange(n_minus + 1), n_minus, q_minus))
    p = pmf[:central + 1].sum()
    return 0. if central >= n_plus + n_minus else min(300., max(0., -np.log10(p)))


def scores(root, coefficients):
    packaged = read(root / 'results/tables/chromatin_site_scores.tsv.gz')
    rng = np.random.default_rng(1)
    checked = calls = 0
    for (library, depth), g in packaged.groupby(['library', 'depth']):
        counts = read(root / f'data/tables/counts.{library}.{depth}.tsv.gz').set_index('site_id')
        nested = read(root / 'data/tables/nested_sites.tsv.gz')
        check(set(g.site_id) == set(nested[nested.cell == library.split('_')[0]].site_id), f'cohort {library} {depth}')
        for arm, h in g.groupby('arm'):
            c = counts.loc[h.site_id].reset_index()
            for col in ['plus_center', 'minus_center', 'plus_total', 'minus_total']:
                check((c[col].to_numpy() == h[col].to_numpy()).all(), f'counts {library} {depth} {arm} {col}')
            qp, qm = probabilities(c, arm, coefficients)
            valid = (c.DWM_plus_valid & c.DWM_minus_valid) if arm == 'DWM_count_tail' else (
                np.ones(len(c), bool) if arm == 'uniform_count_tail' else (c.LOG21_plus_valid & c.LOG21_minus_valid))
            check(np.allclose(np.where(valid, qp, np.nan), h.q_plus, rtol=0, atol=1e-10, equal_nan=True), f'q {library} {depth} {arm}')
            # Exact tail for every call-relevant row and a random 10% of the remainder.
            sample = (h.score.to_numpy() > THRESHOLD - .05) | (rng.random(len(h)) < .1)
            for i in np.flatnonzero(sample & np.asarray(valid, bool)):
                r = c.iloc[i]
                s = round(tail(int(r.plus_total), int(r.minus_total), int(r.plus_center + r.minus_center), qp[i], qm[i]), 10)
                ceil = round(tail(int(r.plus_total), int(r.minus_total), 0, qp[i], qm[i]), 10)
                check(abs(s - h.score.iloc[i]) < 1e-8 and abs(ceil - h.ceiling.iloc[i]) < 1e-8, f'score {library} {depth} {arm} {r.site_id}')
                checked += 1
            expected_call = (h.total > 0) & (h.ceiling > THRESHOLD) & (h.score > THRESHOLD)
            check((expected_call.to_numpy() == h.call.to_numpy()).all(), f'calls {library} {depth} {arm}')
            calls += int(h.call.sum())
    return packaged, checked, calls


def auc_ap(y, s):
    y = np.asarray(y, int); s = np.asarray(s, float); pos = y.sum(); neg = len(y) - pos
    r = rankdata(s)
    auc = (r[y == 1].sum() - pos * (pos + 1) / 2) / (pos * neg)
    order = np.argsort(-s, kind='stable'); ys = y[order]; ss = s[order]
    ends = np.r_[np.flatnonzero(ss[:-1] != ss[1:]), len(y) - 1]
    tp = np.cumsum(ys)[ends]
    return auc, np.sum(np.diff(np.r_[0, tp]) / pos * tp / (ends + 1))


def metrics(root, packaged):
    folds = read(root / 'results/tables/fold_metrics.tsv')
    mine = []
    for (library, depth, tf, arm, chrom), g in packaged.groupby(['library', 'depth', 'tf', 'arm', 'TFBS_chr']):
        if g.chip_label.nunique() == 2 and np.isfinite(g.score).all():
            a, p = auc_ap(g.chip_label, g.score)
            mine.append(dict(library=library, depth=depth, tf=tf, arm=arm, fold=chrom, auroc_v=a, auprc_v=p))
    joined = folds.merge(pd.DataFrame(mine), on=['library', 'depth', 'tf', 'arm', 'fold'], how='inner')
    check(len(joined) == int(folds[folds.arm != 'covariates_only'].status.eq('AVAILABLE').sum()), 'available fold count')
    check(np.allclose(joined.auroc, joined.auroc_v, atol=1e-12) and np.allclose(joined.auprc, joined.auprc_v, atol=1e-12), 'fold metrics')
    summary = read(root / 'results/tables/screen_summary.tsv')
    ledger = read(root / 'data/tables/task_ledger.tsv')
    weights = np.random.default_rng(2026).multinomial(15, np.repeat(1 / 15, 15), size=2000)
    for s in summary.itertuples():
        f = folds[folds.depth == s.depth]
        points, boots = [], []
        for t in ledger.itertuples():
            a = f[(f.cell == t.cell) & (f.tf == t.tf) & (f.arm == s.method)].set_index(['library', 'fold'])
            b = f[(f.cell == t.cell) & (f.tf == t.tf) & (f.arm == s.reference)].set_index(['library', 'fold'])
            d = pd.DataFrame({'roc': a.auroc - b.auroc, 'ap': a.auprc - b.auprc, 'base': b.auprc}).dropna()
            per = d.groupby(level=0).mean()
            if len(per) != 3:
                continue
            points.append((per.roc.mean(), per.ap.mean() / per.base.mean()))
            draws = []
            for lib in per.index:
                delta = d.loc[lib].roc.reindex(CHROMS).to_numpy(); v = np.isfinite(delta)
                den = weights @ v.astype(float)
                draws.append(np.divide(weights @ np.nan_to_num(delta), den, out=np.full(2000, np.nan), where=den > 0))
            boots.append(np.mean(draws, 0))
        draw = np.mean(boots, 0); v = np.isfinite(draw)
        lo, hi = np.quantile(draw[v], [.025, .975])
        check(len(points) == s.available_tasks and s.original_tasks == 21, f'task denominator {s.method} {s.reference}')
        check(np.allclose([np.mean([p[0] for p in points]), np.mean([p[1] for p in points]), lo, hi],
                          [s.mean_auroc_gain, s.mean_relative_auprc_gain, s.bootstrap_low, s.bootstrap_high], atol=1e-12),
              f'screen {s.depth} {s.method} {s.reference}')
    return len(joined), summary


def controls_check(root, controls, coefficients):
    means = read(root / 'results/tables/control_task_means.tsv')
    support = controls.eligible & controls.LOG21_plus_valid & controls.LOG21_minus_valid & controls.DWM_plus_valid & controls.DWM_minus_valid
    for m in means.itertuples():
        g = controls[support & (controls.cell == m.cell) & (controls.tf == m.tf)]
        qp, qm = probabilities(g, m.arm, coefficients)
        o = (g.plus_center + g.minus_center).sum(); e = (g.plus_total * qp + g.minus_total * qm).sum()
        check(abs(o / e - m.ratio) < 1e-10, f'control ratio {m.arm} {m.cell} {m.tf}')
    return means


def decision(root, means, summary):
    screen = read(root / 'results/tables/advancement_screen.tsv').set_index('gate')
    ll = read(root / 'results/tables/control_heldout_loglik.tsv').set_index('arm').heldout_loglik
    disc = means.groupby('arm').discrepant.sum()
    a_pass = ll['LOG21NC_count_tail'] > ll['LOG21_count_tail'] and disc['LOG21NC_count_tail'] <= 3 \
        and disc['LOG21NC_count_tail'] < disc['LOG21_count_tail']
    check(screen.loc['A_calibration', 'status'] == ('PASS' if a_pass else 'FAIL'), 'gate A')
    m = summary[(summary.depth == 'full') & (summary.method == 'LOG21NC_count_tail') & (summary.reference == 'DWM_count_tail')].iloc[0]
    b_pass = m.mean_auroc_gain >= .03 and m.mean_relative_auprc_gain >= .10 and m.bootstrap_low > 0
    check(screen.loc['B_accuracy', 'status'] == ('PASS' if b_pass else 'FAIL'), 'gate B')
    check(not a_pass and 'rejected' in screen.loc['DECISION', 'evidence'], 'decision')
    return dict(A=a_pass, B=b_pass, discrepant=disc.to_dict(), loglik=ll.to_dict())


def ledgers(root):
    t = read(root / 'results/tables/all_task_ledger.tsv')
    check(t[['cell', 'tf']].drop_duplicates().shape[0] == 21, '21 tasks')
    h = t[(t.cell == 'HepG2') & (t.tf == 'HMGA1')]
    check(h.chromatin_status.eq('UNAVAILABLE').all(), 'HMGA1 unavailable not imputed')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('package', type=Path)
    root = parser.parse_args().package.resolve()
    files = checksums(root)
    controls = read(root / 'data/tables/control_counts.full.tsv.gz')
    coefficients, fits = recalibration(root, controls)
    means = controls_check(root, controls, coefficients)
    packaged, checked, calls = scores(root, coefficients)
    folds, summary = metrics(root, packaged)
    gates = decision(root, means, summary)
    ledgers(root)
    head = summary[(summary.depth == 'full') & (summary.reference == 'DWM_count_tail')][['method', 'mean_auroc_gain', 'mean_relative_auprc_gain']]
    print(json.dumps(dict(files=files, recalibration_fits=fits, exact_tail_rows=checked, calls=calls, fold_metrics=folds,
                          gates=gates, headline=head.to_dict('records'), failures=len(FAILURES)), indent=1, default=float))
    sys.exit(1 if FAILURES else 0)


if __name__ == '__main__':
    main()
