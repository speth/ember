#!/usr/bin/env python
"""
Plot error-vs-N grid-convergence curves from the JSON files produced by
run_convergence.py (Task 2.2's convergence study; see spec §6.4).

For each requested case, loads every ``results/<case>_<scheme>_rung*.json``
file (excluding any ``--damp-const`` trials, which are a separate one-off
comparison, not part of the main convergence ladder), takes the finest
(largest-N) ``secondOrderLimited`` run as the reference solution, and plots
normalized error vs. N (log-log) for both schemes on one figure per
(case, metric) pair -- metrics: ``consumption_speed`` and ``peak_T``.

Usage:
    pixi run python test/convergence/plot_convergence.py \\
        [--cases strained twin cylindrical] \\
        [--resultsdir test/convergence/results] \\
        [--outdir test/convergence/results/plots]

Requires at least one secondOrderLimited run per requested case to serve as
the reference; requires at least two runs total (any scheme) to plot a
curve. Cases/metrics without enough data are skipped with a warning, not a
hard failure, since this script may be run against a partial results
directory during Task 2.2.
"""

import argparse
import glob
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = ('strained', 'twin', 'cylindrical')
SCHEMES = ('firstOrderUpwind', 'secondOrderLimited')
METRICS = ('consumption_speed', 'peak_T')

# Categorical slots 1 (blue) and 2 (aqua) from the project's validated
# palette, assigned in fixed order (new default scheme first).
COLORS = {
    'secondOrderLimited': '#2a78d6',
    'firstOrderUpwind': '#1baf7a',
}
MARKERS = {
    'secondOrderLimited': 'o',
    'firstOrderUpwind': 's',
}


def load_runs(resultsdir, case):
    """Load all non-damp-const runs for a case, grouped by scheme."""
    runs = {scheme: [] for scheme in SCHEMES}
    for path in sorted(glob.glob(os.path.join(resultsdir, '%s_*_rung*.json' % case))):
        with open(path) as f:
            data = json.load(f)
        if data.get('case') != case:
            continue
        if data.get('damp_const') is not None:
            continue  # separate trial, not part of the main ladder
        scheme = data.get('scheme')
        if scheme not in runs:
            continue
        runs[scheme].append(data)
    for scheme in runs:
        runs[scheme].sort(key=lambda d: d['N'])
    return runs


def plot_case_metric(case, metric, runs, outdir):
    ref_candidates = runs.get('secondOrderLimited', [])
    if not ref_candidates:
        print('  [skip] %s / %s: no secondOrderLimited runs found for reference' %
              (case, metric))
        return
    reference_run = ref_candidates[-1]  # largest N
    ref_value = reference_run['scalars'].get(metric)
    if ref_value is None:
        print('  [skip] %s / %s: reference run has null %s' % (case, metric, metric))
        return

    fig, ax = plt.subplots(figsize=(6, 4.5))
    any_plotted = False
    for scheme in SCHEMES:
        Ns = []
        errs = []
        for run in runs.get(scheme, []):
            value = run['scalars'].get(metric)
            if value is None:
                continue
            # Skip the reference run itself (zero error, breaks log scale).
            if run is reference_run:
                continue
            err = abs(value - ref_value) / abs(ref_value)
            if err <= 0:
                continue
            Ns.append(run['N'])
            errs.append(err)
        if len(Ns) < 2:
            continue
        any_plotted = True
        ax.loglog(Ns, errs, marker=MARKERS[scheme], markersize=8, linewidth=2,
                   color=COLORS[scheme], label=scheme)

    if not any_plotted:
        print('  [skip] %s / %s: fewer than 2 comparable points for any scheme' %
              (case, metric))
        plt.close(fig)
        return

    ax.set_xlabel('N (grid points)')
    ax.set_ylabel('relative error vs. finest secondOrderLimited run')
    ax.set_title('%s: %s convergence' % (case, metric))
    ax.grid(True, which='both', linewidth=0.5, color='#c3c2b7', alpha=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, '%s_%s.png' % (case, metric))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print('  wrote %s' % out_path)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--resultsdir', default='test/convergence/results')
    parser.add_argument('--outdir', default=None,
                         help='Default: <resultsdir>/plots')
    args = parser.parse_args()
    outdir = args.outdir or os.path.join(args.resultsdir, 'plots')

    for case in args.cases:
        print('=== %s ===' % case)
        runs = load_runs(args.resultsdir, case)
        n_found = sum(len(v) for v in runs.values())
        if n_found == 0:
            print('  no result JSON files found for case %r in %s' % (case, args.resultsdir))
            continue
        for metric in METRICS:
            plot_case_metric(case, metric, runs, outdir)


if __name__ == '__main__':
    main()
