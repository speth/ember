#!/usr/bin/env python
"""
Plot grid-convergence curves from the JSON files produced by
run_convergence.py (the errTol-ladder study; see spec §6.4 and
docs/superpowers/specs/2026-07-05-error-based-grid-adaptation-design.md §4).

For each requested case, loads every ``results/<case>_<scheme>_rung*.json``
file (excluding any ``--damp-const`` trials, which are a separate one-off
comparison, not part of the main convergence ladder), takes the finest
(largest-N) ``secondOrderLimited`` run as the reference solution, and
produces, per (case, metric):

- error vs. N (log-log) -- the original Task 2.2 view.
- error vs. errTol (log-log) -- the parity view (Task G6): same-errTol
  points from both schemes should track within a bounded band.

It also writes, per case, an sol/fou N-ratio table at matched errTol (rung
index), and a combined errTol -> measured-QoI-error markdown table across
all cases/schemes for use in user docs.

Note: this script has always kept its metadata in scheme-and-tolerance-
agnostic terms (``grid_tolerances`` is a free-form dict; runs are grouped
purely by case/scheme/rung) -- there is no ``vtol``/``dvtol`` naming to
migrate here. Task G6 only adds the errTol-specific views above.

Usage:
    pixi run python test/convergence/plot_convergence.py \\
        [--cases strained twin cylindrical] \\
        [--resultsdir test/convergence/results] \\
        [--outdir test/convergence/results/plots]

Requires at least one secondOrderLimited run per requested case to serve as
the reference; requires at least two runs total (any scheme) to plot a
curve. Cases/metrics without enough data are skipped with a warning, not a
hard failure, since this script may be run against a partial results
directory.
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


def reference_run_and_value(case, metric, runs):
    """Pick the finest (largest-N) secondOrderLimited run as the reference
    solution and return (reference_run, ref_value), or (None, None) if
    unavailable. Shared by all plot/table functions below so every view
    uses the same reference convention (matches the Task 2.2 precedent in
    the spec addendum §P2.2 / calibration-notes.md)."""
    ref_candidates = runs.get('secondOrderLimited', [])
    if not ref_candidates:
        print('  [skip] %s / %s: no secondOrderLimited runs found for reference' %
              (case, metric))
        return None, None
    reference_run = ref_candidates[-1]  # largest N
    ref_value = reference_run['scalars'].get(metric)
    if ref_value is None:
        print('  [skip] %s / %s: reference run has null %s' % (case, metric, metric))
        return None, None
    return reference_run, ref_value


def scheme_errors(metric, runs, reference_run, ref_value):
    """Per scheme, list of (N, errTol, rung, err) tuples vs. (reference_run,
    ref_value), excluding the reference run itself and non-positive/None
    errors (zero error breaks the log scale; a real duplicate is not
    informative on a log-log plot either)."""
    out = {}
    for scheme in SCHEMES:
        points = []
        for run in runs.get(scheme, []):
            value = run['scalars'].get(metric)
            if value is None or run is reference_run:
                continue
            err = abs(value - ref_value) / abs(ref_value)
            if err <= 0:
                continue
            errtol = run.get('grid_tolerances', {}).get('errTol')
            points.append((run['N'], errtol, run.get('rung'), err))
        out[scheme] = points
    return out


def plot_case_metric(case, metric, runs, outdir):
    reference_run, ref_value = reference_run_and_value(case, metric, runs)
    if reference_run is None:
        return
    errors = scheme_errors(metric, runs, reference_run, ref_value)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    any_plotted = False
    for scheme in SCHEMES:
        points = [(N, err) for N, errtol, rung, err in errors[scheme]]
        if len(points) < 2:
            continue
        any_plotted = True
        Ns, errs = zip(*points)
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


def plot_case_metric_vs_errtol(case, metric, runs, outdir):
    """The parity view (Task G6, spec §4): error vs. errTol (log-log), one
    line per scheme. At matched errTol the two schemes' errors should track
    within a bounded band (the ~2x parity acceptance / errTol^(1/6) drift
    documented in spec §2)."""
    reference_run, ref_value = reference_run_and_value(case, metric, runs)
    if reference_run is None:
        return
    errors = scheme_errors(metric, runs, reference_run, ref_value)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    any_plotted = False
    for scheme in SCHEMES:
        points = [(errtol, err) for N, errtol, rung, err in errors[scheme]
                  if errtol is not None]
        if len(points) < 2:
            continue
        any_plotted = True
        points.sort()
        errtols, errs = zip(*points)
        ax.loglog(errtols, errs, marker=MARKERS[scheme], markersize=8, linewidth=2,
                   color=COLORS[scheme], label=scheme)

    if not any_plotted:
        print('  [skip] %s / %s (vs errTol): fewer than 2 comparable points for any scheme' %
              (case, metric))
        plt.close(fig)
        return

    ax.invert_xaxis()  # tighter tolerance (smaller errTol) to the right
    ax.set_xlabel('errTol')
    ax.set_ylabel('relative error vs. finest secondOrderLimited run')
    ax.set_title('%s: %s vs. errTol (parity view)' % (case, metric))
    ax.grid(True, which='both', linewidth=0.5, color='#c3c2b7', alpha=0.5)
    ax.legend(frameon=False)
    fig.tight_layout()

    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, '%s_%s_vs_errtol.png' % (case, metric))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print('  wrote %s' % out_path)


def n_ratio_table(case, runs):
    """sol/fou N-ratio at matched errTol (matched by rung index, since both
    schemes share the same nominal errTol per rung -- see RUNGS in
    run_convergence.py). Returns a list of
    (rung, errTol_fou, errTol_sol, N_fou, N_sol, ratio) sorted by rung."""
    by_rung = {'firstOrderUpwind': {}, 'secondOrderLimited': {}}
    for scheme in SCHEMES:
        for run in runs.get(scheme, []):
            rung = run.get('rung')
            if rung is not None:
                by_rung[scheme][rung] = run
    rows = []
    for rung in sorted(set(by_rung['firstOrderUpwind']) & set(by_rung['secondOrderLimited'])):
        fou_run = by_rung['firstOrderUpwind'][rung]
        sol_run = by_rung['secondOrderLimited'][rung]
        n_fou = fou_run['N']
        n_sol = sol_run['N']
        rows.append((rung,
                      fou_run.get('grid_tolerances', {}).get('errTol'),
                      sol_run.get('grid_tolerances', {}).get('errTol'),
                      n_fou, n_sol,
                      n_fou / n_sol if n_sol else float('nan')))
    return rows


def write_errtol_error_table(cases, resultsdir, outdir):
    """Combined errTol -> measured-QoI-error markdown table across all
    cases/schemes/metrics (deliverable: docs artifact per spec §4), plus the
    sol/fou N-ratio at matched errTol per case."""
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, 'errtol_error_table.md')
    lines = ['# errTol -> measured QoI error',
             '',
             'Reference (per case/metric): finest (largest-N) '
             '`secondOrderLimited` run loaded from `%s`. Generated by '
             '`plot_convergence.py`; do not hand-edit.' % resultsdir,
             '']
    for case in cases:
        runs = load_runs(resultsdir, case)
        n_found = sum(len(v) for v in runs.values())
        if n_found == 0:
            continue
        lines.append('## %s' % case)
        lines.append('')
        for metric in METRICS:
            reference_run, ref_value = reference_run_and_value(case, metric, runs)
            if reference_run is None:
                continue
            errors = scheme_errors(metric, runs, reference_run, ref_value)
            lines.append('### %s (ref N=%d, value=%.6g)' %
                          (metric, reference_run['N'], ref_value))
            lines.append('')
            lines.append('| scheme | rung | errTol | N | rel. error |')
            lines.append('|---|---|---|---|---|')
            for scheme in SCHEMES:
                for N, errtol, rung, err in sorted(
                        errors[scheme], key=lambda p: (p[2] if p[2] is not None else -1)):
                    lines.append('| %s | %s | %s | %d | %.3e |' %
                                  (scheme, rung, errtol, N, err))
            lines.append('')
        ratio_rows = n_ratio_table(case, runs)
        if ratio_rows:
            lines.append('### N-ratio (fou/sol) at matched errTol')
            lines.append('')
            lines.append('| rung | errTol (fou) | errTol (sol) | N (fou) | N (sol) | N ratio |')
            lines.append('|---|---|---|---|---|---|')
            for rung, et_fou, et_sol, n_fou, n_sol, ratio in ratio_rows:
                lines.append('| %d | %s | %s | %d | %d | %.2f |' %
                              (rung, et_fou, et_sol, n_fou, n_sol, ratio))
            lines.append('')
    with open(out_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
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
            plot_case_metric_vs_errtol(case, metric, runs, outdir)

    write_errtol_error_table(args.cases, args.resultsdir, outdir)


if __name__ == '__main__':
    main()
