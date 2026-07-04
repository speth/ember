#!/usr/bin/env python
"""
Compare two baseline JSON files produced by run_baselines.py, reporting
relative differences in the captured scalars and L2/L-infinity norms of the
profile differences (T and major species, interpolated to a common grid).

Usage:
    pixi run python test/convergence/compare_baselines.py \\
        test/convergence/baselines/example_single.json \\
        test/convergence/baselines/example_single_run2.json \\
        [--threshold 0.02]

Exits with status 0 if every relative difference is within --threshold,
and status 1 otherwise (so this can be used as a pass/fail gate in later
tasks, e.g. Task 1.5 comparing the new convection scheme against these
Phase-0 baselines).
"""

import argparse
import json
import sys

import numpy as np


def load(path):
    with open(path) as f:
        return json.load(f)


def relative_diff(a, b):
    """Relative difference of two scalars, robust to values near zero."""
    if a is None or b is None:
        return None
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def compare_scalars(a, b):
    rows = []
    keys = sorted(set(a['scalars']) | set(b['scalars']))
    for key in keys:
        va = a['scalars'].get(key)
        vb = b['scalars'].get(key)
        rows.append((key, va, vb, relative_diff(va, vb)))
    return rows


def interpolate_to_common_grid(xa, ya, xb, yb, n=1000):
    """Interpolate two profiles (xa,ya) and (xb,yb) onto a shared uniform
    grid spanning the overlap of their domains."""
    xa = np.asarray(xa)
    ya = np.asarray(ya)
    xb = np.asarray(xb)
    yb = np.asarray(yb)

    lo = max(xa.min(), xb.min())
    hi = min(xa.max(), xb.max())
    if hi <= lo:
        raise ValueError('Profiles do not overlap in x')

    grid = np.linspace(lo, hi, n)
    ia = np.interp(grid, xa, ya)
    ib = np.interp(grid, xb, yb)
    return ia, ib


def profile_norms(xa, ya, xb, yb):
    ia, ib = interpolate_to_common_grid(xa, ya, xb, yb)
    diff = ia - ib
    scale = max(np.max(np.abs(ia)), np.max(np.abs(ib)), 1e-300)
    l2 = np.sqrt(np.mean(diff**2)) / scale
    linf = np.max(np.abs(diff)) / scale
    return l2, linf


def compare_profiles(a, b):
    rows = []
    xa, xb = a['profiles']['x'], b['profiles']['x']

    l2, linf = profile_norms(xa, a['profiles']['T'], xb, b['profiles']['T'])
    rows.append(('T', l2, linf))

    species = sorted(set(a['profiles']['species']) & set(b['profiles']['species']))
    for sp in species:
        l2, linf = profile_norms(xa, a['profiles']['species'][sp],
                                  xb, b['profiles']['species'][sp])
        rows.append(('Y_' + sp, l2, linf))

    only_a = sorted(set(a['profiles']['species']) - set(b['profiles']['species']))
    only_b = sorted(set(b['profiles']['species']) - set(a['profiles']['species']))
    return rows, only_a, only_b


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('baselineA')
    parser.add_argument('baselineB')
    parser.add_argument('--threshold', type=float, default=0.02,
                         help='Relative-difference pass/fail threshold, applied to '
                              'both scalars and normalized profile norms (default: 0.02)')
    args = parser.parse_args()

    a = load(args.baselineA)
    b = load(args.baselineB)

    print('Comparing:')
    print('  A: %s (case=%s, commit=%s)' % (args.baselineA, a.get('case'), a.get('commit')))
    print('  B: %s (case=%s, commit=%s)' % (args.baselineB, b.get('case'), b.get('commit')))
    if a.get('case') != b.get('case'):
        print('WARNING: comparing different cases (%r vs %r)' % (a.get('case'), b.get('case')))
    print()

    worst = 0.0

    print('Scalars (name: A, B, relative diff):')
    for key, va, vb, rel in compare_scalars(a, b):
        if rel is None:
            print('  %-28s %r vs %r (skipped: null value)' % (key, va, vb))
            continue
        worst = max(worst, rel)
        flag = '' if rel <= args.threshold else '  <-- EXCEEDS THRESHOLD'
        print('  %-28s %.6g vs %.6g   rel_diff=%.4g%s' % (key, va, vb, rel, flag))
    print()

    print('Profiles (normalized L2 / Linf norms of the difference, common grid):')
    rows, only_a, only_b = compare_profiles(a, b)
    for name, l2, linf in rows:
        worst = max(worst, l2, linf)
        flag = '' if max(l2, linf) <= args.threshold else '  <-- EXCEEDS THRESHOLD'
        print('  %-28s L2=%.4g  Linf=%.4g%s' % (name, l2, linf, flag))
    if only_a:
        print('  Species only in A (not compared): %s' % ', '.join(only_a))
    if only_b:
        print('  Species only in B (not compared): %s' % ', '.join(only_b))
    print()

    passed = worst <= args.threshold
    print('Result: %s (worst relative difference/norm = %.4g, threshold = %.4g)'
          % ('PASS' if passed else 'FAIL', worst, args.threshold))
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
