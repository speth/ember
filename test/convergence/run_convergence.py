#!/usr/bin/env python
"""
Grid-convergence study harness for the convection-scheme comparison (spec
§6.4, Task 2.2). This script (Task 2.1) defines the resolution ladder and
the three case configurations and can run any single (case, scheme, rung)
combination, writing one JSON file per run into ``--outdir``.

This script only *defines* the study; it does not execute the full matrix
(that is Task 2.2). Task 2.1 only smoke-tests the coarsest rung of each
case.

Cases (fixed strain rate and composition per case; see spec §6.4 and the
convection-scheme design doc):

- ``strained``: an "unbounded" strained premixed H2/O2/Ar counterflow flame.
  Built with ``General(unburnedLeft=False, fixedBurnedVal=False)`` on the
  default planar geometry. Tracing ``FlameSolver::updateBC()``
  (src/flameSolver.cpp): with ``unburnedLeft=False`` the burned/product
  index is ``jb=0`` and the unburned/reactant index is ``ju=jj``; with
  ``fixedBurnedVal=False`` and this being neither a twin nor a cylindrical
  flame, the left boundary falls through to
  ``BoundaryCondition::ZeroGradient`` (an open, "unbounded" product-side
  boundary that can extend outward as the grid adapts) while the right
  boundary (the physical reactant inlet) is ``BoundaryCondition::FixedValue``.
  This is the literal ``FixedValue``/``ZeroGradient`` pairing named in spec
  §6.4 for Case A, and is distinct from Case B's ``ControlVolume`` path.
  Composition/strain follow ``TestPremixedStrained``
  (test/python/test_flame_configs.py) and smoke_continuity_bc.py.

- ``twin``: planar twin flame (``General(twinFlame=True,
  unburnedLeft=False)``), which resolves to the ``ControlVolume`` left BC
  (symmetry/stagnation plane at x=0). Config follows
  ``example_twin`` / ``case_twin`` in run_baselines.py.

- ``cylindrical``: outwardly-propagating cylindrical flame
  (``General(flameGeometry='cylindrical', unburnedLeft=False,
  fixedLeftLocation=True)``, i.e. curvature parameter alpha=1 in
  ``OneDimGrid``), which also resolves to the ``ControlVolume`` left BC but
  exercises the curved-geometry (``rphalf`` weighting) path instead of the
  planar-twin path. Config follows ``example_cylindrical_outward`` /
  ``case_cylindrical_outward`` in run_baselines.py.

Resolution ladder: a 6-rung sequence of ``(vtol, dvtol, gridMax)`` tuples
(see RUNGS below), index 0 = coarsest, index 5 = finest, geometrically
spaced by roughly a factor of 1.5x in vtol/dvtol per rung (~7x from coarsest
to finest). This is a *nominal* ladder sized to produce "at least 5 rungs
spanning roughly 4x in N" per spec §6.4; the actual N achieved by each rung
is configuration-dependent (thin, highly-strained flames need more points
per unit vtol than slow ones) and should be confirmed/retuned in Task 2.2
against the observed grid sizes, without changing the rung *count* or
*intent* documented here.

Usage:
    pixi run python test/convergence/run_convergence.py \\
        --case {strained,twin,cylindrical} \\
        --scheme {firstOrderUpwind,secondOrderLimited} \\
        [--rung N | --rungs N [N ...]] \\
        [--damp-const X] \\
        [--outdir test/convergence/results] \\
        [--workdir build/test/convergence-work] \\
        [--retries 3] \\
        [--list-rungs]

With neither ``--rung`` nor ``--rungs`` given, all 6 rungs are run in
sequence (the Task 2.2 usage). Pass ``--rung 0`` to run only the coarsest
rung (used for the Task 2.1 smoke test).

See test/convergence/README.md for the full protocol and output schema.
"""

import argparse
import datetime
import os
import sys
import time
import json

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python'))
# Reuse the git-commit tag and the CVODE convection step-count log parser
# from the Task 0.1/1.5 baseline harness instead of duplicating them.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ember import (Config, Paths, General, Chemistry, Grid, InitialCondition,
                    StrainParameters, Times, TerminationCondition)

from run_baselines import git_commit, total_convection_steps

# ---------------------------------------------------------------------------
# Resolution ladder: grid-tolerance rungs, coarsest (index 0) to finest.
# See module docstring for the rationale/caveats.
# ---------------------------------------------------------------------------
RUNGS = [
    dict(vtol=0.24,  dvtol=0.40,  gridMax=4.0e-4),   # 0: coarsest
    dict(vtol=0.16,  dvtol=0.27,  gridMax=2.5e-4),   # 1
    dict(vtol=0.11,  dvtol=0.18,  gridMax=1.6e-4),   # 2
    dict(vtol=0.075, dvtol=0.12,  gridMax=1.0e-4),   # 3
    dict(vtol=0.050, dvtol=0.080, gridMax=6.3e-5),   # 4
    dict(vtol=0.033, dvtol=0.055, gridMax=4.0e-5),   # 5: finest
]

CASES = ('strained', 'twin', 'cylindrical')
SCHEMES = ('firstOrderUpwind', 'secondOrderLimited')


def _grid_kwargs(rung_idx, damp_const):
    r = RUNGS[rung_idx]
    kwargs = dict(vtol=r['vtol'], dvtol=r['dvtol'], gridMax=r['gridMax'])
    if damp_const is not None:
        kwargs['dampConst'] = damp_const
    return kwargs


def build_strained(rung_idx, scheme, damp_const, work_dir):
    tag = 'strained_%s_rung%d' % (scheme, rung_idx)
    output = os.path.join(work_dir, tag)
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        General(nThreads=1,
                unburnedLeft=False,
                fixedBurnedVal=False,
                convectionScheme=scheme),
        Chemistry(mechanismFile='h2o2.yaml'),
        InitialCondition(fuel='H2:1.0',
                          oxidizer='O2:1.0, AR:4.0',
                          equivalenceRatio=0.3),
        StrainParameters(initial=800, final=800),
        Grid(**_grid_kwargs(rung_idx, damp_const)),
        Times(regridStepInterval=10),
        TerminationCondition(tEnd=2.0, measurement='dTdt'))
    return conf


def build_twin(rung_idx, scheme, damp_const, work_dir):
    tag = 'twin_%s_rung%d' % (scheme, rung_idx)
    output = os.path.join(work_dir, tag)
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        General(nThreads=1,
                twinFlame=True,
                unburnedLeft=False,
                convectionScheme=scheme),
        InitialCondition(fuel='CH4:1.0',
                          equivalenceRatio=0.70,
                          xLeft=0.0,
                          xRight=0.01),
        StrainParameters(initial=100, final=100),
        Grid(**_grid_kwargs(rung_idx, damp_const)),
        Times(regridStepInterval=10),
        TerminationCondition(tEnd=10))
    return conf


def build_cylindrical(rung_idx, scheme, damp_const, work_dir):
    tag = 'cylindrical_%s_rung%d' % (scheme, rung_idx)
    output = os.path.join(work_dir, tag)
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        Chemistry(mechanismFile='gri30.yaml'),
        General(nThreads=1,
                flameGeometry='cylindrical',
                unburnedLeft=False,
                fixedLeftLocation=True,
                convectionScheme=scheme),
        InitialCondition(fuel='CH4:0.5, H2:0.5',
                          equivalenceRatio=0.60,
                          xLeft=0.0,
                          xRight=0.005),
        StrainParameters(initial=500, final=500),
        Grid(**_grid_kwargs(rung_idx, damp_const)),
        TerminationCondition(tEnd=10, measurement='dTdt'),
        Times(profileStepInterval=10, regridStepInterval=10))
    return conf


BUILDERS = {
    'strained': build_strained,
    'twin': build_twin,
    'cylindrical': build_cylindrical,
}


def safe_consumption_speed(solver):
    """Return (value_or_None, note_or_None). Mirrors the ill-conditioning
    guard in run_baselines.py's run_case(): the consumption-speed formula
    has a near-zero denominator whenever the two domain-boundary
    temperatures are nearly equal."""
    raw = float(solver.consumptionSpeed)
    T = np.asarray(solver.T)
    boundary_dT = abs(float(T[0]) - float(T[-1]))
    if not np.isfinite(raw) or boundary_dT < 5.0:
        return None, ('consumption_speed not physically meaningful: '
                       'boundary temperatures nearly equal '
                       '(T[0]=%.2f K, T[-1]=%.2f K); raw value %r discarded'
                       % (T[0], T[-1], raw))
    return raw, None


def run_once(case, scheme, rung_idx, damp_const, work_dir):
    conf = BUILDERS[case](rung_idx, scheme, damp_const, work_dir)
    log_path = conf.paths.logFile.value

    conf.validate()
    concrete = conf.evaluate()

    t0 = time.time()
    solver = concrete.run()
    runtime = time.time() - t0

    T = np.asarray(solver.T)
    Y = np.asarray(solver.Y)
    x = np.asarray(solver.x)
    if not np.all(np.isfinite(T)) or not np.all(np.isfinite(Y)):
        raise RuntimeError('Non-finite values encountered in T or Y '
                            '(case=%r scheme=%r rung=%d)' % (case, scheme, rung_idx))

    consumption_speed, consumption_speed_note = safe_consumption_speed(solver)

    result = {
        'case': case,
        'scheme': scheme,
        'rung': rung_idx,
        'damp_const': damp_const,
        'commit': git_commit(),
        'generated_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'config_summary': conf.stringify(),
        'grid_tolerances': {
            'vtol': concrete.grid.vtol,
            'dvtol': concrete.grid.dvtol,
            'gridMax': concrete.grid.gridMax,
            'gridMin': concrete.grid.gridMin,
            'dampConst': concrete.grid.dampConst,
        },
        'N': int(len(x)),
        'scalars': {
            'consumption_speed': consumption_speed,
            'peak_T': float(np.max(T)),
        },
        'scalar_notes': ({'consumption_speed': consumption_speed_note}
                          if consumption_speed_note else {}),
        'total_convection_steps': total_convection_steps(log_path),
        'runtime_seconds': runtime,
        'final_time': float(solver.tNow),
    }
    return result


def run_with_retries(case, scheme, rung_idx, damp_const, work_dir, retries):
    attempt = 0
    while True:
        attempt += 1
        try:
            result = run_once(case, scheme, rung_idx, damp_const, work_dir)
            break
        except Exception as exc:
            print('    attempt %d/%d failed: %r' % (attempt, retries, exc))
            if attempt >= retries:
                raise
    result['attempts'] = attempt
    return result


def output_filename(case, scheme, rung_idx, damp_const):
    name = '%s_%s_rung%d' % (case, scheme, rung_idx)
    if damp_const is not None:
        name += '_damp%g' % damp_const
    return name + '.json'


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', choices=CASES,
                         help='Case to run (required unless --list-rungs)')
    parser.add_argument('--scheme', choices=SCHEMES,
                         help='Convection scheme (required unless --list-rungs)')
    parser.add_argument('--damp-const', type=float, default=None, dest='damp_const',
                         help='Override Grid.dampConst for this run (spec §6.4\'s '
                              '"relaxed dampConst" trial). Default: use each rung\'s '
                              'own default (Grid.dampConst=7).')
    parser.add_argument('--rung', type=int, default=None, choices=range(len(RUNGS)),
                         help='Run a single rung index (0=coarsest, %d=finest).'
                              % (len(RUNGS) - 1))
    parser.add_argument('--rungs', type=int, nargs='+', default=None,
                         choices=range(len(RUNGS)),
                         help='Run an explicit list of rung indices. Overrides --rung. '
                              'Default (neither given): run all rungs.')
    parser.add_argument('--outdir', default='test/convergence/results',
                         help='Directory to write per-run JSON files into')
    parser.add_argument('--workdir', default='build/test/convergence-work',
                         help='Scratch directory for Ember run outputs (HDF5/log files)')
    parser.add_argument('--retries', type=int, default=3,
                         help='Max attempts per run before giving up (known flake: '
                              'occasional CVODE integrator error under multi-threaded '
                              'or stiff runs; retrying is usually sufficient). Default: 3')
    parser.add_argument('--list-rungs', action='store_true',
                         help='Print the rung ladder (vtol/dvtol/gridMax per index) and exit')
    args = parser.parse_args()

    if args.list_rungs:
        print('%-5s %-8s %-8s %-10s' % ('rung', 'vtol', 'dvtol', 'gridMax'))
        for i, r in enumerate(RUNGS):
            print('%-5d %-8g %-8g %-10g' % (i, r['vtol'], r['dvtol'], r['gridMax']))
        return

    if args.case is None or args.scheme is None:
        parser.error('--case and --scheme are required (unless --list-rungs)')

    if args.rungs is not None:
        rungs_to_run = list(args.rungs)
    elif args.rung is not None:
        rungs_to_run = [args.rung]
    else:
        rungs_to_run = list(range(len(RUNGS)))

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    total_t0 = time.time()
    for rung_idx in rungs_to_run:
        print('=== case=%s scheme=%s rung=%d damp_const=%s ===' %
              (args.case, args.scheme, rung_idx, args.damp_const))
        result = run_with_retries(args.case, args.scheme, rung_idx,
                                   args.damp_const, args.workdir, args.retries)
        print('    finished in %.1f s (attempts=%d, N=%d, final_time=%.5g, '
              'peak_T=%.1f, consumption_speed=%s, convection_steps=%s)'
              % (result['runtime_seconds'], result['attempts'], result['N'],
                 result['final_time'], result['scalars']['peak_T'],
                 result['scalars']['consumption_speed'],
                 result['total_convection_steps']))
        out_path = os.path.join(
            args.outdir, output_filename(args.case, args.scheme, rung_idx, args.damp_const))
        with open(out_path, 'w') as f:
            json.dump(result, f, indent=2)
        print('    wrote %s' % out_path)

    print('=== Total wall time: %.1f s ===' % (time.time() - total_t0))


if __name__ == '__main__':
    main()
