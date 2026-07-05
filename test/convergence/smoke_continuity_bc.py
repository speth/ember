#!/usr/bin/env python
"""
Smoke test for the never-exercised continuity-BC branches of the
trapezoidal (SecondOrderLimited) continuity march in ConvectionSystemUTW::f().

Context: every case in run_baselines.py / test_flame_configs.py uses the
default General.continuityBC='fixedLeft'. The `stagnationPoint` and
`fixedTemperature` branches (ContinuityBoundaryCondition::Zero / ::Temp in
convectionSystem.cpp) run a structurally different rV march (forward from
jContBC AND backward from jContBC, vs. a single forward march from rV[0] for
fixedLeft) and had never actually executed under the new `secondOrderLimited`
trapezoidal scheme before this smoke test. This script is a does-it-execute-
sanely check (short runs), not a physics/convergence study -- see
test/convergence/README.md "Phase 1 review fixes" section for the results.

Matrix: continuityBC in {stagnationPoint, fixedTemperature} x
        convectionScheme in {secondOrderLimited, firstOrderUpwind}
      = 4 short runs.

Base config: a lightly modified version of TestPremixedStrained from
test/python/test_flame_configs.py (H2/O2/Ar premixed counterflow flame,
h2o2.yaml mechanism -- small and fast), with tEnd shortened only enough to
keep each run in the "few hundred timesteps" range while still letting the
solver reach a well-developed profile with an interior stagnation point and
temperature gradient (needed for the BCs under test to have something
meaningful to lock onto). fixedTemperature requires
General.splittingMethod='balanced', which is already Ember's default, so no
extra options are needed beyond continuityBC itself.

Usage:
    pixi run python test/convergence/smoke_continuity_bc.py
"""

import os
import sys
import traceback

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python'))

from ember import (Config, Paths, General, Chemistry, InitialCondition,
                    StrainParameters, Grid, Times, TerminationCondition)

WORK_DIR = os.path.join(REPO_ROOT, 'build', 'test', 'smoke-continuity-bc')

MATRIX = [
    ('stagnationPoint', 'secondOrderLimited'),
    ('stagnationPoint', 'firstOrderUpwind'),
    ('fixedTemperature', 'secondOrderLimited'),
    ('fixedTemperature', 'firstOrderUpwind'),
]


def build_config(continuity_bc, scheme):
    tag = '%s_%s' % (continuity_bc, scheme)
    output = os.path.join(WORK_DIR, tag)
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        General(nThreads=1,
                continuityBC=continuity_bc,
                convectionScheme=scheme),
        Chemistry(mechanismFile='h2o2.yaml'),
        InitialCondition(fuel='H2:1.0',
                          oxidizer='O2:1.0, AR:4.0',
                          equivalenceRatio=0.3),
        StrainParameters(initial=800, final=800),
        Grid(vtol=0.2, dvtol=0.3),
        Times(regridStepInterval=10),
        TerminationCondition(tEnd=0.006, measurement=None))
    return conf


def run_one(continuity_bc, scheme):
    conf = build_config(continuity_bc, scheme)
    result = {
        'continuityBC': continuity_bc,
        'scheme': scheme,
        'completed': False,
        'error': None,
        'final_time': None,
        'grid_size': None,
        'peak_T': None,
        'min_T': None,
        'finite': None,
    }
    try:
        ok = conf.validate()
        if not ok:
            result['error'] = 'conf.validate() returned False'
            return result
        concrete = conf.evaluate()
        solver = concrete.run()
        T = np.asarray(solver.T)
        Y = np.asarray(solver.Y)
        finite = bool(np.all(np.isfinite(T)) and np.all(np.isfinite(Y)))
        result['completed'] = True
        result['final_time'] = float(solver.tNow)
        result['grid_size'] = int(len(np.asarray(solver.x)))
        result['peak_T'] = float(np.max(T)) if finite else float('nan')
        result['min_T'] = float(np.min(T)) if finite else float('nan')
        result['finite'] = finite
    except Exception as exc:
        result['error'] = '%r\n%s' % (exc, traceback.format_exc())
    return result


def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    results = []
    for continuity_bc, scheme in MATRIX:
        print('=== continuityBC=%s scheme=%s ===' % (continuity_bc, scheme))
        result = run_one(continuity_bc, scheme)
        results.append(result)
        if result['completed']:
            print('    completed: final_time=%.5g grid_size=%d peak_T=%.1f '
                  'min_T=%.1f finite=%s'
                  % (result['final_time'], result['grid_size'],
                     result['peak_T'], result['min_T'], result['finite']))
        else:
            print('    FAILED: %s' % result['error'])

    print()
    print('%-20s %-18s %-10s %-10s %-8s %-8s' %
          ('continuityBC', 'scheme', 'completed', 'final_t', 'peak_T', 'finite'))
    any_fail = False
    for r in results:
        ok = r['completed'] and r['finite'] and 250 <= (r['peak_T'] or -1) <= 3000
        if not ok:
            any_fail = True
        print('%-20s %-18s %-10s %-10s %-8s %-8s' %
              (r['continuityBC'], r['scheme'], r['completed'],
               '%.4g' % r['final_time'] if r['final_time'] is not None else 'n/a',
               '%.1f' % r['peak_T'] if r['peak_T'] is not None else 'n/a',
               r['finite']))
    sys.exit(1 if any_fail else 0)


if __name__ == '__main__':
    main()
