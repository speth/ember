#!/usr/bin/env python
"""
Capture baseline outputs from a curated subset of Ember's example
configurations, for later regression comparison against the modified
convection solver (see the convection-scheme design doc / implementation
plan, Phase 0-1).

Each case below is a faithful transcription of the corresponding
``python/ember/examples/example_*.py`` script, with only the output
paths redirected into a scratch directory (``--outdir``). Any other
deviation from the stock example configuration is recorded per-case in
the ``deviations_from_stock`` list and written into the output JSON.

Usage:
    pixi run python test/convergence/run_baselines.py \\
        --outdir test/convergence/baselines \\
        [--cases example_single example_diffusion ...]

See test/convergence/README.md for details on regenerating baselines
and on the run-to-run reproducibility tolerance floor.
"""

import argparse
import datetime
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'python'))

from ember import (Config, Paths, General, Chemistry, Grid, InitialCondition,
                    StrainParameters, PositionControl, Times,
                    TerminationCondition)

# Minimum peak mass fraction for a species to be considered "major" and
# included in the captured profiles.
MAJOR_SPECIES_MIN_Y = 1e-3
# Upper bound on the number of major species profiles captured, so JSON
# files stay a reasonable size regardless of mechanism size.
MAJOR_SPECIES_MAX_COUNT = 8


def git_commit():
    """Return 'sha' or 'sha-dirty' for the current checkout, or None.

    "Dirty" only reflects changes to tracked files (--untracked-files=no):
    untracked scratch/lock files in the working tree don't affect what code
    actually ran, so they shouldn't taint the recorded commit identity.
    """
    try:
        sha = subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=REPO_ROOT).decode().strip()
        status = subprocess.check_output(
            ['git', 'status', '--porcelain', '--untracked-files=no'],
            cwd=REPO_ROOT).decode().strip()
        return sha + ('-dirty' if status else '')
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Case definitions. Each function builds a Config object equivalent to the
# stock example (see the docstring/comment for the source file) and returns
# (conf, deviations) where deviations is a list of human-readable strings
# describing any difference from the stock example configuration (other than
# the output path redirection, which is common to all cases and noted once).
# ---------------------------------------------------------------------------

def case_single(work_dir):
    # Source: python/ember/examples/example_single.py
    output = os.path.join(work_dir, 'ex_single')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        General(twinFlame=False, flameGeometry='disc', nThreads=4),
        InitialCondition(fuel='CH4:1.0',
                          oxidizer='N2:3.76, O2:1.0',
                          equivalenceRatio=1.0,
                          counterflow='N2:1.0',
                          Tcounterflow=300.0,
                          xLeft=-0.01,
                          xRight=0.01,
                          centerWidth=0.005,
                          slopeWidth=0.001),
        StrainParameters(initial=300.0, final=300.0),
    )
    return conf, []


def case_diffusion(work_dir):
    # Source: python/ember/examples/example_diffusion.py
    output = os.path.join(work_dir, 'ex_diffusion')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        InitialCondition(flameType='diffusion',
                          fuel='CH4:1.0, N2:2.0',
                          oxidizer='N2:3.76, O2:1.0',
                          Tfuel=600,
                          Toxidizer=600,
                          xLeft=-0.004,
                          xRight=0.004,
                          centerWidth=0.002,
                          slopeWidth=0.001),
        StrainParameters(initial=100, final=100),
        General(nThreads=2),
        Times(globalTimestep=1e-5, profileStepInterval=20),
        TerminationCondition(tEnd=0.010))
    return conf, []


def case_twin(work_dir):
    # Source: python/ember/examples/example_twin.py
    output = os.path.join(work_dir, 'ex_twin')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        General(twinFlame=True, unburnedLeft=False, nThreads=4),
        InitialCondition(fuel='CH4:1.0',
                          equivalenceRatio=0.70,
                          xLeft=0.0,
                          xRight=0.01),
        StrainParameters(initial=100, final=100),
        TerminationCondition(tEnd=10))
    return conf, []


def case_cylindrical_outward(work_dir):
    # Source: python/ember/examples/example_cylindrical_outward.py
    output = os.path.join(work_dir, 'ex_cylindrical_outward')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        Chemistry(mechanismFile='gri30.yaml'),
        General(flameGeometry='cylindrical',
                unburnedLeft=False,
                fixedLeftLocation=True,
                nThreads=4),
        InitialCondition(fuel='CH4:0.5, H2:0.5',
                          equivalenceRatio=0.60,
                          xLeft=0.0,
                          xRight=0.005),
        StrainParameters(initial=500, final=500),
        TerminationCondition(tEnd=10, measurement='dTdt'),
        Times(profileStepInterval=10, regridStepInterval=10),
    )
    return conf, []


def case_cylindrical_inward(work_dir):
    # Source: python/ember/examples/example_cylindrical_inward.py
    output = os.path.join(work_dir, 'ex_cylindrical_inward')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        Chemistry(mechanismFile='gri30.yaml'),
        General(flameGeometry='cylindrical',
                unburnedLeft=True,
                fixedLeftLocation=True,
                nThreads=4),
        InitialCondition(fuel='CH4:0.5, H2:0.5',
                          equivalenceRatio=0.60,
                          xLeft=0.0,
                          xRight=0.006),
        StrainParameters(initial=200, final=200),
        PositionControl(xInitial=0.002, xFinal=0.002),
        TerminationCondition(tEnd=10, measurement='dTdt'),
        Times(profileStepInterval=10, regridStepInterval=10),
    )
    return conf, []


def case_laminar_flame_speed(work_dir):
    # Source: python/ember/examples/example_laminarFlameSpeed.py
    output = os.path.join(work_dir, 'ex_lfs')
    conf = Config(
        Paths(outputDir=output, logFile=output + '.log'),
        InitialCondition(fuel='CH4:1.0',
                          oxidizer='O2:1, N2:3.76',
                          equivalenceRatio=0.9,
                          xLeft=0.0,
                          xRight=0.01),
        StrainParameters(initial=0, final=0),
        General(fixedLeftLocation=True, fixedBurnedVal=False, nThreads=4),
        Grid(vtol=0.1, dvtol=0.15, gridMin=5e-6, gridMax=0.001),
        PositionControl(proportionalGain=2000, xInitial=0.005, xFinal=0.005),
        TerminationCondition(tolerance=1e-5),
        Times(profileStepInterval=50))
    deviations = [
        'Note (inherited from stock example, not introduced here): '
        'conf.validate() prints "PositionControl can only be used when '
        'either \'twinFlame\' or \'cylindricalFlame\' is set to True. '
        'Validation failed." This example configures PositionControl on a '
        'planar, non-twin flame; the stock example script has this same '
        'validation warning and still runs it. Non-fatal; run proceeds.',
    ]
    return conf, deviations


CASES = {
    'example_single': case_single,
    'example_diffusion': case_diffusion,
    'example_twin': case_twin,
    'example_cylindrical_outward': case_cylindrical_outward,
    'example_cylindrical_inward': case_cylindrical_inward,
    'example_laminarFlameSpeed': case_laminar_flame_speed,
}


def select_major_species(species_names, Y):
    """Return the names of up to MAJOR_SPECIES_MAX_COUNT species whose peak
    mass fraction (over the final profile) is at least MAJOR_SPECIES_MIN_Y,
    ordered from most to least abundant."""
    peak = Y.max(axis=1)
    order = np.argsort(peak)[::-1]
    selected = [i for i in order if peak[i] >= MAJOR_SPECIES_MIN_Y]
    selected = selected[:MAJOR_SPECIES_MAX_COUNT]
    return [species_names[i] for i in selected]



# Matches the per-global-timestep debug log line written by FlameSolver when
# Debug.timesteps is enabled (the default), e.g.:
#   t = 0.000100 (dt = 1.000e-04) [C: 12]
# The trailing integer is ConvectionSystemSplit::getNumSteps(), i.e. the
# number of CVODE steps taken by the split convection sub-integrators
# *since the last reinitialization* (which happens once per global
# timestep in FlameSolver::setState). It is therefore a per-step count, not
# a running total; to get a whole-run total we sum every occurrence in the
# log file.
CONVECTION_STEPS_RE = re.compile(r'\[C:\s*(\d+)\]')


def total_convection_steps(log_path):
    """Sum the per-timestep convection CVODE step counts logged when
    Debug.timesteps is enabled, across an entire run's log file. Returns
    None if the log file is missing or contains no matching lines (e.g.
    Debug.timesteps was disabled)."""
    if not log_path or not os.path.isfile(log_path):
        return None
    total = 0
    found = False
    with open(log_path) as f:
        for line in f:
            m = CONVECTION_STEPS_RE.search(line)
            if m:
                found = True
                total += int(m.group(1))
    return total if found else None


def run_case(name, work_dir, scheme=None):
    build_fn = CASES[name]
    conf, deviations = build_fn(work_dir)
    deviations = (['Paths.outputDir/logFile redirected into scratch --outdir '
                    '(no effect on physics)'] + deviations)

    if scheme is not None:
        conf.general.convectionScheme.value = scheme
        conf.general.convectionScheme.isSet = True
        deviations.append(
            'General.convectionScheme overridden to %r via --scheme '
            '(harness extension for Task 1.5; no other physics changes)'
            % scheme)

    log_path = conf.paths.logFile.value

    conf.validate()
    concrete = conf.evaluate()

    t0 = time.time()
    solver = concrete.run()
    runtime = time.time() - t0

    species_names = list(concrete.gas.species_names)
    T = np.asarray(solver.T)
    x = np.asarray(solver.x)
    Y = np.asarray(solver.Y)

    major = select_major_species(species_names, Y)
    species_profiles = {sp: Y[species_names.index(sp)].tolist() for sp in major}

    peak_T = float(np.max(T))
    if not np.all(np.isfinite(T)) or not np.all(np.isfinite(Y)):
        raise RuntimeError('Non-finite values encountered in T or Y for case %r' % name)

    def scalar_or_none(value):
        value = float(value)
        return value if np.isfinite(value) else None

    # The consumption-speed formula (Q/cp integral, normalized by
    # rhou*(Tb - Tu)) is ill-conditioned when the two domain-boundary
    # temperatures are nearly equal (e.g. a flame sandwiched between two
    # comparably-cold/comparably-warm streams, as in a counterflow
    # diffusion flame or a premixed flame opposed by a cold inert of the
    # same temperature as the reactants). That near-zero denominator can
    # produce NaN/inf (caught by scalar_or_none) or, due to floating-point
    # error, a huge-but-finite garbage value (not caught by isfinite).
    scalar_notes = {}
    raw_consumption_speed = float(solver.consumptionSpeed)
    boundary_dT = abs(float(T[0]) - float(T[-1]))
    if not np.isfinite(raw_consumption_speed) or boundary_dT < 5.0:
        consumption_speed = None
        scalar_notes['consumption_speed'] = (
            'Not physically meaningful for this configuration: the domain '
            'boundary temperatures are nearly equal (T[0]=%.2f K, '
            'T[-1]=%.2f K), so the premixed consumption-speed formula '
            'Q/cp integral / (rho_u*(Tb-Tu)) has a near-zero denominator. '
            'Raw solver value (%r) discarded and replaced with null.'
            % (T[0], T[-1], raw_consumption_speed))
    else:
        consumption_speed = raw_consumption_speed

    result = {
        'case': name,
        'commit': git_commit(),
        'generated_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'config_summary': conf.stringify(),
        'deviations_from_stock': deviations,
        # Convection discretization scheme actually used for this run
        # (General.convectionScheme, post any --scheme override).
        'scheme': scheme if scheme is not None else concrete.general.convectionScheme,
        'runtime_seconds': runtime,
        'final_time': float(solver.tNow),
        'grid_size': int(len(x)),
        # Sum of the per-global-timestep convection CVODE step counts (see
        # total_convection_steps() above); None if Debug.timesteps was off
        # or the log file wasn't found. Not included under 'scalars' so it
        # is not swept into compare_baselines.py's scalar pass/fail check.
        'total_convection_steps': total_convection_steps(log_path),
        'scalars': {
            'peak_T': peak_T,
            'consumption_speed': consumption_speed,
            'heat_release_rate_integral': scalar_or_none(solver.heatReleaseRate),
            'flame_position': scalar_or_none(solver.flamePosition),
        },
        'scalar_notes': scalar_notes,
        'profiles': {
            'x': x.tolist(),
            'T': T.tolist(),
            'species': species_profiles,
        },
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--outdir', default='test/convergence/baselines',
                         help='Directory to write baseline JSON files into')
    parser.add_argument('--workdir', default='build/test/baselines-work',
                         help='Scratch directory for Ember run outputs (HDF5 profiles, logs)')
    parser.add_argument('--cases', nargs='+', choices=sorted(CASES), default=sorted(CASES),
                         help='Subset of cases to run (default: all curated cases)')
    parser.add_argument('--suffix', default='',
                         help='Optional suffix appended to output JSON filenames, '
                              'e.g. "_run2" for a repeatability check')
    parser.add_argument('--retries', type=int, default=3,
                         help='Max attempts per case. Some cases occasionally fail '
                              'with a CVODE/integrator error under multi-threaded '
                              'execution (observed thread-scheduling nondeterminism, '
                              'see README.md); retrying is usually sufficient. '
                              '(default: 3)')
    parser.add_argument('--scheme', choices=['firstOrderUpwind', 'secondOrderLimited'],
                         default=None,
                         help='Override General.convectionScheme for every case in this '
                              'run (harness extension added for Task 1.5). Default: no '
                              'override, i.e. use each case\'s configured default '
                              '(currently secondOrderLimited).')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    total_t0 = time.time()
    for name in args.cases:
        print('=== Running %s ===' % name)
        attempt = 0
        while True:
            attempt += 1
            try:
                result = run_case(name, args.workdir, scheme=args.scheme)
                break
            except Exception as exc:
                print('    attempt %d/%d failed: %r' % (attempt, args.retries, exc))
                if attempt >= args.retries:
                    raise
        result['attempts'] = attempt
        print('    finished in %.1f s (attempts=%d, final_time=%.5g, grid_size=%d, '
              'peak_T=%.1f, scheme=%s, convection_steps=%s)'
              % (result['runtime_seconds'], attempt, result['final_time'],
                 result['grid_size'], result['scalars']['peak_T'],
                 result['scheme'], result['total_convection_steps']))
        out_path = os.path.join(args.outdir, name + args.suffix + '.json')
        with open(out_path, 'w') as f:
            json.dump(result, f, indent=2)
        print('    wrote %s' % out_path)

    print('=== Total wall time: %.1f s ===' % (time.time() - total_t0))


if __name__ == '__main__':
    main()
