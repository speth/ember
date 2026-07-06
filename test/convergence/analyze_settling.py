"""
Grid-settling analysis for errTol-ladder runs: extracts the grid-size
trajectory from the numbered profNNNNNN.h5 outputs of a case work directory
(profNow.h5, a continuously-rewritten current-state file, is excluded) and
tests whether the grid reached a plateau (the P2.4 failure signature is a
monotonically growing, never-settling grid).

Usage: python analyze_settling.py <workdir> [<workdir> ...]
"""
import sys
import os
import glob
import math
from statistics import median
import h5py


def n_trajectory(case_dir):
    """Grid point count for each numbered profile output, in time order."""
    files = sorted(glob.glob(os.path.join(case_dir, 'prof[0-9]*.h5')))
    traj = []
    for f in files:
        with h5py.File(f, 'r') as h:
            traj.append(int(h['x'].shape[0]))
    return traj


def grid_settled(traj, tail_frac=0.25, tol_pts=None):
    """True if the trajectory ends in a plateau: N varies by <= tol_pts over
    the tail window (the last tail_frac of outputs, at least 3 samples).
    Fewer than 4 samples is not judgeable and returns False; callers should
    treat that case as "inspect the run log", not as a settling failure
    (main() reports it as TOO FEW OUTPUTS).

    ``tol_pts=None`` (the default) uses a tolerance relative to the tail's
    median grid size, ``max(3, ceil(REL_TOL * median(tail)))``, rather than a
    fixed absolute point count: an absolute tol_pts=3 falsely flags runs
    with N in the hundreds as "not settled" even when they've plateaued
    (observed on 14/36 runs of the errTol ladder at N > ~300; see spec
    addendum §A.5). Pass an explicit numeric tol_pts to restore the old
    fixed-tolerance behavior.

    REL_TOL=0.05 (5%), not the 1% first floated in §A.5's "e.g." aside:
    verifying against the actual errTol-ladder data
    (build/test/convergence-work-errtol-final/strained_secondOrderLimited_rung5,
    tail spread 7 pts at tail-median N=218, i.e. ~3.2%) showed 1% still
    under-tolerates (tol_pts would round to 3, same as the old absolute
    floor). 5% clears that case with margin while staying far below the
    ~32% relative spread of the genuine P2.4 monotonic-growth failure
    (build/test/convergence-work/strained_secondOrderLimited_rung5, spread
    94 pts at tail-median N=291), so it still correctly reports that one as
    NOT SETTLED.
    """
    if len(traj) < 4:
        return False
    ntail = max(3, int(math.ceil(len(traj) * tail_frac)))
    tail = traj[-ntail:]
    if tol_pts is None:
        REL_TOL = 0.05
        tol_pts = max(3, math.ceil(REL_TOL * median(tail)))
    return max(tail) - min(tail) <= tol_pts


def main():
    for case_dir in sys.argv[1:]:
        traj = n_trajectory(case_dir)
        if len(traj) < 4:
            status = 'TOO FEW OUTPUTS (%d)' % len(traj)
        else:
            status = 'SETTLED' if grid_settled(traj) else 'NOT SETTLED'
        print('%-60s N: %s -> %s (%d outputs)  %s' %
              (case_dir, traj[0] if traj else '-',
               traj[-1] if traj else '-', len(traj), status))


if __name__ == '__main__':
    main()
