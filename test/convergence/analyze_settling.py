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
import h5py


def n_trajectory(case_dir):
    """Grid point count for each numbered profile output, in time order."""
    files = sorted(glob.glob(os.path.join(case_dir, 'prof[0-9]*.h5')))
    traj = []
    for f in files:
        with h5py.File(f, 'r') as h:
            traj.append(int(h['x'].shape[0]))
    return traj


def grid_settled(traj, tail_frac=0.25, tol_pts=3):
    """True if the trajectory ends in a plateau: N varies by <= tol_pts over
    the tail window (the last tail_frac of outputs, at least 3 samples).
    Fewer than 4 samples is not judgeable and returns False; callers should
    treat that case as "inspect the run log", not as a settling failure
    (main() reports it as TOO FEW OUTPUTS).
    """
    if len(traj) < 4:
        return False
    ntail = max(3, int(math.ceil(len(traj) * tail_frac)))
    tail = traj[-ntail:]
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
