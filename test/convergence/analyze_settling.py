"""
Grid-settling analysis for errTol-ladder runs: extracts the grid-size
trajectory from the profNNNNNN.h5 outputs of a case work directory and
tests whether the grid reached a plateau (the P2.4 failure signature is a
monotonically growing, never-settling grid).

Usage: python analyze_settling.py <workdir> [<workdir> ...]
"""
import sys
import os
import glob
import h5py


def n_trajectory(case_dir):
    """Grid point count for each profile output, in time order."""
    files = sorted(glob.glob(os.path.join(case_dir, 'prof*.h5')))
    traj = []
    for f in files:
        with h5py.File(f, 'r') as h:
            traj.append(int(h['x'].shape[0]))
    return traj


def grid_settled(traj, tail_frac=0.25, tol_pts=3):
    """True if N varies by <= tol_pts over the last tail_frac of outputs."""
    if len(traj) < 8:
        return False
    tail = traj[int(len(traj) * (1 - tail_frac)):]
    return max(tail) - min(tail) <= tol_pts


def main():
    for case_dir in sys.argv[1:]:
        traj = n_trajectory(case_dir)
        status = 'SETTLED' if grid_settled(traj) else 'NOT SETTLED'
        print('%-60s N: %s -> %s (%d outputs)  %s' %
              (case_dir, traj[0] if traj else '-',
               traj[-1] if traj else '-', len(traj), status))


if __name__ == '__main__':
    main()
