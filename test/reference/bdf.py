#!/usr/bin/env python
"""
Simple implementations of BDF integration to generate comparisons
for unit testing.
"""

import numpy as np
from scipy import linalg

class BDFIntegrator(object):
    def __init__(self, h, y, A, k):
        """
        Integrator for y' = Ay + k.
        Uses second-order BDF with step size h.
        First timestep is taken as first-order BDF with step size h/8
        """

        self.stepCount = 0
        self.h = h
        self.y = y
        self.A = A
        self.k = k
        self.N = len(self.y)
        self.nSub = 8 # number of substeps

    def step(self):
        if self.stepCount == 0:
            self.yprev = self.y.copy()
            for j in range(self.nSub):
                # Apply first-order BDF: y_n = y_(n-1) + h * y'_n
                # by solving (I-A) y_n = y_(n-1) + h*k
                h = self.h / self.nSub
                M = np.eye(self.N) - self.A * h
                b = self.y +  h * self.k
                self.y = np.linalg.solve(M, b)

        else:
            # Apply second-order BDF: y_n = 4/3*y_(n-1) - 1/3*y_(n-2) + 2/3*h * y'_n
            # by solving (I-2/3*h*A) y_n = 4/3*y_(n-1) - 1/3*y_(n-2) + (2/3)*h*k
            self.yprev2 = self.yprev.copy()
            self.yprev = self.y.copy()
            M = np.eye(self.N) - 2.0/3.0 * self.h * self.A
            b = 4.0/3.0 * self.yprev - 1.0/3.0 * self.yprev2 + 2.0/3.0 * self.h * self.k
            self.y = np.linalg.solve(M, b)

        self.stepCount += 1

class ExactIntegrator(object):
    def __init__(self, y, A, k):
        """
        Integrator for y' = Ay + k.
        Uses matrix exponential to give exact solution.
        """

        self.y = y
        self.A = A
        self.k = k
        self.N = len(self.y)

        self.Ainv_k = np.dot(linalg.inv(A), k)

    def __call__(self, t):
        return np.dot(linalg.expm(self.A * t), self.y + self.Ainv_k) - self.Ainv_k


def getBdfSolutions(dt, tf, y0, A, k):
    sys = BDFIntegrator(dt, y0, A, k)
    nSteps = int(round(tf/dt))
    Y = [y0.copy()]
    for i in range(nSteps):
        sys.step()
        Y.append(sys.y.copy())
    return Y


def main():
    y0 = np.array([0, 0.5, 2.0, 1.0, 0])
    k = np.array([0, 0, 0, 0.2, 0.4])
    A = np.array([[-2, 1, 0, 0, 0],
                  [1, -2, 1, 0, 0],
                  [0, 1, -2, 1, 0],
                  [0, 0, 1, -2, 1],
                  [0, 0, 0, 1, -2]], dtype=float)

    Y1 = getBdfSolutions(0.20, 1.0, y0, A, k)
    Y2 = getBdfSolutions(0.05, 1.0, y0, A, k)
    Y3 = getBdfSolutions(0.0125, 1.0, y0, A, k)
    exact = ExactIntegrator(y0, A, k)

    print('Temporal BDF2 solutions, dt = 0.20:')
    for y in Y1:
        print('{' + ', '.join(['%16.14f' % yi for yi in y]) + '}')
    print()

    print('BDF2 solution at t = 1.0 (dt = 0.05)')
    print('{' + ', '.join(['%16.14f' % yi for yi in Y2[-1]]) + '}\n')

    print('BDF2 solution at t = 1.0 (dt = 0.0125):')
    print('{' + ', '.join(['%16.14f' % yi for yi in Y3[-1]]) + '}\n')

    print('Exact solution at t = 1.0:')
    print('{' + ', '.join(['%16.14f' % yi for yi in exact(1.0)]) + '}\n')

    print('Apparent order of convergence: dt = 0.05 vs dt 0.0125:')
    print(np.log((np.linalg.norm(Y2[-1] - exact(1.0))
                  / np.linalg.norm(Y3[-1] - exact(1.0)))) / np.log(0.05/0.0125))


if __name__ == '__main__':
    main()
