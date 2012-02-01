"""
Simple implementations of BDF integration to generate comparisons
for unit testing.
"""

import numpy as np

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

def main():
    dt = 0.2
    y0 = np.array([0, 0.5, 2.0, 1.0, 0])
    k = np.array([0, 0, 0, 0.2, 0.4])
    A = np.array([[-2, 1, 0, 0, 0],
                  [1, -2, 1, 0, 0],
                  [0, 1, -2, 1, 0],
                  [0, 0, 1, -2, 1],
                  [0, 0, 0, 1, -2]], dtype=float)
    sys = BDFIntegrator(dt, y0, A, k)

    Y = [y0.copy()]
    for i in range(6):
        print '{' + ', '.join(['%16.14f' % yi for yi in sys.y]) + '}'
        sys.step()
        #print 't = %14.14f' % ((i+1) * dt)
        Y.append(sys.y.copy())

    return Y


if __name__ == '__main__':
    main()
