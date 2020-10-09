import os

from ._ember import writelog
from . import utils
import numpy as np

class OutputFile(object):
    def __init__(self, filename):
        self.filename = filename

    def __enter__(self):
        if self.filename.endswith('.h5'):
            import h5py
            self.data = h5py.File(self.filename, mode='a')
            return self.data
        elif self.filename.endswith('.npz'):
            self.data = {}
            return self.data
        else:
            raise Exception("Unknown output file format for file '{0}'."
                " Expected one of: ('h5', 'npz')".format(filename))

    def __exit__(self, exc_type, exc_val, exc_tb):
        dirname = os.path.dirname(self.filename)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        if self.filename.endswith('h5'):
            self.data.close()
        elif self.filename.endswith('npz'):
            np.savez_compressed(self.filename, **self.data)


class TimeSeriesWriter(object):
    def __init__(self, solver, options):
        self.solver = solver
        self.options = options

        self.t = []
        self.dt = []
        self.Q = []
        self.Sc = []
        self.xFlame = []
        self.a = []
        self.dadt = []

    def __call__(self, name, flag=1):
        if not self.t or self.solver.tNow != self.t[-1]:
            self.t.append(self.solver.tNow)
            self.dt.append(self.solver.dt)
            self.Q.append(self.solver.heatReleaseRate)
            self.Sc.append(self.solver.consumptionSpeed)
            self.xFlame.append(self.solver.flamePosition)
            self.a.append(self.solver.a)
            self.dadt.append(self.solver.dadt)

        if flag:
            filename = '{}/{}.{}'.format(self.options.paths.outputDir, name,
                                         self.options.outputFiles.fileExtension)
            if os.path.exists(filename):
                os.remove(filename)

            with OutputFile(filename) as data:
                data['t'] = self.t
                data['dt'] = self.dt
                data['Q'] = self.Q
                data['Sc'] = self.Sc
                data['xFlame'] = self.xFlame
                data['a'] = self.a
                data['dadt'] = self.dadt


class StateWriter(object):
    def __init__(self, solver, options):
        self.solver = solver
        self.options = options
        self.fileNumber = self.options.outputFiles.firstFileNumber - 1

    def write(self, datafile, keys):
        for k in keys:
            datafile[k] = getattr(self.solver, k)

    def __call__(self, name, errorFile=False):
        if not name:
            # Determine the name of the output file (e.g. profXXXXXX.h5)
            self.fileNumber += 1
            filename = '{}/{}{:06d}.{}'.format(
                self.options.paths.outputDir,
                'prof' if not errorFile else 'error',
                self.fileNumber,
                self.options.outputFiles.fileExtension)
        else:
            filename = '{}/{}.{}'.format(self.options.paths.outputDir, name,
                                         self.options.outputFiles.fileExtension)

        if errorFile:
            writelog('Writing error output file: {}'.format(filename))
        else:
            writelog('Writing output file: {}'.format(filename))

        if os.path.exists(filename):
            os.remove(filename)

        with OutputFile(filename) as data:
            # Basic state information
            data['t'] = self.solver.tNow
            data['P'] = float(self.options.initialCondition.pressure)
            data['fileNumber'] = self.fileNumber
            self.write(data, ['x', 'T', 'U', 'Y', 'V', 'gridAlpha', 'a', 'dadt'])

            # extended information
            if self.options.outputFiles.heatReleaseRate or errorFile:
                data['q'] = self.solver.qDot
                data['rho'] = self.solver.rho

            # Individual terms in the governing equations
            if self.options.outputFiles.timeDerivatives or errorFile:
                self.write(data, ['dUdtDiff', 'dUdtConv', 'dUdtProd',
                                  'dTdtDiff', 'dTdtConv', 'dTdtProd', 'dTdtCross',
                                  'dYdtDiff', 'dYdtConv', 'dYdtProd', 'dYdtCross',
                                  'dWdt', 'drhodt'])

            if self.options.outputFiles.auxiliaryVariables or errorFile:
                self.write(data, ['sumcpj', 'Tleft', 'Yleft', 'dWdx', 'dTdx',
                    'splitConstDiff', 'splitConstConv', 'splitConstProd'])

            if self.options.outputFiles.extraVariables or errorFile:
                # These variables can be recomputed from the state variables
                self.write(data, ['wdot', 'rhoD', 'cp', 'mu', 'k', 'Wmx', 'W',
                                  'cfp', 'cf', 'cfm', 'hh', 'rphalf',
                                  'jFick', 'jSoret', 'jCorr'])
