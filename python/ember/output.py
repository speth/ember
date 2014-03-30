import os

import h5py

from _ember import writelog

class TimeSeriesWriter(object):
    def __init__(self, solver, options):
        self.solver = solver
        self.options = options

    def __call__(self, name, flag):
        pass


class StateWriter(object):
    def __init__(self, solver, options):
        self.solver = solver
        self.options = options
        self.fileNumber = self.options.outputFiles.firstFileNumber - 1

    def write(self, datafile, keys):
        for k in keys:
            datafile[k] = getattr(self.solver, k)

    def __call__(self, name, errorFile):
        if not name:
            # Determine the name of the output file (e.g. profXXXXXX.h5)
            self.fileNumber += 1
            filename = '{}/{}{:06d}.h5'.format(
                self.options.paths.outputDir,
                'prof' if not errorFile else 'error',
                self.fileNumber)
        else:
            filename = '{}/{}.h5'.format(self.options.paths.outputDir, name)

        if errorFile:
            writelog('Writing error output file: {}'.format(filename))
        else:
            writelog('Writing output file: {}'.format(filename))

        if os.path.exists(filename):
            os.remove(filename)

        with h5py.File(filename) as data:
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
                self.write(data, ['sumcpj', 'Tleft', 'Yleft', 'dWdx', 'dTdx'])

            if self.options.outputFiles.extraVariables or errorFile:
                # These variables can be recomputed from the state variables
                self.write(data, ['wdot', 'rhoD', 'cp', 'mu', 'k', 'Wmx', 'W',
                                  'cfp', 'cf', 'cfm', 'hh', 'rphalf',
                                  'jFick', 'jSoret', 'jCorr'])
