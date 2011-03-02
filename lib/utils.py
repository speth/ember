import numpy as np
import Cantera as ct
import h5py
import os
import _pyro
import time

class Struct(object):
    """
    data structure with fields accessible as
    attributes and dictionary keys
    """
    def __init__(self, *args, **kwargs):
        for arg in args:
            if hasattr(arg,'items'):
                for k,v in arg.items():
                    self[k] = v

        for k,v in kwargs.items():
            self[k] = v

    def __getitem__(self, key):
        return getattr(self, key, None)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __delitem__(self, key):
        delattr(self, key)

    def __contains__(self, key):
        return (key in self.__dict__)

    def __repr__(self):
        return repr(self.__dict__.keys())

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self):
        return self.__dict__.items()


class HDFStruct(Struct):
    """
    Like Struct, but converts HDF5 structs to numpy arrays
    """
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception("File not found: " + filename)
        data = h5py.File(filename)
        for key in data:
            self[key] = data[key].value


def load(filename):
    return HDFStruct(filename)


def get_qdot(gas, profile, pressure=101325):
    q = []
    for i in range(len(profile.T)):
        gas.set(P=pressure, T=profile.T[i], Y=profile.Y[:,i])
        hk = gas.enthalpies_RT() * profile.T[i] * ct.GasConstant
        wDot = gas.netProductionRates()
        q.append(-np.dot(wDot,hk))

    return np.array(q)


def multirun(conf):
    strainRates = conf.strainParameters.rates
    if not strainRates:
        print 'No strain rate list specified'
        return

    _logFile = file(conf.paths.logFile, 'w')
    def log(message):
        _logFile.write(message)
        _logFile.write('\n')
        _logFile.flush()

    conf.initialCondition.relativeRestartPath = False
    if not os.path.exists(conf.paths.outputDir):
        os.mkdir(conf.paths.outputDir, 0755)

    Q = []
    Sc = []
    xFlame = []

    for a in strainRates:
        restartFile = 'prof_eps%04i' % a
        historyFile = 'out_eps%04i' % a

        restartPath = os.path.join(conf.paths.outputDir, restartFile)
        historyPath = os.path.join(conf.paths.outputDir, historyFile)

        if os.path.exists(restartPath) and os.path.exists(historyPath):
            # If the output files already exist, we simply retrieve the
            # integral flame properties from the existing profiles and
            # advance to the next strain rate.

            log('Skipping run at strain rate a = %g'
                ' because the output file "%s" already exists.' % (a, restartFile))

            # Compute integral properties using points from the last half
            # of the termination-check period
            data = HDFStruct(historyPath)
            mask = data.t > data.t[-1] - 0.5*conf.terminationCondition.steadyPeriod
            if not any(mask):
                log('Warning: old data file did not contain data'
                    ' spanning the requested period.')
                mask = data.t > 0.5*data.t[-1]

            Q.append(np.mean(data.Q[mask]))
            Sc.append(np.mean(data.Sc[mask]))
            xFlame.append(np.mean(data.xFlame[mask]))
            del data

        else:
            # Data is not already present, so run the flame solver for this strain rate

            log('Beginning run at strain rate a = %g s^-1' % a)

            conf.strainParameters.initial = a
            conf.strainParameters.final = a
            conf.paths.logFile = os.path.join(conf.paths.outputDir, 'log-eps%04i.txt' % a)
            solver = _pyro.FlameSolver(conf)
            t1 = time.time()
            solver.initialize()
            solver.run()
            t2 = time.time()

            log('Completed run at strain rate a = %g s^-1' % a)
            log('Integration took %.1f seconds.' % (t2-t1))

            solver.writeStateFile(restartFile, False, False)
            tRun = np.array(solver.timeVector)
            QRun = np.array(solver.heatReleaseRate)
            ScRun = np.array(solver.consumptionSpeed)
            xFlameRun = np.array(solver.flamePosition)

            # Compute integral properties using points from the last half
            # of the termination-check period
            mask = tRun > tRun[-1] - 0.5*conf.terminationCondition.steadyPeriod
            Q.append(np.mean(QRun[mask]))
            Sc.append(np.mean(ScRun[mask]))
            xFlame.append(np.mean(xFlameRun[mask]))

        conf.initialCondition.restartFile = restartPath

        # Sort by strain rate:
        strainRates, Q, Sc, xFlame = zip(*sorted(zip(strainRates, Q, Sc, xFlame)))
        data = h5py.File(os.path.join(conf.paths.outputDir, "integral.h5"))
        data['a'] = strainRates
        data['Q'] = Q
        data['Sc'] = Sc
        data['xFlame'] = xFlame
        data.close()
