import numpy as np
import Cantera as ct
import h5py
import os
import sys
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
        data = h5py.File(filename, mode='r') # We don't need to write to the file
        for key in data:
            self[key] = data[key].value
        data.close() # I don't know if this is necessary, but it can't hurt


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


def expandProfile(prof, gas):
    N = len(prof.x)

    # Grid properties
    try:
        gridAlpha = prof.gridAlpha
    except AttributeError:
        gridAlpha = 0

    prof.hh = np.zeros(N)
    prof.cfp = np.zeros(N)
    prof.cf = np.zeros(N)
    prof.cfm = np.zeros(N)
    prof.rphalf = np.zeros(N)
    prof.dlj = np.zeros(N)

    for j in range(N-1):
        prof.hh[j] = prof.x[j+1] - prof.x[j]
        prof.rphalf[j] = (0.5 * (prof.x[j]+prof.x[j+1]))**prof.gridAlpha

    hh = prof.hh
    for j in range(1, N-1):
        prof.cfp[j] = hh[j-1]/(hh[j]*(hh[j]+hh[j-1]))
        prof.cf[j] = (hh[j]-hh[j-1])/(hh[j]*hh[j-1])
        prof.cfm[j] = -hh[j]/(hh[j-1]*(hh[j]+hh[j-1]))
        prof.dlj[j] = 0.5 * (prof.x[j+1] - prof.x[j-1])

    # Thermodynamic / Transport / Kinetic properties
    K = gas.nSpecies()

    try:
        P = prof.P
    except AttributeError:
        P = 101325

    prof.rho = np.zeros(N)
    prof.wdot = np.zeros((K,N))
    prof.q = np.zeros(N)
    prof.rhoD = np.zeros((K,N))
    prof.k = np.zeros(N)
    prof.cp = np.zeros(N)
    prof.mu = np.zeros(N)
    prof.Wmx = np.zeros(N)
    prof.Dkt = np.zeros((K,N))
    prof.jFick = np.zeros((K,N))
    prof.jSoret = np.zeros((K,N))
    prof.jCorr = np.zeros(N)

    for j in range(N):
        gas.set(T=prof.T[j], Y=prof.Y[:,j], P=P)
        prof.rho[j] = gas.density()
        hk = gas.enthalpies_RT() * prof.T[j] * ct.GasConstant
        wdot = gas.netProductionRates()
        prof.wdot[:,j] = wdot
        prof.q[j] = -np.dot(wdot,hk)

        Dbin = gas.binaryDiffCoeffs()
        prof.Dkt[:,j] = gas.thermalDiffCoeffs()

        eps = 1e-15;
        for k in range(K):
            X = gas.moleFractions()
            Y = gas.massFractions()
            sum1 = sum(X[i]/Dbin[k,i] for i in range(K) if i != k)
            sum2 = sum((Y[i]+eps/K)/Dbin[k,i] for i in range(K) if i != k)
            prof.rhoD[k,j] = prof.rho[j]/(sum1 + X[k]/(1+eps-Y[k])*sum2)

        prof.k[j] = gas.thermalConductivity()
        prof.cp[j] = gas.cp_mass()
        prof.mu[j] = gas.viscosity()
        prof.Wmx[j] = gas.meanMolecularWeight()

    for j in range(1, N-1):
        for k in range(K):
            prof.jFick[k,j] = -0.5 * ((prof.rhoD[k,j] + prof.rhoD[k,j+1]) *
                                      ((prof.Y[k,j+1]-prof.Y[k,j])/prof.hh[j]))
            prof.jSoret[k,j] = -0.5 * ((prof.Dkt[k,j]/prof.T[j] +
                                        prof.Dkt[k,j+1]/prof.T[j+1]) *
                                        (prof.T[j+1]-prof.T[j])/prof.hh[j])
            prof.jCorr[j] -= prof.jFick[k,j] + prof.jSoret[k,j]

    prof.W = gas.molecularWeights()

def run(conf):
    # Validate the configuration and exit
    if len(sys.argv) > 1 and sys.argv[1].lower() == 'validate':
        conf.validate()
        return

    if not os.path.isdir(conf.paths.outputDir):
        os.makedirs(conf.paths.outputDir, 0755)
    confOutPath = os.path.join(conf.paths.outputDir, 'config')
    if (os.path.exists(confOutPath)):
        os.unlink(confOutPath)
    confOut = file(confOutPath, 'w')
    confOut.write(conf.stringify())

    solver = _pyro.FlameSolver(conf)
    solver.initialize()
    solver.run()


def multirun(conf):
    # Validate the configuration and exit
    if len(sys.argv) > 1 and sys.argv[1].lower() == 'validate':
        conf.validate()
        return

    strainRates = conf.strainParameters.rates
    if not strainRates:
        print 'No strain rate list specified'
        return

    conf.strainParameters.rates = None
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
        configFile = 'conf_eps%04i' % a

        restartPath = os.path.join(conf.paths.outputDir, restartFile)
        historyPath = os.path.join(conf.paths.outputDir, historyFile)
        configPath = os.path.join(conf.paths.outputDir, configFile)

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
            confOut = file(configPath, 'w')
            confOut.write(conf.stringify())

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
            solver.writeTimeseriesFile(historyFile)
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
        strainRates, Q, Sc, xFlame = map(list, zip(*sorted(zip(strainRates, Q, Sc, xFlame))))

        integralFile = os.path.join(conf.paths.outputDir, "integral.h5")
        if os.path.exists(integralFile):
            os.unlink(integralFile)
        data = h5py.File(integralFile)
        data['a'] = strainRates
        data['Q'] = Q
        data['Sc'] = Sc
        data['xFlame'] = xFlame
        data.close()
