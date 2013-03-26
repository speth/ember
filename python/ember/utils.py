import numpy as np
import Cantera as ct
import h5py
import os
import sys
import _ember
import time

class Struct(object):
    """
    A dictionary-like data structure where fields are accessible as both
    attributes and dictionary keys::

        >>> s = Struct()
        >>> s['foo'] = 6
        >>> s.foo
        6
        >>> s.bar = 'x'
        >>> 'bar' in s:
        True
        >>> s['bar']
        'x'
        >>> s.keys()
        ['foo', 'bar']

    Valid methods of initialization, equivalent to the above:

        >>> s = Struct(foo=6, bar='x')
        >>> s = Struct({'foo': 6, 'bar': 'x'})

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
    Like :class:`Struct`, but converts HDF5 structs to numpy arrays.
    """
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception("File not found: " + filename)
        data = h5py.File(filename, mode='r') # We don't need to write to the file
        for key in data:
            self[key] = data[key].value
        data.close() # I don't know if this is necessary, but it can't hurt


def load(filename):
    """
    Generate an :class:`HDFStruct` object from a saved profile.
    """
    return HDFStruct(filename)


def get_qdot(gas, profile, pressure=101325):
    """
    Calculate the heat release rate along the flame coordinate.

    :param gas:
        A Cantera `Solution` object made using the same mechanism file used
        for the simulation.
    :param profile:
        An :class:`.HDFStruct` object created by loading a `profNNNNNN.h5` file.
    :param pressure:
        The pressure at which the simulation was run.
    """
    q = []
    for i in range(len(profile.T)):
        gas.set(P=pressure, T=profile.T[i], Y=profile.Y[:,i])
        hk = gas.enthalpies_RT() * profile.T[i] * ct.GasConstant
        wDot = gas.netProductionRates()
        q.append(-np.dot(wDot,hk))

    return np.array(q)


def expandProfile(prof, gas):
    """
    Reconstruct derived data associated with a flame profile.

    :param prof:
        An :class:`.HDFStruct` object created by loading a `profNNNNNN.h5` file
    :param gas:
        A Cantera `Solution` object made using the same mechanism file used
        for the simulation.

    Arrays which are reconstructed:

        * grid properties: *hh*, *cfp*, *cf*, *cfm*, *rphalf*, *dlj*
        * thermodynamic properties: *rho*, *cp*, *Wmx*, *W*
        * kinetic properties: *wdot*, *q*
        * transport properties: *rhoD*, *k*, *mu*, *Dkt*, *jFick*, *jSoret*, *jCorr*
        * other: *X* (mole fractions)
    """
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
    prof.X = np.zeros((K,N))

    for j in range(N):
        gas.set(T=prof.T[j], Y=prof.Y[:,j], P=P)
        prof.rho[j] = gas.density()
        hk = gas.enthalpies_RT() * prof.T[j] * ct.GasConstant
        wdot = gas.netProductionRates()
        prof.wdot[:,j] = wdot
        prof.q[j] = -np.dot(wdot,hk)
        prof.X[:,j] = gas.moleFractions()

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


def calculateReactantMixture(gas, fuel, oxidizer, equivalenceRatio):
    gas.set(X=fuel, T=300.0, P=101325)
    Xf = gas.moleFractions()
    gas.set(X=oxidizer, T=300.0, P=101325)
    Xo = gas.moleFractions()

    nC = np.array([gas.nAtoms(k, 'C') for k in range(gas.nSpecies())])
    nO = np.array([gas.nAtoms(k, 'O') for k in range(gas.nSpecies())])
    nH = np.array([gas.nAtoms(k, 'H') for k in range(gas.nSpecies())])

    Cf = (nC * Xf).sum()
    Co = (nC * Xo).sum()
    Of = (nO * Xf).sum()
    Oo = (nO * Xo).sum()
    Hf = (nH * Xf).sum()
    Ho = (nH * Xo).sum()

    stoichAirFuelRatio = - (Of - 2*Cf - Hf/2.0) / (Oo - 2*Co - Ho/2.0)
    Xr = Xf * equivalenceRatio + stoichAirFuelRatio * Xo
    Xr /= Xr.sum()

    return Xr


def smooth(v):
    v[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]
    return v
