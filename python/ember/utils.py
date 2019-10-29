import numpy as np
import cantera as ct
import os
import sys
from . import _ember
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
        import h5py
        if not os.path.exists(filename):
            raise Exception("File not found: " + filename)
        data = h5py.File(filename, mode='r') # We don't need to write to the file
        for key in data:
            self[key] = data[key][()]
        data.close() # I don't know if this is necessary, but it can't hurt


class NpzStruct(Struct):
    """
    Like :class:`Struct`, but loads data from NumPy 'npz' data files
    """
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise Exception("File not found: " + filename)
        data = np.load(filename)
        for k,v in data.items():
            self[k] = v


def load(filename):
    """
    Generate an :class:`Struct` object from a saved profile.
    """
    if filename.endswith('h5'):
        return HDFStruct(filename)
    elif filename.endswith('npz'):
        return NpzStruct(filename)


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
        gas.TPY = profile.T[i], pressure, profile.Y[:,i]
        q.append(-np.dot(gas.net_production_rates, gas.partial_molar_enthalpies))

    return np.array(q)


def expandProfile(prof, gas, diffusion=True, reaction_rates=True):
    """
    Reconstruct derived data associated with a flame profile.

    :param prof:
        An :class:`.HDFStruct` object created by loading a `profNNNNNN.h5` file
    :param gas:
        A Cantera `Solution` object made using the same mechanism file used
        for the simulation.
    :param diffusion:
        Set to `False` to disable calculating diffusion properties (which can
        be slow for large mechanisms)
    :param reaction_rates:
        Set to `False` to disable detailed reaction rates (creation /
        destruction / rates-of-progress)

    Arrays which are reconstructed:

        * grid properties: *hh*, *cfp*, *cf*, *cfm*, *rphalf*, *dlj*
        * thermodynamic properties: *rho*, *cp*, *Wmx*, *W*
        * kinetic properties: *wdot*, *q*, *creation_rates*, *destruction_rates*,
          *forward_rates_of_progress*, *reverse_rates_of_progress*,
          *net_rates_of_progress*
        * transport properties: *rhoD*, *k*, *mu*, *Dkt*, *jFick*, *jSoret*, *jCorr*
        * other: *X* (mole fractions)
    """
    N = len(prof.x)
    I = gas.n_reactions

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
    K = gas.n_species

    try:
        P = prof.P
    except AttributeError:
        P = 101325

    if diffusion:
        prof.rhoD = np.zeros((K,N))
        prof.Dkt = np.zeros((K,N))
        prof.jFick = np.zeros((K,N))
        prof.jSoret = np.zeros((K,N))
        prof.jCorr = np.zeros(N)

    prof.rho = np.zeros(N)
    prof.wdot = np.zeros((K,N))
    prof.q = np.zeros(N)
    prof.k = np.zeros(N)
    prof.cp = np.zeros(N)
    prof.mu = np.zeros(N)
    prof.Wmx = np.zeros(N)
    prof.X = np.zeros((K,N))

    if reaction_rates:
        prof.creation_rates = np.zeros((K,N))
        prof.destruction_rates = np.zeros((K,N))
        prof.forward_rates_of_progress = np.zeros((I,N))
        prof.reverse_rates_of_progress = np.zeros((I,N))
        prof.net_rates_of_progress = np.zeros((I,N))

    for j in range(N):
        gas.TPY = prof.T[j], P, prof.Y[:,j]
        prof.rho[j] = gas.density
        wdot = gas.net_production_rates
        prof.wdot[:,j] = wdot
        prof.q[j] = -np.dot(wdot, gas.partial_molar_enthalpies)
        prof.creation_rates[:,j] = gas.creation_rates
        prof.destruction_rates[:,j] = gas.destruction_rates
        prof.forward_rates_of_progress[:,j] = gas.forward_rates_of_progress
        prof.reverse_rates_of_progress[:,j] = gas.reverse_rates_of_progress
        prof.net_rates_of_progress[:,j] = gas.net_rates_of_progress
        prof.X[:,j] = gas.X

        prof.k[j] = gas.thermal_conductivity
        prof.cp[j] = gas.cp_mass
        prof.mu[j] = gas.viscosity
        prof.Wmx[j] = gas.mean_molecular_weight

        if diffusion:
            Dbin = gas.binary_diff_coeffs
            prof.Dkt[:,j] = gas.thermal_diff_coeffs

            eps = 1e-15;
            for k in range(K):
                X = gas.X
                Y = gas.Y
                sum1 = sum(X[i]/Dbin[k,i] for i in range(K) if i != k)
                sum2 = sum((Y[i]+eps/K)/Dbin[k,i] for i in range(K) if i != k)
                prof.rhoD[k,j] = prof.rho[j]/(sum1 + X[k]/(1+eps-Y[k])*sum2)

    if diffusion:
        for j in range(1, N-1):
            for k in range(K):
                prof.jFick[k,j] = -0.5 * ((prof.rhoD[k,j] + prof.rhoD[k,j+1]) *
                                          ((prof.Y[k,j+1]-prof.Y[k,j])/prof.hh[j]))
                prof.jSoret[k,j] = -0.5 * ((prof.Dkt[k,j]/prof.T[j] +
                                            prof.Dkt[k,j+1]/prof.T[j+1]) *
                                            (prof.T[j+1]-prof.T[j])/prof.hh[j])
                prof.jCorr[j] -= prof.jFick[k,j] + prof.jSoret[k,j]

    prof.W = gas.molecular_weights


def smooth(v):
    v[1:-1] = 0.25 * v[:-2] + 0.5 * v[1:-1] + 0.25 * v[2:]
    return v
