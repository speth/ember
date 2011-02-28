import numpy as np
import Cantera as ct
import h5py

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
