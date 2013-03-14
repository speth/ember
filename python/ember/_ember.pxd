from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool as cbool
from cpython import bool as pybool

import numpy as np
cimport numpy as np

cdef extern from "mathUtils.h":
    cdef cppclass CxxEigenVec "dvec":
        pass

    cdef cppclass CxxEigenMatrix "dmatrix":
        pass


cdef extern from "readConfig.h":
    cdef cppclass CxxConfigOptions "ConfigOptions":
        pass

cdef extern from "flameSolver.h":
    cdef cppclass CxxFlameSolver "FlameSolver":
        pass


cdef class ConfigOptions:
    cdef CxxConfigOptions* opts

cdef class FlameSolver:
    cdef CxxFlameSolver* solver
