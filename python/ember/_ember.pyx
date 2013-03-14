#cython: embedsignature=True
#distutils: language = c++

import numpy as np
cimport numpy as np

from _ember cimport *

cdef class ConfigOptions:
    def __cinit__(self, *args, **kwargs):
        self.opts = new CxxConfigOptions()

    def __dealloc__(self):
        del self.opts

    def apply_options(self):
        self.opts.outputDir = self.paths.outputDir


cdef class FlameSolver:
    def __cinit__(self, *args, **kwargs):
        self.solver = new CxxFlameSolver()

    def __dealloc__(self):
        del self.solver
