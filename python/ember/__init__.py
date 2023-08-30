import os.path as _path

import cantera
from ._ember import *
from . import _ember
from .input import *
from .output import *
from . import utils

__version__ = '1.5.0b1'

# Add Ember's data file directory to Cantera's search path. Because the Python
# module is statically linked to Cantera, this needs to be done separately for
# each of the two copies of the Cantera library that have been loaded.
_datapath = _path.join(_path.dirname(_path.abspath(__file__)), 'data')
_ember.addCanteraDirectory(_datapath)
cantera.add_directory(_datapath)
