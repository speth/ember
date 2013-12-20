import os.path as _path

import cantera
from _ember import *
import _ember
from input import *
import utils

__version__ = '1.3.dev'

# Add Ember's data file directory to Cantera's search path. Because the Python
# module is statically linked to Cantera, this needs to be done separately for
# each of the two copies of the Cantera library that have been loaded.
_datapath = _path.join(_path.dirname(_path.abspath(__file__)), 'data')
_ember.addCanteraDirectory(_datapath)
cantera.add_directory(_datapath)
