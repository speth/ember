from _ember import *
import _ember

class FlameSolver(_ember.FlameSolver):
    def __init__(self, options):
        _ember.FlameSolver.__init__(self)
        self.setOptions(_ember._ConfigOptions(options))
