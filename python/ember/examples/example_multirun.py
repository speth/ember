#!/usr/bin/python
"""
Lean methane mixture opposing hot equilibrium products in the planar opposed jet
configuration. Demonstration of option to run a select set of strain rates in list
order, using the converged profile from the previous strain rate as the initial
profile for the following strain rate.
"""

from ember import *
import multiprocessing

conf = Config(
    General(nThreads=multiprocessing.cpu_count()),
    Chemistry(kineticsModel='interp',
              transportModel='Approx'),
    Grid(addPointCount=6),
    Paths(outputDir='run/ex_multirun'),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.75),
    StrainParameters(rates=[4800,2400,1200,600,300]))

if __name__ == '__main__':
    conf.run()
