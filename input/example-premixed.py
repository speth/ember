#!/usr/bin/python
from pyro import *

conf = Config(
    Paths(outputDir='run/premixedTest',
          logFile='out-premixedTest.txt'),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.70),
    StrainParameters(initial=400,
                     final=400),
    TerminationCondition(tEnd=0.010))

run(conf)
