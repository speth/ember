#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run/premixedTest'),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.70),
    StrainParameters(initial=400,
                     final=400),
    TerminationCondition(tEnd=0.010))

if __name__ == '__main__':
    conf.run()
