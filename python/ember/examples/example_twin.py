#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run/twintest'),
    Chemistry(mechanismFile='gri30.xml'),
    General(twinFlame=True,
            unburnedLeft=False,
            nThreads=4),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.70,
                     xLeft=0.0,
                     xRight=0.01),
    StrainParameters(initial=100,
                     final=100),
    TerminationCondition(tEnd=10))

if __name__ == '__main__':
    conf.run()
