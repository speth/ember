#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run/laminar-methane'),
    InitialCondition(fuel='CH4:1.0',
                     oxidizer='O2:1, N2:3.76',
                     equivalenceRatio=0.9,
                     xLeft=0.0,
                     xRight=0.01),
    StrainParameters(initial=0,
                     final=0),
    General(fixedLeftLocation=True,
            fixedBurnedVal=False,
            nThreads=4),
    Grid(vtol=0.1, dvtol=0.15, gridMin=5e-6, gridMax=0.001),
    PositionControl(proportionalGain=2000, xInitial=0.005, xFinal=0.005),
    TerminationCondition(tolerance=1e-5),
    Times(profileStepInterval=50))

run(conf)
