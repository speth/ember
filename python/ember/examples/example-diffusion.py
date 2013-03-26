#!/usr/bin/python
from ember import *

conf = Config(
    Paths(outputDir='run/diffusionTest',
          logFile='out-diffusionTest.txt'),
    InitialCondition(flameType='diffusion',
                     fuel='CH4:1.0, N2:2.0',
                     Tfuel=600,
                     Toxidizer=600,
                     xLeft=-0.004,
                     xRight=0.004,
                     centerWidth=0.002,
                     slopeWidth=0.001),
    StrainParameters(initial=100,
                     final=100),
    Times(globalTimestep=2e-6,
          profileStepInterval=20),
    TerminationCondition(tEnd=0.010))

if __name__ == '__main__':
    conf.run()
