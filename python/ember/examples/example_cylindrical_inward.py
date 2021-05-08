#!/usr/bin/env python
"""
An inwardly-propagating, cylindrical, strained, lean flame.
The flame radius, defined by the centroid of the heat release
rate, is specified. The stagnation surface is a cylinder located
at a greater radius than the flame (on the burned side).

In this configuration, the curvature and stretching of the
flame can be varied independently.
"""

from ember import *
output = 'run/ex_cylindrical_inward'

conf = Config(
    Paths(outputDir=output,
          # logFile='ex_cylindrical_inward.log'
          ),
    Chemistry(mechanismFile='gri30.yaml'),
    General(flameGeometry='cylindrical',
            unburnedLeft=True,
            fixedLeftLocation=True,
            nThreads=4),
    InitialCondition(fuel='CH4:0.5, H2:0.5',
                     equivalenceRatio=0.60,
                     xLeft=0.0,
                     xRight=0.006),
    StrainParameters(initial=200,
                     final=200),
    PositionControl(xInitial=0.002,
                    xFinal=0.002),
    TerminationCondition(tEnd=10,
                         measurement='dTdt'),
    Times(profileStepInterval=10,
          regridStepInterval=10),
)

if __name__ == '__main__':
    conf.run()
