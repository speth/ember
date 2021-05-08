#!/usr/bin/env python
"""
Outwardly-propagating cylindrical geometry for a strained lean methane flame.
The converged axial velocity profile is plotted. The stagnation point is
located at r=0.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_cylindrical_outward'

conf = Config(
    Paths(outputDir=output,
          # logFile='ex_cylindrical_outward.log'
          ),
    Chemistry(mechanismFile='gri30.yaml'),
    General(flameGeometry='cylindrical',
            unburnedLeft=False,
            fixedLeftLocation=True,
            nThreads=4),
    InitialCondition(fuel='CH4:0.5, H2:0.5',
                     equivalenceRatio=0.60,
                     xLeft=0.0,
                     xRight=0.005),
    StrainParameters(initial=500,
                     final=500),
    TerminationCondition(tEnd=10,
                         measurement='dTdt'),
    Times(profileStepInterval=10,
          regridStepInterval=10),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.V / struct.rho)
    plt.xlabel('Position [m]')
    plt.ylabel('Axial Velocity [m/s]')
    plt.savefig(output + '/FinalAxialVelocity.png')
    plt.close()
