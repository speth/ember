#!/usr/bin/python
"""
Planar opposed twin flame geometry for strained lean methane flames. The converged
axial velocity profile is plotted. Note that symmetry of the domain is leveraged
to simplify the computation. The stagnation plane is located as x=0.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_twin'

conf = Config(
    Paths(outputDir=output,
          # logFile='ex_twin.log'
          ),
    # Chemistry(mechanismFile='gri30.yaml'),
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

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.V / struct.rho)
    plt.xlabel('Position [m]')
    plt.ylabel('Axial Velocity [m/s]')
    plt.savefig(output + '/FinalAxialVelocity.png')
    plt.close()
