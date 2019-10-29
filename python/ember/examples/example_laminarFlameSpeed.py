#!/usr/bin/python
"""
An unstrained flame converges to a stable propagation rate for a slightly
lean methane/air mixture.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_lfs'

conf = Config(
    Paths(outputDir=output),
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

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.V / struct.rho)
    plt.xlabel('Position [m]')
    plt.ylabel('Axial Velocity [m/s]')
    plt.twinx()
    plt.plot(struct.x, struct.T, 'r--')
    plt.tick_params( colors='r')
    plt.ylabel('Temperature (K)', color='r')
    plt.savefig(output+'/FinalProfiles.png')
    plt.close()
