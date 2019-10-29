#!/usr/bin/python
"""
Hot methane diluted by nitrogen opposing hot air in the planar opposed jet geometry
for a low strain rate. The converged axial velocity profile is plotted.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_diffusion'

conf = Config(
    Paths(outputDir=output,
          # logFile='ex_diffusion.log'
          ),
    InitialCondition(flameType='diffusion',
                     fuel='CH4:1.0, N2:2.0',
                     oxidizer='N2:3.76, O2:1.0',
                     Tfuel=600,
                     Toxidizer=600,
                     xLeft=-0.004,
                     xRight=0.004,
                     centerWidth=0.002,
                     slopeWidth=0.001),
    StrainParameters(initial=100,
                     final=100),
    General(nThreads=2),
    Times(globalTimestep=1e-5,
          profileStepInterval=20),
    TerminationCondition(tEnd=0.010))

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    f, ax = plt.subplots()
    ax.plot(struct.x, struct.V / struct.rho)
    ax.set_xlabel('Position [m]')
    ax.set_ylabel('Axial Velocity [m/s]')
    f.savefig(output+'/FinalAxialVelocity.png')
