#!/usr/bin/python
"""
Similar to as is done experimentally (https://doi.org/10.1016/j.combustflame.2009.06.011),
a single disc flame opposing a cold inert is established at a given strain rate. Infinite
burner separation distance assumed. The converged axial velocity profile is plotted.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_single'
a = 300.0

conf = Config(
    Paths(outputDir=output),
    # Chemistry(mechanismFile='gri30.yaml'),
    General(twinFlame=False,
            flameGeometry= 'disc',
            nThreads=4),
    InitialCondition(fuel='CH4:1.0',
                     oxidizer = 'N2:3.76, O2:1.0',
                     equivalenceRatio=1.0,
                     counterflow = 'N2:1.0',
                     Tcounterflow = 300.0,
                     xLeft=-0.01,
                     xRight=0.01,
                     centerWidth=0.005,
                     slopeWidth=0.001,
                     ),
    StrainParameters(initial=a,
                     final=a),
)

if __name__ == '__main__':
    conf.run()

    struct = utils.load(output + '/profNow.h5')

    plt.figure()
    plt.plot(struct.x, struct.V / struct.rho)
    plt.xlabel('Position [m]')
    plt.ylabel('Axial Velocity [m/s]')
    plt.savefig(output+'/FinalAxialVelocity.png')
    plt.close()
