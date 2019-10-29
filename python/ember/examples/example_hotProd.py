#!/usr/bin/python
"""
A lean methane jet opposes a jet of hot equilibrium products. Planar opposed
jet geometry. Final simulation velocity and temperature profiles plotted.
"""

from ember import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

output = 'run/ex_hotProd'
a = 400.0

conf = Config(
    Paths(outputDir=output),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.70),
    StrainParameters(initial=a,
                     final=a),
    TerminationCondition(tEnd=0.010))

if __name__ == '__main__':
    # conf.run()

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
