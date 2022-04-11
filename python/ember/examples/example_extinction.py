#!/usr/bin/python
"""
A steady flame is established at a starting strain rate far from exctinction.
Then the strain rate parameter (a) is systematically increased until a steady
flame can no longer be attained. The progression to extinction is summarized
in a plot of maximum temperature vs strain rate. Twin disc flame geometry.
"""

from ember import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing
import shutil

output = 'run/ex_extinction'

conf = Config(
    Paths(outputDir=output),
    # Chemistry(mechanismFile='gri30.yaml'),
    General(twinFlame=True,
            flameGeometry='disc',
            unburnedLeft=False,
            fixedLeftLocation=True,
            nThreads=multiprocessing.cpu_count()//2,
            chemistryIntegrator='qss'  # default = 'qss' other option = 'cvode'
            ),
    InitialCondition(fuel='CH4:1.0',
                     oxidizer='N2:7.52, O2:2.0',
                     equivalenceRatio=0.7,
                     Tu=298.0,
                     xLeft=0.0,
                     xRight=0.020,
                     ),
    Grid(
        gridMin=1E-8,  # Default = 5=-7
        centerGridMin=1E-8,  # Default = 1E-4
    ),
    Extinction(
        method='step',  # default = 'step' other option = 'factor'
        initialStep=75.0,
        minStep=0.5,
        # initialFactor=1.05,
        # minFactor=1.0001,
        reductionFactor=0.4,
        cutoffTemp=1500.0,
        initialStrainRate=500,
    )
)

if __name__ == '__main__':

    if os.path.exists(output):
        shutil.rmtree(output)
    conf.runESR()

    data = np.genfromtxt(output+'/extProfile.csv', skip_header=1, delimiter=',')

    plt.figure()
    plt.xlabel('Strain Rate [1/s]')
    plt.ylabel('Max. Temp. [K]')
    plt.semilogx(data[:,0], data[:,1])
    plt.savefig(output+'/progression.png')
    # plt.show()
    plt.close()
