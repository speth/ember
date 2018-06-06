#!/usr/bin/python
#$ -pe singlenode 8
#$ -l normal
#$ -l h_rt=120:00:00
#$ -j y
#$ -o job-output-example-parallel.txt

"""
Example of running flames multiple configuration conditions (equivalence ratios,
temperatures, dilutions, etc). Default flame: Lean methane opposing hot products
in planar opposed jet configuration.
"""

from ember import *

outputDir = 'run/ex_batch/ex-parallel-phi%4.2f'
logFile = 'run/ex_batch/out-ex-parallel-phi%4.2f'

def start(phi):
    conf = Config(
        General(nThreads=8),
        Paths(outputDir=outputDir % phi,
              logFile=logFile % phi
              ),
        InitialCondition(equivalenceRatio=phi),
        StrainParameters(initial=50,
                         final=50),
        Times(profileStepInterval=1000000,
              profileTimeInterval=5e-3))

    conf.run()

if __name__ == '__main__':
    for phi in [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95]:
        start(phi)
    print('All done!')
