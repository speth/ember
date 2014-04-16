#!/usr/bin/python
#$ -pe singlenode 8
#$ -l normal
#$ -l h_rt=120:00:00
#$ -j y
#$ -o job-output-example-parallel.txt

"""
This example demonstrates some features that are useful when running Ember in a
batch environment.

Working with Grid Engine
------------------------

The lines above starting with '#$' are directives for the Grid Engine batch
queueing system. Since this script is itself executable, it is not necessary to
write a separate job submission script. The various resource requests need to be
specified based on the configuration of the cluster you are using.

Running in parallel
-------------------

The option 'nThreads=8' means that Ember will use 8 processors simultaneously.
This mirrors the resource request to the queue manager for an 8 processor node,
'-pe singlenode 8'. Ember is written as a multithreaded parallel code, so all of
the processors allocated to it need to be on the same node.

Running several cases sequentially
----------------------------------

In this script, we define a function that solves a case at a particular
equivalence ratio, then call that function for each equivalence ratio of
interest. Python's string interpolation (using the '%' operator) is used to
generate consistent names for the output files. This method can be used to set
up parametric studies of various input variables.
"""

from ember import *

outputDir = 'run/example-parallel-phi%4.2f'
logFile = 'out-example-parallel-phi%4.2f'

def start(phi):
    conf = Config(
        General(nThreads=8),
        Paths(outputDir=outputDir % phi,
              logFile=logFile % phi),
        InitialCondition(equivalenceRatio=phi),
        StrainParameters(initial=50,
                         final=50),
        Times(profileStepInterval=1000000,
              profileTimeInterval=5e-3))

    conf.run()

if __name__ == '__main__':
    for phi in [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95]:
        start(phi)
    print 'All done!'
