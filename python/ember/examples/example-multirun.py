#!/usr/bin/python
from ember import *
import multiprocessing

conf = Config(
    General(nThreads=multiprocessing.cpu_count()),
    Chemistry(kineticsModel='interp',
              transportModel='Approx'),
    Grid(addPointCount=6),
    Paths(outputDir='run/multirunTest',
          logFile='out-multirunTest.txt'),
    InitialCondition(fuel='CH4:1.0',
                     equivalenceRatio=0.75),
    StrainParameters(rates=[4800,2400,1200,600,300]))

if __name__ == '__main__':
    conf.run()
