import unittest

from ember import *
import cantera as ct
from unittest import TestCase
import numpy as np

class TestPremixedStrained(TestCase):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/premixedStrained',
                  logFile='build/test/work/premixedStrained.txt'),
            General(nThreads=1),
            Chemistry(mechanismFile='h2o2.yaml'),
            InitialCondition(fuel='H2:1.0',
                             oxidizer='O2:1.0, AR:4.0',
                             equivalenceRatio=0.3),
            StrainParameters(initial=800, final=800),
            Grid(vtol=0.2, dvtol=0.3),
            Times(regridStepInterval=10),
            TerminationCondition(tEnd=0.01, measurement=None))

    def test_validates(self):
        self.assertTrue(self.conf.validate())

    def test_run(self):
        solver = self.conf.run()
        np.testing.assert_allclose(solver.tNow, 0.01, atol=3e-5)

        # Reactants side
        gas = ct.Solution('h2o2.yaml')
        gas.TPX = 300, 101325, 'H2:0.6, O2:1.0, AR:4.0'
        np.testing.assert_approx_equal(gas.T, solver.T[0])
        np.testing.assert_allclose(gas.Y, solver.Y[:,0], atol=1e-8)

        # Products side
        gas.equilibrate('HP')
        np.testing.assert_allclose(gas.Y, solver.Y[:,-1], atol=1e-8)
        np.testing.assert_approx_equal(gas.T, solver.T[-1], 5)

        # Since this is Le < 1
        self.assertGreater(max(solver.T), gas.T)


class TestTwinPremixedStrained(TestCase):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/twinPremixedStrained',
                  logFile='build/test/work/twinPremixedStrained.txt'),
            Chemistry(mechanismFile='h2o2.yaml'),
            InitialCondition(fuel='H2:1.0',
                             oxidizer='O2:1.0, AR:4.0',
                             xRight=0.003,
                             equivalenceRatio=0.3),
            General(nThreads=1, twinFlame=True, unburnedLeft=False),
            StrainParameters(initial=800, final=800),
            Grid(vtol=0.2, dvtol=0.3),
            Times(regridStepInterval=10),
            TerminationCondition(tEnd=0.01, measurement=None))

    def test_validates(self):
        self.assertTrue(self.conf.validate())

    def test_run(self):
        solver = self.conf.run()
        np.testing.assert_allclose(solver.tNow, 0.01, atol=3e-5)

        products = ct.Solution('h2o2.yaml')
        equil_prod = ct.Solution('h2o2.yaml')
        reactants = ct.Solution('h2o2.yaml')

        # Reactants side
        reactants.TPX = 300, 101325, 'H2:0.6, O2:1.0, AR:4.0'
        np.testing.assert_allclose(reactants.T, solver.T[-1])
        np.testing.assert_allclose(reactants.Y, solver.Y[:,-1], atol=1e-8)

        # Products side. Not necessarily equilibrium
        equil_prod.TPX = reactants.TPX
        equil_prod.equilibrate('HP')
        products.TPY = solver.T[0], 101325, solver.Y[:,0]
        self.assertTrue(equil_prod.T - 300 < products.T < equil_prod.T + 300)
        self.assertLess(products['H2'].X, 0.1*reactants['H2'].X)
        self.assertGreater(products['H2O'].X, 0.5*equil_prod['H2O'].X)


class TestDiffusion(TestCase):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/h2Diffusion',
                  logFile='build/test/work/h2Diffusion.txt'),
            General(nThreads=1),
            Chemistry(mechanismFile='h2o2.yaml'),
            InitialCondition(flameType='diffusion',
                             fuel='H2:1.0, AR:1.0',
                             oxidizer='O2:1.0, AR:4.0'),
            StrainParameters(initial=400, final=400),
            Grid(vtol=0.2, dvtol=0.3),
            Times(regridStepInterval=10),
            TerminationCondition(tEnd=0.01, measurement=None))

    def test_validates(self):
        self.assertTrue(self.conf.validate())

    def test_run(self):
        solver = self.conf.run()
        np.testing.assert_allclose(solver.tNow, 0.01, atol=3e-5)

        fuel = ct.Solution('h2o2.yaml')
        prod = ct.Solution('h2o2.yaml')
        ox = ct.Solution('h2o2.yaml')

        fuel.TPX = 300, 101325, 'H2:1.0, AR:1.0'
        ox.TPX = 300, 101325, 'O2:1.0, AR:4.0'
        prod.TPX = 300, 101325, 'H2:1.0, O2:0.5, AR:3.0'
        prod.equilibrate('HP')

        # fuel side (left)
        np.testing.assert_allclose(fuel.T, solver.T[0])
        np.testing.assert_allclose(fuel.Y, solver.Y[:,0], atol=1e-9)

        # oxidizer side (right)
        np.testing.assert_allclose(ox.T, solver.T[-1])
        np.testing.assert_allclose(ox.Y, solver.Y[:,-1], atol=1e-9)

        # products
        self.assertTrue(prod.T - 300 < max(solver.T) < prod.T)
