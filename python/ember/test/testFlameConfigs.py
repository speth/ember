import unittest

from ember import *
import cantera as ct
from cantera.test.utilities import CanteraTest

class TestPremixedStrained(CanteraTest):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/premixedStrained',
                  logFile='build/test/work/premixedStrained.txt'),
            General(nThreads=2),
            Chemistry(mechanismFile='h2o2.cti'),
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
        self.assertNear(solver.tNow, 0.01, atol=3e-5)

        # Reactants side
        gas = ct.Solution('h2o2.cti')
        gas.TPX = 300, 101325, 'H2:0.6, O2:1.0, AR:4.0'
        self.assertNear(gas.T, solver.T[0])
        for i in range(gas.n_species):
            self.assertNear(gas.Y[i], solver.Y[i,0])

        # Products side
        gas.equilibrate('HP')
        for i in range(gas.n_species):
            self.assertNear(gas.Y[i], solver.Y[i,-1], atol=1e-8)
        self.assertNear(gas.T, solver.T[-1], 5)

        # Since this is Le < 1
        self.assertGreater(max(solver.T), gas.T)


class TestTwinPremixedStrained(CanteraTest):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/twinPremixedStrained',
                  logFile='build/test/work/twinPremixedStrained.txt'),
            Chemistry(mechanismFile='h2o2.cti'),
            InitialCondition(fuel='H2:1.0',
                             oxidizer='O2:1.0, AR:4.0',
                             xRight=0.003,
                             equivalenceRatio=0.3),
            General(nThreads=2, twinFlame=True, unburnedLeft=False),
            StrainParameters(initial=800, final=800),
            Grid(vtol=0.2, dvtol=0.3),
            Times(regridStepInterval=10),
            TerminationCondition(tEnd=0.01, measurement=None))

    def test_validates(self):
        self.assertTrue(self.conf.validate())

    def test_run(self):
        solver = self.conf.run()
        self.assertNear(solver.tNow, 0.01, atol=3e-5)

        products = ct.Solution('h2o2.cti')
        equil_prod = ct.Solution('h2o2.cti')
        reactants = ct.Solution('h2o2.cti')

        # Reactants side
        reactants.TPX = 300, 101325, 'H2:0.6, O2:1.0, AR:4.0'
        self.assertNear(reactants.T, solver.T[-1])
        for i in range(reactants.n_species):
            self.assertNear(reactants.Y[i], solver.Y[i,-1])

        # Products side. Not necessarily equilibrium
        equil_prod.TPX = reactants.TPX
        equil_prod.equilibrate('HP')
        products.TPY = solver.T[0], 101325, solver.Y[:,0]
        self.assertTrue(equil_prod.T - 300 < products.T < equil_prod.T + 300)
        self.assertLess(products['H2'].X, 0.1*reactants['H2'].X)
        self.assertGreater(products['H2O'].X, 0.5*equil_prod['H2O'].X)


class TestDiffusion(CanteraTest):
    def setUp(self):
        self.conf = Config(
            Paths(outputDir='build/test/work/h2Diffusion',
                  logFile='build/test/work/h2Diffusion.txt'),
            General(nThreads=2),
            Chemistry(mechanismFile='h2o2.cti'),
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
        self.assertNear(solver.tNow, 0.01, atol=3e-5)

        fuel = ct.Solution('h2o2.cti')
        prod = ct.Solution('h2o2.cti')
        ox = ct.Solution('h2o2.cti')

        fuel.TPX = 300, 101325, 'H2:1.0, AR:1.0'
        ox.TPX = 300, 101325, 'O2:1.0, AR:4.0'
        prod.TPX = 300, 101325, 'H2:1.0, O2:0.5, AR:3.0'
        prod.equilibrate('HP')

        # fuel side (left)
        self.assertNear(fuel.T, solver.T[0])
        for i in range(fuel.n_species):
            self.assertNear(fuel.Y[i], solver.Y[i,0])

        # oxidizer side (right)
        self.assertNear(ox.T, solver.T[-1])
        for i in range(ox.n_species):
            self.assertNear(ox.Y[i], solver.Y[i,-1])

        # products
        print prod.T, solver.T
        self.assertTrue(prod.T - 300 < max(solver.T) < prod.T)
