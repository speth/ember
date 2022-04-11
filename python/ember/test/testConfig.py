import unittest

from ember import *


class TestOptions(unittest.TestCase):
    def test_set_options(self):
        opts = Paths(outputDir='dirname',
                     logFile='filename')
        self.assertEqual(opts.outputDir, 'dirname')
        self.assertEqual(opts.logFile, 'filename')

    def test_bad_option_group(self):
        with self.assertRaises(TypeError):
            c = Config('this is not an option group')

    def test_bad_option_name(self):
        with self.assertRaises(KeyError):
            c = Paths(doesNotExist='invalid')

    def test_bad_option_value(self):
        with self.assertRaises(ValueError):
            c = InitialCondition(equivalenceRatio='foo')

        with self.assertRaises(ValueError):
            c = InitialCondition(equivalenceRatio=-2)


class TestConfig(unittest.TestCase):
    def test_set_options(self):
        c = Config(Paths(outputDir='dirname'),
                   StrainParameters(initial=1234))
        self.assertEqual(c.paths.outputDir, 'dirname')
        self.assertEqual(c.strainParameters.initial, 1234)

        concrete = c.evaluate()
        self.assertEqual(c.paths.outputDir, 'dirname')
        self.assertEqual(c.strainParameters.initial, 1234)

    def test_validate_baseline(self):
        c = Config()
        self.assertTrue(c.validate())

    def test_validate_badmech(self):
        c = Config(Chemistry(mechanismFile='wxyz.yaml'))
        with self.assertRaises(Exception):
            c.validate()

    def test_validate_bad_species(self):
        c = Config(InitialCondition(fuel='CH5:1.0'))
        with self.assertRaises(Exception):
            c.validate()

    def test_validate_bad_restart(self):
        c = Config(InitialCondition(restartFile='nosuchfile.h5'))
        self.assertFalse(c.validate())
