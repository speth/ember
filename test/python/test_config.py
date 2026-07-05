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


def test_deprecated_grid_tolerances_warn(capsys):
    conf = Config(Grid(vtol=0.1, dvtol=0.15))
    conf._warnDeprecated()
    out = capsys.readouterr().out
    assert 'deprecated' in out
    assert 'errTol' in out


def test_errtol_produces_no_warning(capsys):
    conf = Config(Grid(errTol=1e-3))
    conf._warnDeprecated()
    assert capsys.readouterr().out == ''


def test_absurd_errtol_warns(capsys):
    conf = Config(Grid(errTol=0.9))
    conf._warnDeprecated()
    assert 'errTol' in capsys.readouterr().out
