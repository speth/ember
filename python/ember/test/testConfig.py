import unittest

from ember import *

class TestConfig(unittest.TestCase):
    def test_paths(self):
        c = Config(Paths(outputDir='dirname',
                         logFile='filename'))
        self.assertEqual(c.paths.outputDir, 'dirname')
        self.assertEqual(c.paths.logFile, 'filename')

    def test_bad_option_name(self):
        with self.assertRaises(KeyError):
            c = Config(Paths(doesNotExist='invalid'))
