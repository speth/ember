from distutils.core import setup, Extension
from distutils.sysconfig import get_config_var
import os

dataFiles = ['_ember%s' % get_config_var('SO')]

if os.name == 'nt':
    for f in os.listdir('ember'):
        # Expect to find: tbb, boost_filesystem, boost_python, boost_system
        if f.lower().endswith('.dll'):
            dataFiles.append(f)

setup(name='Ember',
      version='1.0',
      description='Multipurpose 1D Flame Solver',
      long_description='',
      author='Raymond L. Speth',
      author_email='speth@mit.edu',
      url='http://github.com/speth/ember',
      packages = ['ember'],
      package_data = {'ember': dataFiles})

# From the directory containing this script:
# Install locally:
#     $ python setup.py install
# Install for your user:
#     $ python setup.py install --user
# Create a Windows .msi installer:
#     > python setup.py bdist_msi --target-version=2.7 --dist-dir=.
