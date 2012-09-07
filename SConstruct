"""
SCons build script for Ember

Basic usage:

    'scons help' - print a description of user-specifiable options

    'scons build' - Compile Ember

    'scons test' - Run the test suite

    'scons install' - Install the Ember Python module
"""

import os
import platform
from distutils.sysconfig import get_config_var
import numpy as np
from buildutils import *

if not COMMAND_LINE_TARGETS:
    print __doc__
    sys.exit(0)

VariantDir('build/core','src', duplicate=0)
VariantDir('build/test','test', duplicate=0)
VariantDir('build/python','src/python', duplicate=0)

mode = ARGUMENTS.get('mode','release')
assert mode in ('release', 'debug'), mode

#try:
#    import multiprocessing
#    SetOption('num_jobs', multiprocessing.cpu_count())
#except:
#    pass

extraEnvArgs = {}

if os.name == 'nt':
    # On Windows, use the same version of Visual Studio that was used
    # to compile Python, and target the same architecture.
    pycomp = platform.python_compiler()
    if pycomp.startswith('MSC v.1400'):
        extraEnvArgs['MSVC_VERSION'] = '8.0' # Visual Studio 2005
        tbbCompiler = 'vc8'
    elif pycomp.startswith('MSC v.1500'):
        extraEnvArgs['MSVC_VERSION'] = '9.0' # Visual Studio 2008
        tbbCompiler = 'vc9'
    elif pycomp.startswith('MSC v.1600'):
        extraEnvArgs['MSVC_VERSION'] = '10.0' # Visual Studio 2010
        tbbCompiler = 'vc10'

    if '64 bit' in pycomp:
        extraEnvArgs['TARGET_ARCH'] = 'amd64'
        tbbArch = 'intel64'
    else:
        extraEnvArgs['TARGET_ARCH'] = 'x86'
        tbbArch = 'ia32'

env = Environment(**extraEnvArgs)

opts = Variables('ember.conf')
opts.AddVariables(
    PathVariable(
        'cantera',
        'Location of the Cantera header and library files.',
        '', PathVariable.PathAccept),
    PathVariable(
        'sundials',
        'Location of the Sundials header and library files.',
        '', PathVariable.PathAccept),
    PathVariable(
        'eigen',
        'Location of the Eigen header files',
        '', PathVariable.PathAccept),
    PathVariable(
        'boost',
        'Location of the Boost header and library files.',
        '', PathVariable.PathAccept),
    PathVariable(
        'hdf5',
        'Location of the HDF5 header and library files.',
        '', PathVariable.PathAccept),
    PathVariable(
        'tbb',
        'Location of the Thread Building Blocks (TBB) header and library files',
        '', PathVariable.PathAccept),
    ('blas_lapack',
     'Comma-separated list of libraries to include for BLAS/LAPACK support',
     'blas,lapack'),
    ('install_args',
     'Command-line arguments passed to "python setup.py install"',
     '--user')
    )

opts.Update(env)
opts.Save('ember.conf', env)

if 'help' in COMMAND_LINE_TARGETS:
    # Print help about configuration options and exit
    print """
        ************************************************
        *   Configuration options for building Ember   *
        ************************************************

The following options can be passed to SCons to customize the Ember build
process. They should be given in the form:

    scons build option1=value1 option2=value2

Variables set in this way will then be stored in the 'ember.conf' file and
reusd automatically on subsequent invocations of scons. Alternatively, the
configuration variables can be entered directly into 'ember.conf' before
running 'scons build'. The format of this file is:

    option1 = 'value1'
    option2 = 'value2'

        ************************************************
"""
    for opt in opts.options:
        print '\n'.join(formatOption(env, opt))
    sys.exit(0)

cantera = ['cantera']
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
if os.name == 'nt':
    hdf5 = ['hdf5','libzlib', 'libszip']
    tbbLibDir = env['tbb']+'/lib/%s/%s' % (tbbArch, tbbCompiler)
else:
    hdf5 = ['hdf5']
    tbbLibDir = env['tbb']+'/lib'

python = ['python%s' % get_config_var('VERSION')]
lastlibs = ['tbb'] + python + hdf5 + env['blas_lapack'].split(',')

include_dirs = ['ext/boost.numpy/include']
library_dirs = []

if env['cantera']:
    include_dirs.append(env['cantera'] + '/include')
    library_dirs.append(env['cantera'] + '/lib')

if env['sundials']:
    include_dirs.append(env['sundials'] + '/include')
    library_dirs.append(env['sundials'] + '/lib')

if env['eigen']:
    include_dirs.append(env['eigen'])

if env['boost']:
    include_dirs.append(env['boost'] + '/include')
    library_dirs.append(env['boost'] + '/lib')

if env['hdf5']:
    include_dirs.append(env['hdf5'] + '/include')
    library_dirs.append(env['hdf5'] + '/lib')

if env['tbb']:
    include_dirs.append(env['tbb'] + '/include')
    library_dirs.append(tbbLibDir)

include_dirs.extend([get_config_var('INCLUDEPY'),
                     np.get_include()])

env.Append(CPPPATH=include_dirs,
           LIBPATH=library_dirs,                  
           LIBS=sundials + cantera + lastlibs)

if env['CC'] == 'gcc':
    common = ['-ftemplate-depth-128', '-fPIC', '-g', '-Wall']
    if mode == 'debug':
        env.Append(CXXFLAGS=common + ['-O0','-fno-inline'])
    else:
        env.Append(CXXFLAGS=common + ['-O3', '-finline-functions', '-Wno-inline'],
                   CPPDEFINES=['NDEBUG'])
    boost_libs = ['boost_%s' % lib
                  for lib in ('python', 'filesystem', 'system')]

elif env['CC'] == 'cl':
    common_flags = ['/nologo', '/Zi', '/W3', '/Zc:wchar_t', '/Zc:forScope', '/EHsc']
    common_defines = ['_SCL_SECURE_NO_WARNINGS', '_CRT_SECURE_NO_WARNINGS',
                      'BOOST_ALL_DYN_LINK']
    if mode == 'debug':
        env.Append(CXXFLAGS=common_flags + ['/Od', '/Ob0', '/MD'],
                   CPPDEFINES=common_defines,
                   LINKFLAGS='/DEBUG')
    else:
        env.Append(CXXFLAGS=common_flags + ['/O2', '/MD'],
                   CPPDEFINES=['NDEBUG'] + common_defines)
    boost_libs = []
    env.Append(LIBPATH=get_config_var('prefix') + '/libs')

else:
    print 'error: unknown c++ compiler: "%s"' % env['CC']
    sys.exit(0)

env.Append(LIBS=boost_libs)

def CheckMemberFunction(context, function, includes=""):
    context.Message('Checking for %s... ' % function)
    src = """
%(include)s
int main(int argc, char** argv) {
    &%(func)s;
    return 0;
}
""" % {'func':function,
       'include':includes}
    result = context.TryLink(src, '.cpp')
    context.Result(result)
    return result

tests = {}
conf = Configure(env, custom_tests={'CheckMemberFunction': CheckMemberFunction})
tests['CanteraExtendedTransport'] = conf.CheckMemberFunction(
    "Cantera::MixTransport::getMixDiffCoeffsMass",
    includes='#include "cantera/transport.h"')

if not tests['CanteraExtendedTransport']:
    raise EnvironmentError("Missing required Cantera method 'getMixDiffCoeffsMass'.")

common = Glob('build/core/*.cpp')

# The Python module
pyenv = env.Clone()
pyenv.Append(LIBS=['boost-numpy'],
             LIBPATH=['lib'])

pylib = pyenv.SharedLibrary('python/ember/_ember',
                            common + Glob('build/python/*.cpp'),
                            SHLIBPREFIX='',
                            SHLIBSUFFIX=get_config_var('SO'))

env.Alias('build', pylib)

if os.name == 'nt':
    tbb = env.Command('python/ember/TBB.dll',
                      env['tbb']+'/bin/%s/%s/TBB.dll' % (tbbArch, tbbCompiler),
                      Copy('$TARGET', '$SOURCE'))
    env.Alias('build', tbb)

# Installation of the Python module:
py_install = env.Command('dummy_target', pylib,
                         ('cd python && '
                          'python setup.py install %s' % env['install_args']))
env.Alias('install', py_install)

# GoogleTest tests
python_lib = 'python%s' % get_config_var('VERSION')
testenv = env.Clone()
testenv.Append(LIBS=[pylib, 'gtest', python_lib],
               CPPPATH=['ext/gtest/include'],
               LIBPATH=['lib'])

if os.name == 'nt':
    testenv.Append(LIBPATH=get_config_var('LIBDEST'))

test_program = testenv.Program('bin/unittest',
                               Glob('build/test/*.cpp'))
test_alias = testenv.Alias('test', [test_program], test_program[0].abspath)
# AlwaysBuild(test_alias)

# Google Test
buildTargets = ['build', test_alias]

Export('env', 'buildTargets')
VariantDir('ext/build', 'ext', duplicate=0)
SConscript('ext/build/SConscript')
