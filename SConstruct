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

env = Environment(tools=['default', 'subst'], **extraEnvArgs)

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
    PathVariable(
        'python_cmd',
        """Path to the Python interpreter to be used for building the Python module,
           if different from the interpreter being used by SCons.""",
        sys.executable, PathVariable.PathAccept),
    BoolVariable(
        'debug_symbols',
        'Include debug symbols in the compiled module',
        True),
    BoolVariable(
        'debug',
        """Enable asserts, bounds checking and other debugging code; disable
           compiler optimizations.""",
        False),
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

lastlibs = ['tbb'] + hdf5 + env['blas_lapack'].split(',')

include_dirs = []
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

if env['CC'] == 'gcc':
    flags = ['-ftemplate-depth-128', '-fPIC', '-g', '-Wall']
    defines = []
    if env['debug_symbols']:
        flags.append('-g')

    if env['debug']:
        flags.extend(['-O0','-fno-inline'])
    else:
        flags.extend(['-O3', '-finline-functions', '-Wno-inline'])
        defines.append('NDEBUG')

    boost_libs = ['boost_%s' % lib
                  for lib in ('filesystem', 'system')]

elif env['CC'] == 'cl':
    flags = ['/nologo', '/W3', '/Zc:wchar_t', '/Zc:forScope', '/EHsc', '/MD']
    defines = ['_SCL_SECURE_NO_WARNINGS', '_CRT_SECURE_NO_WARNINGS',
                      'BOOST_ALL_DYN_LINK']
    if env['debug_symbols']:
        env.Append(LINKFLAGS='/DEBUG')
        flags.append('/Zi')

    if env['debug']:
        flags.extend(['/Od', '/Ob0', '/MD'])
    else:
        flags.append('/O2')
        defines.append('NDEBUG')

    library_dirs.append(get_config_var('prefix') + '/libs')
    boost_libs = []
    env['ENV']['MSSdk'] = 1
    env['ENV']['DISTUTILS_USE_SDK'] = 1

else:
    print 'error: unknown c++ compiler: "%s"' % env['CC']
    sys.exit(0)

env.Append(CPPPATH=include_dirs,
           LIBPATH=library_dirs,
           CXXFLAGS=flags,
           CPPDEFINES=defines,
           LIBS=sundials + cantera + lastlibs + boost_libs)

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

common_objects = env.SharedObject(Glob('build/core/*.cpp'))

corelib = env.Library('build/core/ember', common_objects)
env.Alias('build', corelib)

if os.name == 'nt':
    for dest in ['python/ember/TBB.dll', 'bin/TBB.dll']:
        tbb = env.Command(dest,
                          env['tbb']+'/bin/%s/%s/TBB.dll' % (tbbArch, tbbCompiler),
                          Copy('$TARGET', '$SOURCE'))
    env.Alias('build', tbb)
    py_gui1 = env.Command('python/scripts/ember-script.pyw', 'python/scripts/ember',
                          Copy('$TARGET', '$SOURCE'))
    script = ('''"from pkg_resources import resource_string; open('$TARGET', 'wb').write(resource_string('setuptools', 'gui.exe'))"''')
    py_gui2 = env.Command('python/scripts/ember.exe', py_gui1, '$python_cmd -c ' + script)
    env.Alias('build', [py_gui1, py_gui2])

# The Python module
env['py_include_dirs'] = repr(['../src'] + include_dirs)
env['py_libraries'] = repr([L for L in (['ember'] + env['LIBS']) if L])
env['py_libdirs'] = repr(['../build/core'] + env['LIBPATH'])
make_setup = env.SubstFile('#python/setup.py',
                           '#python/setup.py.in')
script = ('from distutils.sysconfig import *\n'
          'print(get_config_var("SO"))')
env['module_ext'] = getCommandOutput(env['python_cmd'], '-c', script).strip()

setup_cmd = 'cd python && $python_cmd setup.py '
build_args = '--build-lib=../build/python --build-temp=../build/tmp-python'
py_build_ext = env.Command('build/python/ember/_ember${module_ext}', make_setup,
                       setup_cmd + 'build_ext ' + build_args)
if os.name == 'nt':
    env.Depends(py_build_ext, [py_gui1, py_gui2])

env.AddPreAction(py_build_ext, Delete('build/python/ember/_ember${module_ext}'))
env.Depends(py_build_ext, corelib)
env.Depends(py_build_ext, mglob(env, 'python/ember', 'h', 'pxd', 'pyx'))

py_build = env.Command('build/python/ember/__init__.py', py_build_ext,
                       setup_cmd + 'build ' + build_args)
env.Alias('build', py_build)
env.Depends(py_build, mglob(env, 'python/ember', 'py'))

py_install = env.Command('dummy_target', py_build, setup_cmd + 'install $install_args')
env.Alias('install', py_install)

# GoogleTest tests
testenv = env.Clone()
testenv.Append(LIBS=['gtest'],
               CPPPATH=['ext/gtest/include'],
               LIBPATH=['lib'])

if os.name == 'nt':
    testenv.Append(LIBPATH=get_config_var('LIBDEST'))
    testenv['ENV']['PATH'] += ';' + Dir('lib').abspath

test_program = testenv.Program('bin/unittest',
                               Glob('build/test/*.cpp') + common_objects)
run_test = testenv.Command('test_dummy', test_program[0].abspath, '$SOURCE $TARGET')
test_alias = Alias('test', run_test)
# AlwaysBuild(test_alias)

# Google Test
buildTargets = ['build', test_alias]

Export('env', 'buildTargets')
VariantDir('ext/build', 'ext', duplicate=0)
SConscript('ext/build/SConscript')
