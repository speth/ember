"""
SCons build script for Ember

Basic usage:

    'scons help' - print a description of user-specifiable options

    'scons build' - Compile Ember

    'scons test' - Run the test suite

    'scons install' - Install the Ember Python module

    'scons msi' - Create a MSI installer for Windows
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
class defaults: pass

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

    defaults.blas_lapack = ''
elif platform.system() == 'Darwin':
    defaults.blas_lapack = ''
else:
    defaults.blas_lapack = 'lapack,blas'

env = Environment(tools=['default', 'subst'], **extraEnvArgs)

opts = Variables('ember.conf')
opts.AddVariables(
    PathVariable(
        'cantera',
        'Location of the Cantera header and library files.',
        '', PathVariable.PathAccept),
    ('cantera_libs',
     """Comma-separated list of extra libraries required to link with
        Cantera. If Cantera was not built with the "single_library=y" option,
        this may include some of 'execstream', 'ctmath', 'ctlapack', 'ctblas',
        and 'ctf2c'""",
     ''),
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
        'Location of the Boost header files.',
        '', PathVariable.PathAccept),
    ('boost_libs',
     """Comma-separated list of extra libraries required to link using Boost
        Thread. If unspecified, SCons will attempt to automatically determine
        the correct libraries. This option should be set if that fails.
     """,
     ''),
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
     defaults.blas_lapack),
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

cantera = ['cantera'] + env['cantera_libs'].split(',')
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
if os.name == 'nt':
    hdf5 = ['hdf5','libzlib', 'libszip']
    tbbLibDir = env['tbb']+'/lib/%s/%s' % (tbbArch, tbbCompiler)
else:
    hdf5 = ['hdf5']
    tbbLibDir = env['tbb']+'/lib'

if platform.system() == 'Darwin':
    env.Append(FRAMEWORKS=['Accelerate'])

lastlibs = ['tbb'] + hdf5 + env['blas_lapack'].split(',')

include_dirs = []
library_dirs = [env.Dir('lib').abspath]

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
    flags = ['-ftemplate-depth-128', '-fPIC', '-g', '-Wall', '-pthread']
    linkflags = ['-pthread']
    defines = []
    if env['debug_symbols']:
        flags.append('-g')

    if env['debug']:
        flags.extend(['-O0','-fno-inline'])
    else:
        flags.extend(['-O3', '-finline-functions', '-Wno-inline'])
        defines.append('NDEBUG')

elif env['CC'] == 'cl':
    flags = ['/nologo', '/W3', '/Zc:wchar_t', '/Zc:forScope', '/EHsc', '/MD']
    linkflags=[]
    defines = ['_SCL_SECURE_NO_WARNINGS', '_CRT_SECURE_NO_WARNINGS']
    if env['debug_symbols']:
        env.Append(LINKFLAGS='/DEBUG')
        flags.append('/Zi')

    if env['debug']:
        flags.extend(['/Od', '/Ob0', '/MD'])
    else:
        flags.append('/O2')
        defines.append('NDEBUG')

    library_dirs.append(get_config_var('prefix') + '/libs')
    env['ENV']['MSSdk'] = 1
    env['ENV']['DISTUTILS_USE_SDK'] = 1

else:
    print 'error: unknown c++ compiler: "%s"' % env['CC']
    sys.exit(0)

env.Append(CPPPATH=include_dirs,
           LIBPATH=library_dirs,
           CXXFLAGS=flags,
           CPPDEFINES=defines,
           LINKFLAGS=linkflags,
           LIBS=sundials + cantera + lastlibs)

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

def get_expression_value(includes, expression):
    s = ['#include ' + i for i in includes]
    s.extend(('#include <iostream>',
              'int main(int argc, char** argv) {',
              '    std::cout << %s << std::endl;' % expression,
              '    return 0;',
              '}\n'))
    return '\n'.join(s)

configInfo = {}

import SCons.Conftest, SCons.SConf
tests = {}
conf = Configure(env, custom_tests={'CheckMemberFunction': CheckMemberFunction})
context = SCons.SConf.CheckContext(conf)

# Check for required headers
fail = False
for header, quotes in [('cantera/thermo/IdealGasPhase.h', '""'),
                       ('cvode/cvode.h', '<>'),
                       ('Eigen/Dense', '<>'),
                       ('tbb/parallel_for.h', '""'),
                       ('hdf5.h', '<>')]:
    fail |= SCons.Conftest.CheckHeader(context, header, language='C++',
                                     include_quotes=quotes)
if fail:
    raise EnvironmentError("Failed to a required header file. "
                           "See config.log for details.")

# Check for required libraries
src = get_expression_value(["<cmath>"], "sin(3.14)")
retcode, retval = conf.TryRun(src, '.cpp')
if not retcode:
    raise EnvironmentError("Failed to a required library."
                           "See config.log for details.")

# Header optionally used by gtest
env['HAS_TR1_TUPLE'] = conf.CheckCXXHeader('tr1/tuple', '<>')

# Figure out what needs to be linked for Boost Thread support (used by Cantera)
boost_ok = False
if env['boost_libs']:
    boost_lib_choices = [env['boost_libs'].split(',')]
else:
    boost_lib_choices = [[''], ['boost_system'], ['boost_thread', 'boost_system']]
for bt in boost_lib_choices:
    header= "#define BOOST_ALL_NO_LIB\n#include <boost/thread/thread.hpp>"
    call = "boost::mutex foo; boost::mutex::scoped_lock bar(foo);"

    ans = SCons.Conftest.CheckLib(context,
                                  bt,
                                  header=header,
                                  language='C++',
                                  call=call,
                                  autoadd=False)
    if not ans:
        boost_ok = True
        print 'Linking Boost thread with the following libraries: {0}'.format(
            ', '.join(bt) or '<none>')
        if bt[0] and os.name != 'nt':
            env.Append(LIBS=bt)
        break

if not boost_ok:
    raise EnvironmentError("Couldn't determine correct libraries to link for Boost Thread.")

# Determine Sundials version
sundials_version_source = get_expression_value(['"sundials/sundials_config.h"'],
                                                   'SUNDIALS_PACKAGE_VERSION')
retcode, sundials_version = conf.TryRun(sundials_version_source, '.cpp')
if retcode == 0:
    print "Failed to determine Sundials version."
    print "See 'config.log' for details."
    sys.exit(1)

# Ignore the minor version and convert to integer, e.g. 2.4.x -> 24
configInfo['SUNDIALS_VERSION'] = ''.join(sundials_version.strip().split('.')[:2])
print """INFO: Using Sundials version %s""" % sundials_version.strip()

tests['CanteraExtendedTransport'] = conf.CheckMemberFunction(
    "Cantera::MixTransport::getMixDiffCoeffsMass",
    includes='#include "cantera/transport.h"')
tests['CanteraThreadSafe'] = conf.CheckDeclaration('THREAD_SAFE_CANTERA',
                                                   '#include "cantera/base/config.h"',
                                                   'C++')

src = """
#include "cantera/transport/MultiTransport.h"
class Foo : public Cantera::MultiTransport {
    virtual void solveLMatrixEquation() {
        m_Lmatrix; // must not be declared private
    }
};
"""
tests['CanteraExtendedMultiTransport'] = conf.TryCompile(src, '.cpp')

if not tests['CanteraThreadSafe']:
    raise EnvironmentError("Ember requires a thread-safe version of Cantera. "
                           "See 'config.log' for details.")

if not tests['CanteraExtendedTransport']:
    raise EnvironmentError("Missing required Cantera method 'getMixDiffCoeffsMass'. "
                           "See 'config.log' for details.")

if tests['CanteraExtendedMultiTransport']:
    configInfo['EMBER_EXTENDED_MULTITRANSPORT'] = 1
else:
    print ('WARNING: Extended MultiTransport class unavailable. Multicomponent'
           'transport model may not work correctly with multiple threads')

config_h = env.Command('src/config.h',
                       'src/config.h.in',
                       ConfigBuilder(configInfo))
env.AlwaysBuild(config_h)

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

py_msi = env.Command('dummy_target2', py_build,
                     setup_cmd + 'bdist_msi --dist-dir=..')
env.Alias('msi', py_msi)

# GoogleTest tests
testenv = env.Clone()
testenv.Append(LIBS=['gtest'],
               CPPPATH=['ext/gtest/include'])

if os.name == 'nt':
    testenv.Append(LIBPATH=get_config_var('LIBDEST'))
    testenv['ENV']['PATH'] += ';' + Dir('lib').abspath

if not env['HAS_TR1_TUPLE']:
    testenv.Append(CPPDEFINES='GTEST_USE_OWN_TR1_TUPLE')

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
