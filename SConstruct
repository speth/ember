import os
import platform
from distutils.sysconfig import get_config_var

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
    elif pycomp.startswith('MSC v.1500'):
        extraEnvArgs['MSVC_VERSION'] = '9.0' # Visual Studio 2008
    elif pycomp.startswith('MSC v.1600'):
        extraEnvArgs['MSVC_VERSION'] = '10.0' # Visual Studio 2010

    if '64 bit' in pycomp:
        extraEnvArgs['TARGET_ARCH'] = 'amd64'
    else:
        extraEnvArgs['TARGET_ARCH'] = 'x86'

env = Environment(**extraEnvArgs)

opts = Variables('pyro.conf')
opts.AddVariables(
    PathVariable(
        'cantera',
        'Location of the Cantera header and library files.',
        '/usr/local', PathVariable.PathIsDir),
    PathVariable(
        'sundials',
        'Location of the Sundials header and library files.',
        '/usr/local', PathVariable.PathIsDir),
    PathVariable(
        'eigen',
        'Location of the Eigen header files',
        '/usr/local/include', PathVariable.PathIsDir),
    PathVariable(
        'boost',
        'Location of the Boost header and library files.',
        '/usr', PathVariable.PathIsDir),
    PathVariable(
        'hdf5',
        'Location of the HDF5 header and library files.',
        '/usr', PathVariable.PathIsDir),
    ('blas_lapack',
     'Comma-separated list of libraries to include for BLAS/LAPACK support',
     'blas,lapack')
    )

opts.Update(env)
opts.Save('pyro.conf', env)

cantera = ['cantera']
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
if os.name == 'nt':
    hdf5 = ['hdf5','libzlib', 'libszip']
else:
    hdf5 = ['hdf5']

python = ['python%s' % get_config_var('VERSION')]
lastlibs = ['tbb'] + python + hdf5 + env['blas_lapack'].split(',')

env.Append(CPPPATH=[env['cantera']+'/include',
                    env['sundials']+'/include',
                    env['eigen'],
                    env['boost']+'/include',
                    env['hdf5']+'/include',
                    get_config_var('INCLUDEPY')],
           LIBPATH=[env['cantera']+'/lib',
                    env['sundials']+'/lib',
                    env['boost']+'/lib',
                    env['hdf5']+'/lib'],
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

if tests['CanteraExtendedTransport']:
    env.Append(CPPDEFINES=['CANTERA_EXTENDED_TRANSPORT'])

common = [f for f in Glob('build/core/*.cpp')
          if 'strainedFlame.cpp' not in f.name]

# The Python module
pyenv = env.Clone()

pylib = pyenv.SharedLibrary('python/pyro/_pyro',
                            common + Glob('build/python/*.cpp'),
                            SHLIBPREFIX='',
                            SHLIBSUFFIX=get_config_var('SO'))

env.Alias('pylib', pylib)


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
buildTargets = ['pylib', test_alias]

Export('env', 'buildTargets')
VariantDir('ext/build', 'ext', duplicate=0)
SConscript('ext/build/SConscript')

Default(['pylib'])
env.Alias('all', buildTargets)
