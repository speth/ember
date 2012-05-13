import os
import platform
from distutils.sysconfig import get_config_var

VariantDir('src/build','src', duplicate=0)
VariantDir('test/build','test', duplicate=0)
VariantDir('python/build','python', duplicate=0)

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
        'cantera_include',
        'Location of the Cantera header files.',
        '/usr/local/include', PathVariable.PathIsDir),
    PathVariable(
        'cantera_libs',
        'Location of the Cantera library files.',
        '/usr/local/lib', PathVariable.PathIsDir),
    PathVariable(
        'sundials_include',
        'Location of the Sundials header files.',
        '/usr/local/include', PathVariable.PathIsDir),
    PathVariable(
        'sundials_libs',
        'Location of the Sundials library files',
        '/usr/local/lib', PathVariable.PathIsDir),
    PathVariable(
        'eigen_include',
        'Location of the Eigen header files',
        '/usr/local/include', PathVariable.PathIsDir),
    PathVariable(
        'boost_include',
        'Location of the Boost header files.',
        '/usr/include', PathVariable.PathIsDir),
    PathVariable(
        'boost_libs',
        'Location of the Boost library files',
        '/usr/lib', PathVariable.PathIsDir),
    PathVariable(
        'hdf5_include',
        'Location of the HDF5 header files.',
        '/usr/include', PathVariable.PathIsDir),
    PathVariable(
        'hdf5_libs',
        'Location of the HDF5 library files',
        '/usr/lib', PathVariable.PathIsDir),
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

lastlibs = hdf5 + env['blas_lapack'].split(',')

env.Append(CPPPATH=[env['cantera_include'],
                    env['sundials_include'],
                    env['eigen_include'],
                    env['boost_include'],
                    env['hdf5_include'],
                    get_config_var('INCLUDEPY')],
           LIBPATH=[env['cantera_libs'],
                    env['sundials_libs'],
                    env['boost_libs'],
                    env['hdf5_libs'],
                    'lib'],
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

common = [f for f in Glob('src/build/*.cpp')
          if 'strainedFlame.cpp' not in f.name]

# The Python module
pyenv = env.Clone()

env.Alias('pylib',
          pyenv.SharedLibrary('lib/_pyro',
                              common + Glob('python/build/*.cpp'),
                              SHLIBPREFIX='',
                              SHLIBSUFFIX=get_config_var('SO')))

# GoogleTest tests
python_lib = 'python%s' % get_config_var('VERSION')
testenv = env.Clone()
testenv.Append(LIBS=['pylib', 'gtest', python_lib],
               CPPPATH=['ext/gtest/include'],
               LIBPATH=['lib'])

if os.name == 'nt':
    testenv.Append(LIBPATH=get_config_var('LIBDEST'))

test_program = testenv.Program('bin/unittest',
                               Glob('test/build/*.cpp'))
test_alias = testenv.Alias('test', [test_program], test_program[0].abspath)
# AlwaysBuild(test_alias)

# Google Test
buildTargets = ['pylib', test_alias]

Export('env', 'buildTargets')
VariantDir('ext/build', 'ext', duplicate=0)
SConscript('ext/build/SConscript')

Default(['pylib'])
env.Alias('all', buildTargets)
