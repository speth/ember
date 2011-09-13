VariantDir('src/build','src', duplicate=0)
VariantDir('test/build','test', duplicate=0)
VariantDir('python/build','python', duplicate=0)
VariantDir('adapchem/build','adapchem', duplicate=0)

try:
    import multiprocessing
    SetOption('num_jobs', multiprocessing.cpu_count())
except:
    pass

debugCppFlags = ['-O0','-ftemplate-depth-128','-fno-inline','-Wall','-g','-fPIC']
releaseCppFlags = ['-ftemplate-depth-128', '-O3', '-finline-functions', '-Wno-inline',
                   '-fPIC', '-Wall', '-g', '-DNDEBUG']

debugFortFlags = ['-O0', '-g','-fPIC']
releaseFortFlags = ['-O3','-fPIC','-g']

mode = ARGUMENTS.get('mode','release')
cppFlags = releaseCppFlags if mode == 'release' else debugCppFlags
fortFlags = releaseFortFlags if mode == 'release' else debugFortFlags

cantera = '''thermo transport kinetics equil tpx ctnumerics
             ctmath ctf2c ctcxx ctbase clib'''.split()
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
lastlibs = 'gfortran hdf5 blas lapack boost_filesystem'.split()
pythonlibs = 'boost_python python2.6'.split()

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

env = Environment(CPPPATH=['/opt/cantera-gcc/include',
                           '/opt/sundials-2.4.0-gcc/include',
                           '/usr/include/python2.6'],
                  LIBPATH=['/opt/cantera-gcc/lib',
                           '/opt/sundials-2.4.0-gcc/lib',
                           'lib'],
                  CPPFLAGS=cppFlags,
                  FORTRANFLAGS=fortFlags,
                  F90FLAGS=fortFlags,
                  LIBS=sundials + cantera + pythonlibs + lastlibs)

tests = {}
conf = Configure(env, custom_tests={'CheckMemberFunction': CheckMemberFunction})
tests['CanteraExtendedTransport'] = conf.CheckMemberFunction(
    "Cantera::MixTransport::getMixDiffCoeffsMass",
    includes="#include <cantera/Cantera.h>\n#include <cantera/transport.h>")

if tests['CanteraExtendedTransport']:
    env.Append(CPPFLAGS=['-DCANTERA_EXTENDED_TRANSPORT'])

cklib = env.Library('lib/libcklib.a', 'adapchem/build/cklib.f')

adapchem = env.Library('lib/libadapchem',
            ['adapchem/build/adapchem.f',
             'adapchem/build/box.f',
             'adapchem/build/ckcompat.cpp',
             'adapchem/build/wrappers.f90'])

env['LIBS'] = cklib + adapchem + env['LIBS']

common = [f for f in Glob('src/build/*.cpp')
          if 'strainedFlame.cpp' not in f.name]

# The Python module
pyenv = env.Clone()
pyenv.Append(LIBS=pythonlibs)
env.Alias('pylib',
          pyenv.SharedLibrary('lib/_pyro.so',
                              common + Glob('python/build/*.cpp'),
                              SHLIBPREFIX=''))

# UnitTest++ tests
testenv = env.Clone()
testenv.Append(LIBS=['UnitTest++'],
               CPPPATH=['../UnitTest++/src/'],
               LIBPATH=['../UnitTest++'])
test_program = testenv.Program('bin/unittest',
                               common + Glob('test/build/*.cpp'))
test_alias = testenv.Alias('test', [test_program], test_program[0].abspath)
AlwaysBuild(test_alias)

Default(['pylib'])
env.Alias('all', ['pylib','test',adapchem,cklib])
