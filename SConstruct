
VariantDir('src/build','src', duplicate=0)
VariantDir('test/build','test', duplicate=0)

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

startlibs = 'adapchem cklib'.split()
cantera = '''thermo transport kinetics equil tpx ctnumerics 
             ctmath ctf2c ctcxx ctbase clib'''.split()
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
lastlibs = 'gfortran hdf5 config++ blas lapack boost_filesystem'.split()

env = Environment(CPPPATH=['/opt/cantera-gcc/include',
                           '/opt/sundials-2.4.0-gcc/include'],
                  LIBPATH=['/opt/cantera-gcc/lib',
                           '/opt/sundials-2.4.0-gcc/lib',
                           'lib'],
                  CPPFLAGS=cppFlags,
                  FORTRANFLAGS=fortFlags,
                  F90FLAGS=fortFlags,
                  LIBS=startlibs + sundials + cantera + lastlibs)

env.Library('lib/libcklib.a', 'adapchem/cklib.f')

env.Library('lib/libadapchem.a',
            ['adapchem/adapchem.f',
             'adapchem/box.f',
             'adapchem/ckcompat.cpp',
             'adapchem/wrappers.f90'])

common = [f for f in Glob('src/build/*.cpp')
          if 'strainedFlame.cpp' not in f.name]

# The main 1dflame program
pyro = env.Program('bin/1dflameV2', common + ['src/build/strainedFlame.cpp'])
Default(pyro)
env.Alias('pyro', pyro)

env.Alias('qsstest', env.Program('bin/qsstest', common+['test/build/test_qssintegrator.cpp']))
env.Alias('all', ['pyro','qsstest'])
