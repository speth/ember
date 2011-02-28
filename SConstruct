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

startlibs = 'adapchem cklib'.split()
cantera = '''thermo transport kinetics equil tpx ctnumerics 
             ctmath ctf2c ctcxx ctbase clib'''.split()
sundials = 'sundials_nvecserial sundials_ida sundials_cvode'.split()
lastlibs = 'gfortran hdf5 config++ blas lapack boost_filesystem'.split()
pythonlibs = 'boost_python python2.6'.split()

env = Environment(CPPPATH=['/opt/cantera-gcc/include',
                           '/opt/sundials-2.4.0-gcc/include',
                           '/usr/include/python2.6'],
                  LIBPATH=['/opt/cantera-gcc/lib',
                           '/opt/sundials-2.4.0-gcc/lib',
                           'lib'],
                  CPPFLAGS=cppFlags,
                  FORTRANFLAGS=fortFlags,
                  F90FLAGS=fortFlags,
                  LIBS=startlibs + sundials + cantera + pythonlibs + lastlibs)

env.Library('lib/libcklib.a', 'adapchem/build/cklib.f')

env.Library('lib/libadapchem.a',
            ['adapchem/build/adapchem.f',
             'adapchem/build/box.f',
             'adapchem/build/ckcompat.cpp',
             'adapchem/build/wrappers.f90'])

common = [f for f in Glob('src/build/*.cpp')
          if 'strainedFlame.cpp' not in f.name]

# The Python module
pyenv = env.Clone()
pyenv.Append(LIBS=pythonlibs)
env.Alias('pylib',
          pyenv.SharedLibrary('lib/_pyro.so',
                              common + Glob('python/build/*.cpp'),
                              SHLIBPREFIX=''))

# Test programs
env.Alias('qsstest',
          env.Program('bin/qsstest',
                      common+['test/build/test_qssintegrator.cpp']))

Default(['pylib'])
env.Alias('all', ['qsstest','pylib'])
