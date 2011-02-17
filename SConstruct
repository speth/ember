VariantDir('build','src', duplicate=0)

debugCppFlags = ['-O0','-ftemplate-depth-128','-fno-inline','-Wall','-g','-fPIC']
releaseCppFlags = ['-ftemplate-depth-128', '-O3', '-finline-functions', '-Wno-inline',
                   '-fPIC', '-Wall', '-g', '-DNDEBUG']

mode = ARGUMENTS.get('mode','release')
cppFlags = releaseCppFlags if mode == 'release' else debugCppFlags

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
                  LIBS=startlibs + sundials + cantera + lastlibs)


Library('lib/libcklib.a', 'adapchem/cklib.f')

Library('lib/libadapchem.a', 
        ['adapchem/adapchem.f', 
         'adapchem/box.f', 
         'adapchem/ckcompat.cpp', 
         'adapchem/wrappers.f90'])

env.Program('bin/1dflameV2', Glob('build/*.cpp'))
