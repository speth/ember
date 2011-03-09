Getting Started with the 1dflameV2 Code

0. Dependencies: These directions assume you are using a Linux system
with the following installed:
   * Python (>=2.6)
   * Boost
   * SCons
   * Cantera
   * Sundials (2.4)
   * HDF5

1. The 1dflameV2 code is stored in a Git repository on Pharos. To
check out a copy of the code, run:

    $ cd
    $ mkdir -p src
    $ cd src
    $ git clone pharos.mit.edu:/var/cache/git/1dflameV2.git

2. Edit the file "1dflameV2/SConstruct" to point to the correct include
and library directories for your system.

3. Compiling the Code:

    $ cd 1dflameV2
    $ scons

This produces the python module "lib/_pyro.so".

4. Add the "Pyro" Python module to your Python search path. The easiest way
to do this is by creating a symlink into your user Python module path:

    $ mkdir -p ~/.local/lib/python2.6/site-packages
    $ ln -s ~/src/1dflameV2/lib ~/.local/lib/python2.6/site-packages/pyro

And similarly for Cantera (optional):

    $ ln -s /path/to/cantera/lib/python2.6/site-packages/Cantera ~/.local/lib/python2.6/site-packages/

5. Prepare your kinetics mechanism. If you're starting with a
mechanism that's in the Chemkin format, it will need to be converted:

    $ source /wherever/cantera/is/installed/setup_cantera
    $ ck2cti -i <mech> -t <thermo> -tr <transport>

This will produce mech.cti, which needs to be further converted:

    $ cti2ctml mech.cti

to produce mech.xml.

6. Prepare the input file, based on the one specified in

    1dflameV2/input/example-premixed.py

A complete list of all available input parameters may be found in
"lib/input.py". Specify the path to your mechanism file (which needs to
be in whatever directory is specified as the "input" directory in input file)
as the mechanismFile and "gas" as the phaseID. The other parameters you
may want to change are in the Paths, InitialCondition, and StrainParameters sections.

7a. Check the configuration for errors

    python myInputFile.py validate

If this prints "Validation completed successfully.", you're all set. Otherwise,
try to correct the indicated error and try again.

7b. Run the code:

     python myInputFile.py &

This may take a while. You can watch the solver's progress as it is written to
the file specified by "Paths.logFile" in the input file.

8. Examine the output files. The files are HDF5 data files, which can be read
using the Python h5py module or Matlab.

   * outNow.h5 contains integral flame properties (e.g. flame speed) as a function of time
   * profNNNNNN.h5 contain the temperature & species profiles output periodically.
   * profNow.h5 contains the most recently saved profiles.

Using the pyro utils module:

    $ ipython -pylab
    >>> import pyro
    >>> prof = pyro.load('run/test/profNow.h5')
    >>> plot(prof.x, prof.T)
