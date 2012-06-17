Getting Started with Pyro
=========================

1. **Dependencies:**
    These directions assume you are using a Linux or Windows system with the
    following installed:

    * A C++ compiler (g++ or Visual Studio 2008)
    * Python (2.6 or 2.7)
    * Boost, including boost-python and boost-filesystem
    * SCons (2.1.0 recommended)
    * Eigen (>= 3.0)
    * Cantera (>= 2.0.0)
    * Sundials (== 2.4)
    * HDF5 (>= 1.8.0)
    * numpy
    * h5py
    * Intel Thread Building Blocks (>= 4.0)

    Additional Dependencies for processing documentation:

    * Doxygen (>= 1.8.0)
    * Sphinx
    * Pygments
    * pyparsing
    * doxylink (see http://pypi.python.org/pypi/sphinxcontrib-doxylink)

2. The code is stored in a Git repository on Pharos. To check out a copy of the
   code, run:

        $ cd
        $ mkdir -p src
        $ cd src
        $ git clone pharos.mit.edu:/var/cache/git/1dflameV2.git

3. Edit the file ``1dflameV2/SConstruct`` to point to the correct include
   and library directories for your system.

4. Compiling the Code:

        $ cd ~/src/1dflameV2
        $ scons

    This produces the python extension ``python/pyro/_pyro.so``.

5. Add the "Pyro" Python module to your Python search path. The easiest way to
   do this is by creating a symlink into your user Python module path:

        $ mkdir -p ~/.local/lib/pythonX.Y/site-packages
        $ ln -s ~/src/1dflameV2/lib ~/.local/lib/pythonX.Y/site-packages/pyro

    where *X.Y* is your your installed version of Python, e.g 2.7.

    And similarly for Cantera, if it is not already on the Python path:

        $ ln -s /path/to/cantera/lib/pythonX.Y/site-packages/Cantera ~/.local/lib/pythonX.Y/site-packages/

6. Prepare the documentation (optional)

        $ doxygen
        $ cd doc/sphinx
        $ make

    To view the HTML docs, open ``doc/sphix/html/index.html`` in your web browser.

7. Prepare your kinetics mechanism. If you're starting with a mechanism that's
   in the Chemkin format, it will need to be converted:

        $ source /wherever/cantera/is/installed/setup_cantera
        $ ck2cti2 -input=<mech> --thermo=<thermo> --transport=<transport>

   This will produce a ``.cti`` file, which needs to be further converted:

        $ cti2ctml <mechname>.cti

   to produce ``<mechname>.xml``.

8. Prepare the input file, based on the one specified in

        1dflameV2/input/example-premixed.py

   A complete list of all available input parameters may be found in the HTML
   documentation. Alternatively, you look at the definitions in
   ``python/pyro/input.py``. Specify the path to your mechanism file (which
   needs to be in whatever directory is specified as the "input" directory in
   input file) as the ``mechanismFile`` and ``"gas"`` as the ``phaseID``. The
   other parameters you may want to change are in the ``Paths``,
   ``InitialCondition``, and ``StrainParameters`` sections.

9. Check the configuration for errors:

        $ python myInputFile.py validate

   If this prints "Validation completed successfully.", you're all set.
   Otherwise, try to correct the indicated error and try again.

10. Run the code:

        $ python myInputFile.py &

    This may take a while. You can watch the solver's progress as it is written
    to the file specified by ``Paths.logFile`` in the input file.

11. Examine the output files. The files are HDF5 data files, which can be read
    using the Python ``h5py`` module or Matlab.

    * ``out.h5`` contains integral flame properties (e.g. flame speed) as a
      function of time
    * ``profNNNNNN.h5`` contain the temperature & species profiles output
      periodically.
    * ``profNow.h5`` contains the most recently saved profiles.

    Using the ``pyro.utils`` module:

        $ ipython -pylab
        >>> import pyro
        >>> prof = pyro.load('run/test/profNow.h5')
        >>> plot(prof.x, prof.T)
