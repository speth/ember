Getting Started with Ember
==========================

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
    * Intel Thread Building Blocks (>= 4.0)
    * numpy
    * h5py
    * git

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
        $ git clone git://github.com/speth/ember.git

3. Edit the file ``ember/SConstruct`` to point to the correct include
   and library directories for your system.

4. Compiling the Code:

        $ cd ~/src/ember
        $ scons

    This produces the python extension ``python/ember/_ember.so``.

5. Install the "Ember" Python module:

        $ cd python
        $ python setup.py install --user

   Alternatively, if you're doing active development, you can just create a
   symlink into your user Python module path:

        $ mkdir -p ~/.local/lib/pythonX.Y/site-packages
        $ ln -s ~/src/ember/python/ember ~/.local/lib/pythonX.Y/site-packages/ember

    where *X.Y* is your your installed version of Python, e.g 2.7.

    And similarly for Cantera, if it is not already on the Python path:

        $ ln -s /path/to/cantera/lib/pythonX.Y/site-packages/Cantera ~/.local/lib/pythonX.Y/site-packages/

    Or add parent directories of each of these modules to your ``PYTHONPATH``.

6. Prepare the documentation (optional)

        $ doxygen
        $ cd doc/sphinx
        $ make

    To view the HTML docs, open ``doc/sphinx/html/index.html`` in your web browser.

7. Prepare your kinetics mechanism. If you're starting with a mechanism that's
   in the Chemkin format, it will need to be converted:

        $ source /wherever/cantera/is/installed/setup_cantera
        $ ck2cti2 -input=<mech> --thermo=<thermo> --transport=<transport>

   This will produce a ``.cti`` file, which needs to be further converted:

        $ cti2ctml <mechname>.cti

   to produce ``<mechname>.xml``.

8. Prepare the input file, based on the one specified in

        ember/input/example-premixed.py

   A complete list of all available input parameters may be found in the HTML
   documentation. Alternatively, you look at the definitions in
   ``python/ember/input.py``. Specify the path to your mechanism file (which
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

    Using the ``ember.utils`` module:

        $ ipython --pylab
        In [1]: import ember
        In [2]: prof = ember.utils.load('run/test/profNow.h5')
        In [3]: plot(prof.x, prof.T)
