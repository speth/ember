Getting Started with Ember
==========================

Ember is a quasi-one-dimensional, unsteady reacting flow solver. It can be used
to simulate a number of fundamental flame configurations, including premixed
laminar flames, opposed flow strained flames (premixed or diffusion),
axisymmetric (tubular) flames with positive or negative curvature, and steady 2D
flames in a prescribed velocity field (using the method of lines).

Ember integrates the governing equations using a variant on the standard Strang
splitting method which eliminates steady-state errors.

Online documentation for Ember is located at `<http://speth.github.com/ember-doc>`_.

Compiling Ember
---------------

* Build Dependencies: These directions assume you are using a Linux or Windows
  system with the following installed:

  * A C++ compiler (g++ or Visual Studio 2008)
  * Python (2.6 or 2.7)
  * Cython (>= 0.18)
  * Boost (>= 1.40, including boost-filesystem)
  * SCons (2.1.0 recommended)
  * Eigen (>= 3.0)
  * Cantera (>= 2.0.0)
  * Sundials (== 2.4)
  * HDF5 (>= 1.8.0)
  * Intel Threading Building Blocks (>= 4.0)
  * numpy (>= 1.3.0)
  * h5py (>= 1.2.1)
  * git

  Additional dependencies for processing documentation:

  * Doxygen (>= 1.8.0)
  * Sphinx
  * Pygments
  * pyparsing
  * doxylink (see http://pypi.python.org/pypi/sphinxcontrib-doxylink)

* The code is stored in a repository on Github. To check out a copy of the
  code, run::

     $ cd
     $ mkdir -p src
     $ cd src
     $ git clone git://github.com/speth/ember.git
     $ cd ember

* Compiling the code. If all the necessary libraries are installed in system
  default locations, then the following should be sufficient::

    $ scons build

  Otherwise, specify the necessary paths using the corresponding command-line
  options: "cantera", "sundials", "eigen", "boost", "hdf5", "tbb". For example::

    $ scons cantera=/home/$USER/.local sundials=/opt/sundials-2.4.0

  This produces the Python extension module ``build/python/ember/_ember.so``.

* Install the "Ember" Python module::

    $ scons install

  By default, Ember will be installed to your *user* site-packages directory,
  e.g. `~/.local/lib/python2.7/site-packages/ember`. If you wish install to a
  different directory, use the `install_args` option to SCons to specify the
  appropriate command to be passed to `python setup.py install`.

  Alternatively, if you're doing active development, you can just create a
  symlink into your user Python module path::

    $ mkdir -p ~/.local/lib/pythonX.Y/site-packages
    $ ln -s ~/src/ember/build/python/ember ~/.local/lib/pythonX.Y/site-packages/ember

  where *X.Y* is your your installed version of Python, e.g 2.7.

  And similarly for Cantera, if it is not already on the Python path::

    $ ln -s /path/to/cantera/lib/pythonX.Y/site-packages/Cantera ~/.local/lib/pythonX.Y/site-packages/

  Or add parent directories of each of these modules to your ``PYTHONPATH``::

    $ export PYTHONPATH=/path/to/cantera/lib/python2.7/site-packages:~/src/ember/build/python

* Prepare the documentation (optional)::

    $ doxygen
    $ cd doc/sphinx
    $ make html

  To view the HTML docs, open ``doc/sphinx/html/index.html`` in your web browser.

Running Ember
-------------

* Prepare your kinetics mechanism. If you're starting with a mechanism that's
  in the Chemkin format, it will need to be converted::

    $ source /wherever/cantera/is/installed/setup_cantera
    $ ck2cti2 -input=<mech> --thermo=<thermo> --transport=<transport>

  This will produce a ``.cti`` file, which needs to be further converted::

    $ cti2ctml <mechname>.cti

  to produce ``<mechname>.xml``.

* Prepare the input file, based on the one specified in::

    ember/input/example-premixed.py

  A complete list of all available input parameters may be found in the `HTML
  documentation. <http://speth.github.com/ember-doc/sphinx/html/input.html>`_
  Alternatively, you look at the definitions in
  ``python/ember/input.py``. Specify the path to your mechanism file (which
  needs to be in whatever directory is specified as the "input" directory in
  input file) as the ``mechanismFile`` and ``"gas"`` as the ``phaseID``. The
  other parameters you may want to change are in the ``Paths``,
  ``InitialCondition``, and ``StrainParameters`` sections.

* Check the configuration for errors::

    $ python myInputFile.py validate

  If this prints "Validation completed successfully.", you're all set.
  Otherwise, try to correct the indicated error and try again.

* Run the code::

    $ python myInputFile.py &

  This may take a while. You can watch the solver's progress as it is written to
  the file specified by ``Paths.logFile`` in the input file.

* Examine the output files. The files are HDF5 data files, which can be read
  using the Python ``h5py`` module or Matlab.

  * ``out.h5`` contains integral flame properties (e.g. flame speed) as a
    function of time
  * ``profNNNNNN.h5`` contain the temperature & species profiles output
    periodically.
  * ``profNow.h5`` contains the most recently saved profiles.

  Using the ``ember.utils`` module (and assuming you have IPython and
  matplotlib installed)::

    $ ipython --pylab
    In [1]: import ember
    In [2]: prof = ember.utils.load('run/test/profNow.h5')
    In [3]: plot(prof.x, prof.T)
