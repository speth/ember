=====================================
Ember: unsteady strained flame solver
=====================================

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1004753.svg
   :target: https://doi.org/10.5281/zenodo.1004753

Introduction
------------

Ember is a quasi-one-dimensional, unsteady reacting flow solver. It can be used
to simulate a number of fundamental flame configurations, including premixed
laminar flames, opposed flow strained flames (premixed or diffusion),
axisymmetric (tubular) flames with positive or negative curvature, and steady 2D
flames in a prescribed velocity field (using the method of lines).

Ember integrates the governing equations using a variant on the standard Strang
splitting method which eliminates steady-state errors.

Online documentation for Ember is located at `<https://speth.github.io/ember-doc>`_.

Installation from source
------------------------

The easiest way to compile and install Ember is to use `conda-forge` to provide all
dependencies beyond the base C++ compiler. After installing a Conda distribution such as
`Miniforge <https://github.com/conda-forge/miniforge/releases>`_, you can install these
dependencies in a new environment named `ember-build` by running::

    $ mamba env create -n ember-build -f environment.yaml

from the Ember source directory (that is, the directory containing this README).

Then, compile, test, and install Ember by running::

    $ scons build
    $ scons test
    $ scons install

Further `Installation Instructions <https://speth.github.io/ember-doc/sphinx/html/installation.html>`_
are available in the Ember documentation.

Running Ember
-------------

* Prepare your kinetics mechanism. If you're starting with a mechanism that's
  in the Chemkin format, it will need to be converted::

    $ ck2yaml --input=<mech> --thermo=<thermo> --transport=<transport>

  This will produce a``<mech>.yaml`` file that can be used by Cantera and Ember.

* Prepare the input file, based on the one specified in::

    python/ember/examples/example-premixed.py

  A complete list of all available input parameters may be found in the `HTML
  documentation. <https://speth.github.io/ember-doc/sphinx/html/input.html>`_
  Alternatively, you look at the definitions in
  ``python/ember/input.py``. Specify the path to your mechanism file as the
  ``mechanismFile``. Other parameters you may want to change are in the
  ``Paths``, ``InitialCondition``, and ``StrainParameters`` sections.

* Check the configuration for errors::

    $ python myInputFile.py validate

  If this prints "Validation completed successfully.", you're all set.
  Otherwise, try to correct the indicated error and try again.

* Run the code::

    $ python myInputFile.py &

  This may take a while. You can watch the solver's progress as it is written to
  the file specified by ``Paths.logFile`` in the input file.

* Examine the output files. The files are HDF5 data files (file extension
  ``.h5``), which can be read using the Python ``h5py`` module or Matlab, or
  compressed NumPy data files (file extension ``.npz``), which can be read using
  NumPy.

  * ``out.h5`` or ``out.npz`` contains integral flame properties (e.g. flame
    speed) as a function of time
  * ``profNNNNNN.h5`` or ``profNNNNNN.npz`` contain the temperature & species
    profiles output periodically.
  * ``profNow.h5`` or ``profNow.npz`` contains the most recently saved profiles.

  Using the ``ember.utils`` module (and assuming you have IPython and matplotlib
  installed), either of these file types can be read using the ``load``
  function::

    $ ipython --pylab
    In [1]: import ember
    In [2]: prof = ember.utils.load('run/test/profNow.h5')
    In [3]: plot(prof.x, prof.T)
