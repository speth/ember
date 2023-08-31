************
Installation
************

Download
========

Releases
--------

Source code for stable releases of Ember may be downloaded from GitHub:
`<https://github.com/speth/ember/releases>`_

Git Repository
--------------

The current development version of Ember is available on Github:
`<https://github.com/speth/ember>`_. You can check out the source code as
follows on Unix-like systems::

    $ mkdir -p ~/src
    $ cd ~/src
    $ git clone http://github.com/speth/ember.git
    $ cd ember

Build Dependencies
==================

The following software needs to be installed in order to compile Ember. Listed versions
are what has been used for Ember development. Newer or older versions may not work. The
easiest way to install these dependencies is from the ``conda-forge`` Conda channel, using
the `environment.yaml <https://github.com/speth/ember/blob/main/environment.yaml>`_ file
included in the Ember repository.

* A C++17-compliant compiler (g++ >= 9.0, clang++ >= 5.0, Visual Studio 2017 or newer)
* Cantera (== 3.0.0)
* Python (>= 3.8, < 3.12)
* Cython (3.0.0)
* Boost (>= 1.70)

  .. note::

     Only the Boost headers are required. None of the compiled Boost
     libraries are used.

* SCons (4.5.2)
* Eigen (3.4.0)
* Sundials (>= 6.0)

  .. warning::

     Be sure to use the same version of Boost as the one that was used to
     compile Cantera.

* Intel Threading Building Blocks (2021.10.0). TBB is optional, but
  strongly recommended. Without TBB, Ember will only run on a single processor. To
  compile Ember without TBB use the SCons option ``use_tbb=n``.

* numpy (1.25.2)

* h5py (3.9.0, optional)

In order to run the unit tests:

* pytest (>= 7.0 recommended)

In order to process the documentation, the following are also required:

* Doxygen (>= 1.8.0)
* graphviz
* Sphinx
* Pygments
* pyparsing
* doxylink (see http://pypi.python.org/pypi/sphinxcontrib-doxylink)


Compilation
===========

First, check out a copy of the code as described above, or extract a release version of
the code from the downloaded archive. If all the necessary libraries are installed in
system-default locations, or in an active Conda environment, then the following should
be sufficient::

    $ scons build

Otherwise, specify the necessary paths using the corresponding command-line
options: "cantera", "sundials", "eigen", "boost", "tbb". For example::

    $ scons build sundials=/opt/sundials-6.6.0

If compilation succeeds, install the Ember Python module by running::

    $ scons install

By default, Ember will be installed into the active Python environment. Using an an
environment management tool such as ``venv`` or Conda is strongly recommended. If you wish
install to a different directory, use the ``install_args`` option to SCons to specify the
appropriate command to be passed to ``pip install``.

Alternatively, if you're doing active development, you can just add the module directory
to your ``PYTHONPATH``, for example:

    $ export PYTHONPATH=$HOME/src/ember/python

Documentation
-------------

To build the Ember documentation, starting from the ``ember`` directory, run::

    $ doxygen
    $ cd doc/sphinx
    $ make html

To view the HTML docs, open ``doc/sphinx/html/index.html`` in your web browser.
