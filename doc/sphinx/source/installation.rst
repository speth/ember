************
Installation
************

Download
========

Releases
--------

Source code and binaries for stable releases of Ember may be downloaded from
Github: `<https://github.com/speth/ember/releases>`_

Git Repository
--------------

The current development version of Ember is available on Github:
`<https://github.com/speth/ember>`_. You can check out the source code as
follows on Unix-like systems::

    $ mkdir -p ~/src
    $ cd ~/src
    $ git clone http://github.com/speth/ember.git
    $ cd ember


Binary Dependencies
===================

Binary releases of Ember are currently available for 64-bit Python 2.7 on
Windows. To use the binary version of Ember, you must have the following
installed:

    * `Python 2.7 (64-bit) <https://www.python.org/downloads/>`_
    * `Numpy <http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy>`_ - version 1.8 or higher
    * `Cantera 2.1 Python module <https://sourceforge.net/projects/cantera/files/cantera/>`_

        - For Ember 1.2.x, install Cantera 2.1 and the "Legacy" Cantera Python
          module
        - For Ember 1.3 and higher, install the new Python module

    * `h5py <http://www.lfd.uci.edu/~gohlke/pythonlibs/#h5py>`_


Build Dependencies
==================

The following software needs to be installed in order to compile Ember:

    * A C++ compiler (g++, Clang, Intel C++, Visual Studio 2008)
    * Python 2.6 or higher (2.7 recommended)
    * Cython (>= 0.18)
    * Boost (>= 1.40, including the compiled ``boost_thread`` library)
    * SCons (>= 2.1.0 recommended)
    * Eigen (>= 3.0)
    * Cantera (>= 2.1.0)

        - For Ember 1.2.x, install Cantera with the 'legacy' Python module
        - For Ember 1.3 and higher, install the 'new' Python module

    * Sundials (2.4.0 or 2.5.0)
    * Intel Threading Building Blocks (>= 4.0)
    * numpy (>= 1.3.0)
    * h5py (>= 1.2.1)

In order to process the documentation, the following are also required:

    * Doxygen (>= 1.8.0)
    * Sphinx
    * Pygments
    * pyparsing
    * doxylink (see http://pypi.python.org/pypi/sphinxcontrib-doxylink)


Compilation
===========

First, check out a copy of the code as described above, or extract a release
version of the code from the downloaded archive. If all the necessary libraries
are installed in system   default locations, then the following should be
sufficient::

    $ scons build

Otherwise, specify the necessary paths using the corresponding command-line
options: "cantera", "sundials", "eigen", "boost", "tbb". For example::

    $ scons build cantera=/home/$USER/.local sundials=/opt/sundials-2.5.0

If compilation succeeds, install the Ember Python module by running::

    $ scons install

By default, Ember will be installed to your *user* site-packages directory, e.g.
`~/.local/lib/python2.7/site-packages/ember`. If you wish install to a different
directory, use the `install_args` option to SCons to specify the appropriate
command to be passed to `python setup.py install`.

Alternatively, if you're doing active development, you can just create a symlink
into your user Python module path::

    $ mkdir -p ~/.local/lib/pythonX.Y/site-packages
    $ ln -s ~/src/ember/build/python/ember ~/.local/lib/pythonX.Y/site-packages/ember

where *X.Y* is your your installed version of Python, e.g 2.7.

And similarly for Cantera, if it is not already on the Python path::

    $ ln -s /path/to/cantera/lib/pythonX.Y/site-packages/cantera ~/.local/lib/pythonX.Y/site-packages/

Or add parent directories of each of these modules to your ``PYTHONPATH``::

    $ export PYTHONPATH=/path/to/cantera/lib/python2.7/site-packages:~/src/ember/build/python


Documentation
-------------

To build the Ember documentation, starting from the ``ember`` directory, run::

    $ doxygen
    $ cd doc/sphinx
    $ make html

To view the HTML docs, open ``doc/sphinx/html/index.html`` in your web browser.
