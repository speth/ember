.. Ember documentation master file, created by
   sphinx-quickstart on Tue Jun  5 17:49:07 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*****
Ember
*****

Ember is a quasi-one-dimensional, unsteady reacting flow solver. It can be used
to simulate a number of fundamental flame configurations, including premixed
laminar flames, opposed flow strained flames (premixed or diffusion),
axisymmetric (tubular) flames with positive or negative curvature, and steady 2D
flames in a prescribed velocity field (using the method of lines).

Ember integrates the governing equations using a variant on the standard Strang
splitting method called *rebalanced splitting*, which eliminates steady-state
errors.

Ember's solver is implemented primarily in C++, which is made accessible to the
user as a Python module to enable script-driven input files and provide
extensibility through user-supplied functions. Ember is parallelized using the
Intel TBB library to take advantage of modern multi-core processors.

This site contains documentation for the current development version of Ember,
which may differ slightly from the most recent stable release.

Contents
========

.. toctree::
   :maxdepth: 2

   installation
   flame-model
   examples/index
   input
   run-sim
   postprocessing
   cxx-implementation

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

