C++ Implementation Docs
=======================

The one-dimensonal flame model is primarily implemented by the following C++
classes:

* :cxx:`FlameSolver` - Main integrator driver

* Operator split systems of equations:

  * :cxx:`SourceSystem` - chemical source (reaction) term

    * :cxx:`SourceSystemCVODE` - used with *sundialsCVODE*
    * :cxx:`SourceSystemQSS` - used with *QssIntegrator*

  * :cxx:`DiffusionSystem` - diffusion term for a single component

  * :cxx:`ConvectionSystemSplit` - Wrapper that integrates convection terms for
    all components

    * :cxx:`ConvectionSystemY` - convection term for a single species
    * :cxx:`ConvectionSystemUTW` - convection term for velocity, temperature,
      and molecular weight

* Integrators:

  * :cxx:`sundialsCVODE` - Integrator used for convection systems and
    *SourceSystemCVODE*

  * :cxx:`QssIntegrator` - Integrator used with *SourceSystemQSS*

* Other classes of interest:

  * :cxx:`OneDimGrid` - an adaptive, non-uniform one-dimensional grid

  * :cxx:`GridBased` - base class used with several other classes to provide
    access to the grid parameters

  * :cxx:`configOptions` - options used to set up the *FlameSolver*

  * :cxx:`CanteraGas` - Wrapper around a set of Cantera objects used to compute
    thermodynamic, kinetic, and transport properties of an ideal gas mixture.

For the full list of classes, see the `C++ Class List <../../doxygen/html/annotated.html>`_.
