.. currentmodule:: ember.input

***************************
Counterflow Diffusion Flame
***************************

This example shows how to simulate an opposed flow diffusion flame.

``example_diffusion.py``:

.. literalinclude:: ../../../../python/ember/examples/example_diffusion.py

The given configuration parameters override default values. The full set of
configuration parameters is described in :ref:`sec-conf-options`.

* :class:`Paths`

  - ``outputDir`` is the directory where periodic periodic
    output profiles (`prof000000.h5`, `prof0000001.h5`, etc.) and time series
    output files (`out.h5`) will be saved. If a relative path is given (as in
    this example), the path is relative to the directory from which the script
    is run.

  - ``logFile`` is used to write information about the progress of the solver.
    The path is relative to the directory from which the script is run, i.e. it
    is not necessarily placed in ``outputDir``. If no value is specified for
    this option, the log output is written to the screen via ``stdout``.

* :class:`InitialCondition` - For the diffusion flame, properties need to be
  for the separate fuel and oxidizer streams.

  - ``flameType='diffusion'`` specifies that this is a diffusion flame and that
    independent properties for the fuel and oxidizer mixture should be used.

  - ``fuel`` and ``oxidizer`` define the composition of the fuel and oxidizer
    streams, respectively.

  - ``Tfuel`` and ``Toxidizer`` define the temperature of the fuel and oxidizer
    streams, respectively.

  - ``xLeft`` and ``xRight`` define the extent of the initial domain. Depending
    on the grid adaptation parameters, the extent of the domain may change as
    the flame evolves in time.

  - ``centerWidth`` defines the width in the initial profile of the fully-mixed
    region between the two inlet streams. If this region is not wide enough,
    the flame may extinguish unexpectedly. Otherwise, the steady-state solution
    will be independent of this value.

  - ``slopeWidth`` defines the initial width of the transition between the inlet
    streams and the center zone. Making this too narrow wide will increase the
    time it takes to reach a steady-state flame. Making this too wide may result
    in integration difficulties.

* :class:`StrainParameters` specifies how the strain rate imposed on the flame
  varies as a function of time. In this case we specify a constant strain rate.

* :class:`General`

  - ``nThreads`` - Set the number of parallel threads to use when solving. This
    can be set to any value up to the number of processors in the computer.

* :class:`Times`

  - ``globalTimestep`` - sets the timestep size used by the operator-split
    integrator. Decreasing this value may help in cases where numerical
    oscillations are observed. Typically, a value of `2e-5` can be used without
    trouble.

  - ``profileStepInterval`` - specifies that the output profiles
    (`prof000001.h5` etc.) should be produced every 20 global timesteps.

* :class:`TerminationCondition` specifies the conditions under which the time
  integration will be terminated. Here, we specify integration until a specific
  time. Another option is to integrate until a steady-state solution is reached.

To try this example, run the following command::

    $ python -m ember.examples.example_diffusion

Or copy the ``example_diffusion.py`` file to another directory, make any
modifications you desire, then run::

    $ python example_diffusion.py
