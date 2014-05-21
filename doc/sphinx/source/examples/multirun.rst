.. currentmodule:: ember.input

*********************************
Sequential Simulations (multirun)
*********************************

This example shows how to simulate a sequence of steady flames at a range of
strain rate values.

``example_multirun.py``:

.. literalinclude:: ../../../../python/ember/examples/example_multirun.py

The given configuration parameters override default values. The full set of
configuration parameters is described in :ref:`sec-conf-options`. Specific
parameters of interest in this example:

* :class:`General`

  - ``nThreads`` - By using the `multiprocessing.cpu_count())` function, we
    always use the maximum number of processors available. Note that for systems
    with hyperthreading enabled, this will use the total number of virtual
    cores, which may be less efficient than using the number of physical cores.

* :class:`Chemistry` - these alternative kinetic and transport models achieve
  increases in speed at a small cost in accuracy.

* :class:`StrainParameters`

  - ``rates`` - A list of strain rates to step through. For the supported
    opposed flow flame (this example), stepping through strain rates from high
    to low is usually more computationally efficient. For the unsupported or
    axisymmetric flame, the strain rates should be specified in increasing
    order, as the flame will extinguish at a sufficiently strain rate. For each
    strain rate in this list, the solver will generate a steady-state solution.

The output files for this run, saved in the directory specified by
``Paths.outputDir``, will contain the string ``epsNNNN`` where NNNN is the
strain rate, e.g. a run at a strain rate of 480 1/s will produce
``conf_eps0480``, ``log-eps0480.txt``, ``out_eps0480.h5``, and
``prof_eps0480.h5``.

Integral output data (strain rate, total heat release rate, flame consumption
speed, and flame position) for each strain rate will be written to ``integral.h5``.

To try this example, run the following command::

    $ python -m ember.examples.example_multirun

Or copy the ``example_multirun.py`` file to another directory, make any
modifications you desire, then run::

    $ python example_multirun.py
