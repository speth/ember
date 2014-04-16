.. currentmodule:: ember.input

***********************
Premixed Strained Flame
***********************

This example shows how to simulate a steady premixed strained flame. This is
essentially the "default" flame configuration for Ember, and most other flame
configurations are described by how they differ from this one.

``example_premixed.py``:

.. literalinclude:: ../../../../python/ember/examples/example_premixed.py

The given configuration parameters override default values. The full set of
configuration parameters is described in :ref:`sec-conf-options`.

* :class:`Paths`

  - ``outputDir`` is the directory where periodic periodic
    output profiles (`prof000000.h5`, `prof0000001.h5`, etc.) and time series
    output files (`out.h5`) will be saved. If a relative path is given (as in
    this example), the path is relative to the directory from which the script
    is run.

* :class:`InitialCondition` In this example, the molar fuel composition is
  specified using the string format supported by Cantera. Alternatively, a
  vector of mole fractions can be supplied. Here, the oxidizer composition is
  left with at the default value, which corresponds to air, i.e.
  `"N2:3.76, O2:1.0"`.

* :class:`StrainParameters` specifies how the strain rate imposed on the flame
  varies as a function of time. In this case we specify a constant strain rate.

* :class:`TerminationCondition` specifies the conditions under which the time
  integration will be terminated. Here, we specify integration until a specific
  time. Another option is to integrate until a steady-state solution is reached.

To try this example, run the following command::

    $ python -m ember.examples.example_premixed

Or copy the ``example_premixed.py`` file to another directory, make any
modifications you desire, then run::

    $ python example_premixed.py
