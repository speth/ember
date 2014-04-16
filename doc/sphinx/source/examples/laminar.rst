.. currentmodule:: ember.input

*************
Laminar Flame
*************

This example shows how to simulate a freely-propagating, unstrained laminar
flame.

.. literalinclude:: ../../../../python/ember/examples/example_laminar.py

The following configuration parameters are particularly important for unstrained
flames:

  * :class:`StrainParameters` should have both `initial` and `final` set to 0.

  * :class:`General` should have `fixedLeftLocation` set to `True`. Since the
    boundary condition on the burned side of the flame is an outflow condition,
    a zero-gradient boundary condition is appropriate, so we set
    `fixedBurnedVal` to `False`.

  * :class:`PositionControl` - Since an unstrained flame propagating into a stagnant
    mixture has no steady-state solution in the fixed reference frame, we use
    this option to vary the inlet velocity to the domain so that we end up in a
    reference frame that follows the flame. The values of `xInitial` and
    `xFinal` should be equal, and set to a point within the initially specified
    domain. The `proportionalGain` value needs to be set high enough so that the
    flame doesn't move too far within the domain, but low enough to avoid
    numerical instabilities.
