.. |x| replace:: :math:`x`
.. |z| replace:: :math:`z`
.. |r| replace:: :math:`r`
.. |u| replace:: :math:`u`
.. |v| replace:: :math:`v`

.. _sec-flame-model:

***********
Flame Model
***********

The elemental flame model consists of a laminar flame stabilized in a
stagnation flow. The parameters of the stagnation flow may be imposed
explicitly, or provided through coupling with an outer flow within which the
elemental flame is embedded. The kinematics of the coupling between the
elemental flame and the outer flow have been described elsewhere [#Marzouk1999].
Here, the flame is parameterized by the strain rate and position (radius) of
the stagnation point in the non-reacting flow. The stagnation flow models are
described in :ref:`sec-strained-flame-configurations`.

The general 3D conservation equations for mass, momentum, energy and chemical
species are reduced to one dimension by using a boundary layer approximation
across the flame. The derivation of these governing equations is contained in
Section :ref:`sec-governing-equations`. The boundary conditions for the
elemental flame model are discussed in :ref:`sec-boundary-conditions`.

.. _sec-strained-flame-configurations:

Strained Flame Configurations
=============================

These are the models of strained flames implemented in Ember.

.. _sec-opposed-jet:

Single Opposed-Jet Flame
------------------------

.. figure:: /_static/images/planar-single.png
   :align: center

   **Schematic of the single opposed-jet flame**

The simplest configuration for the elemental flame is that of the single
opposed-jet flame [#Marzouk1999], [#Marzouk2003]. Fluid flows from
:math:`x=\pm\infty` toward the stagnation point at :math:`x=0`. The potential
flow velocity field for the non-reacting flow is

.. math:: u = a(t)z \quad\quad v = -a(t) x

where |u| and |v| are the velocity components in the |z| and |x| directions,
respectively, and :math:`a(t)` is the time varying strain rate parameter. If we
take the mixture coming from :math:`x=-\infty` to be premixed reactants and the
mixture coming from :math:`x=+\infty` to be the reaction products, a flame will
be established at some point :math:`$x<0` where the consumption speed balances
the local fluid velocity. Note that a flame in this configuration can never be
completely extinguished because there will always be a region of contact
between reactants and products.

.. _sec-twin-flame:

Twin Opposed-Jet Flame
----------------------

.. figure:: /_static/images/planar-twin.png
   :align: center

   **Schematic of the twin opposed-jet flame**

The potential flow velocity field for the twin opposed-jet flame is identical
to that of the single opposed jet flame, but with both inlet streams consisting
of premixed reactants. As depicted above, two flames are established, one on
each side of the stagnation point, propagating away from one another.
Combustion products exist in the region between the two flames. Because the
composition of the products mixture is not imposed externally, this
configuration allows flame extinction to occur when the flames are pushed
close together (at high strain rates). The domain is symmetric about the
stagnation plane :math:`x=0`, allowing this configuration to be simulated
numerically as a single flame with a symmetry boundary condition applied at
the stagnation point.

.. _sec-tubular0:

Tubular Flame at Zero Stagnation Radius
---------------------------------------

.. figure:: /_static/images/curved0.png
   :align: center

   **Schematic of the tubular flame at zero stagnation radius**

The curved flame at zero stagnation radius is quite similar to the twin
opposed-jet flame. Instead of the planar stagnation flow, we have an
axisymmetric stagnation flow, with reactants flowing inward from
:math:`r=\infty` toward the stagnation point at :math:`r=0`, establishing a
cylindrical flame propagating away from the centerline.  The potential flow
velocity field is

.. math:: u = a(t) z \quad\quad v = -\frac{1}{2} a(t) r

where |u| and |v| are the velocity components in the |z| and |r| directions,
respectively. As in the twin opposed jet flame, a symmetry boundary condition
applies at :math:`r=0`, and flame extinction occurs when the flame is pushed
too close to the centerline.

.. _sec-tubular-finite:

Tubular Flame at Finite Stagnation Radius
-----------------------------------------

.. figure:: /_static/images/curved-finite.png

The tubular flame at zero stagnation radius can be modified by introducing a
line source at :math:`r=0`. This produces a potential flow where the stagnation
surface is a cylinder of radius :math:`R`, as shown above. Since the composition
of the flow coming from :math:`r=0` may be freely specified, this configuration
permits both positively and negatively curved flames. The resulting potential
flow velocity field is:

.. math:: u=az\quad\quad v=\frac{a}{2}\left(\frac{R^{2}-r^{2}}{r}\right)

The boundary condition at :math:`r=0` requires special attention in order to
capture both the symmetry condition and the effect of the mass flux at the
centerline.

.. _sec-unified-formulation:

Unified Formulation
-------------------

The stagnation point flow configurations discussed above can be rewritten in a
unified form, permitting both planar and curved geometry through the
introduction of the parameter :math:`\alpha` where :math:`\alpha=0` for planar
geometry and :math:`\alpha=1` for curved geometry. We introduce the coordinate
system :math:`(r, z)` where |r| is the coordinate normal to the flame
(regardless of whether the coordinate system is curved or planar) and the
|z| coordinate is tangential to the flame. The potential flow velocity field
is then

.. math:: u = az \quad\quad
          v=\frac{a}{\alpha+1} \left(\frac{\left|R\right|R^{\alpha}
                                     -r^{\alpha+1}}{r^{\alpha}}\right)

This formulation permits the use of a fictitious negative stagnation point
radius :math:`R`, with the centerline acting as a sink rather than a source.

The stretch rate :math:`\kappa` for a flame at radius :math:`R_{f}` is:

.. math:: \kappa=a+\frac{\alpha}{R_{f}}\frac{dR_{f}}{dt}

When the flame is stationary, the stretch rate reduces to :math:`\kappa=a` and
thus curvature does not contribute to flame stretch for stationary flames in
this configuration.

.. _sec-governing-equations:

Governing Equations
===================

.. _sec-boundary-conditions:

Boundary Conditions
===================

