.. default-role:: math

.. _sec-boundary-conditions:

Boundary Conditions
===================

Establishing the boundary conditions for the elemental flame requires separate
consideration of each of the flame configurations previously described.

Single Opposed Jet
------------------

For the single opposed jet flame, the boundary conditions for the species and
energy equations consist of defining the temperature and mass fraction of the
two incoming streams. For the momentum equation, spatial gradients in `U` must
vanish on both sides of the flame, and we solve the simplified form of the
momentum equation:

.. math:: \rho\frac{\partial U}{\partial t}+
          \rho U^{2}-\rho_{\infty}\left(\frac{da}{dt}+a^{2}\right)=0

for each boundary point. For steady strain rates, this equation requires that
`U=a` on the reactants side and `U=a\sqrt{\rho_{u}/\rho_{b}}` on the products
side.

For a premixed flame with reactants supplied from `-\infty`, the boundary
conditions for the species and energy equations are:

.. math::
   \begin{array}{ccc}
   r=-\infty: & Y_{k}=Y_{k,u} & T=T_{u}\\
   r=+\infty: & Y_{k}=Y_{k,b} & T=T_{b}
   \end{array}

The burned gas temperature and composition correspond to the unburned mixture
brought to equilibrum at constant enthalpy and pressure. The boundary condition
for the continuity equationis taken at the left side of the computational
domain, where `V` is held fixed.

Central Control Volume
----------------------

For curved flames at finite stagnation radius, or flames specified in terms of
the unified formulation given in :ref:`sec-unified-formulation`, the `r=\infty`
boundary condition is the same as for the curved flame at zero stagnation
radius. The boundary conditions at `r=0`, however, require special attention.
First, the mass flux at the center must be expressed in terms of `rV` because
`V\rightarrow\infty` as `r\rightarrow0` in the potential flow solution. If we
specify the non-reacting stagnation point radius `R`, then the boundary
condition for the continuity equation is:

.. math:: r=0:\quad rV=\frac{1}{\alpha+1}\rho a\left|R\right|R

If the boundary mass flux `\left(rV\right)_{0}` is positive, indicating the
presence of source at the boundary, special care must be taken in specifying
the boundary conditions for the energy, species and momentum equations. While
the zero-gradient condition must still hold because of the symmetry at that
boundary, it is important to retain the effect of the mixture being introduced,
which may not be at the same state as the mixture in the vicinity of `r=0`. To
this end, we consider an integral, control volume approach to the `r=0`
boundary condition. Beginning with the species conservation equation, we
multiply through by `r^{\alpha}` and integrate from 0 to some small radius
`R_{i}`:

.. math:: \int_{0}^{R_{i}}r^{\alpha}\rho\frac{\partial Y_{k}}{\partial t}
          +r^{\alpha}V\frac{\partial Y_{k}}{\partial r}+
          \frac{\partial}{\partial r}\left[r^{\alpha}j_{k}\right]
          -r^{\alpha}\dot{\omega}_{k}W_{k}\, dr=0

Because this volume is small, we assume that variations of `Y_{k}`, `\rho` and
`\dot{\omega}_{k}` are negligible, so the unsteady term and the production term
may be taken out of the integral:

.. math:: \frac{R_{i}^{\alpha+1}}{\alpha+1}\left(\rho\frac{\partial Y_{k}}{\partial t}-\dot{\omega}_{k}W_{k}\right)+
          \int_{0}^{R_{i}}r^{\alpha}V\frac{\partial Y_{k}}{\partial r}+\frac{\partial}{\partial r}\left[r^{\alpha}j_{k}\right]\, dr=0

The convection term may be integrated by noting that variations in
`r^{\alpha}V` are negligible across this small distance. Furthermore, we
recognize that `Y_{k}|_{r=0}=Y_{k,left}` is the mass fraction corresponding to
the inlet mixture (either reactants or products). The diffusion term may also
be integrated, noting that `j_{k}|_{r=0}=0` by the symmetry condition. We now
have an ODE for the mass fraction of species `k` in the vicinity of the
symmetry boundary:

.. math:: \frac{R_{i}^{\alpha+1}}{\alpha+1}\left(\rho\frac{\partial Y_{k}}{\partial t}-\dot{\omega}_{k}W_{k}\right)+
          \left(r^{\alpha}V\right)_{0}\left(Y_{k}-Y_{k,left}\right)+R_{i}^{\alpha}j_{k}=0

Finally, we divide by the leading coefficient so that this equation scales
similarly to the species equation in the rest of the domain:

.. math:: \rho\frac{\partial Y_{k}}{\partial t}
          -\dot{\omega}_{k}W_{k}+\frac{\alpha+1}{R_{i}^{\alpha+1}}\left(r^{\alpha}V\right)_{0}\left(Y_{k}-Y_{k,left}\right)
          +\frac{\alpha+1}{R_{i}}j_{k}=0

A similar analysis for the energy equation yields:

.. math:: \rho\frac{\partial T}{\partial t}
          +\frac{1}{c_{p}}\sum_{k=1}^{K}\hat{h}_{k}\dot{\omega}_{k}
          +\frac{\alpha+1}{R_{i}^{\alpha+1}}\left(r^{\alpha}V\right)_{0}\left(T-T_{left}\right)
          -\frac{\alpha+1}{R_{i}c_{p}}\left(\lambda\frac{\partial T}{\partial r}\right)=0

In the case of the momentum equation, we neglect the term associated with the
boundary mass flux, and obtain:

.. math:: \rho\frac{\partial U}{\partial t}+\rho U^{2}
          -\rho_{\infty}\left(\frac{\partial a}{\partial t}+a^{2}\right)
          -\frac{\alpha+1}{R_{i}}\mu\frac{\partial U}{\partial r}=0

If the boundary mass flux `\left(rV\right)_{0}` is negative, the term including
it in each of these boundary conditions is eliminated.
