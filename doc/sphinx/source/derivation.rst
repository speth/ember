.. default-role:: math

.. _sec-governing-equations:

Governing Equations
===================

The one-dimensional governing equations for the elemental flame are derived
here. We begin with the general 3D governing equations for reacting flow as
given by Kee (Kee2003) and reduce them to a single dimension normal to the
flame using a boundary layer approximation and solving along the stagnation
streamline `z=0`.

Consider the coordinate system `(z, r, \theta)` with velocity components
`(u,v,w)` and let `\hat{r}` be the flame normal. We may allow this coordinate
system to be either Cartesian or cylindrical through the introduction of a
parameter `\alpha`, where `\alpha=1` for the cylindrical case and `\alpha=0`
for the Cartesian case. In making the boundary layer approximation, the
tangential variations of all quantities except the pressure `p` and tangential
velocity `u` are neglected. The velocity and all variations in the `\theta`
direction are zero. For an arbitrary scalar `F`, the gradient and substantial
derivative are defined as:

.. math::
   :label: gradient

        \nabla F\equiv\hat{z}\frac{\partial F}{\partial z}
        + \hat{r}\frac{\partial F}{\partial r}

.. math::
   :label: subst-deriv

          \frac{DF}{Dt}=\frac{\partial F}{\partial t}+{\bf v}\cdot\nabla F
              =\frac{\partial F}{\partial t}+u\frac{\partial F}{\partial z}
              +v\frac{\partial F}{\partial r}

The divergence of a vector `\mathbf{F}` is defined as:

.. math::
   :label: divergence

   \nabla\cdot\mathbf{F}\equiv\frac{\partial F_{z}}{\partial z}
   + \frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left(r^{\alpha}F_{r}\right)

where `F_{z}` and `F_{r}` are respectively the `z` and `r` components of
`\mathbf{F}`.

Momentum Equation
-----------------

First, we will focus on the momentum equation. In general three-dimensional
form, it is:

.. math:: \rho\frac{D\mathbf{v}}{Dt}={\bf f}+\nabla\cdot\mathsf{T}

where `\mathsf{T}` is the stress tensor:

.. math:: \mathsf{T}=\left(\begin{array}{ccc}
             -p+2\mu\frac{\partial u}{\partial z}+\kappa\nabla\cdot\mathbf{v} &
             \mu\left(\frac{\partial u}{\partial r}+\frac{\partial v}{\partial z}\right) &
             \mu\left(\frac{1}{r^{\alpha}}\frac{\partial u}{\partial\theta}+\frac{\partial w}{\partial z}\right)\\
             \mu\left(\frac{\partial u}{\partial r}+\frac{\partial v}{\partial z}\right) &
             -p+2\mu\frac{\partial v}{\partial r} +\kappa\nabla\cdot\mathbf{v} &
             \mu\left(\frac{\partial w}{\partial r}-\alpha\frac{w}{r}+\frac{1}{r^{\alpha}}\frac{\partial v}{\partial\theta}\right)\\
             \mu\left(\frac{1}{r^{\alpha}}\frac{\partial u}{\partial\theta}+\frac{\partial w}{\partial z}\right) &
             \mu\left(\frac{\partial w}{\partial r}-\alpha\frac{w}{r}+\frac{1}{r^{\alpha}}\frac{\partial v}{\partial\theta}\right) &
             -p+2\mu\left(\frac{1}{r^{\alpha}}\frac{\partial w}{\partial\theta}+\alpha\frac{v}{r}\right)+\kappa\nabla\cdot\mathbf{v}
          \end{array}\right)

and the body force `\mathbf{f}=\mathbf{0}`. Here, `p` is the pressure and `\mu`
is the dynamic viscosity of the mixture. Setting velocities and derivatives in
the `\theta` direction to zero, the `z`-momentum equation is:

.. math:: \rho\frac{\partial u}{\partial t} +
          \rho u\frac{\partial u}{\partial z}+\rho v\frac{\partial u}{\partial r} =
          \frac{\partial}{\partial z}\left[-p+2\mu\frac{\partial u}{\partial z}
          +\kappa\nabla\cdot\mathbf{v}\right]+
          \frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\mu\left(\frac{\partial u}{\partial r}
          +\frac{\partial v}{\partial z}\right)\right]

The term containing the second coefficient of viscosity, `\kappa\nabla\cdot\mathbf{v}`,
is taken to be zero. Additionally, `\partial v/\partial z` and `\partial^{2}u/\partial z^{2}`
are neglected by the boundary layer approximation, giving the simplified
`z`-momentum equation:

.. math:: \rho\frac{\partial u}{\partial t} +
          \rho u\frac{\partial u}{\partial z} +
          \rho v\frac{\partial u}{\partial r}=
          -\frac{\partial p}{\partial z}
          +\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\mu\frac{\partial u}{\partial r}\right]

The pressure gradient outside the boundary layer may be obtained by substituting
the potential flow solution for `u` and `v` into the momentum equation:

.. math:: \rho_{\infty}z\frac{da}{dt}+\rho_{\infty}a^{2}z=
          -\frac{\partial p}{\partial z}

where `\rho_{\infty}` is the density of the reactants mixture. By the boundary
layer approximation, this must be the pressure gradient inside the boundary
layer as well.

The `z` dependence of `u` may be found similarly:

.. math:: \rho\frac{\partial u}{\partial t}+\rho u\frac{\partial u}{\partial z}+
          \rho v\frac{\partial u}{\partial r}=\rho_{\infty}z\frac{da}{dt}+
          \rho_{\infty}a^{2}z+\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\mu\frac{\partial u}{\partial r}\right]

Now, introduce the notation `U\equiv ua/u_{\infty}`. Using the potential flow
velocity field, we can then write:

.. math:: u=\frac{Uu_{\infty}}{a}=\frac{Uaz}{a}=Uz

Substituting this into the momentum equation gives:

.. math:: \rho z\frac{\partial U}{\partial t}
          +\rho Uz\left(z\frac{\partial U}{\partial z}
          +U\right)+\rho vz\frac{\partial U}{\partial r}
          =\rho_{\infty}z\frac{da}{dt}+\rho_{\infty}a^{2}z
          +\frac{z}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\mu\frac{\partial U}{\partial r}\right]

Dividing by `z` and solving along the stagnation streamline `z=0`, the momentum
equation simplifies to:

.. math:: \rho\frac{\partial U}{\partial t}+\rho U^{2}
          +\rho v\frac{\partial U}{\partial r}=
          \rho_{\infty}\frac{da}{dt}+\rho_{\infty}a^{2}
          +\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\mu\frac{\partial U}{\partial r}\right]

Continuity Equation
-------------------

Now consider the mass conservation equation:

.. math:: \frac{\partial\rho}{\partial t}+\nabla\cdot\left(\rho{\bf v}\right)=0

Expanding the divergence using Equation :eq:`divergence` gives:

.. math:: \frac{\partial\rho}{\partial t}
          +\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left(r^{\alpha}\rho v\right)
          +\rho\frac{\partial u}{\partial z}=0

Making the substitution for the similarity variable `U`, the mass conservation
equation becomes:

.. math:: \frac{\partial\rho}{\partial t}
          +\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left(r^{\alpha}\rho v\right)+\rho U=0

Species Equation
----------------

The general form of the species continuity equation is:

.. math:: \rho\frac{DY_{k}}{Dt}=-\nabla\cdot\mathbf{j}_{k}+\dot{\omega}_{k}W_{k}

where `Y_{k}` is the mass fraction of species `k`, `\dot{\omega}_{k}` is the
molar production rate of species `k`, `W_{k}` is the molecular weight of species
`k`, and the diffusion mass flux `\mathbf{j}_{k}` is defined as:

.. math:: \mathbf{j}_{k} =  -\rho D_{km} \nabla Y_{k}
                            -\frac{D_{k}^{T}}{T}\nabla T+Y_{k}\mathbf{j}'

Here, `T` is the temperature and `D_{km}` and `D_{k}^{T}` are respectively
the mass-based mixture-averaged diffusion coefficient and the thermal diffusion
coefficient of species `k`. Note that this definition of the diffusion mass flux
includes the thermal diffusion (Soret) effect. The final term introduces a
correction, `\mathbf{j}'`, which corrects for the inaccuracy of mixture-averaged
model so that the requirement `\sum\mathbf{j}_{k}=0` is satisfied. In order to
calculate `\mathbf{j}'`, we first calculate the diffusion mass fluxes ignoring
its contribution:

.. math:: \mathbf{j}_{k}^{*}=-\rho D_{km} \nabla Y_{k}
          -\frac{D_{k}^{T}}{T}\nabla T

Then, `\mathbf{j}'` is calculated as:

.. math::  \mathbf{j}'=-\sum_{k}\mathbf{j}_{k}^{*}

Substituting the gradient, substantial derivative and divergence as defined in
equations :eq:`gradient`, :eq:`subst-deriv` and :eq:`divergence`, respectively,
and setting the `z`-derivatives to zero, the species equation becomes:

.. math::
   :label: species

     \rho\frac{\partial Y_{k}}{\partial t}+V\frac{\partial Y_{k}}{\partial r}=
     -\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}j_{k}\right]
     +\dot{\omega}_{k}W_{k}

and the diffusion mass flux is:

.. math:: j_{k}=-\rho D_{km}\frac{\partial Y_{k}}{\partial r} -
          \frac{D_{k}^{T}}{T}\frac{\partial T}{\partial r}+Y_{k}j'

Energy Equation
---------------

Finally we turn our attention to the energy conservation equation. The general
form of the energy equation, expressed in terms of the enthalpy, is:

.. math:: \rho\frac{Dh}{Dt}=\frac{Dp}{Dt}-\nabla\cdot\mathbf{q}+\Phi

where `h` is the enthalpy, `q` is the heat flux and `\Phi` is the viscous work.
By the zero-Mach-number assumption, we neglect the effect of pressure variations,
so `Dp/Dt=0`. The viscous work term is also assumed to be much smaller than the
energy released by chemical reactions and is therefore neglected. Expanding the
substantial derivative of the enthalpy in terms of `T`, `Y_{k}` and specific
heat capacity `c_{p}` gives:

.. math:: \frac{Dh}{Dt}=\sum_{k=1}^{K}\left(Y_{k}\frac{Dh_{k}}{Dt}+h_{k}\frac{DY_{k}}{Dt}\right)=
          c_{p}\frac{DT}{Dt}+\sum_{k=1}^{K}h_{k}\frac{DY_{k}}{Dt}

The total number of species is `K`. The substantial derivative of `Y_{k}` is
replaced using the species equation :eq:`species`. With these substitutions,
the energy equation then becomes:

.. math:: \rho c_{p}\frac{DT}{Dt}
          +\sum_{k=1}^{K}h_{k}\left(-\nabla\cdot{\bf j}_{k}
          +\dot{\omega}_{k}W_{k}\right)=-\nabla\cdot\mathbf{q}

The heat flux vector `\mathbf{q}` is:

.. math:: \mathbf{q}=-\lambda\nabla T+\sum_{k=1}^{K}h_{k}j_{k}+\mathbf{q}_{r}

where `\lambda` is the thermal conductivity of the mixture and the radiation
heat flux `\mathbf{q}_{r}` is neglected in the present model. Also neglected is
the Dufour heat flux, proportional to the gradients in species concentrations,
which is typically three orders of magnitude smaller than the Fourier heat flux,
`\lambda\nabla T`. Substituting the heat flux vector into the energy equation,
and noting that `\hat{h}_{k}=h_{k}W_{k}` where `\hat{h}_{k}` is the molar
enthalpy, we obtain:

.. math:: \rho c_{p}\frac{DT}{Dt}+\sum_{k=1}^{K}\hat{h}_{k}\dot{\omega}_{k}-
          \sum_{k=1}^{K}h_{k}\nabla\cdot\mathbf{j}_{k}=
          -\nabla\cdot\left(-\lambda\nabla T+\sum_{k=1}^{K}h_{k}\mathbf{j}_{k}\right)

By expanding the term `\nabla\cdot\left(\sum_{k=1}^{K}h_{k}j_{k}\right)` and
performing appropriate cancellations, the energy equation simplifies slightly
to:

.. math:: \rho c_{p}\frac{DT}{Dt}+\sum_{k=1}^{K}\hat{h}_{k}\dot{\omega}_{k}+
          \sum_{k=1}^{K}\mathbf{j}_{k}\nabla\cdot h_{k}=
          \nabla\cdot\left(\lambda\nabla T\right)

Making the substitutions for the gradient, substantial derivative and
divergence, and setting the `z`-derivatives to zero, the one-dimensional form
of the energy equation is:

.. math:: \rho\frac{\partial T}{\partial t}+
          \rho v\frac{\partial T}{\partial r}+
          \frac{1}{c_{p}}\left(\sum_{k=1}^{K}\hat{h}_{k}\dot{\omega}_{k}+
          \sum_{k=1}^{K}j_{k}c_{p,k}\frac{\partial T}{\partial r}\right)=
          \frac{1}{c_{p}}\frac{1}{r^{\alpha}}\frac{\partial}{\partial r}\left[r^{\alpha}\lambda\frac{\partial T}{\partial r}\right]
