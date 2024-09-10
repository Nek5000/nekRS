.. _theory:

Theory
======

This page introduces the governing equations and the numerical discretization
at a high level. nekRS includes models for incompressible flow, a partially
compressible low-Mach formulation, the Stokes equations, and the :math:`k`-:math:`\tau`
:term:`RANS` equations.

.. _ins_model:

Incompressible Flow Model
-------------------------

The governing equations for conservation of mass, momentum, and energy for
an incompressible fluid are

.. math::

  \nabla\cdot\mathbf u=0

.. math::

  \rho\left(\frac{\partial\mathbf u}{\partial t}+\mathbf u\cdot\nabla\mathbf u\right)=-\nabla P+\nabla\cdot\tau+\rho\ \mathbf f

.. math::

  \rho C_p\left(\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T\right)=\nabla\cdot\left(k\nabla T\right)+\dot{q}

In these equations, :math:`t` is time,
:math:`\mathbf u` is the velocity vector, :math:`\rho` is the density, :math:`P` is the pressure,
:math:`\tau` is the viscous stress tensor, :math:`\mathbf f` is a force vector, :math:`C_p` is the
isobaric specific heat capacity, :math:`T` is the temperature, :math:`k` is the thermal conductivity,
and :math:`\dot{q}` is a volumetric heat source. If the viscosity is constant, the viscous stress tensor
can be contracted to give

.. math::

  \nabla\cdot\tau=\nabla^2\mathbf u

This is referred to as the "no-stress" formulation. In the general case for non-constant viscosity,
the viscous stress tensor is given by the Navier-Stokes closure as

.. math::

  \nabla\cdot\tau=\nabla\cdot\left\lbrack\mu\left(\nabla\mathbf u+\nabla\mathbf u^T\right)\right\rbrack

.. _nondimensional_eqs:

Non-Dimensional Formulation
"""""""""""""""""""""""""""

It is often advantageous to solve these equations in non-dimensional form. Here, we introduce
the non-dimensional form in as general a manner as possible, assuming the use of variable
material properties for density, viscosity, specific heat capacity, and thermal conductivity
that are functions of temperature, :math:`\rho=\rho(T)`, :math:`\mu=\mu(T)`,
:math:`C_p=C_p(T)`, and :math:`k=k(T)`. For simplicity, the functional notation is
omitted throughout.

Introduce
the non-dimensional variables :math:`\mathbf x^\dagger=\frac{\mathbf x}{L}`,
:math:`\mathbf u^\dagger=\frac{\mathbf u}{U}`, :math:`t^\dagger=\frac{tU}{L}`,
and :math:`\mathbf f^\dagger=\frac{\mathbf f L}{U^2}`. For the material properties,
we non-dimensionalize based on some reference temperature :math:`T_0`, such that
:math:`\rho^\dagger=\frac{\rho}{\rho_0}`, :math:`\mu^\dagger=\frac{\mu}{\mu_0}`,
:math:`k^\dagger=\frac{k}{k_0}`, and :math:`C_p^\dagger=\frac{C_p}{C_{p,0}}`. Here,
a subscript of :math:`0` is shorthand notation that indicates that the property
is evaluated at :math:`T_0`, such that :math:`k_0\equiv k(T_0)`.
Finally, for convection-dominated flows,
the pressure is non-dimensionalized in terms of the dynamic pressure as
:math:`P^\dagger=\frac{P}{\rho_0 U^2}`.

Inserting these non-dimensional variables
into the conservation of mass and momentum equations gives

.. math::

  \frac{\partial u_i^\dagger}{\partial x_i^\dagger}=0

.. math::

  \rho^\dagger\left(\frac{\partial u_i^\dagger}{\partial t^\dagger}+u_j^\dagger\frac{\partial u_i^\dagger}{\partial x_j^\dagger}\right)=-\frac{\partial P^\dagger}{\partial x_i^\dagger}+\frac{1}{Re}\frac{\partial \tau_{ij}^\dagger}{\partial x_j^\dagger}+\rho^\dagger f_i^\dagger

In these equations, the :math:`\nabla` are expanded to explicitly show that all derivatives
are taken with respect to the nondimensional space variable :math:`\mathbf x^\dagger`. :math:`Re`
is the Reynolds number

.. math::

  Re\equiv\frac{\rho_0 UL}{\mu_0}

To non-dimensionalize the energy conservation equation, use the previous non-dimensional
variables in addition to a non-dimensional temperature, :math:`T^\dagger=\frac{T-T_0}{\Delta T}`,
where :math:`\Delta T` is a reference temperature rise relative to a baseline temperature
:math:`T_0`. The heat source is non-dimensionalized as :math:`\dot{q}^\dagger=\frac{\dot{q}}{\rho_0 C_{p,0} U\Delta T/L}`,
which arises naturally from the simple formulation of bulk energy conservation of
:math:`Q=\dot{m}C_p\Delta T`, where :math:`Q` is a heat source (units of Watts) and
:math:`\dot{m}` is a mass flowrate.

Inserting these non-dimensional variables into the energy conservation equation gives

.. math::

  \rho^\dagger C_p^\dagger\left(\frac{\partial T^\dagger}{\partial t^\dagger}+u_i^\dagger\frac{\partial T^\dagger}{\partial x_i^\dagger}\right)=\frac{1}{Pe}\frac{\partial}{\partial x_i^\dagger}\left(k^\dagger\frac{\partial T^\dagger}{\partial x_i^\dagger}\right)+\dot{q}^\dagger

where :math:`Pe` is the Peclet number,

.. math::

  Pe\equiv\frac{LU}{\alpha}

and :math:`\alpha` is the thermal diffusivity,

.. math::

  \alpha\equiv\frac{k_0}{\rho_0 C_{p,0}}

Low-Mach Partially-Compressible Model
-------------------------------------

Stokes Equations
----------------

.. _rans_models:

RANS Models
-----------

The :term:`RANS` equations are derived from the conservation of mass, momentum, and energy
equations by expressing each term in the equation as the sum of a mean and a fluctuation.
Because nekRS is based on the incompressible flow model, all such averages (even for the
energy equation) are based on the notion of *Reynolds averaging*, where
each field :math:`f` is expressed as the sum of a time mean :math:`\overline{f}` and a time fluctuation,
:math:`f'`,

.. math::

  f(\mathbf x, t)=\overline{f}(\mathbf x)+f'(\mathbf x,t)

where the time averaged is defined as

.. math::

  \overline{f}=\lim_{S\rightarrow\infty}\frac{1}{S}\int_{t}^{t+S}f(\mathbf x,t)dt

For compressible flows in which energy conservation affects density, the :term:`RANS`
equations are instead derived with *Favre averaging*, where each field :math:`f_i`
is expressed as the sum of a density-weighted time average :math:`\tilde{f}_i`
and a fluctuation :math:`f_i^{''}`. It is therefore an important distinction here that
we only consider *Reynolds averaging*, which leads to a simpler formulation of
the :term:`RANS` energy conservation equation that the compressible flow case.

Inserting the above "Reynolds decomposition" for :math:`\mathbf u`, :math:`P`, and :math:`T`
into the governing equations and averaging in time
then gives the :term:`RANS` equations. For
the incompressible flow equations in :ref:`ins_model`,
the :term:`RANS` mass, momentum, and energy equations are

.. math::

  \frac{\partial\overline{u_i}}{\partial x_i} = 0

.. math::

  \rho\left(\frac{\partial\overline{u_i}}{\partial t}+\overline{u_j}\frac{\partial\overline{u_i}}{\partial x_j}+\frac{\partial}{\partial x_j}\overline{u_i'u_j'}\right)=-\frac{\partial \overline{P}}{\partial x_i}+\frac{\partial}{\partial x_j}\left(2\mu \overline{S_{ij}}\right)+\rho\overline{\mathbf f}

.. math::

  \rho C_p\left(\frac{\partial\overline{T}}{\partial t}+\overline{u_i}\frac{\partial\overline{T}}{\partial x_i}+\frac{\partial\overline{u_i'T'}}{\partial x_i}\right)=\frac{\partial}{\partial x_i}\left(k\frac{\partial\overline{T}}{\partial x_i}\right)+\overline{\dot{q}}

where :math:`\overline{S_{ij}}` is the mean strain rate tensor,

.. math::

  \overline{S_{ij}}=\frac{1}{2}\left(\frac{\partial \overline{u_i}}{\partial x_j}+\frac{\partial\overline{u_j}}{\partial x_i}\right)

The mass, momentum, and energy conservation equations have the same form as the instantaneous flow
equations in :ref:`ins_model` except for the addition of another
stress tensor to the momentum equation - :math:`\rho \overline{u_i'u_j'}`, and the addition
of another heat flux vector to the energy equation - :math:`\rho C_p\overline{u_i'T'}`.
The stress term in the momentum equation
is referred to as the Reynolds stress tensor;
:math:`\rho\ \partial(\overline{u_i'u_j'})/\partial x_j` represents the time-averaged rate
of momentum transfer due to turbulence. The heat flux term in the energy equation is
referred to as the turbulent heat flux; :math:`\rho C_p\partial\overline{u_i'T'}/\partial x_i`
represents the time-averaged rate of energy addition due to turbulence.
The objective of :term:`RANS` models is to provide
closures for the Reynolds stress tensor and turbulent heat flux vector in terms of the mean
properties such that the time-averaged
equations can be solved for the mean flow.

Boussinesq Approximation
""""""""""""""""""""""""

The :term:`RANS` models in nekRS are based on the Boussinesq eddy viscosity approximation,
which assumes that the momentum flux that induces the Reynolds stresses shares the same
functional form as the momentum flux that induces the molecular stresses. In other words,
the Navier-Stokes closure that was used to relate the deviatoric stress tensor
:math:`\tau_{ij}` to the strain rate tensor,

.. math::

  \tau_{ij}=\mu\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right)-\underbrace{\frac{2}{3}\mu\frac{\partial u_i}{\partial x_i}\delta_{ij}}_\text{$=\ 0$ if incompressible}

is assumed applicable to the Reynolds stress tensor, but with instantaneous velocities replaced by
mean velocities and the molecular viscosity replaced by the turbulent eddy viscosity
:math:`\mu_T`,

.. math::

  \rho\overline{u_i'u_j'}=\mu_T\left(\frac{\partial\overline{u_i}}{\partial x_j}+\frac{\partial\overline{u_j}}{\partial x_i}\right)-\underbrace{\frac{2}{3}\mu\frac{\partial \overline{u_i}}{\partial x_i}\delta_{ij}}_\text{$=\ 0$ if incompressible}-\frac{2}{3}\rho k\delta_{ij}

Here, :math:`k` is the turbulent kinetic energy,

.. math::

  k\equiv\frac{1}{2}\left(\overline{u_1'u_1'}+\overline{u_2'u_2'}+\overline{u_3'u_3'}\right)

The final term in the Boussinesq approximation for the Reynolds stress tensor simply ensures that
the trace of the Reynolds stress tensor equals :math:`2k`, because otherwise, for incompressible flows,
the trace of the Reynolds stress tensor would be zero. Inserting the Boussinesq eddy viscosity
model for the Reynolds stress tensor into the incompressible flow mean momentum equation then gives

.. math::

  \rho\left(\frac{\partial\overline{u_i}}{\partial t}+\overline{u_j}\frac{\partial\overline{u_i}}{\partial x_j}\right)=-\frac{\partial \overline{P}}{\partial x_i}+\frac{\partial}{\partial x_j}\left\lbrack 2\left(\mu+\mu_T\right) \overline{S_{ij}}-\frac{2}{3}\rho k\delta_{ij}\right\rbrack+\rho\overline{\mathbf f}

In nekRS, as well as many other :term:`RANS` codes, it is commonplace to combine
the gradient of :math:`P+\frac{2}{3}\rho k` terms together into a single reduced pressure,

.. math::

  P_r\equiv\overline{P}+\frac{2}{3}\rho k

such that the term proportional to :math:`\rho k` can be bundled into a single pressure gradient kernel,

.. math::

  -\frac{\partial\overline{P}}{\partial x_i}-\frac{\partial}{\partial x_j}\left(-\frac{2}{3}\rho k\delta_{ij}\right)\rightarrow-\frac{\partial P_r}{\partial x_i}

.. warning::

  The pressure solution, available on the ``nrs->P`` object and written to output under
  the name "pressure," represents this *reduced* pressure. To obtain :math:`\overline{P}`,
  you should subtract :math:`\frac{2}{3}\rho k` from ``nrs->P``.

Turbulent Prandtl Number
""""""""""""""""""""""""

Closure for the turbulent heat flux is typically motivated from considerations
of the analogy between momentum and energy transfer; while the Boussinesq approximation
was used to introduce a relationship between the Reynolds stress tensor :math:`\rho\overline{u_i'u_j'}`
in terms of the mean strain rate, the turbulent heat flux is assumed proportional to
the mean temperature gradient via a gradient diffusion approximation,

.. math::

  \rho C_p\overline{u_i'T'}=k_T\frac{\partial \overline{T}}{\partial x_i}

where :math:`k_T` is the turbulent conductivity. :math:`k_T` is related
to :math:`\mu_T`, the turbulent momentum diffusivity, by the turbulent Prandtl number
:math:`Pr_T`,

.. math::

  Pr_T\equiv\frac{\nu_T}{\alpha_T}

where :math:`\nu_T\equiv\mu_T/rho` and :math:`\alpha_T` is the turbulent thermal diffusivity,

.. math::

  \alpha_T\equiv\frac{k_T}{\rho C_p}

Inserting this gradient diffusion approximation into the incompressible flow
mean energy equation then gives

.. math::

  \rho C_p\left(\frac{\partial\overline{T}}{\partial t}+\overline{u_i}\frac{\partial\overline{T}}{\partial x_i}\right)=\frac{\partial}{\partial x_i}\left\lbrack\left(k+\frac{\mu_T}{Pr_T}C_p\right)\frac{\partial\overline{T}}{\partial x_i}\right\rbrack+\overline{\dot{q}}

.. _ktau:

The :math:`k`-:math:`\tau` Model
""""""""""""""""""""""""""""""""

nekRS uses the :math:`k`-:math:`\tau` turbulence model to close the mean flow equations [Kalitzin]_.
Because the :math:`k`-:math:`\epsilon`, :math:`k`-:math:`\omega`, and :math:`k`-:math:`\omega`
:term:`SST` models tend to dominate the :term:`RANS` space, extra discussion is devoted here
to motivating the use of this particular model. Because :math:`Pr_T` is typically taken as a
constant, often 0.90 [Wilcox]_, the objective of incompressible flow :term:`RANS` models is to compute
the eddy viscosity and :math:`k` needed to close the
mean momentum and energy equations.

.. note::

  Take care not to confuse the inverse of the specific dissipation rate, :math:`\tau`, with
  the deviatoric molecular stress tensor, which is also represented here as :math:`\tau` due to convention.

The :math:`k`-:math:`\tau` model is a modification of the standard :math:`k`-:math:`\omega`
turbulence model that bases the second transport equation on the *inverse* of the specific
dissipation rate :math:`\omega`,

.. math::

  \tau\equiv\frac{1}{\omega}

rather than the on :math:`\omega`.
The :math:`k`-:math:`\tau` model attempts to retain two important
features of the :math:`k`-:math:`\omega` model -

  1. Good predictions for flows with adverse pressure gradients and separation, and
  2. Reasonable prediction of boundary layers and near-wall behavior without wall functions
     or special low-:math:`Re_t` treatments.

These two aspects contribute to better predictions of complex flows with reduced
numerical complexity associated with wall functions or
damping functions that can cause stiff behavior [Kok]_ and inaccurate flow predictions. By introducing the
definition of :math:`\tau\equiv 1/\omega`, the :math:`k`-:math:`\tau` model attemps to improve upon
the :math:`k`-:math:`\omega` model in two main ways -

  1. Simplify wall boundary conditions for the second transport equation, and
  2. Bound the source terms in the second transport equation in near-wall regions.

As :math:`y\rightarrow 0`, :math:`\omega\rightarrow y^{-2}`, while
:math:`k\rightarrow 0` [Kok]_. Therefore, while :math:`\omega` is infinite
at walls, :math:`\tau` is zero. Traditionally, this singular behavior in :math:`\omega`
was treated by applying "rough wall" boundary conditions to :math:`\omega`
with the wall roughness set to a "small enough" value to simulate a hydraulically
smooth wall [Kok]_. However, this ad hoc approach retains a strong dependence
on the near-wall mesh resolution, often requiring prohibitively fine elements to
accurately predict boundary layer properties [Kalitzin]_. And,
such an approach retains near-singular behavior in the first and second derivatives of
:math:`\omega`. Applying a zero boundary condition to :math:`\tau`
on solid walls is comparatively trivial.

With regards to the second point, the :math:`\omega` transport equation contains a source term
propotional to :math:`\omega^2`; because :math:`\omega\rightarrow y^{-2}` as :math:`y\rightarrow 0`,
this source term displays singular behavior as :math:`y\rightarrow 0`. Singular behavior
of the source terms can result in large numerical errors and stiffness that negatively
affects the convergence of the computational solution. Conversely, all source terms in
the :math:`\tau` equation are bounded near walls [Kok]_.

With this motivation, the :math:`k` and :math:`\tau` equations are described next.
A slightly lengthier description is provided for each in order to give greater context
to the genesis of this model.

**The** :math:`k` **Equation**

The :math:`k`
equation is a *model* version of the *true* :math:`k` equation. The *true* :math:`k`
equation is derived by taking the trace of the Reynolds stress equation, a process that
is itself motivated by recognition that the trace of the Reynolds stress tensor is equal
to :math:`2k`,

.. math::

  \overline{u_i'u_i'}=2k

The *true* :math:`k` equation contains terms that depend on the mean flow velocity,
the turbulent kinetic energy, and the dissipation, in addition to more exotic terms such as
:math:`\overline{u_i'u_i'u_j'}` and :math:`\overline{P'u_j'}`. These additional fluctuating
terms do not bring the *true* :math:`k` equation any closer to a tractable solution,
so Prandtl introduced a :math:`\partial k/\partial x_j`
gradient diffusion approximation for the turbulent transport and
pressure diffusion terms (:math:`\frac{1}{2}\rho\overline{u_i'u_i'u_j'}+\overline{P'u_j'}`)
with a diffusion coefficient of :math:`\mu_T/\sigma_k`, where :math:`\sigma_k`
is a constant [Wilcox]_. With this gradient diffusion model, the *true*
:math:`k` equation is simplified to a tractable *model* :math:`k` equation [Launder]_,

.. math::

  \frac{\partial(\rho k )}{\partial t}+\nabla\cdot\left(\rho k\overline{\mathbf u}\right)=\nabla\cdot\left\lbrack\left(\mu+\frac{\mu_T}{\sigma_k}\right)\nabla k\right\rbrack+\mathscr{P}-\rho\epsilon

where :math:`\mathscr{P}` is the production of turbulent kinetic energy by velocity shear,

.. math::

  \mathscr{P}\equiv\rho\overline{u_i'u_j'}\frac{\partial\overline{u_i}}{\partial x_j}

and :math:`\epsilon` is the dissipation per unit mass,

.. math::

  \epsilon\equiv\nu\overline{\frac{\partial u_i'}{\partial x_j}\frac{\partial u_i'}{\partial x_j}}

The production term represents the rate at which energy is transferred from the mean
flow to the turbulent flow, while the dissipation term represents the rate at which
turbulent kinetic energy is converted to heat. Note that the only difference between this
*model* :math:`k` equation and the *true* :math:`k` equation is the introduction of the
gradient diffusion approximation for the turbulent transport and pressure diffusion terms.

The :math:`k` equation used in the
:math:`k`-:math:`\tau` model is then
obtained as a simple transformation of
the standard :math:`k` equation by the following
relationship [Kok]_,

.. math::

  \omega\equiv\frac{\epsilon}{\beta^*k}

where :math:`\beta^*` is a constant. Inserting :math:`\omega\beta^*k`
for :math:`\epsilon` in the dissipation term :math:`\rho\epsilon` gives
the :math:`k` equation used in nekRS,

.. math::

  \frac{\partial(\rho k )}{\partial t}+\nabla\cdot\left(\rho k\overline{\mathbf u}\right)=\nabla\cdot\left\lbrack\left(\mu+\frac{\mu_T}{\sigma_k}\right)\nabla k\right\rbrack+\mathscr{P}-\rho\beta^*\frac{k}{\tau}

**The** :math:`\tau` **Equation**

In two-equation :term:`RANS` turbulence modeling, the greatest source of uncertainty is
the proper choice of the second transport variable. While a *true* :math:`k` equation
is often used as the starting point for developing the *model* :math:`k` equation,
it is commonplace to start immediately from an ad hoc, "fabricated," model equation
for the second turbulence variable. Of course, "exact" equations can always be
derived for the second turbulence variable through various operations on the mean Navier-Stokes
equation or the Reynolds stress equation, but the exact equations for :math:`\epsilon`,
:math:`\omega`, or other turbulence quantities tend to be far more complex than
the exact equation for :math:`k` shown earlier.

In 1942, Kolmogorov was the first to
propose the :math:`k`-:math:`\omega` model [Wilcox]_. His formulation was
very heuristic - from the Boussinesq approximation, it is likely that
:math:`\nu_T\propto k`, which requires another variable with dimensions inverse time.
Based on the work of Kolmogorov and many subsequent researchers of the
:math:`k`-:math:`\omega` model, inserting :math:`\tau\equiv 1/\omega` into the
:math:`\omega` equation gives the :math:`\tau` equation - this approach is very similar
to that used to obtained the :math:`k` equation.
The :math:`\tau` equation used in nekRS is [Kok]_

.. math::

  \frac{\partial(\rho\tau)}{\partial t}+\nabla\cdot\left(\rho\tau\overline{\mathbf u}\right)=\nabla\cdot\left\lbrack\left(\mu+\frac{\mu_T}{\sigma_\tau}\right)\nabla \tau\right\rbrack-\alpha\frac{\tau}{k}\mathscr{P}+\rho\beta-2\frac{\mu}{\tau}\nabla\tau\cdot\nabla\tau

..
   TODO:
   The Kok version has
   mu_t/sigma_tau added to the viscosity on the last term.

where :math:`\sigma_\tau`, :math:`\alpha`, and :math:`\beta` are constants. The
last term on the right-hand side of the :math:`\tau` equation is in practice
implemented in the form

.. math::

  \frac{2}{\tau}\nabla\tau\cdot\nabla\tau\rightarrow 8\nabla\sqrt{\tau}\cdot\nabla\sqrt{\tau}

in order to reduce the discretization error associated with the computation
of gradients of a term that scales as :math:`y^2` as :math:`y\rightarrow 0` [Kok]_.

**The Eddy Viscosity**

The objective of :term:`RANS` models is to estimate the eddy viscosity :math:`\mu_T`
that appears in the Boussinesq approximation. The particular form for :math:`\mu_T`
can be understood here in terms of the standard :math:`k`-:math:`\epsilon`
model [Launder]_, for which :math:`\mu_T` is given as


.. math::

  \mu_T=C_\mu\rho\frac{k^2}{\epsilon}

where :math:`C_\mu` is a constant. Inserting :math:`\tau\equiv 1/\omega` and
:math:`\epsilon=\beta^*\omega k` gives

.. math::

  \mu_T=\rho k\tau

which presumes that :math:`C_\mu` and :math:`\beta^*` are really the same constant,
but with different notation developed separately by the :math:`k`-:math:`\epsilon`
researchers and the :math:`k`-:math:`\tau` researchers [Kok]_.

..
  TODO: is this the correct explanation for why there's no coefficient in the mu_t equation?

**Closure Coefficients and Other Details**

Table :ref:`RANS Coefficients <rans_coeffs>` shows the values for the various
constants used in nekRS's :math:`k`-:math:`\tau` model.

.. _rans_coeffs:

.. table:: RANS Coefficients

  ==================== =================== ======
  Coefficient          Value               Source
  ==================== =================== ======
  :math:`\sigma_k`     :math:`5/3`
  :math:`\sigma_\tau`  :math:`2.0`
  :math:`Pr_T`         user-selected       ---
  ==================== =================== ======

A limiter is applied to both :math:`k` and :math:`\tau` to prevent negative values
of either :math:`k` or :math:`\tau`,

.. math::

  k = \max{\left(k, 0.01|k|\right)}

.. math::

  \tau = \max{\left(\tau, 0.01|\tau|\right)}

.. warning::

  nekRS's :math:`k`-:math:`\tau` implementation currently requires that
  the laminar dynamic viscosity and the density are constant, because the setup
  routines can only accept constant values. See :ref:`RANS Plugin <rans_plugin>`
  for more information.

.. note::

  Even if the molecular viscosity is constant, you must set ``stressFormulation = true`` in
  the input file because the total viscosity (molecular plus turbulent) will not be constant.

**Boundary Conditions**

On walls, because the asymptotic behavior of :math:`\omega` is :math:`\omega\propto y^{-2}`
as :math:`y\rightarrow0`, and because the instantaneous velocity
:math:`u_i\equiv \overline{u_i}+u_i'` must be zero due to the no-slip condition, both
:math:`k` and :math:`\tau` should be set to zero on no-slip boundaries.

On turbulent inlets, however, both :math:`k` and :math:`\tau` are generally nonzero, and
estimates for :math:`k` and :math:`\tau` must be provided. Turbulent inflow conditions are
usually unknown unless the modeler is fortunate enough to have experimental data - therefore,
the specification of inlet conditions on :math:`k` and :math:`\tau` tends to be fairly ad hoc.
An inlet condition for :math:`k` can be estimated by prescribing the turbulent intensity at
the inlet. The turbulent intensity :math:`I` is defined as the root mean square
fluctuating velocity normalized by the magnitude of the mean velocity, or

.. math::

  I\equiv\frac{\sqrt{\frac{1}{3}\left(u_i'u_i'\right)}}{\sqrt{\overline{u_j}\ \overline{u_j}}}

which can equivalently be written in terms of the turbulent kinetic energy as

.. math::

  I\equiv\frac{\sqrt{\frac{2}{3}k}}{\sqrt{\overline{u_j}\ \overline{u_j}}}

Therefore, if an inlet turbulent intensity can be prescribed, the inlet
turbulent kinetic energy is

.. math::

  k\equiv\frac{3}{2}I^2\overline{u_j}\ \overline{u_j}

For instance, it is common to assume a uniform turbulent intensity over an inlet of
between 1 and 5% for pipe flows [Russo]_. Experiments have also quantified the scaling
of turbulent intensity with Reynolds number.
From simulations and experiments of both incompressible
and compressible flows, the turbulent intensities
on the axis of a smooth circular pipe are [Russo]_:

.. math::

  I_\text{axis}=\begin{cases}0.0853Re^{-0.0727} & \text{incompressible}\\
  0.0550Re^{-0.0407} & \text{compressible}\end{cases}

With a slightly different definition based on the turbulent intensity averaged over
the cross-sectional area of a circular pipe, the turbulent intensity instead
scales as [Russo]_:

.. math::

  I_\text{area}=\begin{cases}0.140Re^{-0.0790} & \text{incompressible}\\
  0.227Re^{-0.1} & \text{compressible}\end{cases}

Similar correlations have also been developed for rough pipes [Basse]_. It is
also commonplace to apply conditions on :math:`k` in terms of the friction
velocity :math:`u_\tau`, defined as

.. math::

  u_\tau\equiv\sqrt{\frac{\tau_w}{\rho}}

where :math:`\tau_w` is the wall shear stress. The wall shear stress can then be
estimated using a friction factor correlation for :math:`f_D`, the Darcy friction factor,
which is defined as

.. math::

  \frac{\Delta P}{l}\equiv f_D\frac{1}{D}\frac{\rho \overline{u_j}\overline{u_j}}{2}

where :math:`\Delta P/l` is a pressure gradient and :math:`D` is a hydraulic diameter.
For circular pipes, the friction factor is related to :math:`\tau_w` as

.. math::

  f_D=8\frac{\tau_w}{\rho\overline{u_j}\overline{u_j}}

Typically, a simple friction factor correlation such as the Blasius model for pipes,
is selected

.. math::

  f_D=0.316Re^{-0.25}

Combining :math:`u_\tau`, :math:`f_D`, and one of the previous estimates for :math:`I_\text{area}`,
such as :math:`I_\text{area}=0.317Re^{-0.11}` [Basse]_ gives the following relationship between
:math:`I_\text{area}` and :math:`f_D`,

.. math::

  I_\text{area}=0.526f_D^{0.44}

Then, inserting the relationship between :math:`u_\tau` and :math:`f_D` gives,
after some manipulation,

.. math::

  k=2.56 u_\tau^2\left(\overline{u_j}\ \overline{u_j}\right)^{0.24}

While this algebraic exercise hasn't actually introduced any new information or closures
beyond what has already been discussed, it is information to see :math:`k` expressed
in terms of the friction velocity, because some nekRS inputs set :math:`k` on inlets
by first computing the :math:`u_\tau` from a friction correlation and then using an
expression similar to above.

In any case, one of these correlations for turbulent intensity,
or simply a fixed turbulent intensity of, say, 5%,
can be used to prescribe a uniform value of :math:`k` on an inlet. However, including
some spatial variation in :math:`k` on the inlet may reduce Gibbs phenomena if
the inlet turbulent intensity enforces the physical zero wall value. Spatial fits,
such as those developed for circular pipes [Russo]_, may improve the numerical stability
and accuracy of your simulation.

Finally, on outlets, "free-stream" boundary conditions are usually applied to
:math:`k` and :math:`\tau`.

**Non-Dimensional Formulation**

Now that the :math:`k`-:math:`\tau` model has been presented in its full form,
the non-dimensional formulation is provided. Because nekRS's :math:`k`-:math:`\tau` model
is currently limited to constant densities and laminar viscosities, the non-dimensional
formulation is slightly simpler than the more general case shown in
:ref:`Non-Dimensional Formulation <nondimensional_eqs>` that derived the non-dimensional
form of the instantaneous Navier-Stokes equations.

Introduce
the non-dimensional variables :math:`\mathbf x^\dagger=\frac{\mathbf x}{L}`,
:math:`\overline{\mathbf u}^\dagger=\frac{\overline{\mathbf u}}{U}`, :math:`t^\dagger=\frac{tU}{L}`,
and :math:`\overline{\mathbf f}^\dagger=\frac{\overline{\mathbf f} L}{U^2}`.
For convection-dominated flows,
the pressure is non-dimensionalized in terms of the dynamic pressure as
:math:`P^\dagger=\frac{P}{\rho U^2}`. Finally, the turbulent kinetic energy
and inverse dissipation rate are non-dimensionalized as
:math:`k^\dagger=\frac{k}{U^2}` and :math:`\tau^\dagger=\frac{\tau U}{L}`.

Inserting these non-dimensional variables
into the mean flow mass and momentum conservation equations gives

.. math::

  \frac{\partial \overline{u_i}^\dagger}{\partial x_i^\dagger}=0

.. math::

  \frac{\partial \overline{u_i}^\dagger}{\partial t^\dagger}+\overline{u_j}^\dagger\frac{\partial \overline{u_i}^\dagger}{\partial x_j^\dagger}=-\frac{\partial \overline{P}^\dagger}{\partial x_i^\dagger}+\frac{1}{Re}\left(1+\frac{\mu_T}{\mu}\right)\frac{\partial\overline{\tau_{ij}}^\dagger}{\partial x_j^\dagger}-\frac{\partial}{\partial x_j^\dagger}\left(\frac{2}{3}k^\dagger\delta_{ij}\right)+f_i^\dagger

To non-dimensionalize the energy conservation equation, use the previous non-dimensional
variables in addition to a non-dimensional temperature, :math:`T^\dagger=\frac{T-T_0}{\Delta T}`,
where :math:`\Delta T` is a reference temperature rise relative to a baseline temperature
:math:`T_0`. The heat source is non-dimensionalized as :math:`\dot{q}^\dagger=\frac{\dot{q}}{\rho C_{p} U\Delta T/L}`,
which arises naturally from the simple formulation of bulk energy conservation of
:math:`Q=\dot{m}C_p\Delta T`, where :math:`Q` is a heat source (units of Watts) and
:math:`\dot{m}` is a mass flowrate.
Inserting these non-dimensional variables into the energy conservation equation gives

.. math::

  \frac{\partial \overline{T}^\dagger}{\partial t^\dagger}+\overline{u_i}^\dagger\frac{\partial \overline{T}^\dagger}{\partial x_i^\dagger}=\frac{1}{Pe}\left(1+\frac{\mu_T/Pr_T}{k}C_p\right)\frac{\partial}{\partial x_i^\dagger}\frac{\partial \overline{T}^\dagger}{\partial x_i^\dagger}+\dot{q}^\dagger

To non-dimensionalize the :math:`k` and :math:`\tau` equations,
define :math:`\mathscr{P}^\dagger=\frac{\mathscr{P}}{\rho U^3/L}`
and :math:`\epsilon^\dagger=\frac{\epsilon}{U^3/L}`.
With previous non-dimensional variables already defined, the non-dimensional :math:`k`
and :math:`\tau` equations become

.. math::

  \frac{\partial k^\dagger}{\partial t^\dagger}+\overline{u_i}^\dagger\frac{\partial k^\dagger}{\partial x_i^\dagger}=\frac{1}{Re}\left(1+\frac{\mu_T/\sigma_k}{\mu}\right)\frac{\partial}{\partial x_i^\dagger}\frac{\partial k^\dagger}{\partial x_i^\dagger}+\mathscr{P}^\dagger-\beta^*\frac{k^\dagger}{\tau^\dagger}

.. math::

  \frac{\partial\tau^\dagger}{\partial t^\dagger}+\overline{u_i}^\dagger\frac{\partial\tau^\dagger}{\partial x_i^\dagger}=\frac{1}{Re}\left(1+\frac{\mu_T/\sigma_\tau}{\mu}\right)\frac{\partial}{\partial x_i^\dagger}\frac{\partial\tau^\dagger}{\partial x_i^\dagger}-\alpha\frac{\tau^\dagger}{k^\dagger}\mathscr{P}^\dagger+\beta-\frac{2}{Re}\frac{1}{\tau^\dagger}\frac{\partial\tau^\dagger}{\partial x_i^\dagger}\frac{\partial\tau^\dagger}{\partial x_i^\dagger}

Finally, the eddy viscosity is computed as

.. math::

  \mu_T=\rho k\tau

which after inserting the non-dimensional variables becomes

.. math::

  \mu_T=\rho ULk^\dagger\tau^\dagger

such that the non-dimensional eddy viscosity can be written as :math:`\mu_T^\dagger=k^\dagger\tau^\dagger`.
Therefore, the overall diffusion coefficients in the mean momentum equation,
mean energy equation, :math:`k` equation, and :math:`\tau` equation are, respectively

.. math::

  \frac{1}{Re}\left(1+\frac{\mu_T}{\mu}\right)\rightarrow\frac{1}{Re}+\mu_T^\dagger

.. math::

  \frac{1}{Pe}\left(1+\frac{\mu_T/Pr_T}{k}C_p\right)\rightarrow\frac{1}{Pe}+\frac{\mu_T^\dagger}{Pr_T}

.. math::

  \frac{1}{Re}\left(1+\frac{\mu_T/\sigma_k}{\mu}\right)\rightarrow\frac{1}{Re}+\frac{\mu_T^\dagger}{\sigma_k}

.. math::

  \frac{1}{Re}\left(1+\frac{\mu_T/\sigma_\tau}{\mu}\right)\rightarrow\frac{1}{Re}+\frac{\mu_T^\dagger}{\sigma_\tau}