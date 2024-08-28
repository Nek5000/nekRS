.. _detailed:

Detailed Usage
===================

This page describes how to perform a wide variety of user interactions with nekRS
for setting boundary conditions, converting between mesh formats, defining and
running device kernels, writing output files, and much more. Please first consult
the :ref:`Input File Syntax <input>` page for an overview of the purpose of each
nekRS input file to provide context on where the following instructions fit into
the overall code structure. Throughout this section, variables and data structures
in the nekRS source code are referenced - a list defining these variables and structures
is available on the :ref:`Commonly Used Variables <commonly_used_variables>` page
for reference.

.. _defining_variables_for_device:

Defining Variables to Access in Device Kernels
----------------------------------------------

The customization of a nekRS problem to a specific case is one with both the host-side
user functions in the ``.udf`` file, as well as device-side user functions in the ``.oudf``
file. For convenience purposes, nekRS supports setting non-pointer-type variables in the
``.udf`` file that are accessible in the device kernels in the ``.oudf`` file. This section
shows an example of this usage.

Suppose that a device kernel requires a parameter representing a pressure gradient, which
is then used to determine a forcing kernel. One option would be to pass that pressure gradient
to the device kernel through its function parameters. The kernel in the ``.oudf`` file
would look something like the following.

.. code-block::

    @kernel void myForcingKernel(const dfloat dp_dx, /* more parameters */)
    {
      double foo = 2 * dp_dx;

      // do something
    }

Alternatively, we can define a variable, ``p_dp_dx``, that we set from the ``.udf`` file.
While this variable propagation can be done in any of the user-defined functions that
has ``nrs`` as an input parameter, for consistency purposes we will use the ``UDF_LoadKernels``
function for this purpose.

.. note:

  The convention is to precede any of these host-side kernel variable
  definitions with a ``p_``.

To set ``p_dp_dx`` to 5.5 from the ``.udf`` file, write to the ``kernelInfo`` object
on the ``nrs`` object. The ``defines/<p_name>`` syntax indicates that a variable on
the device is being declared with a name ``p_name`` that will be accessible simply as
``p_name`` in the device kernels.

.. code-block::

   void UDF_LoadKernels(nrs_t * nrs)
   {
     occa::properties & kernelInfo = *nrs->kernelInfo;

     kernelInfo["defines/p_dp_dx"] = 5.5;

     // other stuff related to loading the kernels
   }

Then, the kernel would be simplified to the following. You will note that nothing needs
to be passed through the kernel function arguments - ``p_dp_dx`` is simply available as
if it were a local variable to the function.

.. code-block:: cpp

   @kernel void myForcingKernel(/* more parameters */)
   {
     double foo = 2 * p_dp_dx;

     // do something
   }

If you grep for ``kernelInfo["defines`` in the nekRS source code, you will see that
this variable propagation features is also used extensively throughout a normal problem
setup. For instance, the number of velocity fields to solve for is propagated to the device
in the ``nrsSetup`` function.

.. code-block:: cpp

   nrs_t* nrsSetup(MPI_Comm comm, occa::device device, setupAide &options, int buildOnly)
   {
     // ...

     kernelInfo["defines/p_NVfields"] = nrs->NVfields;

     // ...
   }

Again, the convention is to precede all such propagated variables with the ``p_`` prefix.
No list of all such variables propagated automatically within a nekRS simulation is
maintained, so always check if the information you'd like to propagate is perhaps
already automatically propagated.

.. _boundary_conditions:

Setting Boundary Conditions with Device Kernels
-----------------------------------------------

Because all nekRS solves are performed on the device, boundary conditions on the
solution (which may change from time step to time step and be arbitrary functions
of the solution itself) are also applied on the device. The types of boundary conditions
on each solution field are specified in the ``.par`` file with the ``boundaryTypeMap``
key. 

.. _custom_properties:

Setting Custom Properties
-------------------------

Custom material properties can be set for the flow and passive scalar equations
by assigning the ``udf.properties`` function pointer to a function with a signature
that takes the ``nrs`` pointer to the nekRS solution object, the simulation time
``time``, the velocity solution on the device ``o_U``, the passive scalar solution
on the device ``o_S``, the flow material properties on the device ``o_UProp``,
and the passive scalar material properties on the device ``o_SProp``.

This section provides an example of setting :math:`\mu` and :math:`\rho` for the flow
equations and :math:`k` and :math:`\rho C_p` for two passive scalars. Suppose our problem
contains velocity, pressure, temperature, and two passive scalars. The ``[VELOCITY]``,
``[PRESSURE]``, ``[TEMPERATURE]``, ``[SCALAR01]``, and ``[SCALAR02]`` sections of the
``.par`` file would be as follows. Because we will be setting custom properties for
the pressure, velocity, and first two passive scalars (temperature and ``SCALAR01``),
we can let nekRS assign the default values of unity to all properties for those
governing equations until we override them in our custom property function. We still
need to define the material properties for ``SCALAR02``, however, because we will not
be overriding those properties in our function.

.. code-block::

  [PRESSURE]
  residualTol = 1e-6

  [VELOCITY]
  boundaryTypeMap = v, O, W
  residualTol = 1e-8

  [TEMPERATURE]
  boundaryTypeMap = t, O, I
  residualTol = 1e-8

  [SCALAR01]
  boundaryTypeMap = t, O, I
  residualTol = 1e-8

  [SCALAR02]
  boundaryTypeMap = t, O, t
  residualTol = 1e-7
  conductivity = 3.5
  rhoCp = 2e5

Also suppose that our problem contains conjugate heat transfer, such that some of
the mesh is fluid while some of the mesh is solid.

In ``UDF_Setup``, we next need to assign an address to the ``udf.properties`` function
pointer to a function with the correct signature where we eventually assign our custom
properties. Our ``UDF_Setup`` function would be as follows.

.. code-block:: cpp

   void UDF_Setup(nrs_t* nrs)
   {
     udf.properties = &material_props;
   }

Here, ``material_props`` is our name for a function in the ``.udf`` file that sets the
material properties. Its name is arbitrary, but it must have the following signature.

.. code-block:: cpp

   void material_props(nrs_t* nrs, dfloat time, occa::memory o_U, occa::memory o_S,
     occa::memory o_UProp, occa::memory o_SProp)
   {
     // set the material properties
   }

This function is called *after* the solve has been performed on each time step, so the
material properties are lagged by one time step with respect to the simulation.

.. note::

  You must place the ``material_props`` function *before* ``UDF_Setup`` (and before any other
  function that uses ``material_props``) in the ``.udf`` file in order for the just-in-time
  compilation to succeed.

Suppose we would like to set :math:`\rho=1000.0` and :math:`\mu=2.1e-5 e^{-\phi_0/500}(1+z)` for
the flow equations; because only the fluid domain has flow, we do not need to set
these properties on the solid part of the domain. For the first passive scalar
:math:`\phi_0`, we would
like to set :math:`(\rho C_p)_f=2e3(1000+PV_x)` and :math:`k_f=2.5` in the fluid
domain, and :math:`(\rho C_p)_s=2e3(1000+PV_x)` and :math:`k_s=3.5` in the solid domain.
Here, :math:`P` is the thermodynamic pressure and :math:`V_x` is the :math:`x`-component velocity.
For the second passive scalar :math:`\phi_1`, we would like to set
:math:`\rho C_p=0` and :math:`k=5+\phi_0` in both the fluid and solid domains.
Our material property function would be as follows. Note that these boundary conditions
are selected just to be comprehensive and show all possible options for setting
constant and non-constant properties with dependencies on properties - they do not
necessarily represent any realistic physical case.

.. code-block:: cpp

   // declare all the kernels we will be writing
   static occa::kernel viscosityKernel;
   static occa::kernel constantFillKernel;
   static occa::kernel heatCapacityKernel;
   static occa::kernel conductivityKernel;

   void material_props(nrs_t* nrs, dfloat time, occa::memory o_U, occa::memory o_S,
     occa::memory o_UProp, occa::memory o_SProp)
   {
     mesh_t* mesh = nrs->mesh;

     // viscosity and density for the flow equations
     const occa::memory o_mue = o_UProp.slice(0 * nrs->fieldOffset * sizeof(dfloat));
     const occa::memory first_scalar = o_S.slice(0 * cds->fieldOffset * sizeof(dfloat));
     viscosityKernel(mesh->Nelements, first_scalar, mesh->o_z, o_mue);

     const occa::memory o_rho = o_UProp.slice(1 * nrs->fieldOffset * sizeof(dfloat));    
     constantFillKernel(nrs->mesh->Nelements, 1000.0, 0.0 /* dummy */, nrs->o_elementInfo, o_rho);

     // conductivity and rhoCp for the first passive scalar
     int scalar_number = 0;
     const occa::memory o_con = o_SProp.slice((0 + 2 * scalar_number) *
       cds->fieldOffset * sizeof(dfloat));
     constantFillKernel(mesh->Nelements, 2.5, 3.5, nrs->o_elementInfo, o_con);

     const occa::memory o_rhocp = o_SProp.slice((1 + 2 * scalar_number) *
       cds->fieldOffset * sizeof(dfloat));
     heatCapacityKernel(mesh->Nelements, o_U, nrs->o_P, o_rhocp);

     // conductivity and rhoCp for the second passive scalar
     scalar_number = 1;
     const occa::memory o_con_2 = o_SProp.slice((0 + 2 * scalar_number) *
       cds->fieldOffset * sizeof(dfloat));
     conductivityKernel(mesh->Nelements, first_scalar, o_con_2);

     const occa::memory o_rhocp_2 = o_SProp.slice((1 + 2 * scalar_number) *
       cds->fieldOffset * sizeof(dfloat));
     constantFillKernel(mesh->Nelements, 0.0, 0.0, nrs->o_elementInfo, o_rhocp_2);
   }

The ``o_UProp`` and ``o_SProp`` arrays hold all material
property information for the flow equations and passive scalar equations, respectively.
In this function, you see six "slice" operations performed on ``o_UProp`` and ``o_SProp``
in order to access the two individual properties (diffusive constant and time derivative constant)
for the three equations (momentum, scalar 0, and scalar 1). The diffusive constant
(:math:`\mu` for the momentum equations and :math:`k` for the passive scalar equations)
is always listed first in these arrays, while the coefficient on the time derivative
(:math:`\rho C_p` for the momentum equations and :math:`\rho C_p` for the passive scalar
equations) is always listed second in these arrays.

To further elaborate, :math:`\mu` and :math:`\rho` are accessed as slices on ``o_UProp``.
Because viscosity is listed before density, the offset in the ``o_UProp`` array to get
the viscosity is zero, while the offset to get the density is ``nrs->fieldOffset``.
:math:`k` and :math:`\rho C_p` are accessed as slices in ``o_SProp``. Because the
passive scalars are listed in order and the conductivity is listed first for each user,
the offset in the ``o_SProp`` array to get the conductivity for the first passive scalar
is zero, while the offset to get the heat capacity for the first passive scalar 
is ``cds->fieldOffset``. Finally, the offset in the ``o_SProp`` array to get the conductivity
for the second passive scalar is ``2 * cds->fieldOffset``, while the offset to get the
heat capacity for the second passive scalar is ``3 * cds->fieldOffset``.

The ``viscosityKernel``, ``constantFillKernel``, ``heatCapacityKernel``,
and ``conductivityKernel`` functions are all user-defined device kernels. These
functions must be defined in the ``.oudf`` file, and the names are arbitrary. For each
of these kernels, we declare them at the top of the ``.udf`` file. In order to link
against our device kernels, we must instruct nekRS to use its just-in-time compilation
to build those kernels. We do this in ``UDF_LoadKernels`` by calling the
``udfBuildKernel`` function for each kernel. The second argument to the ``udfBuildKernel``
function is the name of the kernel, which appears as the actual function name of
the desired kernel in the ``.oudf`` file.

.. code-block:: cpp

  void UDF_LoadKernels(nrs_t* nrs)
  {
    viscosityKernel = udfBuildKernel(nrs, "viscosity");
    constantFillKernel = udfBuildKernel(nrs, "constantFill");
    heatCapacityKernel = udfBuildKernel(nrs, "heatCapacity");
    conductivityKernel = udfBuildKernel(nrs, "conductivity");
  }

In order to write these device kernels, you will need some background in programming
with :term:`OCCA`. Please consult the `OCCA documentation <https://libocca.org/#/>`__
before proceeding [#f1]_.

First, let's look at the ``constantFill`` kernel. Here, we want to write a device kernel
that assigns a constant value to a material property. So that we can have a general
function, we will write this such that it can be used to set constant (but potentially
different) properties in the fluid and solid phases for conjugate heat transfer
applications.

.. note::
  
  Material properties for the flow equations (i.e. viscosity and density) do not
  *need* to be specified in the solid phase. If you define flow properties in solid
  regions, they are simply not used.

The ``constantFill`` kernel is defined in the ``.oudf`` file as follows [#f2]_. :term:`OCCA`
kernels operate on the device. As input parameters, they can take non-pointer objects
on the host (such as ``Nelements``, ``fluid_val``, and ``solid_val`` in this example),
as well as pointers to objects of type ``occa::memory``, or device-side memory. The
device-side objects are indicated with the ``@restrict`` tag. 

.. note::

  Device-side memory in nekRS is by convention preceded with a ``o_`` prefix in order
  to differentiate from the host-side objects. In the initialization of nekRS, most of
  the simulation data is copied over to the device. All calculations are done on the
  device. The device-side solution is then only copied back onto the host for the
  purpose of writing output files.

.. warning::

  Because nekRS by default only copies the device-side solution back to the host for
  the purpose of writing output files, if you touch any host-side objects in your
  user-defined functions, such as in ``UDF_ExecuteStep``, you must ensure
  that you only use the host-side objects after they have been copied from device back
  to the host. Otherwise, they would not be "up to date." You can ensure that the host-
  side objects reflect the real-time nekRS solution by either (a) only touching the
  host-side solution on output writing steps (which can be determined based on the
  ``nrs->isOutputStep`` variable), or (b) calling the appropriate routines in nekRS
  to force data to be copied from the device back to the host. For the latter option,
  please refer to the :ref:`Copying From Device to Host <copy_device_to_host>` section.

For this example, we
loop over all the elements. The ``eInfo`` parameter represents a mask, and takes a value
of zero for solid elements and a value of unity for fluid elements. Next, we loop over
all of the :term:`GLL` points on the element, or ``p_Np``. This variable is set within
nekRS to be the same as ``mesh->Np`` using the device variable feature described in
the :ref:`Defining Variables to Access in Device Kernels <defining_variables_for_device>`
section. This particular variable is always available, and you do not need to pass it
explicitly into device functions. Finally, we set the value of the ``property`` to the
value specified in the function parameters.

.. code-block:: cpp

   @kernel void constantFill(const dlong Nelements, const dfloat fluid_val,
             const dfloat solid_val, @restrict const dlong* eInfo, @restrict dfloat* property)
   {
     for (dlong e = 0; e < Nelements; ++e ; @outer(0))
     {
       const bool is_solid = eInfo[e];

       for (int n = 0; n < p_Np; ++n ; @inner(0))
       {
         const int id = e * p_Np + n;

         property[id] = fluid_val;

         if (is_solid)
           property[id] = solid_val;
       }
     }
   }

Now, let's look at the slightly more complex ``conductivity`` kernel. Here, our function
signature is very different from that of the ``constantFill`` kernel. While we still
pass the number of elements, we no longer need to check whether we are in a fluid element
or a solid element, since the conductivity for the second passive scalar is going to be
the same in both phases. All that we need to pass in is the coupled scalar ``scalar``, 
or :math:`\phi_0` in our material property correlation :math:`k=5+\phi_0` that we listed
earlier. The ``property`` passed in then should represent the conductivity we are setting.

.. code-block:: cpp

  @kernel void conductivity(const dlong Nelements, @restrict const dfloat* scalar,
            @restrict dfloat* property)
  {
     for (dlong e = 0; e < Nelements; ++e ; @outer(0))
     {
       for (int n = 0; n < p_Np; ++n ; @inner(0))
       {
         const int id = e * p_Np + n;
         const dfloat scalar = scalar[id];

         property[id] = 5.0 + scalar;
       }
     }
  }

A key aspect of writing device kernels is that the device kernel can only operate on
non-pointer objects or pointers to device memory. Whatever the form of your material properties,
you just need to be sure to pass in all necessary information. Now, let's look at the even
more complex ``viscosity`` kernel. Here, we need to pass in the scalar :math:`\phi_0` and the
:math:`z`-coordinate that appear in the viscosity model.

.. code-block:: cpp

  @kernel void viscosity(const dlong Nelements, @restrict const dfloat* scalar,
            @restrict const dfloat* z, @restrict dfloat* property)
  {
     for (dlong e = 0; e < Nelements; ++e ; @outer(0))
     {
       for (int n = 0; n < p_Np; ++n ; @inner(0))
       {
         const int id = e * p_Np + n;
         const dfloat scalar = scalar[id];
         const dfloat z = z[id];

         property[id] = 2.1E-5 * exp(-scalar / 500.0) * (1.0 + z);
       }
     }
  }

The final kernel that wraps up this example is the ``heatCapacity`` kernel.

.. _custom_sources:

Setting Custom Source Terms
---------------------------

Custom source terms can be added to the momentum conservation equation and/or the
energy conservation equation by assigning the ``udf.uEqnSource`` and
``udf.sEqnSource`` function pointers, respectively, to functions with the appropriate signature.
Each of these cases are described separately next. The process is conceptually very similar
to the process for declaring custom properties in :ref:`Setting Custom Properties <custom_properties>`,
so you may find it useful to first review that section.

The Momentum Equation
"""""""""""""""""""""

To set a custom source term for the momentum equation, you must assign the
``udf.uEqnSource`` function pointer to a function with a signature that takes the ``nrs`` pointer
to the nekRS solution object, the simulation time ``time``, the velocity solution on the device
``o_U``, and the momentum source term on the device ``o_FU``. In ``UDF_Setup``,
we need to assign an address to the ``udf.uEqnSource`` function pointer to a function
with the correct signature where we will eventually compute a momentum source. Our
``UDF_Setup`` function would be as follows.

.. code-block:: cpp

  void UDF_Setup(nrs_t * nrs)
  {
    udf.uEqnSource = &custom_source;
  }

Here, ``custom_source`` is our name for a function in the ``.udf`` file that computes the
momentum source. Its name is arbitrary, but it must have the following signature.

.. code-block:: cpp

  void custom_source(nrs_t * nrs, dfloat time, occa::memory o_U, occa::memory o_FU)
  {
    // compute the momentum source
  }

.. note::

  You must place the ``custom_source`` function _before_ ``UDF_Setup`` (and before any other
  function that uses ``custom_source``) in the ``.udf`` file in order for the just-in-time
  compilation to success.

Suppose we would like to add a gravitational force to the :math:`z` momentum equation, of form
:math:`-\rho_fg`. For the momentum equation, the source term is defined on a
per-mass basis; in other words, we must provide the vector :math:`\vec{f}` for a source with
strong form :math:`\rho\vec{f}`. Our custom source function would be as follows.

.. code-block:: cpp

  // declare all kernels we will be writing
  static occa::kernel constantFillKernel;

  void custom_source(nrs_t * nrs, dfloat time, occa::memory o_U, occa::memory o_FU)
  {
    mesh_t * mesh = nrs->mesh;

    // what momentum equation we want to add gravity to
    int component = 2;

    constantFillKernel(nrs->mesh->Nelements, -9.81, component * nrs->fieldOffset, o_FU);
  }

The ``constantFillKernel`` is a user-defined device kernel. This function must now be defined
in the ``.oudf`` file; the name is arbitrary. In order to link against our device kernels, we must
also instruct nekRS to use its just-in-time compilation to build those kernels. We do this in
``UDF_LoadKernels`` by calling the ``udfBuildKernel`` function for the kernel. The second argument
to the ``udfBuildKernel`` function is the name of the kernel, which appears as the actual function
name of the desired kernel in the ``.oudf`` file.

.. code-block:: cpp

  void UDF_LoadKernels(nrs_t * nrs)
  {
    constantFillKernel = udfBuildKernel(nrs, "constantFill");
  }

The ``constantFill`` kernel is now defined in the ``.oudf`` file as follows.

.. code-block:: cpp

  @kernel void constantFill(const dlong Nelements, const dfloat value,
    const int offset, @restrict dfloat * source)
  {
    for (dlong e = 0; e < Nelements; ++e ; @outer(0))
    {
      for (int n = 0; n < p_Np; ++n ; @inner(0))
      {
        const int id = e * p_Np + n + offset;
        source[id] = value;
      }
    }
  }

The Energy Equation
"""""""""""""""""""

.. _nondimensional:

Solving in Non-Dimensional Form
-------------------------------

nekRS can solve its governing equations in either dimensional or non-dimensional form
with careful attention to the specification of the material properties. To solve in
*dimensional* form, the ``density``, ``viscosity``, ``rhoCp``, ``conductivity``, and
``diffusivity`` parameters in the ``.par`` file simply take dimensional forms. Solving
in *non-dimensional* form requires only small changes from the dimensional approach.
For the case of constant properties, the transformation to non-dimensional form is
trivial, but slightly more care is required to solve in non-dimensional form with
variable properties. These two approaches are described next with reference to
the incompressible Navier-Stokes model described in :ref:`Incompressible Flow Model <ins_model>`.

.. _constant_p:

Constant Properties
"""""""""""""""""""

For the case of constant properties for :math:`\rho`, :math:`\mu`, :math:`C_p`,
and :math:`k`, solution in non-dimensional form is achieved by simply specifying
the non-dimensionalized version of these properties in the ``.par`` file. To be explicit,
for the momentum and energy conservation equations, the input parameters should be specified as:

  * ``rho``:math:`\rightarrow` :math:`\rho^\dagger\equiv\frac{\rho}{\rho_0}`
  * ``viscosity``:math:`\rightarrow` :math:`\frac{1}{Re}\mu^\dagger\equiv\frac{\mu_0}{\rho_0UL}\frac{\mu}{\mu_0}`
  * ``rhoCp``:math:`\rightarrow` :math:`\rho^\dagger C_p^\dagger\equiv\frac{\rho}{\rho_0}\frac{C_p}{C_{p,0}}`
  * ``conductivity``:math:`\rightarrow` :math:`\frac{1}{Pe}k^\dagger\equiv\frac{k_0}{\rho_0C_{p,0}UL}\frac{k}{k_0}`

For the :math:`k` and :math:`\tau` equations, if present, the input parameters for
*both* the :math:`k` equation should be specified as:

  * ``rho``:math:`\rightarrow`:math:`1.0`
  * ``diffusivity``:math:`\rightarrow`:math:`\frac{1}{Re}`

Notice that these non-dimensional forms for the :math:`k` and :math:`\tau` equations
are slightly simpler than the forms for the mean momentum and energy equations - this
occurs because nekRS's :math:`k`-:math:`\tau` model is restricted to constant-property
flows, so we do not need to consider :math:`\rho^\dagger\neq 1` or
:math:`\mu^\dagger\neq 1`.

If a volumetric heat source is present, it must also be specified in non-dimensional form
as

.. math::

  \dot{q}^\dagger=\frac{\dot{q}}{\rho_0C_{p,0}U\Delta T/L}

If a source term is present in the momentum conservation equation, that source term
must also be specified in non-dimensional form as

.. math::

   \mathbf s^\dagger=\frac{\mathbf s}{\rho_0U^2/L}

where :math:`\mathbf s` is the source term in the dimensional equation, with dimensions
of mass / square length / square time.

In addition, all boundary conditions must also be non-dimensionalized appropriately.
Some of the more common boundary conditions and their non-dimensionalizations are:

  * fixed velocity: :math:`u_i^\dagger=\frac{u_i}{U}`, i.e. divide all dimensional
    velocity boundary values by :math:`U`
  * fixed temperature: :math:`T^\dagger=\frac{T-T_0}{\Delta T}`, i.e. from all dimensional temperature
    boundary values, first subtract :math:`T_0` and then divide by :math:`\Delta T`
  * fixed pressure: :math:`P^\dagger=\frac{P}{\rho_0U^2}`, i.e. divide all dimensional
    pressure boundary values by :math:`\rho_0U^2`
  * heat flux: :math:`q^\dagger=\frac{q}{\rho_0C_{p,0}U\Delta T}`, i.e. divide all
    dimensional heat flux boundary values by :math:`\rho_0C_{p,0}U\Delta T`
  * turbulent kinetic energy: :math:`k^\dagger=\frac{k}{U^2}`, i.e. divide the dimensional
    turbulent kinetic energy by :math:`U^2`
  * inverse specific dissipation rate: :math:`\tau^\dagger=\frac{\tau}{L/U}`, i.e.
    divide the dimensional inverse specific dissipation rate by :math:`L/U`

If the Prandtl number is unity, then because :math:`Pe\equiv Re\ Pr`, the coefficient on the
diffusion kernel in both the momentum and energy conservation equations will be the same
(for the case of constant properties).

.. note::

  Several of the nekRS input files use syntax inherited from Nek5000 that allows shorthand
  expressions that are often convenient for the Reynolds and Peclet numbers, which appear
  as inverses in the non-dimensional equations. Specifying ``conductivity = -1000`` is
  shorthand for ``conductivity = 1/1000``.

Variable Properties
"""""""""""""""""""

For the case of variable properties, the procedure is similar to the case for constant
properties, except that the properties must be specified in the ``.oudf`` kernels.
It is best practice to simply omit the ``rho``, ``viscosity``, ``rhoCp``, and
``conductivity`` fields from the ``.par`` file entirely. Then, in the ``.oudf`` kernels,
you must include kernels that apply the variable properties in the same manner as in
:ref:`Constant Properties <constant_p>`. See
:ref:`Setting Custom Properties <custom_properties>` for more
information on the kernel setup.

.. _copy_device_to_host:

Copying From Device to Host
---------------------------

All solutions take place on the host, and data transfer of the solution back to the host
must be manually performed by the user if you would like to access ``nrs->U``, ``nrs->p``,
``nrs->cds->S``, or other solution objects, in host-side functions. To copy the solution
from the device to the host, use the ``nek_ocopyFrom(double time, int tstep)`` routine in the
``nekInterfaceAdapter.cpp`` file. This function performs the following actions:

1. Copy the nekRS solution from the nekRS device arrays to the nekRS host arrays - that is,
``nrs->o_U`` is copied to ``nrs->U``, and so on. This
allows you to access the solution on the host as ``nrs->U``, ``nrs->p``, ``nrs->S``, etc.

1. Copy the nekRS solution from the nekRS host arrays to the Nek5000 backend arrays.

If you only want to access the nekRS host side arays such as ``nrs->U``, you can skip the
second part by directly using :term:`OCCA` memory copy functions like the following, which
copies from the device array ``nrs->o_U`` to the host array ``nrs->U``.

.. code-block::

  nrs->o_U.copyTo(nrs->U);

Calculating the Distance to a Wall
----------------------------------

nekRS allows users to access many Nek5000 "backends" through the (optional)
``<case>.usr`` file. A common use case is to calculate the distance from each
:term:`GLL` point to a boundary, such as for setting initial conditions for turbulent quantities
or other closures. The procedure to compute and then use these values is as follows.

First, in the ``usrdat2`` subroutine, make sure that all boundaries for which
you want to compute the distance for are marked as "wall" boundaries in the ``cbc`` array.
In the example shown below, we assume that
the mesh already has sidesets defined in it (assigned through Cubit/gmsh/however else
the mesh was created). We then loop over all the :term:`GLL` points and determine
if the point is on the boundary of interest by checking if the boundary ID is
equal to the sideset of interest. This is done by checking the absolute difference
between the ``bc`` array and the sideset value of interest (in this example, the sideset
is 7). If the boundary ID matches the sideset of interest, then we set the ``cbc`` array
to ``W``, or the character that indicates a no-slip wall boundary.

.. code-block::

  subroutine usrdat2
  include 'SIZE'
  include 'TOTAL'
  integer e,f

  n = lx1*ly1*lz1*nelv
  nxz = nx1*nz1
  nface = 2*ldim

  do iel=1,nelv
  do ifc=1,2*ndim
     if (abs((bc(5,ifc,iel,1)-7.0)).lt.1e-4) cbc(ifc,iel,1)= 'W  '
  enddo
  enddo

  return
  end

In other words,
if your wall boundaries were instead boundaries 3 and 4, the ``if (abs...)`` lines
in the above example would become:

.. code-block::

  if (abs((bc(5,ifc,iel,1)-3.0)).lt.1e-4) cbc(ifc,iel,1)= 'W  '
  if (abs((bc(5,ifc,iel,1)-4.0)).lt.1e-4) cbc(ifc,iel,1)= 'W  '

Next, in the ``usrdat3`` subroutine, you simply need to call the
``dist`` function, which loops over all boundaries with ``W`` type
and determines the distance of all :term:`GLL` points to those boundaries.
The result of the calculation should be stored into the ``nrs_scptr(1)`` pointer,
which is then what we will access in the ``.udf`` file.

.. code-block::

  subroutine usrdat3
  include 'SIZE'
  include 'TOTAL'

  common /scrach_o1/
   w1(lx1*ly1*lz1*lelv)
  ,w2(lx1*ly1*lz1*lelv)
  ,w3(lx1*ly1*lz1*lelv)
  ,w4(lx1*ly1*lz1*lelv)
  ,w5(lx1*ly1*lz1*lelv)

  common /scrach_o2/
   ywd(lx1,ly1,lz1,lelv)

  COMMON /NRSSCPTR/ nrs_scptr(1)
  integer*8         nrs_scptr

  call distf(ywd,7,'W  ',w1,w2,w3,w4,w5)

  nrs_scptr(1) = loc(ywd)

  return
  end

In other words, if your wall boundaries were instead boundaries 3 and 4, the
``call distf...`` lines in the above example would become:

.. code-block::

  call distf(ywd,3,'W  ',w1,w2,w3,w4,w5)
  call distf(ywd,4,'W  ',w1,w2,w3,w4,w5)

Then, you can access the results of the distance-to-wall calculation in the ``.udf``
by assigning a pointer to the ``nek::scPtr(1)`` array. Note that this call must be
within ``UDF_ExecuteStep`` so that the Nek5000 backend will have been called first.

.. code-block::

  void UDF_ExecuteStep(nrs_t * nrs, dfloat time, int tstep)
  {
    double * wall_distance = (double *) nek::scPtr(1);

    // then, you can copy it into some device-side memory so you can use it in
    // BCs if you want
    auto mesh = nrs->meshV;
    int n_gll_points = mesh->Np * mesh->Nelements;
    int write_location = 2; // "slice" into which you want to write, in case nrs->o_usrwrk holds other info
    nrs->o_usrwrk.copyFrom(wall_distance, n_gll_points * sizeof(dfloat), write_location * nrs->fieldOffset * sizeof(dfloat));
  }


.. rubric:: Footnotes

.. [#f1] There are many different ways to write :term:`OCCA` kernels. The examples shown here are by no means the most optimal form, and are only intended for illustration.
.. [#f2] :term:`OCCA` kernels are programmed in OKL, a thin extension to C++. Unfortunately, the ``pygmentize`` Python syntax highlighter does not recognize OKL syntax, so these examples here lack syntax highlighting.

