.. _properties:

Physical properties
===================

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


Constant terms
--------------

Variable terms
--------------
