.. _intro:

Introduction
============

This section will introduce some main concepts that are needed to setup cases in
nekRS.

.. _nondimensional:

Solving in Dimensional vs Non-Dimensional Form
----------------------------------------------

nekRS can solve its governing equations in either dimensional or non-dimensional form
with careful attention to the specification of the material properties. To solve in
*dimensional* form, the ``density``, ``viscosity``, ``rhoCp``, ``conductivity``, and
``diffusivity`` parameters in the ``.par`` file simply take dimensional forms. Solving
in *non-dimensional* form requires only small changes from the dimensional approach.
For the case of constant properties, the transformation to non-dimensional form is
trivial, but slightly more care is required to solve in non-dimensional form with
variable properties. These two approaches are described next with reference to
the incompressible Navier-Stokes model described in :ref:`ins_model`.

It is recommended to use non-dimensional solves and the other sections of the
documentation will use this as a default.

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
:ref:`constant_p`. See
:ref:`custom_properties` for more
information on the kernel setup.

.. _compute_backend_abstraction:

Compute Backend Abstraction (OCCA)
----------------------------------

One important overarching feature of nekRS is the use of :term:`OCCA` to provide a layer
of abstraction of the potential compute backends (E.G. CPU, GPU's and Intel XPU's)
so a universal language can be used to program the compute intensive areas of a case.
The two main elements of this abstraction is to provide mechanisms to transfer 
relevant data into the memory of the compute target, and a way to write functions 
that can be executed on the compute target.

Here we introduce these elements in the most relevant way to nekRS, but further
information can be found in the `OCCA documentation <https://libocca.org/>`_. The
sections below all refer to code that will be present within the ``.udf`` file 
(see :ref:`udf_functions` for more details)

.. _occa_memory:

Memory
""""""

Memory management is done through the C++ API which allows the user to make data
available on the compute backend device (sometimes referred to as the device) and
copy data into this for future use. 

**TODO** - explanation of any automatic copying??

Typically, relevant fields should be created and initialised in the 
`UDF_loadKernels` function:

.. code-block::

  void UDF_LoadKernels(deviceKernelProperties& kernelInfo)
  { 
    kernelInfo.define("<p_VARIABLE>") = <VALUE>;
  }

.. tip::
  p_ and o_ prefixing

.. _occa_functions:

Functions
"""""""""

The :term:`OKL` language extends C with keywords allowing functions to be written 
in a consistent language which are translated to device specific code (E.G. CUDA).
These functions should typically be in the ``.udf`` file within a ``#ifdef __okl__``
block, and are preceded with a ``@kernel`` keyword. This block would also have
any standard functions that would be required for relevant boundary conditions
(see :ref:`boundary_conditions`). Below is an example showing both of these 
types of function.

.. code-block::

  #ifdef __okl__
  @kernel void sample_function()
  {
    // some code
  }

  void velocityDirichletConditions(bcData *bc)
  {
    // some code
    bc->u = u;
    bc->v = v;
    bc->w = w;
  }
  #endif // __okl__

.. _data_structures:

Data Structures
---------------

UDF Only??

To become a proficient user of nekRS requires some knowledge of the data structures
used to store the mesh, solution fields, and simulation settings. While many
commercial :term:`CFD<CFD>` codes have developed user interfaces that allow most user
code interactions to occur through a :term:`GUI<GUI>` or even a text-based format, nekRS
very much remains a research tool. As such, even "routine" actions such as setting
boundary and initial conditions requires an understanding of the source code structure in
nekRS. This requirement is advantageous from a flexibility perspective, however, because
almost any user action that can be written in C++ ``.udf`` or :term:`OKL<OKL>` in ``.oudf``
files can be incorporated into a nekRS simulation.

This page contains a summary of some of the most commonly-used variables and structures
used to interact with nekRS. For array-type variables, the size of the array is also listed
in terms of the length of each dimension of that array. For instance, if the size of an array
is ``Nelements * Np``, then the data is stored first by each element, and second by each
quadrature point. If the variable is not an array type, the size is shown as ``1``.

Some variables have an equivalent form that is stored on the device that can be accessed
in device kernels. All such device variables and
arrays that live on the device by convention are prefixed with ``o_``. That is, ``mesh->x``
represents all the :math:`x`-coordinates of the quadrature points, and is stored on the host.
The same data, but accessible on the device, is ``mesh->o_x``. Not all variables and arrays
are automatically available on both the host and device, but those that are available are
indicated with a :math:`\checkmark` in the "Device?" table column.

Platform
""""""""

.. _fig:platform_class:

.. figure:: ../doxygen/doxygen_html/structplatform__t__coll__graph.png
   :align: center
   :figclass: align-center
   :alt: Class diagram of the major elements of the platform class


Mesh
""""
.. _fig:mesh_class:

.. figure:: ../doxygen/doxygen_html/classnrs__t__coll__graph.png
   :align: center
   :figclass: align-center
   :alt: Class diagram of the major elements of the Mesh class

This section describes commonly-used variables related to the mesh, which are all stored
on data structures of type ``mesh_t``. nekRS uses an archaic approach for conjugate heat
transfer applications, i.e. problems with separate fluid and solid domains. For problems
without conjugate heat transfer, all mesh information is stored on the ``nrs->mesh`` object,
while for problems with conjugate heat transfer, all mesh information is stored on the
``nrs->cds->mesh`` object. More information is available in the
:ref:`cht_mesh` section. To keep the following
summary table general, the variable names are referred to simply as living on the ``mesh``
object, without any differentiation between whether that ``mesh`` object is the object on
``nrs`` or ``nrs->cds``.

Some notable points of interest that require additional comment:

* The :term:`MPI<MPI>` communicator is stored on the mesh, since domain decomposition
  is used to divide the mesh among processes. *Most* information stored on the ``mesh`` object
  strictly refers to the portion of the mesh "owned" by the current process. For instance,
  ``mesh->Nelements`` only refers to the number of elements "owned" by the current process
  (``mesh->rank``), not the total number of elements in the simulation mesh. Any exceptions
  to this process-local information is noted as applicable.

================== ============================ ================== =================================================
Variable Name      Size                         Device?            Meaning
================== ============================ ================== =================================================
``comm``           1                                               MPI communicator
``device``         1                                               backend device
``dim``            1                                               spatial dimension of mesh
``elementInfo``    ``Nelements``                                   phase of element (0 = fluid, 1 = solid)
``EToB``           ``Nelements * Nfaces``       :math:`\checkmark` boundary ID for each face
``N``              1                                               polynomial order for each dimension
``NboundaryFaces`` 1                                               *total* number of faces on a boundary (rank sum)
``Nelements``      1                                               number of elements
``Nfaces``         1                                               number of faces per element
``Nfp``            1                                               number of quadrature points per face
``Np``             1                                               number of quadrature points per element
``rank``           1                                               parallel process rank
``size``           1                                               size of MPI communicator
``vmapM``          ``Nelements * Nfaces * Nfp`` :math:`\checkmark` quadrature point index for faces on boundaries
``x``              ``Nelements * Np``           :math:`\checkmark` :math:`x`-coordinates of quadrature points
``y``              ``Nelements * Np``           :math:`\checkmark` :math:`y`-coordinates of quadrature points
``z``              ``Nelements * Np``           :math:`\checkmark` :math:`z`-coordinates of quadrature points
================== ============================ ================== =================================================

.. _flow_vars:

Flow Solution Fields and Simulation Settings
""""""""""""""""""""""""""""""""""""""""""""

This section describes the members on the ``nrs`` object, which consist of user settings as well as the flow
solution. Some of this information is simply assigned a value also stored on the ``nrs->mesh`` object.
Some notable points that require additional comment:

* Like the mesh object, the solution fields are stored only on a per-rank basis. That is, ``nrs->U`` only
  contains the velocity solution for the elements "owned" by the current process.
* Solution arrays with more than one component (such as velocity, in ``nrs->U``) are indexed according
  to a ``fieldOffset``. This offset is chosen to be larger than the *actual* length of the velocity
  solution (which is the total number of quadrature points on that rank, or ``nrs->Nlocal``) due to
  performance reasons. That is, you should use the ``fieldOffset`` to index between components, but
  within a single component, you should not attempt to access entries with indices between
  ``i * (fieldOffset - Nlocal)``, where ``i`` is the component number, because those values are not actually
  used to store the solution (they are the end of a storage buffer).

Some members only exist on the device - in this case, the variable name shown in the first column
explicitly shows the ``o_`` prefix to differentiate that this member is not available in this form
on the host. For instance, the ``o_mue`` member is only available on the device - there is no
corresponding array ``nrs->mue`` member.

================== ================================= ================== ======================================================================================================
Variable Name      Size                              Device?            Meaning
================== ================================= ================== ======================================================================================================
``cds``            1                                                    convection-diffusion solution object
``cht``            1                                                    whether the problem contains conjugate heat transfer
``dim``            1                                                    spatial dimension of ``nrs->mesh``
``dt``             3                                                    time step for previous 3 time steps
``fieldOffset``    1                                                    offset in flow solution arrays to access new component
``FU``             ``NVfields * nEXT * fieldOffset`` :math:`\checkmark` source term for each momentum equation for each step in the time stencil
``isOutputStep``   1                                                    if an output file is written on this time step
``lastStep``       1                                                    if this time step is the last time step of the run
``mesh``           1                                                    mesh used for the flow simulation
``nEXT``           1                                                    number of time steps in the time derivative stencil
``NiterU``         1                                                    number of iterations taken in last velocity solve
``NiterP``         1                                                    number of iterations taken in last pressure solve
``Nlocal``         1                                                    number of quadrature points local to this process
``Nscalar``        1                                                    number of passive scalars to solve for
``NTfields``       1                                                    number of flow-related fields to solve for (:math:`\vec{V}` plus :math:`T`)
``NVfields``       1                                                    number of velocity fields to solve for
``o_mue``          ``fieldOffset``                   :math:`\checkmark` total dynamic viscosity (laminar plus turbulent) for the momentum equation
``options``        1                                                    object containing user settings from ``.par`` file
``o_rho``          ``fieldOffset``                   :math:`\checkmark` density for the momentum equation
``P``              ``fieldOffset``                   :math:`\checkmark` pressure solution for most recent time step
``prop``           ``2 * fieldOffset``               :math:`\checkmark` total dynamic viscosity (laminar plus turbulent) and density (in this order) for the momentum equation
``U``              ``NVfields * fieldOffset``        :math:`\checkmark` velocity solution for all components for most recent time step
================== ================================= ================== ======================================================================================================

Passive Scalar Solution Fields and Simulation Settings
""""""""""""""""""""""""""""""""""""""""""""""""""""""

This section describes the members on the ``cds`` object, which consist of user settings as well as the
passive scalar solution. Note that, from :ref:`flow_vars`,
the ``cds`` object is itself stored on the ``nrs`` flow solution object. Many of these members are
copied from the analogous variable on the ``nrs`` object. For instance, ``cds->fieldOffset`` is simply
set equal to ``nrs->fieldOffset``. In a few cases, however, the names on the ``cds`` object differ
from the analogous names on the ``nrs`` object, such as for ``cds->NSfields`` and ``nrs->Nscalar``, which
contain identical information.

================== ============================== ================== ======================================================================================================
Variable Name      Size                           Device?            Meaning
================== ============================== ================== ======================================================================================================
``fieldOffset``    1                                                 offset in passive scalar solution arrays to access new component
``NSfields``       1                                                 number of passive scalars to solve for
``o_diff``         ``NSfields * fieldOffset``     :math:`\checkmark` diffusion coefficient (laminar plus turbulent) for the passive scalar equations
``o_rho``          ``NSfields * fieldOffset``     :math:`\checkmark` coefficient on the time derivative for the passive scalar equations
``prop``           ``2 * NSfields * fieldOffset`` :math:`\checkmark` diffusion coefficient (laminar plus turbulent) and coefficient on the time derivative (in this order) for the passive scalar equations
================== ============================== ================== ======================================================================================================

