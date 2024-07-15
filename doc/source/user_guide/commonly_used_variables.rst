.. _commonly_used_variables:

Commonly-Used Variables in nekRS
================================

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


Mesh
----

This section describes commonly-used variables related to the mesh, which are all stored
on data structures of type ``mesh_t``. nekRS uses an archaic approach for conjugate heat
transfer applications, i.e. problems with separate fluid and solid domains. For problems
without conjugate heat transfer, all mesh information is stored on the ``nrs->mesh`` object,
while for problems with conjugate heat transfer, all mesh information is stored on the
``nrs->cds->mesh`` object. More information is available in the
:ref:`Creating a Mesh for Conjugate Heat Transfer <cht_mesh>` section. To keep the following
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
--------------------------------------------

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
------------------------------------------------------

This section describes the members on the ``cds`` object, which consist of user settings as well as the
passive scalar solution. Note that, from :ref:`Flow Solution Fields and Simulation Settings <flow_vars>`,
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

