.. _case:

Case files
==========

A nekRS simulation is referred to as a "case," which utilises a number of files
which are described in this page. An overview of these are presented in the the 
image below .

.. _fig:case_overview:

.. figure:: ../_static/img/overview.svg
   :align: center
   :figclass: align-center
   :alt: An overview of nekRS case files

There are a minimum of three files required to run a case:

* Parameter file, with ``.par`` extension. This sets parameters used by the case
  and can be modified between runs.
* Mesh file, with ``.re2`` extension. This file defines the geometry of the case.
* User-defined host file, with ``.udf`` extension. This is used to set specific
  equations of the case, initial/boundary conditions, data outputs and other user 
  definable behaviour.

With one optional file

* Trigger file, with ``.upd`` extension. This file allows modifications to the 
  simulation during execution.

The "case name" is then the common prefix applied to these files - for instance,
a complete input description with a case name of "eddy" would be given by the files
``eddy.par``, ``eddy.re2``, ``eddy.udf``, and ``eddy.upd``.
The only restrictions on the case name are:

* It must be used as the prefix on all simulation files, and
* Typical restrictions for naming files for your operating system

The following sections describe the structure and syntax for each of these files
for a general case. Because the :term:`Nek5000` code is a predecessor to
nekRS, some aspects of the current nekRS input file design are selected to enable faster translation of
Nek5000 input files into nekRS input files. Because these
Nek5000-based approaches require proficiency in Fortran, the inclusion of several additional input
files, and in some cases, careful usage of fixed-format text inputs, all
Nek5000-based methods for case setup are referred to here as "legacy" approaches.
All new users are encouraged to adopt the nekRS-based problem setup.

.. _parameter_file:

Parameter File (.par)
---------------------

Most information about the problem setup is defined in the parameter file. This file is organized
in a number of sections, each with a number of keys. Values are assigned to these keys in order to
control the simulation settings.

The general structure of the ``.par`` file is as
follows, where ``FOO`` and ``BAR`` are both section names, with a number of (key, value) pairs.

.. code-block:: ini

  [FOO]
    key = value
    baz = bat

  [BAR]
    alpha = beta
    gamma = delta + keyword=value + ... 

The valid sections for setting up a nekRS simulation are:

* ``BOOMERAMG``: settings for the (optional) :term:`AMG` solver
* ``GENERAL``: generic settings for the simulation
* ``MESH``: settings for the mesh
* ``OCCA``: backend :term:`OCCA` device settings
* ``PRESSURE``: settings for the pressure solution
* ``PROBLEMTYPE``: settings for the governing equations
* ``SCALARXX``: settings for the ``XX``-th scalar
* ``TEMPERATURE``: settings for the temperature solution
* ``VELOCITY``: settings for the velocity solution
* ``CASEDATA``: custom settings

Each of the keys and value types are described below. 

**TODO** Description of keys in table

.. code-block:: cpp

  std::string user_occa_backend;
  options.getArgs("THREAD MODEL", user_occa_backend);

In other words, if you have ``backend = CUDA`` in the ``.par`` file, then
``user_occa_backend`` would be set to ``CUDA`` in the above code.

Generally, most ``.par`` settings are not saved to a data structure, so throughout the code
base, whenever information from the ``.par`` file is needed, it is simply
extracted on-the-fly via the ``options`` structure.

nekRS performs validation of the par file. Invalid sections, invalid keys or values,
invalid value combinations, missing values etc. will terminate the NekRS run with a
clear error message. Deprecated attributes will be highlighted. 

nekRS uses just-in-time compilation to allow the incorporation of user-defined functions
into program execution. These functions can be written to allow ultimate flexibility on
the part of the user to affect the simulation, such as to define custom fluid properties,
specify spatially-dependent boundary and initial conditions, and apply post-processing
operations. Some of the parameters in the sections can be overridden through the use of
user-defined functions - see, for example, the ``viscosity`` key in
the ``VELOCITY`` section. This parameter is used to set a constant viscosity, whereas
for variable-property simulations, a user-defined function will override the ``viscosity``
input parameter. A full description of these user-defined functions on the host and
device are described in Sections :ref:`udf_functions`. So, the description of valid (key, value)
pairs here does not necessarily imply that these parameters reflect the full capabilities
of nekRS.

.. literalinclude:: ../../parHelp.txt
   :language: none

Mesh File (.re2)
----------------

The nekRS mesh file is provided in a binary format with a nekRS-specific
``.re2`` extension. This format can be produced by either:

* Converting a mesh made with commercial meshing software to ``.re2`` format, or
* Directly creating an ``.re2``-format mesh with nekRS-specific scripts

There are three main limitations for the nekRS mesh:

* nekRS is restricted to 3-D hexahedral meshes.
* The numeric IDs for the mesh boundaries must be ordered contiguously beginning from 1.
* The ``.re2`` format only supports HEX8 and HEX 20 (eight- and twenty-node) hexahedral elements.

Lower-dimensional problems can be accommodated on these 3-D meshes by applying zero gradient
boundary conditions to all solution variables in directions perpendicular to the
simulation plane or line, respectively. All source terms and material properties in the
governing equations must therefore also be fixed in the off-interest directions.

For cases with conjugate heat transfer, nekRS uses an archaic process
for differentiating between fluid and solid regions. Rather than block-restricting variables to
particular regions of the same mesh, nekRS retains two independent mesh representations
for the same problem. One of these meshes represents the flow domain, while the other
represents the heat transfer domain. The ``nrs_t`` struct, which encapsulates all of
the nekRS simulation data related to the flow solution, represents the flow mesh as
``nrs_t.mesh``. Similarly,
the ``cds_t`` struct, which encapsulates all of the nekRS simulation data related to the
convection-diffusion passive scalar solution, has one mesh for each passive scalar. That is,
``cds_t.mesh[0]`` is the mesh for the first passive scalar, ``cds_t.mesh[1]`` is the mesh
for the second passive scalar, and so on.
Note that only the temperature passive scalar uses the conjugate heat transfer mesh,
even though the ``cds_t`` struct encapsulates information related to all other
passive scalars (such as chemical concentration, or turbulent kinetic energy). All
non-temperature scalars are only solved on the flow mesh.

.. warning::

  When writing user-defined functions that rely on mesh information (such as boundary
  IDs and spatial coordinates), you must take care to use the correct mesh representation
  for your problem. For instance, to apply initial conditions to a flow variable, you
  would need to loop over the number of quadrature points known on the ``nrs_t`` meshes,
  rather than the ``cds_t`` meshes for the passive scalars (unless the meshes are the same,
  such as if you have heat transfer in a fluid-only domain).
  Also note that the ``cds_t * cds`` object will not exist if your problem
  does not have any passive scalars.

nekRS requires that the flow mesh be a subset of the heat transfer mesh. In other words,
the flow mesh always has less than (or equal to, for cases without conjugate heat transfer)
the number of elements in the heat transfer mesh. Creating a mesh for conjugate heat
transfer problems requires additional pre-processing steps that are described in the
:ref:`Creating a Mesh for Conjugate Heat Tranfser <cht_mesh>` section. The remainder
of this section describes how to generate a mesh in ``.re2`` format, assuming
any pre-processing steps have been done for the special cases of conjugate heat transfer.

.. _udf_functions:

User-Defined Host File (.udf)
-----------------------------

The ``.udf`` file is a :term:`OKL` and C++ mixed language source file where user code 
used to formulate the case is placed. This code is placed in various functions
and these can be used to perform virtually any action that can be programmed in
C++. Some of the more common examples are setting initial conditions, querying
the solution at regular intervals, and defining custom material properties and
source terms. The available functions that you may define in the ``.udf`` file
are as follows.

OKL block
"""""""""

The ``.udf`` typically has a ``#ifdef __okl__`` block near the start which is 
where all OKL code will be placed that will run on the compute backed specified to
:term:`OCCA`. The most frequent use of this block is to provide the functions 
for boundary conditions that require additional information, such as a value to
impose for a Dirichlet velocity condition, or a flux to impose for a Neumann
temperature condition. Additional user functions may be placed in this block to
allow advanced modification of the simulation or functionality such as calculating
exact values at a specified time point.

.. tip::

  If the user generated functions are sufficiently large, or in older nekRS examples 
  you may see a ``.oudf`` file which is included within the ``ifdef`` block 
  instead of the functions being in the ``.udf`` file.

.. code-block::
  
  #ifdef __okl__

  @kernel void computeexact(const dlong Ntotal)
  {
    for (dlong n = 0; n < Ntotal; ++n; @tile(p_blockSize, @outer, @inner)) {
      if (n < Ntotal) {
        // some code
      }
    }
  }

  void velocityDirichletConditions(bcData *bc)
  {
    // some code
    bc->u = u;
    bc->v = v;
    bc->w = w;
  }

  void scalarDirichletConditions(bcData *bc)
  {
    // some code
    bc->s = s
  }

  void scalarNeumannConditions(bcData *bc)
  {
    bc->flux = tflux;
  }

.. _udf_setup0:

UDF_Setup0
""""""""""

This user-defined function is passed the nekRS :term:`MPI` communicator ``comm`` and a data
structure containing all of the user-specified simulation options, ``options``. This function is
called once at the beginning of the simulation *before* initializing the nekRS internals
such as the mesh, solvers, and solution data arrays. Because virtually no aspects of
the nekRS simulation have been initialized at the point when this function is called,
this function is primarily used to modify the user settings. For the typical user,
all relevant settings are already exposed through the ``.par`` file; any desired
changes to settings should therefore be performed by modifying the ``.par`` file.

This function is intended for developers or advanced users to overwrite any user
settings that may not be exposed to the ``.par`` file. For instance, setting
``timeStepper = tombo2`` in the ``GENERAL`` section triggers a number of other internal
settings in nekRS that do not need to be exposed to the typical user, but that perhaps
a developer may want to modify for testing purposes.

UDF_Setup
"""""""""

This user-defined function is passed the nekRS simulation object ``nrs``. This function
is called once at the beginning of the simulation *after* initializing the mesh, solution
arrays, material property arrays, and boundary field mappings. This function is most
commonly used to:

* Apply initial conditions to the solution
* Assign function pointers to user-defined source terms and material properties

Any other additional setup actions that depend on initialization of the solution arrays
and mesh can of course also be placed in this function.

UDF_ExecuteStep
"""""""""""""""

This user-defined function is probably the most flexible of the nekRS user-defined
functions. This function is called once at the start of the simulation just before
beginning the time stepping, and then once per time step after running each step.

.. _trigger_file:

Trigger Files (.upd)
--------------------

**TODO** Full description

Allows modifications to the simulation during execution. Can be edited and then
notify of changes through sending a signal MPI rank 0.

.. rubric:: Footnotes

.. [#f1] While the heading for ``Mesh File (.re2)`` seems to suggest that the contents refer only to the ``.re2`` format, the actual text description still points to the legacy ``.rea`` format.
