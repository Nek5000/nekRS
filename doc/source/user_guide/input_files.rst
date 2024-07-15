.. _input:

The nekRS Input Files
=====================

This page describes the input file structure and syntax needed to run a nekRS simulation.
A nekRS simulation is referred to as a "case," and at a minimum requires four files to run -

* Parameter file, with ``.par`` extension
* Mesh file, with ``.re2`` extension
* User-defined functions for the host, with ``.udf`` extension
* User-defined functions for the device, with ``.oudf`` extension

The "case name" is then the common prefix applied to these files - for instance,
a complete input description with a case name of "eddy" would be given by the files
``eddy.par``, ``eddy.re2``, ``eddy.udf``, and ``eddy.oudf``.
The only restrictions on the case name are:

* It must be used as the prefix on all simulation files, and
* Typical restrictions for naming files for your operating system

The next four sections describe the structure and syntax for each of these four files
for a general case.
Because the :term:`Nek5000` code is a predecessor to
nekRS, some aspects of the current nekRS input file design are selected to enable faster translation of
Nek5000 input files into nekRS input files. Because these
Nek5000-based approaches require proficiency in Fortran, the inclusion of several additional input
files, and in some cases, careful usage of fixed-format text inputs, all
Nek5000-based methods for case setup are referred to here as "legacy" approaches.
All new users are encouraged to adopt the nekRS-based problem setup.

The scope of this page is merely to introduce the format and purpose of the four
files needed to set up a nekRS simulation. Much more detailed instructions are provided
on the :ref:`FAQs <detailed>` page.

.. _parameterFile:

Parameter File (.par)
_____________________

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

Each of the keys and value types are now described for these sections. The
formatting used here to describe valid key, value combinations is as follows.
Take the ``backend`` key in the ``OCCA`` section as an example:

**backend** *(CUDA), CPU, HIP, OPENCL, OPENMP, SERIAL* [``THREAD MODEL``]

Here, ``backend`` is the key, and ``CUDA``, ``CPU``, ``HIP``, ``OPENCL``, ``OPENMP``,
and ``SERIAL`` are all valid values. Defaults are indicated in parentheses - therefore,
if you do not explicitly give the ``backend`` in the ``.par`` file,
the :term:`CUDA` backend is used. Similar conventions are used to describe non-character
type values; for instance, *(3), <int>* indicates that the default value for the indicated
key is 3, but any integer value can be provided.

Most of the values associated with the various keys in the ``.par`` file are read by nekRS
and then saved to various arguments in the ``options`` data structure. The argument
is indicated in this section within square brackets. For example,
the value set by the ``backend`` key is stored in the ``THREAD MODEL`` argument
to ``options``. In other words, if you wanted to grab the value set by the user for the
``backend`` key, and save it in a local variable named ``user_occa_backend``,
you can use the ``getArgs`` function on the ``options`` data structure as follows.

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

.. warning::

  This user guide may quickly become out of date unless developers are careful to keep 
  the keys listed here up to date. A list of possible values is also given in 
  ``doc/parHelp.txt``

nekRS uses just-in-time compilation to allow the incorporation of user-defined functions
into program execution. These functions can be written to allow ultimate flexibility on
the part of the user to affect the simulation, such as to define custom fluid properties,
specify spatially-dependent boundary and initial conditions, and apply post-processing
operations. Some of the parameters in the sections can be overridden through the use of
user-defined functions - see, for example, the ``viscosity`` key in
the ``VELOCITY`` section. This parameter is used to set a constant viscosity, whereas
for variable-property simulations, a user-defined function will override the ``viscosity``
input parameter. A full description of these user-defined functions on the host and
device are described in Sections :ref:`UDF Functions <udf_functions>` and
:ref:`OUDF Functions <oudf_functions>`. So, the description of valid (key, value)
pairs here does not necessarily imply that these parameters reflect the full capabilities
of nekRS.

``BOOMERAMG`` section
^^^^^^^^^^^^^^^^^^^^^

This section is used to describe settings for the (optional) :term:`AMG` solver.

 * **coarsenType** [``BOOMERAMG COARSEN TYPE``]

 * **interpolationType** [``BOOMERAMG INTERPOLATION TYPE``]

 * **iterations** *<int>* [``BOOMERAMG ITERATIONS``]

 * **nonGalerkinTol** [``BOOMERAMG NONGALERKIN TOLERANCE``]

 * **smootherType** [``BOOMERAMG SMOOTHER TYPE``]

 * **strongThreshold** *<double>* [``BOOMERAMG NONGALERKIN TOLERANCE``]

``GENERAL`` section
^^^^^^^^^^^^^^^^^^^

This section is used to describe generic settings for the simulation such as time steppers,
solution order, and file writing control.

* **constFlowRate** *<string>* [``"CONSTANT FLOW RATE = [value is provided]``]
 
  Set a constant flow rate in a given direction. Either ``meanVelocity`` or 
  ``meanVolumetricFlow`` must be provided to set the flow rate,
  and either ``bid`` or ``direction`` must be provided to set the direction.
  The following options are valid:

  * **meanVelocity** *<float>* [``CONSTANT FLOW RATE TYPE = BULK``, ``FLOW RATE``]

    Sets the mean velocity.
  
  * **meanVelocity** *<float>* [``CONSTANT FLOW RATE TYPE = VOLUMETRIC``, ``FLOW RATE``]

    Sets the mean volumetric flow rate.
  
  * **bid** *<int>, <int>* [``CONSTANT FLOW FROM BID``, ``CONSTANT FLOW TO BID``]

    Sets the flow direction based on two boundary IDs.
  
  * **direction** *x, y, z*  [``CONSTANT FLOW DIRECTION``]

    Sets a flow direction parallel to the global coordinate axis.

* **cubaturePolynomialOrder** *<int>* [``CUBATURE POLYNOMIAL DEGREE``]

  Polynomial order for the cubature. If not specified, this defaults to the integer
  closest to :math:`\frac{3}{2}(N + 1)` minus one, where :math:`N` is the polynomial
  order.

  .. TODO: need better description of what cubature is

* **dealiasing** *(true), false*

  If dealiasing is turned on, [``ADVECTION TYPE``] is set to ``CUBATURE+CONVECTIVE``,
  whereas if dealiasing is turned off, [``ADVECTION TYPE``] is set to ``CUBATURE``.

  .. TODO: need better description of what dealiasing is
* **dt** *<string>* [``DT``]

  Time step size. If any of the keyword options ``targetCFL``, ``max`` or ``initial``
  are specified (separated by ``+``), a variable timestep [``VARIABLE DT = TRUE``] 
  is used. Otherwise, ``dt`` is parsed as ``float`` and indicates the time step size.
  
  The following keywords may be given:

  * **targetCFL** *(0.5), <float>* [``TARGET CFL``]: The target :term:`CFL` is also 
    used to set a default for the ``subCyclingSteps``. If not specified, it is given 
    by `max(subcyclingSteps*2, 0.5)``. 
  
  * **max** *(0), <float>* [``MAX DT``]: Largest allowed timestep. If 0 or unset, the 
    option is ignored.

  * **initial** *(0), <float>* [initially written to ``DT``]: initial timestep.

* **elapsedTime** *<double>* [``STOP AT ELAPSED TIME``]

  Elapsed time at which to end the simulation, if using ``stopAt = elapsedTime``.

* **endTime** *<double>* [``END TIME``]

  Final time at which to end the simulation, if using ``stopAt = endTime``.

* **numSteps** *(0), <int>* [``NUMBER TIMESTEPS``]

  Number of time steps to perform, if using ``stopAt = numSteps``. By default, if not
  specified, then it is assumed that no time steps are performed.


* **oudf** *[casename].oudf* [``UDF OKL FILE``]

  File name (including extension) of the ``*.oudf`` file, relative to the current directory.
  By default, the stem of the ``*.par`` file is used as ``casename``.

* **polynomialOrder** *<int>* [``POLYNOMIAL DEGREE``]

  Polynomial order for the spectral element solution. An order of :math:`N` will result
  in :math:`N+1` basis functions for each spatial dimension. The polynomial order is
  currently limited to :math:`N < 10`.

* **startFrom** *<string>* [``RESTART FILE NAME``]

  Absolute or relative path to a nekRS output file from which to start the simulation from.
  When used, the [``RESTART FROM FILE``] option argument is also set to true.
  If the solution in the restart file was obtained with a different polynomial order,
  interpolation is performed to the current simulation settings. To only read select fields
  from the restart file (such as if you wanted to only apply the temperature solution from the
  restart file to the present simulation), append ``+U`` (to read velocity), ``+P`` (to read pressure),
  or ``+T`` (to read temperature) to the end of the restart file name. For instance, if the restart
  file is named ``restart.fld``, using ``restart.fld+T`` will only read the temperature solution.
  If ``startFrom`` is omitted, the simulation is assumed to start based on the user-defined initial conditions at time zero.

* **stopAt** *(numSteps), elapsedTime, endTime*

  When to stop the simulation, either based on a number of time steps *numSteps*, a simulated
  end time *endTime*, or a total elapsed wall time *elapsedTime*. If ``stopAt = numSteps``,
  the ``numSteps`` parameter must be provided. If ``stopAt = endTime``, the ``endTime``
  parameter must be provided. If ``stopAt = elapsedTime``, the ``elapsedTime`` parameter
  must be provided.

* **subCyclingSteps** *(0), <int>, auto* [``SUBCYCLING STEPS``]

  Number of subcycling steps; if ``dt: targetCFL`` is specified, the number of subcycling 
  steps is taken as the integer nearest to half the target :term:`CFL` as given by
  the ``dt: targetCFL`` parameter. In this case, ``auto`` ensures that an error is raised
  if ``dt: targetCFL`` is not specified.

  .. TODO: better description of what subcycling is

* **timeStepper** *(tombo2), bdf1, bdf2, bdf3, tombo1, tombo3* [``TIME INTEGRATOR``]

  The method to use for time stepping. Note that
  if you select any of the :term:`BDF` options, the time integrator is internally set to
  the :term:`TOMBO` time integrator of equivalent order.

* **udf** *[casename].udf* [``UDF FILE``]

  File name (including extension) of the ``*.udf`` file, relative to the current directory.
  By default, the stem of the ``*.par`` file is used as ``casename``.

* **usr** *[casename].usr* [``NEK USR FILE``]

  File name (including extension) of the ``*.usr`` file, relative to the current directory.
  By default, the stem of the ``*.par`` file is used as ``casename``.

* **verbose** *(false), true* [``VERBOSE``]

  Whether to print the simulation results in verbose format to the screen.

* **writeControl** *(timeStep), runTime* [``SOLUTION OUTPUT COTROL``]

  Method to use for the writing of output files, either based on a time step interval with
  *timeStep* (in which case ``SOLUTION OUTPUT CONTROL`` is set to ``STEPS``)
  or a simulated time interval with *runTime* (in which case ``SOLUTION OUTPUT CONTROL``
  is set to ``RUNTIME``).

* **writeInterval** *<double>* [``SOLUTION OUTPUT INTERVAL``]

  Output writing frequency, either in units of time steps for ``writeControl = timeStep`` or
  in units of simulation time for ``writeControl = runTime``. If a runtime step control is
  used that does not perfectly align with the time steps of the simulation, nekRS will write
  an output file on the timestep that most closely matches the desired write interval.

Common keys
^^^^^^^^^^^

These parameters may be specified in any of the ``GENERAL``, ``VELOCITY``, ``TEMPERATURE`` and 
``SCALARXX``  sections. If the parameter is not specified in any given ``VELOCITY``, 
``TEMPERATURE`` or ``SCALARXX`` section, its values are usually inherited from the ``GENERAL``
section.

The key for the ``options`` structure listed here is the ``GENERAL`` key; in the other sections, 
the key is prefixed with the section name.

* **regularization** *("none"), <string>* [``REGULARIZATION METHOD``]
  
  Filtering settings., options are separated by ``+``. This parameter is mutually exclusive
  with the (deprecated) ``filtering`` parameter. The parameter may be specified in any of
  the ``GENERAL``, ``VELOCITY``, ``TEMPERATURE`` and ``SCALARXX``  sections. If the parameter
  is no specified in any given ``VELOCITY``, ``TEMPERATURE`` or ``SCALARXX`` section,
  its values are inherited from the ``GENERAL`` section.

  Filtering is analogous to Nek5000; the ``hpfrt`` filter is described  further in the 
  `Nek5000 documentation <http://nek5000.github.io/NekDoc/problem_setup/filter.html#high-pass-filter>`__.

  The following examples for ``regularization`` are given in ``examples``:

  .. code-block:: ini

    # examples/turbPipePeriodic
    regularization = hpfrt + nModes=1 + scalingCoeff=10

    # examples/double_shear
    regularization=avm+c0+highestModalDecay+scalingCoeff=0.5+rampconstant=1

  * ``none``: regularization is disabled.
  
  * ``hpfrt``: High-pass filter. The following settings apply to this mode:

    * ``nmodes`` *(1), <int>* [``HPFRT MODES``] 
      
      Number of filtered modes :math:`(N-N')`, where :math:`(N)` is the
      polynomial degree and :math:`(N')` the number of fully resolved modes.

    * ``cutoffratio`` *<float>* 
      
      Alternatively, the number of filtered modes can be given 
      by the cutoff ratio, where :math:`\frac{N'+1}{N+1} = {\tt filterCutoffRatio}`.

    * ``scalingcoeff`` *(1.0), <expression>* (required): [``HPFRT STRENGTH``]
      
      filter weight
  
       .. TODO: need better description of what filter weight is
      
  * ``avm+hpfResidual``: use HPF Residual :term:`AVM<AVM>`, or ``avm+highestModalDecay``: 
    use Persson's highest modal decay AVM.
    The AVM is described in [Persson]_, and only allowed for scalars. 
    If specified in ``GENERAL``, the ``regularization`` parameter must be overwritten in the 
    ``VELOCITY`` section. The following settings apply to these modes:

    * ``scalingcoeff`` *(1.0), <expression>* (required) [``REGULARIZATION SCALING COEFF``]
      
      filter weight
  
    * the ``nmodes``, ``cutoffratio`` and ``scalingcoeff`` parameters described above. With 
      ``HighestModalDecay`` mode, ``scalingcoeff`` is interpreted (and overwrites) as 
      ``vismaxcoeff``.

    * ``vismaxcoeff`` *(0.5), <float>* [``REGULARIZATION VISMAX COEFF``]: 
      
      controls maximum artificial viscosity
  
    * ``c0`` [``REGULARIZATION AVM C0``]:
    
      if provided, make viscosity C0 continous across elements
  
    * ``rampconstant`` *(1.0), <float>* [``REGULARIZATION RAMP CONSTANT``]: 
      
      controls ramp to maximum artificial viscosity


``MESH`` section
^^^^^^^^^^^^^^^^

This section is used to describe mesh settings and set up various mesh solvers
for mesh motion.

**partitioner** [``MESH PARTITIONER``]

**solver** *elasticity, none, user*

If ``solver = none``, the mesh does not move and [``MOVING MESH``] is set to false.
Otherwise, the solver is stored in [``MESH SOLVER``]. When ``solver = user``, the
mesh moves according to a user-specified velocity. Alternatively, if
``solver = elasticity``, then the mesh motion is solved with an :term:`ALE` formulation.

``OCCA`` section
^^^^^^^^^^^^^^^^

This section is used to specify the :term:`OCCA` backend for parallelization.

**backend** *(CUDA), CPU, HIP, OPENCL, OPENMP, SERIAL* [``THREAD MODEL``]

OCCA backend; ``CPU`` is the same as ``SERIAL``, and means that parallelism is achieved with
:term:`MPI`.

**deviceNumber** *(LOCAL-RANK), <int>* [``DEVICE NUMBER``]

``PRESSURE`` section
^^^^^^^^^^^^^^^^^^^^

The ``PRESSURE`` section describes solve settings for the pressure equation. Note that
this block is only read if the ``VELOCITY`` block is also present.

.. TODO: This section needs a lot more work describing all the parameters

**downwardSmoother** *ASM, jacobi, RAS* [``PRESSURE MULTIGRID DOWNWARD SMOOTHER``]

**galerkinCoarseOperator** *<bool>* [``GALERKIN COARSE OPERATOR``]

**maxIterations** *<int>* [``PRESSURE MAXIMUM ITERATIONS``]

**pMultigridCoarsening** [``PRESSURE MULTIGRID COARSENING``]

**preconditioner** *jacobi, multigrid, none, semfem, semg* [``PRESSURE PRECONDITIONER``]

The pressure preconditioner to use; ``semg`` and ``multigrid`` both result
in a multigrid preconditioner.

**residualProj** *(true), false* [``PRESSURE RESIDUAL PROJECTION``]

**residualProjectionStart** *<int>* [``PRESSURE RESIDUAL PROJECTION START``]

**residualProjectionVectors** *<int>* [``PRESSURE RESIDUAL PROJECTION VECTORS``]

**residualTol** *<double>* [``PRESSURE SOLVER TOLERANCE``]

Absolute residual tolerance for the pressure solution

**smootherType** *additive, asm, chebyshev, chebyshev+ras, chebyshev+asm, ras* [``PRESSURE MULTIGRID SMOOTHER``]

**solver**

**upwardSmoother** *ASM, JACOBI, RAS* [``PRESSURE MULTIGRID UPWARD SMOOTHER``]

``PROBLEMTYPE`` section
^^^^^^^^^^^^^^^^^^^^^^^

This section is used to control the form of the governing equations used in nekRS.
While individual equations can be turned on/off in the ``VELOCITY``, ``TEMPERATURE``,
and ``SCALARXX`` sections, this block is used for higher-level control of the forms
of those equations themselves.

**equation** *stokes*

Whether to omit the advection term in the conservation of momentum equation, therefore
solving for the Stokes equations. If ``equation = stokes``, then
[``ADVECTION``] is set to false.

**stressFormulation** *(false), true* [``STRESSFORMULATION``]

Whether the viscosity (molecular plus turbulent) is not constant, therefore requiring
use of the full form of the viscous stress tensor :math:`\tau`. By setting ``stressFormulation = false``,
:math:`\nabla\cdot\tau` is represented as :math:`\nabla\cdot\tau=\mu\nabla^2\mathbf u`.
Even if the molecular viscosity is constant, this parameter must be set to ``true``
when using a :term:`RANS` model because the turbulent viscosity portion of the overall
viscosity is not constant.

``SCALARXX`` section
^^^^^^^^^^^^^^^^^^^^

This section is used to define the transport parameters and solver settings for each
passive scalar. For instance, in a simulation with two passive scalars, you would have
two sections - ``SCALAR01`` and ``SCALAR02``, each of which represents a passive scalar.

**boundaryTypeMap** *<string[]>*

Array of strings describing the boundary condition to be applied to each sideset, ordered
by sideset ID. The valid characters/strings are shown in Table
:ref:`Passive Scalar Boundary Conditions <scalar_bcs>`.

**diffusivity** *<double>*

Although this is named ``diffusivity``, this parameter doubly represents the conductivity
governing diffusion of the passive scalar. In other words, the analogue from the
``TEMPERATURE`` section (a passive scalar in its internal representation) is the
``conductivity`` parameter. If a negative value is provided, the
conductivity is internally set to :math:`1/|k|`, where :math:`k` is the value of the
``conductivity`` key. If not specified, this defaults to :math:`1.0`.

**residualProjection** *<bool>*

**residualProjectionStart** *<int>*

**residualProjectionVectors** *<int>*

**residualTol** *<double>*

Absolute residual tolerance for the passive scalar solution

**rho** *<double>*

Although this is name ``rho``, this parameter doubly represents the coefficient on the
total derivative of the passive scalar. In other words, the analogue from the
``TEMPERATURE`` section (a passive scalar in its internal representation) is the
``rhoCp`` parameter. If not specified, this defaults to :math:`1.0`.

``TEMPERATURE`` section
^^^^^^^^^^^^^^^^^^^^^^^

This section is used to define the transport parameters and solver settings for the
temperature passive scalar.

**boundaryTypeMap** *<string[]>*

Array of strings describing the boundary condition to be applied to each sideset, ordered
by sideset ID. The valid characters/strings are shown in Table
:ref:`Passive Scalar Boundary Conditions <scalar_bcs>`.

**conductivity** *<double>* [``SCALAR00 DIFFUSIVITY``]

Constant thermal conductivity; if a negative value is provided, the thermal conductivity
is internally set to :math:`1/|k|`, where :math:`k` is the value of the ``conductivity``
key. If not specified, this defaults to :math:`1.0`.

**residualProj** *<bool>* [``SCALAR00 RESIDUAL PROJECTION``]

**residualProjectionStart** *<int>* [``SCALAR00 RESIDUAL PROJECTION START``]

**residualProjectionVectors** *<int>* [``SCALAR00 RESIDUAL PROJECTION VECTORS``]

**residualTol** *<double>* [``SCALAR00 SOLVER TOLERANCE``]

**rhoCp** *<double>* [``SCALAR00 DENSITY``]

Constant volumetric isobaric specific heat. If not specified, this defaults to :math:`1.0`.

**solver** *none*

You can turn off the solution of temperature by setting the solver to ``none``.

``VELOCITY`` section
^^^^^^^^^^^^^^^^^^^^

This section is used to define the transport properties and solver settings for the
velocity.

**boundaryTypeMap** *<string[]>*

Array of strings describing the boundary condition to be applied to each sideset, ordered
by sideset ID. The valid characters/strings are shown in Table
:ref:`Flow Boundary Conditions <flow_bcs>`. Note that no boundary conditions need to be
specified in the ``PRESSURE`` section, since the form of the pressure conditions are
specified in tandem with the velocity conditions with this parameter.

**density** *<double>* [``DENSITY``]

Constant fluid density. If not specified, this defaults to :math:`1.0`.

**maxIterations** *(200), <int>* [``VELOCITY MAXIMUM ITERATIONS``]

Maximum number of iterations for the velocity solve

**residualProj** *<bool>* [``VELOCITY RESIDUAL PROJECTION``]

**residualProjectionStart** *<int>* [``VELOCITY RESIDUAL PROJECTION START``]

**residualProjectionVectors** *<int>* [``VELOCITY RESIDUAL PROJECTION VECTORS``]

**residualTol** *<double>* [``VELOCITY SOLVER TOLERANCE``]

Absolute tolerance used for the velocity solve.

**solver** *none* [``VELOCITY SOLVER``]

You can turn off the solution of the flow (velocity and pressure) by setting the solver
to ``none``. Otherwise, if you omit ``solver`` entirely, the velocity solve will be turned on.
If you turn the velocity solve off, then you automatically also turn off the pressure solve.

**viscosity** *<double>* [``VISCOSITY``]

Constant dynamic viscosity; if a negative value is provided, the dynamic viscosity is
internally set to :math:`1/|\mu|`, where :math:`\mu` is the value of the ``viscosity`` key.
If not specified, this defaults to :math:`1.0`.

``CASEDATA`` section
^^^^^^^^^^^^^^^^^^^^^

This section may be used to provide custom parameters in the ``.par`` file that are to be read
in the ``.udf`` file. For example, you may specify 

.. code-block:: ini

  [CASEDATA]
    Re_tau = 550

in the ``.par`` file; the parameters should be read in the :ref:`UDF_Setup0 <udf_setup0>` 
function, e.g.

.. code-block:: cpp

  static dfloat Re_tau;
  platform->par->extract("casedata", "re_tau",Re_tau);

NekRS does not check the contents of the ``CASEDATA`` section; such checks may be added in the
``UDF_Setup0`` function as well.

Deprecated parameters
^^^^^^^^^^^^^^^^^^^^^

``GENERAL`` section
"""""""""""""""""""

* **filterCutoffRatio** *<double>* [deprecated, see **regularization**]

  .. TODO: need better description of what filter cutoff ratio is

* **filtering** *hpfrt* [deprecated, see **regularization**]

  If ``filtering = hpfrt``, [``FILTER STABILIZATION``] is set to ``RELAXATION``,
  and ``filterWeight`` must be specified. If ``filtering`` is not specified,
  [``FILTER STABILIZATION``] is set to ``NONE`` by default.

  .. TODO: need better description of what filtering is

*  **filterModes** *<int>* [``HPFRT MODES``] [deprecated, see **regularization**]

  Number of filter modes; minimum value is 1. If not specified, the number of modes
  is set by default to the nearest integer to :math:`(N+1)(1-f_c)`, where :math:`f_c`
  is the filter cutoff ratio.

  .. TODO: need better description of what filter modes is

*  **filterWeight** *<double>* [``HPFRT STRENGTH``] [deprecated, see **regularization**]

  .. TODO: need better description of what filter weight is

Legacy Option (.rea)
^^^^^^^^^^^^^^^^^^^^

An alternative to the use of the ``.par`` file is to use the legacy Nek5000-based ``.rea`` file
to set up the case parameters.
See the ``Mesh File (.re2)`` section of the :term:`Nek5000`
`documentation <http://nek5000.github.io/NekDoc/problem_setup/case_files.html>`__ [#f1]_
for further details on the format for the ``.rea`` file.

The ``.rea`` file contains both simulation parameters (now covered by the ``.par`` file) as well
as mesh information (now covered by the ``.re2`` file). This section
here only describes the legacy approach to setting simulation parameters via the ``.rea`` file.

.. TODO: describe the .rea file approach

Mesh File (.re2)
________________

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

Converting an Existing Commercial Mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most general and flexible approach for creating a mesh is to use commercial meshing software
such as Cubit or Gmsh. After creating the mesh, it must be converted to the ``.re2`` binary format. Depending
on the mesh format (such as Exodus II format or Gmsh format), a conversion script is used to
convert the mesh to ``.re2`` format. See the
:ref:`Converting a Mesh to .re2 Format <converting_mesh>` section for examples demonstrating
conversion of Exodus and Gmsh meshes into ``.re2`` format.

.. _nek5000_mesh:

Nek5000 Script-Based Meshing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A number of meshing scripts ship with the :term:`Nek5000` dependency, which allow
you to directly create ``.re2`` format meshes without the need of commercial meshing
tools. These scripts, such as ``genbox``, take user input related to the desired
grid spacing to generate meshes for fairly simple geometries. Please consult the
`Nek5000 documentation <http://nek5000.github.io/NekDoc/index.html>`__
for more information on the use of these scripts.

Legacy Option (.rea)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An alternative to the use of the ``.re2`` mesh file is to use the legacy Nek5000-based ``.rea`` file
to set up the mesh.
See the ``Mesh File (.re2)`` section of the :term:`Nek5000`
`documentation <http://nek5000.github.io/NekDoc/problem_setup/case_files.html>`__ [#f1]_
for further details on the format for the ``.rea`` file.

The ``.rea`` file contains both simulation parameters (now covered by the ``.par`` file) as well as
mesh information (now covered by the ``.re2`` file). This section
here only describes the legacy approach to setting mesh information via the ``.rea`` file.

The mesh section of the ``.rea`` file can be generated in two different manners -
either by specifying all the element nodes by hand, or with the :term:`Nek5000` mesh
generation scripts introduced in Section :ref:`Nek5000 Script-Based Meshing <nek5000_mesh>`.
Because the binary ``.re2`` format is preferred for very large meshes where memory may be
a concern, the ``.rea`` file approach is considered to be a legacy option.
The mesh portion of the legacy ``.rea``
file can be converted to the ``.re2`` format with the ``reatore2`` script, which also
ships with the :term:`Nek5000` dependency.

.. _udf_functions:

User-Defined Host Functions (.udf)
__________________________________

User-defined functions for the host are specified in the ``.udf`` file. These
functions can be used to perform virtually any action that can be programmed in C++.
Some of the more common examples are setting initial conditions, querying the solution
at regular intervals, and defining custom material properties and source terms. The
available functions that you may define in the ``.udf`` file are as follows. From the
examples shown on the :ref:`Detailed Usage <detailed>` page, you will see that usage
of these functions requires some proficiency in the C++
language as well as some knowledge of the nekRS source code internals.

``UDF_ExecuteStep(nrs_t* nrs, dfloat time, int tstep)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This user-defined function is probably the most flexible of the nekRS user-defined
functions. This function is called once at the start of the simulation just before
beginning the time stepping, and then once per time step after running each step.

``UDF_LoadKernels(nrs_t*  nrs)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This user-defined function is used to load case-specific device kernels that are
used in other UDF functions. For instance, if you add a custom forcing term to the
momentum equations, you need to tell nekRS to compile that kernel by loading it in
this function. The custom material property example shown in the
:ref:`Setting Custom Properties with UDF_Setup <custom_properties>` section
demonstrates how to load kernels with this function. The process is quite simple,
and only involves:

* Declaring all kernels as ``static occa::kernel`` at the top of the ``.udf`` file
* Loading those kernels in ``UDF_LoadKernels``
* Defining those kernels in the device user file (the ``.oudf`` file)

The only kernels in the ``.oudf`` file that don't need to be explicitly loaded are
the boundary condition kernels that ship with nekRS. During the ``.oudf`` just-in-time
compilation, nekRS will search the ``.oudf`` file for any functions that match the
nekRS boundary condition functions, and automatically create and load a kernel based
on the function internals set by the user. For instance, in the ``setOUDF`` function
in the nekRS source code,
the ``.oudf`` file is scanned for a string matching ``scalarDirichletConditions`` (one
of the boundary condition functions in Table :ref:`Passive Scalar Boundary Conditions <scalar_bcs>`).
If this string is found, then the function internals written by the user are cast
into a generic :term:`OCCA` kernel that is then written into a just-in-time compiled
:term:`OKL`-language file at ``.cache/udf/udf.okl``.

.. code-block:: cpp

   found = buffer.str().find("void scalarDirichletConditions");
   if (found == std::string::npos)
     out << "void scalarDirichletConditions(bcData *bc){}\n";

   out <<
     "@kernel void __dummy__(int N) {"
     "  for (int i = 0; i < N; ++i; @tile(16, @outer, @inner)) {}"
     "}";

The ``UDF_LoadKernels`` function is passed the nekRS simulation object ``nrs`` to provide optional
access to the ``occa::properties`` object on the ``nrs->kernelInfo`` object. In
addition to loading kernels, this function can also be used to propagate user-defined
variables to the kernels. See
the :ref:`Defining Variables to Access in Device Kernels <defining_variables_for_device>`
section for a description of this feature.

.. _udf_setup0:

``UDF_Setup0(MPI_Comm comm, setupAide & options)``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

``UDF_Setup(nrs_t* nrs)``
^^^^^^^^^^^^^^^^^^^^^^^^^

This user-defined function is passed the nekRS simulation object ``nrs``. This function
is called once at the beginning of the simulation *after* initializing the mesh, solution
arrays, material property arrays, and boundary field mappings. This function is most
commonly used to:

* Apply initial conditions to the solution
* Assign function pointers to user-defined source terms and material properties

Any other additional setup actions that depend on initialization of the solution arrays
and mesh can of course also be placed in this function.

Other Functions for Custom Sources on the ``udf`` Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the ``UDF_Setup0``, ``UDF_Setup``, ``UDF_ExecuteStep``, and ``UDF_LoadKernels``,
there are other user-defined functions. These functions
are handled in a slightly different manner - rather than be tied to a specific function name
like ``UDF_Setup0``, these functions are provided in terms of generic function pointers to
*any* function (provided the function parameters match those of the pointer). The four
function pointers are named as follows in nekRS:

================== ======================================================== ===================
Function pointer   Function signature                                       Purpose
================== ======================================================== ===================
``udf.converged``  ``f(nrs_t* nrs, int stage)``
``udf.uEqnSource`` ``f(nrs_t* nrs, float t, m o_U, m o_FU)``                momentum source
``udf.sEqnSource`` ``f(nrs_t* nrs, float t, m o_S, m o_SU)``                scalar source
``udf.properties`` ``f(nrs_t* nrs, float t, m o_U, m o_S, m o_Up, m o_Sp)`` material properties
``udf.div``        ``f(nrs_t* nrs, float t, m o_div)``                      thermal divergence
================== ======================================================== ===================

To shorten the syntax above, the type ``m`` is shorthand for ``occa::memory``, and ``f`` is the
name of the function, which can be *any* user-defined name. Other parameters that appear in the
function signatures are as follows:

* ``nrs`` is a pointer to the nekRS simulation object
* ``stage``
* ``t`` is the current simulation time
* ``o_U`` is the velocity solution on the device
* ``o_S`` is the scalar solution on the device
* ``o_FU`` is the forcing term in the momentum equation
* ``o_SU`` is the forcing term in the scalar equation(s)
* ``o_Up`` is the material properties (:math:`\mu` and :math:`\rho`) for the momentum equation
* ``o_Sp`` is the material properties (:math:`k` and :math:`\rho C_p`) for the scalar equation(s)
* ``o_div``

The ``udf.uEqnSource`` allows specification of a momentum source, such as a gravitational force, or
a friction form loss. The ``udf.sEqnSource`` allows specification of a source term for the passive
scalars. For a temperature passive scalar, this source term might represent a volumetric heat source,
while for a chemical concentration passive scalar, this source term could represent a mass
source. See the :ref:`Setting Custom Source Terms <custom_sources>` section for an example
of setting custom source terms.

The ``udf.properties`` allows specification of custom material properties for the flow
and passive scalar equations,
which can be a function of the solution as well as position and time. See the
:ref:`Setting Custom Properties <custom_properties>` section for an example of setting custom
properties.

.. TODO: describe what ``udf.converged`` is

Finally, ``udf.div``
allows specification of the thermal divergence term needed for the low Mach formulation.


Legacy Option (.usr)
^^^^^^^^^^^^^^^^^^^^

The legacy alternative to user-defined functions in the ``.udf`` file is to write
Fortran routines in a ``.usr`` file based on Nek5000 code internals.

.. TODO: describe how to use the ``.usr`` file

.. _oudf_functions:

User-Defined Device Functions (.oudf)
_____________________________________

This file contains all user-defined functions that are to run on the device. These functions include
all functions used to apply boundary conditions that are built in to nekRS, as well as any other
problem-specific device functions.

Boundary Condition Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The type of condition to apply for each boundary is specified by the ``boundaryTypeMap`` parameter
in the ``.par`` file. A character or longer-form word is used to indicate each boundary condition, where the
entries in ``boundaryTypeMap`` are listed in increasing boundary ID order.
However, this single line only specifies the *type* of boundary condition.
If that boundary condition requires additional information, such as a value to impose for
a Dirichlet velocity condition, or a flux to impose for a Neumann temperature condition, then
a device function must be provided in the ``.oudf`` file. A list of all possible boundary
conditions is as follows. For boundary conditions that require additional input from the user,
a device function is also listed. For other boundary conditions that are fully specified simply
by the type of condition (such as a wall boundary condition for velocity, which sets all
three components of velocity to zero without additional user input), no device function is
needed to apply that condition.

.. _flow_bcs:

.. table:: Flow Boundary Conditions

  =========================================== ============================== =============================
  Function                                    Character Map                  Purpose
  =========================================== ============================== =============================
  ``pressureDirichletConditions(bcData* bc)``                                Dirichlet pressure condition
  ``velocityDirichletConditions(bcData* bc)`` ``v``, ``inlet``               Dirichlet velocity condition
  ``velocityNeumannConditions(bcData* bc)``                                  Neumann velocity condition
  N/A                                         ``p``                          Periodic
  N/A                                         ``w``, ``wall``                No-slip wall for velocity
  N/A                                         ``o``, ``outlet``, ``outflow`` Zero-gradient velocity
  N/A                                         ``slipx``                      ?
  N/A                                         ``slipy``                      ?
  N/A                                         ``slipz``                      ?
  N/A                                         ``symx``                       ?
  N/A                                         ``symy``                       ?
  N/A                                         ``symz``                       ?
  =========================================== ============================== =============================

.. _scalar_bcs:

.. table:: Passive Scalar Boundary Conditions

  =========================================== ============================== ===================
  Function                                    Character Map                  Purpose
  =========================================== ============================== ===================
  ``scalarDirichletConditions(bcData* bc)``   ``t``, ``inlet``               Dirichlet condition
  ``scalarNeumannConditions(bcData* bc)``     ``f``, ``flux``                Neumann condition
  N/A                                         ``p``                          Periodic
  N/A                                         ``i``, ``zeroflux``            Zero-gradient
  N/A                                         ``o``, ``outlet``, ``outflow`` Zero-gradient
  =========================================== ============================== ===================

Each function has the same signature, and takes as input the ``bc`` object. This object contains
all information needed to apply a boundary condition - the position, unit normals, and solution
components. The "character map" refers to the character in the ``boundaryTypeMap`` key in the
``.par`` file that will trigger this boundary condition. The character map can either be a single
letter, or a more verbose (and equivalent) string.

The ``scalar``-type boundary conditions
are called for boundary conditions on passive scalars, while the ``pressure``- and ``velocity``-type
conditions are called for the boundary conditions on the flow.

Each of these functions is *only* called on boundaries that contain that boundary. For instance,
if only boundaries 3 and 4 are primitive conditions on velocity, then ``velocityDirichletConditions``
is only called on boundaries 3 and 4. See the :ref:`Setting Boundary Conditions <boundary_conditions>`
section for several examples on how to set boundary conditions with device functions.

.. rubric:: Footnotes

.. [#f1] While the heading for ``Mesh File (.re2)`` seems to suggest that the contents refer only to the ``.re2`` format, the actual text description still points to the legacy ``.rea`` format.
