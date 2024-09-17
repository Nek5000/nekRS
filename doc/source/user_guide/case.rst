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
"""""""""""""""""""""

This section is used to describe settings for the (optional) :term:`AMG` solver.

 * **coarsenType** [``BOOMERAMG COARSEN TYPE``]

 * **interpolationType** [``BOOMERAMG INTERPOLATION TYPE``]

 * **iterations** *<int>* [``BOOMERAMG ITERATIONS``]

 * **nonGalerkinTol** [``BOOMERAMG NONGALERKIN TOLERANCE``]

 * **smootherType** [``BOOMERAMG SMOOTHER TYPE``]

 * **strongThreshold** *<double>* [``BOOMERAMG NONGALERKIN TOLERANCE``]

``GENERAL`` section
"""""""""""""""""""

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
"""""""""""

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
""""""""""""""""

This section is used to describe mesh settings and set up various mesh solvers
for mesh motion.

**partitioner** [``MESH PARTITIONER``]

**solver** *elasticity, none, user*

If ``solver = none``, the mesh does not move and [``MOVING MESH``] is set to false.
Otherwise, the solver is stored in [``MESH SOLVER``]. When ``solver = user``, the
mesh moves according to a user-specified velocity. Alternatively, if
``solver = elasticity``, then the mesh motion is solved with an :term:`ALE` formulation.

``OCCA`` section
""""""""""""""""

This section is used to specify the :term:`OCCA` backend for parallelization.

**backend** *(CUDA), CPU, HIP, OPENCL, OPENMP, SERIAL* [``THREAD MODEL``]

OCCA backend; ``CPU`` is the same as ``SERIAL``, and means that parallelism is achieved with
:term:`MPI`.

**deviceNumber** *(LOCAL-RANK), <int>* [``DEVICE NUMBER``]

``PRESSURE`` section
""""""""""""""""""""

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
"""""""""""""""""""""""

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
""""""""""""""""""""

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
"""""""""""""""""""""""

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
""""""""""""""""""""

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
""""""""""""""""""""

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
"""""""""""""""""""""

``GENERAL`` section
^^^^^^^^^^^^^^^^^^^

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
""""""""""""""""""""

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

TODO Full description

Allows modifications to the simulation during execution. Can be edited and then
notify of changes through sending a signal MPI rank 0.

.. rubric:: Footnotes

.. [#f1] While the heading for ``Mesh File (.re2)`` seems to suggest that the contents refer only to the ``.re2`` format, the actual text description still points to the legacy ``.rea`` format.
