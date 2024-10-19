# Release v24.0

## What is new? 

* FP32 solver mode
* DPCPP backend to support Intel GPUs
* Interpolation based velocity recycling
* [Ascent](https://ascent.readthedocs.io/en/latest/) in situ visualisation plugin
* iofld class reading/writing field files including [ADIOS2](https://adios2.readthedocs.io/) support 
* Addtional output options (element filter and interpolation on uniform grid / different polynomial-order)
* Multi session nek-nek including multi-rate time stepping
* Improved memory management
* CHT nek-nek support
* nek-nek support for nrsqsub scripts
* Improved JIT compilation performance
* HIP support for SEMFEM  
* Aero forces
* opSEM class
* Mesh surface ops
* Linear implicit velocity source term
* Various bug fixes

## Good to know

* after fixing a bug in the linear solver residual norm, iteration counts have increased compared to previous versions
* GPU aware MPI is disabled by default (`NEKRS_GPU_MPI=0`)
* [reproducibility] variable time step controller restricts dt to 5 significant digits
* nrsman <par, env>  can be used to display the par file or environment settings
* HYPRE replaces AmgX
* Field file extension starts with 0-index

## Breaking Changes

This list provides an overview of the most significant changes in this release, although it may not encompass all modifications. We acknowledge that this release introduces several breaking changes. These adjustments were essential to enhance the stability of the user interface in future iterations. We apologize for any inconvenience this may cause.

* run `build.sh` instead of `nrsconfig` to build the code
* change par section `SCALAR00` to `TEMPERATURE` in case it represent indeed a physical temperature
* `nek::userchk` is no longer called automatically during the setup phase 
* host mirrored variables including `nrs->U, cds->S, mesh->x, nrs->usrwrk` have been removed 
* send signal (defined in env-var `NEKRS_SIGNUM_UPD`) to process trigger file `nekrs.upd`
* use `auto foo = platform->deviceMemoryPool.reserve<T>(nWords)` instead of pre-allocated dfloat slices like `platform->o_mempool.slice0`
* change count argument of `occa::memory::slice, occa::memory::copyFrom, occa::memory::copyTo` to number of words instead of bytes 
* use `nekrs_registerPtr` instead of common blocks NRSSCPTR / SCNRS in usr file and access them using `nek::ptr` in udf (see examples)

### Name Changes
* `velocityDirichletConditions` -> `codedFixedValueVelocity` (same for scalars)
* `velocityNeumannConditions` -> `codedFixedGradientVelocity` (same for scalars)
* `nrs->_mesh` -> `cds->mesh[0]`
* `nek::ocopyToNek` -> `nrs->copyToNek`
* `nek::ocopyFromNek` -> `nrs->copyFromNek`
* `nrs->o_FU` -> `nrs->o_NLT`
* `cds->o_FS` -> `cds->o_NLT`
* `occaKernel` -> `deviceKernel`
* `occaProperties` > `deviceKernelProperties`
* `occa::memory` -> `deviceMemory` 
* `nrs->isOutputStep` -> `nrs->checkpointStep`

### Interface Changes 
* define `time` as double (instead of defloat) in all UDF functions
* remove `nrs_t` argument from all UDF functions (nrs object is now globally accessible within udf if the Navier Stokes solver is enabled)
* `nrs_t::userProperties = std::function<void(double)>` -> `udf::properties = std::function<void(nrs_t *, dfloat, occa::memory, occa::memory, occa::memory, occa::memory)>`
* `nrs_t::userVelocitySource = std::function<void(double)>` -> `udf::uEqnSource = std::function<void(nrs_t *, dfloat, occa::memory, occa::memory)>`
* `nrs_t::userScalarSource = std::function<void(double)>` -> `udf::sEqnSource = std::function<void(nrs_t *, dfloat, occa::memory, occa::memory)>`
* `nrs_t::userConvergenceCheck = std::function<bool(int)>` -> `udf::udfconv = std::function<int(nrs_t *, int)>`
* `nrs_t::userDivergence = std::function<void(double)>` -> `udf::udfdif = std::function<void(nrs_t *, dfloat, occa::memory)>`
* `tavg::setup(dlong fieldOffset, const fields& fields)` -> `tavg::setup(nrs_t*)`
* `planarAvg(mesh_t*, const std::string&, int, int, int, int, dlong, occa::memory o_avg)` -> `postProcessing::planarAvg(nrs_t*, const std::string&, int, int, int, int, occa::memory)`
* `::postProcessing` functions are now members of `nrs_t` (except planarAvg)
* remove `nrs_t` argument from `<plugin>::setup`
* `pointInterpolation_t::setPoints(int, dfloat*, dfloat*, dfloat*)` -> `pointInterpolation_t::setPoints(const std::vector<dfloat>&, const std::vector<dfloat>&, const std::vector<dfloat>&)`
* use `iofld` instead of `writeFld`

## Known Bugs / Restrictions

* Code is not fully optimized on CPUs and Intel GPUs
* [507](https://github.com/Nek5000/nekRS/issues/507)
* [485](https://github.com/Nek5000/nekRS/issues/485)
* [258](https://github.com/Nek5000/nekRS/issues/258)

## Thanks to our Contributors

@kris-rowe, @tcew, @yslan, @MalachiTimothyPhillips, @thilinarmtb

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.


# Release v23.0

## What is new? 

* Lagrangian phase model (one-way coupling) 
* Overset grids (neknek)
* Particle tracking 
* Single source udf+oudf
* Device support BoomerAMG
* Improved runtime statistics
* 4th-kind Chebyshev smoothers
* Configureable time averaging 
* Extrapolation initialGuess method
* Scaleable JIT compilation
* Real gas support for lowMach
* More examples
* Various bug fixes 

## Good to know

* [udf] Changes in include files do not trigger a rebuild automatically 
* [udf] Plugins kernels will be loaded automatically (call in `UDF_LoadKernels` no longer required)

## Breaking Changes
* [nrsconfig] Ensure env-vars `CC`, `CXX` and `FC` point to the correct MPI compiler wrappers (see README.md for an example)
* [udf] Plugin header files need to be included explicitly
* [udf] Rename `bc->wrk` => `bc->usrwrk`
* [udf] Update to new API of lowMach plugin (see lowMach example)
* Time step was added to `nekRS::outfld(..., int step, ...)`
* [par] Use `pMGSchedule` instead of `pMultigridCoarsening` (see help for more details)
* [par] Rename writeControl value `runTime` => `simulationTime`
* [par] Remove multigrid qualifier `coarse`
* [par] Remove SEMFEM solver specification from key `preconditioner`, use `semfemSolver` instead
* [par] Replace `stressFormulation = true` by `equation = navierStokes+variableViscosity` 
* [par] Replace bcType `fixedValue` by `codedFixedValue`
* [par] Replace `elasticity` by `pcg+block` for mesh solver
* [okl] Replace `@barrier("local")` by `@barrier()` 
* [oudf] `bc` struct member `trn` was removed
* Use occa::memory mesh_t objects for vgeo, cubvgeo, ggeom, sgeom, LMM, invLMM (no longer mirrored on host)
* All `boundaryIDs` need to be assigned in  `boundaryTypeMap` (use `none` for an internal boundary)

## Known Bugs / Restrictions

* Code is not fully optimized on CPUs in general and Intel GPUs
* [485](https://github.com/Nek5000/Nek5000/issues/485)
* [729](https://github.com/Nek5000/Nek5000/issues/759)
* [258](https://github.com/Nek5000/nekRS/issues/258)

## Thanks to our Contributors

@neil-lindquist, @kris-rowe, @pwang234, @nandu90, @yhaomin2007

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.


# Release v22.0

## What is new? 

* Multi-session (uncoupled) support
* Support unaligned symmetry boundary condition
* Support (unaligned) traction boundary condition
* Better performance on AMD MI-GPUs
* FLOP counters
* Various bug fixes 

## Good to know

* OpenCL support is now disabled by default 

## Breaking Changes

* [udf] Rename `udfBuildKernel` => `oudfBuildKernel`
* [par] Separate details of coarse grid discretization from coarse grid solver
        e.g., `coarseSolver = SEMFEM+AmgX` is replaced by
        `coarseSolver = AmgX` and `coarseGridDiscretization = SEMFEM`
* [par] Remove `preconditioner=semg` and `preconditioner=semg_amg`
* [udf] Rename plug-in name `avg`  => `tavg`
* [udf] Rename `udf.converged` => `udf.timeStepConverged`
* [nrsconfig] Rename env-var `AMGX_ENABLE` => `ENABLE_AMGX`

## Known Bugs / Restrictions

* Mesh solver does not support CHT and unaligned sym/shl BCs
* [729](https://github.com/Nek5000/Nek5000/issues/759)
* [300](https://github.com/Nek5000/nekRS/issues/300)
* [258](https://github.com/Nek5000/nekRS/issues/258)

## Thanks to our Contributors

@tcew, @kris-rowe, @aprilnovak

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

A special shout out to Tim Warburton at VT for tuning some critical kernels. 


# Release v21.1

## What is new? 

* Flexible GMRES
* Constant flow rate
* Time step controller for targetCFL
* Improved runtime statistics
* Support for ROCm version > v4.0
* AVM for scalars
* FEMSEM preconditioner
* Update file (nekrs.upd) for runtime modifications
* Validate key/value input in par
* Various bug fixes 

## Good to know 
* [par] `preconditioner = multigrid` was replaced by `preconditioner = multigrid+coarse`
* [par] Only valid `key/value` pairs will be accepted 
* [par] Default smootherType is `ASM+Chebyshev+degree=2` (instead of degree=1)
* [fld] Only first checkpoint will contain mesh coordinates 
* GMRES is now the default linear solver for pressure (higher memory usage)

## Breaking Changes 

* [udf] Use std namespace qualifier e.g. `std::cout` instead of `cout`
* [udf] Rename `UDF_LoadKernels(nrs_t *nrs)` => `UDF_LoadKernels(occa::properties& kernelInfo)`
* [udf] Replace argument `nrs_t *nrs` by `occa::properties& kernelInfo` in `udfBuildKernel()`, `(plugin)::buildKernel()`
* [udf] `UDF_LoadKernels(occa::properties& kernelInfo)` is no longer optional
* Code crashes (Segmentation fault: invalid permissions) if MPI installation is not GPU aware unless you specify `NEKRS_GPU_MPI=0` in `$NEKRS_HOME/nekrs.conf`

## Known Bugs / Restrictions

* [383](https://github.com/Nek5000/nekRS/issues/383)
* [300](https://github.com/Nek5000/nekRS/issues/300)
* [258](https://github.com/Nek5000/nekRS/issues/258)
* [201](https://github.com/Nek5000/nekRS/issues/201)

## Thanks to our Contributors

@RonRahaman, @aprilnovak, @yslan

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

# Hofix Release v21.0.1

* Update to latest parRSB version
* Fix restart issue if restart time is non-zero
* Fix io-frequency issue
* Fix JIT issue for lowMach
* Disable gs-timers to prevent performance regression

# Release v21.0

## What is new? 

* ASM and RAS smoother + Chebyshev acceleration
* Improved gs performance
* Initial guess projection
* Runtime averages
* Stress formulation
* ALE formulation to support moving meshes
* Linear algebra helpers 
* Various bug fixes 

## What you may have to change to be compatible 

* common block SCRNS was replaced by pointer array NRSSCPTR (see ethier example) 
* boundary device functions and bc struct members in oudf were renamed
* manually copying nek's IC in UDF_Setup() is no longer required 
* nrs->Nlocal was replaced by mesh->Nlocal
* nrs->options was replaced by platform->options
* nrs->linAlg was replaced by platform->linAlg
* nek_copyFrom() was renamed to nek::copyToNek()
* nek_copyTo() was renamed to nek::copyFromNek()
* cds->fieldOffset was replaced by cds->fieldOffset[i] 
* nrs->mesh was replaced by nrs->meshV
* cds->mesh was replaced by cds->mesh[i] 
* nrs->meshT was replaced by cds->mesh[0]
* mesh->rank was replaced by platform->comm.mpiRank
* mesh->comm was replaced by platform->comm.mpiComm
* mesh->device was replaced by platform->device

## Known Bugs 

* [201](https://github.com/Nek5000/nekRS/issues/201)
* [199](https://github.com/Nek5000/nekRS/issues/199)
* [166](https://github.com/Nek5000/nekRS/issues/166)
* [2](https://github.com/Nek5000/nekRS/issues/2)

## Thanks to our Contributors

@RonRahaman, @aprilnovak, @roystgnr, @yslan, @pwang234

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

A special thanks goes to the CAPS Lab at ETHZ who helped to develop the moving mesh support. 

# Release v20.0

## What is new? 

* Initial release

## What you may have to change to be compatible 

* n/a 

## Known Bugs 

[80](https://github.com/Nek5000/nekRS/issues/80)

## Thanks to our Contributors

@AliKarakus, @thilinarmtb, @noelchalmers and @tcew for helping 

We are grateful to all who added new features, filed issues or helped resolve them, 
asked and answered questions, and were part of inspiring discussions.

