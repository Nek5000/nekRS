# Release v22.0

## What is new? 

* Multi-session (uncoupled) support
* Support unaligned symmetry boundary condition
* Support (unaligned) traction boundary condition
* Better performance on AMD MI-GPUs
* FLOP counters
* Various bug fixes 

## Good to know

* OpenCL support is disabled by default 

## Breaking Changes

* [udf] Rename `udfBuildKernel` => `oudfBuildKernel`
* [par] Separate details of coarse grid discretization from coarse grid solver
        e.g., `coarseSolver = SEMFEM+AmgX` is replaced by
        `coarseSolver = AmgX` and `coarseGridDiscretization = SEMFEM`
* [par] Remove `preconditioner=semg` and `preconditioner=semg_amg`
* [udf] Rename plug-in name `avg`  => `tavg`
* [udf] Rename `udf.converged` => `udf.timeStepConverged`

## Known Bugs / Restrictions

* Mesh solver does not support CHT and unaligned sym/shl BCs
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

