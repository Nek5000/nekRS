# nekRS .par File Reference

Complete reference for all valid `.par` file sections and keys.

## [GENERAL] Section

| Key | Values | Default | Notes |
|-----|--------|---------|-------|
| `polynomialOrder` | int | **required** | Spectral element degree (typically 5-9) |
| `verbose` | true/false | false | |
| `startFrom` | string | - | Restart file. Modifiers: `+time=<f>`, `+x`, `+u`, `+s`, `+int` |
| `timeStepper` | bdf1/tombo1, bdf2/tombo2, bdf3/tombo3 | tombo2 | |
| `stopAt` | numSteps, endTime, elapsedTime | numSteps | |
| `numSteps` | int | - | |
| `endTime` | float | - | |
| `dt` | float or expr | - | Modifiers: `+targetCFL=<f>`, `+max=<f>`, `+initial=<f>` |
| `advectionSubCyclingSteps` | int, auto | 0 | OIFS subcycling |
| `checkpointControl` | steps, simulationTime | steps | |
| `checkpointInterval` | int/float | 0 | -1 to disable |
| `checkpointPrecision` | 32, 64 | 32 | |
| `dealiasing` | true/false | true | |
| `cubaturePolynomialOrder` | int | 3/2*(P+1)-1 | |
| `constFlowRate` | expression | - | e.g., `meanVelocity=1.0+direction=X` |
| `scalars` | name1, name2, ... | - | Scalar field names |
| `regularization` | none, hpfrt, gjp, avm | - | See regularization section |
| `userSections` | string, ... | - | Custom .par sections to parse |
| `udf` / `oudf` / `usr` | string | - | Override default file names |
| `redirectOutputTo` | string | - | Log file name |

## [OCCA] Section

| Key | Values | Default |
|-----|--------|---------|
| `backend` | SERIAL, CUDA, HIP, DPCPP | from nekrs.conf |
| `deviceNumber` | int, LOCAL-RANK | 0 |

## [MESH] Section

| Key | Values | Default |
|-----|--------|---------|
| `file` | string | `<caseName>.re2` |
| `partitioner` | rcb, rcb+rsb | rcb+rsb |
| `connectivityTol` | float | 0.2 |
| `boundaryIDMap` | int, int, ... | - |
| `hRefine` | int, int, ... | - |

## [PROBLEMTYPE] Section

| Key | Values | Default |
|-----|--------|---------|
| `equation` | stokes, navierStokes | navierStokes |
| `stressFormulation` | true/false | false |

Modifier: `+variableViscosity` on equation enables stress formulation.

## [FLUID VELOCITY] Section

| Key | Values | Default |
|-----|--------|---------|
| `boundaryTypeMap` | BC types (comma-separated) | **required** |
| `density` / `rho` | float or expr | 1.0 |
| `viscosity` / `mu` | float or expr | 1.0 |
| `solver` | CG, CG+BLOCK, GMRES, ... | CG+BLOCK |
| `preconditioner` | jacobi, multigrid | jacobi |
| `residualTol` | float | 1e-4 |
| `initialGuess` | extrapolation, projection | extrapolation |
| `regularization` | hpfrt, gjp, none | - |

## [FLUID PRESSURE] Section

| Key | Values | Default |
|-----|--------|---------|
| `solver` | GMRES+FLEXIBLE+NVECTOR=15 | GMRES+FLEXIBLE |
| `preconditioner` | multigrid | multigrid |
| `residualTol` | float | 1e-4 |
| `initialGuess` | projection, projectionAconj | projectionAconj |
| `boundaryTypeMap` | BC types | - |

## [SCALAR name] Section

| Key | Values | Default |
|-----|--------|---------|
| `boundaryTypeMap` | BC types | - |
| `diffusionCoeff` / `conductivity` | float or expr | - |
| `transportCoeff` / `rhoCp` | float or expr | - |
| `diffusionCoeffSolid` | float or expr | - |
| `transportCoeffSolid` | float or expr | - |
| `mesh` | fluid, fluid+solid | fluid |
| `solver` | CG, GMRES, CVODE, NONE | CG |
| `preconditioner` | jacobi, multigrid | jacobi |
| `residualTol` | float | 1e-4 |
| `checkpointing` | true/false | true |

## [GEOM] Section (Moving Mesh)

| Key | Values | Default |
|-----|--------|---------|
| `solver` | CG+BLOCK, NONE, USER | CG+BLOCK |
| `boundaryTypeMap` | BC types | - |
| `preconditioner` | jacobi, multigrid | jacobi |

## Regularization Options

Format: `regularization = method+param1=val1+param2=val2`

- **hpfrt**: `+nModes=<int>+scalingCoeff=<float>` (both required)
- **gjp**: `+scalingCoeff=<float>`
- **avm**: `+c0+scalingCoeff=<float>` (scalars only)

## Solver Modifiers

- CG: `+combined`, `+block`, `+flexible`, `+maxiter=<int>`
- GMRES: `+flexible`, `+nVector=<int>`, `+maxiter=<int>`, `+ir`

## Preconditioner Modifiers

- multigrid: `+multiplicative`/`+additive`, `+SEMFEM`
- smootherType: ASM, RAS, `+Chebyshev`, `+FourthOptChebyshev`, `+degree=<int>`
- coarseSolver: smoother, jpcg, boomerAMG (`+cpu`/`+device`, `+overlap`)

## Boundary Condition Types

### Velocity BCs
| Type | Aliases | Description |
|------|---------|-------------|
| `zeroDirichlet` | w, wall, inlet | No-slip wall / zero velocity |
| `udfDirichlet` | v | User-defined velocity (calls `udfDirichlet`) |
| `zeroNeumann` | o, O, outlet, outflow | Zero-gradient outflow |
| `zeroDirichletN/zeroNeumann` | slip, sym | Free-slip / symmetry |
| `zeroDirichletX/zeroNeumann` | slipx, symx | Symmetry in X |
| `zeroDirichletY/zeroNeumann` | slipy, symy | Symmetry in Y |
| `zeroDirichletZ/zeroNeumann` | slipz, symz | Symmetry in Z |
| `zeroDirichletN/udfNeumann` | traction, shl | User traction BC |
| `interpolation` | int | NekNek interpolated BC |
| `none` | - | Internal / periodic |

### Scalar BCs
| Type | Description |
|------|-------------|
| `udfDirichlet` | User-defined value (calls `udfDirichlet`) |
| `udfNeumann` | User-defined flux (calls `udfNeumann`) |
| `udfRobin` / `convective` | Mixed BC (calls `udfRobin`) |
| `zeroNeumann` / `insulated` | Zero flux (adiabatic) |
| `interpolation` | NekNek interpolated |
| `none` | Internal / periodic |
