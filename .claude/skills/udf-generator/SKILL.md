---
name: udf-generator
description: >
  Generate nekRS simulation case files (.udf, .par, .oudf) from a user's
  physics description. Use when the user describes a CFD simulation setup
  including flow conditions (Reynolds number, Mach number), boundary
  conditions (inlet, outlet, wall), geometry, or physical models (turbulence,
  heat transfer, multiphysics). Also use when the user asks to create, modify,
  or set up a nekRS case, even if they don't mention "UDF" explicitly.
compatibility: Requires nekRS codebase. Designed for Claude Code.
allowed-tools: Bash(mkdir:*) Read Write Edit Glob Grep
metadata:
  author: nekRS
  version: "1.0"
---

# nekRS UDF Case Generator

Generate complete nekRS simulation cases from user physics descriptions.

## Workflow

### Step 1: Gather Requirements

Extract these from the user's prompt. If any **critical** item is missing, ask before proceeding.

**Critical (must know):**
- Flow type: incompressible, low-Mach, Stokes
- Boundary conditions: what happens at each boundary (wall, inlet, outlet, symmetry)
- Fluid properties: Reynolds number OR viscosity/density directly
- Domain: geometry description or existing `.re2` mesh file

**Important (ask if ambiguous):**
- Scalar transport: temperature, passive scalars, species
- Turbulence modeling: DNS, LES (with regularization), RANS (k-tau)
- Time stepping: steady vs transient, timestep size, duration
- Source terms: body forces, buoyancy, Coriolis

**Optional (use sensible defaults):**
- Polynomial order (default: 7)
- Time stepper (default: tombo2)
- Solver tolerances (default: 1e-6 velocity, 1e-4 pressure)
- Checkpointing (default: every 500 steps)

### Step 2: Determine Case Type

Match the user's description to one of these physics categories to choose the right pattern:

| Category | Key Indicators | Callbacks Needed | Plugins |
|----------|---------------|-----------------|---------|
| Laminar flow | Low Re, no turbulence model | Minimal or none | - |
| Turbulent flow (DNS/LES) | High Re, no model | `userSource` (forcing) | tavg |
| RANS | Mentions k-tau, SST | `userSource`, `userProperties` | RANSktau, tavg |
| Heat transfer | Temperature, Nusselt, Prandtl | `userProperties` (if variable) | - |
| Conjugate HT | Solid+fluid thermal | `userProperties` | - |
| Low-Mach | Variable density, large dT | `userDivergence`, `userSource` | lowMach |
| Moving mesh | ALE, piston, oscillating | mesh velocity kernel | - |
| Particle tracking | Lagrangian, droplets | LPM callbacks | lpm |

### Step 3: Create Case Directory

Create the case directory under `examples/`:

```bash
mkdir -p examples/<caseName>
```

Use the name the user provides, or infer a descriptive snake_case name from the physics (e.g., `heated_channel`, `turbulent_pipe`, `lid_driven_cavity`).

### Step 4: Generate Files

Generate these files in order. Read reference files as needed:

1. **`<caseName>.par`** — Read [PAR_REFERENCE.md](references/PAR_REFERENCE.md) for all valid keys
2. **`<caseName>.udf`** — Read [UDF_PATTERNS.md](references/UDF_PATTERNS.md) for callback patterns
3. **`<caseName>.oudf`** (if needed) — Read [BOUNDARY_CONDITIONS.md](references/BOUNDARY_CONDITIONS.md) for BC function signatures

### Step 5: Validate

Check the generated files for consistency:

- [ ] Every boundary ID in `boundaryTypeMap` has a corresponding handler if `udfDirichlet`/`udfNeumann`/`udfRobin`
- [ ] Scalar names in `[GENERAL] scalars = ...` match `[SCALAR name]` sections
- [ ] `udfDirichlet` in `.oudf` handles all fields that use it (check with `isField()`)
- [ ] Properties (viscosity, density) are consistent with stated Reynolds number
- [ ] If RANS: scalars include `K` and `TAU`, and `RANSktau` plugin is initialized
- [ ] If lowMach: `userDivergence` is registered and `lowMach::setup()` is called
- [ ] `[CASEDATA]` parameters match what `UDF_Setup0` extracts

## Gotchas

- Boundary IDs in `boundaryTypeMap` are positional: first type = boundary ID 1, second = ID 2, etc.
- `isField()` comparisons in `.oudf` are **case-insensitive** but must match the field name exactly (e.g., `"fluid velocity"`, `"scalar temperature"`, `"scalar00"`)
- Scalar names from `.par` map to field names: `[SCALAR TEMPERATURE]` becomes field `"scalar temperature"` in `isField()`. For unnamed scalars, use `"scalar00"`, `"scalar01"`, etc.
- `platform->par->extract()` keys are **case-insensitive**
- OKL kernels in `#ifdef __okl__` blocks use `@kernel`, `@outer`, `@inner`, `@tile`, `@shared`, `@exclusive` decorators — NOT CUDA/HIP syntax
- `dfloat` is the nekRS floating-point type (double or float depending on build) — always use it instead of `double`/`float` in UDF code
- `dlong` is the nekRS integer type for global indices
- Device memory uses `o_` prefix convention (e.g., `o_U` is velocity on GPU)
- Field offset: velocity components are stored contiguously as `[ux | uy | uz]` each of size `fieldOffset`
- The `usrwrk` array in `bcData` is a **read-only** pointer to a global workspace set from the UDF side
- `UDF_Setup0` runs **before** the mesh and solvers exist — only `platform->par` and MPI are available
- `UDF_Setup` runs **after** solvers are created — this is where you access `nrs`, `mesh`, register callbacks
- For constant pressure-gradient driven flows, prefer `constFlowRate` in `[GENERAL]` over manual forcing
- `regularization = hpfrt` requires `+nModes=<int>+scalingCoeff=<float>` — omitting either causes errors
- Never include `ci.inc` or CI-related code in generated cases — that's for regression tests only

## File Templates

### Minimal `.par` Template

```ini
[GENERAL]
polynomialOrder = 7
numSteps = 1000
dt = 1e-03
timeStepper = tombo2
checkpointInterval = 500

[FLUID VELOCITY]
boundaryTypeMap = <types matching mesh boundary IDs>
residualTol = 1e-06
viscosity = <value>
density = 1.0

[FLUID PRESSURE]
residualTol = 1e-04
```

### Minimal `.udf` Template

```cpp
#include "nrs.hpp"
#include "udf.hpp"

#ifdef __okl__
// OKL kernels here (boundary conditions, forcing, etc.)
#endif

void UDF_Setup0(MPI_Comm comm, setupAide &options)
{
  // Extract [CASEDATA] parameters
  // platform->par->extract("casedata", "key", variable);
}

void UDF_Setup()
{
  // Register callbacks
  // nrs->userSource = &myForcing;
  // nrs->userProperties = &myProperties;
}

void UDF_ExecuteStep(double time, int tstep)
{
  // Per-timestep operations (diagnostics, I/O)
}
```

### Minimal `.oudf` Template (when custom BCs are needed)

```okl
void udfDirichlet(bcData *bc)
{
  if (isField("fluid velocity")) {
    bc->uxFluid = 0.0;
    bc->uyFluid = 0.0;
    bc->uzFluid = 0.0;
  }
}
```
