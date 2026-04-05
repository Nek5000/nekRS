# nekRS UDF Patterns by Physics Type

## Required Headers

```cpp
#include "nrs.hpp"
#include "udf.hpp"
```

Additional headers by feature (use short form, matching existing examples):
- `#include "tavg.hpp"` — time averaging
- `#include "RANSktau.hpp"` — k-tau RANS
- `#include "lowMach.hpp"` — low-Mach
- `#include "lpm.hpp"` — Lagrangian particles

## Core Callback Signatures

```cpp
void UDF_Setup0(MPI_Comm comm, setupAide &options);
void UDF_Setup();
void UDF_LoadKernels(deviceKernelProperties& kernelInfo);
void UDF_ExecuteStep(double time, int tstep);
```

## Available User Callbacks (registered in UDF_Setup)

```cpp
nrs->userSource = &myForcing;           // Body forces
nrs->userProperties = &myProperties;    // Variable rho, mu
nrs->userDivergence = &myDivergence;    // Custom divergence (lowMach)
nrs->userConvergenceCheck = &myCheck;   // Outer iteration control
nrs->fluid->userImplicitLinearTerm = &myImplicit;  // Implicit forcing
```

## Key Data Access

```cpp
auto mesh = nrs->meshV;
dlong fieldOffset = nrs->fieldOffset;
dlong Nlocal = mesh->Nlocal;

// Velocity (device)
auto o_U = nrs->fluid->o_U;  // [ux | uy | uz], each fieldOffset long

// Scalars
auto& o_S = nrs->scalar->o_solution("temperature");

// Properties
auto& o_mue = nrs->fluid->o_mue;   // viscosity
auto& o_rho = nrs->fluid->o_rho;   // density
auto& o_diff = nrs->scalar->o_diffusionCoeff("temperature");

// Mesh coordinates (host)
auto [x, y, z] = mesh->xyzHost();

// Parameter extraction
dfloat val;
platform->par->extract("casedata", "key", val);
platform->options.getArgs("FLUID VISCOSITY", val);
```

## Linear Algebra Operations

```cpp
platform->linAlg->fill(N, value, o_array);
platform->linAlg->add(N, value, o_array);
platform->linAlg->scale(N, factor, o_array);
platform->linAlg->axpby(N, a, o_x, b, o_y);  // y = a*x + b*y
platform->linAlg->axpbyz(N, a, o_x, b, o_y, o_z);  // z = a*x + b*y
dfloat dot = platform->linAlg->innerProd(N, o_x, o_y, comm);
dfloat vmin = platform->linAlg->min(N, o_x, comm);
dfloat vmax = platform->linAlg->max(N, o_x, comm);
```

---

## Pattern: Laminar Flow

Minimal UDF. Boundary conditions handled in `.oudf` or via `zeroDirichlet`/`zeroNeumann`.

```cpp
#include "nrs.hpp"
#include "udf.hpp"

void UDF_Setup0(MPI_Comm comm, setupAide &options) {}
void UDF_Setup() {}
void UDF_ExecuteStep(double time, int tstep) {}
```

---

## Pattern: Constant Pressure-Gradient Driven Flow

For channel/pipe flows. Prefer `constFlowRate` in `.par` when possible.

**Manual forcing approach:**

```cpp
static dfloat P_DPDX;

#ifdef __okl__
@kernel void constantForcing(const dlong N,
                             const dlong offset,
                             const dfloat dpdx,
                             @restrict dfloat *EXT)
{
  for (dlong n = 0; n < N; ++n; @tile(p_blockSize, @outer, @inner)) {
    EXT[n + 0 * offset] += dpdx;  // x-direction force
  }
}
#endif

void userSource(double time)
{
  auto mesh = nrs->meshV;
  constantForcing(mesh->Nlocal, nrs->fieldOffset, P_DPDX, nrs->fluid->o_EXT);
}

void UDF_Setup0(MPI_Comm comm, setupAide &options)
{
  platform->par->extract("casedata", "p_dpdx", P_DPDX);
}

void UDF_Setup()
{
  nrs->userSource = &userSource;
}
```

---

## Pattern: Temperature-Dependent Properties

```cpp
static dfloat P_BETA;  // thermal expansion coefficient

#ifdef __okl__
@kernel void setProperties(const dlong N,
                           const dfloat rho0,
                           const dfloat mu0,
                           const dfloat beta,
                           @restrict const dfloat *TEMP,
                           @restrict dfloat *RHO,
                           @restrict dfloat *MUE)
{
  for (dlong n = 0; n < N; ++n; @tile(p_blockSize, @outer, @inner)) {
    RHO[n] = rho0 * (1.0 - beta * TEMP[n]);
    MUE[n] = mu0;
  }
}
#endif

void uservp(double time)
{
  auto mesh = nrs->meshV;
  dfloat rho0, mu0;
  platform->options.getArgs("FLUID DENSITY", rho0);
  platform->options.getArgs("FLUID VISCOSITY", mu0);

  setProperties(mesh->Nlocal, rho0, mu0, P_BETA,
                nrs->scalar->o_solution("temperature"),
                nrs->fluid->o_rho, nrs->fluid->o_mue);
}

void UDF_Setup()
{
  nrs->userProperties = &uservp;
}
```

---

## Pattern: Buoyancy (Boussinesq Approximation)

```cpp
static dfloat P_BETA, P_T_REF;
static int gravityDir = 1;  // y-direction

#ifdef __okl__
@kernel void buoyancyForcing(const dlong N,
                             const dlong offset,
                             const dfloat beta,
                             const dfloat T_ref,
                             const int gDir,
                             @restrict const dfloat *TEMP,
                             @restrict dfloat *EXT)
{
  for (dlong n = 0; n < N; ++n; @tile(p_blockSize, @outer, @inner)) {
    const dfloat dT = TEMP[n] - T_ref;
    EXT[n + gDir * offset] += -beta * dT;  // -beta*(T-Tref)*g_hat
  }
}
#endif

void userSource(double time)
{
  auto mesh = nrs->meshV;
  buoyancyForcing(mesh->Nlocal, nrs->fieldOffset, P_BETA, P_T_REF, gravityDir,
                  nrs->scalar->o_solution("temperature"), nrs->fluid->o_EXT);
}
```

---

## Pattern: RANS k-tau SST

Requires scalars `K` and `TAU` in `.par`, and `equation = navierStokes+variableViscosity`
in `[PROBLEMTYPE]`.

### Minimal RANS (no temperature)

```cpp
#include "RANSktau.hpp"

#ifdef __okl__

void udfDirichlet(bcData *bc)
{
  // Wall BC for k and tau
  if (isField("scalar k")) {
    bc->sScalar = 0.0;
  }
  if (isField("scalar tau")) {
    bc->sScalar = 1e-12;
  }
}

#endif

void userq(double time)
{
  RANSktau::updateSourceTerms();
}

void uservp(double time)
{
  RANSktau::updateProperties();
}

void UDF_Setup0(MPI_Comm comm, setupAide &options) {}

void UDF_Setup()
{
  nrs->userProperties = &uservp;
  nrs->userSource = &userq;

  // Set IC if not restarting
  if (platform->options.getArgs("RESTART FILE NAME").empty()) {
    auto mesh = nrs->meshV;
    std::vector<dfloat> U(nrs->fluid->fieldOffsetSum, 0.0);
    std::vector<dfloat> k(mesh->Nlocal, 0.01);
    std::vector<dfloat> tau(mesh->Nlocal, 0.1);
    for (int n = 0; n < mesh->Nlocal; n++) {
      U[n + 0 * nrs->fieldOffset] = 1.0;
    }
    nrs->fluid->o_U.copyFrom(U.data(), U.size());
    nrs->scalar->o_solution("k").copyFrom(k.data(), k.size());
    nrs->scalar->o_solution("tau").copyFrom(tau.data(), tau.size());
  }

  RANSktau::setup(nrs->scalar->nameToIndex.find("k")->second);
}

void UDF_ExecuteStep(double time, int tstep) {}
```

### RANS with Inlet Turbulence (inflow/outflow cases)

For cases with inflow BCs, pass inlet k/tau as compile-time constants:

```cpp
#include "RANSktau.hpp"

static dfloat P_U_INLET, P_K_INLET, P_TAU_INLET;

#ifdef __okl__

void udfDirichlet(bcData *bc)
{
  if (isField("fluid velocity")) {
    if (bc->id == 1) {  // inlet
      bc->uxFluid = p_U_INLET;
      bc->uyFluid = 0.0;
      bc->uzFluid = 0.0;
    }
  }
  if (isField("scalar k")) {
    if (bc->id == 1) {
      bc->sScalar = p_K_INLET;
    } else if (bc->id == 3) {  // wall
      bc->sScalar = 0.0;
    }
  }
  if (isField("scalar tau")) {
    if (bc->id == 1) {
      bc->sScalar = p_TAU_INLET;
    } else if (bc->id == 3) {  // wall
      bc->sScalar = 1e-12;
    }
  }
}

#endif

void userq(double time)
{
  RANSktau::updateSourceTerms();
}

void uservp(double time)
{
  RANSktau::updateProperties();
}

void UDF_LoadKernels(deviceKernelProperties &kernelInfo)
{
  kernelInfo.define("p_U_INLET") = P_U_INLET;
  kernelInfo.define("p_K_INLET") = P_K_INLET;
  kernelInfo.define("p_TAU_INLET") = P_TAU_INLET;
}

void UDF_Setup0(MPI_Comm comm, setupAide &options)
{
  dfloat TI;
  platform->par->extract("casedata", "u_inlet", P_U_INLET);
  platform->par->extract("casedata", "ti", TI);

  // k = 1.5 * (U * TI)^2
  P_K_INLET = 1.5 * (P_U_INLET * TI) * (P_U_INLET * TI);

  // tau = 1/omega, omega = k^0.5 / (Cmu^0.25 * l_t), l_t = 0.07 * D
  dfloat D;
  platform->par->extract("casedata", "diameter", D);
  const dfloat l_t = 0.07 * D;
  const dfloat Cmu = 0.09;
  const dfloat omega = sqrt(P_K_INLET) / (pow(Cmu, 0.25) * l_t);
  P_TAU_INLET = 1.0 / omega;
}

void UDF_Setup()
{
  nrs->userProperties = &uservp;
  nrs->userSource = &userq;

  if (platform->options.getArgs("RESTART FILE NAME").empty()) {
    auto mesh = nrs->meshV;
    std::vector<dfloat> U(nrs->fluid->fieldOffsetSum, 0.0);
    std::vector<dfloat> k(mesh->Nlocal, P_K_INLET);
    std::vector<dfloat> tau(mesh->Nlocal, P_TAU_INLET);
    for (int n = 0; n < mesh->Nlocal; n++) {
      U[n + 0 * nrs->fieldOffset] = P_U_INLET;
    }
    nrs->fluid->o_U.copyFrom(U.data(), U.size());
    nrs->scalar->o_solution("k").copyFrom(k.data(), k.size());
    nrs->scalar->o_solution("tau").copyFrom(tau.data(), tau.size());
  }

  RANSktau::setup(nrs->scalar->nameToIndex.find("k")->second);
}

void UDF_ExecuteStep(double time, int tstep) {}
```

### RANS with Temperature (turbulent heat transfer)

Add turbulent diffusion to the temperature scalar's diffusion coefficient:

```cpp
// In uservp, AFTER RANSktau::updateProperties():
dfloat conduct;
platform->options.getArgs("SCALAR00 DIFFUSIVITY", conduct);
const dfloat Pr_t = 0.85;
// Need a scalarScaledAdd kernel for this:
// Y[n] = conduct + (1/Pr_t) * mue_t[n]
scalarScaledAdd(mesh->Nlocal, conduct, 1.0 / Pr_t,
                RANSktau::o_mue_t(),
                nrs->scalar->o_diffusionCoeff("temperature"));
```

Only include the `scalarScaledAdd` OKL kernel when the case has a temperature scalar.

### RANS .par template

```ini
[GENERAL]
polynomialOrder = 7
dt = targetCFL=0.5+initial=1e-3
timeStepper = tombo2
scalars = k, tau

[PROBLEMTYPE]
equation = navierStokes+variableViscosity

[SCALAR K]
boundaryTypeMap = udfDirichlet, zeroNeumann, udfDirichlet
residualTol = 1e-08

[SCALAR TAU]
boundaryTypeMap = udfDirichlet, zeroNeumann, udfDirichlet
residualTol = 1e-08
```
```

---

## Pattern: Time Averaging (tavg)

```cpp
#include "tavg.hpp"

static std::unique_ptr<tavg> avg;

void UDF_Setup()
{
  auto mesh = nrs->meshV;
  dlong fieldOffset = nrs->fieldOffset;

  std::vector<tavg::field> tavgFields;
  deviceMemory<dfloat> o_ux(nrs->fluid->o_solution("x"));
  deviceMemory<dfloat> o_uy(nrs->fluid->o_solution("y"));
  deviceMemory<dfloat> o_uz(nrs->fluid->o_solution("z"));

  // Mean fields
  tavgFields.push_back({"ux", std::vector{o_ux}});
  tavgFields.push_back({"uy", std::vector{o_uy}});
  tavgFields.push_back({"uz", std::vector{o_uz}});

  // Reynolds stresses
  tavgFields.push_back({"uxux", std::vector{o_ux, o_ux}});
  tavgFields.push_back({"uyuy", std::vector{o_uy, o_uy}});
  tavgFields.push_back({"uzuz", std::vector{o_uz, o_uz}});
  tavgFields.push_back({"uxuy", std::vector{o_ux, o_uy}});

  avg = std::make_unique<tavg>(fieldOffset, tavgFields);
}

void UDF_ExecuteStep(double time, int tstep)
{
  if (nrs->timeStepConverged) avg->run(time);
  if (nrs->checkpointStep) avg->writeToFile(nrs->meshV);
}
```

---

## Pattern: Low-Mach Variable Density

```cpp
#include "lowMach.hpp"

static dfloat P_GAMMA;
static deviceMemory<dfloat> o_beta, o_kappa;

void qtl(double time)
{
  lowMach::qThermalSingleComponent(time);
}

void uservp(double time)
{
  // Update density from temperature: rho = p0 / (R * T)
  // Fill o_beta and o_kappa based on temperature field
}

void UDF_Setup()
{
  dfloat alphaRef = (P_GAMMA - 1.0) / P_GAMMA;
  lowMach::setup(alphaRef, o_beta, o_kappa);
  nrs->userDivergence = &qtl;
  nrs->userProperties = &uservp;
}
```

---

## Pattern: Custom Initial Conditions

Set IC in `UDF_Setup()` after solvers are initialized. Use `fieldOffsetSum` for velocity:

```cpp
void UDF_Setup()
{
  auto mesh = nrs->meshV;
  auto [x, y, z] = mesh->xyzHost();

  if (platform->options.getArgs("RESTART FILE NAME").empty()) {
    std::vector<dfloat> U(nrs->fluid->fieldOffsetSum, 0.0);
    for (dlong n = 0; n < mesh->Nlocal; n++) {
      U[n + 0 * nrs->fieldOffset] = 1.0;  // ux
      U[n + 1 * nrs->fieldOffset] = 0.0;  // uy
      U[n + 2 * nrs->fieldOffset] = 0.0;  // uz
    }
    nrs->fluid->o_U.copyFrom(U.data(), U.size());
  }
}
```

---

## Pattern: Checkpoint Custom Fields

```cpp
void UDF_Setup()
{
  // Add Q-criterion as checkpoint field
  deviceMemory<dfloat> o_Qcriterion(nrs->fieldOffset);
  nrs->addUserCheckpointField("Q", o_Qcriterion);
}
```
