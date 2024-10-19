#include "platform.hpp"
#include "linAlg.hpp"

namespace
{

static oogs_t *gsh = nullptr;
static mesh_t *mesh = nullptr;
static mesh_t *meshV = nullptr;
static occa::kernel opKernel;
static dlong cubatureOffset;
static dlong fieldOffset;
static dlong meshOffset;

static occa::kernel subCycleInitU0Kernel;
static occa::kernel nStagesSum3Kernel;
static occa::kernel subCycleRKKernel;

static void flops(mesh_t *mesh, int Nfields)
{
  const auto cubNq = meshV->cubNq;
  const auto cubNp = meshV->cubNp;
  const auto Nq = meshV->Nq;
  const auto Np = meshV->Np;
  const auto nEXT = 3;
  const auto Nelements = meshV->Nelements;
  double flopCount = 0.0; // per elem basis
  if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
    flopCount += 6. * cubNp * nEXT;            // extrapolate U(r,s,t) to current time
    flopCount += 6. * cubNp * cubNq * Nfields; // apply Dcub
    flopCount += 3. * Np * Nfields;            // compute NU
    flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq) * Nfields; // interpolation
  } else {
    flopCount = Nq * Nq * Nq * (6. * Nq + 6. * nEXT + 8.) * Nfields;
  }
  flopCount *= Nelements;

  platform->flopCounter->add("subcycling", flopCount);
}

static dfloat *
extCoeffs(int nEXT, double time, dfloat tstage, dfloat sdt, dfloat *dt, dfloat *nodesRK, int rk)
{
  const double t = tstage + sdt * nodesRK[rk];
  const double tn0 = time;
  const double tn1 = time - dt[1];
  const double tn2 = time - (dt[1] + dt[2]);
  static dfloat extC[3] = {0, 0, 0};
  switch (nEXT) {
  case 1:
    extC[0] = 1;
    extC[1] = 0;
    extC[2] = 0;
    break;
  case 2:
    extC[0] = (t - tn1) / (tn0 - tn1);
    extC[1] = (t - tn0) / (tn1 - tn0);
    extC[2] = 0;
    break;
  case 3:
    extC[0] = (t - tn1) * (t - tn2) / ((tn0 - tn1) * (tn0 - tn2));
    extC[1] = (t - tn0) * (t - tn2) / ((tn1 - tn0) * (tn1 - tn2));
    extC[2] = (t - tn0) * (t - tn1) / ((tn2 - tn0) * (tn2 - tn1));
    break;
  }
  return extC;
}

static void
applyOperator(int nFields, dfloat *extC, occa::memory &o_Urst, occa::memory &o_u1, occa::memory &o_rhs)
{
  if (meshV->NglobalGatherElements) {
    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      opKernel(meshV->NglobalGatherElements,
               meshV->o_globalGatherElementList,
               meshV->o_cubDiffInterpT,
               meshV->o_cubInterpT,
               fieldOffset,
               cubatureOffset,
               meshOffset,
               mesh->o_invLMM,
               meshV->o_divU,
               extC[0],
               extC[1],
               extC[2],
               o_Urst,
               o_u1,
               o_rhs);
    } else {
      opKernel(meshV->NglobalGatherElements,
               meshV->o_globalGatherElementList,
               meshV->o_D,
               fieldOffset,
               meshOffset,
               mesh->o_invLMM,
               meshV->o_divU,
               extC[0],
               extC[1],
               extC[2],
               o_Urst,
               o_u1,
               o_rhs);
    }
  }

  oogs::start(o_rhs, nFields, fieldOffset, ogsDfloat, ogsAdd, gsh);

  if (meshV->NlocalGatherElements) {
    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      opKernel(meshV->NlocalGatherElements,
               meshV->o_localGatherElementList,
               meshV->o_cubDiffInterpT,
               meshV->o_cubInterpT,
               fieldOffset,
               cubatureOffset,
               meshOffset,
               mesh->o_invLMM,
               meshV->o_divU,
               extC[0],
               extC[1],
               extC[2],
               o_Urst,
               o_u1,
               o_rhs);
    } else {
      opKernel(meshV->NlocalGatherElements,
               meshV->o_localGatherElementList,
               meshV->o_D,
               fieldOffset,
               meshOffset,
               mesh->o_invLMM,
               meshV->o_divU,
               extC[0],
               extC[1],
               extC[2],
               o_Urst,
               o_u1,
               o_rhs);
    }
  }

  oogs::finish(o_rhs, nFields, fieldOffset, ogsDfloat, ogsAdd, gsh);

  flops(meshV, nFields);
}

static void rk44(int nFields,
                 int nEXT,
                 double time,
                 dfloat tstage,
                 dfloat sdt,
                 dfloat *dt,
                 occa::memory &o_Urst,
                 occa::memory &o_u0)
{
  constexpr int nRK = 4;
  dfloat nodes[nRK] = {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
  dfloat weights[nRK] = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
  static occa::memory o_weights;
  if (!o_weights.isInitialized()) {
    o_weights = platform->device.malloc<dfloat>(nRK);
    o_weights.copyFrom(weights);
  }

  linAlg_t *linAlg = platform->linAlg;

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  occa::memory o_u1 = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);
  o_u1.copyFrom(o_u0);

  std::vector<occa::memory> o_rhs(4);
  o_rhs[0] = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[1] = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[2] = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[3] = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);

  occa::memory o_LMMe = (movingMesh) ? platform->deviceMemoryPool.reserve<dfloat>(fieldOffset) : nullptr;

  for (int rk = 0; rk < nRK; ++rk) {
    auto extC = extCoeffs(nEXT, time, tstage, sdt, dt, nodes, rk);

    if (movingMesh) {
      nStagesSum3Kernel(meshV->Nlocal, meshOffset, nEXT, extC[0], extC[1], extC[2], meshV->o_LMM, o_LMMe);
      linAlg->aydxMany(meshV->Nlocal, nFields, fieldOffset, 0, 1.0, o_LMMe, o_u1);
    }

    applyOperator(nFields, extC, o_Urst, o_u1, o_rhs[rk]);

    if (movingMesh) {
      linAlg->axmyMany(meshV->Nlocal, nFields, fieldOffset, 0, 1.0, o_LMMe, o_rhs[rk]);
    }

    if (rk != nRK - 1) {
      linAlg
          ->axpbyzMany(meshV->Nlocal, nFields, fieldOffset, 1.0, o_u0, -sdt * nodes[rk + 1], o_rhs[rk], o_u1);
    } else {
      subCycleRKKernel(meshV->Nlocal,
                       nFields,
                       fieldOffset,
                       sdt,
                       o_weights,
                       o_rhs[0],
                       o_rhs[1],
                       o_rhs[2],
                       o_rhs[3],
                       o_u0);
    }
  }
}

} // namespace

occa::memory advectionSubcyclingRK(mesh_t *_meshT,
                                   mesh_t *_meshV,
                                   double time,
                                   dfloat *dt,
                                   int Nsubsteps,
                                   dfloat *coeffBDF,
                                   int nEXT,
                                   int nFields,
                                   occa::kernel _opKernel,
                                   oogs_t *_gsh,
                                   dlong _meshOffset,
                                   dlong _fieldOffset,
                                   dlong _cubatureOffset,
                                   dlong fieldOffsetSum,
                                   occa::memory o_Urst,
                                   occa::memory o_U)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  mesh = _meshT;
  meshV = _meshV;

  gsh = _gsh;
  opKernel = _opKernel;
  meshOffset = _meshOffset;
  fieldOffset = _fieldOffset;
  cubatureOffset = _cubatureOffset;

  subCycleInitU0Kernel = platform->kernelRequests.load("nrs-subCycleInitU0");
  nStagesSum3Kernel = platform->kernelRequests.load("core-nStagesSum3");
  subCycleRKKernel = platform->kernelRequests.load("nrs-subCycleRK");

  occa::memory o_u0 = platform->deviceMemoryPool.reserve<dfloat>(nFields * fieldOffset);
  platform->linAlg->fill(o_u0.length(), 0.0, o_u0);

  for (int torder = nEXT - 1; torder >= 0; torder--) {
    // Initialize u0 = U^(t-torder*dt)
    subCycleInitU0Kernel(mesh->Nlocal,
                         nFields,
                         fieldOffset,
                         torder,
                         nEXT,
                         torder * fieldOffsetSum,
                         (movingMesh) ? torder * meshOffset : 0, // offset to lagged LMM
                         coeffBDF[torder],
                         mesh->o_LMM,
                         o_U,
                         o_u0);

    // Advance sub-problem from t^(n-torder) to t^(n-torder+1)
    const double dtSubStep = dt[torder] / Nsubsteps;

    auto t0 = [&](int step) {
      double sum = 0;
      for (int i = torder; i > 0; i--) {
        sum += dt[i];
      }
      return (time - sum) + step * dtSubStep;
    };

    for (int tSubStep = 0; tSubStep < Nsubsteps; ++tSubStep) {
      rk44(nFields, nEXT, time, t0(dtSubStep), dtSubStep, dt, o_Urst, o_u0);
    }
  }

  if (!movingMesh) {
    platform->linAlg->axmyMany(mesh->Nlocal, nFields, fieldOffset, 0, 1.0, mesh->o_LMM, o_u0);
  }

  return o_u0;
}
