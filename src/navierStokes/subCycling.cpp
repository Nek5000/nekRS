#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "linAlg.hpp"
#include "nrs.hpp"

namespace {

static oogs_t *gsh = nullptr;
static mesh_t *mesh = nullptr;
static mesh_t *meshV = nullptr;
static occa::kernel opKernel;
static dlong cubatureOffset;
static dlong fieldOffset;
static dlong meshOffset;

static void flops(mesh_t *mesh, int Nfields)
{
  const auto cubNq = mesh->cubNq;
  const auto cubNp = mesh->cubNp;
  const auto Nq = mesh->Nq;
  const auto Np = mesh->Np;
  const auto nEXT = 3;
  const auto Nelements = mesh->Nelements;
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

static dfloat *coeffs(nrs_t *nrs, int nEXT, double time, dfloat tstage, dfloat sdt, int rk)
{
  const double t = tstage + sdt * nrs->nodesRK[rk];
  const double tn0 = time;
  const double tn1 = time - nrs->dt[1];
  const double tn2 = time - (nrs->dt[1] + nrs->dt[2]);
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

static void applyOperator(int nFields,
                          dfloat *extC,
                          occa::memory &o_Urst,
                          occa::memory &o_u1,
                          occa::memory &o_rhs)
{
  if (meshV->NglobalGatherElements) {
    if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
      opKernel(meshV->NglobalGatherElements,
             meshV->o_globalGatherElementList,
             mesh->o_cubDiffInterpT,
             mesh->o_cubInterpT,
             fieldOffset,
             cubatureOffset,
             meshOffset,
             mesh->o_invLMM,
             mesh->o_divU,
             extC[0],
             extC[1],
             extC[2],
             o_Urst,
             o_u1,
             o_rhs);
    } else {
      opKernel(meshV->NglobalGatherElements,
             meshV->o_globalGatherElementList,
             mesh->o_D,
             fieldOffset,
             meshOffset,
             mesh->o_invLMM,
             mesh->o_divU,
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
             mesh->o_cubDiffInterpT,
             mesh->o_cubInterpT,
             fieldOffset,
             cubatureOffset,
             meshOffset,
             mesh->o_invLMM,
             mesh->o_divU,
             extC[0],
             extC[1],
             extC[2],
             o_Urst,
             o_u1,
             o_rhs);
    } else {
      opKernel(meshV->NlocalGatherElements,
             meshV->o_localGatherElementList,
             mesh->o_D,
             fieldOffset,
             meshOffset,
             mesh->o_invLMM,
             mesh->o_divU,
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

static void rk(nrs_t *nrs,
               int nFields,
               int nEXT,
               double time,
               dfloat tstage,
               dfloat sdt,
               occa::memory &o_Urst,
               occa::memory &o_u0)
{
  linAlg_t *linAlg = platform->linAlg;

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  occa::memory o_u1 = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);
  o_u1.copyFrom(o_u0);

  std::vector<occa::memory> o_rhs(4);
  o_rhs[0] = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[1] = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[2] = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);
  o_rhs[3] = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);

  occa::memory o_LMMe = (movingMesh) ? platform->o_memPool.reserve<dfloat>(fieldOffset) : nullptr;

  for (int rk = 0; rk < nrs->nRK; ++rk) {
      auto extC = coeffs(nrs, nEXT, time, tstage, sdt, rk);

      if (movingMesh) {
        nrs->nStagesSum3Kernel(meshV->Nlocal,
                               meshOffset,
                               nEXT,
                               extC[0],
                               extC[1],
                               extC[2],
                               mesh->o_LMM,
                               o_LMMe);
        linAlg->aydxMany(meshV->Nlocal, nFields, fieldOffset, 0, 1.0, o_LMMe, o_u1);
      }

      applyOperator(nFields, extC, o_Urst, o_u1, o_rhs[rk]);

      if (movingMesh) {
        linAlg->axmyMany(meshV->Nlocal, nFields, fieldOffset, 0, 1.0, o_LMMe, o_rhs[rk]);
      }

      if (rk != nrs->nRK - 1) {
        linAlg->axpbyzMany(meshV->Nlocal,
                           nFields,
                           fieldOffset,
                           1.0,
                           o_u0,
                           -sdt * nrs->coeffsfRK[rk + 1],
                           o_rhs[rk],
                           o_u1);
      } else {
        nrs->subCycleRKKernel(meshV->Nlocal,
                              nFields,
                              fieldOffset,
                              sdt,
                              nrs->o_weightsRK,
                              o_rhs[0],
                              o_rhs[1],
                              o_rhs[2],
                              o_rhs[3],
                              o_u0);
      }
  }
}

static occa::memory _subCycle(nrs_t *nrs, mesh_t *_mesh, oogs_t *_gsh, occa::kernel _opKernel,
                              int nEXT, int nFields, 
                              dlong _meshOffset, dlong _fieldOffset, dlong fieldOffsetSum,
                              double time, occa::memory o_Urst, occa::memory o_U)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  mesh = _mesh;
  gsh =  _gsh;
  opKernel = _opKernel;
  meshOffset = _meshOffset;
  fieldOffset = _fieldOffset;

  meshV = nrs->meshV;
  cubatureOffset = nrs->cubatureOffset; // fine grid convection velocity

  occa::memory o_u0 = platform->o_memPool.reserve<dfloat>(nFields * fieldOffset);
  platform->linAlg->fill(o_u0.length(), 0.0, o_u0);

  for (int torder = nEXT - 1; torder >= 0; torder--) {
     // Initialize u0 = U^(t-torder*dt)
     nrs->subCycleInitU0Kernel(mesh->Nlocal,
                               nFields,
                               fieldOffset,
                               torder,
                               nEXT,
                               torder * fieldOffsetSum,
                               (movingMesh) ? torder * meshOffset : 0, // offset to lagged LMM 
                               nrs->coeffBDF[torder],
                               mesh->o_LMM,
                               o_U,
                               o_u0);

     // Advance sub-problem from t^(n-torder) to t^(n-torder+1)
     const double dtSubStep = nrs->dt[torder] / nrs->Nsubsteps;

     auto t0 = [&](int step) {
       double sum = 0;
       for (int i = torder; i > 0; i--) {
         sum += nrs->dt[i];
       }
       return (time - sum) + step * dtSubStep;
     }; 

     for (int tSubStep = 0; tSubStep < nrs->Nsubsteps; ++tSubStep) {
       rk(nrs,
          nFields, 
          nEXT, 
          time,
          t0(dtSubStep), 
          dtSubStep, 
          o_Urst, 
          o_u0);
     }
  }

  if (!movingMesh) {
     platform->linAlg->axmyMany(mesh->Nlocal, 
                      nFields, 
                      fieldOffset, 
                      0,
                      1.0, 
                      mesh->o_LMM, 
                      o_u0);
  }

  return o_u0;
}

} // namespace

occa::memory velocitySubCycle(nrs_t *nrs, int nEXT, double time)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const auto mesh = nrs->meshV;
  const auto gsh = nrs->gsh;

  const auto nFields = nrs->NVfields;
  const auto fieldOffset = nrs->fieldOffset;
  const auto fieldOffsetSum = nFields * fieldOffset;
  const auto meshOffset = nrs->fieldOffset; // used for mesh->o_invLMM and mesh->o_divU in case of moving mesh

  auto o_U = nrs->o_U;

  auto o_Urst = (movingMesh) ? nrs->o_relUrst : nrs->o_Urst;
  auto opKernel = (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) ? 
                  nrs->subCycleStrongCubatureVolumeKernel :
                  nrs->subCycleStrongVolumeKernel;

  return _subCycle(nrs, mesh, gsh, opKernel, nEXT, nFields, meshOffset, fieldOffset, fieldOffsetSum, time, o_Urst, o_U);
}

occa::memory scalarSubCycle(nrs_t *nrs, int nEXT, double time, int scalarIdx)
{
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  const auto mesh = (scalarIdx == 0) ? nrs->cds->mesh[0] : nrs->cds->meshV;
  const auto gsh = nrs->cds->gsh;

  const auto nFields = 1;
  const auto fieldOffset = nrs->cds->fieldOffset[scalarIdx]; 
  const auto fieldOffsetSum = nrs->cds->fieldOffsetSum;
  const auto meshOffset = nrs->fieldOffset; // used for mesh->o_invLMM and mesh->o_divU in case of moving mesh

  auto o_U = nrs->cds->o_S.slice(nrs->cds->fieldOffsetScan[scalarIdx], 
                                fieldOffset);

  auto o_Urst = (movingMesh) ? nrs->cds->o_relUrst : nrs->cds->o_Urst;
  auto opKernel = (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE")) ?
                  nrs->cds->subCycleStrongCubatureVolumeKernel :  
                  nrs->cds->subCycleStrongVolumeKernel;

  return _subCycle(nrs, mesh, gsh, opKernel, nEXT, nFields, meshOffset, fieldOffset, fieldOffsetSum, time, o_Urst, o_U);
}
