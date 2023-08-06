#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "linAlg.hpp"
#include "nrs.hpp"

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
    flopCount += 6. * cubNp * nEXT;  // extrapolate U(r,s,t) to current time
    flopCount += 6. * cubNp * cubNq * Nfields;                                       // apply Dcub
    flopCount += 3. * Np * Nfields;                                                  // compute NU
    flopCount += 4. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq) * Nfields; // interpolation
  }
  else {
    flopCount = Nq * Nq * Nq * (6. * Nq + 6. * nEXT + 8.) * Nfields;
  }
  flopCount *= Nelements;

  platform->flopCounter->add("subcycling", flopCount);
}

occa::memory
scalarSubCycleMovingMesh(cds_t *cds, int nEXT, dfloat time, int is, occa::memory o_U, occa::memory o_S)
{
  std::string sid = scalarDigitStr(is);

  linAlg_t *linAlg = platform->linAlg;

  occa::memory &o_r1 = platform->o_mempool.slice2;
  occa::memory &o_r2 = platform->o_mempool.slice3;
  occa::memory &o_r3 = platform->o_mempool.slice4;
  occa::memory &o_r4 = platform->o_mempool.slice5;

  occa::memory &o_p0 = platform->o_mempool.slice0;
  occa::memory &o_u1 = platform->o_mempool.slice6;

  occa::memory &o_LMMe = platform->o_mempool.slice1;

  // Solve for Each SubProblem
  for (int torder = (nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset =
        cds->fieldOffsetScan[is] + torder * cds->fieldOffsetSum;
    const dlong offset = torder * cds->fieldOffset[is];
    cds->subCycleInitU0Kernel(cds->mesh[0]->Nlocal,
        1,
        cds->fieldOffset[is],
        torder,
        nEXT,
        toffset,
        offset,
        cds->coeffBDF[torder],
        cds->mesh[0]->o_LMM,
        o_S,
        o_p0);

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0; i--)
      tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder] / cds->Nsubsteps;

    for (int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;
      o_u1.copyFrom(o_p0, cds->mesh[0]->Nlocal * sizeof(dfloat));
      for (int rk = 0; rk < cds->nRK; ++rk) {
        occa::memory o_rhs;
        if (rk == 0)
          o_rhs = o_r1;
        if (rk == 1)
          o_rhs = o_r2;
        if (rk == 2)
          o_rhs = o_r3;
        if (rk == 3)
          o_rhs = o_r4;

        // Extrapolate velocity to subProblem stage time
        const dfloat t = tstage + sdt * cds->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        dfloat extC[3] = {0., 0., 0.};
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
        cds->nStagesSum3Kernel(cds->mesh[0]->Nlocal,
            cds->vFieldOffset,
            nEXT,
            extC[0],
            extC[1],
            extC[2],
            cds->mesh[0]->o_LMM,
            o_LMMe);
        linAlg->aydx(cds->mesh[0]->Nlocal, 1.0, o_LMMe, o_u1);

        if (cds->meshV->NglobalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(cds->meshV->NglobalGatherElements,
                                                    cds->meshV->o_globalGatherElementList,
                                                    cds->meshV->o_cubDiffInterpT,
                                                    cds->meshV->o_cubInterpT,
                                                    cds->vFieldOffset,
                                                    cds->vCubatureOffset,
                                                    0,
                                                    cds->mesh[0]->o_invLMM,
                                                    cds->mesh[0]->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    cds->o_relUrst,
                                                    o_u1,
                                                    o_rhs);
          else
            cds->subCycleStrongVolumeKernel(cds->meshV->NglobalGatherElements,
                cds->meshV->o_globalGatherElementList,
                cds->meshV->o_D,
                cds->vFieldOffset,
                0,
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_relUrst,
                o_u1,
                o_rhs);
        }

        oogs::start(
            o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        if (cds->meshV->NlocalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(cds->meshV->NlocalGatherElements,
                                                    cds->meshV->o_localGatherElementList,
                                                    cds->meshV->o_cubDiffInterpT,
                                                    cds->meshV->o_cubInterpT,
                                                    cds->vFieldOffset,
                                                    cds->vCubatureOffset,
                                                    0,
                                                    cds->mesh[0]->o_invLMM,
                                                    cds->mesh[0]->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    cds->o_relUrst,
                                                    o_u1,
                                                    o_rhs);
          else
            cds->subCycleStrongVolumeKernel(cds->meshV->NlocalGatherElements,
                cds->meshV->o_localGatherElementList,
                cds->meshV->o_D,
                cds->vFieldOffset,
                0,
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_relUrst,
                o_u1,
                o_rhs);
        }

        oogs::finish(
            o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        flops(cds->mesh[0], 1);

        linAlg->axmy(cds->mesh[0]->Nlocal, 1.0, o_LMMe, o_rhs);
        if (rk != 3)
          linAlg->axpbyz(cds->mesh[0]->Nlocal,
              1.0,
              o_p0,
              -sdt * cds->coeffsfRK[rk + 1],
              o_rhs,
              o_u1);
        else
          cds->subCycleRKKernel(cds->mesh[0]->Nlocal,
              sdt,
              cds->o_weightsRK,
              o_r1,
              o_r2,
              o_r3,
              o_r4,
              o_p0);
      }
    }
  }
  return o_p0;
}

occa::memory scalarSubCycle(cds_t *cds, int nEXT, dfloat time, int is, occa::memory o_U, occa::memory o_S)
{
  std::string sid = scalarDigitStr(is);

  linAlg_t *linAlg = platform->linAlg;

  // Solve for Each SubProblem
  for (int torder = (nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset =
        cds->fieldOffsetScan[is] + torder * cds->fieldOffsetSum;
    cds->subCycleInitU0Kernel(cds->mesh[0]->Nlocal,
        1,
        cds->fieldOffset[is],
        torder,
        nEXT,
        toffset,
        0,
        cds->coeffBDF[torder],
        cds->mesh[0]->o_LMM,
        o_S,
        platform->o_mempool.slice0);

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0; i--)
      tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder] / cds->Nsubsteps;

    for (int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_mempool.slice0.copyFrom(platform->o_mempool.slice0,
          cds->fieldOffset[is] * sizeof(dfloat),
          cds->fieldOffset[is] * sizeof(dfloat),
          0);

      for (int rk = 0; rk < cds->nRK; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t = tstage + sdt * cds->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        dfloat extC[3] = {0., 0., 0.};
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

        if (cds->meshV->NglobalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
                cds->meshV->NglobalGatherElements,
                cds->meshV->o_globalGatherElementList,
                cds->meshV->o_cubDiffInterpT,
                cds->meshV->o_cubInterpT,
                cds->vFieldOffset,
                cds->vCubatureOffset,
                rk * cds->fieldOffset[is],
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice2);
          else
            cds->subCycleStrongVolumeKernel(cds->meshV->NglobalGatherElements,
                cds->meshV->o_globalGatherElementList,
                cds->meshV->o_D,
                cds->vFieldOffset,
                rk * cds->fieldOffset[is],
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice2);
        }

        occa::memory o_rhs;
        if (rk == 0)
          o_rhs = platform->o_mempool.slice2;
        if (rk == 1)
          o_rhs = platform->o_mempool.slice3;
        if (rk == 2)
          o_rhs = platform->o_mempool.slice4;
        if (rk == 3)
          o_rhs = platform->o_mempool.slice5;

        oogs::start(
            o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        if (cds->meshV->NlocalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
                cds->meshV->NlocalGatherElements,
                cds->meshV->o_localGatherElementList,
                cds->meshV->o_cubDiffInterpT,
                cds->meshV->o_cubInterpT,
                cds->vFieldOffset,
                cds->vCubatureOffset,
                rk * cds->fieldOffset[is],
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice2);
          else
            cds->subCycleStrongVolumeKernel(cds->meshV->NlocalGatherElements,
                cds->meshV->o_localGatherElementList,
                cds->meshV->o_D,
                cds->vFieldOffset,
                rk * cds->fieldOffset[is],
                cds->mesh[0]->o_invLMM,
                cds->mesh[0]->o_divU,
                extC[0],
                extC[1],
                extC[2],
                cds->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice2);
        }

        oogs::finish(
            o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        flops(cds->mesh[0], 1);

        cds->subCycleRKUpdateKernel(cds->meshV->Nlocal,
            rk,
            sdt,
            cds->fieldOffset[is],
            cds->o_coeffsfRK,
            cds->o_weightsRK,
            platform->o_mempool.slice1,
            platform->o_mempool.slice2,
            platform->o_mempool.slice0);
      }
    }
  }
  linAlg->axmy(cds->mesh[0]->Nlocal,
      1.0,
      cds->mesh[0]->o_LMM,
      platform->o_mempool.slice0);
  return platform->o_mempool.slice0;
}
