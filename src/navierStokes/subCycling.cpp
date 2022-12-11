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

occa::memory velocitySubCycleMovingMesh(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  occa::memory &o_p0 = platform->o_mempool.slice0;
  occa::memory &o_u1 = platform->o_mempool.slice3;

  occa::memory &o_r1 = platform->o_mempool.slice6;
  occa::memory &o_r2 = platform->o_mempool.slice9;
  occa::memory &o_r3 = platform->o_mempool.slice12;
  occa::memory &o_r4 = platform->o_mempool.slice15;

  occa::memory &o_LMMe = platform->o_mempool.slice18;

  // Solve for Each SubProblem
  for (int torder = nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    const dlong offset = torder * nrs->fieldOffset;
    nrs->subCycleInitU0Kernel(mesh->Nlocal,
        nrs->NVfields,
        nrs->fieldOffset,
        torder,
        nEXT,
        toffset,
        offset,
        nrs->coeffBDF[torder],
        mesh->o_LMM,
        o_U,
        o_p0);

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0; i--)
      tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder] / nrs->Nsubsteps;

    for (int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      o_u1.copyFrom(o_p0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

      for (int rk = 0; rk < nrs->nRK; ++rk) {
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
        const dfloat t = tstage + sdt * nrs->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
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

        nrs->nStagesSum3Kernel(mesh->Nlocal,
            nrs->fieldOffset,
            nEXT,
            extC[0],
            extC[1],
            extC[2],
            mesh->o_LMM,
            o_LMMe);
        linAlg->aydxMany(mesh->Nlocal,
            nrs->NVfields,
            nrs->fieldOffset,
            0,
            1.0,
            o_LMMe,
            o_u1);

        if (mesh->NglobalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(mesh->NglobalGatherElements,
                                                    mesh->o_globalGatherElementList,
                                                    mesh->o_cubDiffInterpT,
                                                    mesh->o_cubInterpT,
                                                    nrs->fieldOffset,
                                                    nrs->cubatureOffset,
                                                    0,
                                                    mesh->o_invLMM,
                                                    mesh->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    nrs->o_relUrst,
                                                    o_u1,
                                                    o_rhs);
          else
            nrs->subCycleStrongVolumeKernel(mesh->NglobalGatherElements,
                mesh->o_globalGatherElementList,
                mesh->o_D,
                nrs->fieldOffset,
                0,
                mesh->o_invLMM,
                mesh->o_divU,
                extC[0],
                extC[1],
                extC[2],
                nrs->o_relUrst,
                o_u1,
                o_rhs);
        }

        oogs::start(o_rhs,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsAdd,
            nrs->gsh);

        if (mesh->NlocalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(mesh->NlocalGatherElements,
                                                    mesh->o_localGatherElementList,
                                                    mesh->o_cubDiffInterpT,
                                                    mesh->o_cubInterpT,
                                                    nrs->fieldOffset,
                                                    nrs->cubatureOffset,
                                                    0,
                                                    mesh->o_invLMM,
                                                    mesh->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    nrs->o_relUrst,
                                                    o_u1,
                                                    o_rhs);
          else
            nrs->subCycleStrongVolumeKernel(mesh->NlocalGatherElements,
                mesh->o_localGatherElementList,
                mesh->o_D,
                nrs->fieldOffset,
                0,
                mesh->o_invLMM,
                mesh->o_divU,
                extC[0],
                extC[1],
                extC[2],
                nrs->o_relUrst,
                o_u1,
                o_rhs);
        }

        oogs::finish(o_rhs,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsAdd,
            nrs->gsh);

        flops(nrs->meshV, nrs->NVfields);

        linAlg->axmyMany(mesh->Nlocal,
            nrs->NVfields,
            nrs->fieldOffset,
            0,
            1.0,
            o_LMMe,
            o_rhs);

        if (rk != 3)
          linAlg->axpbyzMany(mesh->Nlocal,
              nrs->NVfields,
              nrs->fieldOffset,
              1.0,
              o_p0,
              -sdt * nrs->coeffsfRK[rk + 1],
              o_rhs,
              o_u1);
        else
          nrs->subCycleRKKernel(mesh->Nlocal,
              nrs->NVfields,
              nrs->fieldOffset,
              sdt,
              nrs->o_weightsRK,
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
occa::memory velocitySubCycle(
    nrs_t *nrs, int nEXT, dfloat time, occa::memory o_U) {
  mesh_t *mesh = nrs->meshV;
  linAlg_t *linAlg = platform->linAlg;

  // Solve for Each SubProblem
  for (int torder = nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    nrs->subCycleInitU0Kernel(mesh->Nlocal,
        nrs->NVfields,
        nrs->fieldOffset,
        torder,
        nEXT,
        toffset,
        0,
        nrs->coeffBDF[torder],
        mesh->o_LMM,
        o_U,
        platform->o_mempool.slice0);

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0; i--)
      tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder] / nrs->Nsubsteps;

    for (int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_mempool.slice0.copyFrom(platform->o_mempool.slice0,
          nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),
          nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),
          0);

      for (int rk = 0; rk < nrs->nRK; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t = tstage + sdt * nrs->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
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

        if (mesh->NglobalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(mesh->NglobalGatherElements,
                                                    mesh->o_globalGatherElementList,
                                                    mesh->o_cubDiffInterpT,
                                                    mesh->o_cubInterpT,
                                                    nrs->fieldOffset,
                                                    nrs->cubatureOffset,
                                                    rk * nrs->NVfields * nrs->fieldOffset,
                                                    mesh->o_invLMM,
                                                    mesh->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    nrs->o_Urst,
                                                    platform->o_mempool.slice0,
                                                    platform->o_mempool.slice6);
          else
            nrs->subCycleStrongVolumeKernel(mesh->NglobalGatherElements,
                mesh->o_globalGatherElementList,
                mesh->o_D,
                nrs->fieldOffset,
                rk * nrs->NVfields * nrs->fieldOffset,
                mesh->o_invLMM,
                mesh->o_divU,
                extC[0],
                extC[1],
                extC[2],
                nrs->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice6);
        }

        occa::memory o_rhs;
        if (rk == 0)
          o_rhs = platform->o_mempool.slice6;
        if (rk == 1)
          o_rhs = platform->o_mempool.slice9;
        if (rk == 2)
          o_rhs = platform->o_mempool.slice12;
        if (rk == 3)
          o_rhs = platform->o_mempool.slice15;

        oogs::start(o_rhs,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsAdd,
            nrs->gsh);

        if (mesh->NlocalGatherElements) {
          if (platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(mesh->NlocalGatherElements,
                                                    mesh->o_localGatherElementList,
                                                    mesh->o_cubDiffInterpT,
                                                    mesh->o_cubInterpT,
                                                    nrs->fieldOffset,
                                                    nrs->cubatureOffset,
                                                    rk * nrs->NVfields * nrs->fieldOffset,
                                                    mesh->o_invLMM,
                                                    mesh->o_divU,
                                                    extC[0],
                                                    extC[1],
                                                    extC[2],
                                                    nrs->o_Urst,
                                                    platform->o_mempool.slice0,
                                                    platform->o_mempool.slice6);
          else
            nrs->subCycleStrongVolumeKernel(mesh->NlocalGatherElements,
                mesh->o_localGatherElementList,
                mesh->o_D,
                nrs->fieldOffset,
                rk * nrs->NVfields * nrs->fieldOffset,
                mesh->o_invLMM,
                mesh->o_divU,
                extC[0],
                extC[1],
                extC[2],
                nrs->o_Urst,
                platform->o_mempool.slice0,
                platform->o_mempool.slice6);
        }

        oogs::finish(o_rhs,
            nrs->NVfields,
            nrs->fieldOffset,
            ogsDfloat,
            ogsAdd,
            nrs->gsh);

        flops(nrs->meshV, nrs->NVfields);

        nrs->subCycleRKUpdateKernel(mesh->Nlocal,
            rk,
            sdt,
            nrs->fieldOffset,
            nrs->o_coeffsfRK,
            nrs->o_weightsRK,
            platform->o_mempool.slice3,
            platform->o_mempool.slice6,
            platform->o_mempool.slice0);
      }
    }
  }
  linAlg->axmyMany(mesh->Nlocal,
      3,
      nrs->fieldOffset,
      0,
      1.0,
      mesh->o_LMM,
      platform->o_mempool.slice0);
  return platform->o_mempool.slice0;
}
