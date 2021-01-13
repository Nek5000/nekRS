#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp"
#include "cfl.hpp"

void extbdfCoefficents(nrs_t* nrs, int order);

void makef(nrs_t* nrs, dfloat time, occa::memory o_FU, occa::memory o_BF);
occa::memory velocityStrongSubCycle(nrs_t* nrs, dfloat time,
                                    occa::memory o_U);
void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U);

void makeq(nrs_t* nrs, dfloat time, occa::memory o_FS, occa::memory o_BF);
occa::memory scalarStrongSubCycle(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
void scalarSolve(nrs_t* nrs, dfloat time, occa::memory o_S);

double tElapsed = 0;

void runStep(nrs_t* nrs, dfloat time, dfloat dt, int tstep)
{
  mesh_t* mesh = nrs->mesh;
  cds_t* cds = nrs->cds;

  mesh->device.finish();
  MPI_Barrier(mesh->comm);
  double tStart = MPI_Wtime();

  nrs->dt[0] = dt;

  nrs->idt = 1/nrs->dt[0];
  if(nrs->Nscalar) cds->idt = 1/cds->dt[0]; 
  extbdfCoefficents(nrs, mymin(tstep, nrs->temporalOrder));

  if(nrs->flow) 
    nrs->extrapolateKernel(mesh->Nelements,
                           nrs->NVfields,
                           nrs->ExplicitOrder,
                           nrs->fieldOffset,
                           nrs->o_extbdfA,
                           nrs->o_U,
                           nrs->o_Ue);
  if(nrs->Nscalar) 
    nrs->extrapolateKernel(cds->mesh->Nelements,
                           cds->NSfields,
                           cds->ExplicitOrder,
                           cds->fieldOffset,
                           cds->o_extbdfA,
                           cds->o_S,
                           cds->o_Se);

  if(nrs->Nscalar)
    scalarSolve(nrs, time, cds->o_S);

  if(udf.properties) {
    timer::tic("udfProperties", 1);
    occa::memory o_S = nrs->o_wrk0;
    occa::memory o_SProp = nrs->o_wrk0;
    if(nrs->Nscalar) {
      o_S = cds->o_S;
      o_SProp = cds->o_prop;
    }
    udf.properties(nrs, time + nrs->dt[0], nrs->o_U, o_S, nrs->o_prop, o_SProp);
    timer::toc("udfProperties");
  }

  if(udf.div) udf.div(nrs, time + nrs->dt[0], nrs->o_div);

  if(nrs->flow) fluidSolve(nrs, time, nrs->o_U); 

  nrs->dt[2] = nrs->dt[1];
  nrs->dt[1] = nrs->dt[0];

  mesh->device.finish();
  MPI_Barrier(mesh->comm);
  const double tElapsedStep = MPI_Wtime() - tStart;
  tElapsed += tElapsedStep;
  timer::set("solve", tElapsed);

  // print some diagnostics
  const dfloat cfl = computeCFL(nrs);
  if(mesh->rank == 0) {
    printf("step= %d  t= %.8e  dt=%.1e  C= %.2f",
           tstep, time + nrs->dt[0], nrs->dt[0], cfl);

    if(nrs->flow) {
      if(nrs->uvwSolver)
        printf("  UVW: %d  P: %d", nrs->NiterU, nrs->NiterP);
      else
        printf("  U: %d  V: %d  W: %d  P: %d", nrs->NiterU, nrs->NiterV, nrs->NiterW, nrs->NiterP);
    }

    for(int is = 0; is < nrs->Nscalar; is++)
      if(cds->compute[is]) printf("  S: %d", cds->Niter[is]);

    printf("  eTime= %.2e, %.5e s\n", tElapsedStep, tElapsed);
  }

  if(cfl > 30 || std::isnan(cfl)) {
    if(mesh->rank == 0) cout << "Unreasonable CFL! Dying ...\n" << endl;
    ABORT(1);
  }

  if(tstep % 10 == 0) fflush(stdout);
}

void extbdfCoefficents(nrs_t* nrs, int order)
{
  if(order == 1) {
    nrs->g0 = 1.0;
    nrs->extbdfB[0] = 1.0;
    nrs->extbdfB[1] = 0.0;
    nrs->extbdfB[2] = 0.0;
    nrs->extbdfA[0] = 1.0;
    nrs->extbdfA[1] = 0.0;
    nrs->extbdfA[2] = 0.0;
    nrs->ExplicitOrder = 1;
  } else if(order == 2) {
    nek_bdfCoeff(&nrs->g0, nrs->extbdfB, nrs->dt, order);
    nrs->extbdfB[2] = 0.0;
    nrs->ExplicitOrder = 2;
    nek_extCoeff(nrs->extbdfA, nrs->dt, nrs->ExplicitOrder);
    nrs->extbdfA[2] = 0.0;
  } else if(order == 3) {
    nek_bdfCoeff(&nrs->g0, nrs->extbdfB, nrs->dt, order);
    nrs->ExplicitOrder = 3;
    nek_extCoeff(nrs->extbdfA, nrs->dt, nrs->ExplicitOrder);
  }

  nrs->ig0 = 1.0 / nrs->g0;
  nrs->o_extbdfB.copyFrom(nrs->extbdfB);
  nrs->o_extbdfA.copyFrom(nrs->extbdfA);

#if 0
  if (nrs->mesh->rank == 0) {
    cout << "DT:" << nrs->dt[0] << "," << nrs->dt[1] << "," << nrs->dt[2] << "\n";
    cout << "BDF:" << nrs->g0 << "," << nrs->extbdfB[0] << "," << nrs->extbdfB[1] << "," << nrs->extbdfB[2] << "\n";
    cout << "EXT:" << nrs->extbdfA[0] << "," << nrs->extbdfA[1] << "," << nrs->extbdfA[2] << "\n";
  }
#endif

  if (nrs->Nscalar) {
    nrs->cds->ExplicitOrder = nrs->ExplicitOrder;
    nrs->cds->g0 = nrs->g0;
    nrs->cds->ig0 = nrs->ig0;
  }
}

void makeq(nrs_t* nrs, dfloat time, occa::memory o_FS, occa::memory o_BF)
{
  cds_t* cds   = nrs->cds;
  mesh_t* mesh = cds->mesh;

  if(udf.sEqnSource) {
    timer::tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_FS);
    timer::toc("udfSEqnSource");
  }

  for(int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;
    const dlong isOffset = is * cds->fieldOffset;
    occa::memory o_adv = cds->o_wrk0;

    if(cds->options[is].compareArgs("FILTER STABILIZATION", "RELAXATION"))
      cds->filterRTKernel(
        cds->meshV->Nelements,
        nrs->o_filterMT,
        nrs->filterS,
        isOffset,
        cds->o_rho,
        cds->o_S,
        o_FS);

    if(cds->options[is].compareArgs("ADVECTION", "TRUE")) {
      if(cds->Nsubsteps) {
        o_adv = scalarStrongSubCycle(cds, time, is, cds->o_U, cds->o_S);
      } else {
        if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
          cds->advectionStrongCubatureVolumeKernel(
            cds->meshV->Nelements,
            mesh->o_vgeo,
            mesh->o_cubvgeo,
            mesh->o_cubDiffInterpT,
            mesh->o_cubInterpT,
            mesh->o_cubProjectT,
            cds->vFieldOffset,
            isOffset,
            cds->o_U,
            cds->o_S,
            cds->o_rho,
            cds->o_wrk0);
        else
          cds->advectionStrongVolumeKernel(
            cds->meshV->Nelements,
            mesh->o_vgeo,
            mesh->o_Dmatrices,
            cds->vFieldOffset,
            isOffset,
            cds->o_U,
            cds->o_S,
            cds->o_rho,
            cds->o_wrk0);
 
        nrs->scaledAddKernel(
          cds->meshV->Nelements * cds->meshV->Np,
          -1.0,
          0 * cds->fieldOffset,
          cds->o_wrk0,
          1.0,
          isOffset,
          o_FS);
      }
    } else {
      cds->fillKernel(cds->fieldOffset * cds->NVfields, 0.0, o_adv);
    } 

    cds->sumMakefKernel(
      mesh->Nelements,
      mesh->o_vgeo,
      cds->idt,
      cds->o_extbdfA,
      cds->o_extbdfB,
      cds->fieldOffset * cds->NSfields,
      isOffset,
      cds->o_S,
      o_adv,
      o_FS,
      cds->o_rho,
      o_BF);
  }
}

void scalarSolve(nrs_t* nrs, dfloat time, occa::memory o_S)
{
  cds_t* cds   = nrs->cds;

  timer::tic("makeq", 1);
  cds->fillKernel(cds->fieldOffset * cds->NSfields, 0.0, cds->o_FS);
  makeq(nrs, time, cds->o_FS, cds->o_BF);
  timer::toc("makeq");

  for (int s = cds->Nstages; s > 1; s--) {
    const dlong Nbyte = cds->fieldOffset * cds->NSfields * sizeof(dfloat);
    cds->o_FS.copyFrom(cds->o_FS, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    cds->o_S.copyFrom (cds->o_S , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }

  timer::tic("scalarSolve", 1);
  for (int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;

    cds->setEllipticCoeffKernel(
      cds->Nlocal,
      cds->g0 * cds->idt,
      is * cds->fieldOffset,
      cds->fieldOffset,
      cds->o_diff,
      cds->o_rho,
      cds->o_ellipticCoeff);

    if(cds->o_BFDiag.ptr())
      cds->scaledAddKernel(
        cds->Nlocal,
        1.0,
        is * cds->fieldOffset,
        cds->o_BFDiag,
        1.0,
        cds->fieldOffset,
        cds->o_ellipticCoeff);

    occa::memory o_Snew = cdsSolve(is, cds, time + cds->dt[0]);
    o_Snew.copyTo(o_S, cds->Ntotal * sizeof(dfloat), is * cds->fieldOffset * sizeof(dfloat));
  }
  timer::toc("scalarSolve");
}

void makef(nrs_t* nrs, dfloat time, occa::memory o_FU, occa::memory o_BF)
{
  mesh_t* mesh = nrs->mesh;

  if(udf.uEqnSource) {
    timer::tic("udfUEqnSource", 1);
    udf.uEqnSource(nrs, time, nrs->o_U, o_FU);
    timer::toc("udfUEqnSource");
  }

  if(nrs->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    nrs->filterRTKernel(
      mesh->Nelements,
      nrs->o_filterMT,
      nrs->filterS,
      nrs->fieldOffset,
      nrs->o_U,
      o_FU);

  occa::memory o_adv = nrs->o_wrk0;
  if(nrs->options.compareArgs("ADVECTION", "TRUE")) {
    if(nrs->Nsubsteps) {
      o_adv = velocityStrongSubCycle(nrs, time, nrs->o_U);
    } else {
      if(nrs->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
        nrs->advectionStrongCubatureVolumeKernel(
          mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_cubvgeo,
          mesh->o_cubDiffInterpT,
          mesh->o_cubInterpT,
          mesh->o_cubProjectT,
          nrs->fieldOffset,
          nrs->o_U,
          nrs->o_wrk0);
      else
        nrs->advectionStrongVolumeKernel(
          mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_Dmatrices,
          nrs->fieldOffset,
          nrs->o_U,
          nrs->o_wrk0);
 
      nrs->scaledAddKernel(
        nrs->NVfields * nrs->fieldOffset,
        -1.0,
        0,
        nrs->o_wrk0,
        1.0,
        0,
        o_FU);
    }
  } else {
    if(nrs->Nsubsteps) nrs->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, o_adv);
  }

  nrs->sumMakefKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    nrs->idt,
    nrs->o_extbdfA,
    nrs->o_extbdfB,
    nrs->fieldOffset,
    nrs->o_U,
    o_adv,
    o_FU,
    o_BF);
}

void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->mesh;

  timer::tic("makef", 1);
  nrs->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
  makef(nrs, time, nrs->o_FU, nrs->o_BF);
  timer::toc("makef");

  for (int s = nrs->Nstages; s > 1; s--) {
    const dlong Nbyte = nrs->fieldOffset * nrs->NVfields * sizeof(dfloat);
    nrs->o_FU.copyFrom(nrs->o_FU, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    nrs->o_U.copyFrom (nrs->o_U , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }

  timer::tic("pressureSolve", 1);
  nrs->setEllipticCoeffPressureKernel(
    nrs->Nlocal,
    nrs->fieldOffset,
    nrs->o_rho,
    nrs->o_ellipticCoeff);
  occa::memory o_Pnew = tombo::pressureSolve(nrs, time + nrs->dt[0]);
  nrs->o_P.copyFrom(o_Pnew, nrs->Ntotal * sizeof(dfloat));
  timer::toc("pressureSolve");

  timer::tic("velocitySolve", 1);
  nrs->setEllipticCoeffKernel(
    nrs->Nlocal,
    nrs->g0 * nrs->idt,
    0 * nrs->fieldOffset,
    nrs->fieldOffset,
    nrs->o_mue,
    nrs->o_rho,
    nrs->o_ellipticCoeff);

  occa::memory o_Unew = tombo::velocitySolve(nrs, time + nrs->dt[0]);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  timer::toc("velocitySolve");
}

occa::memory velocityStrongSubCycle(nrs_t* nrs, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->mesh;

  // Solve for Each SubProblem
  for (int torder = nrs->ExplicitOrder - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    const dfloat b = nrs->extbdfB[torder];
    if (torder == nrs->ExplicitOrder - 1)
      nrs->scaledAddKernel(nrs->NVfields * nrs->fieldOffset, b, toffset,
                           o_U, 0.0, 0, nrs->o_wrk0);
    else
      nrs->scaledAddKernel(nrs->NVfields * nrs->fieldOffset, b, toffset,
                           o_U, 1.0, 0, nrs->o_wrk0);

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder]/nrs->Nsubsteps;

    for(int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      nrs->o_wrk0.copyFrom(nrs->o_wrk0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),
                           nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),0);

      for(int rk = 0; rk < nrs->SNrk; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * nrs->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
        switch(nrs->ExplicitOrder) {
        case 1:
          nrs->extC[0] = 1;
          nrs->extC[1] = 0;
          nrs->extC[2] = 0;
          break;
        case 2:
          nrs->extC[0] = (t - tn1) / (tn0 - tn1);
          nrs->extC[1] = (t - tn0) / (tn1 - tn0);
          nrs->extC[2] = 0;
          break;
        case 3:
          nrs->extC[0] = (t - tn1) * (t - tn2) / ((tn0 - tn1) * (tn0 - tn2));
          nrs->extC[1] = (t - tn0) * (t - tn2) / ((tn1 - tn0) * (tn1 - tn2));
          nrs->extC[2] = (t - tn0) * (t - tn1) / ((tn2 - tn0) * (tn2 - tn1));
          break;
        }
        nrs->o_extC.copyFrom(nrs->extC);

        if(mesh->NglobalGatherElements) {
          if(nrs->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              nrs->o_extC,
              o_U,
              nrs->o_wrk0,
              nrs->o_wrk6);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_vgeo,
              mesh->o_Dmatrices,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              nrs->o_extC,
              o_U,
              nrs->o_wrk0,
              nrs->o_wrk6);
        }

        occa::memory o_rhs;
        if(rk == 0) o_rhs = nrs->o_wrk6;
        if(rk == 1) o_rhs = nrs->o_wrk9;
        if(rk == 2) o_rhs = nrs->o_wrk12;
        if(rk == 3) o_rhs = nrs->o_wrk15;

        oogs::start(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        if(mesh->NlocalGatherElements) {
          if(nrs->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              nrs->o_extC,
              o_U,
              nrs->o_wrk0,
              nrs->o_wrk6);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_vgeo,
              mesh->o_Dmatrices,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              nrs->o_extC,
              o_U,
              nrs->o_wrk0,
              nrs->o_wrk6);
        }

        oogs::finish(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        nrs->subCycleRKUpdateKernel(
          mesh->Nelements,
          rk,
          sdt,
          nrs->fieldOffset,
          nrs->o_Srka,
          nrs->o_Srkb,
          nrs->o_wrk3,
          nrs->o_wrk6,
          nrs->o_wrk0);
      }
    }
  }
  return nrs->o_wrk0;
}

occa::memory scalarStrongSubCycle(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S)
{
  mesh_t* mesh = cds->meshV;

  // Solve for Each SubProblem
  for (int torder = (cds->ExplicitOrder - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = is * cds->fieldOffset +
                          torder * cds->NSfields * cds->fieldOffset;
    if (torder == cds->ExplicitOrder - 1)
      cds->scaledAddKernel(cds->fieldOffset, cds->extbdfB[torder],
                           toffset, o_S, 0.0, 0, cds->o_wrk0);
    else
      cds->scaledAddKernel(cds->fieldOffset, cds->extbdfB[torder],
                           toffset, o_S, 1.0, 0, cds->o_wrk0);

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder]/cds->Nsubsteps;

    for(int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      cds->o_wrk0.copyFrom(cds->o_wrk0, cds->fieldOffset * sizeof(dfloat),
                           cds->fieldOffset * sizeof(dfloat), 0);

      for(int rk = 0; rk < cds->SNrk; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * cds->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        switch(cds->ExplicitOrder) {
        case 1:
          cds->extC[0] = 1;
          cds->extC[1] = 0;
          cds->extC[2] = 0;
          break;
        case 2:
          cds->extC[0] = (t - tn1) / (tn0 - tn1);
          cds->extC[1] = (t - tn0) / (tn1 - tn0);
          cds->extC[2] = 0;
          break;
        case 3:
          cds->extC[0] = (t - tn1) * (t - tn2) / ((tn0 - tn1) * (tn0 - tn2));
          cds->extC[1] = (t - tn0) * (t - tn2) / ((tn1 - tn0) * (tn1 - tn2));
          cds->extC[2] = (t - tn0) * (t - tn1) / ((tn2 - tn0) * (tn2 - tn1));
          break;
        }
        cds->o_extC.copyFrom(cds->extC);

        if(mesh->NglobalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              mesh->o_invLMM,
              cds->o_extC,
              o_U,
              cds->o_wrk0,
              cds->o_wrk2);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_Dmatrices,
              mesh->o_invLMM,
              cds->o_extC,
              o_U,
              cds->o_wrk0,
              cds->o_wrk2);
        }

        occa::memory o_rhs;
        if(rk == 0) o_rhs = cds->o_wrk2;
        if(rk == 1) o_rhs = cds->o_wrk3;
        if(rk == 2) o_rhs = cds->o_wrk4;
        if(rk == 3) o_rhs = cds->o_wrk5;

        oogs::start(o_rhs, 1, cds->fieldOffset, ogsDfloat, ogsAdd, cds->gsh);

        if(mesh->NlocalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              mesh->o_invLMM,
              cds->o_extC, 
              o_U,
              cds->o_wrk0,
              cds->o_wrk2);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_Dmatrices,
              mesh->o_invLMM,
              cds->o_extC,
              o_U,
              cds->o_wrk0,
              cds->o_wrk2);
        }

        oogs::finish(o_rhs, 1, cds->fieldOffset, ogsDfloat, ogsAdd, cds->gsh);

        cds->subCycleRKUpdateKernel(
          mesh->Nelements,
          rk,
          sdt,
          cds->fieldOffset,
          cds->o_Srka,
          cds->o_Srkb,
          cds->o_wrk1,
          cds->o_wrk2,
          cds->o_wrk0);
      }
    }
  }
  return cds->o_wrk0;
}
