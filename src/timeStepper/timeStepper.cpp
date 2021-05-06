#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp"
#include "cfl.hpp"
#include "linAlg.hpp"

#include "timeStepper.hpp"

namespace timeStepper {

double tElapsed = 0;

void step(nrs_t* nrs, dfloat time, dfloat dt, int tstep)
{
  const double tStart = MPI_Wtime();
      
  mesh_t* mesh = nrs->meshV;
  
  cds_t* cds = nrs->cds;

  coeffs(nrs, dt, tstep);

  const bool movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  if(nrs->flow) 
    nrs->extrapolateKernel(mesh->Nelements,
                           nrs->NVfields,
                           nrs->nEXT,
                           nrs->fieldOffset,
                           nrs->o_coeffEXT,
                           nrs->o_U,
                           nrs->o_Ue);

  if(nrs->Nscalar) 
    nrs->extrapolateKernel(cds->mesh[0]->Nelements,
                           cds->NSfields,
                           cds->nEXT,
                           cds->fieldOffset[0],
                           cds->o_coeffEXT,
                           cds->o_S,
                           cds->o_Se);

  dlong cubatureOffset;
  if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
    cubatureOffset = std::max(nrs->fieldOffset, mesh->Nelements * mesh->cubNp);
  else
    cubatureOffset = nrs->fieldOffset;

  if(nrs->Nsubsteps) {
    mesh_t* mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    const dlong NbyteCubature = nrs->NVfields * cubatureOffset * sizeof(dfloat);
    for (int s = nrs->nEXT; s > 1; s--) {
      const dlong Nbyte = nrs->fieldOffset * sizeof(dfloat);
      if(movingMesh){
        mesh->o_divU.copyFrom(mesh->o_divU, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
        nrs->o_relUrst.copyFrom(nrs->o_relUrst , NbyteCubature, (s - 1)*NbyteCubature, (s - 2)*NbyteCubature);
      } else {
        nrs->o_Urst.copyFrom(nrs->o_Urst , NbyteCubature, (s - 1)*NbyteCubature, (s - 2)*NbyteCubature);
      }
    }
    if(movingMesh)
      nrs->divergenceVolumeKernel(
        mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_D,
        nrs->fieldOffset,
        mesh->o_U,
        mesh->o_divU);
  }

  const bool relative = movingMesh && nrs->Nsubsteps;
  occa::memory& o_Urst = relative ? nrs->o_relUrst : nrs->o_Urst;
  mesh = nrs->meshV;
  if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
    nrs->UrstCubatureKernel(
      mesh->Nelements,
      mesh->o_cubvgeo,
      mesh->o_cubDiffInterpT,
      mesh->o_cubInterpT,
      mesh->o_cubProjectT,
      nrs->fieldOffset,
      cubatureOffset,
      nrs->o_U,
      mesh->o_U,
      o_Urst
    );
  else
    nrs->UrstKernel(
      mesh->Nelements,
      mesh->o_vgeo,
      mesh->o_D,
      nrs->fieldOffset,
      nrs->o_U,
      mesh->o_U,
      o_Urst
    );

  if(nrs->Nscalar) {
    platform->timer.tic("makeq", 1);
    platform->linAlg->fillKernel(cds->fieldOffsetSum, 0.0, cds->o_FS);
    makeq(nrs, time, tstep, cds->o_FS, cds->o_BF);
    platform->timer.toc("makeq");
  }

  if(nrs->flow) {
    platform->timer.tic("makef", 1);
    platform->linAlg->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
    makef(nrs, time, tstep, nrs->o_FU, nrs->o_BF);
    platform->timer.toc("makef");
  }

  if(movingMesh) {
    mesh_t *mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
      const dlong NbyteScalar = nrs->fieldOffset * sizeof(dfloat);
      mesh->o_LMM.copyFrom(mesh->o_LMM , NbyteScalar, (s - 1)*NbyteScalar, (s - 2)*NbyteScalar);
      mesh->o_invLMM.copyFrom(mesh->o_invLMM , NbyteScalar, (s - 1)*NbyteScalar, (s - 2)*NbyteScalar);

    }
    mesh->move();
    if(nrs->cht) nrs->meshV->computeInvLMM();
    for (int s = std::max(nrs->nEXT, mesh->nAB); s > 1; s--) {
      const dlong Nbyte = nrs->fieldOffset * nrs->NVfields * sizeof(dfloat);
      mesh->o_U.copyFrom(mesh->o_U , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    }
  } 

  platform->device.finish();
  MPI_Barrier(platform->comm.mpiComm);
  const double tPreStep = MPI_Wtime() - tStart;

  const int isOutputStep = nrs->isOutputStep;
  bool converged;
  int stage = 0;
  do {
    stage++;
    platform->device.finish();
    MPI_Barrier(platform->comm.mpiComm);
    const double tStartStep = MPI_Wtime();
     
    const dfloat timeNew = time + nrs->dt[0]; 

    if(nrs->Nscalar)
      scalarSolve(nrs, timeNew, cds->o_S, stage); 

    if(udf.properties) {
      platform->timer.tic("udfProperties", 1);
      occa::memory o_S = platform->o_mempool.slice0;
      occa::memory o_SProp = platform->o_mempool.slice0;
      if(nrs->Nscalar) {
        o_S = cds->o_S;
        o_SProp = cds->o_prop;
      }
      udf.properties(nrs, timeNew, nrs->o_U, o_S, nrs->o_prop, o_SProp);
      platform->timer.toc("udfProperties");
    }

    if(udf.div){
      linAlg_t* linAlg = platform->linAlg;
      linAlg->fill(mesh->Nlocal, 0.0, nrs->o_div);
      udf.div(nrs, timeNew, nrs->o_div);
    }


    if(nrs->flow)
      fluidSolve(nrs, timeNew, nrs->o_U, stage);

    if(platform->options.compareArgs("MESH SOLVER", "ELASTICITY"))
      meshSolve(nrs, timeNew, nrs->meshV->o_U, stage);

    platform->device.finish();
    MPI_Barrier(platform->comm.mpiComm);
    double tElapsedStep = MPI_Wtime() - tStartStep;
    if(stage == 1) tElapsedStep += tPreStep;
    tElapsed += tElapsedStep;

    converged = true;
    if(udf.converged) converged = udf.converged(nrs, stage);

    printInfo(nrs, timeNew, tstep, tElapsedStep, tElapsed);

    platform->timer.tic("udfExecuteStep", 1);
    if(isOutputStep && converged) {
      nek::ifoutfld(1);
      nrs->isOutputStep = 1;
    } 
    if(udf.executeStep) udf.executeStep(nrs, timeNew, tstep);
    nek::ifoutfld(0);
    nrs->isOutputStep = 0;
    platform->timer.toc("udfExecuteStep");
  }
  while(!converged);

  nrs->dt[2] = nrs->dt[1];
  nrs->dt[1] = nrs->dt[0];
  platform->timer.set("solve", tElapsed);
}

void coeffs(nrs_t* nrs, double dt, int tstep)
{
  nrs->dt[0] = dt;
  nrs->idt = 1/nrs->dt[0];

  const int bdfOrder = mymin(tstep, nrs->nBDF);
  const int extOrder = mymin(tstep, nrs->nEXT);

  nek::bdfCoeff(&nrs->g0, nrs->coeffBDF, nrs->dt, bdfOrder);
  nek::extCoeff(nrs->coeffEXT, nrs->dt, extOrder, bdfOrder);

  for(int i = nrs->nEXT; i > extOrder; i--) nrs->coeffEXT[i-1] = 0.0;
  for(int i = nrs->nBDF; i > bdfOrder; i--) nrs->coeffBDF[i-1] = 0.0;

  if(platform->options.compareArgs("MOVING MESH", "TRUE")) {
    mesh_t* mesh = nrs->meshV;
    if(nrs->cht) mesh = nrs->cds->mesh[0];
    const int meshOrder = mymin(tstep, mesh->nAB);
    nek::coeffAB(mesh->coeffAB, nrs->dt, meshOrder);
    for(int i = 0 ; i < meshOrder; ++i) mesh->coeffAB[i] *= nrs->dt[0];
    for(int i = mesh->nAB; i > meshOrder; i--) mesh->coeffAB[i-1] = 0.0;
    mesh->o_coeffAB.copyFrom(mesh->coeffAB, mesh->nAB * sizeof(dfloat));
  }

  nrs->ig0 = 1.0 / nrs->g0;
  nrs->o_coeffBDF.copyFrom(nrs->coeffBDF);
  nrs->o_coeffEXT.copyFrom(nrs->coeffEXT);

  if (nrs->Nscalar) {
    nrs->cds->idt = 1/nrs->cds->dt[0];
    nrs->cds->g0 = nrs->g0;
    nrs->cds->ig0 = nrs->ig0;
  }
}
void makeq(nrs_t* nrs, dfloat time, int tstep, occa::memory o_FS, occa::memory o_BF)
{
  cds_t* cds   = nrs->cds;
  mesh_t* mesh = cds->mesh[0];
  

  if(udf.sEqnSource) {
    platform->timer.tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_FS);
    platform->timer.toc("udfSEqnSource");
  }
  
  const dlong cubatureOffset = std::max(cds->vFieldOffset, cds->meshV->Nelements * cds->meshV->cubNp);

  for(int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];
    const dlong isOffset = cds->fieldOffsetScan[is];

    if(cds->options[is].compareArgs("FILTER STABILIZATION", "RELAXATION"))
      cds->filterRTKernel(
        cds->meshV->Nelements,
        nrs->o_filterMT,
        nrs->filterS,
        isOffset,
        cds->o_rho,
        cds->o_S,
        o_FS);
    const int movingMesh = cds->options[is].compareArgs("MOVING MESH", "TRUE");
    if(movingMesh && !cds->Nsubsteps){
      cds->advectMeshVelocityKernel(
        cds->meshV->Nelements,
        mesh->o_vgeo,
        mesh->o_D,
        isOffset,
        cds->vFieldOffset,
        cds->o_rho,
        mesh->o_U,
        cds->o_S,
        o_FS
      );
    }

    occa::memory o_Usubcycling = platform->o_mempool.slice0;
    if(cds->options[is].compareArgs("ADVECTION", "TRUE")) {
      if(cds->Nsubsteps) {
        if(movingMesh)
          o_Usubcycling = scalarStrongSubCycleMovingMesh(cds, mymin(tstep, cds->nEXT), time, is, cds->o_U, cds->o_S);
        else
          o_Usubcycling = scalarStrongSubCycle(cds, mymin(tstep, cds->nEXT), time, is, cds->o_U, cds->o_S);
      } else {
        if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
          cds->advectionStrongCubatureVolumeKernel(
            cds->meshV->Nelements,
            mesh->o_vgeo,
            mesh->o_cubDiffInterpT,
            mesh->o_cubInterpT,
            mesh->o_cubProjectT,
            cds->vFieldOffset,
            isOffset,
            cubatureOffset,
            cds->o_S,
            cds->o_Urst,
            cds->o_rho,
            platform->o_mempool.slice0);
        else
          cds->advectionStrongVolumeKernel(
            cds->meshV->Nelements,
            mesh->o_D,
            cds->vFieldOffset,
            isOffset,
            cds->o_S,
            cds->o_Urst,
            cds->o_rho,
            platform->o_mempool.slice0);
        platform->linAlg->axpby(
          cds->meshV->Nelements * cds->meshV->Np,
          -1.0,
          platform->o_mempool.slice0,
          1.0,
          o_FS,
          0, isOffset
        );
      }
    } else {
      platform->linAlg->fill(cds->fieldOffsetSum, 0.0, o_Usubcycling);
    } 

    cds->sumMakefKernel(
      mesh->Nlocal,
      mesh->o_LMM,
      cds->idt,
      cds->o_coeffEXT,
      cds->o_coeffBDF,
      cds->fieldOffsetSum,
      cds->fieldOffset[is],
      isOffset,
      cds->o_S,
      o_Usubcycling,
      o_FS,
      cds->o_rho,
      o_BF);
  }

  for (int s = std::max(cds->nBDF, cds->nEXT); s > 1; s--) {
    const dlong Nbyte = cds->fieldOffsetSum * sizeof(dfloat);
    cds->o_S.copyFrom (cds->o_S , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    o_FS.copyFrom(o_FS, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }
}

void scalarSolve(nrs_t* nrs, dfloat time, occa::memory o_S, int stage)
{
  cds_t* cds   = nrs->cds;
  
  platform->timer.tic("scalarSolve", 1);
  for (int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh[0];

    cds->setEllipticCoeffKernel(
      mesh->Nlocal,
      cds->g0 * cds->idt,
      cds->fieldOffsetScan[is],
      cds->fieldOffset[is],
      cds->o_diff,
      cds->o_rho,
      cds->o_ellipticCoeff);

    if(cds->o_BFDiag.ptr())
      platform->linAlg->axpby(
        mesh->Nlocal,
        1.0,
        cds->o_BFDiag,
        1.0,
        cds->o_ellipticCoeff,
        cds->fieldOffsetScan[is],
        cds->fieldOffset[is]
      );

    occa::memory o_Snew = cdsSolve(is, cds, time, stage);
    o_Snew.copyTo(o_S, cds->fieldOffset[is] * sizeof(dfloat), cds->fieldOffsetScan[is] * sizeof(dfloat));
  }
  platform->timer.toc("scalarSolve");
}

void makef(nrs_t* nrs, dfloat time, int tstep, occa::memory o_FU, occa::memory o_BF)
{
  mesh_t* mesh = nrs->meshV;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE"); 
  const int movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  if(udf.uEqnSource) {
    platform->timer.tic("udfUEqnSource", 1);
    udf.uEqnSource(nrs, time, nrs->o_U, o_FU);
    platform->timer.toc("udfUEqnSource");
  }

  if(platform->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    nrs->filterRTKernel(
      mesh->Nelements,
      nrs->o_filterMT,
      nrs->filterS,
      nrs->fieldOffset,
      nrs->o_U,
      o_FU);

  if(movingMesh && !nrs->Nsubsteps){
    nrs->advectMeshVelocityKernel(
      mesh->Nelements,
      mesh->o_vgeo,
      mesh->o_D,
      nrs->fieldOffset,
      mesh->o_U,
      nrs->o_U,
      o_FU
    );
  }

  occa::memory o_Usubcycling = platform->o_mempool.slice0;
  if(platform->options.compareArgs("ADVECTION", "TRUE")) {
    if(nrs->Nsubsteps) {
      if(movingMesh)     
        o_Usubcycling = velocityStrongSubCycleMovingMesh(nrs, mymin(tstep, nrs->nEXT), time, nrs->o_U);
      else 
        o_Usubcycling = velocityStrongSubCycle(nrs, mymin(tstep, nrs->nEXT), time, nrs->o_U);
    } else {
      if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
        nrs->advectionStrongCubatureVolumeKernel(
          mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_cubDiffInterpT,
          mesh->o_cubInterpT,
          mesh->o_cubProjectT,
          nrs->fieldOffset,
          std::max(nrs->fieldOffset, mesh->Nelements * mesh->cubNp),
          nrs->o_U,
          nrs->o_Urst,
          platform->o_mempool.slice0);
      else
        nrs->advectionStrongVolumeKernel(
          mesh->Nelements,
          mesh->o_D,
          nrs->fieldOffset,
          nrs->o_U,
          nrs->o_Urst,
          platform->o_mempool.slice0);

      platform->linAlg->axpby(
        nrs->NVfields * nrs->fieldOffset,
        -1.0,
        platform->o_mempool.slice0,
        1.0,
        o_FU
      );
    }
  } else {
    if(nrs->Nsubsteps) platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_Usubcycling);
  }

  nrs->sumMakefKernel(
    mesh->Nlocal,
    mesh->o_LMM,
    nrs->idt,
    nrs->o_coeffEXT,
    nrs->o_coeffBDF,
    nrs->fieldOffset,
    nrs->o_U,
    o_Usubcycling,
    o_FU,
    o_BF);

  if(verbose) {
    const dfloat debugNorm =
      platform->linAlg->weightedNorm2Many(
        mesh->Nlocal,
        nrs->NVfields,
        nrs->fieldOffset,
        mesh->ogs->o_invDegree,
        o_BF,
        platform->comm.mpiComm
      );
    if(platform->comm.mpiRank == 0) printf("BF norm: %.15e\n", debugNorm);
  }

  for (int s = std::max(nrs->nBDF, nrs->nEXT); s > 1; s--) {
    const dlong Nbyte = nrs->fieldOffset * nrs->NVfields * sizeof(dfloat);
    nrs->o_U.copyFrom (nrs->o_U , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    o_FU.copyFrom(o_FU, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }
}

void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U, int stage)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  platform->timer.tic("pressureSolve", 1);
  nrs->setEllipticCoeffPressureKernel(
    mesh->Nlocal,
    nrs->fieldOffset,
    nrs->o_rho,
    nrs->o_ellipticCoeff);
  occa::memory o_Pnew = tombo::pressureSolve(nrs, time, stage);
  nrs->o_P.copyFrom(o_Pnew, nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("pressureSolve");

  platform->timer.tic("velocitySolve", 1);
  nrs->setEllipticCoeffKernel(
    mesh->Nlocal,
    nrs->g0 * nrs->idt,
    0 * nrs->fieldOffset,
    nrs->fieldOffset,
    nrs->o_mue,
    nrs->o_rho,
    nrs->o_ellipticCoeff);

  occa::memory o_Unew = tombo::velocitySolve(nrs, time, stage);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("velocitySolve");
}

void meshSolve(nrs_t* nrs, dfloat time, occa::memory o_U, int stage)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  platform->timer.tic("meshSolve", 1);
  nrs->setEllipticCoeffKernel(
    mesh->Nlocal,
    1.0,
    0 * nrs->fieldOffset,
    nrs->fieldOffset,
    nrs->o_meshMue,
    nrs->o_meshRho,
    nrs->o_ellipticCoeff);
  occa::memory o_Unew = tombo::meshSolve(nrs, time, stage);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("meshSolve");
}

occa::memory velocityStrongSubCycleMovingMesh(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  occa::memory& o_p0 = platform->o_mempool.slice0;
  occa::memory& o_u1 = platform->o_mempool.slice3;

  occa::memory& o_r1 = platform->o_mempool.slice6;
  occa::memory& o_r2 = platform->o_mempool.slice9;
  occa::memory& o_r3 = platform->o_mempool.slice12;
  occa::memory& o_r4 = platform->o_mempool.slice15;

  occa::memory& o_LMMe = platform->o_mempool.slice18;

  const dlong cubatureOffset = std::max(nrs->fieldOffset, mesh->cubNp * mesh->Nelements);

  // Solve for Each SubProblem
  for (int torder = nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    const dlong offset = torder * nrs->fieldOffset;
    nrs->subCycleInitU0Kernel(
      mesh->Nlocal,
      nrs->NVfields,
      nrs->fieldOffset,
      torder,
      nEXT,
      toffset,
      offset,
      nrs->coeffBDF[torder],
      mesh->o_LMM,
      o_U,
      o_p0
    );

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder]/nrs->Nsubsteps;

    for(int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      o_u1.copyFrom(o_p0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

      for(int rk = 0; rk < nrs->nRK; ++rk) {
        occa::memory o_rhs;
        if(rk == 0) o_rhs = o_r1;
        if(rk == 1) o_rhs = o_r2;
        if(rk == 2) o_rhs = o_r3;
        if(rk == 3) o_rhs = o_r4;
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * nrs->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
        dfloat extC[3] = {0., 0., 0.};
        switch(nEXT) {
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

        nrs->nStagesSum3Kernel(
          mesh->Nlocal,
          nrs->fieldOffset,
          nEXT,
          extC[0],
          extC[1],
          extC[2],
          mesh->o_LMM,
          o_LMMe
        );
        linAlg->aydxMany(
          mesh->Nlocal,
          nrs->NVfields,
          nrs->fieldOffset,
          0,
          1.0,
          o_LMMe,
          o_u1
        );

        if(mesh->NglobalGatherElements) {
          if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              cubatureOffset,
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
            nrs->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
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

        oogs::start(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        if(mesh->NlocalGatherElements) {
          if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              cubatureOffset,
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
            nrs->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
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

        oogs::finish(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     
        linAlg->axmyMany(
          mesh->Nlocal,
          nrs->NVfields,
          nrs->fieldOffset,
          0,
          1.0,
          o_LMMe,
          o_rhs
        );

        if(rk != 3 ) 
          linAlg->axpbyzMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 1.0, o_p0, -sdt * nrs->coeffsfRK[rk+1], o_rhs, o_u1);
        else
          nrs->subCycleRKKernel(
            mesh->Nlocal,
            nrs->NVfields,
            nrs->fieldOffset,
            sdt,
            nrs->o_weightsRK,
            o_r1,
            o_r2,
            o_r3,
            o_r4,
            o_p0
          );
      }
    }
  }
  return o_p0;
}
occa::memory velocityStrongSubCycle(nrs_t* nrs, int nEXT, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = platform->linAlg;

  dlong cubatureOffset;
  if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
    cubatureOffset = std::max(nrs->fieldOffset, mesh->Nelements * mesh->cubNp);
  else
    cubatureOffset = nrs->fieldOffset;

  // Solve for Each SubProblem
  for (int torder = nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    nrs->subCycleInitU0Kernel(
      mesh->Nlocal,
      nrs->NVfields,
      nrs->fieldOffset,
      torder,
      nEXT,
      toffset,
      0,
      nrs->coeffBDF[torder],
      mesh->o_LMM,
      o_U,
      platform->o_mempool.slice0
    );

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder]/nrs->Nsubsteps;

    for(int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_mempool.slice0.copyFrom(platform->o_mempool.slice0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),
                           nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),0);

      for(int rk = 0; rk < nrs->nRK; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * nrs->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
        dfloat extC[3] = {0., 0., 0.};
        switch(nEXT) {
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

        if(mesh->NglobalGatherElements) {
          if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              cubatureOffset,
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
            nrs->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
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
        if(rk == 0) o_rhs = platform->o_mempool.slice6;
        if(rk == 1) o_rhs = platform->o_mempool.slice9;
        if(rk == 2) o_rhs = platform->o_mempool.slice12;
        if(rk == 3) o_rhs = platform->o_mempool.slice15;

        oogs::start(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        if(mesh->NlocalGatherElements) {
          if(platform->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
            nrs->subCycleStrongCubatureVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              nrs->fieldOffset,
              cubatureOffset,
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
            nrs->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
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

        oogs::finish(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        nrs->subCycleRKUpdateKernel(
          mesh->Nelements,
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
  linAlg->axmyMany(mesh->Nlocal, 3, nrs->fieldOffset, 0, 1.0, mesh->o_LMM, platform->o_mempool.slice0);
  return platform->o_mempool.slice0;
}

occa::memory scalarStrongSubCycleMovingMesh(cds_t* cds, int nEXT, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S)
{

  linAlg_t* linAlg = platform->linAlg;

  occa::memory& o_r1 = platform->o_mempool.slice2;
  occa::memory& o_r2 = platform->o_mempool.slice3;
  occa::memory& o_r3 = platform->o_mempool.slice4;
  occa::memory& o_r4 = platform->o_mempool.slice5;

  occa::memory& o_p0 = platform->o_mempool.slice0;
  occa::memory& o_u1 = platform->o_mempool.slice6;

  occa::memory& o_LMMe = platform->o_mempool.slice1;
  
  dlong cubatureOffset = std::max(cds->vFieldOffset, cds->meshV->Nelements * cds->meshV->cubNp);

  // Solve for Each SubProblem
  for (int torder = (nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = cds->fieldOffsetScan[is] +
                          torder * cds->fieldOffsetSum;
    const dlong offset = torder * cds->fieldOffset[is];
    cds->subCycleInitU0Kernel(
      cds->mesh[0]->Nlocal,
      1,
      cds->fieldOffset[is],
      torder,
      nEXT,
      toffset,
      offset,
      cds->coeffBDF[torder],
      cds->mesh[0]->o_LMM,
      o_S,
      o_p0
    );

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder]/cds->Nsubsteps;

    for(int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;
      o_u1.copyFrom(o_p0, cds->mesh[0]->Nlocal * sizeof(dfloat));
      for(int rk = 0; rk < cds->nRK; ++rk)
      {
        occa::memory o_rhs;
        if(rk == 0) o_rhs = o_r1;
        if(rk == 1) o_rhs = o_r2;
        if(rk == 2) o_rhs = o_r3;
        if(rk == 3) o_rhs = o_r4;

        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * cds->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        dfloat extC[3] = {0., 0., 0.};
        switch(nEXT) {
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
        cds->nStagesSum3Kernel(
          cds->mesh[0]->Nlocal,
          cds->vFieldOffset,
          nEXT,
          extC[0],
          extC[1],
          extC[2],
          cds->mesh[0]->o_LMM,
          o_LMMe
        );
        linAlg->aydx(cds->mesh[0]->Nlocal, 1.0, o_LMMe, o_u1);

        if(cds->meshV->NglobalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              cds->meshV->NglobalGatherElements,
              cds->meshV->o_globalGatherElementList,
              cds->meshV->o_cubDiffInterpT,
              cds->meshV->o_cubInterpT,
              cds->meshV->o_cubProjectT,
              cds->vFieldOffset,
              cubatureOffset,
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
            cds->subCycleStrongVolumeKernel(
              cds->meshV->NglobalGatherElements,
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

        oogs::start(o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        if(cds->meshV->NlocalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              cds->meshV->NlocalGatherElements,
              cds->meshV->o_localGatherElementList,
              cds->meshV->o_cubDiffInterpT,
              cds->meshV->o_cubInterpT,
              cds->meshV->o_cubProjectT,
              cds->vFieldOffset,
              cubatureOffset,
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
            cds->subCycleStrongVolumeKernel(
              cds->meshV->NlocalGatherElements,
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

        oogs::finish(o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        linAlg->axmy(cds->mesh[0]->Nlocal, 1.0, o_LMMe, o_rhs);
        if(rk != 3)
          linAlg->axpbyz(cds->mesh[0]->Nlocal, 1.0, o_p0, -sdt * cds->coeffsfRK[rk+1], o_rhs, o_u1);
        else
          cds->subCycleRKKernel(
            cds->mesh[0]->Nlocal,
            sdt,
            cds->o_weightsRK,
            o_r1,
            o_r2,
            o_r3,
            o_r4,
            o_p0
          );
      }
    }
  }
  return o_p0;
}
occa::memory scalarStrongSubCycle(cds_t* cds, int nEXT, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S)
{
  linAlg_t* linAlg= platform->linAlg;
  dlong offset = std::max(cds->vFieldOffset, cds->meshV->Nelements * cds->meshV->cubNp);

  // Solve for Each SubProblem
  for (int torder = (nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = cds->fieldOffsetScan[is] +
                          torder * cds->fieldOffsetSum;
    cds->subCycleInitU0Kernel(
      cds->mesh[0]->Nlocal,
      1,
      cds->fieldOffset[is],
      torder,
      nEXT,
      toffset,
      0,
      cds->coeffBDF[torder],
      cds->mesh[0]->o_LMM,
      o_S,
      platform->o_mempool.slice0
    );

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder]/cds->Nsubsteps;

    for(int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_mempool.slice0.copyFrom(platform->o_mempool.slice0, cds->fieldOffset[is] * sizeof(dfloat),
                           cds->fieldOffset[is] * sizeof(dfloat), 0);

      for(int rk = 0; rk < cds->nRK; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * cds->nodesRK[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        dfloat extC[3] = {0.,0.,0.};
        switch(nEXT) {
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

        if(cds->meshV->NglobalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              cds->meshV->NglobalGatherElements,
              cds->meshV->o_globalGatherElementList,
              cds->meshV->o_cubDiffInterpT,
              cds->meshV->o_cubInterpT,
              cds->meshV->o_cubProjectT,
              cds->vFieldOffset,
              offset,
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
            cds->subCycleStrongVolumeKernel(
              cds->meshV->NglobalGatherElements,
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
        if(rk == 0) o_rhs = platform->o_mempool.slice2;
        if(rk == 1) o_rhs = platform->o_mempool.slice3;
        if(rk == 2) o_rhs = platform->o_mempool.slice4;
        if(rk == 3) o_rhs = platform->o_mempool.slice5;

        oogs::start(o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        if(cds->meshV->NlocalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              cds->meshV->NlocalGatherElements,
              cds->meshV->o_localGatherElementList,
              cds->meshV->o_cubDiffInterpT,
              cds->meshV->o_cubInterpT,
              cds->meshV->o_cubProjectT,
              cds->vFieldOffset,
              offset,
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
            cds->subCycleStrongVolumeKernel(
              cds->meshV->NlocalGatherElements,
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

        oogs::finish(o_rhs, 1, cds->fieldOffset[is], ogsDfloat, ogsAdd, cds->gsh);

        cds->subCycleRKUpdateKernel(
          cds->meshV->Nelements,
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
  linAlg->axmy(cds->mesh[0]->Nlocal, 1.0, cds->mesh[0]->o_LMM, platform->o_mempool.slice0);
  return platform->o_mempool.slice0;
}

void printInfo(nrs_t *nrs, dfloat time, int tstep, double tElapsedStep, double tElapsed)
{
  cds_t *cds = nrs->cds;
      
  const int enforceVerbose = tstep < 101;
  const dfloat cfl = computeCFL(nrs);
  if(platform->comm.mpiRank == 0) {
    if(platform->options.compareArgs("VERBOSE SOLVER INFO", "TRUE") || enforceVerbose) {
      if(nrs->flow) {
        elliptic_t *solver = nrs->pSolver;
        printf("  P  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n", 
               solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
 
        if(nrs->uvwSolver) {
          solver = nrs->uvwSolver;
          printf("  UVW: iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
               solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
        } else {
          solver = nrs->uSolver;
          printf("  U  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
                 solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
          solver = nrs->vSolver;
          printf("  V  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
                 solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
          solver = nrs->wSolver;
          printf("  W  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
                 solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
        }
      }

      if(nrs->meshSolver)
      {
        elliptic_t* solver = nrs->meshSolver;
        printf("  M  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
               solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
      }
 
      for(int is = 0; is < nrs->Nscalar; is++) {
        elliptic_t * solver = cds->solver[is];
        printf("  S%02d: iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n", is,
               solver->Niter, solver->res00Norm, solver->res0Norm, solver->resNorm);
      }	
    }
    printf("step= %d  t= %.8e  dt=%.1e  C= %.2f",
           tstep, time, nrs->dt[0], cfl);

    if(nrs->flow) {
      if(nrs->uvwSolver)
        printf("  UVW: %d  P: %d", nrs->uvwSolver->Niter, nrs->pSolver->Niter);
      else
        printf("  U: %d  V: %d  W: %d  P: %d", 
       	       nrs->uSolver->Niter, nrs->vSolver->Niter, nrs->wSolver->Niter, nrs->pSolver->Niter);
    }
    if(nrs->meshSolver)
      printf("  M: %d", 
       	nrs->meshSolver->Niter);
      
    for(int is = 0; is < nrs->Nscalar; is++)
      if(cds->compute[is]) printf("  S: %d", cds->solver[is]->Niter);
 
    printf("  eTimeStep= %.2es eTime= %.5es\n", tElapsedStep, tElapsed);
  }

  if(cfl > 30 || std::isnan(cfl)) {
    if(platform->comm.mpiRank == 0) cout << "Unreasonable CFL! Dying ...\n" << endl;
    ABORT(1);
  }

  if(tstep % 10 == 0) fflush(stdout);
}

} // namespace
