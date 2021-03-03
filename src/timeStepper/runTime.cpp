#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp"
#include "cfl.hpp"
#include "linAlg.hpp"

void computeCoefficients(nrs_t* nrs, int order, int meshOrder);

void makef(nrs_t* nrs, dfloat time, occa::memory o_FU, occa::memory o_BF);
occa::memory velocityStrongSubCycle(nrs_t* nrs, dfloat time,
                                    occa::memory o_U);
occa::memory velocityStrongSubCycleMovingMesh(nrs_t* nrs, dfloat time,
                                    occa::memory o_U);
void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U);

void makeq(nrs_t* nrs, dfloat time, occa::memory o_FS, occa::memory o_BF);
occa::memory scalarStrongSubCycleMovingMesh(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
occa::memory scalarStrongSubCycle(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S);
void scalarSolve(nrs_t* nrs, dfloat time, occa::memory o_S);
void meshUpdate(nrs_t* nrs);
void meshUpdateBdivW(nrs_t* nrs);

double tElapsed = 0;

void runStep(nrs_t* nrs, dfloat time, dfloat dt, int tstep)
{
  mesh_t* mesh = nrs->meshV;
  
  cds_t* cds = nrs->cds;

  platform->device.finish();
  MPI_Barrier(platform->comm.mpiComm);
  double tStart = MPI_Wtime();

  nrs->dt[0] = dt;

  nrs->idt = 1/nrs->dt[0];
  if(nrs->Nscalar) cds->idt = 1/cds->dt[0]; 
  computeCoefficients(nrs, mymin(tstep, nrs->nBDF), mymin(tstep, mesh->nAB));

  for(int geom = 0; geom < 2; geom++) {
    if(nrs->flow && geom == 0) 
      nrs->extrapolateKernel(mesh->Nelements,
                             nrs->NVfields,
                             nrs->nEXT,
                             nrs->fieldOffset,
                             nrs->o_coeffEXT,
                             nrs->o_U,
                             nrs->o_Ue);
    if(nrs->Nscalar && geom == 0) 
      nrs->extrapolateKernel(cds->meshT[0]->Nelements,
                             cds->NSfields,
                             cds->nEXT,
                             cds->fieldOffset,
                             cds->o_coeffEXT,
                             cds->o_S,
                             cds->o_Se);

    const bool movingMesh = nrs->options.compareArgs("MOVING MESH", "TRUE");
    if(movingMesh && geom == 0) meshUpdateBdivW(nrs);
    if(movingMesh && geom == 1) meshUpdate(nrs); 

    if(nrs->Nscalar) {
      if(geom == 0) {
        platform->timer.tic("makeq", 1);
        platform->linAlg->fillKernel(cds->fieldOffset * cds->NSfields, 0.0, cds->o_FS);
        makeq(nrs, time, cds->o_FS, cds->o_BF);
        platform->timer.toc("makeq");
      }
      if(geom == 1) scalarSolve(nrs, time, cds->o_S); 
    }

    if(udf.properties && geom == 1) {
      platform->timer.tic("udfProperties", 1);
      occa::memory o_S = platform->o_slice0;
      occa::memory o_SProp = platform->o_slice0;
      if(nrs->Nscalar) {
        o_S = cds->o_S;
        o_SProp = cds->o_prop;
      }
      udf.properties(nrs, time + nrs->dt[0], nrs->o_U, o_S, nrs->o_prop, o_SProp);
      platform->timer.toc("udfProperties");
    }

    if(udf.div && geom == 1){
      linAlg_t* linAlg = nrs->linAlg;
      linAlg->fill(mesh->Nlocal, 0.0, nrs->o_div);
      udf.div(nrs, time + nrs->dt[0], nrs->o_div);
    }

    if(nrs->flow) {
      if(geom == 0) {
        platform->timer.tic("makef", 1);
        platform->linAlg->fillKernel(nrs->fieldOffset * nrs->NVfields, 0.0, nrs->o_FU);
        makef(nrs, time, nrs->o_FU, nrs->o_BF);
        platform->timer.toc("makef");
      }
      if(geom == 1) fluidSolve(nrs, time, nrs->o_U); 
      if(geom == 0) nrs->meshT->solve(); 
    }
  }

  nrs->dt[2] = nrs->dt[1];
  nrs->dt[1] = nrs->dt[0];

  platform->device.finish();
  MPI_Barrier(platform->comm.mpiComm);
  const double tElapsedStep = MPI_Wtime() - tStart;
  tElapsed += tElapsedStep;
  platform->timer.set("solve", tElapsed);

  // print some diagnostics
  const dfloat cfl = computeCFL(nrs);
  if(platform->comm.mpiRank == 0) {
    printf("step= %d  t= %.8e  dt=%.1e  C= %.2f",
           tstep, time + nrs->dt[0], nrs->dt[0], cfl);

    if(nrs->options.compareArgs("VERBOSE SOLVER INFO", "TRUE") || tstep < 101) {
      printf("  eTime= %.2e, %.5e s\n", tElapsedStep, tElapsed);
      if(nrs->flow) {
        elliptic_t *solver = nrs->pSolver;
        printf("  P  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n", 
	       solver->Niter, solver->res00, solver->res0, solver->res);
 
        if(nrs->uvwSolver) {
          solver = nrs->uvwSolver;
          printf("  UVW: iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",

	       solver->Niter, solver->res00, solver->res0, solver->res);
        } else {
          solver = nrs->uSolver;
          printf("  U  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
	         solver->Niter, solver->res00, solver->res0, solver->res);
          solver = nrs->vSolver;
          printf("  V  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
	         solver->Niter, solver->res00, solver->res0, solver->res);
          solver = nrs->wSolver;
          printf("  W  : iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n",
	         solver->Niter, solver->res00, solver->res0, solver->res);
        }
      }

      for(int is = 0; is < nrs->Nscalar; is++) {
	elliptic_t * solver = cds->solver[is];
        printf("  S%02d: iter %03d  resNorm00 %e  resNorm0 %e  resNorm %e\n", is,
	       solver->Niter, solver->res00, solver->res0, solver->res);
      }	
    }  else {
      if(nrs->flow) {
        if(nrs->uvwSolver)
          printf("  UVW: %d  P: %d", nrs->uvwSolver->Niter, nrs->pSolver->Niter);
        else
          printf("  U: %d  V: %d  W: %d  P: %d", 
         	       nrs->uSolver->Niter, nrs->vSolver->Niter, nrs->wSolver->Niter, nrs->pSolver->Niter);
      }
      for(int is = 0; is < nrs->Nscalar; is++)
        if(cds->compute[is]) printf("  S: %d", cds->solver[is]->Niter);

      printf("  eTime= %.2e, %.5e s\n", tElapsedStep, tElapsed);
    }
  }

  if(cfl > 30 || std::isnan(cfl)) {
    if(platform->comm.mpiRank == 0) cout << "Unreasonable CFL! Dying ...\n" << endl;
    ABORT(1);
  }

  if(tstep % 10 == 0) fflush(stdout);
}

void computeCoefficients(nrs_t* nrs, int order, int meshOrder)
{
  mesh_t* mesh = nrs->meshV;
  if(order == 1) {
    nrs->g0 = 1.0;
    nrs->coeffBDF[0] = 1.0;
    nrs->coeffBDF[1] = 0.0;
    nrs->coeffBDF[2] = 0.0;
    nrs->coeffEXT[0] = 1.0;
    nrs->coeffEXT[1] = 0.0;
    nrs->coeffEXT[2] = 0.0;
    nrs->nEXT = 1;
  } else if(order == 2) {
    nek::bdfCoeff(&nrs->g0, nrs->coeffBDF, nrs->dt, order);
    nrs->coeffBDF[2] = 0.0;

    nrs->nEXT = 2;
    nek::extCoeff(nrs->coeffEXT, nrs->dt, nrs->nEXT);
    nrs->coeffEXT[2] = 0.0;
  } else if(order == 3) {
    nek::bdfCoeff(&nrs->g0, nrs->coeffBDF, nrs->dt, order);
    nrs->nEXT = 3;
    nek::extCoeff(nrs->coeffEXT, nrs->dt, nrs->nEXT);
  }

  if(nrs->options.compareArgs("MOVING MESH", "TRUE"))
  {
    const int maxIntegrationOrder = 3;
    for(int i = 0 ; i < maxIntegrationOrder; ++i){
      mesh->coeffAB[i] = 0.0;
    }
    nek::coeffAB(mesh->coeffAB, nrs->dt, meshOrder);
    for(int i = 0 ; i < maxIntegrationOrder; ++i){
      mesh->coeffAB[i] *= nrs->dt[0];
    }
    mesh->o_coeffAB.copyFrom(mesh->coeffAB, maxIntegrationOrder * sizeof(dfloat));
  }

  nrs->ig0 = 1.0 / nrs->g0;
  nrs->o_coeffBDF.copyFrom(nrs->coeffBDF);
  nrs->o_coeffEXT.copyFrom(nrs->coeffEXT);

#if 0
  if (platform->comm.mpiRank == 0) {
    cout << "DT:" << nrs->dt[0] << "," << nrs->dt[1] << "," << nrs->dt[2] << "\n";
    cout << "BDF:" << nrs->g0 << "," << nrs->coeffBDF[0] << "," << nrs->coeffBDF[1] << "," << nrs->coeffBDF[2] << "\n";
    cout << "EXT:" << nrs->coeffEXT[0] << "," << nrs->coeffEXT[1] << "," << nrs->coeffEXT[2] << "\n";
  }
#endif

  if (nrs->Nscalar) {
    nrs->cds->nEXT = nrs->nEXT;
    nrs->cds->g0 = nrs->g0;
    nrs->cds->ig0 = nrs->ig0;
  }
}
void makeq(nrs_t* nrs, dfloat time, occa::memory o_FS, occa::memory o_BF)
{
  cds_t* cds   = nrs->cds;
  mesh_t* mesh = cds->meshT[0];
  

  if(udf.sEqnSource) {
    platform->timer.tic("udfSEqnSource", 1);
    udf.sEqnSource(nrs, time, cds->o_S, o_FS);
    platform->timer.toc("udfSEqnSource");
  }

  for(int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->meshT[0];
    const dlong isOffset = is * cds->fieldOffset;
    occa::memory o_adv = platform->o_slice0;

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
        mesh->Nelements,
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

    if(cds->options[is].compareArgs("ADVECTION", "TRUE")) {
      if(cds->Nsubsteps) {
        o_adv = movingMesh ? scalarStrongSubCycleMovingMesh(cds, time, is, cds->o_U, cds->o_S) :
          scalarStrongSubCycle(cds, time, is, cds->o_U, cds->o_S);
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
            platform->o_slice0);
        else
          cds->advectionStrongVolumeKernel(
            cds->meshV->Nelements,
            mesh->o_vgeo,
            mesh->o_D,
            cds->vFieldOffset,
            isOffset,
            cds->o_U,
            cds->o_S,
            cds->o_rho,
            platform->o_slice0);
        platform->linAlg->axpby(
          cds->meshV->Nelements * cds->meshV->Np,
          -1.0,
          platform->o_slice0,
          1.0,
          o_FS,
          0, isOffset
        );
      }
    } else {
      platform->linAlg->fill(cds->fieldOffset * cds->NVfields, 0.0, o_adv);
    } 

    cds->sumMakefKernel(
      mesh->Nelements,
      mesh->o_LMM,
      cds->idt,
      cds->o_coeffEXT,
      cds->o_coeffBDF,
      cds->fieldOffset * cds->NSfields,
      cds->fieldOffset,
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
  

  for (int s = cds->nBDF; s > 1; s--) {
    const dlong Nbyte = cds->fieldOffset * cds->NSfields * sizeof(dfloat);
    cds->o_FS.copyFrom(cds->o_FS, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    cds->o_S.copyFrom (cds->o_S , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }

  platform->timer.tic("scalarSolve", 1);
  for (int is = 0; is < cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t* mesh;
    (is) ? mesh = cds->meshV : mesh = cds->meshT[0];

    cds->setEllipticCoeffKernel(
      mesh->Nlocal,
      cds->g0 * cds->idt,
      is * cds->fieldOffset,
      cds->fieldOffset,
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
        is * cds->fieldOffset,
        cds->fieldOffset
      );

    occa::memory o_Snew = cdsSolve(is, cds, time + cds->dt[0]);
    o_Snew.copyTo(o_S, cds->fieldOffset * sizeof(dfloat), is * cds->fieldOffset * sizeof(dfloat));
  }
  platform->timer.toc("scalarSolve");
}

void makef(nrs_t* nrs, dfloat time, occa::memory o_FU, occa::memory o_BF)
{
  mesh_t* mesh = nrs->meshV;
  

  if(udf.uEqnSource) {
    platform->timer.tic("udfUEqnSource", 1);
    udf.uEqnSource(nrs, time, nrs->o_U, o_FU);
    platform->timer.toc("udfUEqnSource");
  }

  if(nrs->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    nrs->filterRTKernel(
      mesh->Nelements,
      nrs->o_filterMT,
      nrs->filterS,
      nrs->fieldOffset,
      nrs->o_U,
      o_FU);
  const int movingMesh = nrs->options.compareArgs("MOVING MESH", "TRUE");
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

  occa::memory o_adv = platform->o_slice0;
  if(nrs->options.compareArgs("ADVECTION", "TRUE")) {
    if(nrs->Nsubsteps) {
      o_adv = movingMesh ? velocityStrongSubCycleMovingMesh(nrs, time, nrs->o_U) :
              velocityStrongSubCycle(nrs, time, nrs->o_U);
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
          platform->o_slice0);
      else
        nrs->advectionStrongVolumeKernel(
          mesh->Nelements,
          mesh->o_vgeo,
          mesh->o_D,
          nrs->fieldOffset,
          nrs->o_U,
          platform->o_slice0);
      platform->linAlg->axpby(
        nrs->NVfields * nrs->fieldOffset,
        -1.0,
        platform->o_slice0,
        1.0,
        o_FU
      );
    }
  } else {
    if(nrs->Nsubsteps) platform->linAlg->fill(nrs->fieldOffset * nrs->NVfields, 0.0, o_adv);
  }

  nrs->sumMakefKernel(
    mesh->Nelements,
    mesh->o_LMM,
    nrs->idt,
    nrs->o_coeffEXT,
    nrs->o_coeffBDF,
    nrs->fieldOffset,
    nrs->o_U,
    o_adv,
    o_FU,
    o_BF);
}

void fluidSolve(nrs_t* nrs, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = nrs->linAlg;

  for (int s = nrs->nBDF; s > 1; s--) {
    const dlong Nbyte = nrs->fieldOffset * nrs->NVfields * sizeof(dfloat);
    nrs->o_FU.copyFrom(nrs->o_FU, Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
    nrs->o_U.copyFrom (nrs->o_U , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }

  platform->timer.tic("pressureSolve", 1);
  nrs->setEllipticCoeffPressureKernel(
    mesh->Nlocal,
    nrs->fieldOffset,
    nrs->o_rho,
    nrs->o_ellipticCoeff);
  occa::memory o_Pnew = tombo::pressureSolve(nrs, time + nrs->dt[0]);
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

  occa::memory o_Unew = tombo::velocitySolve(nrs, time + nrs->dt[0]);
  o_U.copyFrom(o_Unew, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));
  platform->timer.toc("velocitySolve");
}

occa::memory velocityStrongSubCycleMovingMesh(nrs_t* nrs, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = mesh->linAlg;

  occa::memory& o_p0 = platform->o_slice0;
  occa::memory& o_u1 = platform->o_slice3;

  occa::memory& o_r1 = platform->o_slice6;
  occa::memory& o_r2 = platform->o_slice9;
  occa::memory& o_r3 = platform->o_slice12;
  occa::memory& o_r4 = platform->o_slice15;

  occa::memory o_bmst = platform->o_slice18;

  // Solve for Each SubProblem
  for (int torder = nrs->nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    const dlong offset = torder * nrs->fieldOffset;
    const dfloat b = nrs->coeffBDF[torder];
    if (torder == nrs->nEXT - 1){
      linAlg->axpby(nrs->NVfields * nrs->fieldOffset, b, o_U, 0.0, o_p0, toffset, 0);
      occa::memory o_LMM_slice = mesh->o_LMM + offset * sizeof(dfloat);
      linAlg->axmyMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 0, 1.0, o_LMM_slice, o_p0);
    }
    else{
      nrs->subCycleExtrapolateFieldKernel(
        mesh->Nlocal,
        nrs->NVfields,
        nrs->fieldOffset,
        toffset,
        offset,
        b,
        mesh->o_LMM,
        o_U,
        o_p0
      );
    }

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder]/nrs->Nsubsteps;

    for(int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      o_u1.copyFrom(o_p0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat));

      for(int rk = 0; rk < nrs->SNrk; ++rk) {
        occa::memory o_rhs;
        if(rk == 0) o_rhs = o_r1;
        if(rk == 1) o_rhs = o_r2;
        if(rk == 2) o_rhs = o_r3;
        if(rk == 3) o_rhs = o_r4;
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * nrs->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
        switch(nrs->nEXT) {
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
        nrs->subCycleExtrapolateScalarKernel(mesh->Nlocal, nrs->nEXT, nrs->fieldOffset, nrs->o_extC, mesh->o_LMM, o_bmst);
        linAlg->aydxMany(
          mesh->Nlocal,
          nrs->NVfields,
          nrs->fieldOffset,
          0,
          1.0,
          o_bmst,
          o_u1
        );

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
              0,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              o_u1,
              o_rhs);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_vgeo,
              mesh->o_D,
              nrs->fieldOffset,
              0,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              o_u1,
              o_rhs);
        }

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
              0,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              o_u1,
              o_rhs);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_vgeo,
              mesh->o_D,
              nrs->fieldOffset,
              0,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
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
          o_bmst,
          o_rhs
        );

        if(rk != 3 ) linAlg->axpbyzMany(mesh->Nlocal, nrs->NVfields, nrs->fieldOffset, 1.0, o_p0, -sdt * nrs->Srka[rk+1], o_rhs, o_u1);
        else{
          nrs->subCycleRKKernel(
            mesh->Nlocal,
            nrs->NVfields,
            nrs->fieldOffset,
            sdt,
            nrs->o_Srkb,
            o_r1,
            o_r2,
            o_r3,
            o_r4,
            o_p0
          );
        }
      }
    }
  }
  return o_p0;
}
occa::memory velocityStrongSubCycle(nrs_t* nrs, dfloat time, occa::memory o_U)
{
  mesh_t* mesh = nrs->meshV;
  linAlg_t* linAlg = mesh->linAlg;

  // Solve for Each SubProblem
  for (int torder = nrs->nEXT - 1; torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder * nrs->NVfields * nrs->fieldOffset;
    const dfloat b = nrs->coeffBDF[torder];
    if (torder == nrs->nEXT - 1){
      platform->linAlg->axpby(
        nrs->NVfields * nrs->fieldOffset,
        b,
        o_U,
        0.0,
        platform->o_slice0,
        toffset,
        0
      );
    }
    else{
      platform->linAlg->axpby(
        nrs->NVfields * nrs->fieldOffset,
        b,
        o_U,
        1.0,
        platform->o_slice0,
        toffset,
        0
      );
    }

    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= nrs->dt[i];
    const dfloat sdt = nrs->dt[torder]/nrs->Nsubsteps;

    for(int ststep = 0; ststep < nrs->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_slice0.copyFrom(platform->o_slice0, nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),
                           nrs->NVfields * nrs->fieldOffset * sizeof(dfloat),0);

      for(int rk = 0; rk < nrs->SNrk; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * nrs->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - nrs->dt[1];
        const dfloat tn2 = time - (nrs->dt[1] + nrs->dt[2]);
        switch(nrs->nEXT) {
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
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice6);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              mesh->o_vgeo,
              mesh->o_D,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice6);
        }

        occa::memory o_rhs;
        if(rk == 0) o_rhs = platform->o_slice6;
        if(rk == 1) o_rhs = platform->o_slice9;
        if(rk == 2) o_rhs = platform->o_slice12;
        if(rk == 3) o_rhs = platform->o_slice15;

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
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice6);
          else
            nrs->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              mesh->o_vgeo,
              mesh->o_D,
              nrs->fieldOffset,
              rk * nrs->NVfields * nrs->fieldOffset,
              mesh->o_invLMM,
              mesh->o_BdivW,
              nrs->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice6);
        }

        oogs::finish(o_rhs, nrs->NVfields, nrs->fieldOffset,ogsDfloat, ogsAdd, nrs->gsh);                     

        nrs->subCycleRKUpdateKernel(
          mesh->Nelements,
          rk,
          sdt,
          nrs->fieldOffset,
          nrs->o_Srka,
          nrs->o_Srkb,
          platform->o_slice3,
          platform->o_slice6,
          platform->o_slice0);
      }
    }
  }
  linAlg->axmyMany(mesh->Nlocal, 3, nrs->fieldOffset, 0, 1.0, mesh->o_LMM, platform->o_slice0);
  return platform->o_slice0;
}

occa::memory scalarStrongSubCycleMovingMesh(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S)
{

  mesh_t* mesh = cds->meshV;
  linAlg_t* linAlg = mesh->linAlg;

  occa::memory& o_r1 = platform->o_slice2;
  occa::memory& o_r2 = platform->o_slice3;
  occa::memory& o_r3 = platform->o_slice4;
  occa::memory& o_r4 = platform->o_slice5;

  occa::memory& o_p0 = platform->o_slice0;
  occa::memory& o_u1 = platform->o_slice6;

  occa::memory& o_bmst = platform->o_slice1;

  // Solve for Each SubProblem
  for (int torder = (cds->nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = is * cds->fieldOffset +
                          torder * cds->NSfields * cds->fieldOffset;
    const dlong offset = torder * cds->fieldOffset;
    if (torder == cds->nEXT - 1){
      linAlg->axpby(cds->fieldOffset, cds->coeffBDF[torder],
                           o_S, 0.0, o_p0, toffset, 0);
      occa::memory o_LMM_slice = mesh->o_LMM + offset * sizeof(dfloat);
      linAlg->axmy(mesh->Nlocal, 1.0, o_LMM_slice, o_p0);
    }
    else{
      cds->subCycleExtrapolateFieldKernel(mesh->Nlocal, toffset, offset, cds->coeffBDF[torder], mesh->o_LMM,
        o_S, o_p0);
    }

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder]/cds->Nsubsteps;

    for(int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;
      o_u1.copyFrom(o_p0, mesh->Nlocal * sizeof(dfloat));
      for(int rk = 0; rk < cds->SNrk; ++rk)
      {
        occa::memory o_rhs;
        if(rk == 0) o_rhs = o_r1;
        if(rk == 1) o_rhs = o_r2;
        if(rk == 2) o_rhs = o_r3;
        if(rk == 3) o_rhs = o_r4;

        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * cds->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        switch(cds->nEXT) {
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
        cds->subCycleExtrapolateScalarKernel(mesh->Nlocal, cds->nEXT, cds->fieldOffset, cds->o_extC, mesh->o_LMM, o_bmst);
        linAlg->aydx(mesh->Nlocal, 1.0, o_bmst, o_u1);

        if(mesh->NglobalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              cds->vFieldOffset,
              0,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              o_u1,
              o_rhs);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              cds->vFieldOffset,
              0,
              mesh->o_vgeo,
              mesh->o_D,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              o_u1,
              o_rhs);
        }

        oogs::start(o_rhs, 1, cds->fieldOffset, ogsDfloat, ogsAdd, cds->gsh);

        if(mesh->NlocalGatherElements) {
          if(cds->options[is].compareArgs("ADVECTION TYPE", "CUBATURE"))
            cds->subCycleStrongCubatureVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              cds->vFieldOffset,
              0,
              mesh->o_vgeo,
              mesh->o_cubvgeo,
              mesh->o_cubDiffInterpT,
              mesh->o_cubInterpT,
              mesh->o_cubProjectT,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC, 
              o_U,
              o_u1,
              o_rhs);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              cds->vFieldOffset,
              0,
              mesh->o_vgeo,
              mesh->o_D,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              o_u1,
              o_rhs);
        }

        oogs::finish(o_rhs, 1, cds->fieldOffset, ogsDfloat, ogsAdd, cds->gsh);

        linAlg->axmy(mesh->Nlocal, 1.0, o_bmst, o_rhs);
        if(rk != 3) {
          linAlg->axpbyz(mesh->Nlocal, 1.0, o_p0, -sdt * cds->Srka[rk+1], o_rhs, o_u1);
        }
        else{
          cds->subCycleRKKernel(
            mesh->Nlocal,
            sdt,
            cds->o_Srkb,
            o_r1,
            o_r2,
            o_r3,
            o_r4,
            o_p0
          );
        }
      }
    }
  }
  return o_p0;
}
occa::memory scalarStrongSubCycle(cds_t* cds, dfloat time, int is,
                                  occa::memory o_U, occa::memory o_S)
{
  mesh_t* mesh = cds->meshV;
  linAlg_t* linAlg= mesh->linAlg;

  // Solve for Each SubProblem
  for (int torder = (cds->nEXT - 1); torder >= 0; torder--) {
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = is * cds->fieldOffset +
                          torder * cds->NSfields * cds->fieldOffset;
    if (torder == cds->nEXT - 1){
      platform->linAlg->axpby(
        cds->fieldOffset,
        cds->coeffBDF[torder],
        o_S,
        0.0,
        platform->o_slice0,
        toffset, 0
      );
    }
    else{
      platform->linAlg->axpby(
        cds->fieldOffset,
        cds->coeffBDF[torder],
        o_S,
        1.0,
        platform->o_slice0,
        toffset, 0
      );
    }

    // Advance SubProblem to t^(n-torder+1)
    dfloat tsub = time;
    for (int i = torder; i > 0 ; i--) tsub -= cds->dt[i];
    const dfloat sdt = cds->dt[torder]/cds->Nsubsteps;

    for(int ststep = 0; ststep < cds->Nsubsteps; ++ststep) {
      const dfloat tstage = tsub + ststep * sdt;

      platform->o_slice0.copyFrom(platform->o_slice0, cds->fieldOffset * sizeof(dfloat),
                           cds->fieldOffset * sizeof(dfloat), 0);

      for(int rk = 0; rk < cds->SNrk; ++rk) {
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  sdt * cds->Srkc[rk];
        const dfloat tn0 = time;
        const dfloat tn1 = time - cds->dt[1];
        const dfloat tn2 = time - (cds->dt[1] + cds->dt[2]);
        switch(cds->nEXT) {
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
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice2);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NglobalGatherElements,
              mesh->o_globalGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_D,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice2);
        }

        occa::memory o_rhs;
        if(rk == 0) o_rhs = platform->o_slice2;
        if(rk == 1) o_rhs = platform->o_slice3;
        if(rk == 2) o_rhs = platform->o_slice4;
        if(rk == 3) o_rhs = platform->o_slice5;

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
              mesh->o_BdivW,
              cds->o_extC, 
              o_U,
              platform->o_slice0,
              platform->o_slice2);
          else
            cds->subCycleStrongVolumeKernel(
              mesh->NlocalGatherElements,
              mesh->o_localGatherElementList,
              cds->vFieldOffset,
              rk * cds->fieldOffset,
              mesh->o_vgeo,
              mesh->o_D,
              mesh->o_invLMM,
              mesh->o_BdivW,
              cds->o_extC,
              o_U,
              platform->o_slice0,
              platform->o_slice2);
        }

        oogs::finish(o_rhs, 1, cds->fieldOffset, ogsDfloat, ogsAdd, cds->gsh);

        cds->subCycleRKUpdateKernel(
          mesh->Nelements,
          rk,
          sdt,
          cds->fieldOffset,
          cds->o_Srka,
          cds->o_Srkb,
          platform->o_slice1,
          platform->o_slice2,
          platform->o_slice0);
      }
    }
  }
  linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_LMM, platform->o_slice0);
  return platform->o_slice0;
}
void meshUpdate(nrs_t* nrs)
{
  mesh_t * mesh = nrs->meshT;
  cds_t * cds = nrs->cds;

  for (int s = nrs->nBDF; s > 1; s--) {
    const dlong NbyteScalar = nrs->fieldOffset * sizeof(dfloat);
    mesh->o_LMM.copyFrom (mesh->o_LMM , NbyteScalar, (s - 1)*NbyteScalar, (s - 2)*NbyteScalar);
    mesh->o_invLMM.copyFrom (mesh->o_invLMM , NbyteScalar, (s - 1)*NbyteScalar, (s - 2)*NbyteScalar);
  }

  mesh->move();

  if(nrs->meshV!= nrs->meshT) nrs->meshV->computeInvLMM();

  // lag mesh velocities
  for (int s = nrs->nBDF; s > 1; s--) {
    const dlong Nbyte = nrs->fieldOffset * nrs->NVfields * sizeof(dfloat);
    mesh->o_U.copyFrom (mesh->o_U , Nbyte, (s - 1)*Nbyte, (s - 2)*Nbyte);
  }
}

void meshUpdateBdivW(nrs_t* nrs)
{
  mesh_t * mesh = nrs->meshT;
  cds_t * cds = nrs->cds;

  for (int s = nrs->nBDF; s > 1; s--) {
    const dlong NbyteScalar = nrs->fieldOffset * sizeof(dfloat);
    mesh->o_BdivW.copyFrom (mesh->o_BdivW , NbyteScalar, (s - 1)*NbyteScalar, (s - 2)*NbyteScalar);
  }

  mesh->computeBdivW();

  if(nrs->meshV!= nrs->meshT) nrs->meshV->computeBdivW();
}
