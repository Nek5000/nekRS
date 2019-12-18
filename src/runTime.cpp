#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp" 
#include "cfl.hpp"

void extbdfCoefficents(ins_t *ins, int order);

void makef(ins_t *ins, dfloat time, occa::memory o_wrk, occa::memory o_BF);
void makeq(ins_t *ins, dfloat time, occa::memory o_wrk, occa::memory o_BF);

void fluidSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U);
void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S);

void velocityStrongSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_wrk,
                            occa::memory o_U, occa::memory o_Ud);
void scalarStrongSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_wrk, 
                          occa::memory o_U, occa::memory o_S, occa::memory o_Sd);

void qthermal(ins_t *ins, dfloat time, occa::memory o_qtl);

double etime0 = 0;


void runStep(ins_t *ins, dfloat time, dfloat dt, int tstep)
{
  mesh_t *mesh = ins->mesh;
  cds_t *cds = ins->cds;

  ins->dt = dt;
  if(tstep<=1){
    etime0 = MPI_Wtime();
    extbdfCoefficents(ins,tstep);
  } else if(tstep<=2 && ins->temporalOrder>=2){ 
    extbdfCoefficents(ins,tstep);
  } else if(tstep<=3 && ins->temporalOrder>=3){ 
    extbdfCoefficents(ins,tstep);
  }

  // First extrapolate velocity to t^(n+1)
  ins->velocityExtKernel(mesh->Nelements,
                         ins->ExplicitOrder,
                         ins->fieldOffset,
                         ins->o_extbdfC,
                         ins->o_U,
                         ins->o_Ue);

  timer::tic("scalarSolve"); 
  if(ins->Nscalar)
    scalarSolve(ins, time, dt, cds->o_S);
  timer::toc("scalarSolve"); 

  timer::tic("udfProperties"); 
  if(udf.properties)
    udf.properties(ins, time+dt, ins->o_U, cds->o_S, ins->o_prop, cds->o_prop);
  timer::toc("udfProperties"); 

  if(ins->lowMach){
    if(udf.qtl)
      udf.qtl(ins, time+dt, ins->o_qtl);
    else
      qthermal(ins, time+dt, ins->o_qtl);
  }

  timer::tic("fluidSolve"); 
  fluidSolve(ins, time, dt, ins->o_U);
  timer::toc("fluidSolve"); 

  const dfloat cfl = computeCFL(ins, time+dt, tstep);

  if(mesh->rank==0) {
    printf("step= %d  t= %.5e  dt=%.1e  C= %.2f  U: %d  V: %d  W: %d  P: %d",
           tstep, time+dt, dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, 
           ins->NiterP);

    for(int is=0; is<ins->Nscalar; is++) {
      if(cds->compute[is]) printf("  S: %d", cds->Niter[is]);
    }

    printf("  tElapsed= %.5e s\n", MPI_Wtime()-etime0);
  }

  if(cfl > 20) {
    if(mesh->rank==0) cout << "CFL too high! Dying ...\n" << endl; 
    EXIT(1);
  }

  if(tstep%5==0) fflush(stdout);
}

void extbdfCoefficents(ins_t *ins, int order) {

  if(order==1) {

    ins->g0 =  1.0f;
    dfloat extbdfB[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfA[3] = {1.0f, 0.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder = 1;

  } else if(order==2) {

    ins->g0 =  1.5f;
    dfloat extbdfB[3] = {2.0f,-0.5f, 0.0f};
    dfloat extbdfA[3] = {2.0f,-1.0f, 0.0f};
    dfloat extbdfC[3] = {1.0f, 0.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=2;

  } else if(order==3) {

    ins->g0 =  11.f/6.f;
    dfloat extbdfB[3] = {3.0f,-1.5f, 1.0f/3.0f};
    dfloat extbdfA[3] = {3.0f,-3.0f, 1.0f};
    dfloat extbdfC[3] = {2.0f,-1.0f, 0.0f};

    memcpy(ins->extbdfB, extbdfB, 3*sizeof(dfloat));
    memcpy(ins->extbdfA, extbdfA, 3*sizeof(dfloat));
    memcpy(ins->extbdfC, extbdfC, 3*sizeof(dfloat));

    ins->o_extbdfB.copyFrom(extbdfB);
    ins->o_extbdfA.copyFrom(extbdfA);
    ins->o_extbdfC.copyFrom(extbdfC);

    ins->ExplicitOrder=3;

  }

  ins->ig0 = 1.0/ins->g0;

  if (ins->Nscalar) {
    ins->cds->ExplicitOrder = ins->ExplicitOrder;  
    ins->cds->g0 = ins->g0;    
    ins->cds->ig0 = ins->ig0; 
  }
}

void makeq(ins_t *ins, dfloat time, occa::memory o_scratch6, occa::memory o_BF)
{
  timer::tic("makeq");

  cds_t *cds   = ins->cds;

  occa::memory o_adv = o_scratch6.slice(0*1*cds->fieldOffset*sizeof(dfloat));
  occa::memory o_wrk5 = o_scratch6.slice(1*cds->fieldOffset*sizeof(dfloat));

  cds->setScalarKernel(cds->fieldOffset*cds->NSfields, 0.0, cds->o_FS);

  if(udf.sEqnSource) udf.sEqnSource(ins, time, cds->o_S, cds->o_FS);
    
  for(int is=0; is<cds->NSfields; is++) {
    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;

    if(!cds->compute[is]) continue;

    occa::memory o_rho  = cds->o_rho.slice (is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_diff = cds->o_diff.slice(is*cds->fieldOffset*sizeof(dfloat));

    occa::memory o_Si  = cds->o_S.slice (is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_FSi = cds->o_FS.slice(is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_BFi = cds->o_BF.slice(is*cds->fieldOffset*sizeof(dfloat));

    if(cds->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
      cds->filterRTKernel(
           cds->meshV->Nelements,
           ins->o_filterMT,
           ins->filterS,
           0*cds->fieldOffset,
           o_rho,
           o_Si,
           o_FSi);

    cds->setScalarKernel(cds->fieldOffset, 0.0, o_adv);
    if(cds->Nsubsteps) {
      scalarStrongSubCycle(cds, time, cds->Nstages, o_wrk5, cds->o_U, o_Si, o_adv);
    } else {
      if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        cds->advectionStrongCubatureVolumeKernel(
               cds->meshV->Nelements,
               mesh->o_vgeo,
               mesh->o_cubvgeo,
               mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
               mesh->o_cubInterpT,
               mesh->o_cubProjectT,
               cds->vFieldOffset,
               cds->fieldOffset,
               cds->o_Ue,
               o_Si,
               o_rho,
               o_adv);
      } else {
        cds->advectionStrongVolumeKernel(
               cds->meshV->Nelements,
               mesh->o_vgeo,
               mesh->o_Dmatrices,
               cds->vFieldOffset,
               cds->fieldOffset,
               cds->o_Ue,
               o_Si,
               o_rho,
               o_adv);
      }
      ins->scaledAddKernel(
           cds->meshV->Nelements*cds->meshV->Np, 
           -1.0, 
           0*cds->fieldOffset, 
           o_adv, 
           1.0, 
           0*cds->fieldOffset, 
           o_FSi);
    }
    cds->sumMakefKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_MM,
         cds->idt,
         cds->o_extbdfA,
         cds->o_extbdfB,
         cds->o_extbdfC,
         cds->fieldOffset,
         o_Si,
         o_adv,
         o_FSi,
         o_rho,
         o_BFi);
  }
  timer::toc("makeq");
}

void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S)
{
  cds_t *cds   = ins->cds;

  occa::memory o_scratch = ins->o_scratch.slice(ins->ellipticWrkSize*sizeof(dfloat)); 
  occa::memory o_scratch6 = o_scratch.slice(0*6*cds->fieldOffset*sizeof(dfloat));

  for (int s=cds->Nstages;s>1;s--) {
    cds->o_FS.copyFrom(
           cds->o_FS, 
           cds->fieldOffset*cds->NSfields*sizeof(dfloat), 
           (s-1)*cds->fieldOffset*cds->NSfields*sizeof(dfloat), 
           (s-2)*cds->fieldOffset*cds->NSfields*sizeof(dfloat));
  }

  makeq(ins, time, o_scratch6, cds->o_BF);

  for (int is=0; is<cds->NSfields; is++) {
    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;

    if(!cds->compute[is]) continue; 

    occa::memory o_Si   = cds->o_S.slice   (is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_rho  = cds->o_rho.slice (is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_diff = cds->o_diff.slice(is*cds->fieldOffset*sizeof(dfloat));
    occa::memory o_Snew = o_scratch.slice(6*cds->fieldOffset*sizeof(dfloat));

    cds->setEllipticCoeffKernel(
         cds->Nlocal,
         cds->g0*cds->idt,
         cds->solver[is]->Ntotal, // offset required by elliptic
         o_diff,
         o_rho,
         cds->o_ellipticCoeff);

    cdsSolve(is, cds, time+dt, o_scratch6, o_Snew);

    for (int s=cds->Nstages;s>1;s--) {
      o_S.copyFrom(
        o_S, 
        cds->fieldOffset*sizeof(dfloat), 
        ((s-1)*(cds->fieldOffset*cds->NSfields) + is*cds->fieldOffset)*sizeof(dfloat), 
        ((s-2)*(cds->fieldOffset*cds->NSfields) + is*cds->fieldOffset)*sizeof(dfloat));
    }
    o_Si.copyFrom(o_Snew, cds->Ntotal*sizeof(dfloat));
  }
}

void makef(ins_t *ins, dfloat time, occa::memory o_scratch18, occa::memory o_BF)
{
  mesh_t *mesh = ins->mesh;

  for (int s=ins->Nstages;s>1;s--) {
    ins->o_FU.copyFrom(
      ins->o_FU, 
      ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
      (s-1)*ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
      (s-2)*ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  }

  ins->setScalarKernel(ins->fieldOffset*ins->NVfields, 0.0, ins->o_FU);
  if(udf.uEqnSource) udf.uEqnSource(ins, time, ins->o_U, ins->o_FU);
  
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    ins->filterRTKernel(
         mesh->Nelements,
         ins->o_filterMT,
         ins->filterS,
         ins->fieldOffset,
         ins->o_U,
         ins->o_FU);


  occa::memory o_adv   = o_scratch18.slice(0*3*ins->fieldOffset*sizeof(dfloat));;
  occa::memory o_wrk15 = o_scratch18.slice(3*ins->fieldOffset*sizeof(dfloat));

  if(ins->Nsubsteps) {
    velocityStrongSubCycle(ins, time, ins->Nstages, o_wrk15, ins->o_U, o_adv);
  } else {
    if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
       ins->advectionStrongCubatureVolumeKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_cubvgeo,
         mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
         mesh->o_cubInterpT,
         mesh->o_cubProjectT,
         ins->fieldOffset,
         ins->o_U,
         o_adv);
     else
       ins->advectionStrongVolumeKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_Dmatrices,
         ins->fieldOffset,
         ins->o_U,
         o_adv);

    ins->scaledAddKernel(
       ins->NVfields*ins->fieldOffset, 
       -1.0, 
       0, 
       o_adv, 
       1.0, 
       0, 
       ins->o_FU);
  }

  ins->sumMakefKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_MM,
       ins->idt,
       ins->o_extbdfA,
       ins->o_extbdfB,
       ins->fieldOffset,
       ins->o_U,
       o_adv,
       ins->o_FU,
       o_BF);
}

void fluidSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U)
{
  mesh_t *mesh = ins->mesh;

  timer::tic("makef");
  makef(ins, time, ins->o_scratch, ins->o_BF);
  timer::toc("makef");

  {
    ins->setEllipticCoeffPressureKernel(
      ins->Nlocal,
      ins->pSolver->Ntotal, // offset required by elliptic
      ins->o_rho,
      ins->o_ellipticCoeff);

    occa::memory o_scratch = ins->o_scratch.slice(ins->ellipticWrkSize*sizeof(dfloat)); 
    occa::memory o_Pnew = o_scratch.slice(0*1*ins->fieldOffset*sizeof(dfloat));
    occa::memory o_wrk6 = o_scratch.slice(1*ins->fieldOffset*sizeof(dfloat));
 
    timer::tic("pressureSolve");
    tombo::pressureSolve(ins, time+dt, o_wrk6, o_Pnew);
    timer::toc("pressureSolve");
    ins->o_P.copyFrom(o_Pnew, ins->Ntotal*sizeof(dfloat));
  }

  {
    ins->setEllipticCoeffKernel(
      ins->Nlocal,
      ins->g0*ins->idt,
      ins->uSolver->Ntotal, // offset required by elliptic
      ins->o_mue,
      ins->o_rho,
      ins->o_ellipticCoeff);

    occa::memory o_scratch = ins->o_scratch.slice(ins->ellipticWrkSize*sizeof(dfloat)); 
    occa::memory o_Unew = o_scratch.slice(0*3*ins->fieldOffset*sizeof(dfloat));
    occa::memory o_wrk6 = o_scratch.slice(3*ins->fieldOffset*sizeof(dfloat));

    timer::tic("velocitySolve");
    tombo::velocitySolve(ins, time+dt, o_wrk6, o_Unew);
    timer::toc("velocitySolve");
    for (int s=ins->Nstages;s>1;s--) {
      o_U.copyFrom(
        o_U, 
        ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
        (s-1)*ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
        (s-2)*ins->fieldOffset*ins->NVfields*sizeof(dfloat));
    }
    o_U.copyFrom(o_Unew, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
  } 
}

void velocityStrongSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_scratch15,
                            occa::memory o_U, occa::memory o_Ud)
{
  mesh_t *mesh = ins->mesh;

  occa::memory o_tmp = o_scratch15.slice(0*ins->NVfields*ins->fieldOffset * sizeof(dfloat));

  const dfloat tn0 = time - 0*ins->dt;
  const dfloat tn1 = time - 1*ins->dt;
  const dfloat tn2 = time - 2*ins->dt;

  dfloat zero = 0.0, one = 1.0;
  int izero = 0;

  dfloat b, bScale=0;

  // Solve for Each SubProblem
  for (int torder=ins->ExplicitOrder-1; torder>=0; torder--){
    
    b=ins->extbdfB[torder];
    bScale += b;

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*ins->NVfields*ins->fieldOffset;
    if (torder==ins->ExplicitOrder-1) { //first substep
      ins->scaledAddKernel(ins->NVfields*ins->fieldOffset, b, toffset, o_U, zero, izero, o_Ud);
    } else { //add the next field
      ins->scaledAddKernel(ins->NVfields*ins->fieldOffset, b, toffset, o_U,  one, izero, o_Ud);
    }   
    
    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    const dfloat tsub = time - torder*ins->dt;
    for(int ststep = 0; ststep<ins->Nsubsteps;++ststep){
      const dfloat tstage = tsub + ststep*ins->sdt;     

      if(ins->SNrk==4) o_tmp.copyFrom(o_Ud, ins->NVfields*ins->fieldOffset*sizeof(dfloat));   

      for(int rk=0;rk<ins->SNrk;++rk){// LSERK4 stages
        // Extrapolate velocity to subProblem stage time
        dfloat t = tstage +  ins->sdt*ins->Srkc[rk]; 

        switch(ins->ExplicitOrder){
	case 1:
	  ins->extC[0] = 1.f; ins->extC[1] = 0.f; ins->extC[2] = 0.f;
	  break;
	case 2:
	  ins->extC[0] = (t-tn1)/(tn0-tn1);
	  ins->extC[1] = (t-tn0)/(tn1-tn0);
	  ins->extC[2] = 0.f; 
	  break;
	case 3:
	  ins->extC[0] = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
	  ins->extC[1] = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
	  ins->extC[2] = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
	  break;
        }
        ins->o_extC.copyFrom(ins->extC);

        //compute advective velocity fields at time t
        ins->velocityExtKernel(mesh->Nelements,
                               Nstages,
                               ins->fieldOffset,
                               ins->o_extC,
                               o_U,
                               ins->o_Ue);

        occa::memory o_rhsU = o_scratch15.slice(ins->NVfields*ins->fieldOffset * sizeof(dfloat));
        if(ins->SNrk==4)
          o_rhsU = o_scratch15.slice((1 + rk)*ins->NVfields*ins->fieldOffset * sizeof(dfloat));

        // Compute Volume Contribution
	if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){        
	  ins->subCycleStrongCubatureVolumeKernel(
            mesh->Nelements,
            mesh->o_vgeo,
            mesh->o_cubvgeo,
            mesh->o_cubDiffInterpT, //mesh->o_cubDWmatrices,
            mesh->o_cubInterpT,
            mesh->o_cubProjectT,
            ins->o_InvM,
            ins->fieldOffset,
            ins->o_Ue,
            o_Ud,
            o_rhsU);
	}else{
	  ins->subCycleStrongVolumeKernel(mesh->Nelements,
	    mesh->o_vgeo,
	    mesh->o_Dmatrices,
	    ins->fieldOffset,
	    ins->o_Ue,
	    o_Ud,
	    o_rhsU);
	}
	
        ogsGatherScatterMany(o_rhsU, ins->NVfields, ins->fieldOffset,
                             ogsDfloat, ogsAdd, mesh->ogs);
  
        ins->invMassMatrixKernel(mesh->Nelements,
				 ins->fieldOffset,
				 ins->NVfields,
				 mesh->o_vgeo,
				 ins->o_InvM, // mesh->o_MM, 
				 o_rhsU);
        
        if(ins->SNrk==5){ // LSERK
          ins->subCycleRKUpdateKernel(
               mesh->Nelements,
  	       ins->sdt,
  	       ins->Srka[rk],
  	       ins->Srkb[rk],
  	       ins->fieldOffset,
               o_rhsU,
  	       o_tmp, 
  	       o_Ud);
        }else{ // ERK
          if(rk==ins->SNrk-1) o_rhsU = o_scratch15.slice(ins->NVfields*ins->fieldOffset * sizeof(dfloat));
          ins->subCycleRKUpdateKernel(
               mesh->Nelements,
               rk,
               ins->sdt,
               ins->fieldOffset,
               ins->o_Srka,
               ins->o_Srkb,
               o_tmp,
               o_rhsU, 
               o_Ud);
        }
      }
    }
  }
}

void scalarStrongSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_scratch5, 
                          occa::memory o_U, occa::memory o_S, occa::memory o_Sd)
{
  mesh_t *mesh = cds->meshV;

  occa::memory o_tmp = o_scratch5.slice(0*cds->fieldOffset * sizeof(dfloat));

  const dfloat tn0 = time - 0*cds->dt;
  const dfloat tn1 = time - 1*cds->dt;
  const dfloat tn2 = time - 2*cds->dt;

  dfloat zero = 0.0, one = 1.0;
  int izero = 0;
  dfloat b, bScale=0;

   // Solve for Each SubProblem
  for (int torder=(cds->ExplicitOrder-1); torder>=0; torder--){
    
    b=cds->extbdfB[torder];
    bScale += b; 

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*cds->NSfields*cds->fieldOffset;

    if (torder==cds->ExplicitOrder-1) { //first substep
      cds->scaledAddKernel(cds->fieldOffset, b, toffset, o_S, zero, izero, o_Sd);
    } else { //add the next field
      cds->scaledAddKernel(cds->fieldOffset, b, toffset, o_S,  one, izero, o_Sd);
    } 

    // SubProblem  starts from here from t^(n-torder)
    const dfloat tsub = time - torder*cds->dt;
    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<cds->Nsubsteps;++ststep){

      if(cds->SNrk==4) o_tmp.copyFrom(o_Sd, cds->fieldOffset*sizeof(dfloat));   
     
      const dfloat tstage = tsub + ststep*cds->sdt;     
     
      for(int rk=0;rk<cds->SNrk;++rk){// loop over stages
        // Extrapolate velocity to subProblem stage time
        dfloat t = tstage +  cds->sdt*cds->Srkc[rk]; 

        switch(cds->ExplicitOrder){
          case 1:
            cds->extC[0] = 1.f; cds->extC[1] = 0.f; cds->extC[2] = 0.f;
            break;
          case 2:
            cds->extC[0] = (t-tn1)/(tn0-tn1);
            cds->extC[1] = (t-tn0)/(tn1-tn0);
            cds->extC[2] = 0.f; 
            break;
          case 3:
            cds->extC[0] = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
            cds->extC[1] = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
            cds->extC[2] = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
            break;
        }
        cds->o_extC.copyFrom(cds->extC);

        //compute advective velocity fields at time t
        cds->velocityExtKernel(mesh->Nelements,
                               Nstages,
                               cds->vFieldOffset,
                               cds->o_extC,
                               o_U,
                               cds->o_Ue);

        occa::memory o_rhsS = o_scratch5.slice(1*cds->fieldOffset * sizeof(dfloat));
        if(cds->SNrk==4)
          o_rhsS = o_scratch5.slice((1+rk)*cds->fieldOffset * sizeof(dfloat));

        // Compute Volume Contribution
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          cds->subCycleStrongCubatureVolumeKernel(
                                            mesh->Nelements,
                                            mesh->o_vgeo,
                                            mesh->o_cubvgeo,
                                            mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
                                            mesh->o_cubInterpT,
                                            mesh->o_cubProjectT,
                                            cds->vFieldOffset,
                                            cds->fieldOffset,             
                                            cds->o_Ue,
                                            o_Sd,
                                            o_rhsS);
        } else{
          cds->subCycleStrongVolumeKernel(mesh->Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_Dmatrices,
                                    cds->vFieldOffset,
                                    cds->fieldOffset,           
                                    cds->o_Ue,
                                    o_Sd,
                                    o_rhsS);

        }

        ogsGatherScatter(o_rhsS, ogsDfloat, ogsAdd, mesh->ogs);

        // int nfield = ins->dim==2 ? 2:3; 
        cds->invMassMatrixKernel(mesh->Nelements,
                                 cds->fieldOffset,
                                 1,
                                 mesh->o_vgeo,
                                 cds->o_InvMV,
                                 o_rhsS);

        // Update Kernel
        if(cds->SNrk==5){
          cds->subCycleRKUpdateKernel(
               mesh->Nelements,
               cds->sdt,
               cds->Srka[rk],
               cds->Srkb[rk],
               cds->fieldOffset,
               o_rhsS,
               o_tmp, 
               o_Sd);
        }else{
          if(rk==cds->SNrk-1) o_rhsS = o_scratch5.slice(1*cds->fieldOffset * sizeof(dfloat));
          cds->subCycleRKUpdateKernel(
               mesh->Nelements,
               rk,
               cds->sdt,
               cds->fieldOffset,
               cds->o_Srka,
               cds->o_Srkb,
               o_tmp,
               o_rhsS, 
               o_Sd);
        }
      }
    }
  }
}


// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void qthermal(ins_t *ins, dfloat time, occa::memory o_qtl)
{
  cds_t *cds = ins->cds;
  mesh_t *mesh = ins->mesh;

  occa::memory o_gradS = ins->o_scratch.slice(0*3*cds->fieldOffset*sizeof(dfloat));
  occa::memory o_src   = ins->o_scratch.slice(3*cds->fieldOffset*sizeof(dfloat));

  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       cds->o_S,
       o_gradS);

  ogsGatherScatterMany(o_gradS, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        ins->NVfields,
        mesh->o_vgeo,
        ins->o_InvM, 
        o_gradS);

   if(udf.sEqnSource)
     udf.sEqnSource(ins, time, cds->o_S, o_src);
   else
     ins->setScalarKernel(mesh->Nelements*mesh->Np, 0.0, o_src);

   ins->qtlKernel(
        mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_Dmatrices,
        ins->fieldOffset,
        o_gradS,
        cds->o_S,
        cds->o_diff,
        cds->o_rho,
        o_src,
        o_qtl);

   ogsGatherScatterMany(o_qtl, 1, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        1,
        mesh->o_vgeo,
        ins->o_InvM, 
        o_qtl);
}

