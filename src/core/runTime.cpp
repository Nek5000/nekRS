#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp" 
#include "cfl.hpp"

void extbdfCoefficents(ins_t *ins, int order);

void makef(ins_t *ins, dfloat time, occa::memory &o_BF);
occa::memory velocityStrongSubCycle(ins_t *ins, dfloat time, 
                                    occa::memory o_U);
void fluidSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U);

void makeq(ins_t *ins, dfloat time, occa::memory o_BF);
occa::memory scalarStrongSubCycle(cds_t *cds, dfloat time, int is, 
                                  occa::memory o_U, occa::memory o_S); 
void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S);

void qthermal(ins_t *ins, dfloat time, occa::memory o_div);

double tElapsed = 0;


void runStep(ins_t *ins, dfloat time, dfloat dt, int tstep)
{
  mesh_t *mesh = ins->mesh;
  cds_t *cds = ins->cds;

  double tStart = MPI_Wtime();

  ins->dt = dt;
  if(tstep<=1){
    extbdfCoefficents(ins,tstep);
  } else if(tstep<=2 && ins->temporalOrder>=2){ 
    extbdfCoefficents(ins,tstep);
  } else if(tstep<=3 && ins->temporalOrder>=3){ 
    extbdfCoefficents(ins,tstep);
  }

  // First extrapolate velocity to t^(n+1)
  int velocityExtrapolationOrder = ins->ExplicitOrder;
  if(!ins->flow) velocityExtrapolationOrder = 1; 
  ins->velocityExtKernel(mesh->Nelements,
                         velocityExtrapolationOrder,
                         ins->fieldOffset,
                         ins->o_extbdfC,
                         ins->o_U,
                         ins->o_Ue);

  if(ins->Nscalar)
    scalarSolve(ins, time, dt, cds->o_S);

  timer::tic("udfProperties"); 
  if(udf.properties)
    udf.properties(ins, time+dt, ins->o_U, cds->o_S, ins->o_prop, cds->o_prop);
  timer::toc("udfProperties"); 

  if(ins->lowMach){
    if(udf.div)
      udf.div(ins, time+dt, ins->o_div);
    else
      qthermal(ins, time+dt, ins->o_div);
  }

  if(ins->flow) fluidSolve(ins, time, dt, ins->o_U);
    
  const dfloat cfl = computeCFL(ins, time+dt, tstep);

  if(mesh->rank==0) {
    printf("step= %d  t= %.8e  dt=%.1e  C= %.2f",
           tstep, time+dt, dt, cfl); 

    if(ins->flow) {
      if(ins->uvwSolver)
        printf("  UVW: %d  P: %d", ins->NiterU, ins->NiterP); 
      else
        printf("  U: %d  V: %d  W: %d  P: %d", ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP); 
    }

    for(int is=0; is<ins->Nscalar; is++) {
      if(cds->compute[is]) printf("  S: %d", cds->Niter[is]);
    }
    
    const double tElapsedStep = MPI_Wtime()-tStart;
    tElapsed += tElapsedStep;
    printf("  eTime= %.2e, %.5e s\n", tElapsedStep, tElapsed);
  }

  if(cfl > 30 || std::isnan(cfl)) {
    if(mesh->rank==0) cout << "Unreasonable CFL! Dying ...\n" << endl; 
    ABORT(1);
  }

  if(tstep%10==0) fflush(stdout);
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

void makeq(ins_t *ins, dfloat time, occa::memory o_BF)
{
  cds_t *cds   = ins->cds;
  mesh_t *mesh = cds->mesh;

  cds->setScalarKernel(cds->fieldOffset*cds->NSfields, 0.0, cds->o_FS);
  if(udf.sEqnSource) udf.sEqnSource(ins, time, cds->o_S, cds->o_FS);

  for(int is=0; is<cds->NSfields; is++) {
    if(!cds->compute[is]) continue;

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;
    const dlong sOffset = is*cds->fieldOffset;
    occa::memory o_adv = cds->o_wrk0;

    if(cds->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
      cds->filterRTKernel(
           cds->meshV->Nelements,
           ins->o_filterMT,
           ins->filterS,
           sOffset,
           cds->o_rho,
           cds->o_S,
           cds->o_FS);

    if(cds->Nsubsteps) {

      o_adv = scalarStrongSubCycle(cds, time, is, cds->o_U, cds->o_S);

    } else {

     if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")) {
        cds->advectionStrongCubatureVolumeKernel(
             cds->meshV->Nelements,
             mesh->o_vgeo,
             mesh->o_cubvgeo,
             mesh->o_cubDiffInterpT,
             mesh->o_cubInterpT,
             mesh->o_cubProjectT,
             cds->vFieldOffset,
             sOffset,
             cds->o_Ue,
             cds->o_S,
             cds->o_rho,
             cds->o_wrk0);
      } else {
        cds->advectionStrongVolumeKernel(
             cds->meshV->Nelements,
             mesh->o_vgeo,
             mesh->o_Dmatrices,
             cds->vFieldOffset,
             sOffset,
             cds->o_Ue,
             cds->o_S,
             cds->o_rho,
             cds->o_wrk0);
      }

      ins->scaledAddKernel(
           cds->meshV->Nelements*cds->meshV->Np, 
           -1.0, 
           0*cds->fieldOffset, 
           cds->o_wrk0, 
           1.0, 
           sOffset, 
           cds->o_FS);

    }

    cds->sumMakefKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_MM,
         cds->idt,
         cds->o_extbdfA,
         cds->o_extbdfB,
         cds->o_extbdfC,
         cds->fieldOffset*cds->NSfields,
         sOffset,
         cds->o_S,
         o_adv,
         cds->o_FS,
         cds->o_rho,
         cds->o_BF);
  }
}

void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S)
{
  cds_t *cds   = ins->cds;

  for (int s=cds->Nstages;s>1;s--) {
    cds->o_FS.copyFrom(
           cds->o_FS, 
           cds->fieldOffset*cds->NSfields*sizeof(dfloat), 
           (s-1)*cds->fieldOffset*cds->NSfields*sizeof(dfloat), 
           (s-2)*cds->fieldOffset*cds->NSfields*sizeof(dfloat));
  }

  timer::tic("makeq");
  makeq(ins, time, cds->o_BF);
  timer::toc("makeq");

  timer::tic("scalarSolve");
  for (int is=0; is<cds->NSfields; is++) {
    if(!cds->compute[is]) continue; 

    mesh_t *mesh;
    (is) ? mesh = cds->meshV : mesh = cds->mesh;

    cds->setEllipticCoeffKernel(
         cds->Nlocal,
         cds->g0*cds->idt,
         is*cds->fieldOffset,
         cds->fieldOffset,
         cds->o_diff,
         cds->o_rho,
         cds->o_ellipticCoeff);

    if(cds->o_BFDiag.ptr()) {
      cds->scaledAddKernel(
           cds->Nlocal, 
           1.0, 
           is*cds->fieldOffset, 
           cds->o_BFDiag, 
           1.0, 
           cds->fieldOffset, 
           cds->o_ellipticCoeff); 
    }

    occa::memory o_Snew = cdsSolve(is, cds, time+dt);

    for (int s=cds->Nstages;s>1;s--) {
      o_S.copyFrom(
        o_S, 
        cds->Ntotal*sizeof(dfloat), 
        ((s-1)*(cds->fieldOffset*cds->NSfields) + is*cds->fieldOffset)*sizeof(dfloat), 
        ((s-2)*(cds->fieldOffset*cds->NSfields) + is*cds->fieldOffset)*sizeof(dfloat));
    }
    o_Snew.copyTo(o_S, cds->Ntotal*sizeof(dfloat), is*cds->fieldOffset*sizeof(dfloat));
  }
  timer::toc("scalarSolve");
}

void makef(ins_t *ins, dfloat time, occa::memory &o_BF)
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

  occa::memory o_adv = ins->o_wrk0;
  if(ins->Nsubsteps) {
    o_adv = velocityStrongSubCycle(ins, time, ins->o_U);
  } else {
    if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
       ins->advectionStrongCubatureVolumeKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_cubvgeo,
         mesh->o_cubDiffInterpT,
         mesh->o_cubInterpT,
         mesh->o_cubProjectT,
         ins->fieldOffset,
         ins->o_U,
         ins->o_wrk0);
     else
       ins->advectionStrongVolumeKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_Dmatrices,
         ins->fieldOffset,
         ins->o_U,
         ins->o_wrk0);

    ins->scaledAddKernel(
       ins->NVfields*ins->fieldOffset, 
       -1.0, 
       0, 
       ins->o_wrk0, 
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
  makef(ins, time, ins->o_BF);
  timer::toc("makef");

  timer::tic("pressureSolve");
  ins->setEllipticCoeffPressureKernel(
    ins->Nlocal,
    ins->fieldOffset,
    ins->o_rho,
    ins->o_ellipticCoeff);
  occa::memory o_Pnew = tombo::pressureSolve(ins, time+dt); 
  ins->o_P.copyFrom(o_Pnew, ins->Ntotal*sizeof(dfloat));
  timer::toc("pressureSolve");

  timer::tic("velocitySolve");
  ins->setEllipticCoeffKernel(
    ins->Nlocal,
    ins->g0*ins->idt,
    0*ins->fieldOffset,
    ins->fieldOffset,
    ins->o_mue,
    ins->o_rho,
    ins->o_ellipticCoeff);

  occa::memory o_Unew = tombo::velocitySolve(ins, time+dt);
  for (int s=ins->Nstages;s>1;s--) {
    o_U.copyFrom(
      o_U, 
      ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
      (s-1)*ins->fieldOffset*ins->NVfields*sizeof(dfloat), 
      (s-2)*ins->fieldOffset*ins->NVfields*sizeof(dfloat));
  }
  o_U.copyFrom(o_Unew, ins->NVfields*ins->fieldOffset*sizeof(dfloat));
  timer::toc("velocitySolve");
}

occa::memory velocityStrongSubCycle(ins_t *ins, dfloat time, occa::memory o_U)
{
  mesh_t *mesh = ins->mesh;

  // Solve for Each SubProblem
  for (int torder=ins->ExplicitOrder-1; torder>=0; torder--){
    
    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*ins->NVfields*ins->fieldOffset;
    const dfloat b=ins->extbdfB[torder];
    if (torder==ins->ExplicitOrder-1) 
      ins->scaledAddKernel(ins->NVfields*ins->fieldOffset, b, toffset, 
                           o_U, 0.0, 0, ins->o_wrk0);
    else
      ins->scaledAddKernel(ins->NVfields*ins->fieldOffset, b, toffset, 
                           o_U, 1.0, 0, ins->o_wrk0);
    
    // Advance subproblem from here from t^(n-torder) to t^(n-torder+1)
    for(int ststep = 0; ststep<ins->Nsubsteps;++ststep){

      const dfloat tsub   = time - torder*ins->dt;
      const dfloat tstage = tsub + ststep*ins->sdt;     

      //ins->o_wrk3.copyFrom(ins->o_wrk0, ins->NVfields*ins->fieldOffset*sizeof(dfloat));   
      ins->o_wrk0.copyFrom(ins->o_wrk0, ins->NVfields*ins->fieldOffset*sizeof(dfloat),
                           ins->NVfields*ins->fieldOffset*sizeof(dfloat),0);   

      for(int rk=0;rk<ins->SNrk;++rk){
        
        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  ins->sdt*ins->Srkc[rk]; 
        const dfloat tn0 = time - 0*ins->dt;
        const dfloat tn1 = time - 1*ins->dt;
        const dfloat tn2 = time - 2*ins->dt;
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
        ins->velocityExtKernel(mesh->Nelements,
                               ins->Nstages,
                               ins->fieldOffset,
                               ins->o_extC,
                               o_U,
                               ins->o_Ue);

        // Compute Volume Contribution
	if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
	  ins->subCycleStrongCubatureVolumeKernel(
            mesh->Nelements,
            mesh->o_vgeo,
            mesh->o_cubvgeo,
            mesh->o_cubDiffInterpT, //mesh->o_cubDWmatrices,
            mesh->o_cubInterpT,
            mesh->o_cubProjectT,
            ins->o_InvM,
            ins->fieldOffset,
            rk*ins->NVfields*ins->fieldOffset,
            ins->o_Ue,
            ins->o_wrk0,
            ins->o_wrk6);
	else
	  ins->subCycleStrongVolumeKernel(mesh->Nelements,
	    mesh->o_vgeo,
	    mesh->o_Dmatrices,
	    ins->fieldOffset,
            rk*ins->NVfields*ins->fieldOffset,
	    ins->o_Ue,
	    ins->o_wrk0,
	    ins->o_wrk6);

        occa::memory o_rhs;
        if(rk==0) o_rhs = ins->o_wrk6;
        if(rk==1) o_rhs = ins->o_wrk9;
        if(rk==2) o_rhs = ins->o_wrk12;
        if(rk==3) o_rhs = ins->o_wrk15;
        ogsGatherScatterMany(o_rhs, ins->NVfields, ins->fieldOffset,
                             ogsDfloat, ogsAdd, mesh->ogs);
        ins->invMassMatrixKernel(
             mesh->Nelements,
	     ins->fieldOffset,
	     ins->NVfields,
	     mesh->o_vgeo,
	     ins->o_InvM,
	     o_rhs);
        
        ins->subCycleRKUpdateKernel(
             mesh->Nelements,
             rk,
             ins->sdt,
             ins->fieldOffset,
             ins->o_Srka,
             ins->o_Srkb,
             ins->o_wrk3,
             ins->o_wrk6, 
             ins->o_wrk0);
      }
    }
  }
  return ins->o_wrk0;
}

occa::memory scalarStrongSubCycle(cds_t *cds, dfloat time, int is, 
                                  occa::memory o_U, occa::memory o_S)
{
  mesh_t *mesh = cds->meshV;

   // Solve for Each SubProblem
  for (int torder=(cds->ExplicitOrder-1); torder>=0; torder--){

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    const dlong toffset = is*cds->fieldOffset + 
                          torder*cds->NSfields*cds->fieldOffset;
    if (torder==cds->ExplicitOrder-1)
      cds->scaledAddKernel(cds->fieldOffset, cds->extbdfB[torder], 
                           toffset, o_S, 0.0, 0, cds->o_wrk0);
    else
      cds->scaledAddKernel(cds->fieldOffset, cds->extbdfB[torder],
                           toffset, o_S, 1.0, 0, cds->o_wrk0);

    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<cds->Nsubsteps;++ststep){

      const dfloat tsub   = time - torder*cds->dt;
      const dfloat tstage = tsub + ststep*cds->sdt;     

      //cds->o_wrk1.copyFrom(cds->o_wrk0, cds->fieldOffset*sizeof(dfloat));
      cds->o_wrk0.copyFrom(cds->o_wrk0, cds->fieldOffset*sizeof(dfloat), 
                           cds->fieldOffset*sizeof(dfloat), 0);  

      for(int rk=0;rk<cds->SNrk;++rk){

        // Extrapolate velocity to subProblem stage time
        const dfloat t   = tstage +  cds->sdt*cds->Srkc[rk]; 
        const dfloat tn0 = time - 0*cds->dt;
        const dfloat tn1 = time - 1*cds->dt;
        const dfloat tn2 = time - 2*cds->dt;
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
        cds->velocityExtKernel(
             mesh->Nelements,
             cds->Nstages,
             cds->vFieldOffset,
             cds->o_extC,
             o_U,
             cds->o_Ue);

        // Compute Volume Contribution
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
          cds->subCycleStrongCubatureVolumeKernel(
               mesh->Nelements,
               cds->vFieldOffset,
               rk*cds->fieldOffset,            
               mesh->o_vgeo,
               mesh->o_cubvgeo,
               mesh->o_cubDiffInterpT,
               mesh->o_cubInterpT,
               mesh->o_cubProjectT,
               cds->o_Ue,
               cds->o_wrk0,
               cds->o_wrk2);
        else
          cds->subCycleStrongVolumeKernel(
               mesh->Nelements,
               cds->vFieldOffset,
               rk*cds->fieldOffset,           
               mesh->o_vgeo,
               mesh->o_Dmatrices,
               cds->o_Ue,
               cds->o_wrk0,
               cds->o_wrk2);

        occa::memory o_rhs;
        if(rk==0) o_rhs = cds->o_wrk2;
        if(rk==1) o_rhs = cds->o_wrk3;
        if(rk==2) o_rhs = cds->o_wrk4;
        if(rk==3) o_rhs = cds->o_wrk5;
        ogsGatherScatter(o_rhs, ogsDfloat, ogsAdd, mesh->ogs);
        cds->invMassMatrixKernel(
             mesh->Nelements,
             cds->fieldOffset,
             1,
             mesh->o_vgeo,
             cds->o_InvMV,
             o_rhs);

        cds->subCycleRKUpdateKernel(
             mesh->Nelements,
             rk,
             cds->sdt,
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

// qtl = 1/(rho*cp*T) * (div[k*grad[T] ] + qvol)
void qthermal(ins_t *ins, dfloat time, occa::memory o_div)
{
  cds_t *cds = ins->cds;
  mesh_t *mesh = ins->mesh;

  ins->gradientVolumeKernel(
       mesh->Nelements,
       mesh->o_vgeo,
       mesh->o_Dmatrices,
       ins->fieldOffset,
       cds->o_S,
       cds->o_wrk0);

  ogsGatherScatterMany(cds->o_wrk0, ins->NVfields, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        ins->NVfields,
        mesh->o_vgeo,
        ins->o_InvM, 
        cds->o_wrk0);

   if(udf.sEqnSource)
     udf.sEqnSource(ins, time, cds->o_S, cds->o_wrk3);
   else
     ins->setScalarKernel(mesh->Nelements*mesh->Np, 0.0, cds->o_wrk3);

   ins->qtlKernel(
        mesh->Nelements,
        mesh->o_vgeo,
        mesh->o_Dmatrices,
        ins->fieldOffset,
        cds->o_wrk0,
        cds->o_S,
        cds->o_diff,
        cds->o_rho,
        cds->o_wrk3,
        o_div);

   ogsGatherScatterMany(o_div, 1, ins->fieldOffset,
                       ogsDfloat, ogsAdd, mesh->ogs);

   ins->invMassMatrixKernel(
        mesh->Nelements,
        ins->fieldOffset,
        1,
        mesh->o_vgeo,
        ins->o_InvM, 
        o_div);
}
