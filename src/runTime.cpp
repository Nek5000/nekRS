#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nekrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp" 
#include "cfl.hpp"

void extbdfCoefficents(ins_t *ins, int order);
void makeq(ins_t *ins, dfloat time, occa::memory o_NS, occa::memory o_FS);
void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S);
void makef(ins_t *ins, dfloat time, occa::memory o_NU, occa::memory o_FU);
void velocitySolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U);

void runTime(ins_t *ins)
{
  mesh_t *mesh = ins->mesh;
  cds_t *cds = ins->cds; 

  double etime0 = MPI_Wtime();
  for(int tstep=0;tstep<ins->NtimeSteps;++tstep){

    if(tstep<1) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<2 && ins->temporalOrder>=2) 
      extbdfCoefficents(ins,tstep+1);
    else if(tstep<3 && ins->temporalOrder>=3) 
      extbdfCoefficents(ins,tstep+1);

    dfloat time = ins->startTime + tstep*ins->dt;

    ins->isOutputStep = 0;
    nek_ifoutfld(0);
    if (ins->outputStep > 0) {
      if (((tstep+1)%(ins->outputStep))==0 ||  tstep+1 == ins->NtimeSteps) {
        ins->isOutputStep = 1;
        nek_ifoutfld(1);
      }
    }
    
    if(ins->Nscalar) scalarSolve(ins, time, ins->dt, cds->o_S); 
    velocitySolve(ins, time, ins->dt, ins->o_U);

    dfloat cfl = computeCFL(ins, time+ins->dt, tstep+1);
    if (cfl > 20) {
      if (mesh->rank==0) cout << "CFL too high! Dying ...\n" << endl; 
      EXIT(0);
    }

    if (mesh->rank==0) {
      if(ins->Nscalar)
        printf("step= %d  t= %.5e  dt=%.1e  C= %.2f  U: %d  V: %d  W: %d  P: %d  S: %d  tElapsed= %.5e s\n",
          tstep+1, time+ins->dt, ins->dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, 
          ins->NiterP, cds->Niter, MPI_Wtime()-etime0);
      else
        printf("step= %d  t= %.5e  dt=%.1e  C= %.2f  U: %d  V: %d  W: %d  P: %d  tElapsed= %.5e s\n",
          tstep+1, time+ins->dt, ins->dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, 
          ins->NiterP, MPI_Wtime()-etime0);

      if ((tstep+1)%5==0) fflush(stdout);
    }

    if (ins->isOutputStep) nek_ocopyFrom(ins, time+ins->dt, tstep+1); 

    if (udf.executeStep) udf.executeStep(ins, time+ins->dt, tstep+1);

    if (ins->isOutputStep) nek_outfld(); 

  }
}

void extbdfCoefficents(ins_t *ins, int order) {

  if(order==1) {
     //advection, first order in time, increment
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

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;

     if(ins->Nscalar){
       ins->cds->ExplicitOrder = ins->ExplicitOrder;  
       ins->cds->g0 = ins->g0;    
       ins->cds->lambda = ins->cds->g0 / (ins->dt * ins->cds->diff);
       ins->cds->ig0 = 1.0/ins->cds->g0; 
    }

  } else if(order==2) {
    //advection, second order in time, increment
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

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;

     if(ins->Nscalar){
      ins->cds->ExplicitOrder = ins->ExplicitOrder=2;  
      ins->cds->g0 = ins->g0;  
      ins->cds->lambda = ins->cds->g0 / (ins->cds->dt * ins->cds->diff);
      ins->cds->ig0 = 1.0/ins->cds->g0; 
    }
  } else if(order==3) {
    //advection, third order in time, increment
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

    ins->lambda = ins->g0 / (ins->dt * ins->nu);
    ins->ig0 = 1.0/ins->g0;

    if(ins->Nscalar){
      ins->cds->ExplicitOrder = ins->ExplicitOrder=3;  
      ins->cds->g0 = ins->g0;  
      ins->cds->lambda = ins->cds->g0 / (ins->cds->dt * ins->cds->diff);
      ins->cds->ig0 = 1.0/ins->cds->g0; 
    }
  }
}

void makeq(ins_t *ins, dfloat time, occa::memory o_NS, occa::memory o_FS){
  cds_t *cds   = ins->cds; 
  mesh_t *mesh = ins->mesh;

  if(cds->Nsubsteps) {
    cdsStrongSubCycle(cds, time, cds->Nstages, ins->o_U, cds->o_S, o_NS);
  } else {
    // First extrapolate velocity to t^(n+1)
     cds->subCycleExtKernel(mesh->Nelements,
                            cds->ExplicitOrder,
                            cds->vOffset,
                            cds->o_extbdfC,
                            ins->o_U,
                            cds->o_Ue);
     
    cdsAdvection(cds, time, cds->o_Ue, cds->o_S, o_NS);
  }

  ins->setScalarKernel(cds->Ntotal*cds->NSfields, 0.0, o_FS);
  
  if(udf.sEqnSource)
    udf.sEqnSource(ins, time, cds->o_S, o_FS);

  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    ins->SFilterKernel(mesh->Nelements,
                       ins->o_filterMT,
                       ins->filterS,
                       cds->sOffset,
                       cds->o_S,
                       o_FS);
}

void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S){

  cds_t *cds   = ins->cds;
  mesh_t *mesh = cds->mesh;
  
  hlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  makeq(ins, time, cds->o_NS, cds->o_FS); 

  cdsHelmholtzRhs(cds, time+dt, cds->Nstages, cds->o_rhsS);

  cdsHelmholtzSolve(cds, time+dt, cds->Nstages, cds->o_rhsS, cds->o_rkS);

  for (int s=cds->Nstages;s>1;s--) {
    o_S.copyFrom( o_S, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
      (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
      (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));
  }

  //copy updated scalar
  o_S.copyFrom(cds->o_rkS, cds->NSfields*cds->Ntotal*sizeof(dfloat)); 
   
  for (int s=cds->Nstages;s>1;s--) {
    cds->o_NS.copyFrom(cds->o_NS, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
           (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
           (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));

    cds->o_FS.copyFrom(cds->o_FS, cds->Ntotal*cds->NSfields*sizeof(dfloat), 
           (s-1)*cds->Ntotal*cds->NSfields*sizeof(dfloat), 
           (s-2)*cds->Ntotal*cds->NSfields*sizeof(dfloat));
  }
}

void makef(ins_t *ins, dfloat time, occa::memory o_NU, occa::memory o_FU)
{
  mesh_t *mesh = ins->mesh;

  if(ins->Nsubsteps) {
    subCycle(ins, time, ins->Nstages, ins->o_U, o_NU);
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
         ins->o_cU,
         o_NU);
     else
       ins->advectionStrongVolumeKernel(
         mesh->Nelements,
         mesh->o_vgeo,
         mesh->o_Dmatrices,
         ins->fieldOffset,
         ins->o_U,
         o_NU);
  }

  ins->setScalarKernel(ins->Ntotal*ins->NVfields, 0.0, o_FU);
  if(udf.uEqnSource) udf.uEqnSource(ins, time, ins->o_U, o_FU);
  
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    ins->VFilterKernel(mesh->Nelements,
                       ins->o_filterMT,
                       ins->filterS,
                       ins->fieldOffset,
                       ins->o_U,
                       o_FU);
}

void velocitySolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U)
{
  makef(ins, time, ins->o_NU, ins->o_FU);

  curlCurl(ins, time, o_U, ins->o_NC);

  // o_rkU = sum_j (alpha_j U_j)- sum_j ( beta_j (FU_J + NU_j + NC_j) )
  pressureRhs  (ins, time+dt);
  pressureSolve(ins, time+dt, ins->o_rkP); 
  ins->o_P.copyFrom(ins->o_rkP, ins->Ntotal*sizeof(dfloat)); 

  // rhsU^s = MM*1/nu*[ -(grad P) + sum_i ( (a_i) U^(n-i)/dt - b_i (NU+FU)^(n-i) )]
  velocityRhs  (ins, time+dt);
  velocitySolve(ins, time+dt, ins->o_rkU);

  for (int s=ins->Nstages;s>1;s--)
    o_U.copyFrom(o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));

  o_U.copyFrom(ins->o_rkU, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

  for (int s=ins->Nstages;s>1;s--) {
    ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    
    ins->o_NC.copyFrom(ins->o_NC, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    
    ins->o_FU.copyFrom(ins->o_FU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
  }
}
