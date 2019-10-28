#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"
#include "udf.hpp"
#include "tombo.hpp" 
#include "cfl.hpp"

void extbdfCoefficents(ins_t *ins, int order);
void makeq(ins_t *ins, dfloat time, occa::memory o_NS, occa::memory o_FS);
void scalarSolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_S);
void makef(ins_t *ins, dfloat time, occa::memory o_NU, occa::memory o_FU);
void velocitySolve(ins_t *ins, dfloat time, dfloat dt, occa::memory o_U);
void velocityStrongSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_Ud);
void scalarStrongSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_U, occa::memory o_S, occa::memory o_Sd);

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
  ins->subCycleExtKernel(mesh->Nelements,
                         ins->ExplicitOrder,
                         ins->fieldOffset,
                         ins->o_extbdfC,
                         ins->o_U,
                         ins->o_Ue);

  if(ins->Nscalar) scalarSolve(ins, time, dt, cds->o_S); 
  velocitySolve(ins, time, dt, ins->o_U);

  dfloat cfl = computeCFL(ins, time+dt, tstep);

  if (mesh->rank==0) {
    if(ins->Nscalar)
      printf("step= %d  t= %.5e  dt=%.1e  C= %.2f  U: %d  V: %d  W: %d  P: %d  S: %d  tElapsed= %.5e s\n",
        tstep, time+dt, dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, 
        ins->NiterP, cds->Niter, MPI_Wtime()-etime0);
    else
      printf("step= %d  t= %.5e  dt=%.1e  C= %.2f  U: %d  V: %d  W: %d  P: %d  tElapsed= %.5e s\n",
        tstep, time+dt, dt, cfl, ins->NiterU, ins->NiterV, ins->NiterW, 
        ins->NiterP, MPI_Wtime()-etime0);
  }

  if (cfl > 20) {
    if (mesh->rank==0) cout << "CFL too high! Dying ...\n" << endl; 
    EXIT(0);
  }

  if (tstep%5==0) fflush(stdout);
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

  if(cds->Nsubsteps)
    scalarStrongSubCycle(cds, time, cds->Nstages, cds->o_U, cds->o_S, o_NS);
  else
    cdsAdvection(cds, time, cds->o_Ue, cds->o_S, o_NS);

  ins->setScalarKernel(cds->Ntotal*cds->NSfields, 0.0, o_FS);
  
  if(udf.sEqnSource)
    udf.sEqnSource(ins, time, cds->o_S, o_FS);

  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
    ins->SFilterKernel(cds->meshV->Nelements,
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
    velocityStrongSubCycle(ins, time, ins->Nstages, ins->o_U, o_NU);
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

  tombo::curlCurl(ins, time, ins->o_Ue, ins->o_NC);

  occa::memory &o_wrk = ins->o_UH;

  tombo::pressureRhs  (ins, time+dt, ins->o_rhsP);
  tombo::pressureSolve(ins, time+dt, ins->o_rhsP, o_wrk); 
  ins->o_P.copyFrom(o_wrk, ins->Ntotal*sizeof(dfloat)); 

  tombo::velocityRhs  (ins, time+dt, ins->o_rhsU);
  tombo::velocitySolve(ins, time+dt, ins->o_rhsU, o_wrk);

  for (int s=ins->Nstages;s>1;s--)
    o_U.copyFrom(o_U, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));

  o_U.copyFrom(o_wrk, ins->NVfields*ins->Ntotal*sizeof(dfloat)); 

  for (int s=ins->Nstages;s>1;s--) {
    ins->o_NU.copyFrom(ins->o_NU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
    
    ins->o_FU.copyFrom(ins->o_FU, ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-1)*ins->Ntotal*ins->NVfields*sizeof(dfloat), 
     (s-2)*ins->Ntotal*ins->NVfields*sizeof(dfloat));
  }
}

void velocityStrongSubCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_Ud)
{
  mesh_t *mesh = ins->mesh;

  const dlong NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);  

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
    dlong toffset = torder*ins->NVfields*ins->Ntotal;

    if (torder==ins->ExplicitOrder-1) { //first substep
      ins->scaledAddKernel(ins->NVfields*ins->Ntotal, b, toffset, o_U, zero, izero, o_Ud);
    } else { //add the next field
      ins->scaledAddKernel(ins->NVfields*ins->Ntotal, b, toffset, o_U,  one, izero, o_Ud);
    }     

    // SubProblem  starts from here from t^(n-torder)
    const dfloat tsub = time - torder*ins->dt;
    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<ins->Nsubsteps;++ststep){
      const dfloat tstage = tsub + ststep*ins->sdt;     

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
        ins->subCycleExtKernel(NtotalElements,
                               Nstages,
                               ins->fieldOffset,
                               ins->o_extC,
                               o_U,
                               ins->o_Ue);

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
            ins->o_rhsUd);
	}else{
	  ins->subCycleStrongVolumeKernel(mesh->Nelements,
	    mesh->o_vgeo,
	    mesh->o_Dmatrices,
	    ins->fieldOffset,
	    ins->o_Ue,
	    o_Ud,
	    ins->o_rhsUd);
	}
	
        ogsGatherScatterMany(ins->o_rhsUd, ins->NVfields, ins->fieldOffset,
                             ogsDfloat, ogsAdd, mesh->ogs);
  

        int nfield = ins->dim==2 ? 2:3; 
        ins->invMassMatrixKernel(mesh->Nelements,
				 ins->fieldOffset,
				 nfield,
				 mesh->o_vgeo,
				 ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
				 ins->o_rhsUd);

        ins->subCycleRKUpdateKernel(mesh->Nelements,
				    ins->sdt,
				    ins->Srka[rk],
				    ins->Srkb[rk],
				    ins->fieldOffset,
				    ins->o_rhsUd,
				    ins->o_resU, 
				    o_Ud);
      }
    }
  }
}

void scalarStrongSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_U, occa::memory o_S, occa::memory o_Sd){

  mesh_t *mesh = cds->mesh;

  const dfloat tn0 = time - 0*cds->dt;
  const dfloat tn1 = time - 1*cds->dt;
  const dfloat tn2 = time - 2*cds->dt;

  const dlong Nelements = cds->meshV->Nelements;

  dfloat zero = 0.0, one = 1.0;
  int izero = 0;
  dfloat b, bScale=0;

   // Solve for Each SubProblem
  for (int torder=(cds->ExplicitOrder-1); torder>=0; torder--){
    
    b=cds->extbdfB[torder];
    bScale += b; 

    // Initialize SubProblem Velocity i.e. Ud = U^(t-torder*dt)
    dlong toffset = torder*cds->NSfields*cds->Ntotal;

    if (torder==cds->ExplicitOrder-1) { //first substep
      cds->scaledAddKernel(cds->NSfields*cds->Ntotal, b, toffset, o_S, zero, izero, o_Sd);
    } else { //add the next field
      cds->scaledAddKernel(cds->NSfields*cds->Ntotal, b, toffset, o_S,  one, izero, o_Sd);
    }     

    // SubProblem  starts from here from t^(n-torder)
    const dfloat tsub = time - torder*cds->dt;
    // Advance SubProblem to t^(n-torder+1) 
    for(int ststep = 0; ststep<cds->Nsubsteps;++ststep){
     
      const dfloat tstage = tsub + ststep*cds->sdt;     
     
      for(int rk=0;rk<cds->SNrk;++rk){// LSERK4 stages
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
        cds->subCycleExtKernel(Nelements,
                               Nstages,
                               cds->vOffset,
                               cds->o_extC,
                               o_U,
                               cds->o_Ue);

     
        // Compute Volume Contribution
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          cds->subCycleStrongCubatureVolumeKernel(
                                            Nelements,
                                            mesh->o_vgeo,
                                            mesh->o_cubvgeo,
                                            mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
                                            mesh->o_cubInterpT,
                                            mesh->o_cubProjectT,
                                            cds->vOffset,
                                            cds->sOffset,             
                                            cds->o_Ue,
                                                 o_Sd,
                                            cds->o_rhsSd);
        } else{
          cds->subCycleStrongVolumeKernel(Nelements,
                                    mesh->o_vgeo,
                                    mesh->o_Dmatrices,
                                    cds->vOffset,
                                    cds->sOffset,           
                                    cds->o_Ue,
                                    o_Sd,
                                    cds->o_rhsSd);

        }

        ogsGatherScatter(cds->o_rhsSd, ogsDfloat, ogsAdd, mesh->ogs);

        // int nfield = ins->dim==2 ? 2:3; 
        cds->invMassMatrixKernel(Nelements,
                                 cds->sOffset,
                                 cds->NSfields,
                                 mesh->o_vgeo,
                                 cds->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
                                 cds->o_rhsSd);

        // Update Kernel
        cds->subCycleRKUpdateKernel(Nelements,
                                    cds->sdt,
                                    cds->Srka[rk],
                                    cds->Srkb[rk],
                                    cds->sOffset,
                                    cds->o_rhsSd,
                                    cds->o_resS, 
                                         o_Sd);
      }
    }
  }
}
