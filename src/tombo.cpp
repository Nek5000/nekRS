#include "nekrs.hpp"

// Compute divergence of the velocity field using physical boundary data at t = time. 
void insDivergence(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_DU){

  mesh_t *mesh = ins->mesh;

  //  if (ins->vOptions.compareArgs("DISCRETIZATION","IPDG")) {
  if(mesh->totalHaloPairs>0){

    // make sure compute device is ready to perform halo extract
    mesh->device.finish();
    
    // switch to data stream
    mesh->device.setStream(mesh->dataStream);
    
    ins->haloGetKernel(mesh->totalHaloPairs,
		       ins->NVfields,
		       ins->fieldOffset,
		       mesh->o_haloElementList,
		       mesh->o_haloGetNodeIds,
		       o_U,
		       ins->o_vHaloBuffer);
    
    dlong Ndata = ins->NVfields*mesh->Nfp*mesh->totalHaloPairs;
    
    // copy extracted halo to HOST 
    ins->o_vHaloBuffer.copyTo(ins->vSendBuffer, Ndata*sizeof(dfloat), 0, "async: true");// zero offset             

    mesh->device.setStream(mesh->defaultStream);
  }

  // computes div u^(n+1) volume term
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             ins->fieldOffset,
                             o_U,
                             o_DU);
  
  if(mesh->totalHaloPairs>0){

    // make sure compute device is ready to perform halo extract
    mesh->device.setStream(mesh->dataStream);
    mesh->device.finish();
    
    // start halo exchange
    meshHaloExchangeStart(mesh,
			  mesh->Nfp*(ins->NVfields)*sizeof(dfloat),
			  ins->vSendBuffer,
			  ins->vRecvBuffer);
        
    meshHaloExchangeFinish(mesh);
    
    dlong Ndata = ins->NVfields*mesh->Nfp*mesh->totalHaloPairs;
    
    ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer, Ndata*sizeof(dfloat), 0);  // zero offset
    
    ins->haloPutKernel(mesh->totalHaloPairs,
		       ins->NVfields,
		       ins->fieldOffset,
		       mesh->o_haloElementList,
		       mesh->o_haloPutNodeIds,
		       ins->o_vHaloBuffer,
		       o_U);

    mesh->device.setStream(mesh->defaultStream);
  }
  
  //computes div u^(n+1) surface term
  const dfloat lambda = -ins->g0*ins->idt; 
  ins->divergenceSurfaceKernel(
    mesh->Nelements,
    mesh->o_vgeo,
    mesh->o_sgeo,
    mesh->o_LIFTT,
    mesh->o_vmapM,
    mesh->o_vmapP,
    ins->o_EToB,
    time,
    lambda, 
    mesh->o_x,
    mesh->o_y,
    mesh->o_z,
    ins->fieldOffset,
    ins->o_Wrk,
    o_U,
    o_DU);
}

void pressureRhs(ins_t *ins, dfloat time)
{
  mesh_t *mesh = ins->mesh;

  ins->pressureRhsKernel(mesh->Nelements,
                         mesh->o_vgeo,
                         mesh->o_MM,
                         ins->idt,
                         ins->nu,
                         ins->g0,
                         ins->o_extbdfA,
                         ins->o_extbdfB,
                         ins->fieldOffset,
                         ins->o_U,
                         ins->o_NU,
                         ins->o_NC,
                         ins->o_FU,
                         ins->o_rkU);

  // weak divergence term (not that no boundary contribution) 
  // -> (F, grad(phi))_sigma - g_0/dt*(n*U^n+1)_del_sigma
  insDivergence(ins, time, ins->o_rkU, ins->o_rhsP);

  // Add  -(grad P, grad phi) to rhsP
  const dfloat lambda = 0.0;
  ins->pressureAxKernel(mesh->Nelements,
                        mesh->o_ggeo,
                        mesh->o_Dmatrices,
                        mesh->o_Smatrices,
                        mesh->o_MM,
                        lambda,
                        ins->o_P,
                        ins->o_rhsP);
}


void pressureSolve(ins_t *ins, dfloat time, occa::memory o_rkP)
{
  mesh_t *mesh = ins->mesh;
  elliptic_t *solver = ins->pSolver;

  ogsGatherScatter(ins->o_rhsP, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsP);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);

  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI);

  ins->pressureAddBCKernel(mesh->Nelements,
                           time,
                           ins->dt,
                           ins->fieldOffset,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           ins->o_PmapB,
                           ins->o_EToB,
                           ins->o_Wrk,
                           ins->o_U,
                           ins->o_PI);

  ins->pressureUpdateKernel(mesh->Nelements,
                            ins->fieldOffset,
                            ins->o_PmapB,
                            ins->o_PI,
                            ins->o_P,
                            o_rkP);
}

void insGradient(ins_t *ins, dfloat time, occa::memory o_P, occa::memory o_GP)
{
  mesh_t *mesh = ins->mesh;

  if(mesh->totalHaloPairs>0){

    // make sure compute device is ready to perform halo extract
    mesh->device.finish();
    
    // switch to data stream
    mesh->device.setStream(mesh->dataStream);

    int one = 1;
    dlong Ndata = one*mesh->Nfp*mesh->totalHaloPairs;

    ins->haloGetKernel(mesh->totalHaloPairs,
		       one,
		       ins->fieldOffset,
		       mesh->o_haloElementList,
		       mesh->o_haloGetNodeIds,
		       o_P,
		       ins->o_pHaloBuffer);   
    
    // copy extracted halo to HOST 
    ins->o_pHaloBuffer.copyTo(ins->pSendBuffer,
			      Ndata*sizeof(dfloat), 0, "async: true");
    
    mesh->device.setStream(mesh->defaultStream);
  }
  
  // Compute Volume Contribution
  ins->gradientVolumeKernel(mesh->Nelements,
			    mesh->o_vgeo,
			    mesh->o_Dmatrices,
                            ins->fieldOffset,
                            o_P,
                            o_GP);

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

    // make sure compute device is ready to perform halo extract
    mesh->device.setStream(mesh->dataStream);
    mesh->device.finish();

    // start halo exchange
    meshHaloExchangeStart(mesh,
			  mesh->Nfp*sizeof(dfloat),
			  ins->pSendBuffer,
			  ins->pRecvBuffer);
    
    meshHaloExchangeFinish(mesh);
    
    int one = 1;
    dlong Ndata = one*mesh->Nfp*mesh->totalHaloPairs;
    
    ins->o_pHaloBuffer.copyFrom(ins->pRecvBuffer, Ndata*sizeof(dfloat), 0);  // zero offset
    
    ins->haloPutKernel(mesh->totalHaloPairs,
		       one,
		       ins->fieldOffset,
		       mesh->o_haloElementList,
		       mesh->o_haloPutNodeIds,
		       ins->o_pHaloBuffer,
		       o_P);

    mesh->device.finish();
    mesh->device.setStream(mesh->defaultStream);
    mesh->device.finish();
  }
}

void velocityRhs(ins_t *ins, dfloat time) 
{
  mesh_t *mesh = ins->mesh;
  insGradient (ins, time, ins->o_P, ins->o_rkGP);

  ins->velocityRhsKernel(mesh->Nelements,
                         mesh->o_vgeo,
                         mesh->o_MM,
                         ins->idt,
                         ins->inu,
                         ins->g0,
                         ins->o_extbdfA,
                         ins->o_extbdfB,
                         ins->fieldOffset,
                         ins->o_U,
                         ins->o_NU,
                         ins->o_FU,
                         ins->o_rkGP,
                         ins->o_rhsU,
                         ins->o_rhsV,
                         ins->o_rhsW);
}

void velocitySolve(ins_t *ins, dfloat time, occa::memory o_Uhat) 
{
  mesh_t *mesh = ins->mesh;
  elliptic_t *usolver = ins->uSolver;
  elliptic_t *vsolver = ins->vSolver;
  elliptic_t *wsolver = ins->wSolver;

  ins->velocityRhsBCKernel(mesh->Nelements,
                           ins->fieldOffset,
                           mesh->o_ggeo,
                           mesh->o_sgeo,
                           mesh->o_Dmatrices,
                           mesh->o_Smatrices,
                           mesh->o_MM,
                           mesh->o_vmapM,
                           ins->o_EToB,
                           mesh->o_sMT,
                           ins->lambda,
                           time,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           ins->o_VmapB,
                           ins->o_Wrk,
                           ins->o_U,
                           ins->o_rhsU,
                           ins->o_rhsV,
                           ins->o_rhsW);

  // TODO: change me to gs_many
  ogsGatherScatter(ins->o_rhsU, ogsDfloat, ogsAdd, mesh->ogs);
  ogsGatherScatter(ins->o_rhsV, ogsDfloat, ogsAdd, mesh->ogs);
  if (ins->dim==3)
    ogsGatherScatter(ins->o_rhsW, ogsDfloat, ogsAdd, mesh->ogs);

  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  // Use old velocity for velocity solver initial condition
  // TODO: fure into one kernel
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,0*ins->fieldOffset*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,1*ins->fieldOffset*sizeof(dfloat));
  if (ins->dim==3)
    ins->o_WH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,2*ins->fieldOffset*sizeof(dfloat));

  // TODO: fuse into single kernel 
  if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_UH);
  if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_VH);
  if (ins->dim==3 && wsolver->Nmasked)
    mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_WH);

  // TODO: fuse into single kernel
  if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_rhsU);
  if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_rhsV);
  if (ins->dim==3 && wsolver->Nmasked)
    mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_rhsW);

  ins->NiterU = ellipticSolve(usolver, ins->lambda, ins->velTOL, ins->o_rhsU, ins->o_UH);
  ins->NiterV = ellipticSolve(vsolver, ins->lambda, ins->velTOL, ins->o_rhsV, ins->o_VH);
  if (ins->dim==3) 
    ins->NiterW = ellipticSolve(wsolver, ins->lambda, ins->velTOL, ins->o_rhsW, ins->o_WH);
  
  ins->velocityAddBCKernel(mesh->Nelements,
                           ins->fieldOffset,
                           time,
                           mesh->o_sgeo,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mesh->o_vmapM,
                           ins->o_VmapB,
                           ins->o_Wrk,
                           ins->o_U,
                           ins->o_UH,
                           ins->o_VH,
                           ins->o_WH);

  //copy into intermediate stage storage
  // TODO: fuse into one copyTo
  ins->o_UH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
  ins->o_VH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);
  if (ins->dim==3)
    ins->o_WH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);
}

void subCycle(ins_t *ins, dfloat time, int Nstages, occa::memory o_U, occa::memory o_Ud)
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
            ins->o_invLumpedMassMatrix,
            ins->fieldOffset,
            ins->o_Ue,
            o_Ud,
            ins->o_cU,     
            ins->o_cUd,     
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
	
        // TODO: change me to gs_many
	for(int k=0;k<ins->dim;++k){
	  ogsGatherScatter(ins->o_rhsUd+k*ins->fieldOffset*sizeof(dfloat), 
                           ogsDfloat, ogsAdd, mesh->ogs);
	} 

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

void curlCurl(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NC)
{
  mesh_t *mesh = ins->mesh;

  // TODO: extrapolte velocity

  // No that UH, VH and WH is used as the temporary arrays to store curl(u)
  // Note that multiplied with Mass Matrix i.e. JW to prepare smoothing !!!!!
  ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  o_U,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz

  // TODO: change me to gs_many
  if(ins->dim==3){    
   ogsGatherScatter(ins->o_UH, ogsDfloat, ogsAdd, mesh->ogs);
   ogsGatherScatter(ins->o_VH, ogsDfloat, ogsAdd, mesh->ogs);
  } 
   ogsGatherScatter(ins->o_WH, ogsDfloat, ogsAdd, mesh->ogs);

  int nfield = 0; 
  if (ins->dim==3){
    nfield = 3; 
    // TODO: fuse into single copyTo
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);    
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);    
  }else{
    nfield =1; 
    // if quad or tri copy Wz to first place
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0); 
  }
  
 
  // Multiply With Inverse Mass Matrix
  ins->invMassMatrixKernel(
    mesh->Nelements,
    ins->fieldOffset,
    nfield,
    mesh->o_vgeo,
    ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
    ins->o_rkU);


  if(ins->dim==2){
  // Second curl on smoothed curl(u)
  ins->curlBKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  ins->o_rkU,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz
  }else{
    ins->curlKernel(mesh->Nelements,
                  mesh->o_vgeo,
                  mesh->o_Dmatrices,
                  mesh->o_MM, 
                  ins->fieldOffset,
                  ins->o_rkU,
                  ins->o_UH, // Wx
                  ins->o_VH, // Wy
                  ins->o_WH);// Wz
  }
 
  // TODO: change me to gs_many
  ogsGatherScatter(ins->o_UH, ogsDfloat, ogsAdd, mesh->ogs);
  ogsGatherScatter(ins->o_VH, ogsDfloat, ogsAdd, mesh->ogs);
  if(ins->dim==3)       
    ogsGatherScatter(ins->o_WH, ogsDfloat, ogsAdd, mesh->ogs);

  nfield = 0; 
  if (ins->dim==3){
    nfield = 3;
    // TODO: fuse into single copyTo 
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);    
    ins->o_WH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);    
  }else{
    nfield =2; 
    // TODO: fuse into single copyTo 
    ins->o_UH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_rkU,ins->Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0); 
  }

  ins->invMassMatrixKernel(
    mesh->Nelements,
    ins->fieldOffset,
    nfield,
    mesh->o_vgeo,
    ins->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
    ins->o_rkU);
   
  ins->o_rkU.copyTo(o_NC,ins->NVfields*ins->Ntotal*sizeof(dfloat),
                    0*ins->fieldOffset*sizeof(dfloat),0);
}



