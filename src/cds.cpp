/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "cds.h"


void cdsHelmholtzSolve(cds_t *cds, dfloat time, int stage,occa::memory o_rhsS,occa::memory o_Shat){
  
  mesh_t     *mesh   = cds->mesh; 
  elliptic_t *solver = cds->solver;
  
  cds->helmholtzRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            cds->o_EToB,
                            mesh->o_sMT,
                            cds->lambda,
                            time,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            cds->o_mapB,
                            o_rhsS);

  ogsGatherScatter(o_rhsS, ogsDfloat, ogsAdd, mesh->ogs);
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_rhsS);

  //copy current solution fields as initial guess? (could use Shat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  o_Shat.copyFrom(cds->o_S,Ntotal*sizeof(dfloat),0,0*cds->sOffset*sizeof(dfloat)); 
 
  if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, o_Shat);
  
  cds->Niter = ellipticSolve(solver, cds->lambda, cds->TOL, o_rhsS, o_Shat);

  cds->helmholtzAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            cds->o_mapB,
                            o_Shat);
}

void cdsHelmholtzRhs(cds_t *cds, dfloat time, int stage, occa::memory o_rhsS){
  
  mesh_t *mesh = cds->mesh; 

  // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
  cds->helmholtzRhsKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_MM,
                          cds->idt,
                          cds->idiff,
                          cds->o_extbdfA,
                          cds->o_extbdfB,
                          cds->o_extbdfC,
                          cds->sOffset,
                          cds->o_S,
                          cds->o_NS,
                          cds->o_FS,
                          o_rhsS);
}

void cdsAdvection(cds_t *cds, dfloat time, occa::memory o_U, occa::memory o_S, occa::memory o_NS){

  mesh_t *mesh = cds->mesh;

  if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE"))
    cds->advectionStrongCubatureVolumeKernel(
           mesh->Nelements,
           mesh->o_vgeo,
           mesh->o_cubvgeo,
           mesh->o_cubDiffInterpT, // mesh->o_cubDWmatrices,
           mesh->o_cubInterpT,
           mesh->o_cubProjectT,
           cds->vOffset,
           cds->sOffset,
           o_U,
           o_S,
           o_NS);
  else
    cds->advectionStrongVolumeKernel(
           mesh->Nelements,
           mesh->o_vgeo,
           mesh->o_Dmatrices,
           cds->vOffset,
           cds->sOffset,
           o_U,
           o_S,
           o_NS);

}

void cdsStrongSubCycle(cds_t *cds, dfloat time, int Nstages, occa::memory o_U, occa::memory o_S, occa::memory o_Sd){

  mesh_t *mesh = cds->mesh;
  const dlong NtotalElements = (mesh->Nelements+mesh->totalHaloPairs);  

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
        cds->subCycleExtKernel(NtotalElements,
                               Nstages,
                               cds->vOffset,
                               cds->o_extC,
                               o_U,
                               cds->o_Ue);

     
        // Compute Volume Contribution
        if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
          cds->subCycleStrongCubatureVolumeKernel(mesh->Nelements,
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
          cds->subCycleStrongVolumeKernel(mesh->Nelements,
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
        cds->invMassMatrixKernel(mesh->Nelements,
                                 cds->sOffset,
                                 cds->NSfields,
                                 mesh->o_vgeo,
                                 cds->o_InvM, // mesh->o_MM, // should be invMM for tri/tet
                                 cds->o_rhsSd);

        // Update Kernel
        cds->subCycleRKUpdateKernel(mesh->Nelements,
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
