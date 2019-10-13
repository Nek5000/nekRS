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


