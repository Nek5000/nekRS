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

#include "elliptic.h"
void ellipticPreconditioner(elliptic_t *elliptic, occa::memory &o_r, occa::memory &o_z){

  mesh_t *mesh = elliptic->mesh;
  precon_t *precon = elliptic->precon;
  setupAide options = elliptic->options;

  const dlong Nlocal = mesh->Np*mesh->Nelements;

  if(options.compareArgs("PRECONDITIONER", "JACOBI")){
    
    if(elliptic->blockSolver){
      elliptic->dotMultiplyKernel(Nlocal, elliptic->Ntotal, o_r, precon->o_invDiagA, o_z);      
    }else{
      elliptic->dotMultiplyKernel(Nlocal, o_r, precon->o_invDiagA, o_z);      
    }
  }
  else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {

    //mesh->device.finish();
    //double t0 = MPI_Wtime();
    parAlmond::Precon(precon->parAlmond, o_z, o_r);
    //ogsGatherScatter(o_z, ogsDfloat, ogsAdd, elliptic->ogs);
    //elliptic->collocateKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_z);
    //mesh->device.finish();
    //printf("tPreconditioner: %g\n", MPI_Wtime()-t0);

  } else if (options.compareArgs("PRECONDITIONER", "FULLALMOND")) {

    if (options.compareArgs("DISCRETIZATION", "IPDG")) {
      parAlmond::Precon(precon->parAlmond, o_z, o_r);
    } else if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
      ogsGather(precon->o_rhsG, o_r, ogsDfloat, ogsAdd, elliptic->ogs);
      elliptic->dotMultiplyKernel(elliptic->ogs->Ngather,
                      elliptic->ogs->o_gatherInvDegree, precon->o_rhsG, precon->o_rhsG);

      parAlmond::Precon(precon->parAlmond, precon->o_xG, precon->o_rhsG);
      ogsScatter(o_z, precon->o_xG, ogsDfloat, ogsAdd, elliptic->ogs);
    }

  } else if(options.compareArgs("PRECONDITIONER", "MASSMATRIX")){

    dfloat invLambda = 1./elliptic->lambda[0];

    if (options.compareArgs("DISCRETIZATION", "IPDG")) {
      precon->blockJacobiKernel(mesh->Nelements, invLambda, mesh->o_vgeo, precon->o_invMM, o_r, o_z);
    } else if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
      ogs_t *ogs = elliptic->ogs;

      elliptic->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_r, elliptic->o_rtmp);

      if(mesh->NglobalGatherElements)
        precon->partialblockJacobiKernel(mesh->NglobalGatherElements,
                                mesh->o_globalGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, elliptic->o_rtmp, o_z);

      ogsGatherScatterStart(o_z, ogsDfloat, ogsAdd, ogs);

      if(mesh->NlocalGatherElements)
        precon->partialblockJacobiKernel(mesh->NlocalGatherElements,
                                mesh->o_localGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, elliptic->o_rtmp, o_z);

      ogsGatherScatterFinish(o_z, ogsDfloat, ogsAdd, ogs);

      elliptic->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_z, o_z);

      //post-mask
      if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_z);
    }

  } else if (options.compareArgs("PRECONDITIONER", "SEMFEM")) {

    if (elliptic->elementType==TRIANGLES||elliptic->elementType==TETRAHEDRA) {
      o_z.copyFrom(o_r);
      elliptic->dotMultiplyKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_z, o_z);
      precon->SEMFEMInterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,o_z,precon->o_rFEM);
      ogsGather(precon->o_GrFEM, precon->o_rFEM, ogsDfloat, ogsAdd, precon->FEMogs);

      parAlmond::Precon(precon->parAlmond, precon->o_GzFEM, precon->o_GrFEM);

      ogsScatter(precon->o_zFEM, precon->o_GzFEM, ogsDfloat, ogsAdd, precon->FEMogs);
      precon->SEMFEMAnterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,precon->o_zFEM,o_z);
      elliptic->dotMultiplyKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_z, o_z);

      ogsGatherScatter(o_z, ogsDfloat, ogsAdd, elliptic->ogs);
      if (elliptic->Nmasked) mesh->maskKernel(elliptic->Nmasked, elliptic->o_maskIds, o_z);
    } else {
      ogsGather(precon->o_rhsG, o_r, ogsDfloat, ogsAdd, precon->FEMogs);
      elliptic->dotMultiplyKernel(precon->FEMogs->Ngather,
                      precon->FEMogs->o_gatherInvDegree, precon->o_rhsG, precon->o_rhsG);

      parAlmond::Precon(precon->parAlmond, precon->o_xG, precon->o_rhsG);

      ogsScatter(o_z, precon->o_xG, ogsDfloat, ogsAdd, precon->FEMogs);
    }

  }else if (options.compareArgs("PRECONDITIONER", "OAS")) {

    //    printf("IN OAS PRECONDITIONER\n");

    //ellipticOasSolve(elliptic, lambda, o_r, o_z);

  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }

  if(elliptic->allNeumann) // zero mean of RHS
    ellipticZeroMean(elliptic, o_z);
}


