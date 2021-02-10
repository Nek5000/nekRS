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

void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs, occa::properties &kernelInfo)
{
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = elliptic->precon;
  setupAide options = elliptic->options;

  if(options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    ellipticMultiGridSetup(elliptic,precon);
  } else if(options.compareArgs("PRECONDITIONER", "SEMFEM")) {
    //ellipticSEMFEMSetup(elliptic,precon);
    printf("ERROR: SEMFEM does not work right now.\n");
    ABORT(EXIT_FAILURE);;
  } else if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    dfloat* invDiagA;
    ellipticBuildJacobi(elliptic,&invDiagA);
    const dlong Nlocal =  mesh->Np * mesh->Nelements;
    int Ntotal = elliptic->blockSolver ? elliptic->Ntotal * elliptic->Nfields: Nlocal;
    precon->o_invDiagA = mesh->device.malloc(Ntotal * sizeof(dfloat), invDiagA);
    free(invDiagA);
  } else {
    printf("ERROR: Unknown preconditioner!\n");
    ABORT(EXIT_FAILURE);
  }
}
