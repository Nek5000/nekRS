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
#include "timer.hpp"

void ellipticPreconditioner(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_z)
{
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = elliptic->precon;
  setupAide options = elliptic->options;

  const dlong Nlocal = mesh->Np * mesh->Nelements;

  if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    if(elliptic->blockSolver)
      elliptic->dotMultiplyKernel(Nlocal, elliptic->Ntotal, o_r, precon->o_invDiagA, o_z);
    else
      elliptic->dotMultiplyKernel(Nlocal, o_r, precon->o_invDiagA, o_z);
  }else if (options.compareArgs("PRECONDITIONER", "MULTIGRID"))  {
    timer::tic("preconditioner", 1);
    parAlmond::Precon(precon->parAlmond, o_z, o_r);
    //ogsGatherScatter(o_z, ogsDfloat, ogsAdd, elliptic->ogs);
    //elliptic->collocateKernel(mesh->Nelements*mesh->Np, elliptic->o_invDegree, o_z);
    timer::toc("preconditioner");
  }else  {
    if(mesh->rank == 0) printf("ERRROR: Unknown preconditioner\n");
    MPI_Abort(mesh->comm, 1);
    //o_z.copyFrom(o_r);
  }

  if(elliptic->allNeumann) // zero mean of RHS
    ellipticZeroMean(elliptic, o_z);
}
