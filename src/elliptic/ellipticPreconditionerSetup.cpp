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
#include "ellipticPrecon.h"
#include "platform.hpp"

void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs)
{
  elliptic->precon = new precon_t();
  precon_t* precon = (precon_t*) elliptic->precon;
  const long long int pre = platform->device.occaDevice().memoryAllocated();
  
  mesh_t* mesh = elliptic->mesh;
  setupAide& options = elliptic->options;

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();

  if(options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    ellipticMultiGridSetup(elliptic,precon);
  } else if(options.compareArgs("PRECONDITIONER", "SEMFEM")) {
    if(platform->comm.mpiRank == 0) printf("building SEMFEM preconditioner ...\n"); fflush(stdout);
    precon->SEMFEMSolver = new SEMFEMSolver_t(elliptic);
  } else if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    if(platform->comm.mpiRank == 0) printf("building Jacobi preconditioner ... "); fflush(stdout);
    precon->o_invDiagA = platform->device.malloc(elliptic->Nfields * elliptic->fieldOffset, sizeof(pfloat));
    ellipticUpdateJacobi(elliptic, precon->o_invDiagA);
  } else if(options.compareArgs("PRECONDITIONER", "NONE")) {
    // nothing 
  } else {
    printf("ERROR: Unknown preconditioner!\n");
    ABORT(EXIT_FAILURE);
  }

  precon->preconBytes = platform->device.occaDevice().memoryAllocated() - pre;

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}

precon_t::~precon_t()
{
  if(SEMFEMSolver) delete SEMFEMSolver;
  if(MGSolver) delete MGSolver;
  o_diagA.free();
  o_invDiagA.free();
}
