/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "parAlmond.hpp"

namespace parAlmond {


solver_t *Init(occa::device device, MPI_Comm comm, setupAide options) {
  solver_t *M = new solver_t(device, comm, options);

  if (Nrefs==0) buildParAlmondKernels(comm, device);
  Nrefs++;

  return M;
}

void AMGSetup(solver_t *MM,
               hlong* globalRowStarts,       //global partition
               dlong nnz,                    //--
               hlong* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
               hlong* Aj,                    //--
               dfloat* Avals,                //--
               bool nullSpace,
               dfloat nullSpacePenalty){

  solver_t *M = (solver_t *) MM;

  int rank, size;
  MPI_Comm_rank(M->comm, &rank);
  MPI_Comm_size(M->comm, &size);

  hlong TotalRows = globalRowStarts[M->size];
  dlong numLocalRows = (dlong) (globalRowStarts[M->rank+1]-globalRowStarts[M->rank]);

  MPI_Barrier(M->comm);
  double startTime = MPI_Wtime();
  if(rank==0) printf("Setting up AMG...");fflush(stdout);

  M->coarseLevel = new coarseSolver(M->options, M->comm);
  M->coarseLevel->setup(numLocalRows, globalRowStarts, nnz, Ai, Aj, Avals, nullSpace);
  M->baseLevel = M->numLevels;
  M->numLevels++;

  MPI_Barrier(M->comm);
  if(rank==0) printf("done (%gs)\n", MPI_Wtime()-startTime);
}

void Precon(solver_t *M, occa::memory o_x, occa::memory o_rhs) {

  M->levels[0]->o_x   = o_x;
  M->levels[0]->o_rhs = o_rhs;

  if       ((M->exact)&&(M->ktype==PCG)){
    M->device_pcg(1000,1e-8);
  } else if((M->exact)&&(M->ktype==GMRES)){
    M->device_pgmres(1000,1e-8);
  } else if(M->ctype==KCYCLE) {
    M->device_kcycle(0);
  } else if(M->ctype==VCYCLE) {
    //if(M->additive){
    if(0){
      M->additiveVcycle();
    } else {
      M->device_vcycle(0);
    }
  }
}

void Report(solver_t *M) {
  M->Report();
}

void Free(solver_t* M) {
  Nrefs--;
  if (Nrefs==0) {
    freeParAlmondKernels();
    freeScratchSpace();
    freePinnedScratchSpace();
  }

  delete M;
}

} //namespace parAlmond
