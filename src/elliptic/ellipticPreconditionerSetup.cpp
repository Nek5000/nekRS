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
#include "platform.hpp"
#include "gslib.h"
extern "C" void fem_amg_setup(const sint *n_x_, const sint *n_y_, const sint *n_z_, 
                   const sint *n_elem_, const sint *n_dim_, 
                   double *x_m_, double *y_m_, double *z_m_, 
                   double *pmask_, const sint *nullspace,
                   double *param,
                   MPI_Comm comm,
                   long long int* gatherGlobalNodes
                   );

void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs, occa::properties &kernelInfo)
{
  
  mesh_t* mesh = elliptic->mesh;
  precon_t* precon = elliptic->precon;
  setupAide options = elliptic->options;

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();

  if(options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    if(platform->comm.mpiRank == 0) printf("building MG preconditioner ... \n"); fflush(stdout);
    ellipticMultiGridSetup(elliptic,precon);
  } else if(options.compareArgs("PRECONDITIONER", "SEMFEM")) {
    // initialize solver here
    const int N = mesh->Nq;
    const int Nelements = mesh->Nelements;
    const int Ndim = 3;
    const int nullspace = (int) elliptic->allNeumann;
    double nullParams[] = {0}; // use default parameters
    double* mask = (double*) malloc(N*N*N*Nelements*sizeof(double));
    for(int i = 0; i < N*N*N*Nelements; ++i) mask[i] = 1.0;
    for(dlong n = 0; n < elliptic->Nmasked; n++){
      mask[elliptic->maskIds[n]] = 0.0;
    }
    
    fem_amg_setup(
      &N,
      &N,
      &N,
      &Nelements,
      &Ndim,
      mesh->x,
      mesh->y,
      mesh->z,
      mask,
      &nullspace,
      &nullParams[0],
      platform->comm.mpiComm,
      mesh->globalIds
    );

  } else if(options.compareArgs("PRECONDITIONER", "JACOBI")) {
    if(platform->comm.mpiRank == 0) printf("building Jacobi preconditioner ... "); fflush(stdout);
    precon->o_invDiagA = platform->device.malloc(elliptic->Nfields * elliptic->Ntotal ,  sizeof(dfloat));
    ellipticUpdateJacobi(elliptic);
  } else if(options.compareArgs("PRECONDITIONER", "NONE")) {
    // nothing 
  } else {
    printf("ERROR: Unknown preconditioner!\n");
    ABORT(EXIT_FAILURE);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}
