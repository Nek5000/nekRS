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

#include "stdio.h"
#include "parAlmond.hpp"

#include "timer.hpp"

#include "omp.h"
#include "limits.h"
#include "boomerAMG.h"
#include "amgx.h"
#include "platform.hpp"

namespace parAlmond {

coarseSolver::coarseSolver(setupAide options_, MPI_Comm comm_) {
  gatherLevel = false;
  options = options_;
  comm = comm_;
}

int coarseSolver::getTargetSize() {
  if(options.compareArgs("AMG SOLVER", "BOOMERAMG"))
    return INT_MAX;
  else 
    return 1000;
}

void coarseSolver::setup(
               dlong Nrows, 
      	       hlong* globalRowStarts,       //global partition
               dlong nnz,                    //--
               hlong* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
               hlong* Aj,                    //--
               dfloat* Avals,                //--
               bool nullSpace)
{
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

   if(options.compareArgs("BUILD ONLY", "TRUE"))
    return; // bail early as this will not get used

  if(options.compareArgs("PARALMOND SMOOTH COARSEST", "TRUE"))
    return; // bail early as this will not get used

  if ((rank==0)&&(options.compareArgs("VERBOSE","TRUE")))
    printf("Setting up coarse solver...");fflush(stdout);

  if (options.compareArgs("AMG SOLVER", "BOOMERAMG")){
    int Nthreads = 1;
 
    double settings[BOOMERAMG_NPARAM+1];
    settings[0]  = 1;    /* custom settings              */
    settings[1]  = 8;    /* coarsening                   */
    settings[2]  = 6;    /* interpolation                */
    settings[3]  = 1;    /* number of cycles             */
    settings[4]  = 16;   /* smoother for crs level       */
    settings[5]  = 3;    /* number of coarse sweeps      */
    settings[6]  = 16;   /* smoother                     */
    settings[7]  = 1;    /* number of sweeps             */
    settings[8]  = 0.25; /* strong threshold             */
    settings[9]  = 0.05; /* non galerkin tol             */
    settings[10] = 0;    /* aggressive coarsening levels */

    options.getArgs("BOOMERAMG COARSEN TYPE", settings[1]);
    options.getArgs("BOOMERAMG INTERPOLATION TYPE", settings[2]);
    options.getArgs("BOOMERAMG SMOOTHER TYPE", settings[6]);
    options.getArgs("BOOMERAMG SMOOTHER SWEEPS", settings[7]);
    options.getArgs("BOOMERAMG ITERATIONS", settings[3]);
    options.getArgs("BOOMERAMG STRONG THRESHOLD", settings[8]);
    options.getArgs("BOOMERAMG NONGALERKIN TOLERANCE" , settings[9]);
    options.getArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS" , settings[10]);

    boomerAMGSetup(Nrows,
                   nnz,
                   Ai,
                   Aj,
                   Avals,
                   (int) nullSpace,
                   comm,
                   Nthreads,
                   -1, /* device ID, if negative run on host */
                   0,  /* useFP32 */
                   settings);
 
    N = (int) Nrows;
    xLocal   = (dfloat*) calloc(N,sizeof(dfloat));
    rhsLocal = (dfloat*) calloc(N,sizeof(dfloat));
  }
  else if (options.compareArgs("AMG SOLVER", "AMGX")){
    const int useFP32 = options.compareArgs("SEMFEM SOLVER PRECISION", "FP32");
    if(useFP32)
    {
      if(platform->comm.mpiRank == 0) printf("FP32 is not supported for the coarse grid solver.\n");
      MPI_Barrier(platform->comm.mpiComm);
      ABORT(1);
    }
    if(platform->device.mode() != "CUDA") {
      if(platform->comm.mpiRank == 0) printf("AmgX only supports CUDA!\n");
      MPI_Barrier(platform->comm.mpiComm);
      ABORT(1);
    } 
    string configFile;
    options.getArgs("AMGX CONFIG FILE", configFile);
    char *cfg = NULL;
    if(configFile.size()) cfg = (char*) configFile.c_str();
    AMGXsetup(
      Nrows,
      nnz,
      Ai,
      Aj,
      Avals,
      (int) nullSpace,
      comm,
      platform->device.id(),
      useFP32,
      std::stoi(getenv("NEKRS_GPU_MPI")),
      cfg);
    N = (int) Nrows;
  } else {
    if(platform->comm.mpiRank == 0){
      std::string amgSolver;
      options.getArgs("SEMFEM SOLVER", amgSolver);
      printf("SEMFEM SOLVER %s is not supported!\n", amgSolver.c_str());
    }
    ABORT(EXIT_FAILURE);
  }
}

void coarseSolver::setup(parCSR *A) {

  comm = A->comm;

  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  if ((rank==0)&&(options.compareArgs("VERBOSE","TRUE")))
    printf("Setting up coarse solver...");fflush(stdout);

  if (options.compareArgs("AMG SOLVER", "BOOMERAMG")){
    int Nthreads = 1;
 
    int totalNNZ = A->diag->nnz + A->offd->nnz;
    hlong *rows;
    hlong *cols;
    dfloat *vals;
 
    if (totalNNZ) {
      rows = (hlong *) calloc(totalNNZ,sizeof(hlong));
      cols = (hlong *) calloc(totalNNZ,sizeof(hlong));
      vals = (dfloat *) calloc(totalNNZ,sizeof(dfloat));
    }
 
    // populate local COO (AIJ) matrix using global (row,col) indices
    int cnt = 0;
    for (int n = 0; n < A->Nrows; n++) {
      for (int m = A->diag->rowStarts[n]; m < A->diag->rowStarts[n+1]; m++) {
        rows[cnt] = n + A->globalRowStarts[rank];
        cols[cnt] = A->diag->cols[m] + A->globalRowStarts[rank];
        vals[cnt] = A->diag->vals[m];
        cnt++;
      }
      for (int m = A->offd->rowStarts[n]; m < A->offd->rowStarts[n+1]; m++) {
        rows[cnt] = n + A->globalRowStarts[rank];
        cols[cnt] = A->colMap[A->offd->cols[m]];
        vals[cnt] = A->offd->vals[m];
        cnt++;
      }
    }

    setup(A->Nrows, A->globalRowStarts, totalNNZ, rows, cols, vals,(int) A->nullSpace); 
 
    if (totalNNZ) {
      free(rows);
      free(cols);
      free(vals);
    }
 
    return;
  }

  //copy the global coarse partition as ints
  coarseOffsets = (int* ) calloc(size+1,sizeof(int));
  for (int r=0;r<size+1;r++) coarseOffsets[r] = (int) A->globalRowStarts[r];

  coarseTotal   = coarseOffsets[size];
  coarseOffset  = coarseOffsets[rank];

  N = (int) A->Nrows;

  coarseCounts = (int*) calloc(size,sizeof(int));

  if(A->nullSpace){
    if(rank==0) printf("Current null space handling not available for parAlmond!\n");
    fflush(stdout);
    exit(1);
  }

  int sendNNZ = (int) (A->diag->nnz+A->offd->nnz);
  int *rows;
  int *cols;
  dfloat *vals;

  // Make the MPI_NONZERO_T data type
  nonzero_t NZ;
  MPI_Datatype MPI_NONZERO_T;
  MPI_Datatype dtype[3] = {MPI_HLONG, MPI_HLONG, MPI_DFLOAT};
  int blength[3] = {1, 1, 1};
  MPI_Aint addr[3], displ[3];
  MPI_Get_address ( &(NZ.row), addr+0);
  MPI_Get_address ( &(NZ.col), addr+1);
  MPI_Get_address ( &(NZ.val), addr+2);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  MPI_Type_create_struct (3, blength, displ, dtype, &MPI_NONZERO_T);
  MPI_Type_commit (&MPI_NONZERO_T);

  nonzero_t *sendNonZeros = (nonzero_t *) calloc(sendNNZ, sizeof(nonzero_t));

  //populate matrix
  int cnt = 0;
  for (int n=0;n<N;n++) {
    int start = (int) A->diag->rowStarts[n];
    int end   = (int) A->diag->rowStarts[n+1];
    for (int m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A->diag->cols[m] + coarseOffset;
      sendNonZeros[cnt].val = A->diag->vals[m];
      cnt++;
    }
    start = (int) A->offd->rowStarts[n];
    end   = (int) A->offd->rowStarts[n+1];
    for (dlong m=start;m<end;m++) {
      sendNonZeros[cnt].row = n + coarseOffset;
      sendNonZeros[cnt].col = A->colMap[A->offd->cols[m]];
      sendNonZeros[cnt].val = A->offd->vals[m];
      cnt++;
    }
  }

  //get the nonzero counts from all ranks
  int *recvNNZ    = (int*) calloc(size,sizeof(int));
  int *NNZoffsets = (int*) calloc(size+1,sizeof(int));
  MPI_Allgather(&sendNNZ, 1, MPI_INT,
                 recvNNZ, 1, MPI_INT, comm);

  int totalNNZ = 0;
  for (int r=0;r<size;r++) {
    totalNNZ += recvNNZ[r];
    NNZoffsets[r+1] = NNZoffsets[r] + recvNNZ[r];
  }

  nonzero_t *recvNonZeros = (nonzero_t *) calloc(totalNNZ, sizeof(nonzero_t));

  MPI_Allgatherv(sendNonZeros, sendNNZ,             MPI_NONZERO_T,
                 recvNonZeros, recvNNZ, NNZoffsets, MPI_NONZERO_T, comm);

  //gather null vector
  dfloat *nullTotal = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  for (int r=0;r<size;r++)
    coarseCounts[r] = coarseOffsets[r+1]-coarseOffsets[r];

  MPI_Allgatherv(  A->null,          N,                MPI_DFLOAT,
                 nullTotal, coarseCounts, coarseOffsets, MPI_DFLOAT,
                 comm);

  //clean up
  MPI_Barrier(comm);
  MPI_Type_free(&MPI_NONZERO_T);
  free(sendNonZeros);
  free(NNZoffsets);
  free(recvNNZ);

  //assemble the full matrix
  dfloat *coarseA = (dfloat *) calloc(coarseTotal*coarseTotal,sizeof(dfloat));
  for (int i=0;i<totalNNZ;i++) {
    int n = recvNonZeros[i].row;
    int m = recvNonZeros[i].col;
    coarseA[n*coarseTotal+m] = recvNonZeros[i].val;
  }

  if (A->nullSpace) { //A is dense due to nullspace augmentation
    for (int n=0;n<coarseTotal;n++) {
      for (int m=0;m<coarseTotal;m++) {
        coarseA[n*coarseTotal+m] += A->nullSpacePenalty*nullTotal[n]*nullTotal[m];
      }
    }
  }

  free(recvNonZeros);
  free(nullTotal);

  matrixInverse(coarseTotal, coarseA);

  //store only the local rows of the full inverse
  invCoarseA = (dfloat *) calloc(N*coarseTotal,sizeof(dfloat));
  for (int n=0;n<N;n++) {
    for (int m=0;m<coarseTotal;m++) {
      invCoarseA[n*coarseTotal+m] = coarseA[(n+coarseOffset)*coarseTotal+m];
    }
  }

  xLocal   = (dfloat*) calloc(N,sizeof(dfloat));
  rhsLocal = (dfloat*) calloc(N,sizeof(dfloat));

  xCoarse   = (dfloat*) calloc(coarseTotal,sizeof(dfloat));
  rhsCoarse = (dfloat*) calloc(coarseTotal,sizeof(dfloat));

  free(coarseA);

  if((rank==0)&&(options.compareArgs("VERBOSE","TRUE"))) printf("done.\n");
}

void coarseSolver::syncToDevice() {}

void coarseSolver::solve(dfloat *rhs, dfloat *x) {
  exit(1);
}

void coarseSolver::gather(occa::memory o_rhs, occa::memory o_x)
{
  if (gatherLevel) {
    //weight
    vectorDotStar(ogs->N, 1.0, ogs->o_invDegree, o_rhs, 0.0, o_Sx);
    ogsGather(o_Gx, o_Sx, ogsDfloat, ogsAdd, ogs);
    if(N)
      o_Gx.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  } else {
    if(N)
      o_rhs.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  }
}
void coarseSolver::scatter(occa::memory o_rhs, occa::memory o_x)
{
  if (gatherLevel) {
    if(N)
      o_Gx.copyFrom(xLocal, N*sizeof(dfloat), 0);
    ogsScatter(o_x, o_Gx, ogsDfloat, ogsAdd, ogs);
  } else {
    if(N)
      o_x.copyFrom(xLocal, N*sizeof(dfloat), 0);
  }
}
void coarseSolver::BoomerAMGSolve() {
  platform->timer.hostTic("BoomerAMGSolve", 1);
  boomerAMGSolve(xLocal, rhsLocal);
  platform->timer.hostToc("BoomerAMGSolve");
}
void coarseSolver::AmgXSolve(occa::memory o_rhs, occa::memory o_x) {
  platform->timer.tic("AmgXSolve", 1);
  AMGXsolve(o_x.ptr(), o_rhs.ptr());
  platform->timer.toc("AmgXSolve");
}
void coarseSolver::solve(occa::memory o_rhs, occa::memory o_x) {

  const bool useDevice = options.compareArgs("AMG SOLVER", "AMGX");
  if (gatherLevel) {
    //weight
    vectorDotStar(ogs->N, 1.0, ogs->o_invDegree, o_rhs, 0.0, o_Sx);
    ogsGather(o_Gx, o_Sx, ogsDfloat, ogsAdd, ogs);
    if(N && !useDevice)
      o_Gx.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  } else {
    if(N && !useDevice)
      o_rhs.copyTo(rhsLocal, N*sizeof(dfloat), 0);
  }

  if (options.compareArgs("AMG SOLVER", "BOOMERAMG")){
    BoomerAMGSolve(); 
  } else if (options.compareArgs("AMG SOLVER", "AMGX")){
    occa::memory o_b = gatherLevel ? o_Gx : o_rhs;
    AmgXSolve(o_b, o_x);
  } else {
    //gather the full vector
    MPI_Allgatherv(rhsLocal,             N,                MPI_DFLOAT,
                   rhsCoarse, coarseCounts, coarseOffsets, MPI_DFLOAT, comm);

    //multiply by local part of the exact matrix inverse
    // #pragma omp parallel for
    for (int n=0;n<N;n++) {
      xLocal[n] = 0.;
      for (int m=0;m<coarseTotal;m++) {
        xLocal[n] += invCoarseA[n*coarseTotal+m]*rhsCoarse[m];
      }
    }
  }

  if (gatherLevel) {
    if(N && !useDevice)
      o_Gx.copyFrom(xLocal, N*sizeof(dfloat));
    if(N && useDevice)
      o_Gx.copyFrom(o_x, N*sizeof(dfloat));
    ogsScatter(o_x, o_Gx, ogsDfloat, ogsAdd, ogs);
  } else {
    if(N && !useDevice)
      o_x.copyFrom(xLocal, N*sizeof(dfloat), 0);
  }
}

} //namespace parAlmond
