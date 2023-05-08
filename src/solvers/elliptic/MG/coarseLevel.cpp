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

#include "limits.h"
#include "stdio.h"
#include "timer.hpp"

#include "AMGX.hpp"

#include "platform.hpp"
#include "linAlg.hpp"
#include "MGSolver.hpp"

static occa::kernel vectorDotStarKernel;

MGSolver_t::coarseLevel_t::coarseLevel_t(setupAide options, MPI_Comm comm)
{
  this->options = options;
  this->comm = comm;
  this->solvePtr = &MGSolver_t::coarseLevel_t::solve;
}

void MGSolver_t::coarseLevel_t::setupSolver(
               hlong* globalRowStarts,
               dlong nnz,                    //--
               hlong* Ai,                    //-- Local A matrix data (globally indexed, COO storage, row sorted)
               hlong* Aj,                    //--
               dfloat* Avals,                //--
               bool nullSpace)
{
  int rank, size;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  MPI_Barrier(comm);
  double startTime = MPI_Wtime();
  if(rank==0) printf("setup FEM solver ...");fflush(stdout);

  N = (dlong) (globalRowStarts[rank+1] - globalRowStarts[rank]);

  const int verbose = (platform->options.compareArgs("VERBOSE", "TRUE")) ? 1 : 0;
  const bool useDevice = options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");
  const int useFP32 = options.compareArgs("COARSE SOLVER PRECISION", "FP32");

  const std::string kernelName = "vectorDotStar";
  if (!vectorDotStarKernel.isInitialized()) 
    vectorDotStarKernel = platform->kernels.get(kernelName);

  o_xBuffer = platform->device.malloc(N * sizeof(pfloat));
  h_xBuffer = platform->device.mallocHost(N * sizeof(pfloat));
  xBuffer = (pfloat*) h_xBuffer.ptr(); 

  if (options.compareArgs("COARSE SOLVER", "BOOMERAMG")){
 
    double settings[hypreWrapperDevice::NPARAM+1];
    settings[0]  = 1;    /* custom settings              */
    settings[1]  = 10;   /* coarsening                   */
    if (useDevice) 
      settings[1]  = 8;  /*  PMIS currently not supported on device */
    settings[2]  = 6;    /* interpolation                */
    settings[3]  = 1;    /* number of cycles             */
    settings[4]  = 16;   /* smoother for crs level       */
    settings[5]  = 3;    /* number of coarse sweeps      */
    settings[6]  = 16;   /* smoother                     */
    settings[7]  = 1;    /* number of sweeps             */
    settings[8]  = 0.25; /* strong threshold             */
    settings[9]  = 0.05; /* non galerkin tol             */
    settings[10] = 0;    /* aggressive coarsening levels */
    settings[11] = 1;    /* chebyRelaxOrder */

    platform->options.getArgs("BOOMERAMG COARSEN TYPE", settings[1]);
    platform->options.getArgs("BOOMERAMG INTERPOLATION TYPE", settings[2]);
    platform->options.getArgs("BOOMERAMG COARSE SMOOTHER TYPE", settings[4]);
    platform->options.getArgs("BOOMERAMG SMOOTHER TYPE", settings[6]);
    platform->options.getArgs("BOOMERAMG SMOOTHER SWEEPS", settings[7]);
    platform->options.getArgs("BOOMERAMG ITERATIONS", settings[3]);
    platform->options.getArgs("BOOMERAMG STRONG THRESHOLD", settings[8]);
    platform->options.getArgs("BOOMERAMG NONGALERKIN TOLERANCE" , settings[9]);
    platform->options.getArgs("BOOMERAMG AGGRESSIVE COARSENING LEVELS" , settings[10]);
    platform->options.getArgs("BOOMERAMG CHEBYSHEV RELAX ORDER" , settings[11]);

    if(useDevice) {
      boomerAMG = new hypreWrapperDevice::boomerAMG_t(
        N,
        nnz,
        Ai,
        Aj,
        Avals,
        (int) nullSpace,
        comm,
        platform->device.occaDevice(),
        useFP32,
        settings,
        verbose);
    } else {
      const int Nthreads = 1;
      boomerAMG = new hypreWrapper::boomerAMG_t(
        N,
        nnz,
        Ai,
        Aj,
        Avals,
        (int) nullSpace,
        comm,
        Nthreads,
        useFP32,
        settings,
        verbose);
    }
  }
  else if (options.compareArgs("COARSE SOLVER", "AMGX")){
    std::string configFile;
    options.getArgs("AMGX CONFIG FILE", configFile);
    char *cfg = NULL;
    if(configFile.size()) cfg = (char*) configFile.c_str();
    AMGX = new AMGX_t(
      N,
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
  } else {
    std::string amgSolver;
    options.getArgs("COARSE SOLVER", amgSolver);
    nrsAbort(platform->comm.mpiComm, EXIT_FAILURE,
             "COARSE SOLVER <%s> is not supported!\n", amgSolver.c_str());
  }

  MPI_Barrier(comm);
  if(rank==0) printf("done (%gs)\n", MPI_Wtime()-startTime);
}

MGSolver_t::coarseLevel_t::~coarseLevel_t()
{
  const auto useDevice = options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");
  if(boomerAMG) {
    if(useDevice) 
      delete (hypreWrapperDevice::boomerAMG_t*) this->boomerAMG;
    else
      delete (hypreWrapper::boomerAMG_t*) this->boomerAMG;
  }
  if(AMGX) delete AMGX;

  h_xBuffer.free();
  o_xBuffer.free();
  h_Sx.free();
  h_Gx.free();
  o_Sx.free();
  o_Gx.free();
}

void MGSolver_t::coarseLevel_t::solve(occa::memory& o_rhs, occa::memory& o_x) 
{
  platform->timer.tic("coarseSolve", 1);

  {
    const bool useDevice = options.compareArgs("COARSE SOLVER LOCATION", "DEVICE");

    const pfloat zero = 0.0;
    platform->linAlg->pfill(N, zero, o_xBuffer);
    if(!useDevice) o_xBuffer.copyTo(xBuffer, N*sizeof(pfloat));

    // E->T
    const pfloat one = 1.0;
    vectorDotStarKernel(ogs->N, one, zero, o_weight, o_rhs, o_Sx); 
    ogsGather(o_Gx, o_Sx, ogsPfloat, ogsAdd, ogs);
    if(!useDevice) o_Gx.copyTo(Gx, N*sizeof(pfloat));

    if (options.compareArgs("COARSE SOLVER", "BOOMERAMG")){
      if(useDevice) {
        auto boomerAMG = (hypreWrapperDevice::boomerAMG_t*) this->boomerAMG;
        boomerAMG->solve(o_Gx, o_xBuffer);
      } else {
        auto boomerAMG = (hypreWrapper::boomerAMG_t*) this->boomerAMG;
        boomerAMG->solve(Gx, xBuffer); 
      }
    } else if (options.compareArgs("COARSE SOLVER", "AMGX")){
        AMGX->solve(o_Gx.ptr(), o_xBuffer.ptr());
    }

    // T->E
    if(useDevice) {
      ogsScatter(o_x, o_xBuffer, ogsPfloat, ogsAdd, ogs);
    } else {
      o_Gx.copyFrom(xBuffer, N*sizeof(pfloat));
      ogsScatter(o_x, o_Gx, ogsPfloat, ogsAdd, ogs);
    }

  }

  platform->timer.toc("coarseSolve");
}
