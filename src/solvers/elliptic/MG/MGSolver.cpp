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

#include "MGSolver.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

namespace {

void coarsenV(MGSolver_t* M)
{
  for(int k = 0 ; k < M->numLevels-1; ++k){
    auto level = M->levels[k];
    auto o_rhs = level->o_rhs;
    auto o_res = level->o_res;
    auto levelC = M->levels[k + 1];
    auto o_rhsC = levelC->o_rhs;

    o_res.copyFrom(o_rhs, level->Nrows * sizeof(pfloat));
    levelC->coarsen(o_res, o_rhsC);
  }

}
void prolongateV(MGSolver_t* M)
{
  for(int k = M->numLevels-2; k >= 0; --k){
    auto level = M->levels[k];
    auto o_x = level->o_x;

    auto levelC = M->levels[k + 1];
    auto o_xC = levelC->o_x;

    // x = x + P xC
    levelC->prolongate(o_xC, o_x);
  }
}
void schwarzSolve(MGSolver_t* M)
{    
  for(int k = 0 ; k < M->numLevels-1; ++k){
    auto level = M->levels[k];
    auto o_rhs = level->o_rhs;
    auto o_x = level->o_x;
    auto o_res = level->o_res;
    auto levelC = M->levels[k + 1];
    auto o_rhsC = levelC->o_rhs;
    auto o_xC = levelC->o_x;

    // apply smoother to x and then compute res = rhs-Ax
    level->smooth(o_rhs, o_x, true);

    o_res.copyFrom(o_rhs, level->Nrows * sizeof(pfloat));

    // rhsC = P^T res
    levelC->coarsen(o_res, o_rhsC);
    }
}

}

MGSolver_t::MGSolver_t(occa::device device_, MPI_Comm comm_,
                       setupAide options_) 
{
  device = device_;
  comm = comm_;
  options = options_;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  levels = (MGSolver_t::multigridLevel **) calloc(MAX_LEVELS,sizeof(MGSolver_t::multigridLevel *));

  coarseLevel = new coarseLevel_t(options, comm);

  numLevels = 0;

  if(options.compareArgs("MGSOLVER CYCLE", "VCYCLE")) {
    ctype = VCYCLE;
    additive = false;
    if(options.compareArgs("MGSOLVER CYCLE", "ADDITIVE")) {
      if (options.compareArgs("MGSOLVER SMOOTHER", "CHEBYSHEV")) {
        if(rank==0) printf("Additive vcycle is not supported for Chebyshev!\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
      additive = true;
      overlapCrsGridSolve = false;
      if(options.compareArgs("MGSOLVER CYCLE", "OVERLAPCRS")){
        if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP"){
          overlapCrsGridSolve = false;
        } else {
          overlapCrsGridSolve = true;
          int provided;
          MPI_Query_thread(&provided);
          if(provided != MPI_THREAD_MULTIPLE) {
            overlapCrsGridSolve = false;
            if(rank ==0 && size > 1) 
              printf("disable overlapping coarse solve as MPI_THREAD_MULTIPLE is not supported!\n");
          }
          if(size == 1) overlapCrsGridSolve = true;
        }
        if(rank ==0 && overlapCrsGridSolve) printf("overlapping coarse grid solve enabled\n");
      }
    } else {
      if (options.compareArgs("MGSOLVER SMOOTHER", "RAS") || 
          options.compareArgs("MGSOLVER SMOOTHER", "ASM")) {
        if(!options.compareArgs("MGSOLVER SMOOTHER", "CHEBYSHEV")){
          if(rank==0) 
            printf("Multiplicative vcycle is not supported for RAS/ASM smoother without Chebyshev!\n");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
        }
      }
    }
  } else {
    if(rank==0) printf("Unknown multigrid cycle type!\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
  }
}

MGSolver_t::~MGSolver_t() {
  for (int n=0;n<numLevels;n++)
    delete levels[n];

  free(levels);

  if(coarseLevel) delete coarseLevel; 
}

void MGSolver_t::Run(occa::memory o_rhs, occa::memory o_x) 
{
  levels[0]->o_x   = o_x;
  levels[0]->o_rhs = o_rhs;

  if(ctype == VCYCLE) {
    if(additive)
      runAdditiveVcycle();
    else
      runVcycle(0);
  }
}

void MGSolver_t::Report() {

}

void MGSolver_t::runVcycle(int k)
{
  MGSolver_t::multigridLevel *level = levels[k];
  occa::memory o_rhs = level->o_rhs;
  occa::memory o_x   = level->o_x;
  occa::memory o_res = level->o_res;

  if(k == baseLevel) {
    coarseLevel->solvePtr(coarseLevel, o_rhs, o_x);
    return;
  }

  MGSolver_t::multigridLevel *levelC = levels[k+1];
  occa::memory o_rhsC = levelC->o_rhs;
  occa::memory o_xC   = levelC->o_x;

  level->smooth(o_rhs, o_x, true);
  level->residual(o_rhs, o_x, o_res);

  levelC->coarsen(o_res, o_rhsC);

  this->runVcycle(k+1);

  levelC->prolongate(o_xC, o_x);

  level->smooth(o_rhs, o_x, false);
}

void MGSolver_t::runAdditiveVcycle()
{
  {
    coarsenV(this);
  }

  const int nThreads = this->overlapCrsGridSolve ? 2 : 1;
  occa::memory o_rhs = levels[baseLevel]->o_rhs;
  occa::memory o_x   = levels[baseLevel]->o_x;

  auto xBuffer = this->coarseLevel->xBuffer;
  auto ogs = this->coarseLevel->ogs;

  auto Gx = this->coarseLevel->Gx;
  auto Sx = this->coarseLevel->Sx;

  // local E vector size
  const auto Nlocal = ogs->N;

  // local T vector size
  const auto NlocalT = this->coarseLevel->N;

  o_rhs.copyTo(Sx, Nlocal*sizeof(pfloat));

  o_x.getDevice().finish();
  #pragma omp parallel proc_bind(close) num_threads(nThreads)
  {
    #pragma omp single
    {
      #pragma omp task
      {
        schwarzSolve(this);
      }
      #pragma omp task
      {
        for(int i = 0; i < Nlocal; i++)
          Sx[i] *= this->coarseLevel->weight[i]; 
        ogsGather(Gx, Sx, ogsPfloat, ogsAdd, ogs);
    
        for(int i = 0; i < NlocalT; i++) {
          xBuffer[i] = 0; 
        }
 
        auto boomerAMG = (hypreWrapper::boomerAMG_t*) coarseLevel->boomerAMG;
        boomerAMG->solve(Gx, xBuffer);

        ogsScatter(Sx, xBuffer, ogsPfloat, ogsAdd, ogs);
      }
    }
  }

  o_x.copyFrom(Sx, Nlocal*sizeof(pfloat));

  {
    prolongateV(this);
  }
}
