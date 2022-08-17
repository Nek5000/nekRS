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
#include "hypreWrapper.hpp"
#include <omp.h>

namespace parAlmond {

void solver_t::device_vcycle(int k){

  multigridLevel *level = levels[k];

  dlong m = level->Nrows;

  occa::memory o_rhs = level->o_rhs;
  occa::memory o_x   = level->o_x;
  occa::memory o_res = level->o_res;

  if(k==baseLevel) {
    if(options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
      level->smooth(o_rhs,o_x,true);
    else
      coarseLevel->solve(o_rhs, o_x);
    
    return;
  }

  multigridLevel *levelC = levels[k+1];
  dlong mCoarse = levelC->Nrows;
  occa::memory o_rhsC = levelC->o_rhs;
  occa::memory o_xC   = levelC->o_x;

  //apply smoother to x and then compute res = rhs-Ax
  level->smooth(o_rhs, o_x, true);
  level->residual(o_rhs, o_x, o_res);

  // rhsC = P^T res
  levelC->coarsen(o_res, o_rhsC);

  this->device_vcycle(k+1);

  // x = x + P xC
  levelC->prolongate(o_xC, o_x);

  level->smooth(o_rhs, o_x, false);
}

namespace {

void coarsenV(solver_t* M)
{
  for(int k = 0 ; k < M->numLevels-1; ++k){
    multigridLevel *level = M->levels[k];
    occa::memory o_rhs = level->o_rhs;
    occa::memory o_x   = level->o_x;
    occa::memory o_res = level->o_res;
    multigridLevel *levelC = M->levels[k+1];
    occa::memory o_rhsC = levelC->o_rhs;
    occa::memory o_xC   = levelC->o_x;
    level->residual(o_rhs, o_x, o_res);
    levelC->coarsen(o_res, o_rhsC);
  }

}
void prolongateV(solver_t* M)
{
  for(int k = M->numLevels-2; k >= 0; --k){
    multigridLevel *level = M->levels[k];
    occa::memory o_rhs = level->o_rhs;
    occa::memory o_x   = level->o_x;
    occa::memory o_res = level->o_res;
    multigridLevel *levelC = M->levels[k+1];
    dlong mCoarse = levelC->Nrows;
    occa::memory o_rhsC = levelC->o_rhs;
    occa::memory o_xC   = levelC->o_x;
    // x = x + P xC
    levelC->prolongate(o_xC, o_x);
    level->smooth(o_rhs, o_x, false);
  }
}
void schwarzSolve(solver_t* M)
{    
  for(int k = 0 ; k < M->numLevels-1; ++k){
      multigridLevel *level = M->levels[k];
      occa::memory o_rhs = level->o_rhs;
      occa::memory o_x   = level->o_x;
      occa::memory o_res = level->o_res;
      multigridLevel *levelC = M->levels[k+1];
      occa::memory o_rhsC = levelC->o_rhs;
      occa::memory o_xC   = levelC->o_x;

      //apply smoother to x and then compute res = rhs-Ax
      level->smooth(o_rhs, o_x, true);
      level->residual(o_rhs, o_x, o_res);

      // rhsC = P^T res
      levelC->coarsen(o_res, o_rhsC);
    }
}
}

void solver_t::additiveVcycle()
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
  const auto N = this->coarseLevel->N;

  o_rhs.copyTo(Sx, Nlocal*sizeof(pfloat));

  o_x.getDevice().finish();
  #pragma omp parallel proc_bind(close) num_threads(nThreads)
  {
    #pragma omp single
    {
      #pragma omp task
      {
        //printf("Schwarz solve omp thread%d\n", omp_get_thread_num());
        schwarzSolve(this);
      }
      #pragma omp task
      {
        //printf("Coarse solve omp thread %d\n", omp_get_thread_num());

        for(int i = 0; i < Nlocal; i++)
          Sx[i] *= this->coarseLevel->weight[i]; 
        ogsGather(Gx, Sx, ogsPfloat, ogsAdd, ogs);
    
        for(int i = 0; i < N; i++) {
          xBuffer[i] = 0; 
        }

        hypreWrapper::BoomerAMGSolve(Gx, xBuffer);

        ogsScatter(Sx, xBuffer, ogsPfloat, ogsAdd, ogs);
      }
    }
  }

  o_x.copyFrom(Sx, Nlocal*sizeof(pfloat));

  {
    prolongateV(this);
  }
}

} //hamespace parAlmond
