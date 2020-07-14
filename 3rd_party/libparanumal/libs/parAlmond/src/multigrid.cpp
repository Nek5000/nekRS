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
#include <omp.h>

namespace parAlmond {

void solver_t::kcycle(int k){

  multigridLevel *level = levels[k];

  dlong m = level->Nrows;

  dfloat* rhs = level->rhs;
  dfloat*   x = level->x;
  dfloat* res = level->res;

  //check for base level
  if(k==baseLevel) {
    if(options.compareArgs("PARALMOND SMOOTH COARSEST", "TRUE") &&
       !options.compareArgs("AMG SOLVER", "AMG"))
      level->smooth(rhs,x,true);
    else
      coarseLevel->solve(rhs, x);

    return;
  }

  multigridLevel *levelC = levels[k+1];
  dlong mCoarse = levelC->Nrows;
  dfloat* rhsC   = levelC->rhs;
  dfloat*   xC   = levelC->x;

  //apply smoother to x and then return res = rhs-Ax
  level->smooth(rhs, x, true);
  level->residual(rhs, x, res);

  // rhsC = P^T res
  levelC->coarsen(res, rhsC);

  if(k+1>NUMKCYCLES) {
    this->vcycle(k+1);
  } else{
    // first inner krylov iteration
    this->kcycle(k+1);

    // ck = x
    // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
    // rhsC = rhsC - (alpha1/rho1)*vkp1
    // norm_rtilde = sqrt(rhsC*rhsC)
    dfloat rho1, alpha1, norm_rhs, norm_rhstilde;
    levelC->kcycleOp1(&alpha1, &rho1, &norm_rhs, &norm_rhstilde);

    if(norm_rhstilde < KCYCLETOL*norm_rhs){
      // xC = (alpha1/rho1)*xC
      vectorScale(mCoarse, alpha1/rho1, xC);
    } else{

      // second inner krylov iteration
      this->kcycle(k+1);

      // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
      // rho2=beta - gamma*gamma/rho1
      // xC = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*xC
      levelC->kcycleOp2(alpha1, rho1);
    }
  }

  // x = x + P xC
  levelC->prolongate(xC, x);

  level->smooth(rhs, x, false);
}


void solver_t::device_kcycle(int k){

  multigridLevel *level = levels[k];

  dlong m = level->Nrows;

  occa::memory o_rhs = level->o_rhs;
  occa::memory o_x   = level->o_x;
  occa::memory o_res = level->o_res;

  //check for device<->host handoff
  if(m < GPU_CPU_SWITCH_SIZE){
    o_rhs.copyTo(level->rhs, m*sizeof(dfloat));
    this->kcycle(k);
    o_x.copyFrom(level->x, m*sizeof(dfloat));
    return;
  }

  //check for base level
  if(k==baseLevel) {
    //    coarseLevel->solve(o_rhs, o_x);

    if(options.compareArgs("PARALMOND SMOOTH COARSEST", "TRUE"))
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

  if(k+1>NUMKCYCLES) {
    this->device_vcycle(k+1);
  } else{
    // first inner krylov iteration
    this->device_kcycle(k+1);

    // alpha1=ck*rhsC, rho1=ck*Ack, norm_rhs=sqrt(rhsC*rhsC)
    // rhsC = rhsC - (alpha1/rho1)*vkp1
    // norm_rtilde = sqrt(rhsC*rhsC)
    dfloat rho1, alpha1, norm_rhs, norm_rhstilde;
    levelC->device_kcycleOp1(&alpha1, &rho1, &norm_rhs, &norm_rhstilde);

    if(norm_rhstilde < KCYCLETOL*norm_rhs){
      // xC = (alpha1/rho1)*xC
      vectorScale(mCoarse, alpha1/rho1, o_xC);
    } else{

      // second inner krylov iteration
      this->device_kcycle(k+1);

      // gamma=xC*Ack, beta=xC*AxC, alpha2=xC*rhsC
      // rho2=beta - gamma*gamma/rho1
      // xC = (alpha1/rho1 - (gam*alpha2)/(rho1*rho2))*ck + (alpha2/rho2)*xC
      levelC->device_kcycleOp2(alpha1, rho1);
    }
  }

  // x = x + P xC
  levelC->prolongate(o_xC, o_x);
  level->smooth(o_rhs, o_x, false);
}



void solver_t::vcycle(int k) {

  multigridLevel *level = levels[k];

  dlong m = level->Nrows;

  dfloat* rhs = level->rhs;
  dfloat*   x = level->x;
  dfloat* res = level->res;

  //check for base level
  if(k==baseLevel) {
    //    coarseLevel->solve(rhs, x);

    if(options.compareArgs("PARALMOND SMOOTH COARSEST", "TRUE"))
      level->smooth(rhs,x,true);
    else
      coarseLevel->solve(rhs, x);
    
    return;
  }

  multigridLevel *levelC = levels[k+1];
  dlong mCoarse = levelC->Nrows;
  dfloat* rhsC   = levelC->rhs;
  dfloat*   xC   = levelC->x;

  //apply smoother to x and then return res = rhs-Ax
  level->smooth(rhs, x, true);
  level->residual(rhs, x, res);

  // rhsC = P^T res
  levelC->coarsen(res, rhsC);

  this->vcycle(k+1);

  // x = x + P xC
  levelC->prolongate(xC, x);

  level->smooth(rhs, x, false);
}


void solver_t::device_vcycle(int k){

  multigridLevel *level = levels[k];

  dlong m = level->Nrows;

  occa::memory o_rhs = level->o_rhs;
  occa::memory o_x   = level->o_x;
  occa::memory o_res = level->o_res;

  //check for device<->host handoff
  if(m < GPU_CPU_SWITCH_SIZE){
    o_rhs.copyTo(level->rhs, m*sizeof(dfloat));
    vcycle(k);
    o_x.copyFrom(level->x, m*sizeof(dfloat));
    return;
  }

  //check for base level
  if(k==baseLevel) {
    //    coarseLevel->solve(o_rhs, o_x);

    if(options.compareArgs("PARALMOND SMOOTH COARSEST", "TRUE")){
      level->smooth(o_rhs,o_x,true);
    }
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
void coarseSolve(solver_t* M)
{   
  M->coarseLevel->BoomerAMGSolve();
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

  coarseLevel->gather(o_rhs, o_x);
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
        coarseSolve(this);
      }
    }
  }
  o_x.getDevice().finish();
  coarseLevel->scatter(o_rhs, o_x);
  o_x.getDevice().finish();

  {
    prolongateV(this);
  }

}

} //hamespace parAlmond
