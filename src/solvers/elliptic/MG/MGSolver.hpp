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

#ifndef MGSOLVER_HPP
#define MGSOLVER_HPP

#include <functional>

#include "nrssys.hpp"
#include "defines.hpp"

#include "parseMultigridSchedule.hpp"
#include "determineMGLevels.hpp"

#include "hypreWrapper.hpp"
#include "hypreWrapperDevice.hpp"
#include "AMGX.hpp"

class MGSolver_t {

public:
  class multigridLevel {
  
  public:
    dlong Nrows, Ncols;
  
    occa::memory o_x, o_rhs, o_res;
  
    SmoothType smootherType;
  
    MPI_Comm comm;
  
    multigridLevel(dlong N, dlong M, MPI_Comm comm);
    ~multigridLevel();
  
    virtual void Ax(occa::memory o_x, occa::memory o_Ax)=0;
  
    virtual void smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero)=0;
  
    virtual void residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res)=0;
  
    virtual void coarsen(occa::memory o_x, occa::memory o_Cx)=0;
  
    virtual void prolongate(occa::memory o_x, occa::memory o_Px)=0;
  
    virtual void Report()=0;
  };

  class coarseLevel_t {
  
  public:
    int N;
  
    occa::memory h_xBuffer;
    occa::memory o_xBuffer;
    pfloat *xBuffer;
  
    ogs_t *ogs;
    pfloat *Gx, *Sx;
    occa::memory h_Sx, h_Gx;
    occa::memory o_Sx, o_Gx;
  
    pfloat *weight = NULL;
    occa::memory o_weight;
  
    MPI_Comm comm;
    occa::device device;
  
    setupAide options;
  
    coarseLevel_t(setupAide options, MPI_Comm comm);
    ~coarseLevel_t();

    void setupSolver(hlong* globalRowStarts, dlong nnz, hlong* Ai, hlong* Aj, dfloat* Avals, bool nullSpace);
    void solve(occa::memory& o_rhs, occa::memory& o_x);
    std::function<void(coarseLevel_t *, occa::memory&, occa::memory&)> solvePtr = nullptr;
 
    void *boomerAMG = nullptr;
    AMGX_t *AMGX = nullptr;
  
  };


  MPI_Comm comm;
  int rank, size;

  occa::device device;
  setupAide options;

  CycleType ctype;
  SmoothType smootherType;

  int numLevels;
  int AMGstartLev, baseLevel;
  multigridLevel **levels = nullptr;

  coarseLevel_t *coarseLevel = nullptr;

  bool additive;
  bool overlapCrsGridSolve;

  MGSolver_t(occa::device otherdevice, MPI_Comm othercomm,
           setupAide otheroptions);

  ~MGSolver_t();

  void Run(occa::memory o_rhs, occa::memory o_x);

  void Report();


private:
  void runVcycle(int k);
  void runAdditiveVcycle();


};

#endif
