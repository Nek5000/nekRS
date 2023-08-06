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

#ifndef SEMFEMSOLVER_HPP
#define SEMFEMSOLVER_HPP

#include "nrssys.hpp"
#include "elliptic.h"

#include "hypreWrapper.hpp"
#include "hypreWrapperDevice.hpp"
#include "AMGX.hpp"

class SEMFEMSolver_t {

public:
  SEMFEMSolver_t(elliptic_t*);
  ~SEMFEMSolver_t();

  void run(occa::memory&, occa::memory&);

private:
  dlong numRowsSEMFEM;
  occa::memory o_dofMap;
  occa::memory o_SEMFEMBuffer1;
  occa::memory o_SEMFEMBuffer2;
  void *SEMFEMBuffer1_h_d;
  void *SEMFEMBuffer2_h_d;
  void *boomerAMG = nullptr;
  AMGX_t *AMGX = nullptr;

  elliptic_t *elliptic;

  struct matrix_t {
    long long *Ai;
    long long *Aj;
    double *Av;
    int nnz;
    long long rowStart;
    long long rowEnd;
    long long *dofMap;
  };
  
  matrix_t *build(const int N_,
                  const int n_elem_,
                  occa::memory _o_x,
                  occa::memory _o_y,
                  occa::memory _o_z,
                  double *pmask_,
                  double lambda,
                  hypreWrapper::IJ_t &hypreIJ,
                  MPI_Comm comm,
                  long long int *gatherGlobalNodes);

};

#endif
