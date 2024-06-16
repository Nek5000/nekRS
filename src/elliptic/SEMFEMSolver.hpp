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

#include "elliptic.h"

#include "hypreWrapper.hpp"
#include "hypreWrapperDevice.hpp"
#include "AMGX.hpp"

class SEMFEMSolver_t
{

public:
  SEMFEMSolver_t(elliptic_t *);
  ~SEMFEMSolver_t();

  void run(const occa::memory &, occa::memory &);

private:
  occa::memory o_dofMap;

  void *boomerAMG = nullptr;
  AMGX_t *AMGX = nullptr;
  elliptic_t *elliptic = nullptr;

  struct matrix_t {
    int nnz;
    long long rowStart;
    long long rowEnd;

    std::vector<long long> Ai;
    std::vector<long long> Aj;
    std::vector<double> Av;
    std::vector<dlong> dofMap;
  };

  matrix_t buildMatrix(const int N_,
                       const int n_elem_,
                       occa::memory _o_x,
                       occa::memory _o_y,
                       occa::memory _o_z,
                       const std::vector<int>& pmask_,
                       double lambda,
                       void *gsh);
};

#endif
