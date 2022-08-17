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

#ifndef PARALMOND_COARSESOLVE_HPP
#define PARALMOND_COARSESOLVE_HPP

#include <functional>

namespace parAlmond {

class coarseSolver {

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
  std::function<void(occa::memory,occa::memory)> semfemSolver = nullptr;

  coarseSolver(setupAide options, MPI_Comm comm);
  ~coarseSolver();

  void setup(dlong Nrows, hlong* globalRowStarts, dlong nnz, hlong* Ai, hlong* Aj, dfloat* Avals, bool nullSpace);

  void solve(occa::memory o_rhs, occa::memory o_x);
};

}

#endif
