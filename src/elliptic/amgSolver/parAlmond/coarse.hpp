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

namespace parAlmond {

class coarseSolver {

public:
  int coarseTotal;
  int coarseOffset;
  int *coarseOffsets=NULL;
  int *coarseCounts=NULL;

  int N;
  dfloat *invCoarseA=NULL;

  dfloat *xLocal=NULL;
  dfloat *rhsLocal=NULL;

  dfloat *xCoarse=NULL;
  dfloat *rhsCoarse=NULL;

  bool gatherLevel;
  ogs_t *ogs;
  dfloat *Gx, *Sx;
  occa::memory o_Sx, o_Gx;

  MPI_Comm comm;
  occa::device device;

  setupAide options;

  coarseSolver(setupAide options, MPI_Comm comm);
  ~coarseSolver();

  int getTargetSize();

  void setup(parCSR *A);
  void setup(dlong Nrows, hlong* globalRowStarts, dlong nnz, hlong* Ai, hlong* Aj, dfloat* Avals, bool nullSpace);

  void syncToDevice();

  void solve(dfloat *rhs, dfloat *x);
  void solve(occa::memory o_rhs, occa::memory o_x);
  void gather(occa::memory o_rhs, occa::memory o_x);
  void scatter(occa::memory o_rhs, occa::memory o_x);
  void BoomerAMGSolve();
  void AmgXSolve(occa::memory o_rhs, occa::memory o_x);
};

}

#endif
