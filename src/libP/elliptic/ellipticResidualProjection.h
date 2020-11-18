/*

   The MIT License (MIT)

   Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#ifndef ELLIPTIC_RESIDUAL_PROJECTION_H
#define ELLIPTIC_RESIDUAL_PROJECTION_H
#include <occa.hpp>
#include <types.h>
#include <vector>
#include <sstream>
#include <elliptic.h>
#include <functional>

class ResidualProjection final
{
public:
  ResidualProjection(elliptic_t& _elliptic,
                     const dlong _maxNumVecsProjection = 8,
                     const dlong _numTimeSteps = 5);
  void pre(occa::memory& o_r);
  void post(occa::memory& o_x);
private:
  void computePreProjection(occa::memory& o_r);
  void computePostProjection(occa::memory& o_x);
  void updateProjectionSpace();
  void matvec(occa::memory& o_Ax, const dlong Ax_offset, occa::memory& o_x, const dlong x_offset);
  void multiWeightedInnerProduct(
                              occa::memory& o_a,
                              const dlong m,
                              occa::memory& o_b,
                              const dlong offset);
  const dlong maxNumVecsProjection;
  const dlong numTimeSteps;
  dlong timestep;
  bool verbose;

  occa::memory o_xbar;
  occa::memory o_xx;
  occa::memory o_bb;
  occa::memory o_alpha;
  // references to memory on elliptic
  occa::memory& o_invDegree;
  occa::memory& o_rtmp;
  occa::memory& o_Ap;
  occa::memory& o_tmp;
  occa::memory& o_tmp2;
  occa::memory& o_wrk;

  occa::kernel scalarMultiplyKernel;
  occa::kernel multiScaledAddwOffsetKernel;
  occa::kernel multiWeightedInnerProduct2Kernel;
  occa::kernel accumulateKernel;
  occa::kernel scaledAddKernel;
  occa::kernel sumKernel;

  dfloat * alpha;
  dfloat * work;
  dfloat * multiwork;
  dfloat* tmp;

  dlong numVecsProjection;
  const dlong Nlocal; // vector size
  const dlong fieldOffset; // offset
  const dlong Nblock;
  const dlong Nblock2;
  const dfloat resNormFactor;
  const int rank;
  const int size;
  MPI_Comm comm;

  std::function<void(occa::memory&,occa::memory&)> matvecOperator;
  std::function<dfloat(occa::memory&)> weightedNorm;


};
#endif