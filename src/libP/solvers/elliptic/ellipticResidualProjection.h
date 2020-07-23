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

class ResidualProjection final
{
public:
  ResidualProjection(elliptic_t& _elliptic,
                     const dlong _maxNumVecsProjection = 8,
                     const dlong _numTimeSteps = 5);
  void preSolveProjection(occa::memory& o_r);
  void postSolveProjection(occa::memory& o_x);
private:
  void computePreProjection(occa::memory& o_r);
  void computePostProjection(occa::memory& o_x);
  void updateProjectionSpace();
  void reOrthogonalize();
  bool checkOrthogonalize();
  void matvec(occa::memory& o_Ax, const dlong Ax_offset, occa::memory& o_x, const dlong x_offset);
  void gop(dfloat*,dfloat*, const dlong);
  dfloat computeInnerProduct(occa::memory& o_a,
                             const dlong a_offset,
                             occa::memory& o_b,
                             const dlong b_offset);
  dfloat weightedInnerProduct(occa::memory& o_w,
                              occa::memory& o_a,
                              const dlong a_offset,
                              occa::memory& o_b,
                              const dlong b_offset);
  elliptic_t& elliptic;
  const dlong maxNumVecsProjection;
  const dlong numTimeSteps;
  dlong timestep;
  bool initialized;
  bool verbose;

  occa::memory o_xbar;
  occa::memory o_xx;
  occa::memory o_bb;
  occa::memory o_alpha;

  occa::kernel scalarMultiplyKernel;
  occa::kernel scaledAddwOffsetTwoVecKernel;
  occa::kernel scaledAddwOffsetKernel;
  occa::kernel placeVectorKernel;
  occa::kernel extractVectorKernel;
  occa::kernel weightedInnerProduct2Kernel;
  occa::kernel innerProductKernel;
  occa::kernel accumulateKernel;

  std::vector<dfloat> alpha; // host shadow
  std::vector<dfloat> work; // O(m) work array

  dlong numVecsProjection;
  dlong Ntotal; // vector size
  bool useWeightedFormulation;

  // logging...
  std::ostringstream log;
};
#endif