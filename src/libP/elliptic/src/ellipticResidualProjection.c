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
#include "elliptic.h"
#include "ellipticResidualProjection.h"
#include <iostream>
#include <timer.hpp>

void ResidualProjection::matvec(occa::memory& o_Ax,
                                const dlong Ax_offset,
                                occa::memory& o_x,
                                const dlong x_offset)
{
  occa::memory o_xtmp = o_x + fieldOffset * x_offset * sizeof(dfloat);
  occa::memory o_Axtmp = o_Ax + fieldOffset * Ax_offset * sizeof(dfloat);
  matvecOperator(o_xtmp, o_Axtmp);
}
void ResidualProjection::updateProjectionSpace()
{
  dlong m = numVecsProjection;
  if(m <= 0) return;

  multiWeightedInnerProduct(o_xx, m, o_bb, m-1);
  const dfloat norm_orig = alpha[m - 1];
  dfloat norm_new = norm_orig;
  const dfloat one = 1.0;
  multiScaledAddwOffsetKernel(Nlocal, m, (m-1)*fieldOffset, fieldOffset, o_alpha, one, o_xx);
  multiScaledAddwOffsetKernel(Nlocal, m, (m-1)*fieldOffset, fieldOffset, o_alpha, one, o_bb);
  for(int k = 0; k < m - 1; ++k) {
    norm_new = norm_new - alpha[k] * alpha[k];
  }
  norm_new = sqrt(norm_new);
  dfloat tol = 1e-7;
  const dfloat test = norm_new / norm_orig;
  if(test > tol) {
    const dfloat scale = 1.0 / norm_new;
    scalarMultiplyKernel(Nlocal, (m-1) * fieldOffset, scale, o_xx);
    scalarMultiplyKernel(Nlocal, (m-1) * fieldOffset, scale, o_bb);
  } else {
    if(verbose && rank == 0){
      std::cout << "Detected rank deficiency: " << test << ".\n";
      std::cout << "Removing column : " << numVecsProjection << ".\n";
    }
    numVecsProjection--;
  }
}
void ResidualProjection::computePreProjection(occa::memory& o_r)
{
  dfloat one = 1.0;
  dfloat zero = 0.0;
  dfloat mone = -1.0;
  const int m = numVecsProjection;
  if(m <= 0) return;
  multiWeightedInnerProduct(o_xx,m,o_r,0);

  accumulateKernel(Nlocal, m, fieldOffset, o_alpha, o_xx, o_xbar);
  accumulateKernel(Nlocal, m, fieldOffset, o_alpha, o_bb, o_rtmp);
  scaledAddKernel(Nlocal, mone, o_rtmp, one, o_r);
}
void ResidualProjection::computePostProjection(occa::memory & o_x)
{
  const dfloat one = 1.0;
  const dfloat zero = 0.0;

  if(numVecsProjection == 0) {
    // reset bases
    numVecsProjection = 1;
    o_xx.copyFrom(o_x, Nlocal * sizeof(dfloat));
  } else if(numVecsProjection == maxNumVecsProjection) {
    numVecsProjection = 1;
    scaledAddKernel(Nlocal, one, o_xbar, one, o_x);
    o_xx.copyFrom(o_x, Nlocal * sizeof(dfloat));
  } else {
    numVecsProjection++;
    // xx[m-1] = x
    o_xx.copyFrom(o_x, Nlocal*sizeof(dfloat), (numVecsProjection - 1)*fieldOffset*sizeof(dfloat), 0);
    // x = x + xbar
    scaledAddKernel(Nlocal, one, o_xbar, one, o_x);
  }
  const dlong m_save = numVecsProjection;
  matvec(o_bb,numVecsProjection - 1,o_xx,numVecsProjection - 1);

  updateProjectionSpace();
  if (numVecsProjection < m_save) { // Last vector was linearly dependent, reset space
    numVecsProjection = 1;
    o_xx.copyFrom(o_x, Nlocal * sizeof(dfloat)); // writes first n words of o_xx, first approximation vector
    matvec(o_bb,0,o_xx,0);
    updateProjectionSpace();
  }
}
ResidualProjection::ResidualProjection(elliptic_t& elliptic,
                                       const dlong _maxNumVecsProjection,
                                       const dlong _numTimeSteps)
  :
  maxNumVecsProjection(_maxNumVecsProjection),
  numTimeSteps(_numTimeSteps),
  Nlocal(elliptic.mesh->Np * elliptic.mesh->Nelements),
  fieldOffset(elliptic.Ntotal),
  Nblock(elliptic.Nblock),
  Nblock2(elliptic.Nblock2),
  resNormFactor(elliptic.resNormFactor),
  rank(elliptic.mesh->rank),
  size(elliptic.mesh->size),
  comm(elliptic.mesh->comm),
  o_tmp(elliptic.o_tmp),
  o_tmp2(elliptic.o_tmp2),
  o_wrk(elliptic.o_wrk),
  o_invDegree(elliptic.mesh->ogs->o_invDegree),
  o_rtmp(elliptic.o_rtmp),
  o_Ap(elliptic.o_Ap)
{
  tmp = elliptic.tmp;
  timestep = 0;
  const dlong Nblock = elliptic.Nblock;
  numVecsProjection = 0;
  verbose = elliptic.options.compareArgs("VERBOSE","TRUE");
  alpha = (dfloat*) calloc(maxNumVecsProjection, sizeof(dfloat));
  work = (dfloat*) calloc(maxNumVecsProjection, sizeof(dfloat));
  multiwork = (dfloat*) calloc(Nblock*maxNumVecsProjection, sizeof(dfloat));
  o_alpha = elliptic.mesh->device.malloc(maxNumVecsProjection * sizeof(dfloat));
  o_xbar = elliptic.mesh->device.malloc(Nlocal * sizeof(dfloat));
  o_xx = elliptic.mesh->device.malloc(fieldOffset * maxNumVecsProjection * sizeof(dfloat));
  o_bb = elliptic.mesh->device.malloc(fieldOffset * maxNumVecsProjection * sizeof(dfloat));

  char fileName[BUFSIZ], kernelName[BUFSIZ];
  for (int r = 0; r < 2; r++) {
    if ((r == 0 && elliptic.mesh->rank == 0) || (r == 1 && elliptic.mesh->rank > 0)) {
      occa::properties properties;
      properties += elliptic.mesh->device.properties();
      properties["defines/p_threadBlockSize"] = blockSize;
      properties["defines/p_blockSize"] = blockSize;
      properties["defines/dfloat"] = dfloatString;
      properties["defines/dlong"] = dlongString;

      sprintf(fileName, DELLIPTIC "/okl/ellipticResidualProjection.okl");
      scalarMultiplyKernel = elliptic.mesh->device.buildKernel(fileName,
                                                               "scalarMultiply",
                                                               properties);
      multiScaledAddwOffsetKernel = elliptic.mesh->device.buildKernel(fileName,
                                                                 "multiScaledAddwOffset",
                                                                 properties);
      multiWeightedInnerProduct2Kernel = elliptic.mesh->device.buildKernel(fileName,
                                                                      "multiWeightedInnerProduct2",
                                                                      properties);
      accumulateKernel = elliptic.mesh->device.buildKernel(fileName, "accumulate", properties);
    }
    MPI_Barrier(elliptic.mesh->comm);
  }
  scaledAddKernel = elliptic.scaledAddKernel;
  sumKernel = elliptic.mesh->sumKernel;
  matvecOperator = [&](occa::memory& o_x, occa::memory & o_Ax)
  {
    ellipticOperator(&elliptic, o_x, o_Ax, dfloatString);
  };
  weightedNorm = [&](occa::memory& o_x)
  {
    return ellipticWeightedNorm2(&elliptic, o_invDegree, o_x);
  };
}
void ResidualProjection::pre(occa::memory& o_r)
{
  ++timestep;
  if(timestep < numTimeSteps)
    return;

  const int m = numVecsProjection;
  if(m <= 0) return;
  dfloat priorResidualNorm = 0.0;
  if(verbose) priorResidualNorm =
      sqrt(weightedNorm(o_r) * resNormFactor);
  computePreProjection(o_r);
  dfloat postResidualNorm = 0.0;
  dfloat ratio = 0.0;
  if(verbose) {
    postResidualNorm = sqrt(weightedNorm(o_r) * resNormFactor);
    ratio = priorResidualNorm / postResidualNorm;
  }
  if(rank == 0 && verbose)
    std::cout << "Residual projection : "
              << std::cout.precision(15)
              << priorResidualNorm << ", "
              << postResidualNorm << ", "
              << ratio << "\n";
}
void ResidualProjection::post(occa::memory& o_x)
{
  if(timestep < numTimeSteps)
    return;
  computePostProjection(o_x);
}
void ResidualProjection::multiWeightedInnerProduct(
                                                occa::memory &o_a,
                                                const dlong m,
                                                occa::memory &o_b,
                                                const dlong offset)
{
#ifdef ELLIPTIC_ENABLE_TIMER
  timer::tic("dotp",1);
#endif
  multiWeightedInnerProduct2Kernel(Nlocal, fieldOffset, Nblock, m, offset*fieldOffset, o_invDegree, o_a, o_b, o_wrk);

  o_wrk.copyTo(multiwork, sizeof(dfloat)*m*Nblock);
  for(dlong k = 0 ; k < m; ++k){
    dfloat accum = 0.0;
    for(dlong n = 0; n < Nblock; ++n){
      accum += multiwork[n+k*Nblock];
    }
    alpha[k] = accum;
  }
  MPI_Allreduce(MPI_IN_PLACE, alpha, m, MPI_DFLOAT, MPI_SUM, comm);
  o_alpha.copyFrom(alpha,sizeof(dfloat)*m);
#ifdef ELLIPTIC_ENABLE_TIMER
  timer::toc("dotp");
#endif
}
