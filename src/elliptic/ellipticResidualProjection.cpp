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
#include "mesh.h"
#include "elliptic.h"
#include "ellipticResidualProjection.h"
#include <iostream>
#include "timer.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

void ResidualProjection::matvec(occa::memory& o_Ax,
                                const dlong Ax_offset,
                                occa::memory& o_x,
                                const dlong x_offset)
{
  occa::memory o_xtmp = o_x + Nfields * fieldOffset * x_offset * sizeof(dfloat);
  occa::memory o_Axtmp = o_Ax + Nfields * fieldOffset * Ax_offset * sizeof(dfloat);
  matvecOperator(o_xtmp, o_Axtmp);
}

void ResidualProjection::updateProjectionSpace()
{
  
  if(numVecsProjection <= 0) return;

  platform->linAlg->weightedInnerProdMulti(
    Nlocal,
    numVecsProjection,
    Nfields,
    fieldOffset,
    o_invDegree,
    o_xx,
    o_bb,
    platform->comm.mpiComm,
    alpha,
    (type == ProjectionType::CLASSIC) ?
      Nfields * (numVecsProjection-1) * fieldOffset :
      0
  );
  o_alpha.copyFrom(alpha,sizeof(dfloat) * numVecsProjection);

  const dfloat norm_orig = alpha[numVecsProjection - 1];
  dfloat norm_new = norm_orig;
  const dfloat one = 1.0;
  multiScaledAddwOffsetKernel(Nlocal, numVecsProjection, Nfields * (numVecsProjection - 1) * fieldOffset, fieldOffset, o_alpha, one, o_xx);
  if(type == ProjectionType::CLASSIC) multiScaledAddwOffsetKernel(Nlocal, numVecsProjection, Nfields * (numVecsProjection - 1) * fieldOffset, fieldOffset, o_alpha, one, o_bb);
  for(int k = 0; k < numVecsProjection - 1; ++k)
    norm_new = norm_new - alpha[k] * alpha[k];
  norm_new = sqrt(norm_new);
  dfloat tol = 1e-7;
  const dfloat test = norm_new / norm_orig;
  if(test > tol) {
    const dfloat scale = 1.0 / norm_new;
    platform->linAlg->scaleMany(Nlocal, Nfields, fieldOffset, scale, o_xx, fieldOffset * Nfields * (numVecsProjection - 1));
    if(type == ProjectionType::CLASSIC) platform->linAlg->scaleMany(Nlocal, Nfields, fieldOffset, scale, o_bb, fieldOffset * Nfields * (numVecsProjection - 1));
  } else {
    if(verbose && platform->comm.mpiRank == 0) {
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
  if(numVecsProjection <= 0) return;
  platform->linAlg->weightedInnerProdMulti(
    Nlocal,
    numVecsProjection,
    Nfields,
    fieldOffset,
    o_invDegree,
    o_xx,
    o_r,
    platform->comm.mpiComm,
    alpha,
    Nfields * 0 * fieldOffset
  );
  o_alpha.copyFrom(alpha,sizeof(dfloat) * numVecsProjection);

  accumulateKernel(Nlocal, numVecsProjection, fieldOffset, o_alpha, o_xx, o_xbar);
  if(type == ProjectionType::CLASSIC){
    accumulateKernel(Nlocal, numVecsProjection, fieldOffset, o_alpha, o_bb, o_rtmp);
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, mone, o_rtmp, one, o_r);
  }
  else if (type == ProjectionType::ACONJ)
  {
    matvec(o_bb, 0, o_xbar, 0);
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, mone, o_bb, one, o_r);
  }
}

void ResidualProjection::computePostProjection(occa::memory & o_x)
{
  
  const dfloat one = 1.0;
  const dfloat zero = 0.0;

  if(numVecsProjection == 0) {
    // reset bases
    numVecsProjection = 1;
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat));
  } else if(numVecsProjection == maxNumVecsProjection) {
    numVecsProjection = 1;
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, one, o_xbar, one, o_x);
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat));
  } else {
    numVecsProjection++;
    // xx[m-1] = x
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat), Nfields * (numVecsProjection - 1) * fieldOffset * sizeof(dfloat), 0);
    // x = x + xbar
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, one, o_xbar, one, o_x);
  }
  const dlong previousNumVecsProjection = numVecsProjection;
  const dlong bOffset = (type == ProjectionType::CLASSIC) ? numVecsProjection - 1 : 0;
  matvec(o_bb,bOffset,o_xx,numVecsProjection - 1);

  updateProjectionSpace();
  if (numVecsProjection < previousNumVecsProjection) { // Last vector was linearly dependent, reset space
    numVecsProjection = 1;
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat)); // writes first n words of o_xx, first approximation vector
    matvec(o_bb,0,o_xx,0);
    updateProjectionSpace();
  }
}

ResidualProjection::ResidualProjection(elliptic_t& elliptic,
                                       const ProjectionType _type,
                                       const dlong _maxNumVecsProjection,
                                       const dlong _numTimeSteps)
  :
  maxNumVecsProjection(_maxNumVecsProjection),
  numTimeSteps(_numTimeSteps),
  type(_type),
  Nlocal(elliptic.mesh->Np * elliptic.mesh->Nelements),
  fieldOffset(elliptic.Ntotal),
  Nfields(elliptic.Nfields),
  o_invDegree(elliptic.mesh->ogs->o_invDegree),
  o_rtmp(elliptic.o_rtmp),
  o_Ap(elliptic.o_Ap)
{
  platform_t* platform = platform_t::getInstance();
  timestep = 0;
  numVecsProjection = 0;
  verbose = elliptic.options.compareArgs("VERBOSE","TRUE");
  alpha = (dfloat*) calloc(maxNumVecsProjection, sizeof(dfloat));
  o_alpha = platform->device.malloc(maxNumVecsProjection, sizeof(dfloat));
  o_xbar = platform->device.malloc(Nfields * fieldOffset, sizeof(dfloat));
  o_xx = platform->device.malloc(Nfields * fieldOffset * maxNumVecsProjection, sizeof(dfloat));
  o_bb = platform->device.malloc(
    (type == ProjectionType::CLASSIC) ? Nfields * fieldOffset * maxNumVecsProjection :
    Nfields * fieldOffset
    , sizeof(dfloat));

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
  const string oklpath = install_dir + "/okl/elliptic/";
  string filename, kernelName;

  {
    occa::properties properties = platform->kernelInfo;
    properties["defines/p_Nfields"] = Nfields;

    filename = oklpath + "ellipticResidualProjection.okl";
    multiScaledAddwOffsetKernel = platform->device.buildKernel(filename,
                                                                    "multiScaledAddwOffset",
                                                                    properties);
    accumulateKernel = platform->device.buildKernel(filename, "accumulate", properties);
  }
  matvecOperator = [&](occa::memory& o_x, occa::memory & o_Ax)
                   {
                     ellipticOperator(&elliptic, o_x, o_Ax, dfloatString);
                   };
}

void ResidualProjection::pre(occa::memory& o_r)
{
  ++timestep;
  if(timestep < numTimeSteps)
    return;

  if(numVecsProjection <= 0) return;
  computePreProjection(o_r);
}

void ResidualProjection::post(occa::memory& o_x)
{
  if(timestep < numTimeSteps)
    return;
  computePostProjection(o_x);
}