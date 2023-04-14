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
#include "ellipticSolutionProjection.h"
#include <iostream>
#include "timer.hpp"
#include "platform.hpp"
#include "linAlg.hpp"

void SolutionProjection::matvec(occa::memory &o_Ax,
                                const dlong Ax_offset,
                                occa::memory &o_x,
                                const dlong x_offset)
{
  occa::memory o_xtmp = o_x + (Nfields * x_offset * sizeof(dfloat)) * fieldOffset;
  occa::memory o_Axtmp = o_Ax + (Nfields * Ax_offset * sizeof(dfloat)) * fieldOffset;
  matvecOperator(o_xtmp, o_Axtmp);
}

void SolutionProjection::updateProjectionSpace()
{

  if (numVecsProjection <= 0)
    return;

  double flopCount = 0.0;

#if USE_WEIGHTED_INNER_PROD_MULTI_DEVICE
  platform->linAlg->weightedInnerProdMulti(
      Nlocal,
      numVecsProjection,
      Nfields,
      fieldOffset,
      o_invDegree,
      o_xx,
      o_bb,
      platform->comm.mpiComm,
      o_alpha,
      (type == ProjectionType::CLASSIC) ? Nfields * (numVecsProjection - 1) * fieldOffset : 0);
  o_alpha.copyTo(alpha, sizeof(dfloat) * numVecsProjection);
#else
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
      (type == ProjectionType::CLASSIC) ? Nfields * (numVecsProjection - 1) * fieldOffset : 0);
  o_alpha.copyFrom(alpha, sizeof(dfloat) * numVecsProjection);
#endif

  const dfloat norm_orig = alpha[numVecsProjection - 1];
  const dfloat one = 1.0;
  multiScaledAddwOffsetKernel(Nlocal,
                              numVecsProjection,
                              Nfields * (numVecsProjection - 1) * fieldOffset,
                              fieldOffset,
                              o_alpha,
                              one,
                              o_xx);
  if (type == ProjectionType::CLASSIC)
    multiScaledAddwOffsetKernel(Nlocal,
                                numVecsProjection,
                                Nfields * (numVecsProjection - 1) * fieldOffset,
                                fieldOffset,
                                o_alpha,
                                one,
                                o_bb);

  flopCount += 3 * static_cast<double>(Nlocal) * Nfields * (numVecsProjection - 1);
  flopCount *= (type == ProjectionType::CLASSIC) ? 2 : 1;

  dfloat sumAlpha = 0;
  for (int k = 0; k < numVecsProjection - 1; ++k)
    sumAlpha += alpha[k] * alpha[k];

  dfloat norm_new = norm_orig - sumAlpha;
  // printf("norm_new:%g norm_orig:%g sumAlpha:%g\n", norm_new, norm_orig, sumAlpha);
  norm_new = sqrt(norm_new);

  dfloat tol = 1e-7;
  const dfloat test = norm_new / norm_orig;
  if (test > tol) {
    const dfloat scale = 1.0 / norm_new;
    platform->linAlg->scaleMany(Nlocal,
                                Nfields,
                                fieldOffset,
                                scale,
                                o_xx,
                                fieldOffset * Nfields * (numVecsProjection - 1));
    if (type == ProjectionType::CLASSIC)
      platform->linAlg->scaleMany(Nlocal,
                                  Nfields,
                                  fieldOffset,
                                  scale,
                                  o_bb,
                                  fieldOffset * Nfields * (numVecsProjection - 1));
    flopCount += static_cast<double>(Nlocal) * Nfields;
    flopCount *= (type == ProjectionType::CLASSIC) ? 2 : 1;
  }
  else {
    if (platform->comm.mpiRank == 0) {
      std::cout << "solutionProjection " << solverName
                << ": Discard new solution as it is linearly dependent!\n";
    }
    numVecsProjection--;
  }

  platform->flopCounter->add(solverName + " SolutionProjection::updateProjectionSpace", flopCount);
}

void SolutionProjection::computePreProjection(occa::memory &o_r)
{

  dfloat flopCount = 0.0;

  dfloat one = 1.0;
  dfloat zero = 0.0;
  dfloat mone = -1.0;
  if (numVecsProjection <= 0)
    return;

#if USE_WEIGHTED_INNER_PROD_MULTI_DEVICE
  platform->linAlg->weightedInnerProdMulti(Nlocal,
                                           numVecsProjection,
                                           Nfields,
                                           fieldOffset,
                                           o_invDegree,
                                           o_xx,
                                           o_r,
                                           platform->comm.mpiComm,
                                           o_alpha,
                                           Nfields * 0 * fieldOffset);
#else
  platform->linAlg->weightedInnerProdMulti(Nlocal,
                                           numVecsProjection,
                                           Nfields,
                                           fieldOffset,
                                           o_invDegree,
                                           o_xx,
                                           o_r,
                                           platform->comm.mpiComm,
                                           alpha,
                                           Nfields * 0 * fieldOffset);
  o_alpha.copyFrom(alpha, sizeof(dfloat) * numVecsProjection);
#endif

  // o_xbar = sum_i alpha_i * o_xx_i
  accumulateKernel(Nlocal, numVecsProjection, fieldOffset, o_alpha, o_xx, o_xbar);

  flopCount += Nfields * (1 + 2 * (numVecsProjection - 1)) * static_cast<double>(Nlocal);
  if (type == ProjectionType::CLASSIC) {
    accumulateKernel(Nlocal, numVecsProjection, fieldOffset, o_alpha, o_bb, o_rtmp);
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, mone, o_rtmp, one, o_r);

    flopCount += Nfields * (1 + 2 * (numVecsProjection - 1)) * static_cast<double>(Nlocal); // accumulation
  }
  else if (type == ProjectionType::ACONJ) {
    matvec(o_bb, 0, o_xbar, 0);
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, mone, o_bb, one, o_r);
  }

  platform->flopCounter->add(solverName + " SolutionProjection::computePreProjection", flopCount);
}

void SolutionProjection::computePostProjection(occa::memory &o_x)
{
  const dfloat one = 1.0;
  const dfloat zero = 0.0;

  if (numVecsProjection == 0) {
    // reset bases
    numVecsProjection = 1;
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat));
  }
  else if (numVecsProjection == maxNumVecsProjection) {
    numVecsProjection = 1;
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, one, o_xbar, one, o_x);
    o_xx.copyFrom(o_x, Nfields * fieldOffset * sizeof(dfloat));
  }
  else {
    numVecsProjection++;
    // xx[m-1] = x
    o_xx.copyFrom(o_x,
                  fieldOffset *  Nfields * sizeof(dfloat),
                  fieldOffset * (Nfields * sizeof(dfloat) * (numVecsProjection - 1)),
                  0);
    // x = x + xbar
    platform->linAlg->axpbyMany(Nlocal, Nfields, fieldOffset, one, o_xbar, one, o_x);
  }
  const dlong previousNumVecsProjection = numVecsProjection;
  const dlong bOffset = (type == ProjectionType::CLASSIC) ? numVecsProjection - 1 : 0;
  matvec(o_bb, bOffset, o_xx, numVecsProjection - 1);

  updateProjectionSpace();
  if (numVecsProjection < previousNumVecsProjection) { // Last vector was linearly dependent, reset space
    numVecsProjection = 1;
    o_xx.copyFrom(o_x,
                  Nfields * fieldOffset * sizeof(dfloat)); // writes first n words of o_xx, first approximation vector
    matvec(o_bb, 0, o_xx, 0);
    updateProjectionSpace();
  }
}

SolutionProjection::SolutionProjection(elliptic_t &elliptic,
                                       const ProjectionType _type,
                                       const dlong _maxNumVecsProjection,
                                       const dlong _numTimeSteps)
    : maxNumVecsProjection(_maxNumVecsProjection), numTimeSteps(_numTimeSteps), type(_type),
      alpha((dfloat *)calloc(maxNumVecsProjection, sizeof(dfloat))), numVecsProjection(0),
      prevNumVecsProjection(0), Nlocal(elliptic.mesh->Np * elliptic.mesh->Nelements),
      fieldOffset(elliptic.fieldOffset), Nfields(elliptic.Nfields), timestep(0),
      verbose(platform->options.compareArgs("VERBOSE", "TRUE")), o_invDegree(elliptic.mesh->ogs->o_invDegree),
      o_rtmp(elliptic.o_z), o_Ap(elliptic.o_Ap)
{
  solverName = elliptic.name;

  platform_t *platform = platform_t::getInstance();

  o_alpha = platform->device.malloc(maxNumVecsProjection * sizeof(dfloat));
  o_xbar = platform->device.malloc((Nfields * sizeof(dfloat)) * fieldOffset);
  o_xx = platform->device.malloc((Nfields * maxNumVecsProjection * sizeof(dfloat)) * fieldOffset);
  o_bb =
      platform->device.malloc((type == ProjectionType::CLASSIC) ? Nfields * fieldOffset * maxNumVecsProjection
                                                                : Nfields * fieldOffset,
                              sizeof(dfloat));

  const std::string sectionIdentifier = std::to_string(Nfields) + "-";

  {
    multiScaledAddwOffsetKernel = platform->kernels.get(sectionIdentifier + "multiScaledAddwOffset");
    accumulateKernel = platform->kernels.get(sectionIdentifier + "accumulate");
  }

  matvecOperator = [&](occa::memory &o_x, occa::memory &o_Ax) {
    ellipticOperator(&elliptic, o_x, o_Ax, dfloatString);
  };

  maskOperator = [&](occa::memory &o_x) { ellipticApplyMask(&elliptic, o_x, dfloatString); };
}

void SolutionProjection::pre(occa::memory &o_r)
{
  ++timestep;
  if (timestep < numTimeSteps)
    return;

  if (numVecsProjection <= 0)
    return;

  prevNumVecsProjection = numVecsProjection;
  computePreProjection(o_r);
}

void SolutionProjection::post(occa::memory &o_x)
{
  if (timestep < numTimeSteps)
    return;
  computePostProjection(o_x);
}
