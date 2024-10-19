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
#include "ellipticPrecon.h"
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticSolve(elliptic_t *elliptic,
                   const occa::memory &o_lambda0,
                   const occa::memory &o_lambda1,
                   const occa::memory &o_rhs,
                   occa::memory o_x)
{
  elliptic->o_lambda0 = o_lambda0;
  elliptic->o_lambda1 = o_lambda1;

  auto &options = elliptic->options;
  auto &precon = elliptic->precon;
  auto &mesh = elliptic->mesh;

  int maxIter = 999;
  options.getArgs("MAXIMUM ITERATIONS", maxIter);

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const auto movingMesh = platform->options.compareArgs("MOVING MESH", "TRUE");

  auto printNorm = [&](const occa::memory &o_u, const std::string &txt) {
    const dfloat norm = platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                                            elliptic->Nfields,
                                                            elliptic->fieldOffset,
                                                            elliptic->o_invDegree,
                                                            o_u,
                                                            platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("%s %s norm: %.15e\n", elliptic->name.c_str(), txt.c_str(), norm);
    }
    nekrsCheck(std::isnan(norm),
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s unreasonable %s!\n",
               elliptic->name.c_str(),
               txt.c_str());
  };

  auto o_x0 = platform->deviceMemoryPool.reserve<dfloat>(
      (elliptic->Nfields > 1) ? elliptic->Nfields * elliptic->fieldOffset : mesh->Nlocal);
  nekrsCheck(o_x.size() < o_x0.size(), MPI_COMM_SELF, EXIT_FAILURE, "%s!\n", "unreasonable size of o_x");
  nekrsCheck(o_rhs.size() < o_x.size(), MPI_COMM_SELF, EXIT_FAILURE, "%s!\n", "unreasonable size of o_rhs");

  auto updateResidualWeight = [&]() {
    if (platform->options.compareArgs("LINEAR SOLVER STOPPING CRITERION TYPE", "LEGACY")) {
      if (!elliptic->o_residualWeight.isInitialized()) {
        elliptic->o_residualWeight = platform->device.malloc<dfloat>(mesh->Nlocal);
      }
      elliptic->o_residualWeight.copyFrom(elliptic->o_invDegree);
      platform->linAlg->scale(mesh->Nlocal, 1 / mesh->volume, elliptic->o_residualWeight);
    } else if (platform->options.compareArgs("LINEAR SOLVER STOPPING CRITERION TYPE", "l2_RESIDUAL")) {
      if (!elliptic->o_residualWeight.isInitialized()) {
        elliptic->o_residualWeight = platform->device.malloc<dfloat>(mesh->Nlocal);
      }
      auto Nglobal = mesh->NelementsGlobal * mesh->Np;
      platform->linAlg->axmyz(mesh->Nlocal,
                              1. / Nglobal,
                              mesh->o_invAJw,
                              mesh->o_invAJw,
                              elliptic->o_residualWeight);
      platform->linAlg->axmy(mesh->Nlocal, 1.0, elliptic->o_invDegree, elliptic->o_residualWeight);
    } else if (platform->options.compareArgs("LINEAR SOLVER STOPPING CRITERION TYPE", "L2_RESIDUAL")) {
      elliptic->o_residualWeight = mesh->o_invAJwTimesInvDegree;
    } else {
      const auto txt = platform->options.getArgs("LINEAR SOLVER STOPPING CRITERION TYPE");
      nekrsAbort(MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s <%s>\n",
                 "Invalid LINEAR SOLVER STOPPING CRITERION TYPE",
                 txt.c_str());
    }
  };

  if (!elliptic->o_residualWeight.isInitialized()) {
    updateResidualWeight();
  } else if (movingMesh) {
    updateResidualWeight();
  }

  std::string timerName = elliptic->name;
  if (timerName.find("scalar") != std::string::npos) {
    timerName = "scalar";
  }

  ellipticAllocateWorkspace(elliptic);

  if (options.compareArgs("ELLIPTIC PRECO COEFF FIELD", "TRUE")) {
    if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
      ellipticMultiGridUpdateLambda(elliptic);
    }

    if (options.compareArgs("PRECONDITIONER", "JACOBI") ||
        options.compareArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI")) {
      ellipticUpdateJacobi(elliptic);
    }
  }

  o_x0.copyFrom(o_x);
  if (platform->verbose) {
    printNorm(o_x0, "o_x0");
    printNorm(o_rhs, "o_rhs");
  }

  // compute initial residual r = rhs - Ax0
  auto o_r = [&]() {
    auto o_r = platform->deviceMemoryPool.reserve<dfloat>(o_x0.size());
    auto &o_Ap = o_x;
    ellipticAx(elliptic, mesh->Nelements, mesh->o_elementList, o_x0, o_Ap, dfloatString);
    platform->linAlg
        ->axpbyzMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, -1.0, o_Ap, 1.0, o_rhs, o_r);

    if (elliptic->nullspace) {
      ellipticZeroMean(elliptic, o_r);
    }
    ellipticApplyMask(elliptic, o_r, dfloatString);
    oogs::startFinish(o_r, elliptic->Nfields, elliptic->fieldOffset, ogsDfloat, ogsAdd, elliptic->oogs);
    return o_r;
  }();

  if (platform->verbose) {
    printNorm(o_r, "o_r");
  }

  const auto rdotr = [&]() {
    return platform->linAlg->weightedNorm2Many(mesh->Nlocal,
                                               elliptic->Nfields,
                                               elliptic->fieldOffset,
                                               elliptic->o_residualWeight,
                                               o_r,
                                               platform->comm.mpiComm);
  };

  if (options.compareArgs("INITIAL GUESS", "PROJECTION") ||
      options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ")) {

    platform->timer.tic(timerName + " proj pre", 1);
    elliptic->res00Norm = rdotr();
    nekrsCheck(std::isnan(elliptic->res00Norm),
               MPI_COMM_SELF,
               EXIT_FAILURE,
               "%s unreasonable res00Norm!\n",
               elliptic->name.c_str());

    elliptic->solutionProjection->pre(o_r);

    platform->timer.toc(timerName + " proj pre");
  }

  elliptic->res0Norm = rdotr();
  nekrsCheck(std::isnan(elliptic->res0Norm),
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "%s unreasonable res00Norm!\n",
             elliptic->name.c_str());

  // linear solve
  {
    // absolute tol
    dfloat tol = 1e-6;
    options.getArgs("SOLVER TOLERANCE", tol);

    // absolute tol + relative
    if (!options.getArgs("SOLVER RELATIVE TOLERANCE").empty()) {
      dfloat relTol;
      options.getArgs("SOLVER RELATIVE TOLERANCE", relTol);
      tol = std::max(relTol * elliptic->res0Norm, tol);
    } else { // relative and absolute tolerance are the same
      if (options.compareArgs("LINEAR SOLVER STOPPING CRITERION", "RELATIVE")) {
        tol *= elliptic->res0Norm;
      }
    }

    elliptic->resNorm = elliptic->res0Norm;
    platform->linAlg->fill(o_x.size(), 0.0, o_x);

    if (options.compareArgs("SOLVER", "PCG")) {
      elliptic->Niter = pcg(elliptic, tol, maxIter, elliptic->resNorm, o_r, o_x);
    } else if (options.compareArgs("SOLVER", "PGMRES")) {
      elliptic->Niter = pgmres(elliptic, tol, maxIter, elliptic->resNorm, o_r, o_x);
    } else {
      nekrsAbort(platform->comm.mpiComm,
                 EXIT_FAILURE,
                 "Unknown linear solver %s!\n",
                 options.getArgs("SOLVER").c_str());
    }

    if (elliptic->Niter == maxIter && platform->comm.mpiRank == 0) {
      printf("iteration limit of %s linear solver reached!\n", elliptic->name.c_str());
    }
  }

  if (options.compareArgs("INITIAL GUESS", "PROJECTION") ||
      options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ")) {
    platform->timer.tic(timerName + " proj post", 1);
    elliptic->solutionProjection->post(o_x);
    platform->timer.toc(timerName + " proj post");
  } else {
    elliptic->res00Norm = elliptic->res0Norm;
  }

  platform->linAlg->axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, 1.0, o_x0, 1.0, o_x);

  if (elliptic->nullspace) {
    ellipticZeroMean(elliptic, o_x);
  }

  elliptic->o_lambda0 = nullptr;
  elliptic->o_lambda1 = nullptr;
  ellipticFreeWorkspace(elliptic);
}
