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

#include "platform.hpp"
#include "elliptic.h"
#include "ellipticPrecon.h"
#include "ellipticMultiGrid.h"
#include "ellipticBuildFEM.hpp"

void ellipticMultiGridSetup(elliptic_t *elliptic_)
{
  if (platform->comm.mpiRank == 0) {
    printf("building MG preconditioner ... \n");
  }
  fflush(stdout);

  elliptic_->precon = new precon_t();
  const auto precon = elliptic_->precon;

  // setup new object from fine grid but with constant coeff
  elliptic_t *elliptic = ellipticBuildMultigridLevelFine(elliptic_);
  setupAide options = elliptic_->options;
  mesh_t *mesh = elliptic->mesh;

  // read all the nodes files and load them in a dummy mesh array
  std::vector<mesh_t *> meshLevels(mesh->N + 1);
  for (int n = 1; n < mesh->N + 1; n++) {
    meshLevels[n] = new mesh_t();
    meshLevels[n]->Nverts = mesh->Nverts;
    meshLevels[n]->Nfaces = mesh->Nfaces;
    meshLevels[n]->Nfields = mesh->Nfields; // TW: ahem

    switch (elliptic->elementType) {
    case HEXAHEDRA:
      meshLoadReferenceNodesHex3D(meshLevels[n], n, 1);
      break;
    }
  }

  // set the number of MG levels and their degree
  int numMGLevels = elliptic->nLevels;
  std::vector<int> levelDegree(numMGLevels);
  for (int i = 0; i < numMGLevels; ++i) {
    levelDegree[i] = elliptic->levels[i];
  }

  int Nmax = levelDegree[0];
  int Nmin = levelDegree[numMGLevels - 1];

  precon->MGSolver = new MGSolver_t(platform->device.occaDevice(), platform->comm.mpiComm, options);
  MGSolver_t::multigridLevel **levels = precon->MGSolver->levels;

  oogs_mode oogsMode = OOGS_AUTO;

  auto autoOverlap = [&](elliptic_t *elliptic) {
    if (!options.compareArgs("MULTIGRID SMOOTHER", "CHEBYSHEV")) {
      return;
    }

    auto o_p = platform->deviceMemoryPool.reserve<pfloat>(mesh->Nlocal);
    auto o_Ap = platform->deviceMemoryPool.reserve<pfloat>(mesh->Nlocal);

    auto timeOperator = [&]() {
      const int Nsamples = 10;
      ellipticOperator(elliptic, o_p, o_Ap, pfloatString);

      platform->device.finish();
      MPI_Barrier(platform->comm.mpiComm);
      const double start = MPI_Wtime();

      for (int test = 0; test < Nsamples; ++test) {
        ellipticOperator(elliptic, o_p, o_Ap, pfloatString);
      }

      platform->device.finish();
      double elapsed = (MPI_Wtime() - start) / Nsamples;
      MPI_Allreduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);

      return elapsed;
    };

    if (platform->options.compareArgs("ENABLE GS COMM OVERLAP", "TRUE")) {
      auto nonOverlappedTime = timeOperator();
      auto callback = [&]() {
        ellipticAx(elliptic,
                   elliptic->mesh->NlocalGatherElements,
                   elliptic->mesh->o_localGatherElementList,
                   o_p,
                   o_Ap,
                   pfloatString);
      };

      elliptic->oogsAx = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, callback, oogsMode);

      auto overlappedTime = timeOperator();
      if (overlappedTime > nonOverlappedTime) {
        elliptic->oogsAx = elliptic->oogs;
      }

      if (platform->comm.mpiRank == 0) {
        printf("testing overlap in ellipticOperator: %.2es %.2es ", nonOverlappedTime, overlappedTime);
        if (elliptic->oogsAx != elliptic->oogs) {
          printf("(overlap enabled)");
        }

        printf("\n");
      }
    }
  };

  // set up the finest level 0
  if (Nmax > Nmin) {
    if (platform->comm.mpiRank == 0) {
      printf("============= BUILDING pMG%d ==================\n", Nmax);
    }

    elliptic->oogs = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    elliptic->oogsAx = elliptic->oogs;

    levels[0] = new pMGLevel(elliptic, Nmax, options, platform->comm.mpiComm);
    precon->MGSolver->numLevels++;

    autoOverlap(elliptic);
  }

  // build intermediate MGLevels
  for (int n = 1; n < numMGLevels - 1; n++) {
    int Nc = levelDegree[n];
    int Nf = levelDegree[n - 1];
    elliptic_t *ellipticFine = ((pMGLevel *)levels[n - 1])->elliptic;
    if (platform->comm.mpiRank == 0) {
      printf("============= BUILDING pMG%d ==================\n", Nc);
    }

    elliptic_t *ellipticC = ellipticBuildMultigridLevel(ellipticFine, Nc, Nf);

    ellipticC->oogs = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticC->oogsAx = ellipticC->oogs;

    levels[n] = new pMGLevel(elliptic,
                             meshLevels.data(),
                             ellipticFine,
                             ellipticC,
                             Nf,
                             Nc,
                             options,
                             platform->comm.mpiComm);
    precon->MGSolver->numLevels++;

    autoOverlap(ellipticC);
  }

  // set up coarse level numMGLevels - 1
  elliptic_t *ellipticCoarse;
  if (platform->comm.mpiRank == 0) {
    printf("============= BUILDING COARSE pMG%d ==================\n", Nmin);
  }

  if (Nmax > Nmin) {
    int Nc = levelDegree[numMGLevels - 1];
    int Nf = levelDegree[numMGLevels - 2];
    elliptic_t *ellipticFine = ((pMGLevel *)levels[numMGLevels - 2])->elliptic;

    ellipticCoarse = ellipticBuildMultigridLevel(ellipticFine, Nc, Nf);

    ellipticCoarse->oogs = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticCoarse->oogsAx = ellipticCoarse->oogs;

    levels[numMGLevels - 1] = new pMGLevel(elliptic,
                                           meshLevels.data(),
                                           ellipticFine,
                                           ellipticCoarse,
                                           Nf,
                                           Nc,
                                           options,
                                           platform->comm.mpiComm,
                                           true);

    if (options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE") ||
        options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE")) {
      autoOverlap(ellipticCoarse);
    }
  } else {
    ellipticCoarse = elliptic;
    levels[numMGLevels - 1] = new pMGLevel(ellipticCoarse, Nmin, options, platform->comm.mpiComm, true);
  }
  precon->MGSolver->baseLevel = precon->MGSolver->numLevels;
  precon->MGSolver->numLevels++;

  if (options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")) {
    if (options.compareArgs("MULTIGRID SEMFEM", "TRUE")) {
      precon->SEMFEMSolver = new SEMFEMSolver_t(ellipticCoarse);
      if (options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE")) {
        auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];

        precon->MGSolver->coarseLevel->solvePtr =
            [elliptic,
             baseLevel](MGSolver_t::coarseLevel_t *coarseLevel, occa::memory &o_rhs, occa::memory &o_x) {
              auto &o_res = baseLevel->o_res;
              baseLevel->smooth(o_rhs, o_x, true);
              baseLevel->residual(o_rhs, o_x, o_res);

              auto o_tmp = platform->deviceMemoryPool.reserve<pfloat>(o_x.size());
              elliptic->precon->SEMFEMSolver->run(o_res, o_tmp);

              platform->linAlg->paxpby(o_x.size(), 1.0, o_tmp, 1.0, o_x);
              baseLevel->smooth(o_rhs, o_x, false);
            };
      } else {
        precon->MGSolver->coarseLevel->solvePtr =
            [elliptic](MGSolver_t::coarseLevel_t *, occa::memory &o_rhs, occa::memory &o_x) {
              elliptic->precon->SEMFEMSolver->run(o_rhs, o_x);
            };
      }
    } else {

      hlong *coarseGlobalStarts = (hlong *)calloc(platform->comm.mpiCommSize + 1, sizeof(hlong));

      nonZero_t *coarseA;
      dlong nnzCoarseA;

      if (options.compareArgs("GALERKIN COARSE OPERATOR", "TRUE")) {
        ellipticBuildFEMGalerkinHex3D(ellipticCoarse, elliptic, &coarseA, &nnzCoarseA, coarseGlobalStarts);
      } else {
        ellipticBuildFEM(ellipticCoarse, &coarseA, &nnzCoarseA, coarseGlobalStarts);
      }

      hlong *Rows = (hlong *)calloc(nnzCoarseA, sizeof(hlong));
      hlong *Cols = (hlong *)calloc(nnzCoarseA, sizeof(hlong));
      dfloat *Vals = (dfloat *)calloc(nnzCoarseA, sizeof(dfloat));

      for (dlong i = 0; i < nnzCoarseA; i++) {
        Rows[i] = coarseA[i].row;
        Cols[i] = coarseA[i].col;
        Vals[i] = coarseA[i].val;

        nekrsCheck(Rows[i] < 0 || Cols[i] < 0 || std::isnan(Vals[i]),
                   MPI_COMM_SELF,
                   EXIT_FAILURE,
                   "invalid {row %lld, col %lld , val %g}\n",
                   Rows[i],
                   Cols[i],
                   Vals[i]);
      }
      free(coarseA);

      precon->MGSolver->coarseLevel
          ->setupSolver(coarseGlobalStarts, nnzCoarseA, Rows, Cols, Vals, elliptic->nullspace);

      free(coarseGlobalStarts);
      free(Rows);
      free(Cols);
      free(Vals);

      MGSolver_t::coarseLevel_t *coarseLevel = precon->MGSolver->coarseLevel;
      coarseLevel->ogs = ellipticCoarse->ogs;

      coarseLevel->o_weight = ellipticCoarse->o_invDegree;
      coarseLevel->weight = (pfloat *)calloc(ellipticCoarse->mesh->Nlocal, sizeof(pfloat));
      coarseLevel->o_weight.copyTo(coarseLevel->weight, ellipticCoarse->mesh->Nlocal);

      coarseLevel->h_Gx = platform->device.mallocHost<pfloat>(coarseLevel->ogs->Ngather);
      coarseLevel->Gx = (pfloat *)coarseLevel->h_Gx.ptr();
      coarseLevel->o_Gx = platform->device.malloc<pfloat>(coarseLevel->ogs->Ngather);

      coarseLevel->h_Sx = platform->device.mallocHost<pfloat>(ellipticCoarse->mesh->Nlocal);
      coarseLevel->Sx = (pfloat *)coarseLevel->h_Sx.ptr();
      coarseLevel->o_Sx = platform->device.malloc<pfloat>(ellipticCoarse->mesh->Nlocal);

      if (options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE")) {
        auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];

        precon->MGSolver->coarseLevel->solvePtr =
            [baseLevel](MGSolver_t::coarseLevel_t *coarseLevel, occa::memory &o_rhs, occa::memory &o_x) {
              occa::memory o_res = baseLevel->o_res;
              baseLevel->smooth(o_rhs, o_x, true);
              baseLevel->residual(o_rhs, o_x, o_res);

              auto o_tmp = platform->deviceMemoryPool.reserve<pfloat>(baseLevel->Nrows);
              coarseLevel->solve(o_res, o_tmp);

              platform->linAlg->paxpby(baseLevel->Nrows, 1.0, o_tmp, 1.0, o_x);
              baseLevel->smooth(o_rhs, o_x, false);
            };
      }
    }
  } else {
    auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];
    precon->MGSolver->coarseLevel->solvePtr =
        [baseLevel](MGSolver_t::coarseLevel_t *, occa::memory &o_rhs, occa::memory &o_x) {
          baseLevel->smooth(o_rhs, o_x, true);
        };
  }

  if (platform->comm.mpiRank == 0) {
    printf("-----------------------------------------------------------------------\n");
    printf("level|    Type    |                 |     Smoother                    |\n");
    printf("     |            |                 |                                 |\n");
    printf("-----------------------------------------------------------------------\n");
  }

  for (int lev = 0; lev < precon->MGSolver->numLevels; lev++) {
    if (platform->comm.mpiRank == 0) {
      printf(" %3d ", lev);
    }
    levels[lev]->Report();
  }

  if (platform->comm.mpiRank == 0) {
    printf("-----------------------------------------------------------------------\n");
  }

  fflush(stdout);
}
