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

void pMGLevelAllocateStorage(pMGLevel *level, int k)
{
  // allocate but reuse finest level
  if (k)
    level->o_x = platform->device.malloc(level->Ncols * sizeof(pfloat));
  if (k)
    level->o_rhs = platform->device.malloc(level->Nrows * sizeof(pfloat));

  level->o_res = platform->device.malloc(level->Ncols * sizeof(pfloat));

  // extra storage for smoother
  const size_t Nbytes = level->Ncols * sizeof(pfloat);

  if (pMGLevel::o_smootherResidual.size() < Nbytes) {
    pMGLevel::o_smootherResidual.free();
    pMGLevel::o_smootherResidual = platform->device.malloc(Nbytes);
  }
  if (pMGLevel::o_smootherResidual2.size() < Nbytes) {
    pMGLevel::o_smootherResidual2.free();
    pMGLevel::o_smootherResidual2 = platform->device.malloc(Nbytes);
  }
  if (pMGLevel::o_smootherUpdate.size() < Nbytes) {
    pMGLevel::o_smootherUpdate.free();
    pMGLevel::o_smootherUpdate = platform->device.malloc(Nbytes);
  }
}

void ellipticMultiGridSetup(elliptic_t *elliptic_, precon_t *precon_)
{
  if (platform->comm.mpiRank == 0)
    printf("building MG preconditioner ... \n");
  fflush(stdout);

  precon_t *precon = precon_;
  // setup new object from fine grid but with constant coeff
  elliptic_t *elliptic = ellipticBuildMultigridLevelFine(elliptic_);
  setupAide options = elliptic_->options;
  mesh_t *mesh = elliptic->mesh;

  // read all the nodes files and load them in a dummy mesh array
  mesh_t **meshLevels = (mesh_t **)calloc(mesh->N + 1, sizeof(mesh_t *));
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
  int *levelDegree = (int *)calloc(numMGLevels, sizeof(int));
  for (int i = 0; i < numMGLevels; ++i)
    levelDegree[i] = elliptic->levels[i];

  int Nmax = levelDegree[0];
  int Nmin = levelDegree[numMGLevels - 1];

  precon->MGSolver = new MGSolver_t(platform->device.occaDevice(), platform->comm.mpiComm, options);
  MGSolver_t::multigridLevel **levels = precon->MGSolver->levels;

  oogs_mode oogsMode = OOGS_AUTO;

  auto autoOverlap = [&](elliptic_t *elliptic) {
    if (!options.compareArgs("MULTIGRID SMOOTHER", "CHEBYSHEV"))
      return;

    auto timeOperator = [&]() {
      const int Nsamples = 10;
      ellipticOperator(elliptic, elliptic->o_p, elliptic->o_Ap, pfloatString);

      platform->device.finish();
      MPI_Barrier(platform->comm.mpiComm);
      const double start = MPI_Wtime();

      for (int test = 0; test < Nsamples; ++test)
        ellipticOperator(elliptic, elliptic->o_p, elliptic->o_Ap, pfloatString);

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
                   elliptic->o_p,
                   elliptic->o_Ap,
                   pfloatString);
      };

      elliptic->oogsAx = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, callback, oogsMode);

      auto overlappedTime = timeOperator();
      if (overlappedTime > nonOverlappedTime)
        elliptic->oogsAx = elliptic->oogs;

      if (platform->comm.mpiRank == 0) {
        printf("autotuning overlap in ellipticOperator: %.2es %.2es ", nonOverlappedTime, overlappedTime);
        if (elliptic->oogsAx != elliptic->oogs)
          printf("(overlap enabled)");

        printf("\n");
      }
    }
  };

  // set up the finest level
  if (Nmax > Nmin) {
    if (platform->comm.mpiRank == 0)
      printf("============= BUILDING pMG%d ==================\n", Nmax);

    elliptic->oogs = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    elliptic->oogsAx = elliptic->oogs;

    levels[0] = new pMGLevel(elliptic, Nmax, options, platform->comm.mpiComm);
    pMGLevelAllocateStorage((pMGLevel *)levels[0], 0);
    precon->MGSolver->numLevels++;

    autoOverlap(elliptic);
  }

  // build intermediate MGLevels
  for (int n = 1; n < numMGLevels - 1; n++) {
    int Nc = levelDegree[n];
    int Nf = levelDegree[n - 1];
    elliptic_t *ellipticFine = ((pMGLevel *)levels[n - 1])->elliptic;
    if (platform->comm.mpiRank == 0)
      printf("============= BUILDING pMG%d ==================\n", Nc);

    elliptic_t *ellipticC = ellipticBuildMultigridLevel(ellipticFine, Nc, Nf);

    ellipticC->oogs = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticC->oogsAx = ellipticC->oogs;

    levels[n] =
        new pMGLevel(elliptic, meshLevels, ellipticFine, ellipticC, Nf, Nc, options, platform->comm.mpiComm);
    pMGLevelAllocateStorage((pMGLevel *)levels[n], n);
    precon->MGSolver->numLevels++;

    autoOverlap(ellipticC);
  }

  // set up coarse level
  elliptic_t *ellipticCoarse;
  if (platform->comm.mpiRank == 0)
    printf("============= BUILDING COARSE pMG%d ==================\n", Nmin);

  if (Nmax > Nmin) {
    int Nc = levelDegree[numMGLevels - 1];
    int Nf = levelDegree[numMGLevels - 2];

    ellipticCoarse = ellipticBuildMultigridLevel(elliptic, Nc, Nf);

    ellipticCoarse->oogs = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticCoarse->oogsAx = ellipticCoarse->oogs;

    elliptic_t *ellipticFine = ((pMGLevel *)levels[numMGLevels - 2])->elliptic;
    levels[numMGLevels - 1] = new pMGLevel(elliptic,
                                           meshLevels,
                                           ellipticFine,
                                           ellipticCoarse,
                                           Nf,
                                           Nc,
                                           options,
                                           platform->comm.mpiComm,
                                           true);

    if (options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE") ||
        options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE"))
      autoOverlap(ellipticCoarse);
  }
  else {
    ellipticCoarse = elliptic;
    levels[numMGLevels - 1] = new pMGLevel(ellipticCoarse, Nmin, options, platform->comm.mpiComm, true);
  }
  pMGLevelAllocateStorage((pMGLevel *)levels[numMGLevels - 1], numMGLevels - 1);
  precon->MGSolver->baseLevel = precon->MGSolver->numLevels;
  precon->MGSolver->numLevels++;

  if (options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")) {
    if (options.compareArgs("MULTIGRID SEMFEM", "TRUE")) {

      precon->SEMFEMSolver = new SEMFEMSolver_t(ellipticCoarse);
      if (options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE")) {
        auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];
        auto &o_tmp = pMGLevel::o_smootherUpdate;
        precon->MGSolver->coarseLevel->solvePtr =
            [elliptic, baseLevel, &o_tmp](MGSolver_t::coarseLevel_t *coarseLevel,
                                          occa::memory &o_rhs,
                                          occa::memory &o_x) {
              occa::memory o_res = baseLevel->o_res;
              baseLevel->smooth(o_rhs, o_x, true);
              baseLevel->residual(o_rhs, o_x, o_res);

              auto precon = elliptic->precon;
              precon->SEMFEMSolver->run(o_res, o_tmp);

              platform->linAlg->paxpby(baseLevel->Nrows, 1.0, o_tmp, 1.0, o_x);
              baseLevel->smooth(o_rhs, o_x, false);
            };
      }
      else {
        precon->MGSolver->coarseLevel->solvePtr =
            [elliptic](MGSolver_t::coarseLevel_t *, occa::memory &o_rhs, occa::memory &o_x) {
              auto precon = elliptic->precon;
              precon->SEMFEMSolver->run(o_rhs, o_x);
            };
      }
    }
    else {

      hlong *coarseGlobalStarts = (hlong *)calloc(platform->comm.mpiCommSize + 1, sizeof(hlong));

      nonZero_t *coarseA;
      dlong nnzCoarseA;

      if (options.compareArgs("GALERKIN COARSE OPERATOR", "TRUE"))
        ellipticBuildFEMGalerkinHex3D(ellipticCoarse, elliptic, &coarseA, &nnzCoarseA, coarseGlobalStarts);
      else
        ellipticBuildFEM(ellipticCoarse, &coarseA, &nnzCoarseA, coarseGlobalStarts);

      hlong *Rows = (hlong *)calloc(nnzCoarseA, sizeof(hlong));
      hlong *Cols = (hlong *)calloc(nnzCoarseA, sizeof(hlong));
      dfloat *Vals = (dfloat *)calloc(nnzCoarseA, sizeof(dfloat));

      for (dlong i = 0; i < nnzCoarseA; i++) {
        Rows[i] = coarseA[i].row;
        Cols[i] = coarseA[i].col;
        Vals[i] = coarseA[i].val;
      }
      free(coarseA);

      precon->MGSolver->coarseLevel
          ->setupSolver(coarseGlobalStarts, nnzCoarseA, Rows, Cols, Vals, elliptic->allNeumann);

      free(coarseGlobalStarts);
      free(Rows);
      free(Cols);
      free(Vals);

      MGSolver_t::coarseLevel_t *coarseLevel = precon->MGSolver->coarseLevel;
      coarseLevel->ogs = ellipticCoarse->ogs;
      coarseLevel->o_weight = ellipticCoarse->o_invDegree;
      coarseLevel->weight = (pfloat *)calloc(ellipticCoarse->mesh->Nlocal, sizeof(pfloat));
      coarseLevel->o_weight.copyTo(coarseLevel->weight, ellipticCoarse->mesh->Nlocal * sizeof(pfloat));
      coarseLevel->h_Gx = platform->device.mallocHost(coarseLevel->ogs->Ngather * sizeof(pfloat));
      coarseLevel->Gx = (pfloat *)coarseLevel->h_Gx.ptr();
      coarseLevel->o_Gx = platform->device.malloc(coarseLevel->ogs->Ngather * sizeof(pfloat));
      coarseLevel->h_Sx = platform->device.mallocHost(ellipticCoarse->mesh->Nlocal * sizeof(pfloat));
      coarseLevel->Sx = (pfloat *)coarseLevel->h_Sx.ptr();
      coarseLevel->o_Sx = platform->device.malloc(ellipticCoarse->mesh->Nlocal * sizeof(pfloat));

      if (options.compareArgs("MULTIGRID COARSE SOLVE AND SMOOTH", "TRUE")) {
        auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];
        auto &o_tmp = pMGLevel::o_smootherUpdate;
        precon->MGSolver->coarseLevel->solvePtr = [baseLevel, &o_tmp](MGSolver_t::coarseLevel_t *coarseLevel,
                                                                      occa::memory &o_rhs,
                                                                      occa::memory &o_x) {
          occa::memory o_res = baseLevel->o_res;
          baseLevel->smooth(o_rhs, o_x, true);
          baseLevel->residual(o_rhs, o_x, o_res);

          coarseLevel->solve(o_res, o_tmp);

          platform->linAlg->paxpby(baseLevel->Nrows, 1.0, o_tmp, 1.0, o_x);
          baseLevel->smooth(o_rhs, o_x, false);
        };
      }
    }
  }
  else {
    auto baseLevel = (pMGLevel *)levels[numMGLevels - 1];
    precon->MGSolver->coarseLevel->solvePtr =
        [baseLevel](MGSolver_t::coarseLevel_t *, occa::memory &o_rhs, occa::memory &o_x) {
          baseLevel->smooth(o_rhs, o_x, true);
        };
  }

  for (int n = 1; n < mesh->N + 1; n++)
    delete meshLevels[n];
  free(meshLevels);

  if (platform->comm.mpiRank == 0) {
    printf("-----------------------------------------------------------------------\n");
    printf("level|    Type    |                 |     Smoother                    |\n");
    printf("     |            |                 |                                 |\n");
    printf("-----------------------------------------------------------------------\n");
  }

  for (int lev = 0; lev < precon->MGSolver->numLevels; lev++) {
    if (platform->comm.mpiRank == 0) {
      printf(" %3d ", lev);
      fflush(stdout);
    }
    levels[lev]->Report();
  }

  if (platform->comm.mpiRank == 0)
    printf("-----------------------------------------------------------------------\n");
}
