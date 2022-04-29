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
#include "platform.hpp"

void ellipticMultiGridSetup(elliptic_t* elliptic_, precon_t* precon)
{
  // setup new object from fine grid but with constant coeff
  elliptic_t* elliptic = ellipticBuildMultigridLevelFine(elliptic_);
  setupAide options = elliptic_->options;
  mesh_t* mesh = elliptic->mesh;

  //read all the nodes files and load them in a dummy mesh array
  mesh_t** meshLevels = (mesh_t**) calloc(mesh->N + 1,sizeof(mesh_t*));
  for (int n = 1; n < mesh->N + 1; n++) {
    meshLevels[n] = new mesh_t();
    meshLevels[n]->Nverts = mesh->Nverts;
    meshLevels[n]->Nfaces = mesh->Nfaces;
    meshLevels[n]->Nfields = mesh->Nfields; // TW: ahem

    switch(elliptic->elementType) {
    case HEXAHEDRA:
      meshLoadReferenceNodesHex3D(meshLevels[n], n, 1);
      break;
    }
  }

  //set the number of MG levels and their degree
  int numMGLevels = elliptic->nLevels;
  int* levelDegree = (int*) calloc(numMGLevels,sizeof(int));
  for(int i = 0; i < numMGLevels; ++i)
    levelDegree[i] = elliptic->levels[i];

  int Nmax = levelDegree[0];
  int Nmin = levelDegree[numMGLevels - 1];

  //initialize parAlmond
  precon->parAlmond = parAlmond::Init(platform->device.occaDevice(), platform->comm.mpiComm, options);
  parAlmond::multigridLevel** levels = precon->parAlmond->levels;

  oogs_mode oogsMode = OOGS_AUTO;
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;

  //set up the finest level
  if (Nmax > Nmin) {
    if(platform->comm.mpiRank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nmax);

    elliptic->o_lambdaPfloat = platform->device.malloc(2 * mesh->Nelements * mesh->Np, sizeof(pfloat));
    elliptic->copyDfloatToPfloatKernel(2 * mesh->Nelements * mesh->Np,
      elliptic->o_lambda,
      elliptic->o_lambdaPfloat);

    auto callback = [&]()
                    {
                      ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
                                 elliptic->o_p, elliptic->o_Ap, pfloatString);
                    };
    elliptic->oogs   = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    elliptic->oogsAx = elliptic->oogs;
    if(options.compareArgs("GS OVERLAP", "TRUE"))
      elliptic->oogsAx = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, callback, oogsMode);

    levels[0] = new MGLevel(elliptic, Nmax, options,
                            precon->parAlmond->ktype, platform->comm.mpiComm);
    MGLevelAllocateStorage((MGLevel*) levels[0], 0,
                           precon->parAlmond->ctype);
    precon->parAlmond->numLevels++;
  }

  //build a MGLevel for every degree (except degree 1)
  for (int n = 1; n < numMGLevels - 1; n++) {
    int Nc = levelDegree[n];
    int Nf = levelDegree[n - 1];
    elliptic_t* ellipticFine = ((MGLevel*) levels[n - 1])->elliptic;
    //build elliptic struct for this degree
    if(platform->comm.mpiRank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nc);

    elliptic_t* ellipticC = ellipticBuildMultigridLevel(ellipticFine,Nc,Nf);

    auto callback = [&]()
                    {
                      ellipticAx(ellipticC,
                                 ellipticC->mesh->NlocalGatherElements,
                                 ellipticC->mesh->o_localGatherElementList,
                                 ellipticC->o_p,
                                 ellipticC->o_Ap,
                                 pfloatString);
                    };
    ellipticC->oogs   = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticC->oogsAx = ellipticC->oogs; 
    if(options.compareArgs("GS OVERLAP", "TRUE"))
      ellipticC->oogsAx = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, callback, oogsMode);

    //add the level manually
    levels[n] = new MGLevel(elliptic,
                            meshLevels,
                            ellipticFine,
                            ellipticC,
                            Nf, Nc,
                            options,
                            precon->parAlmond->ktype, platform->comm.mpiComm);
    MGLevelAllocateStorage((MGLevel*) levels[n], n,
                           precon->parAlmond->ctype);
    precon->parAlmond->numLevels++;
  }

  /* build degree 1 problem and pass to AMG */
  nonZero_t* coarseA;
  dlong nnzCoarseA;
  ogs_t* coarseogs;

  //set up the base level
  elliptic_t* ellipticCoarse;
  if (Nmax > Nmin) {
    int Nc = levelDegree[numMGLevels - 1];
    int Nf = levelDegree[numMGLevels - 2];

    if(platform->comm.mpiRank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nmin);

    ellipticCoarse = ellipticBuildMultigridLevel(elliptic,Nc,Nf);

    auto callback = [&]()
                    {
                      ellipticAx(ellipticCoarse,
                                 ellipticCoarse->mesh->NlocalGatherElements,
                                 ellipticCoarse->mesh->o_localGatherElementList,
                                 ellipticCoarse->o_p,
                                 ellipticCoarse->o_Ap,
                                 pfloatString);
                    };
    ellipticCoarse->oogs   = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticCoarse->oogsAx = ellipticCoarse->oogs;
    if(options.compareArgs("GS OVERLAP", "TRUE") && options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
      ellipticCoarse->oogsAx = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, callback, oogsMode);
  } else {
    ellipticCoarse = elliptic;
  }

  if(options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")){
    if(options.compareArgs("MULTIGRID COARSE SEMFEM", "TRUE")){
      ellipticSEMFEMSetup(ellipticCoarse);
      precon->parAlmond->coarseLevel = new parAlmond::coarseSolver(precon->parAlmond->options, platform->comm.mpiComm);
      precon->parAlmond->coarseLevel->useSEMFEM = true;
      precon->parAlmond->coarseLevel->semfemSolver = [ellipticCoarse](occa::memory o_rhs, occa::memory o_x)
      {
        ellipticSEMFEMSolve(ellipticCoarse, o_rhs, o_x);
      };
      precon->parAlmond->baseLevel = precon->parAlmond->numLevels;
      precon->parAlmond->numLevels++;
    } else {
      int basisNp = ellipticCoarse->mesh->Np;
      dfloat* basis = NULL;

      hlong* coarseGlobalStarts = (hlong*) calloc(platform->comm.mpiCommSize + 1, sizeof(hlong));

      if(options.compareArgs("GALERKIN COARSE OPERATOR","TRUE"))
        ellipticBuildContinuousGalerkinHex3D(ellipticCoarse,elliptic,&coarseA,&nnzCoarseA,
                                             &coarseogs,coarseGlobalStarts);
      else
        ellipticBuildContinuous(ellipticCoarse, &coarseA, &nnzCoarseA,&coarseogs,
                                coarseGlobalStarts);

      hlong* Rows = (hlong*) calloc(nnzCoarseA, sizeof(hlong));
      hlong* Cols = (hlong*) calloc(nnzCoarseA, sizeof(hlong));
      dfloat* Vals = (dfloat*) calloc(nnzCoarseA,sizeof(dfloat));

      for (dlong i = 0; i < nnzCoarseA; i++) {
        Rows[i] = coarseA[i].row;
        Cols[i] = coarseA[i].col;
        Vals[i] = coarseA[i].val;
      }
      free(coarseA);

      // build amg starting at level N=1
      parAlmond::AMGSetup(precon->parAlmond,
                          coarseGlobalStarts,
                          nnzCoarseA,
                          Rows,
                          Cols,
                          Vals,
                          elliptic->allNeumann,
                          0.0);

      free(Rows);
      free(Cols);
      free(Vals);
    }
  } else {
    precon->parAlmond->baseLevel = precon->parAlmond->numLevels;
    precon->parAlmond->numLevels++;
  }

  //overwrite the finest AMG level with the degree 1 matrix free level
  // delete levels[numMGLevels-1];
  if (Nmax > Nmin) {
    int Nc = levelDegree[numMGLevels - 1];
    int Nf = levelDegree[numMGLevels - 2];
    elliptic_t* ellipticFine = ((MGLevel*) levels[numMGLevels - 2])->elliptic;
    levels[numMGLevels - 1] = new MGLevel(elliptic,
                                          meshLevels,
                                          ellipticFine,
                                          ellipticCoarse,
                                          Nf, Nc,
                                          options,
                                          precon->parAlmond->ktype, platform->comm.mpiComm,
                                          true);
  } else {
    levels[numMGLevels - 1] = new MGLevel(ellipticCoarse, Nmin, options,
                                          precon->parAlmond->ktype, platform->comm.mpiComm, true);
  }
  MGLevelAllocateStorage((MGLevel*) levels[numMGLevels - 1], numMGLevels - 1,
                         precon->parAlmond->ctype);

  if(options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")){
    if(options.compareArgs("MULTIGRID COARSE SEMFEM", "TRUE")){
    } else {
      //tell parAlmond to gather when going to the next level
      if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
        if (precon->parAlmond->numLevels > numMGLevels) {
          parAlmond::agmgLevel* nextLevel
            = (parAlmond::agmgLevel*)precon->parAlmond->levels[numMGLevels];

          nextLevel->gatherLevel = true;
          nextLevel->ogs = ellipticCoarse->ogs;
          nextLevel->Gx = (dfloat*) calloc(levels[numMGLevels - 1]->Ncols,sizeof(dfloat));
          nextLevel->Sx = (dfloat*) calloc(ellipticCoarse->mesh->Np * ellipticCoarse->mesh->Nelements,
                                           sizeof(dfloat));
          nextLevel->o_Gx =
            platform->device.malloc(levels[numMGLevels - 1]->Ncols * sizeof(dfloat),
                                                nextLevel->Gx);
          nextLevel->o_Sx = platform->device.malloc(
            ellipticCoarse->mesh->Np * ellipticCoarse->mesh->Nelements * sizeof(dfloat),
            nextLevel->Sx);
        } else {
          //this level is the base
          parAlmond::coarseSolver* coarseLevel = precon->parAlmond->coarseLevel;

          coarseLevel->gatherLevel = true;
          coarseLevel->ogs = ellipticCoarse->ogs;
          coarseLevel->Gx = (dfloat*) calloc(coarseLevel->ogs->Ngather,sizeof(dfloat));
          coarseLevel->Sx = (dfloat*) calloc(ellipticCoarse->mesh->Np * ellipticCoarse->mesh->Nelements,
                                             sizeof(dfloat));
          coarseLevel->o_Gx =
            platform->device.malloc(coarseLevel->ogs->Ngather * sizeof(dfloat),
                                                coarseLevel->Gx);
          coarseLevel->o_Sx = platform->device.malloc(
            ellipticCoarse->mesh->Np * ellipticCoarse->mesh->Nelements * sizeof(dfloat),
            coarseLevel->Sx);
        }
      }
    }
  }

  for (int n = 1; n < mesh->N + 1; n++)
    delete meshLevels[n];
  //  for (int n=1;n<mesh->N+1;n++) delete[] meshLevels[n];
  free(meshLevels);

  if (platform->comm.mpiRank == 0) {
    printf("---------------------------------------------------------\n");
    printf("level|    Type    |                 |     Smoother      |\n");
    printf("     |            |                 |                   |\n");
    printf("---------------------------------------------------------\n");
  }

  for(int lev = 0; lev < precon->parAlmond->numLevels; lev++) {
    if(platform->comm.mpiRank == 0) {
      printf(" %3d ", lev);
      fflush(stdout);
    }
    levels[lev]->Report();
  }

  if (platform->comm.mpiRank == 0)
    printf("---------------------------------------------------------\n");
}

void MGLevelAllocateStorage(MGLevel* level, int k, parAlmond::CycleType ctype)
{
  
  // extra storage for smoothing op
  size_t Nbytes = level->Ncols * sizeof(pfloat);
  if (MGLevel::smootherResidualBytes < Nbytes) {
    if (MGLevel::o_smootherResidual.size()) {
      free(MGLevel::smootherResidual);
      MGLevel::o_smootherResidual.free();
      MGLevel::o_smootherResidual2.free();
      MGLevel::o_smootherUpdate.free();
    }

    MGLevel::smootherResidual = (pfloat*) calloc(level->Ncols,sizeof(pfloat));
    MGLevel::o_smootherResidual = platform->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::o_smootherResidual2 = platform->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::o_smootherUpdate = platform->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::smootherResidualBytes = Nbytes;
  }

  if (k) level->x    = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
  if (k) level->rhs  = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
  if (k) level->o_x   = platform->device.malloc(level->Ncols * sizeof(dfloat),level->x);
  if (k) level->o_rhs = platform->device.malloc(level->Nrows * sizeof(dfloat),level->rhs);

  level->res  = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
  level->o_res = platform->device.malloc(level->Ncols * sizeof(dfloat),level->res);

  //kcycle vectors
  if (ctype == parAlmond::KCYCLE) {
    if ((k > 0) && (k < NUMKCYCLES + 1)) {
      level->ck = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
      level->vk = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
      level->wk = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
      level->o_ck = platform->device.malloc(level->Ncols * sizeof(dfloat),level->ck);
      level->o_vk = platform->device.malloc(level->Nrows * sizeof(dfloat),level->vk);
      level->o_wk = platform->device.malloc(level->Nrows * sizeof(dfloat),level->wk);
    }
  }
}
