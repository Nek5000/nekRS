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

void ellipticMultiGridSetup(elliptic_t* elliptic_, precon_t* precon)
{
  // setup new object with constant coeff
  elliptic_t* elliptic = ellipticBuildMultigridLevelFine(elliptic_);
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  const dfloat lambda = elliptic->lambda[0];

  //read all the nodes files and load them in a dummy mesh array
  mesh_t** meshLevels = (mesh_t**) calloc(mesh->N + 1,sizeof(mesh_t*));
  for (int n = 1; n < mesh->N + 1; n++) {
    meshLevels[n] = new mesh_t();
    meshLevels[n]->Nverts = mesh->Nverts;
    meshLevels[n]->Nfaces = mesh->Nfaces;
    meshLevels[n]->Nfields = mesh->Nfields; // TW: ahem

    switch(elliptic->elementType) {
    case TRIANGLES:
      meshLoadReferenceNodesTri2D(meshLevels[n], n);
      break;
    case QUADRILATERALS:
      meshLoadReferenceNodesQuad2D(meshLevels[n], n);
      break;
    case TETRAHEDRA:
      meshLoadReferenceNodesTet3D(meshLevels[n], n);
      break;
    case HEXAHEDRA:
      meshLoadReferenceNodesHex3D(meshLevels[n], n, 1);
      break;
    }
  }

  //set the number of MG levels and their degree
  int numMGLevels;
  int* levelDegree;

  if (options.compareArgs("MULTIGRID COARSENING","CUSTOM")) {
    numMGLevels = elliptic->nLevels;
    levelDegree = (int*) calloc(numMGLevels,sizeof(int));
    for(int i = 0; i < numMGLevels; ++i)
      levelDegree[i] = elliptic->levels[i];
  }else if (options.compareArgs("MULTIGRID COARSENING","ALLDEGREES")) {
    numMGLevels = mesh->N;
    levelDegree = (int*) calloc(numMGLevels,sizeof(int));
    for (int n = 0; n < numMGLevels; n++) levelDegree[n] = mesh->N - n; //all degrees
  }else if (options.compareArgs("MULTIGRID COARSENING","HALFDEGREES")) {
    numMGLevels = floor(mesh->N / 2.) + 1;
    levelDegree = (int*) calloc(numMGLevels,sizeof(int));
    for (int n = 0; n < numMGLevels; n++) levelDegree[n] = mesh->N - 2 * n; //decrease by two
    levelDegree[numMGLevels - 1] = 1; //ensure the last level is degree 1
  }else {  //default "HALFDOFS"
    // pick the degrees so the dofs of each level halfs (roughly)
    //start by counting the number of levels neccessary
    numMGLevels = 1;
    int degree = mesh->N;
    int dofs = meshLevels[degree]->Np;
    int basedofs = mesh->Nverts;
    while (dofs > basedofs) {
      numMGLevels++;
      for (; degree > 0; degree--)
        if (meshLevels[degree]->Np <= dofs / 2)
          break;
      dofs = meshLevels[degree]->Np;
    }
    levelDegree = (int*) calloc(numMGLevels,sizeof(int));
    degree = mesh->N;
    numMGLevels = 1;
    levelDegree[0] = degree;
    dofs = meshLevels[degree]->Np;
    while (dofs > basedofs) {
      for (; degree > 0; degree--)
        if (meshLevels[degree]->Np <= dofs / 2)
          break;
      dofs = meshLevels[degree]->Np;
      levelDegree[numMGLevels] = degree;
      numMGLevels++;
    }
  }

  int Nmax = levelDegree[0];
  int Nmin = levelDegree[numMGLevels - 1];

  //initialize parAlmond
  precon->parAlmond = parAlmond::Init(mesh->device, mesh->comm, options);
  parAlmond::multigridLevel** levels = precon->parAlmond->levels;

  oogs_mode oogsMode = OOGS_AUTO;
  if(options.compareArgs("THREAD MODEL", "SERIAL")) oogsMode = OOGS_DEFAULT;
  if(options.compareArgs("THREAD MODEL", "OPENMP")) oogsMode = OOGS_DEFAULT;

  //set up the finest level
  if (Nmax > Nmin) {
    if(mesh->rank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nmax);

    auto callback = [&]()
      {
        ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
                   elliptic->o_p, elliptic->o_Ap, pfloatString);
      };
    elliptic->oogs   = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    elliptic->oogsAx = oogs::setup(elliptic->ogs, 1, 0, ogsPfloat, callback, oogsMode);

    levels[0] = new MGLevel(elliptic, lambda, Nmax, options,
                            precon->parAlmond->ktype, mesh->comm);
    MGLevelAllocateStorage((MGLevel*) levels[0], 0,
                           precon->parAlmond->ctype);
    precon->parAlmond->numLevels++;
  }

  //build a MGLevel for every degree (except degree 1)
  for (int n = 1; n < numMGLevels - 1; n++) {
    int Nc = levelDegree[n];
    int Nf = levelDegree[n - 1];
    //build elliptic struct for this degree
    if(mesh->rank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nc);

    elliptic_t* ellipticC = ellipticBuildMultigridLevel(elliptic,Nc,Nf);

    auto callback = [&]()
      {
        ellipticAx(ellipticC, ellipticC->mesh->NlocalGatherElements, ellipticC->mesh->o_localGatherElementList,
                   ellipticC->o_p, ellipticC->o_Ap, pfloatString);
      };
    ellipticC->oogs   = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    ellipticC->oogsAx = oogs::setup(ellipticC->ogs, 1, 0, ogsPfloat, callback, oogsMode);

    //add the level manually
    levels[n] = new MGLevel(elliptic,
                            meshLevels,
                            ((MGLevel*) levels[n - 1])->elliptic,
                            ellipticC,
                            lambda,
                            Nf, Nc,
                            options,
                            precon->parAlmond->ktype, mesh->comm);
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

    if(mesh->rank == 0)
      printf("=============BUILDING MULTIGRID LEVEL OF DEGREE %d==================\n", Nmin);

    ellipticCoarse = ellipticBuildMultigridLevel(elliptic,Nc,Nf);

    auto callback = [&]()
      {
        ellipticAx(ellipticCoarse, ellipticCoarse->mesh->NlocalGatherElements, ellipticCoarse->mesh->o_localGatherElementList,
                   ellipticCoarse->o_p, ellipticCoarse->o_Ap, pfloatString);
      };
    ellipticCoarse->oogs   = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, NULL, oogsMode);
    //ellipticCoarse->oogsAx = oogs::setup(ellipticCoarse->ogs, 1, 0, ogsPfloat, callback, oogsMode);
  } else {
    ellipticCoarse = elliptic;
  }
  int basisNp = ellipticCoarse->mesh->Np;
  dfloat* basis = NULL;

  if (options.compareArgs("BASIS","BERN")) basis = ellipticCoarse->mesh->VB;

  hlong* coarseGlobalStarts = (hlong*) calloc(mesh->size + 1, sizeof(hlong));

  if (options.compareArgs("DISCRETIZATION","IPDG")) {
    ellipticBuildIpdg(ellipticCoarse,
                      basisNp,
                      basis,
                      lambda,
                      &coarseA,
                      &nnzCoarseA,
                      coarseGlobalStarts);
  } else if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    if(options.compareArgs("GALERKIN COARSE OPERATOR","TRUE"))
      ellipticBuildContinuousGalerkinHex3D(ellipticCoarse,elliptic,lambda,&coarseA,&nnzCoarseA,
                                           &coarseogs,coarseGlobalStarts);
    else
      ellipticBuildContinuous(ellipticCoarse, &coarseA, &nnzCoarseA,&coarseogs,
                              coarseGlobalStarts);
  }

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
                                          lambda,
                                          Nf, Nc,
                                          options,
                                          precon->parAlmond->ktype, mesh->comm
                                          );
  } else {
    levels[numMGLevels - 1] = new MGLevel(ellipticCoarse, lambda, Nmin, options,
                                          precon->parAlmond->ktype, mesh->comm);
  }
  MGLevelAllocateStorage((MGLevel*) levels[numMGLevels - 1], numMGLevels - 1,
                         precon->parAlmond->ctype);

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
        ellipticCoarse->mesh->device.malloc(levels[numMGLevels - 1]->Ncols * sizeof(dfloat),
                                            nextLevel->Gx);
      nextLevel->o_Sx = ellipticCoarse->mesh->device.malloc(
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
        ellipticCoarse->mesh->device.malloc(coarseLevel->ogs->Ngather * sizeof(dfloat),
                                            coarseLevel->Gx);
      coarseLevel->o_Sx = ellipticCoarse->mesh->device.malloc(
        ellipticCoarse->mesh->Np * ellipticCoarse->mesh->Nelements * sizeof(dfloat),
        coarseLevel->Sx);
    }
  }

  for (int n = 1; n < mesh->N + 1; n++)
    delete meshLevels[n];
  //  for (int n=1;n<mesh->N+1;n++) delete[] meshLevels[n];
  free(meshLevels);

  //report top levels
  if (mesh->rank == 0) { //report the upper multigrid levels
    printf("--------------------Multigrid Report---------------------\n");
    printf("---------------------------------------------------------\n");
    printf("level|    Type    |                 |     Smoother      |\n");
    printf("     |            |                 |                   |\n");
    printf("---------------------------------------------------------\n");
  }

  for(int lev = 0; lev < precon->parAlmond->numLevels; lev++) {
    if(mesh->rank == 0) {
      printf(" %3d ", lev);
      fflush(stdout);
    }
    levels[lev]->Report();
  }

  if (mesh->rank == 0)
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
    MGLevel::o_smootherResidual = level->mesh->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::o_smootherResidual2 = level->mesh->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::o_smootherUpdate = level->mesh->device.malloc(Nbytes,MGLevel::smootherResidual);
    MGLevel::smootherResidualBytes = Nbytes;
  }

  if (k) level->x    = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
  if (k) level->rhs  = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
  if (k) level->o_x   = level->mesh->device.malloc(level->Ncols * sizeof(dfloat),level->x);
  if (k) level->o_rhs = level->mesh->device.malloc(level->Nrows * sizeof(dfloat),level->rhs);

  level->res  = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
  level->o_res = level->mesh->device.malloc(level->Ncols * sizeof(dfloat),level->res);

  //kcycle vectors
  if (ctype == parAlmond::KCYCLE) {
    if ((k > 0) && (k < NUMKCYCLES + 1)) {
      level->ck = (dfloat*) calloc(level->Ncols,sizeof(dfloat));
      level->vk = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
      level->wk = (dfloat*) calloc(level->Nrows,sizeof(dfloat));
      level->o_ck = level->mesh->device.malloc(level->Ncols * sizeof(dfloat),level->ck);
      level->o_vk = level->mesh->device.malloc(level->Nrows * sizeof(dfloat),level->vk);
      level->o_wk = level->mesh->device.malloc(level->Nrows * sizeof(dfloat),level->wk);
    }
  }
}
