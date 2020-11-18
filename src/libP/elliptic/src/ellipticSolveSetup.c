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
#include <string>

void ellipticSolveSetup(elliptic_t* elliptic, occa::properties kernelInfo)
{
  mesh_t* mesh      = elliptic->mesh;
  setupAide options = elliptic->options;

  const dlong Nlocal = mesh->Np * mesh->Nelements;
  elliptic->resNormFactor = 1 / (elliptic->Nfields*mesh->volume);

  const int serial = options.compareArgs("THREAD MODEL", "SERIAL");

  if (elliptic->blockSolver &&  elliptic->elementType != HEXAHEDRA &&
      !options.compareArgs("DISCRETIZATION",
                           "CONTINUOUS") && !options.compareArgs("PRECONDITIONER","JACOBI") ) {
    if(mesh->rank == 0)
      printf("ERROR: Block solver is implemented for C0-HEXAHEDRA with Jacobi preconditioner only\n");

    MPI_Finalize();
    exit(-1);
  }

  // Sanity check for discretization type
  if (options.compareArgs("COEFFICIENT","VARIABLE") &&  elliptic->elementType != HEXAHEDRA &&
      !options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    if(mesh->rank == 0)
      printf("ERROR: Varibale coefficient solver is implemented for C0-HEXAHEDRA only\n");

    MPI_Finalize();
    exit(-1);
  }

  // Sanity check for preconditioner type
  if (options.compareArgs("COEFFICIENT","VARIABLE")) {
    if(options.compareArgs("PRECONDITIONER",
                           "MULTIGRID") &&
       !options.compareArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE")) {
      if(mesh->rank == 0)
        printf(
          "ERROR: Varibale coefficient solver is implemented for constant multigrid preconditioner only\n");
      MPI_Finalize();
      exit(-1);
    }

    if(!options.compareArgs("PRECONDITIONER",
                            "MULTIGRID") && !options.compareArgs("PRECONDITIONER", "JACOBI")
       && !options.compareArgs("PRECONDITIONER", "NONE")) {
      if(mesh->rank == 0)
        printf(
          "ERROR: Varibale coefficient solver is implemented for multigrid/Jacobi/None preconditioners only\n");

      MPI_Finalize();
      exit(-1);
    }
  }

  //sanity checking
  if (options.compareArgs("BASIS","BERN") && elliptic->elementType != TRIANGLES) {
    printf("ERROR: BERN basis is only available for triangular elements\n");
    MPI_Finalize();
    exit(-1);
  }

  if (options.compareArgs("PRECONDITIONER","MASSMATRIX") && elliptic->elementType != TRIANGLES
      && elliptic->elementType != TETRAHEDRA ) {
    printf(
      "ERROR: MASSMATRIX preconditioner is only available for triangle and tetrhedra elements. Use JACOBI instead.\n");
    MPI_Finalize();
    exit(-1);
  }

  dlong Nblock  = mymax(1,(Nlocal + blockSize - 1) / blockSize);
  dlong Nblock2 = mymax(1,(Nblock + blockSize - 1) / blockSize);

  dlong NthreadsUpdatePCG = 256;
  dlong NblocksUpdatePCG = mymin((Nlocal + NthreadsUpdatePCG - 1) / NthreadsUpdatePCG, 160);

  elliptic->NthreadsUpdatePCG = NthreadsUpdatePCG;
  elliptic->NblocksUpdatePCG = NblocksUpdatePCG;

  //tau
  if (elliptic->elementType == TRIANGLES ||
      elliptic->elementType == QUADRILATERALS) {
    elliptic->tau = 2.0 * (mesh->N + 1) * (mesh->N + 2) / 2.0;
    if(elliptic->dim == 3)
      elliptic->tau *= 1.5;
  }else {
    elliptic->tau = 2.0 * (mesh->N + 1) * (mesh->N + 3);
  }

  // Assumes wrkoffset is set properly, i.e. workoffset = wrkoffset*Nfields
  if (elliptic->wrk) { // user-provided scratch space
    elliptic->p    = elliptic->wrk + 0 * elliptic->Ntotal * elliptic->Nfields;
    elliptic->z    = elliptic->wrk + 1 * elliptic->Ntotal * elliptic->Nfields;
    elliptic->Ap   = elliptic->wrk + 2 * elliptic->Ntotal * elliptic->Nfields;
    elliptic->grad = elliptic->wrk + 4 * elliptic->Ntotal * elliptic->Nfields;

    elliptic->o_p    =
      elliptic->o_wrk.slice(0 * elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
    elliptic->o_z    =
      elliptic->o_wrk.slice(1 * elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
    elliptic->o_Ap   =
      elliptic->o_wrk.slice(2 * elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
    elliptic->o_rtmp =
      elliptic->o_wrk.slice(3 * elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));

    //elliptic->o_grad =
    //  elliptic->o_wrk.slice(4 * elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
  } else {
    elliptic->p    = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));
    elliptic->z    = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));
    elliptic->Ap   = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));
    elliptic->grad = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields * 4, sizeof(dfloat));

    elliptic->o_p    = mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                           elliptic->p);
    elliptic->o_z    = mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                           elliptic->z);
    elliptic->o_Ap   = mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                           elliptic->Ap);
    elliptic->o_rtmp = mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                           elliptic->p);
    elliptic->o_grad = mesh->device.malloc(
      elliptic->Ntotal * elliptic->Nfields * 4 * sizeof(dfloat),
      elliptic->grad);
  }

  elliptic->o_x0 = mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));

  elliptic->tmp = (dfloat*) calloc(Nblock, sizeof(dfloat));
  elliptic->o_tmp = mesh->device.malloc(Nblock * sizeof(dfloat), elliptic->tmp);
  elliptic->o_tmp2 = mesh->device.malloc(Nblock2 * sizeof(dfloat), elliptic->tmp);

  elliptic->tmpNormr = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
  elliptic->o_tmpNormr = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                             elliptic->tmpNormr);

  int useFlexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");

  if(options.compareArgs("KRYLOV SOLVER", "NONBLOCKING")) {
    int Nwork = (useFlexible) ? 9: 6;

    elliptic->o_pcgWork = new occa::memory[Nwork];
    for(int n = 0; n < Nwork; ++n)
      elliptic->o_pcgWork[n]  =
        mesh->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat), elliptic->z);

    elliptic->tmppdots = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmppdots = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmppdots);

    elliptic->tmprdotz = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmprdotz = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmprdotz);

    elliptic->tmpzdotz = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmpzdotz = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmpzdotz);

    elliptic->tmprdotr = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmprdotr = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmprdotr);

    elliptic->tmpudotr = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmpudotr = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmpudotr);

    elliptic->tmpudots = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmpudots = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmpudots);

    elliptic->tmpudotw = (dfloat*) calloc(elliptic->NblocksUpdatePCG,sizeof(dfloat));
    elliptic->o_tmpudotw = mesh->device.malloc(elliptic->NblocksUpdatePCG * sizeof(dfloat),
                                               elliptic->tmpudotw);
  }

  //setup async halo stream
  elliptic->defaultStream = mesh->defaultStream;
  elliptic->dataStream = mesh->dataStream;

  dlong Nbytes = mesh->totalHaloPairs * mesh->Np * sizeof(dfloat);
  if(Nbytes > 0) {
    elliptic->sendBuffer = (dfloat*) occaHostMallocPinned(mesh->device,
                                                          Nbytes * elliptic->Nfields,
                                                          NULL,
                                                          elliptic->o_sendBuffer,
                                                          elliptic->h_sendBuffer);
    elliptic->recvBuffer = (dfloat*) occaHostMallocPinned(mesh->device,
                                                          Nbytes * elliptic->Nfields,
                                                          NULL,
                                                          elliptic->o_recvBuffer,
                                                          elliptic->h_recvBuffer);
    elliptic->gradSendBuffer = (dfloat*) occaHostMallocPinned(mesh->device,
                                                              2 * Nbytes * elliptic->Nfields,
                                                              NULL,
                                                              elliptic->o_gradSendBuffer,
                                                              elliptic->h_gradSendBuffer);
    elliptic->gradRecvBuffer = (dfloat*) occaHostMallocPinned(mesh->device,
                                                              2 * Nbytes * elliptic->Nfields,
                                                              NULL,
                                                              elliptic->o_gradRecvBuffer,
                                                              elliptic->h_gradRecvBuffer);
  }else{
    elliptic->sendBuffer = NULL;
    elliptic->recvBuffer = NULL;
  }
  mesh->device.setStream(elliptic->defaultStream);

  elliptic->type = strdup(dfloatString);

  elliptic->Nblock = Nblock;
  elliptic->Nblock2 = Nblock2;

  //fill geometric factors in halo
  if(mesh->totalHaloPairs) {
    dlong Nlocal = mesh->Nelements;
    dlong Nhalo = mesh->totalHaloPairs;
    size_t Nbytes = mesh->Nvgeo * sizeof(dfloat);

    if (elliptic->elementType == QUADRILATERALS || elliptic->elementType == HEXAHEDRA) {
      Nlocal *= mesh->Np;
      Nhalo *= mesh->Np;
      Nbytes *= mesh->Np;
    }

    dfloat* vgeoSendBuffer = (dfloat*) calloc(Nhalo * mesh->Nvgeo, sizeof(dfloat));

    // import geometric factors from halo elements
    mesh->vgeo = (dfloat*) realloc(mesh->vgeo, (Nlocal + Nhalo) * mesh->Nvgeo * sizeof(dfloat));

    meshHaloExchange(mesh,
                     Nbytes,
                     mesh->vgeo,
                     vgeoSendBuffer,
                     mesh->vgeo + Nlocal * mesh->Nvgeo);

    mesh->o_vgeo =
      mesh->device.malloc((Nlocal + Nhalo) * mesh->Nvgeo * sizeof(dfloat), mesh->vgeo);
    free(vgeoSendBuffer);
  }

  // count total number of elements
  hlong NelementsLocal = mesh->Nelements;
  hlong NelementsGlobal = 0;

  MPI_Allreduce(&NelementsLocal, &NelementsGlobal, 1, MPI_HLONG, MPI_SUM, mesh->comm);

  elliptic->NelementsGlobal = NelementsGlobal;

  elliptic->allNeumannPenalty = 1.;
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  elliptic->allNeumannScale = 1. / sqrt((dfloat)mesh->Np * totalElements);

  elliptic->allNeumannPenalty = 0;
  elliptic->allNeumannScale = 0;

  elliptic->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces * elliptic->Nfields,sizeof(int));
  int* allNeumann = (int*)calloc(elliptic->Nfields, sizeof(int));
  // check based on the coefficient
  for(int fld = 0; fld < elliptic->Nfields; fld++) {
    if(elliptic->var_coeff) {
      int allzero = 1;
      for(int n = 0; n < Nlocal; n++) { // check any non-zero value for each field
        const dfloat lambda = elliptic->lambda[n + elliptic->Ntotal + fld * elliptic->loffset];
        if(lambda) {
          allzero = 0;
          break;
        }
      }
      allNeumann[fld] = allzero;
    }else{
      allNeumann[fld] = (elliptic->lambda[fld] == 0) ? 1 : 0;
    }
  }

  // check based on BC
  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong e = 0; e < mesh->Nelements; e++)
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = mesh->EToB[e * mesh->Nfaces + f];
        if (bc > 0) {
          int BC = elliptic->BCType[bc + elliptic->NBCType * fld];
          elliptic->EToB[f + e * mesh->Nfaces + fld * mesh->Nelements * mesh->Nfaces] = BC; //record it
          if (BC != 2) allNeumann[fld] = 0; //check if its a Dirchlet for each field
        }
      }

  elliptic->allNeumann = 0;
  elliptic->allBlockNeumann = (int*)calloc(elliptic->Nfields, sizeof(int));
  for(int fld = 0; fld < elliptic->Nfields; fld++) {
    int lallNeumann, gallNeumann;
    lallNeumann = allNeumann[fld] ? 0:1;
    MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, mesh->comm);
    elliptic->allBlockNeumann[fld] = (gallNeumann > 0) ? 0: 1;
    // even if there is a single allNeumann activate Null space correction
    if(elliptic->allBlockNeumann[fld])
      elliptic->allNeumann = 1;
  }

  if(mesh->rank == 0)
    printf("allNeumann = %d \n", elliptic->allNeumann);

  //set surface mass matrix for continuous boundary conditions
  mesh->sMT = (dfloat*) calloc(mesh->Np * mesh->Nfaces * mesh->Nfp,sizeof(dfloat));
  for (int n = 0; n < mesh->Np; n++)
    for (int m = 0; m < mesh->Nfp * mesh->Nfaces; m++) {
      dfloat MSnm = 0;
      for (int i = 0; i < mesh->Np; i++)
        MSnm += mesh->MM[n + i * mesh->Np] * mesh->LIFT[m + i * mesh->Nfp * mesh->Nfaces];
      mesh->sMT[n + m * mesh->Np]  = MSnm;
    }
  mesh->o_sMT =
    mesh->device.malloc(mesh->Np * mesh->Nfaces * mesh->Nfp * sizeof(dfloat), mesh->sMT);

  //copy boundary flags
  elliptic->o_EToB = mesh->device.malloc(
    mesh->Nelements * mesh->Nfaces * elliptic->Nfields * sizeof(int),
    elliptic->EToB);

#if 0
  if (mesh->rank == 0 && options.compareArgs("VERBOSE","TRUE"))
    occa::setVerboseCompilation(true);
  else
    occa::setVerboseCompilation(false);
#endif

  //setup an unmasked gs handle
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  if(mesh->ogs == NULL) meshParallelGatherScatterSetup(mesh, Nlocal, mesh->globalIds, mesh->comm, verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  const int mapSize = elliptic->blockSolver ? elliptic->Ntotal * elliptic->Nfields: Nlocal;
  elliptic->mapB = (int*) calloc(mapSize,sizeof(int));
  int largeNumber = 1 << 20;

  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong e = 0; e < mesh->Nelements; e++) {
      for (int n = 0; n < mesh->Np; n++)
        elliptic->mapB[n + e * mesh->Np + fld * elliptic->Ntotal] = largeNumber;
      for (int f = 0; f < mesh->Nfaces; f++) {
        int bc = mesh->EToB[f + e * mesh->Nfaces];
        if (bc > 0) {
          for (int n = 0; n < mesh->Nfp; n++) {
            int BCFlag = elliptic->BCType[bc + elliptic->NBCType * fld];
            int fid = mesh->faceNodes[n + f * mesh->Nfp];
            elliptic->mapB[fid + e * mesh->Np + fld * elliptic->Ntotal] = mymin(BCFlag,
                                                  elliptic->mapB[fid + e *mesh->Np + fld*elliptic->Ntotal]);
          }
        }
      }
    }

  if(elliptic->blockSolver)
    ogsGatherScatterMany(elliptic->mapB,
                         elliptic->Nfields,
                         elliptic->Ntotal,
                         ogsInt,
                         ogsMin,
                         mesh->ogs);
  else
    ogsGatherScatter(elliptic->mapB, ogsInt, ogsMin, mesh->ogs);

  // Create field based
  elliptic->Nmasked      = 0;
  elliptic->fNmasked     = (dlong*)calloc(elliptic->Nfields, sizeof(dlong));
  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++) {
      if (elliptic->mapB[n + fld * elliptic->Ntotal] == largeNumber) {
        elliptic->mapB[n + fld * elliptic->Ntotal] = 0.;
      } else if (elliptic->mapB[n + fld * elliptic->Ntotal] == 1) { //Dirichlet boundary
        elliptic->Nmasked++; // increase global accumulator
        elliptic->fNmasked[fld]++; // increase local accumulator
      }
    }
  elliptic->o_mapB = mesh->device.malloc(mapSize * sizeof(int), elliptic->mapB);

  elliptic->maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
  elliptic->Nmasked = 0; //reset
  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
      if (elliptic->mapB[n + fld * elliptic->Ntotal] == 1)
        elliptic->maskIds[elliptic->Nmasked++] = n + fld * elliptic->Ntotal;
  if (elliptic->Nmasked) elliptic->o_maskIds = mesh->device.malloc(
      elliptic->Nmasked * sizeof(dlong),
      elliptic->maskIds);


  if(elliptic->blockSolver) { // Create a gs handle independent from BC handler
    elliptic->ogs = ogsSetup(Nlocal, mesh->globalIds, mesh->comm, verbose, mesh->device);
    // Create copy of invDegree so that we can accelerate vector form of masking!!!!!!
    elliptic->invDegree = (dfloat*)calloc(elliptic->Ntotal * elliptic->Nfields, sizeof(dfloat));

    for(int n = 0; n < elliptic->Ntotal * elliptic->Nfields; n++) elliptic->invDegree[n] = 1.0;
    for(int fld = 0; fld < elliptic->Nfields; fld++)
      for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
        elliptic->invDegree[n + fld * elliptic->Ntotal] = elliptic->ogs->invDegree[n];

    elliptic->o_invDegree = mesh->device.malloc(
      elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
      elliptic->invDegree);
  }else{ //make a masked version of the global id numbering if scalar
    mesh->maskedGlobalIds = (hlong*) calloc(Nlocal,sizeof(hlong));
    memcpy(mesh->maskedGlobalIds, mesh->globalIds, Nlocal * sizeof(hlong));
    for (dlong n = 0; n < elliptic->Nmasked; n++)
      mesh->maskedGlobalIds[elliptic->maskIds[n]] = 0;

    //use the masked ids to make another gs handle
    elliptic->ogs = ogsSetup(Nlocal, mesh->maskedGlobalIds, mesh->comm, verbose, mesh->device);
    elliptic->o_invDegree = elliptic->ogs->o_invDegree;
  }

  /*preconditioner setup */
  elliptic->precon = new precon_t();

  //  kernelInfo["parser/" "automate-add-barriers"] =  "disabled";
  kernelInfo["defines/pfloat"] = pfloatString;


  // set kernel name suffix
  char* suffix;

  if(elliptic->elementType == TRIANGLES) {
    if(elliptic->dim == 2)
      suffix = strdup("Tri2D");
    else
      suffix = strdup("Tri3D");
  }
  if(elliptic->elementType == QUADRILATERALS) {
    if(elliptic->dim == 2)
      suffix = strdup("Quad2D");
    else
      suffix = strdup("Quad3D");
  }
  if(elliptic->elementType == TETRAHEDRA)
    suffix = strdup("Tet3D");
  if(elliptic->elementType == HEXAHEDRA)
    suffix = strdup("Hex3D");

  char fileName[BUFSIZ], kernelName[BUFSIZ];

  kernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;
  kernelInfo["defines/p_Nalign"] = USE_OCCA_MEM_BYTE_ALIGN;
  kernelInfo["defines/" "p_blockSize"] = blockSize;

  occa::properties pfloatKernelInfo = kernelInfo;
  pfloatKernelInfo["defines/dfloat"] = pfloatString;
  pfloatKernelInfo["defines/pfloat"] = pfloatString;

  occa::properties kernelInfoNoOKL = kernelInfo;
  if(serial) kernelInfoNoOKL["okl/enabled"] = false;

  MPI_Barrier(mesh->comm);
  double tStartLoadKernel = MPI_Wtime(); 
  if(mesh->rank == 0)  printf("loading elliptic kernels ... "); fflush(stdout);

  for (int r = 0; r < 2; r++) {
    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {
      //mesh kernels
      mesh->haloExtractKernel =
        mesh->device.buildKernel(DHOLMES "/okl/meshHaloExtract2D.okl",
                                 "meshHaloExtract2D",
                                 kernelInfo);

      mesh->addScalarKernel =
        mesh->device.buildKernel(DHOLMES "/okl/addScalar.okl",
                                 "addScalar",
                                 kernelInfo);

      mesh->maskKernel =
        mesh->device.buildKernel(DHOLMES "/okl/mask.okl",
                                 "mask",
                                 kernelInfo);
      mesh->maskPfloatKernel =
        mesh->device.buildKernel(DHOLMES "/okl/mask.okl",
                                 "mask",
                                 pfloatKernelInfo);

      mesh->sumKernel =
        mesh->device.buildKernel(DHOLMES "/okl/sum.okl",
                                 "sum",
                                 kernelInfo);

      elliptic->fillKernel =
        mesh->device.buildKernel(DHOLMES "/okl/fill.okl",
                                 "fill",
                                 kernelInfo);

      elliptic->dotMultiplyAddKernel =
        mesh->device.buildKernel(DHOLMES "/okl/dotMultiplyAdd.okl",
                                 "dotMultiplyAdd",
                                 kernelInfo);

      elliptic->dotDivideKernel =
        mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
                                 "dotDivide",
                                 kernelInfo);

      elliptic->scalarDivideKernel =
        mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
                                 "scalarDivide",
                                 kernelInfo);

      if(elliptic->blockSolver) {
        elliptic->sumBlockKernel =
          mesh->device.buildKernel(DHOLMES "/okl/sum.okl", "sumBlock", kernelInfo);

        elliptic->sumBlockFieldKernel =
          mesh->device.buildKernel(DHOLMES "/okl/sum.okl", "sumBlockField", kernelInfo);

        elliptic->addScalarBlockFieldKernel =
          mesh->device.buildKernel(DHOLMES "/okl/addScalar.okl",
                                   "addBlockScalarField",
                                   kernelInfo);

        elliptic->weightedInnerProduct1Kernel =
          mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
                                   "weightedBlockInnerProduct1",
                                   kernelInfo);
        if(serial)
          elliptic->weightedInnerProduct2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialWeightedInnerProduct2.c",
                                     "weightedBlockInnerProduct2",
                                     kernelInfoNoOKL);
        else
          elliptic->weightedInnerProduct2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
                                     "weightedBlockInnerProduct2",
                                     kernelInfo);

        elliptic->innerProductKernel =
          mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
                                   "innerBlockProduct",
                                   kernelInfo);

        elliptic->innerProductFieldKernel =
          mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
                                   "innerBlockProductField",
                                   kernelInfo);

        if(serial)
          elliptic->weightedNorm2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialWeightedNorm2.c",
                                     "weightedBlockNorm2",
                                     kernelInfoNoOKL);
        else
          elliptic->weightedNorm2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
                                     "weightedBlockNorm2",
                                     kernelInfo);

        elliptic->norm2Kernel =
          mesh->device.buildKernel(DHOLMES "/okl/norm2.okl",
                                   "normBlock2",
                                   kernelInfo);
        elliptic->scaledAddPfloatKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                   "scaledBlockAdd",
                                   pfloatKernelInfo);

        if(serial)
          elliptic->scaledAddKernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialScaledAdd.c",
                                     "scaledBlockAdd",
                                     kernelInfoNoOKL);
        else{
          elliptic->scaledAddKernel =
            mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                     "scaledBlockAdd",
                                     kernelInfo);
        }

        elliptic->collocateKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                   "collocate",
                                   kernelInfo);
          elliptic->dotMultiplyPfloatKernel =
            mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                     "dotBlockMultiply",
                                     pfloatKernelInfo);

        if(serial)
          elliptic->dotMultiplyKernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialDotMultiply.c",
                                     "dotBlockMultiply",
                                     kernelInfoNoOKL);
        else
          elliptic->dotMultiplyKernel =
            mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                     "dotBlockMultiply",
                                     kernelInfo);

        elliptic->scalarDivideManyKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotDivide.okl",
                                   "scalarDivideMany",
                                   kernelInfo);
      }else{
        elliptic->weightedInnerProduct1Kernel =
          mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct1.okl",
                                   "weightedInnerProduct1",
                                   kernelInfo);

        if(serial)
          elliptic->weightedInnerProduct2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialWeightedInnerProduct2.c",
                                     "weightedInnerProduct2",
                                     kernelInfoNoOKL);
        else
          elliptic->weightedInnerProduct2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/weightedInnerProduct2.okl",
                                     "weightedInnerProduct2",
                                     kernelInfo);

        elliptic->innerProductKernel =
          mesh->device.buildKernel(DHOLMES "/okl/innerProduct.okl",
                                   "innerProduct",
                                   kernelInfo);

        if(serial)
          elliptic->weightedNorm2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialWeightedNorm2.c",
                                     "weightedNorm2",
                                     kernelInfoNoOKL);
        else
          elliptic->weightedNorm2Kernel =
            mesh->device.buildKernel(DHOLMES "/okl/weightedNorm2.okl",
                                     "weightedNorm2",
                                     kernelInfo);

        elliptic->norm2Kernel =
          mesh->device.buildKernel(DHOLMES "/okl/norm2.okl",
                                   "norm2",
                                   kernelInfo);
        elliptic->scaledAddPfloatKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                   "scaledAdd",
                                   pfloatKernelInfo);
        elliptic->copyDfloatToPfloatKernel =
          mesh->device.buildKernel(DHOLMES "/okl/copyDfloatToPfloat.okl",
                                   "copyDfloatToPfloat",
                                   kernelInfo);
        elliptic->copyPfloatToDPfloatKernel =
          mesh->device.buildKernel(DHOLMES "/okl/copyPfloatToDfloat.okl",
                                   "copyPfloatToDfloat",
                                   kernelInfo);
        elliptic->updateSmoothedSolutionVecKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                   "updateSmoothedSolutionVec",
                                   pfloatKernelInfo);
        elliptic->updateChebyshevSolutionVecKernel =
          mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                   "updateChebyshevSolutionVec",
                                   pfloatKernelInfo);

        if(serial)
          elliptic->scaledAddKernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialScaledAdd.c",
                                     "scaledAdd",
                                     kernelInfoNoOKL);
        else{
          elliptic->scaledAddKernel =
            mesh->device.buildKernel(DHOLMES "/okl/scaledAdd.okl",
                                     "scaledAdd",
                                     kernelInfo);
        }

        elliptic->collocateKernel =
          mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                   "collocate",
                                   kernelInfo);
          elliptic->dotMultiplyKernel =
            mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                     "dotMultiply",
                                     kernelInfo);
          elliptic->dotMultiplyPfloatKernel =
            mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                     "dotMultiply",
                                     pfloatKernelInfo);

        if(serial)
          elliptic->dotMultiplyKernel =
            mesh->device.buildKernel(DHOLMES "/okl/serialDotMultiply.c",
                                     "dotMultiply",
                                     kernelInfoNoOKL);
        else
          elliptic->dotMultiplyKernel =
            mesh->device.buildKernel(DHOLMES "/okl/dotMultiply.okl",
                                     "dotMultiply",
                                     kernelInfo);
      }
    }
    MPI_Barrier(mesh->comm);
  }

  // add custom defines
  kernelInfo["defines/" "p_Nverts"] = mesh->Nverts;

  //sizes for the coarsen and prolongation kernels. degree N to degree 1
  kernelInfo["defines/" "p_NpFine"] = mesh->Np;
  kernelInfo["defines/" "p_NpCoarse"] = mesh->Nverts;

  if (elliptic->elementType == QUADRILATERALS || elliptic->elementType == HEXAHEDRA) {
    kernelInfo["defines/" "p_NqFine"] = mesh->N + 1;
    kernelInfo["defines/" "p_NqCoarse"] = 2;
  }

  int Nmax = mymax(mesh->Np, mesh->Nfaces * mesh->Nfp);
  kernelInfo["defines/" "p_Nmax"] = Nmax;

  int maxNodes = mymax(mesh->Np, (mesh->Nfp * mesh->Nfaces));
  kernelInfo["defines/" "p_maxNodes"] = maxNodes;

  int NblockV = mymax(1,maxNthreads / mesh->Np); // works for CUDA
  int NnodesV = 1; //hard coded for now
  kernelInfo["defines/" "p_NblockV"] = NblockV;
  kernelInfo["defines/" "p_NnodesV"] = NnodesV;
  kernelInfo["defines/" "p_NblockVFine"] = NblockV;
  kernelInfo["defines/" "p_NblockVCoarse"] = NblockV;

  int NblockS = mymax(1,maxNthreads / maxNodes); // works for CUDA
  kernelInfo["defines/" "p_NblockS"] = NblockS;

  int NblockP = mymax(1,maxNthreads / (4 * mesh->Np)); // get close to maxNthreads threads
  kernelInfo["defines/" "p_NblockP"] = NblockP;

  int NblockG;
  if(mesh->Np <= 32) NblockG = ( 32 / mesh->Np );
  else NblockG = maxNthreads / mesh->Np;
  kernelInfo["defines/" "p_NblockG"] = NblockG;

  kernelInfo["defines/" "p_halfC"] = (int)((mesh->cubNq + 1) / 2);
  kernelInfo["defines/" "p_halfN"] = (int)((mesh->Nq + 1) / 2);

  kernelInfo["defines/" "p_NthreadsUpdatePCG"] = (int) NthreadsUpdatePCG; // WARNING SHOULD BE MULTIPLE OF 32
  kernelInfo["defines/" "p_NwarpsUpdatePCG"] = (int) (NthreadsUpdatePCG / 32); // WARNING: CUDA SPECIFIC

  occa::properties dfloatKernelInfo = kernelInfo;
  occa::properties floatKernelInfo = kernelInfo;
  floatKernelInfo["defines/" "pfloat"] = pfloatString;
  floatKernelInfo["defines/" "dfloat"] = pfloatString;
  dfloatKernelInfo["defines/" "pfloat"] = dfloatString;

  occa::properties AxKernelInfo = dfloatKernelInfo;
  occa::properties dfloatKernelInfoNoOKL = kernelInfoNoOKL;
  dfloatKernelInfoNoOKL["defines/" "pfloat"] = dfloatString;
  if(serial) AxKernelInfo = dfloatKernelInfoNoOKL;

  for (int r = 0; r < 2; r++) {
    if ((r == 0 && mesh->rank == 0) || (r == 1 && mesh->rank > 0)) {

      if(elliptic->var_coeff) {
        sprintf(fileName,  DELLIPTIC "/okl/ellipticBuildDiagonal%s.okl", suffix);
        if(elliptic->blockSolver)
          sprintf(kernelName, "ellipticBlockBuildDiagonal%s", suffix);
        else
          sprintf(kernelName, "ellipticBuildDiagonal%s", suffix);
        elliptic->updateDiagonalKernel = mesh->device.buildKernel(fileName,
                                                                  kernelName,
                                                                  dfloatKernelInfo);
      }

      if(elliptic->blockSolver) {
        sprintf(fileName,  DELLIPTIC "/okl/ellipticBlockAx%s.okl", suffix);
        if(serial) sprintf(fileName,  DELLIPTIC "/okl/ellipticSerialAx%s.c", suffix);
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA){
          if(elliptic->stressForm){
            sprintf(kernelName, "ellipticStressAxVar%s", suffix);
          } else {
            sprintf(kernelName, "ellipticBlockAxVar%s_N%d", suffix, elliptic->Nfields);
          }
        }
        else{
          if(elliptic->stressForm){
            sprintf(kernelName, "ellipticStressAx%s", suffix);
          } else {
            sprintf(kernelName, "ellipticBlockAx%s_N%d", suffix, elliptic->Nfields);
          }
        }
      }else{
        sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
        if(serial) sprintf(fileName,  DELLIPTIC "/okl/ellipticSerialAx%s.c", suffix);
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA)
          sprintf(kernelName, "ellipticAxVar%s", suffix);
        else
          sprintf(kernelName, "ellipticAx%s", suffix);
      }
      elliptic->AxStressKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
      if(elliptic->blockSolver) {
        sprintf(fileName,  DELLIPTIC "/okl/ellipticBlockAx%s.okl", suffix);
        if(serial) sprintf(fileName,  DELLIPTIC "/okl/ellipticSerialAx%s.c", suffix);
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA){
            sprintf(kernelName, "ellipticBlockAxVar%s_N%d", suffix, elliptic->Nfields);
        }
        else{
            sprintf(kernelName, "ellipticBlockAx%s_N%d", suffix, elliptic->Nfields);
        }
      }else{
        sprintf(fileName,  DELLIPTIC "/okl/ellipticAx%s.okl", suffix);
        if(serial) sprintf(fileName,  DELLIPTIC "/okl/ellipticSerialAx%s.c", suffix);
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA)
          sprintf(kernelName, "ellipticAxVar%s", suffix);
        else
          sprintf(kernelName, "ellipticAx%s", suffix);
      }
      // Keep other kernel around
      elliptic->AxKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);

      if(!serial) {
        if(elliptic->elementType != HEXAHEDRA) {
          sprintf(kernelName, "ellipticPartialAx%s", suffix);
        }else {
          if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR")) {
            if(elliptic->var_coeff || elliptic->blockSolver) {
              printf(
                "ERROR: TRILINEAR form is not implemented for varibale coefficient and block solver yet \n");
              exit(-1);
            }
            sprintf(kernelName, "ellipticPartialAxTrilinear%s", suffix);
          }else {
            if(elliptic->blockSolver) {
              if(elliptic->var_coeff){
                if(elliptic->stressForm){
                  sprintf(kernelName, "ellipticStressPartialAxVar%s", suffix);
                } else {
                  sprintf(kernelName, "ellipticBlockPartialAxVar%s_N%d", suffix, elliptic->Nfields);
                }
              }
              else{
                if(elliptic->stressForm){
                  sprintf(kernelName, "ellipticStessPartialAx%s", suffix);
                } else {
                  sprintf(kernelName, "ellipticBlockPartialAx%s_N%d", suffix, elliptic->Nfields);
                }
              }
            }else {
              if(elliptic->var_coeff)
                sprintf(kernelName, "ellipticPartialAxVar%s", suffix);
              else
                sprintf(kernelName, "ellipticPartialAx%s", suffix);
            }
          }
        }
        elliptic->partialAxKernel = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
        elliptic->partialAxKernel2 = mesh->device.buildKernel(fileName,kernelName,AxKernelInfo);
      }

/*
      // only for Hex3D - cubature Ax
      if(elliptic->elementType == HEXAHEDRA && !elliptic->var_coeff && !elliptic->blockSolver) {
        sprintf(fileName,  DELLIPTIC "/okl/ellipticCubatureAx%s.okl", suffix);

        sprintf(kernelName, "ellipticCubaturePartialAx%s", suffix);
        elliptic->partialCubatureAxKernel = mesh->device.buildKernel(fileName,
                                                                     kernelName,
                                                                     dfloatKernelInfo);
      }
*/

      // combined PCG update and r.r kernel
      if(elliptic->blockSolver) {
        if(serial)
          elliptic->updatePCGKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSerialUpdatePCG.c",
                                     "ellipticUpdatePCG", dfloatKernelInfoNoOKL);
        else
          elliptic->updatePCGKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdatePCG.okl",
                                     "ellipticBlockUpdatePCG", dfloatKernelInfo);

        elliptic->update1NBPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
                                   "ellipticBlockUpdate1NBPCG", dfloatKernelInfo);

        elliptic->update2NBPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
                                   "ellipticBlockUpdate2NBPCG", dfloatKernelInfo);

        // combined update for Non-blocking flexible PCG
        elliptic->update0NBFPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
                                   "ellipticBlockUpdate0NBFPCG", dfloatKernelInfo);

        elliptic->update1NBFPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
                                   "ellipticBlockUpdate1NBFPCG", dfloatKernelInfo);
      }else{
        if(serial)
          elliptic->updatePCGKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSerialUpdatePCG.c",
                                     "ellipticUpdatePCG", dfloatKernelInfoNoOKL);
        else
          elliptic->updatePCGKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdatePCG.okl",
                                     "ellipticUpdatePCG", dfloatKernelInfo);

        // combined update for Non-blocking PCG
        elliptic->update1NBPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
                                   "ellipticUpdate1NBPCG", dfloatKernelInfo);

        elliptic->update2NBPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBPCG.okl",
                                   "ellipticUpdate2NBPCG", dfloatKernelInfo);

        // combined update for Non-blocking flexible PCG
        elliptic->update0NBFPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
                                   "ellipticUpdate0NBFPCG", dfloatKernelInfo);

        elliptic->update1NBFPCGKernel =
          mesh->device.buildKernel(DELLIPTIC "/okl/ellipticUpdateNBFPCG.okl",
                                   "ellipticUpdate1NBFPCG", dfloatKernelInfo);
      }

      if(!elliptic->blockSolver) {
        // Not implemented for Quad3D !!!!!
        if (options.compareArgs("BASIS","BERN")) {
          sprintf(fileName, DELLIPTIC "/okl/ellipticGradientBB%s.okl", suffix);
          sprintf(kernelName, "ellipticGradientBB%s", suffix);

          elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

          sprintf(kernelName, "ellipticPartialGradientBB%s", suffix);
          elliptic->partialGradientKernel =
            mesh->device.buildKernel(fileName,kernelName,kernelInfo);
/*
          sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdgBB%s.okl", suffix);
          sprintf(kernelName, "ellipticAxIpdgBB%s", suffix);
          elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

          sprintf(kernelName, "ellipticPartialAxIpdgBB%s", suffix);
          elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
*/
        } else if (options.compareArgs("BASIS","NODAL")) {
          sprintf(fileName, DELLIPTIC "/okl/ellipticGradient%s.okl", suffix);
          sprintf(kernelName, "ellipticGradient%s", suffix);

          elliptic->gradientKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

          sprintf(kernelName, "ellipticPartialGradient%s", suffix);
          elliptic->partialGradientKernel =
            mesh->device.buildKernel(fileName,kernelName,kernelInfo);
/*
          sprintf(fileName, DELLIPTIC "/okl/ellipticAxIpdg%s.okl", suffix);
          sprintf(kernelName, "ellipticAxIpdg%s", suffix);
          elliptic->ipdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

          sprintf(kernelName, "ellipticPartialAxIpdg%s", suffix);
          elliptic->partialIpdgKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);
*/
        }

        // Use the same kernel with quads for the following kenels
        if(elliptic->dim == 3) {
          if(elliptic->elementType == QUADRILATERALS)
            suffix = strdup("Quad2D");
          else if(elliptic->elementType == TRIANGLES)
            suffix = strdup("Tri2D");
        }

        sprintf(fileName, DELLIPTIC "/okl/ellipticPreconCoarsen%s.okl", suffix);
        sprintf(kernelName, "ellipticPreconCoarsen%s", suffix);
        elliptic->precon->coarsenKernel = mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticPreconProlongate%s.okl", suffix);
        sprintf(kernelName, "ellipticPreconProlongate%s", suffix);
        elliptic->precon->prolongateKernel =
          mesh->device.buildKernel(fileName,kernelName,kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticBlockJacobiPrecon.okl");
        sprintf(kernelName, "ellipticBlockJacobiPrecon");
        elliptic->precon->blockJacobiKernel = mesh->device.buildKernel(fileName,
                                                                       kernelName,
                                                                       kernelInfo);

        sprintf(kernelName, "ellipticPartialBlockJacobiPrecon");
        elliptic->precon->partialblockJacobiKernel = mesh->device.buildKernel(fileName,
                                                                              kernelName,
                                                                              kernelInfo);

        sprintf(fileName, DELLIPTIC "/okl/ellipticPatchSolver.okl");
        sprintf(kernelName, "ellipticApproxBlockJacobiSolver");
        elliptic->precon->approxBlockJacobiSolverKernel = mesh->device.buildKernel(fileName,
                                                                                   kernelName,
                                                                                   kernelInfo);

        if (   elliptic->elementType == TRIANGLES
               || elliptic->elementType == TETRAHEDRA) {
          elliptic->precon->SEMFEMInterpKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMInterp.okl",
                                     "ellipticSEMFEMInterp",
                                     kernelInfo);

          elliptic->precon->SEMFEMAnterpKernel =
            mesh->device.buildKernel(DELLIPTIC "/okl/ellipticSEMFEMAnterp.okl",
                                     "ellipticSEMFEMAnterp",
                                     kernelInfo);
        }
      }
    }

    MPI_Barrier(mesh->comm);
  }

  MPI_Barrier(mesh->comm);
  if(mesh->rank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel); fflush(stdout);

  if(elliptic->blockSolver) {
    elliptic->nullProjectBlockWeightGlobal = (dfloat*)calloc(elliptic->Nfields, sizeof(dfloat));

    for(int fld = 0; fld < elliptic->Nfields; fld++) {
      elliptic->sumBlockFieldKernel(Nlocal,
                                    fld,
                                    elliptic->Ntotal,
                                    elliptic->o_invDegree,
                                    elliptic->o_tmp);
      elliptic->o_tmp.copyTo(elliptic->tmp);

      dfloat nullProjectWeightLocal = 0;
      dfloat nullProjectWeightGlobal = 0;
      for(dlong n = 0; n < elliptic->Nblock; ++n)
        nullProjectWeightLocal += elliptic->tmp[n];

      MPI_Allreduce(&nullProjectWeightLocal,
                    &nullProjectWeightGlobal,
                    1,
                    MPI_DFLOAT,
                    MPI_SUM,
                    mesh->comm);

      elliptic->nullProjectBlockWeightGlobal[fld] = 1.0 / nullProjectWeightGlobal;
    }
  }else{
    // TW: WARNING C0 appropriate only
    mesh->sumKernel(mesh->Nelements * mesh->Np, elliptic->o_invDegree, elliptic->o_tmp);
    elliptic->o_tmp.copyTo(elliptic->tmp);

    dfloat nullProjectWeightLocal = 0;
    dfloat nullProjectWeightGlobal = 0;
    for(dlong n = 0; n < elliptic->Nblock; ++n)
      nullProjectWeightLocal += elliptic->tmp[n];

    MPI_Allreduce(&nullProjectWeightLocal,
                  &nullProjectWeightGlobal,
                  1,
                  MPI_DFLOAT,
                  MPI_SUM,
                  mesh->comm);

    elliptic->nullProjectWeightGlobal = 1. / nullProjectWeightGlobal;
  }

  long long int pre = mesh->device.memoryAllocated();

  oogs_mode oogsMode = OOGS_AUTO;
  if(options.compareArgs("THREAD MODEL", "SERIAL")) oogsMode = OOGS_DEFAULT;
  if(options.compareArgs("THREAD MODEL", "OPENMP")) oogsMode = OOGS_DEFAULT;
  auto callback = [&]() // hardwired to FP64 variable coeff
    {
      ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
                 elliptic->o_p, elliptic->o_Ap, dfloatString);
    };
  elliptic->oogs = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, NULL, oogsMode);
  elliptic->oogsAx = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, callback, oogsMode);

  ellipticPreconditionerSetup(elliptic, elliptic->ogs, kernelInfo);

  long long int usedBytes = mesh->device.memoryAllocated() - pre;

  elliptic->precon->preconBytes = usedBytes;

  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) {
    dlong nVecsProject = 8;
    try {
      nVecsProject = static_cast < dlong > (std::stoi(options.getArgs(
                                                        "RESIDUAL PROJECTION VECTORS")));
    } catch(std::invalid_argument& e) {
      if(elliptic->mesh->rank == 0) {
        std::cout << "Encountered invalid argument when getting RESIDUAL PROJECTION VECTORS!\n";
        std::cout << e.what();
      }
      exit(-1);
    } catch (std::out_of_range& e) {
      if(elliptic->mesh->rank == 0) {
        std::cout << "Encountered out_of_range error when getting RESIDUAL PROJECTION VECTORS!\n";
        std::cout << e.what();
      }
      exit(-1);
    }
    dlong nStepsStart = 5;
    try {
      nStepsStart = static_cast < dlong > (std::stoi(options.getArgs("RESIDUAL PROJECTION START")));
    } catch(std::invalid_argument& e) {
      if(elliptic->mesh->rank == 0) {
        std::cout << "Encountered invalid argument when getting RESIDUAL PROJECTION START!\n";
        std::cout << e.what();
      }
      exit(-1);
    } catch (std::out_of_range& e) {
      if(elliptic->mesh->rank == 0) {
        std::cout << "Encountered out_of_range error when getting RESIDUAL PROJECTION START!\n";
        std::cout << e.what();
      }
      exit(-1);
    }
    elliptic->residualProjection = new ResidualProjection(* elliptic, nVecsProject, nStepsStart);
  }

}
