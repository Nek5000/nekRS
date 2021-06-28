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
#include "platform.hpp"
#include "linAlg.hpp"

void ellipticSolveSetup(elliptic_t* elliptic, occa::properties kernelInfo)
{
  
  mesh_t* mesh      = elliptic->mesh;
  
  setupAide options = elliptic->options;

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();

  const dlong Nlocal = mesh->Np * mesh->Nelements;
  elliptic->resNormFactor = 1 / (elliptic->Nfields * mesh->volume);

  const int serial = platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP";

  if (elliptic->blockSolver &&  elliptic->elementType != HEXAHEDRA &&
      !options.compareArgs("DISCRETIZATION",
                           "CONTINUOUS") && !options.compareArgs("PRECONDITIONER","JACOBI") ) {
    if(platform->comm.mpiRank == 0)
      printf("ERROR: Block solver is implemented for C0-HEXAHEDRA with Jacobi preconditioner only\n");

    ABORT(EXIT_FAILURE);
  }

  if (options.compareArgs("COEFFICIENT","VARIABLE") &&  elliptic->elementType != HEXAHEDRA &&
      !options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    if(platform->comm.mpiRank == 0)
      printf("ERROR: Varibale coefficient solver is implemented for C0-HEXAHEDRA only\n");

    ABORT(EXIT_FAILURE);
  }

  if (options.compareArgs("COEFFICIENT","VARIABLE")) {
    if(options.compareArgs("PRECONDITIONER",
                           "MULTIGRID") &&
       !options.compareArgs("MULTIGRID VARIABLE COEFFICIENT", "FALSE")) {
      if(platform->comm.mpiRank == 0)
        printf(
          "ERROR: Varibale coefficient solver is implemented for constant multigrid preconditioner only\n");

      ABORT(EXIT_FAILURE);
    }
  }

  if(options.compareArgs("KRYLOV SOLVER", "PGMRES")){
    initializeGmresData(elliptic);
    string install_dir;
    install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
    const string oklpath = install_dir + "/okl/elliptic/";
    string filename;

    occa::properties gmresKernelInfo = platform->kernelInfo;
    gmresKernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;
    gmresKernelInfo["defines/" "p_Nfields"] = elliptic->Nfields;
    if(serial) gmresKernelInfo["okl/enabled"] = false;
    filename = serial ? oklpath + "ellipticGramSchmidtOrthogonalization.c" : oklpath + "ellipticGramSchmidtOrthogonalization.okl";
    elliptic->gramSchmidtOrthogonalizationKernel =
      platform->device.buildKernel(filename,
                               "gramSchmidtOrthogonalization",
                               gmresKernelInfo);
    filename = serial ? oklpath + "ellipticUpdatePGMRES.c" : oklpath + "ellipticUpdatePGMRES.okl";
    elliptic->updatePGMRESSolutionKernel =
      platform->device.buildKernel(filename,
                               "updatePGMRESSolution",
                               gmresKernelInfo);
    filename = serial ? oklpath + "ellipticFusedResidualAndNorm.c" : oklpath + "ellipticFusedResidualAndNorm.okl";
    elliptic->fusedResidualAndNormKernel =
      platform->device.buildKernel(filename,
                               "fusedResidualAndNorm",
                               gmresKernelInfo);
  }

  elliptic->p    = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));
  elliptic->z    = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));
  elliptic->Ap   = (dfloat*) calloc(elliptic->Ntotal * elliptic->Nfields,   sizeof(dfloat));

  elliptic->o_p    = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                         elliptic->p);
  elliptic->o_z    = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                         elliptic->z);
  elliptic->o_Ap   = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                         elliptic->Ap);
  elliptic->o_rtmp = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat),
                                         elliptic->p);

  elliptic->o_x0 = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields ,  sizeof(dfloat));

  dlong Nblocks = (Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
  elliptic->tmpNormr = (dfloat*) calloc(Nblocks,sizeof(dfloat));
  elliptic->o_tmpNormr = platform->device.malloc(Nblocks * sizeof(dfloat),
                                             elliptic->tmpNormr);

  int useFlexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");

  elliptic->type = strdup(dfloatString);

  // count total number of elements
  hlong NelementsLocal = mesh->Nelements;
  hlong NelementsGlobal = 0;

  MPI_Allreduce(&NelementsLocal, &NelementsGlobal, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);

  elliptic->NelementsGlobal = NelementsGlobal;

  elliptic->allNeumannPenalty = 1.;
  hlong localElements = (hlong) mesh->Nelements;
  hlong totalElements = 0;
  MPI_Allreduce(&localElements, &totalElements, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
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
    MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm);
    elliptic->allBlockNeumann[fld] = (gallNeumann > 0) ? 0: 1;
    // even if there is a single allNeumann activate Null space correction
    if(elliptic->allBlockNeumann[fld])
      elliptic->allNeumann = 1;
  }

  if(platform->comm.mpiRank == 0)
    printf("allNeumann = %d \n", elliptic->allNeumann);

  //copy boundary flags
  elliptic->o_EToB = platform->device.malloc(
    mesh->Nelements * mesh->Nfaces * elliptic->Nfields * sizeof(int),
    elliptic->EToB);

  //setup an unmasked gs handle
  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  if(mesh->ogs == NULL) meshParallelGatherScatterSetup(mesh, Nlocal, mesh->globalIds, platform->comm.mpiComm, verbose);

  //make a node-wise bc flag using the gsop (prioritize Dirichlet boundaries over Neumann)
  const int mapSize = elliptic->blockSolver ? elliptic->Ntotal * elliptic->Nfields: Nlocal;
  elliptic->mapB = (int*) calloc(mapSize,sizeof(int));
  const int largeNumber = 1 << 20;
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
            elliptic->mapB[fid + e * mesh->Np + fld * elliptic->Ntotal] = 
              mymin(BCFlag, elliptic->mapB[fid + e * mesh->Np + fld * elliptic->Ntotal]);
          }
        }
      }
    }
  ogsGatherScatterMany(elliptic->mapB,
                       elliptic->Nfields,
                       elliptic->Ntotal,
                       ogsInt,
                       ogsMin,
                       mesh->ogs);

  // Create mask Ids
  elliptic->Nmasked  = 0;
  elliptic->fNmasked = (dlong*)calloc(elliptic->Nfields, sizeof(dlong));
  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++) {
      if (elliptic->mapB[n + fld * elliptic->Ntotal] == largeNumber) {
        elliptic->mapB[n + fld * elliptic->Ntotal] = 0.;
      } else if (elliptic->mapB[n + fld * elliptic->Ntotal] == 1) { //Dirichlet boundary
        elliptic->Nmasked++; // increase global accumulator
        elliptic->fNmasked[fld]++; // increase local accumulator
      }
    }
  elliptic->o_mapB = platform->device.malloc(mapSize * sizeof(int), elliptic->mapB);

  elliptic->maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
  elliptic->Nmasked = 0;
  for(int fld = 0; fld < elliptic->Nfields; fld++)
    for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
      if (elliptic->mapB[n + fld * elliptic->Ntotal] == 1)
        elliptic->maskIds[elliptic->Nmasked++] = n + fld * elliptic->Ntotal;
  if (elliptic->Nmasked) 
    elliptic->o_maskIds = platform->device.malloc(elliptic->Nmasked * sizeof(dlong), elliptic->maskIds);

  if(elliptic->blockSolver) {
    elliptic->ogs = mesh->ogs; // cannot use masked version as mixed BC's possible in each field
  } else {
    hlong* maskedGlobalIds = (hlong*) calloc(Nlocal,sizeof(hlong));
    memcpy(maskedGlobalIds, mesh->globalIds, Nlocal * sizeof(hlong));
    for (dlong n = 0; n < elliptic->Nmasked; n++)
      maskedGlobalIds[elliptic->maskIds[n]] = 0;

    elliptic->ogs = ogsSetup(Nlocal, maskedGlobalIds, platform->comm.mpiComm, verbose, platform->device);
    free(maskedGlobalIds);
  }
  elliptic->o_invDegree = elliptic->ogs->o_invDegree;

  elliptic->precon = new precon_t();

  string suffix = "Hex3D";
  string filename, kernelName;

  kernelInfo["defines/pfloat"] = pfloatString;
  kernelInfo["defines/" "p_eNfields"] = elliptic->Nfields;

  occa::properties pfloatKernelInfo = kernelInfo;
  pfloatKernelInfo["defines/dfloat"] = pfloatString;
  pfloatKernelInfo["defines/pfloat"] = pfloatString;

  occa::properties kernelInfoNoOKL = kernelInfo;
  if(serial) kernelInfoNoOKL["okl/enabled"] = false;

  string install_dir;
  install_dir.assign(getenv("NEKRS_INSTALL_DIR"));

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0) printf("loading elliptic kernels ... ");
  fflush(stdout);

  {
      const string oklpath = install_dir + "/okl/core/";
      string filename;

      filename = oklpath + "mask.okl";
      mesh->maskKernel =
        platform->device.buildKernel(filename,
                                 "mask",
                                 kernelInfo);

      filename = oklpath + "mask.okl";
      mesh->maskPfloatKernel =
        platform->device.buildKernel(filename,
                                 "mask",
                                 pfloatKernelInfo);
        filename = install_dir + "/okl/elliptic/ellipticLinAlg.okl";
        elliptic->fusedCopyDfloatToPfloatKernel =
          platform->device.buildKernel(filename,
                                   "fusedCopyDfloatToPfloat",
                                   kernelInfo);
        elliptic->copyDfloatToPfloatKernel =
          platform->device.buildKernel(filename,
                                   "copyDfloatToPfloat",
                                   kernelInfo);

        elliptic->copyPfloatToDPfloatKernel =
          platform->device.buildKernel(filename,
                                   "copyPfloatToDfloat",
                                   kernelInfo);

        elliptic->scaledAddPfloatKernel =
          platform->device.buildKernel(filename,
                                   "scaledAdd",
                                   kernelInfo);
        elliptic->dotMultiplyPfloatKernel =
          platform->device.buildKernel(filename,
                                   "dotMultiply",
                                   kernelInfo);
        filename = install_dir + "/okl/elliptic/chebyshev.okl";
        elliptic->updateSmoothedSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateSmoothedSolutionVec",
                                   kernelInfo);
        elliptic->updateChebyshevSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateChebyshevSolutionVec",
                                   kernelInfo);

        elliptic->updateIntermediateSolutionVecKernel =
          platform->device.buildKernel(filename,
                                   "updateIntermediateSolutionVec",
                                   kernelInfo);

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

  kernelInfo["defines/" "p_halfC"] = (int)((mesh->cubNq + 1) / 2);
  kernelInfo["defines/" "p_halfN"] = (int)((mesh->Nq + 1) / 2);

  occa::properties dfloatKernelInfo = kernelInfo;
  occa::properties floatKernelInfo = kernelInfo;
  floatKernelInfo["defines/" "pfloat"] = pfloatString;
  floatKernelInfo["defines/" "dfloat"] = pfloatString;
  dfloatKernelInfo["defines/" "pfloat"] = dfloatString;

  occa::properties AxKernelInfo = dfloatKernelInfo;
  occa::properties dfloatKernelInfoNoOKL = kernelInfoNoOKL;
  dfloatKernelInfoNoOKL["defines/" "pfloat"] = dfloatString;
  if(serial) AxKernelInfo = dfloatKernelInfoNoOKL;

  {
      const string oklpath = install_dir + "/okl/elliptic/";
      string filename;

      filename = oklpath + "ellipticBuildDiagonal" + suffix + ".okl";
      kernelName = "ellipticBlockBuildDiagonal" + suffix;
      elliptic->updateDiagonalKernel = platform->device.buildKernel(filename,
                                                                    kernelName,
                                                                    dfloatKernelInfo);
      if(elliptic->blockSolver) {
        filename =  oklpath + "ellipticBlockAx" + suffix + ".okl";
        if(serial) filename = oklpath + "ellipticSerialAx" +  suffix + ".c";
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA) {
          if(elliptic->stressForm)
            kernelName = "ellipticStressAxVar" + suffix;
          else
            kernelName = "ellipticBlockAxVar" + suffix + "_N" + std::to_string(elliptic->Nfields);
        }else {
          if(elliptic->stressForm)
            kernelName = "ellipticStressAx" + suffix;
          else
            kernelName = "ellipticBlockAx", suffix + "_N" + std::to_string(elliptic->Nfields);
        }
      }else{
        filename = oklpath + "ellipticAx" + suffix + ".okl";
        if(serial) filename = oklpath + "ellipticSerialAx" + suffix + ".c";
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA)
          kernelName = "ellipticAxVar" + suffix;
        else
          kernelName =  "ellipticAx" + suffix;
      }
      elliptic->AxStressKernel = platform->device.buildKernel(filename,kernelName,AxKernelInfo);
      if(elliptic->blockSolver) {
        filename = oklpath + "ellipticBlockAx" + suffix + ".okl";
        if(serial) filename = oklpath + "ellipticSerialAx" + suffix + ".c";
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA)
          kernelName = "ellipticBlockAxVar" + suffix + "_N" + std::to_string(elliptic->Nfields);
        else
          kernelName = "ellipticBlockAx" + suffix + "_N" + std::to_string(elliptic->Nfields);
      }else{
        filename = oklpath + "ellipticAx" + suffix + ".okl";
        if(serial) filename = oklpath + "ellipticSerialAx" + suffix + ".c";
        if(elliptic->var_coeff && elliptic->elementType == HEXAHEDRA)
          kernelName = "ellipticAxVar" + suffix;
        else
          kernelName = "ellipticAx" + suffix;
      }
      // Keep other kernel around
      elliptic->AxKernel = platform->device.buildKernel(filename,kernelName,AxKernelInfo);

      if(!serial) {
        if(elliptic->elementType != HEXAHEDRA) {
          kernelName = "ellipticPartialAx" + suffix;
        }else {
          if(elliptic->options.compareArgs("ELEMENT MAP", "TRILINEAR")) {
            if(elliptic->var_coeff || elliptic->blockSolver) {
              printf(
                "ERROR: TRILINEAR form is not implemented for varibale coefficient and block solver yet \n");
              ABORT(EXIT_FAILURE);
            }
            kernelName = "ellipticPartialAxTrilinear" + suffix;
          }else {
            if(elliptic->blockSolver) {
              if(elliptic->var_coeff) {
                if(elliptic->stressForm)
                  kernelName = "ellipticStressPartialAxVar" + suffix;
                else
                  kernelName = "ellipticBlockPartialAxVar" + suffix + "_N" + std::to_string(elliptic->Nfields);
              }else {
                if(elliptic->stressForm)
                  kernelName = "ellipticStessPartialAx" + suffix;
                else
                  kernelName = "ellipticBlockPartialAx" + suffix + "_N" + std::to_string(elliptic->Nfields);
              }
            }else {
              if(elliptic->var_coeff)
                kernelName = "ellipticPartialAxVar" + suffix;
              else
                kernelName = "ellipticPartialAx" + suffix;
            }
          }
        }
        elliptic->partialAxKernel = platform->device.buildKernel(filename,kernelName,AxKernelInfo);
        elliptic->partialAxKernel2 = platform->device.buildKernel(filename,kernelName,AxKernelInfo);
      }

      // combined PCG update and r.r kernel
      if(serial) {
        filename = oklpath + "ellipticSerialUpdatePCG.c";
        elliptic->updatePCGKernel =
          platform->device.buildKernel(filename,
                                   "ellipticUpdatePCG", dfloatKernelInfoNoOKL);
      } else {
        filename = oklpath + "ellipticUpdatePCG.okl";
        elliptic->updatePCGKernel =
          platform->device.buildKernel(filename,
                                   "ellipticBlockUpdatePCG", dfloatKernelInfo);
      }

      if(!elliptic->blockSolver) {
        if(serial){
          filename = oklpath + "ellipticPreconCoarsen" + suffix + ".c";
          kernelName = "ellipticPreconCoarsen" + suffix;
          occa::properties serialPreconKernelInfo = kernelInfo;
          serialPreconKernelInfo["okl/enabled"] = false;
          elliptic->precon->coarsenKernel = platform->device.buildKernel(filename,kernelName,serialPreconKernelInfo);
          filename = oklpath + "ellipticPreconProlongate" + suffix + ".c";
          kernelName = "ellipticPreconProlongate" + suffix;
          elliptic->precon->prolongateKernel =
            platform->device.buildKernel(filename,kernelName,serialPreconKernelInfo);
        } else {
          filename = oklpath + "ellipticPreconCoarsen" + suffix + ".okl";
          kernelName = "ellipticPreconCoarsen" + suffix;
          elliptic->precon->coarsenKernel = platform->device.buildKernel(filename,kernelName,kernelInfo);
          filename = oklpath + "ellipticPreconProlongate" + suffix + ".okl";
          kernelName = "ellipticPreconProlongate" + suffix;
          elliptic->precon->prolongateKernel =
            platform->device.buildKernel(filename,kernelName,kernelInfo);
        }

      }
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel);
  fflush(stdout);

  oogs_mode oogsMode = OOGS_AUTO;
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  auto callback = [&]() // hardwired to FP64 variable coeff
                  {
                    ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
                               elliptic->o_p, elliptic->o_Ap, dfloatString);
                  };
  elliptic->oogs = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, NULL, oogsMode);
  elliptic->oogsAx = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, callback, oogsMode);

  long long int pre = platform->device.memoryAllocated();
  ellipticPreconditionerSetup(elliptic, elliptic->ogs, kernelInfo);

  long long int usedBytes = platform->device.memoryAllocated() - pre;

  elliptic->precon->preconBytes = usedBytes;

  if(options.compareArgs("RESIDUAL PROJECTION","TRUE")) {
    dlong nVecsProject = 8;
    options.getArgs("RESIDUAL PROJECTION VECTORS", nVecsProject);

    dlong nStepsStart = 5;
    options.getArgs("RESIDUAL PROJECTION START", nStepsStart);

    ResidualProjection::ProjectionType type = ResidualProjection::ProjectionType::CLASSIC;
    if(options.compareArgs("RESIDUAL PROJECTION METHOD", "CLASSIC"))
      type = ResidualProjection::ProjectionType::CLASSIC;
    else if (options.compareArgs("RESIDUAL PROJECTION METHOD", "ACONJ"))
      type = ResidualProjection::ProjectionType::ACONJ;
    else {
      if(platform->comm.mpiRank == 0){
        printf("Encountered invalid residual projection method %s!\n",
          options.getArgs("RESIDUAL PROJECTION METHOD").c_str());
      }
      exit(1);
    }

    elliptic->residualProjection = new ResidualProjection(*elliptic, type, nVecsProject, nStepsStart);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
  fflush(stdout);
}
