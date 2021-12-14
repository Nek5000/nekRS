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

void ellipticSolveSetup(elliptic_t* elliptic)
{
  
  mesh_t* mesh      = elliptic->mesh;
  
  setupAide options = elliptic->options;

  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();

  const dlong Nlocal = mesh->Np * mesh->Nelements;
  elliptic->resNormFactor = 1 / (elliptic->Nfields * mesh->volume);

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
    const std::string sectionIdentifier = std::to_string(elliptic->Nfields) + "-";
    elliptic->gramSchmidtOrthogonalizationKernel =
      platform->kernels.get(sectionIdentifier + "gramSchmidtOrthogonalization");
    elliptic->updatePGMRESSolutionKernel =
      platform->kernels.get(sectionIdentifier + "updatePGMRESSolution");
    elliptic->fusedResidualAndNormKernel =
      platform->kernels.get(sectionIdentifier + "fusedResidualAndNorm");
  }

  const size_t offsetBytes = elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat);
  if(elliptic->o_wrk.size() < elliptic_t::NScratchFields * offsetBytes) {
    if(platform->comm.mpiRank == 0) printf("ERROR: mempool assigned for elliptic too small!");
    ABORT(EXIT_FAILURE);
  }

#if 0  
  elliptic->o_p    = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
  elliptic->o_z    = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
  elliptic->o_Ap   = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
  elliptic->o_x0   = platform->device.malloc(elliptic->Ntotal * elliptic->Nfields * sizeof(dfloat));
#else
  elliptic->o_p    = elliptic->o_wrk + 0*offsetBytes;
  elliptic->o_z    = elliptic->o_wrk + 1*offsetBytes; 
  elliptic->o_Ap   = elliptic->o_wrk + 2*offsetBytes; 
  elliptic->o_x0   = elliptic->o_wrk + 3*offsetBytes; 
#endif

  const dlong Nblocks = (Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
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

  elliptic->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces * elliptic->Nfields,sizeof(int));
  int* allNeumann = (int*)calloc(elliptic->Nfields, sizeof(int));

  dfloat* lambda = (dfloat*) calloc(2*elliptic->Ntotal, sizeof(dfloat));
  elliptic->o_lambda.copyTo(lambda, 2*elliptic->Ntotal*sizeof(dfloat));
  // check based on the coefficient
  for(int fld = 0; fld < elliptic->Nfields; fld++) {
    if(elliptic->coeffField) {
      int allzero = 1;
      for(int n = 0; n < Nlocal; n++) { // check any non-zero value for each field
        if(lambda[n + elliptic->Ntotal + fld * elliptic->loffset]) {
          allzero = 0;
          break;
        }
      }
      allNeumann[fld] = allzero;
    }else{
      allNeumann[fld] = (lambda[elliptic->Ntotal + fld * elliptic->loffset] == 0) ? 1 : 0;
    }
  }

  free(lambda);

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
  int* allBlockNeumann = (int*)calloc(elliptic->Nfields, sizeof(int));
  for(int fld = 0; fld < elliptic->Nfields; fld++) {
    int lallNeumann, gallNeumann;
    lallNeumann = allNeumann[fld] ? 0:1;
    MPI_Allreduce(&lallNeumann, &gallNeumann, 1, MPI_INT, MPI_SUM, platform->comm.mpiComm);
    allBlockNeumann[fld] = (gallNeumann > 0) ? 0: 1;
    // even if there is a single allNeumann activate Null space correction
    if(allBlockNeumann[fld])
      elliptic->allNeumann = 1;
  }
  free(allBlockNeumann);

  if(platform->comm.mpiRank == 0)
    printf("allNeumann = %d \n", elliptic->allNeumann);

  int verbose = options.compareArgs("VERBOSE","TRUE") ? 1:0;
  if(mesh->ogs == NULL) meshParallelGatherScatterSetup(mesh, Nlocal, mesh->globalIds, platform->comm.mpiComm, verbose);

  { //setup an unmasked gs handle
    ogs_t *ogs = NULL;
    if(elliptic->blockSolver) ogs = mesh->ogs; 
    ellipticOgs(mesh, elliptic->Ntotal, elliptic->Nfields, elliptic->Ntotal, elliptic->BCType, elliptic->NBCType, 
                elliptic->Nmasked, elliptic->o_mapB, elliptic->o_maskIds, &ogs);
    elliptic->ogs = ogs;
    elliptic->o_invDegree = elliptic->ogs->o_invDegree;
  }

  elliptic->precon = new precon_t();

  std::string suffix = "Hex3D";
  std::string kernelName;

  MPI_Barrier(platform->comm.mpiComm);
  double tStartLoadKernel = MPI_Wtime();
  if(platform->comm.mpiRank == 0) printf("loading elliptic kernels ... ");
  fflush(stdout);

  {
      mesh->maskKernel =
        platform->kernels.get("mask");
  }

  {
      const std::string sectionIdentifier = std::to_string(elliptic->Nfields) + "-";
      kernelName = "ellipticBlockBuildDiagonal" + suffix;
      elliptic->ellipticBlockBuildDiagonalKernel = platform->kernels.get(sectionIdentifier + kernelName);
      elliptic->axmyzManyPfloatKernel = platform->kernels.get("axmyzManyPfloat");
      elliptic->adyManyPfloatKernel = platform->kernels.get("adyManyPfloat");

      std::string kernelNamePrefix = "";
      if(elliptic->poisson) kernelNamePrefix += "poisson-";
      kernelNamePrefix += "elliptic";
      if (elliptic->blockSolver)
        kernelNamePrefix += (elliptic->stressForm) ? "Stress" : "Block";
 
      kernelName = "Ax";
      if (elliptic->coeffField) kernelName += "Coeff";
      if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR")) kernelName += "Trilinear";
      kernelName += suffix; 
      if (elliptic->blockSolver && !elliptic->stressForm) 
        kernelName += "_N" + std::to_string(elliptic->Nfields);

      elliptic->AxKernel = 
        platform->kernels.get(kernelNamePrefix + "Partial" + kernelName);

      elliptic->updatePCGKernel =
        platform->kernels.get(sectionIdentifier + "ellipticBlockUpdatePCG");
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStartLoadKernel);
  fflush(stdout);

  oogs_mode oogsMode = OOGS_AUTO;
  auto callback = [&]() // hardwired to FP64 variable coeff
                  {
                    ellipticAx(elliptic, mesh->NlocalGatherElements, mesh->o_localGatherElementList,
                               elliptic->o_p, elliptic->o_Ap, dfloatString);
                  };
  elliptic->oogs = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, NULL, oogsMode);
  elliptic->oogsAx = oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, callback, oogsMode);

  long long int pre = platform->device.occaDevice().memoryAllocated();
  ellipticPreconditionerSetup(elliptic, elliptic->ogs);

  long long int usedBytes = platform->device.occaDevice().memoryAllocated() - pre;

  elliptic->precon->preconBytes = usedBytes;

  if(options.compareArgs("INITIAL GUESS","PROJECTION") ||
     options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ"))
  {
    dlong nVecsProject = 8;
    options.getArgs("RESIDUAL PROJECTION VECTORS", nVecsProject);

    dlong nStepsStart = 5;
    options.getArgs("RESIDUAL PROJECTION START", nStepsStart);

    SolutionProjection::ProjectionType type = SolutionProjection::ProjectionType::CLASSIC;
    if(options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ"))
      type = SolutionProjection::ProjectionType::ACONJ;
    else if (options.compareArgs("INITIAL GUESS", "PROJECTION"))
      type = SolutionProjection::ProjectionType::CLASSIC;

    elliptic->solutionProjection = new SolutionProjection(*elliptic, type, nVecsProject, nStepsStart);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
  fflush(stdout);
}
