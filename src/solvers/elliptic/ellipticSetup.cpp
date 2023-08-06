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

occa::memory elliptic_t::o_wrk = occa::memory();

void checkConfig(elliptic_t *elliptic)
{
  mesh_t *mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  int err = 0;

  if (!options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
    if (platform->comm.mpiRank == 0)
      printf("solver only supports CG\n");
    err++;
  }

  if (elliptic->elementType != HEXAHEDRA) {
    if (platform->comm.mpiRank == 0)
      printf("solver only supports HEX elements\n");
    err++;
  }

  if (elliptic->blockSolver && options.compareArgs("PRECONDITIONER", "MULTIGRID")) {
    if (platform->comm.mpiRank == 0)
      printf("Block solver does not support multigrid preconditioner\n");
    err++;
  }

  if (!elliptic->poisson && options.compareArgs("PRECONDITIONER", "MULTIGRID") &&
      !options.compareArgs("MULTIGRID SMOOTHER", "DAMPEDJACOBI")) {
    if (platform->comm.mpiRank == 0)
      printf("Non-Poisson type equations require Jacobi multigrid smoother\n");
    err++;
  }

  if (options.compareArgs("PRECONDITIONER", "MULTIGRID") &&
      options.compareArgs("MULTIGRID COARSE SOLVE", "TRUE")) {
    if (elliptic->poisson == 0) {
      if (platform->comm.mpiRank == 0)
        printf("Multigrid + coarse solve only supported for Poisson type equations\n");
      err++;
    }
  }

  if (elliptic->mesh->ogs == NULL) {
    if (platform->comm.mpiRank == 0)
      printf("mesh->ogs == NULL!");
    err++;
  }

  if (elliptic->Nfields < 1 || elliptic->Nfields > 3) {
    if (platform->comm.mpiRank == 0)
      printf("Invalid Nfields = %d!", elliptic->Nfields);
    err++;
  }

  nrsCheck(elliptic->EToB == nullptr,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s",
           "elliptic->EToB not allocated!\n");

  {
    int found = 0;
    for (int fld = 0; fld < elliptic->Nfields; fld++) {
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int offset = fld * mesh->Nelements * mesh->Nfaces;
          const int bc = elliptic->EToB[f + e * mesh->Nfaces + offset];
          if (bc == ZERO_NORMAL || bc == ZERO_TANGENTIAL)
            found = 1;
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    if (found && options.compareArgs("BLOCK SOLVER", "FALSE")) {
      if (platform->comm.mpiRank == 0)
        printf("Unaligned BCs require block solver!\n");
      err++;
    }
  }

  nrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s", "\n");
}

void ellipticSolveSetup(elliptic_t *elliptic)
{
  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();

  nrsCheck(elliptic->name.size() == 0,
           platform->comm.mpiComm,
           EXIT_FAILURE,
           "%s\n",
           "Empty elliptic solver name!");

  elliptic->options.setArgs("DISCRETIZATION", "CONTINUOUS");

  platform->options.getArgs("ELEMENT TYPE", elliptic->elementType);

  elliptic->blockSolver = elliptic->Nfields > 1;

  // create private options based on platform
  for (auto &entry : platform->options.keyWordToDataMap) {
    std::string prefix = elliptic->name;
    upperCase(prefix);
    if (entry.first.find(prefix) == 0) {
      std::string key = entry.first;
      key.erase(0, prefix.size() + 1);
      elliptic->options.setArgs(key, entry.second);
    }
  }

  if (platform->comm.mpiRank == 0 && platform->verbose)
    std::cout << elliptic->options << std::endl;

  if (platform->device.mode() == "Serial")
    elliptic->options.setArgs("COARSE SOLVER LOCATION", "CPU");

  setupAide &options = elliptic->options;
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE") ? 1 : 0;

  mesh_t *mesh = elliptic->mesh;
  const dlong Nlocal = mesh->Np * mesh->Nelements;

  const dlong Nblocks = (Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;

  elliptic->type = strdup(dfloatString);

  hlong NelementsLocal = mesh->Nelements;
  hlong NelementsGlobal = 0;
  MPI_Allreduce(&NelementsLocal, &NelementsGlobal, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
  elliptic->NelementsGlobal = NelementsGlobal;

  elliptic->o_EToB = platform->device.malloc(mesh->Nelements * mesh->Nfaces * elliptic->Nfields * sizeof(int),
                                             elliptic->EToB);

  checkConfig(elliptic);

  if (options.compareArgs("SOLVER", "PGMRES")) {
    initializeGmresData(elliptic);
    const std::string sectionIdentifier = std::to_string(elliptic->Nfields) + "-";
    elliptic->gramSchmidtOrthogonalizationKernel =
        platform->kernels.get(sectionIdentifier + "gramSchmidtOrthogonalization");
    elliptic->updatePGMRESSolutionKernel = platform->kernels.get(sectionIdentifier + "updatePGMRESSolution");
    elliptic->fusedResidualAndNormKernel = platform->kernels.get(sectionIdentifier + "fusedResidualAndNorm");
  }

  mesh->maskKernel = platform->kernels.get("mask");
  mesh->maskPfloatKernel = platform->kernels.get("maskPfloat");
 
  ellipticUpdateWorkspace(elliptic);

  elliptic->tmpNormr = (dfloat *)calloc(Nblocks, sizeof(dfloat));
  elliptic->o_tmpNormr = platform->device.malloc(Nblocks * sizeof(dfloat), elliptic->tmpNormr);

  elliptic->allNeumann = 0;
  if (elliptic->poisson) {
    int allNeumann = 1;

    // check based on BC
    for (int fld = 0; fld < elliptic->Nfields; fld++) {
      for (dlong e = 0; e < mesh->Nelements; e++) {
        for (int f = 0; f < mesh->Nfaces; f++) {
          const int offset = fld * mesh->Nelements * mesh->Nfaces;
          const int bc = elliptic->EToB[f + e * mesh->Nfaces + offset];
          if (bc > 0 && bc != NEUMANN)
            allNeumann = 0;
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &allNeumann, 1, MPI_INT, MPI_MIN, platform->comm.mpiComm);
    elliptic->allNeumann = allNeumann;
    if (platform->comm.mpiRank == 0 && elliptic->allNeumann)
      printf("non-trivial nullSpace detected\n");
  }

  { // setup an masked gs handle
    ogs_t *ogs = (elliptic->blockSolver) ? mesh->ogs : nullptr;
    ellipticOgs(mesh,
                elliptic->fieldOffset,
                elliptic->Nfields,
                elliptic->fieldOffset,
                elliptic->EToB,
                elliptic->Nmasked,
                elliptic->o_maskIds,
                elliptic->NmaskedLocal,
                elliptic->o_maskIdsLocal,
                elliptic->NmaskedGlobal,
                elliptic->o_maskIdsGlobal,
                &ogs);
    elliptic->ogs = ogs;
    elliptic->o_invDegree = elliptic->ogs->o_invDegree;
  }

  {
    std::string kernelName;
    const std::string suffix = "Hex3D";

    const std::string sectionIdentifier = std::to_string(elliptic->Nfields) + "-";
    kernelName = "ellipticBlockBuildDiagonal" + suffix;
    elliptic->ellipticBlockBuildDiagonalKernel = platform->kernels.get(sectionIdentifier + kernelName);

    kernelName = "fusedCopyDfloatToPfloat";
    elliptic->fusedCopyDfloatToPfloatKernel = platform->kernels.get(kernelName);

    std::string kernelNamePrefix = "";
    if (elliptic->poisson)
      kernelNamePrefix += "poisson-";
    kernelNamePrefix += "elliptic";
    if (elliptic->blockSolver)
      kernelNamePrefix += (elliptic->stressForm) ? "Stress" : "Block";

    kernelName = "Ax";
    kernelName += "Coeff";
    if (platform->options.compareArgs("ELEMENT MAP", "TRILINEAR"))
      kernelName += "Trilinear";
    kernelName += suffix;

    elliptic->AxKernel = platform->kernels.get(kernelNamePrefix + "Partial" + kernelName);

    elliptic->updatePCGKernel = platform->kernels.get(sectionIdentifier + "ellipticBlockUpdatePCG");
  }

  auto timeEllipticOperator = [&]() {
    const int Nsamples = 10;
    ellipticOperator(elliptic, elliptic->o_p, elliptic->o_Ap, dfloatString);

    platform->device.finish();
    MPI_Barrier(platform->comm.mpiComm);
    const double start = MPI_Wtime();

    for (int test = 0; test < Nsamples; ++test)
      ellipticOperator(elliptic, elliptic->o_p, elliptic->o_Ap, dfloatString);

    platform->device.finish();
    double elapsed = (MPI_Wtime() - start) / Nsamples;
    MPI_Allreduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);

    return elapsed;
  };

  oogs_mode oogsMode = OOGS_AUTO;
  elliptic->oogs =
      oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->fieldOffset, ogsDfloat, NULL, oogsMode);
  elliptic->oogsAx = elliptic->oogs;

  if (platform->options.compareArgs("ENABLE GS COMM OVERLAP", "TRUE")) {
    auto nonOverlappedTime = timeEllipticOperator();
    auto callback = [&]() {
      ellipticAx(elliptic,
                 mesh->NlocalGatherElements,
                 mesh->o_localGatherElementList,
                 elliptic->o_p,
                 elliptic->o_Ap,
                 dfloatString);
    };
    elliptic->oogsAx =
        oogs::setup(elliptic->ogs, elliptic->Nfields, elliptic->fieldOffset, ogsDfloat, callback, oogsMode);

    auto overlappedTime = timeEllipticOperator();
    if (overlappedTime > nonOverlappedTime)
      elliptic->oogsAx = elliptic->oogs;

    if (platform->comm.mpiRank == 0) {
      printf("testing Ax overlap %.2es %.2es ", nonOverlappedTime, overlappedTime);
      if (elliptic->oogsAx != elliptic->oogs)
        printf("(overlap enabled)");

      printf("\n");
    }
  }

  ellipticPreconditionerSetup(elliptic, elliptic->ogs);

  if (options.compareArgs("INITIAL GUESS", "PROJECTION") ||
      options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ")) {
    dlong nVecsProject = 8;
    options.getArgs("RESIDUAL PROJECTION VECTORS", nVecsProject);

    dlong nStepsStart = 5;
    options.getArgs("RESIDUAL PROJECTION START", nStepsStart);

    SolutionProjection::ProjectionType type = SolutionProjection::ProjectionType::CLASSIC;
    if (options.compareArgs("INITIAL GUESS", "PROJECTION-ACONJ"))
      type = SolutionProjection::ProjectionType::ACONJ;
    else if (options.compareArgs("INITIAL GUESS", "PROJECTION"))
      type = SolutionProjection::ProjectionType::CLASSIC;

    elliptic->solutionProjection = new SolutionProjection(*elliptic, type, nVecsProject, nStepsStart);
  }

  MPI_Barrier(platform->comm.mpiComm);
  if (platform->comm.mpiRank == 0)
    printf("done (%gs)\n", MPI_Wtime() - tStart);
  fflush(stdout);
}

elliptic_t::~elliptic_t()
{
  if (precon)
    delete this->precon;
  free(this->tmpNormr);
  this->o_tmpNormr.free();
  this->o_EToB.free();
}
