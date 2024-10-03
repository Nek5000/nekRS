#include "platform.hpp"
#include "nekInterfaceAdapter.hpp"
#include "meshNekReader.hpp"

static void checkEToB(mesh_t *mesh)
{
  const auto nid = mesh->Nbid;

  if (nid <= 0) {
    return;
  }

  int err = 0;
  int found = 0;

  for (int id = 1; id <= nid; id++) {
    found = 0;
    for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) {
      if (mesh->EToB[f] == id) {
        found = 1;
        break;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
    err += (found ? 0 : 1);
    if (err && platform->comm.mpiRank == 0) {
      printf("Cannot find boundary ID %d in EToB!\n", id);
    }
  }
  nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "");

  found = 0;
  for (int f = 0; f < mesh->Nelements * mesh->Nfaces; f++) {
    if (mesh->EToB[f] < -1 || mesh->EToB[f] == 0 || mesh->EToB[f] > nid) {
      found = 1;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  nekrsCheck(found, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "EToB has invalid entries!");
}

static void meshVOccaSetup3D(mesh_t *mesh, occa::properties &kernelInfo);

mesh_t *createMeshV(MPI_Comm comm, int N, int cubN, mesh_t *meshT, occa::properties &kernelInfo);

occa::properties meshKernelProperties(int N)
{
  occa::properties meshProperties;
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;
  const int Nfp = Nq * Nq;
  constexpr int Nfaces{6};

  constexpr int Nvgeo{12};
  constexpr int Nggeo{7};
  constexpr int Nsgeo{13};

  nekrsCheck(BLOCKSIZE < Nq * Nq,
             MPI_COMM_SELF,
             EXIT_FAILURE,
             "BLOCKSIZE of %d < %d (Nq * Nq)\n!",
             BLOCKSIZE,
             Nq * Nq);

  meshProperties["defines/"
                 "p_dim"] = 3;
  meshProperties["defines/"
                 "p_Nverts"] = 8;
  meshProperties["defines/"
                 "p_Nfields"] = 1;
  meshProperties["defines/"
                 "p_N"] = N;
  meshProperties["defines/"
                 "p_Nq"] = Nq;
  meshProperties["defines/"
                 "p_Nq_g"] = Nq;
  meshProperties["defines/"
                 "p_Np"] = Np;
  meshProperties["defines/"
                 "p_Np_g"] = Np;
  meshProperties["defines/"
                 "p_Nfp"] = Nfp;
  meshProperties["defines/"
                 "p_Nfaces"] = Nfaces;
  meshProperties["defines/"
                 "p_NfacesNfp"] = Nfp * Nfaces;

  meshProperties["defines/"
                 "p_Nvgeo"] = Nvgeo;
  meshProperties["defines/"
                 "p_Nsgeo"] = Nsgeo;
  meshProperties["defines/"
                 "p_Nggeo"] = Nggeo;

  meshProperties["defines/"
                 "p_NXID"] = NXID;
  meshProperties["defines/"
                 "p_NYID"] = NYID;
  meshProperties["defines/"
                 "p_NZID"] = NZID;
  meshProperties["defines/"
                 "p_SJID"] = SJID;
  meshProperties["defines/"
                 "p_IJID"] = IJID;
  meshProperties["defines/"
                 "p_WIJID"] = WIJID;
  meshProperties["defines/"
                 "p_WSJID"] = WSJID;
  meshProperties["defines/"
                 "p_T1XID"] = T1XID;
  meshProperties["defines/"
                 "p_T1YID"] = T1YID;
  meshProperties["defines/"
                 "p_T1ZID"] = T1ZID;
  meshProperties["defines/"
                 "p_T2XID"] = T2XID;
  meshProperties["defines/"
                 "p_T2YID"] = T2YID;
  meshProperties["defines/"
                 "p_T2ZID"] = T2ZID;

  meshProperties["defines/"
                 "p_G00ID"] = G00ID;
  meshProperties["defines/"
                 "p_G01ID"] = G01ID;
  meshProperties["defines/"
                 "p_G02ID"] = G02ID;
  meshProperties["defines/"
                 "p_G11ID"] = G11ID;
  meshProperties["defines/"
                 "p_G12ID"] = G12ID;
  meshProperties["defines/"
                 "p_G22ID"] = G22ID;
  meshProperties["defines/"
                 "p_GWJID"] = GWJID;

  meshProperties["defines/"
                 "p_RXID"] = RXID;
  meshProperties["defines/"
                 "p_SXID"] = SXID;
  meshProperties["defines/"
                 "p_TXID"] = TXID;

  meshProperties["defines/"
                 "p_RYID"] = RYID;
  meshProperties["defines/"
                 "p_SYID"] = SYID;
  meshProperties["defines/"
                 "p_TYID"] = TYID;

  meshProperties["defines/"
                 "p_RZID"] = RZID;
  meshProperties["defines/"
                 "p_SZID"] = SZID;
  meshProperties["defines/"
                 "p_TZID"] = TZID;

  meshProperties["defines/"
                 "p_JID"] = JID;
  meshProperties["defines/"
                 "p_JWID"] = JWID;
  meshProperties["defines/"
                 "p_IJWID"] = IJWID;
  return meshProperties;
}

void meshLoadKernels(mesh_t *mesh)
{
  const std::string meshPrefix = "mesh-";
  const std::string orderSuffix = "_" + std::to_string(mesh->N);

  mesh->geometricFactorsKernel =
      platform->kernelRequests.load(meshPrefix + "geometricFactorsHex3D" + orderSuffix);
  mesh->surfaceGeometricFactorsKernel =
      platform->kernelRequests.load(meshPrefix + "surfaceGeometricFactorsHex3D" + orderSuffix);
  mesh->cubatureGeometricFactorsKernel =
      platform->kernelRequests.load(meshPrefix + "cubatureGeometricFactorsHex3D" + orderSuffix);

  mesh->setBIDKernel = platform->kernelRequests.load(meshPrefix + "setBIDHex3D");
  mesh->distanceKernel = platform->kernelRequests.load(meshPrefix + "distanceHex3D");
  mesh->hlongSumKernel = platform->kernelRequests.load("hlong-" + meshPrefix + "sum");

  mesh->velocityDirichletKernel = platform->kernelRequests.load(meshPrefix + "velocityDirichletBCHex3D");
  mesh->nStagesSumVectorKernel = platform->kernelRequests.load(meshPrefix + "nStagesSumVector");

  { 
    const auto prefix = "coarsenHex3D_Nf_" + std::to_string(mesh->N);
    for (int N = 1; N < mesh->N; N++) {
      mesh->intpKernel[N] = 
        platform->kernelRequests.load(meshPrefix + prefix + std::string("_Nc_") + std::to_string(N));
    }
  }
  { 
    for (int N = mesh->N; N < mesh->maxNqIntp; N++) {
      const auto prefix = "prolongateHex3D_Nf_" + std::to_string(N);
      mesh->intpKernel[N] = 
        platform->kernelRequests.load(meshPrefix + prefix + std::string("_Nc_") + std::to_string(mesh->N));
    }
  }
}

std::pair<mesh_t*, mesh_t*> createMesh(MPI_Comm comm, int N, int cubN, bool cht, occa::properties &kernelInfo)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh_t *mesh = new mesh_t();
  platform->options.getArgs("MESH INTEGRATION ORDER", mesh->nAB);
  mesh->cht = cht;

  if (platform->comm.mpiRank == 0) {
    if (mesh->cht) {
      printf("generating t-mesh ...\n");
    } else {
      printf("generating mesh ...\n");
    }
  }

  meshNekReaderHex3D(N, mesh);

  nekrsCheck(static_cast<size_t>(mesh->Nelements) * mesh->Nvgeo * cubN > std::numeric_limits<int>::max(),
             platform->comm.mpiComm,
             EXIT_FAILURE,
             "%s\n",
             "mesh->Nelements * mesh->Nvgeo * mesh->cubN exceeds int limit!");

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N, cubN);
  if (platform->comm.mpiRank == 0) {
    printf("polynomial order N: %d", mesh->N);
    if (cubN) {
      printf(", over-integration order cubN: %d", mesh->cubNq - 1);
    }
    printf("\n");
  }

  if (platform->comm.mpiRank == 0 && mesh->N < 5) {
    std::cout << std::endl << "    WARNING: N < 5 may degrade performance!\n" << std::endl;
  }

  meshLoadKernels(mesh);

  meshPhysicalNodesHex3D(mesh);

  meshConnectFaceNodes3D(mesh);

  meshOccaSetup3D(mesh, platform->options, kernelInfo);

  checkEToB(mesh);

  meshGlobalIds(mesh);

  meshParallelGatherScatterSetup(mesh,
                                 mesh->Nelements * mesh->Np,
                                 mesh->globalIds,
                                 platform->comm.mpiComm,
                                 OOGS_AUTO,
                                 0);

  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nlocal, ogsDfloat, NULL, OOGS_AUTO);

  mesh->update();

  int err = 0;
  int Nfine;
  platform->options.getArgs("POLYNOMIAL DEGREE", Nfine);
  if (mesh->N == Nfine) {
    dfloat *tmp = (dfloat *)calloc(mesh->Nlocal, sizeof(dfloat));
    mesh->ogs->o_invDegree.copyTo(tmp, mesh->Nlocal);
    double *mult = (cht) ? nek::ptr<double>("tmult") : nek::ptr<double>("vmult");
    dfloat sum1 = 0;
    for (int i = 0; i < mesh->Nlocal; i++) {
      sum1 += std::abs(tmp[i] - static_cast<dfloat>(mult[i]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    if (sum1 > 10 * std::numeric_limits<dfloat>::epsilon()) {
      if (platform->comm.mpiRank == 0) {
        printf("multiplicity test err=%g!\n", sum1);
      }
      fflush(stdout);
      err++;
    }
    free(tmp);
  }
  nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "");

  if (platform->options.compareArgs("MOVING MESH", "TRUE")) {
    const int maxTemporalOrder = 3;
    mesh->coeffAB = (dfloat *)calloc(maxTemporalOrder, sizeof(dfloat));
    mesh->o_coeffAB = platform->device.malloc<dfloat>(maxTemporalOrder, mesh->coeffAB);
  }

  {
    std::vector<dfloat> tmp(mesh->Nlocal);
    mesh->ogs->o_invDegree.copyTo(tmp.data());

    double sum = 0;
    for (int i = 0; i < mesh->Nlocal; i++) {
      sum += tmp[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, platform->comm.mpiComm);
    auto Nglobal = static_cast<hlong>(sum + 0.1);
    if (platform->comm.mpiRank == 0) {
      printf("unique number of gridpoints: : %lld\n", Nglobal);
    }
    mesh->Nglobal = Nglobal;

    hlong NelementsGlobal = mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &NelementsGlobal, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
    mesh->NelementsGlobal = NelementsGlobal;
  }

  {
    double valMin = (double)mesh->NlocalGatherElements / mesh->Nelements;
    double valMax = (double)mesh->NlocalGatherElements / mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &valMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);
    MPI_Allreduce(MPI_IN_PLACE, &valMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0 && platform->comm.mpiCommSize > 1) {
      printf("number of interior elements min/max: %2.0f%%  %2.0f%%\n", 100 * valMin, 100 * valMax);
    }
  }

  mesh_t* meshV = nullptr;
  if (mesh->cht) {
    meshV = createMeshV(comm, N, cubN, mesh, kernelInfo);
  }

  return {mesh, meshV};
}

mesh_t *createMeshMG(mesh_t *_mesh, int Nc)
{
  mesh_t *mesh = new mesh_t();
  memcpy(mesh, _mesh, sizeof(mesh_t));

  const int cubN = 0;
  meshLoadReferenceNodesHex3D(mesh, Nc, cubN);

  const std::string orderSuffix = "_" + std::to_string(mesh->N);

  int p;
  platform->options.getArgs("POLYNOMIAL DEGREE", p);
  const std::string prefix = (mesh->N != p) ? "pMGmesh-" : "mesh-";

  mesh->geometricFactorsKernel =
      platform->kernelRequests.load(prefix + "geometricFactorsHex3D" + orderSuffix);

  mesh->velocityDirichletKernel = nullptr;
  mesh->surfaceGeometricFactorsKernel = nullptr;
  mesh->cubatureGeometricFactorsKernel = nullptr;
  mesh->nStagesSumVectorKernel = nullptr;

  mesh->o_D = platform->device.malloc<dfloat>(mesh->Nq * mesh->Nq, mesh->D);

  dfloat *DT = (dfloat *)calloc(mesh->Nq * mesh->Nq, sizeof(dfloat));
  for (int j = 0; j < mesh->Nq; j++) {
    for (int i = 0; i < mesh->Nq; i++) {
      DT[j * mesh->Nq + i] = mesh->D[i * mesh->Nq + j];
    }
  }
  mesh->o_DT = platform->device.malloc<dfloat>(mesh->Nq * mesh->Nq, DT);
  free(DT);

  mesh->o_gllw = platform->device.malloc<dfloat>(mesh->Nq, mesh->gllw);

  mesh->o_faceNodes = platform->device.malloc<int>(mesh->Nfaces * mesh->Nfp, mesh->faceNodes);

  mesh->o_LMM = platform->device.malloc<dfloat>(mesh->Nlocal);

  mesh->o_vgeo = platform->device.malloc<dfloat>(mesh->Nlocal * mesh->Nvgeo);

  mesh->o_ggeo = platform->device.malloc<dfloat>(mesh->Nlocal * mesh->Nggeo);

  meshPhysicalNodesHex3D(mesh);

  meshConnectFaceNodes3D(mesh);

  meshGlobalIds(mesh);
  meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, OOGS_AUTO, 0);

  mesh->geometricFactors();

  mesh->o_vgeo.free(); // dfloat version not required
  mesh->o_LMM.free();  // dfloat version not required

  {
    const auto length = mesh->o_ggeo.length();
    auto o_tmp = platform->device.malloc<dfloat>(length);
    mesh->o_ggeo.copyTo(o_tmp);
    mesh->o_ggeo.free();
    mesh->o_ggeo = platform->device.malloc<pfloat>(length);
    platform->copyDfloatToPfloatKernel(length, o_tmp, mesh->o_ggeo);
  }

  {
    const auto length = mesh->o_D.length();
    auto o_tmp = platform->device.malloc<dfloat>(length);
    mesh->o_D.copyTo(o_tmp);
    mesh->o_D.free();
    mesh->o_D = platform->device.malloc<pfloat>(length);
    platform->copyDfloatToPfloatKernel(length, o_tmp, mesh->o_D);
  }

  {
    const auto length = mesh->o_DT.length();
    auto o_tmp = platform->device.malloc<dfloat>(length);
    mesh->o_DT.copyTo(o_tmp);
    mesh->o_DT.free();
    mesh->o_DT = platform->device.malloc<pfloat>(length);
    platform->copyDfloatToPfloatKernel(length, o_tmp, mesh->o_DT);
  }

  return mesh;
}

mesh_t *createMeshV(MPI_Comm comm, int N, int cubN, mesh_t *meshT, occa::properties &kernelInfo)
{
  mesh_t *mesh = new mesh_t();

  if (platform->comm.mpiRank == 0) {
    printf("generating v-mesh ...\n");
  }

  // derive from meshT using shallow copy and adjust only what is needed
  memcpy(mesh, meshT, sizeof(*meshT));
  mesh->cht = 0;

  // find EToV and boundaryInfo
  meshNekReaderHex3D(N, mesh);
  mesh->Nlocal = mesh->Nelements * mesh->Np;

  free(mesh->elementInfo);
  mesh->elementInfo = meshT->elementInfo;
  free(mesh->EX);
  mesh->EX = meshT->EX;
  free(mesh->EY);
  mesh->EY = meshT->EY;
  free(mesh->EZ);
  mesh->EZ = meshT->EZ;

  // find mesh->EToP, mesh->EToE and mesh->EToF, required mesh->EToV
  meshParallelConnect(mesh);

  // find vmapM, based on EToE and EToF
  meshConnectFaceNodes3D(mesh);

  checkEToB(mesh);

  meshVOccaSetup3D(mesh, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nlocal, mesh->globalIds, platform->comm.mpiComm, OOGS_AUTO, 0);

  int err = 0;
  int Nfine;
  platform->options.getArgs("POLYNOMIAL DEGREE", Nfine);
  if (mesh->N == Nfine) {
    dfloat *tmp = (dfloat *)calloc(mesh->Nlocal, sizeof(dfloat));
    mesh->ogs->o_invDegree.copyTo(tmp, mesh->Nlocal);
    auto mult = nek::ptr<double>("vmult");
    dfloat sum1 = 0;
    for (int i = 0; i < mesh->Nlocal; i++) {
      sum1 += std::abs(tmp[i] - static_cast<dfloat>(mult[i]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    if (sum1 > 10 * std::numeric_limits<dfloat>::epsilon()) {
      if (platform->comm.mpiRank == 0) {
        printf("matching invDegree test failed, err=%g!\n", sum1);
      }
      fflush(stdout);
      err++;
    }
    free(tmp);
  }
  nekrsCheck(err, platform->comm.mpiComm, EXIT_FAILURE, "%s\n", "");

  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, OOGS_AUTO);

  mesh->computeInvLMM();

  mesh->volume = platform->linAlg->sum(mesh->Nlocal, mesh->o_LMM, platform->comm.mpiComm);

  {
    std::vector<dfloat> tmp(mesh->Nlocal);
    mesh->ogs->o_invDegree.copyTo(tmp.data());

    double sum = 0;
    for (int i = 0; i < mesh->Nlocal; i++) {
      sum += tmp[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, platform->comm.mpiComm);
    if (platform->comm.mpiRank == 0) {
      printf("unique number of gridpoints: : %lld\n", static_cast<hlong>(sum + 0.1));
    }
  }

  {
    double valMin = (double)mesh->NlocalGatherElements / mesh->Nelements;
    double valMax = (double)mesh->NlocalGatherElements / mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &valMin, 1, MPI_DOUBLE, MPI_MIN, platform->comm.mpiComm);
    MPI_Allreduce(MPI_IN_PLACE, &valMax, 1, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm);

    if (platform->comm.mpiRank == 0 && platform->comm.mpiCommSize > 1) {
      printf("number of interior elements min/max: %2.0f%%  %2.0f%%\n", 100 * valMin, 100 * valMax);
    }
  }

  return mesh;
}

static void meshVOccaSetup3D(mesh_t *mesh, occa::properties &kernelInfo)
{
  mesh->o_EToB = platform->device.malloc<int>(mesh->Nelements * mesh->Nfaces, mesh->EToB);
  mesh->o_vmapM = platform->device.malloc<dlong>(mesh->Nelements * mesh->Nfp * mesh->Nfaces, mesh->vmapM);
  mesh->o_invLMM = platform->device.malloc<dfloat>(mesh->Nelements * mesh->Np);
}
