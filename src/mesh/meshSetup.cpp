#include "nrs.hpp"
#include "bcMap.hpp"
#include "meshNekReader.hpp"
#include <string>

void meshVOccaSetup3D(mesh_t* mesh, occa::properties &kernelInfo);
mesh_t *createMeshV(MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    occa::properties& kernelInfo);

occa::properties meshKernelProperties(int N)
{
  occa::properties meshProperties;
  const int Nq = N+1;
  const int Np = Nq * Nq * Nq;
  const int Nfp = Nq * Nq;
  constexpr int Nfaces {6};

  constexpr int Nvgeo{12};
  constexpr int Nggeo{7};
  constexpr int Nsgeo{13};

  meshProperties["defines/" "p_dim"] = 3;
  meshProperties["defines/" "p_Nverts"] = 8;
  meshProperties["defines/" "p_Nfields"] = 1;
  meshProperties["defines/" "p_N"] = N;
  meshProperties["defines/" "p_Nq"] = Nq;
  meshProperties["defines/" "p_Nq_g"] = Nq;
  meshProperties["defines/" "p_Np"] = Np;
  meshProperties["defines/" "p_Np_g"] = Np;
  meshProperties["defines/" "p_Nfp"] = Nfp;
  meshProperties["defines/" "p_Nfaces"] = Nfaces;
  meshProperties["defines/" "p_NfacesNfp"] = Nfp * Nfaces;

  meshProperties["defines/" "p_Nvgeo"] = Nvgeo;
  meshProperties["defines/" "p_Nsgeo"] = Nsgeo;
  meshProperties["defines/"
                 "p_Nggeo"] = Nggeo;

  meshProperties["defines/" "p_NXID"] = NXID;
  meshProperties["defines/" "p_NYID"] = NYID;
  meshProperties["defines/" "p_NZID"] = NZID;
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

  meshProperties["defines/" "p_G00ID"] = G00ID;
  meshProperties["defines/" "p_G01ID"] = G01ID;
  meshProperties["defines/" "p_G02ID"] = G02ID;
  meshProperties["defines/" "p_G11ID"] = G11ID;
  meshProperties["defines/" "p_G12ID"] = G12ID;
  meshProperties["defines/" "p_G22ID"] = G22ID;
  meshProperties["defines/" "p_GWJID"] = GWJID;

  meshProperties["defines/" "p_RXID"] = RXID;
  meshProperties["defines/" "p_SXID"] = SXID;
  meshProperties["defines/" "p_TXID"] = TXID;

  meshProperties["defines/" "p_RYID"] = RYID;
  meshProperties["defines/" "p_SYID"] = SYID;
  meshProperties["defines/" "p_TYID"] = TYID;

  meshProperties["defines/" "p_RZID"] = RZID;
  meshProperties["defines/" "p_SZID"] = SZID;
  meshProperties["defines/" "p_TZID"] = TZID;

  meshProperties["defines/" "p_JID"] = JID;
  meshProperties["defines/" "p_JWID"] = JWID;
  meshProperties["defines/" "p_IJWID"] = IJWID;
  return meshProperties;
}
void loadKernels(mesh_t* mesh);

mesh_t *createMesh(MPI_Comm comm,
                   int N,
                   int cubN,
                   bool cht,
                   occa::properties& kernelInfo)
{
  mesh_t *mesh = new mesh_t();
  platform->options.getArgs("MESH INTEGRATION ORDER", mesh->nAB);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh->cht  = cht;

  if (platform->comm.mpiRank == 0)
    printf("generating t-mesh ...\n");

  // get mesh from nek
  meshNekReaderHex3D(N, mesh);

  if ((hlong) mesh->Nelements * (mesh->Nvgeo * cubN) > std::numeric_limits<int>::max()) {
    if (platform->comm.mpiRank == 0) 
      printf("ERROR: mesh->Nelements * mesh->Nvgeo * cubN exceeds <int> limit!");
    ABORT(EXIT_FAILURE);
  }

  mesh->Nfields = 1; // TW: note this is a temporary patch (halo exchange depends on nfields)

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N, cubN);
  if (platform->comm.mpiRank == 0) {
    printf("N: %d, Nq: %d", mesh->N, mesh->Nq);
    if (cubN) printf(", cubNq: %d", mesh->cubNq);
    printf("\n");
  }

  loadKernels(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  // compute physical (x,y) locations of the element nodes
  meshPhysicalNodesHex3D(mesh);

  meshHaloPhysicalNodes(mesh);

  // compute geometric factors
  meshGeometricFactorsHex3D(mesh);

  // connect face nodes (find trace indices)
  meshConnectFaceNodes3D(mesh);

  // compute surface geofacs (including halo)
  meshSurfaceGeometricFactorsHex3D(mesh);

  // global nodes
  meshGlobalIds(mesh);
  bcMap::check(mesh);
  bcMap::checkBoundaryAlignment(mesh);
  bcMap::remapUnalignedBoundaries(mesh);

  meshOccaSetup3D(mesh, platform->options, kernelInfo);

  meshParallelGatherScatterSetup(mesh,
                                 mesh->Nelements * mesh->Np,
                                 mesh->globalIds,
                                 platform->comm.mpiComm,
                                 OOGS_AUTO,
                                 0);

  int err = 0;
  int Nfine;
  platform->options.getArgs("POLYNOMIAL DEGREE", Nfine);
  if(mesh->N == Nfine) {
    dfloat* tmp = (dfloat*) calloc(mesh->Nlocal, sizeof(dfloat));
    mesh->ogs->o_invDegree.copyTo(tmp, mesh->Nlocal * sizeof(dfloat));
    double* mult = (cht) ? (double*) nek::ptr("tmult") : (double*) nek::ptr("vmult");
    dfloat sum1 = 0;
    for(int i = 0; i < mesh->Nlocal; i++) sum1 += std::abs(tmp[i] - mult[i]);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    if(sum1 > 1e-14) {
      if(platform->comm.mpiRank == 0) printf("multiplicity test err=%g!\n", sum1);
      fflush(stdout);
      err++;
    }
    free(tmp);
  }
  if(err) ABORT(1);

  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, OOGS_AUTO);

  // build mass + inverse mass matrix
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int n = 0; n < mesh->Np; ++n)
      mesh->LMM[e * mesh->Np + n] = mesh->vgeo[e * mesh->Np * mesh->Nvgeo + JWID * mesh->Np + n];
  mesh->o_LMM.copyFrom(mesh->LMM, mesh->Nelements * mesh->Np * sizeof(dfloat));
  mesh->computeInvLMM();

  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    const int maxTemporalOrder = 3;
    mesh->coeffAB = (dfloat*) calloc(maxTemporalOrder, sizeof(dfloat));
    mesh->o_coeffAB = platform->device.malloc(maxTemporalOrder * sizeof(dfloat), mesh->coeffAB);
  }

  mesh->fluid = mesh;
  if(mesh->cht) mesh->fluid = createMeshV(comm, N, cubN, mesh, kernelInfo); 

  return mesh;
}

/*
mesh_t* duplicateMesh(MPI_Comm comm,
                      int N,
                      int cubN,
                      mesh_t* meshT,
                      occa::device device,
                      occa::properties& kernelInfo)
{
  mesh_t* mesh = new mesh_t[1];

  // shallow copy
  memcpy(mesh, meshT, sizeof(*meshT));

  mesh->Nfields = 1; // TW: note this is a temporary patch (halo exchange depends on nfields)

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N, cubN);
  if (platform->comm.mpiRank == 0)
    printf("Nq: %d cubNq: %d \n", mesh->Nq, mesh->cubNq);

  loadKernels(mesh);

  meshHaloSetup(mesh);
  meshPhysicalNodesHex3D(mesh);
  meshHaloPhysicalNodes(mesh);
  meshGeometricFactorsHex3D(mesh);
  meshConnectFaceNodes3D(mesh);
  meshSurfaceGeometricFactorsHex3D(mesh);
  meshGlobalIds(mesh);

  bcMap::check(mesh);
  bcMap::checkBoundaryAlignment(mesh);
  bcMap::remapUnalignedBoundaries(mesh);

  meshOccaSetup3D(mesh, platform->options, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, OOGS_AUTO, 0);
  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, OOGS_AUTO);

  // build mass + inverse mass matrix
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int n = 0; n < mesh->Np; ++n)
      mesh->LMM[e * mesh->Np + n] = mesh->vgeo[e * mesh->Np * mesh->Nvgeo + JWID * mesh->Np + n];
  mesh->o_LMM.copyFrom(mesh->LMM, mesh->Nelements * mesh->Np * sizeof(dfloat));
  mesh->computeInvLMM();

  return mesh;
}
*/

mesh_t *createMeshMG(mesh_t* _mesh,
                     int Nc)
{
  mesh_t* mesh = new mesh_t();
  memcpy(mesh, _mesh, sizeof(mesh_t));

  meshLoadReferenceNodesHex3D(mesh, Nc, 1);
  meshHaloSetup(mesh);
  meshPhysicalNodesHex3D(mesh);
  meshHaloPhysicalNodes(mesh);
  meshGeometricFactorsHex3D(mesh);

  meshConnectFaceNodes3D(mesh);
  meshSurfaceGeometricFactorsHex3D(mesh);

  meshGlobalIds(mesh);
  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, OOGS_AUTO, 0);

  mesh->o_x = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(dfloat), mesh->x);
  mesh->o_y = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(dfloat), mesh->y);
  mesh->o_z = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(dfloat), mesh->z);

  free(mesh->x);
  free(mesh->y);
  free(mesh->z);

  mesh->o_D = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), mesh->D);

  dfloat* DT = (dfloat*) calloc(mesh->Nq * mesh->Nq, sizeof(dfloat));
  for (int j = 0; j < mesh->Nq; j++)
    for (int i = 0; i < mesh->Nq; i++)
      DT[j * mesh->Nq + i] = mesh->D[i * mesh->Nq + j];
  mesh->o_DT = platform->device.malloc(mesh->Nq * mesh->Nq * sizeof(dfloat), DT);
  free(DT);

  mesh->o_ggeo = platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo * sizeof(dfloat),
                                         mesh->ggeo);

  if(!strstr(pfloatString,dfloatString)) {
    mesh->o_ggeoPfloat = platform->device.malloc(mesh->Nelements * mesh->Np * mesh->Nggeo, sizeof(pfloat));
    mesh->o_DPfloat = platform->device.malloc(mesh->Nq * mesh->Nq, sizeof(pfloat));
    mesh->o_DTPfloat = platform->device.malloc(mesh->Nq * mesh->Nq, sizeof(pfloat));
    platform->copyDfloatToPfloatKernel(mesh->Nelements * mesh->Np * mesh->Nggeo,
                                       mesh->o_ggeo,
                                       mesh->o_ggeoPfloat);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       mesh->o_D,
                                       mesh->o_DPfloat);
    platform->copyDfloatToPfloatKernel(mesh->Nq * mesh->Nq,
                                       mesh->o_DT,
                                       mesh->o_DTPfloat);

    // TODO: once full preconditioner is in FP32, uncomment below
    //mesh->o_D.free();
    //mesh->o_DT.free();
    //mesh->o_ggeo.free();
  }

  return mesh;
}

mesh_t *createMeshV(
                    MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    occa::properties& kernelInfo)
{
  mesh_t *mesh = new mesh_t();

  if (platform->comm.mpiRank == 0)
    printf("generating v-mesh ...\n");

  // shallow copy
  memcpy(mesh, meshT, sizeof(*meshT));
  mesh->cht = 0;

  // find EToV and boundaryInfo
  meshNekReaderHex3D(N, mesh);
  free(mesh->elementInfo);
  mesh->elementInfo = meshT->elementInfo;

  mesh->Nlocal = mesh->Nelements * mesh->Np;

  mesh->Nfields = 1; // temporary patch (halo exchange depends on nfields)

  // find mesh->EToP, mesh->EToE and mesh->EToF, required mesh->EToV
  meshParallelConnect(mesh);

  // set up halo exchange info for MPI (do before connect face nodes)
  meshHaloSetup(mesh);

  //meshPhysicalNodesHex3D(mesh);
  mesh->x = meshT->x;
  mesh->y = meshT->y;
  mesh->z = meshT->z;

  meshHaloPhysicalNodes(mesh);

  // meshGeometricFactorsHex3D(mesh);
  mesh->vgeo = meshT->vgeo;
  mesh->cubvgeo = meshT->cubvgeo;
  mesh->ggeo = meshT->ggeo;

  // connect face nodes (find trace indices)
  // find vmapM, vmapP, mapP based on EToE and EToF
  meshConnectFaceNodes3D(mesh);

  // meshGlobalIds(mesh);
  mesh->globalIds = meshT->globalIds;

  bcMap::check(mesh);
  bcMap::checkBoundaryAlignment(mesh);
  bcMap::remapUnalignedBoundaries(mesh);

  meshVOccaSetup3D(mesh, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, OOGS_AUTO, 0);

  int err = 0;
  int Nfine;
  platform->options.getArgs("POLYNOMIAL DEGREE", Nfine);
  if(mesh->N == Nfine) {
    dfloat* tmp = (dfloat*) calloc(mesh->Nlocal, sizeof(dfloat));
    mesh->ogs->o_invDegree.copyTo(tmp, mesh->Nlocal * sizeof(dfloat));
    double* mult = (double*) nek::ptr("vmult");
    dfloat sum1 = 0;
    for(int i = 0; i < mesh->Nlocal; i++) sum1 += std::abs(tmp[i] - mult[i]);
    MPI_Allreduce(MPI_IN_PLACE, &sum1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    if(sum1 > 1e-14) {
      if(platform->comm.mpiRank == 0) printf("multiplicity test err=%g!\n", sum1);
      fflush(stdout);
      err++;
    }
    free(tmp);
  }
  if(err) ABORT(1);

  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, OOGS_AUTO);

  mesh->computeInvLMM();

  // compute V mesh volume
  dfloat volume = 0.0;
  const auto Np = mesh->Np;
  const auto Nggeo = mesh->Nggeo;
  for(dlong e = 0; e < mesh->Nelements; ++e) {
    for(dlong n = 0; n < Np; ++n){
        volume += mesh->ggeo[Nggeo * Np * e + n + Np * GWJID];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
  mesh->volume = volume;


  return mesh;
}

void meshVOccaSetup3D(mesh_t* mesh, occa::properties &kernelInfo)
{
  if(mesh->totalHaloPairs > 0) {
    // copy halo element list to DEVICE
    mesh->o_haloElementList =
      platform->device.malloc(mesh->totalHaloPairs * sizeof(dlong), mesh->haloElementList);

    // temporary DEVICE buffer for halo (maximum size Nfields*Np for dfloat)
    mesh->o_haloBuffer =
      platform->device.malloc(mesh->totalHaloPairs * mesh->Np * mesh->Nfields ,  sizeof(dfloat));

    // node ids
    mesh->o_haloGetNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloGetNodeIds);

    mesh->o_haloPutNodeIds =
      platform->device.malloc(mesh->Nfp * mesh->totalHaloPairs * sizeof(dlong), mesh->haloPutNodeIds);
  }

  mesh->o_EToB =
    platform->device.malloc(mesh->Nelements * mesh->Nfaces * sizeof(int),
                        mesh->EToB);
  mesh->o_vmapM =
    platform->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapM);
  mesh->o_vmapP =
    platform->device.malloc(mesh->Nelements * mesh->Nfp * mesh->Nfaces * sizeof(dlong),
                        mesh->vmapP);
  mesh->o_invLMM =
    platform->device.malloc(mesh->Nelements * mesh->Np ,  sizeof(dfloat));
}

void loadKernels(mesh_t* mesh)
{
  const std::string meshPrefix = "mesh-";
  mesh->avgBIDValueKernel = platform->kernels.get(meshPrefix + "avgBIDValue");
  mesh->velocityDirichletKernel = platform->kernels.get(meshPrefix + "velocityDirichletBCHex3D");
  mesh->geometricFactorsKernel = platform->kernels.get(meshPrefix + "geometricFactorsHex3D");
  mesh->surfaceGeometricFactorsKernel = platform->kernels.get(meshPrefix + "surfaceGeometricFactorsHex3D");
  mesh->cubatureGeometricFactorsKernel = platform->kernels.get(meshPrefix + "cubatureGeometricFactorsHex3D");
  mesh->nStagesSumVectorKernel = platform->kernels.get(meshPrefix + "nStagesSumVector");
}
