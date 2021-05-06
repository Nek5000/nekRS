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

occa::properties populateMeshProperties(mesh_t* mesh)
{
  occa::properties meshProperties = platform->kernelInfo;

  meshProperties["defines/" "p_dim"] = 3;
  meshProperties["defines/" "p_Nfields"] = mesh->Nfields;
  meshProperties["defines/" "p_N"] = mesh->N;
  meshProperties["defines/" "p_Nq"] = mesh->N + 1;
  meshProperties["defines/" "p_Np"] = mesh->Np;
  meshProperties["defines/" "p_Nfp"] = mesh->Nfp;
  meshProperties["defines/" "p_Nfaces"] = mesh->Nfaces;
  meshProperties["defines/" "p_NfacesNfp"] = mesh->Nfp * mesh->Nfaces;

  meshProperties["defines/" "p_Nvgeo"] = mesh->Nvgeo;
  meshProperties["defines/" "p_Nsgeo"] = mesh->Nsgeo;
  meshProperties["defines/" "p_Nggeo"] = mesh->Nggeo;

  meshProperties["defines/" "p_NXID"] = NXID;
  meshProperties["defines/" "p_NYID"] = NYID;
  meshProperties["defines/" "p_NZID"] = NZID;
  meshProperties["defines/" "p_SJID"] = SJID;
  meshProperties["defines/" "p_IJID"] = IJID;
  meshProperties["defines/" "p_IHID"] = IHID;
  meshProperties["defines/" "p_WSJID"] = WSJID;
  meshProperties["defines/" "p_WIJID"] = WIJID;
  meshProperties["defines/" "p_STXID"] = STXID;
  meshProperties["defines/" "p_STYID"] = STYID;
  meshProperties["defines/" "p_STZID"] = STZID;
  meshProperties["defines/" "p_SBXID"] = SBXID;
  meshProperties["defines/" "p_SBYID"] = SBYID;
  meshProperties["defines/" "p_SBZID"] = SBZID;

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
void loadKernels(mesh_t* mesh, occa::properties kernelInfo);

void meshDummyHex3D(int N, mesh_t* mesh)
{
  mesh->cht = 0;
  mesh->Nfields = 1;
  mesh->dim = 3;
  mesh->Nverts = 8; // number of vertices per element
  mesh->Nfaces = 6;
  mesh->NfaceVertices = 4;

  // vertices on each face
  int faceVertices[6][4] =
  {{0,1,2,3},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{4,5,6,7}};

  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices * mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices * mesh->Nfaces * sizeof(int));

  // build an NX x NY x NZ box grid
  hlong NX = 3, NY = 3, NZ = 3;
  dfloat XMIN = -1, XMAX = +1;
  dfloat YMIN = -1, YMAX = +1;
  dfloat ZMIN = -1, ZMAX = +1;

  hlong allNelements = NX * NY * NZ;

  hlong chunkNelements = allNelements / platform->comm.mpiCommSize;

  hlong start = chunkNelements * platform->comm.mpiRank;
  hlong end   = chunkNelements * (platform->comm.mpiRank + 1);

  if(platform->comm.mpiRank == (platform->comm.mpiCommSize - 1))
    end = allNelements;

  mesh->Nnodes = NX * NY * NZ;
  mesh->Nelements = end - start;
  mesh->NboundaryFaces = 0;

  mesh->EToV = (hlong*) calloc(mesh->Nelements * mesh->Nverts, sizeof(hlong));

  mesh->EX = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));

  mesh->elementInfo = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  // [0,NX]
  dfloat dx = (XMAX - XMIN) / NX; // xmin+0*dx, xmin + NX*(XMAX-XMIN)/NX
  dfloat dy = (YMAX - YMIN) / NY;
  dfloat dz = (ZMAX - ZMIN) / NZ;
  for(hlong n = start; n < end; ++n) {
    int i = n % NX;        // [0, NX)
    int j = (n / NY) % NZ; // [0, NY)
    int k = n / (NX * NY); // [0, NZ)

    hlong e = n - start;

    int ip = (i + 1) % NX;
    int jp = (j + 1) % NY;
    int kp = (k + 1) % NZ;

    mesh->EToV[e * mesh->Nverts + 0] = i  +  j * NX + k * NX * NY;
    mesh->EToV[e * mesh->Nverts + 1] = ip +  j * NX + k * NX * NY;
    mesh->EToV[e * mesh->Nverts + 2] = ip + jp * NX + k * NX * NY;
    mesh->EToV[e * mesh->Nverts + 3] = i  + jp * NX + k * NX * NY;

    mesh->EToV[e * mesh->Nverts + 4] = i  +  j * NX + kp * NX * NY;
    mesh->EToV[e * mesh->Nverts + 5] = ip +  j * NX + kp * NX * NY;
    mesh->EToV[e * mesh->Nverts + 6] = ip + jp * NX + kp * NX * NY;
    mesh->EToV[e * mesh->Nverts + 7] = i  + jp * NX + kp * NX * NY;

    dfloat xo = XMIN + dx * i;
    dfloat yo = YMIN + dy * j;
    dfloat zo = ZMIN + dz * k;

    dfloat* ex = mesh->EX + e * mesh->Nverts;
    dfloat* ey = mesh->EY + e * mesh->Nverts;
    dfloat* ez = mesh->EZ + e * mesh->Nverts;

    ex[0] = xo;
    ey[0] = yo;
    ez[0] = zo;
    ex[1] = xo + dx;
    ey[1] = yo;
    ez[1] = zo;
    ex[2] = xo + dx;
    ey[2] = yo + dy;
    ez[2] = zo;
    ex[3] = xo;
    ey[3] = yo + dy;
    ez[3] = zo;

    ex[4] = xo;
    ey[4] = yo;
    ez[4] = zo + dz;
    ex[5] = xo + dx;
    ey[5] = yo;
    ez[5] = zo + dz;
    ex[6] = xo + dx;
    ey[6] = yo + dy;
    ez[6] = zo + dz;
    ex[7] = xo;
    ey[7] = yo + dy;
    ez[7] = zo + dz;

    mesh->elementInfo[e] = 0;
  }

  mesh->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
  mesh->boundaryInfo = NULL; // no boundaries
}

mesh_t *createMesh(MPI_Comm comm,
                   int N,
                   int cubN,
                   bool cht,
                   occa::properties& kernelInfo)
{
  mesh_t *mesh = new mesh_t();
  int buildOnly = 0;
  if(platform->options.compareArgs("BUILD ONLY", "TRUE")) buildOnly = 1;
  platform->options.getArgs("MESH INTEGRATION ORDER", mesh->nAB);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  mesh->cht  = cht;

  // get mesh from nek
  if(buildOnly)
    meshDummyHex3D(N, mesh);
  else
    meshNekReaderHex3D(N, mesh);

  if (platform->comm.mpiRank == 0)
    printf("generating mesh ... ");
      
  mesh->Nfields = 1; // TW: note this is a temporary patch (halo exchange depends on nfields)

  // connect elements using parallel sort
  meshParallelConnect(mesh);

  // connect elements to boundary faces
  if(!buildOnly) meshConnectBoundary(mesh);

  // load reference (r,s,t) element nodes
  meshLoadReferenceNodesHex3D(mesh, N, cubN);
  if (platform->comm.mpiRank == 0)
    printf("Nq: %d cubNq: %d\n", mesh->Nq, mesh->cubNq);

  mesh->Nlocal = mesh->Nelements * mesh->Np;

  occa::properties meshKernelInfo = populateMeshProperties(mesh);
  loadKernels(mesh, meshKernelInfo);

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
  if(!buildOnly) bcMap::check(mesh);

  meshOccaSetup3D(mesh, platform->options, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, oogsMode);

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

  occa::properties meshKernelInfo = populateMeshProperties(mesh);
  loadKernels(mesh, meshKernelInfo);

  meshHaloSetup(mesh);
  meshPhysicalNodesHex3D(mesh);
  meshHaloPhysicalNodes(mesh);
  meshGeometricFactorsHex3D(mesh);
  meshConnectFaceNodes3D(mesh);
  meshSurfaceGeometricFactorsHex3D(mesh);
  meshGlobalIds(mesh);

  bcMap::check(mesh);
  meshOccaSetup3D(mesh, platform->options, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, oogsMode);

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

  return mesh;
}
*/

mesh_t *createMeshV(
                    MPI_Comm comm,
                    int N,
                    int cubN,
                    mesh_t* meshT,
                    occa::properties& kernelInfo)
{
  mesh_t *mesh = new mesh_t();

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

  // find mesh->EToB, required mesh->EToV and mesh->boundaryInfo
  meshConnectBoundary(mesh);

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

  meshVOccaSetup3D(mesh, kernelInfo);

  meshParallelGatherScatterSetup(mesh, mesh->Nelements * mesh->Np, mesh->globalIds, platform->comm.mpiComm, 0);
  oogs_mode oogsMode = OOGS_AUTO; 
  //if(platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP") oogsMode = OOGS_DEFAULT;
  mesh->oogs = oogs::setup(mesh->ogs, 1, mesh->Nelements * mesh->Np, ogsDfloat, NULL, oogsMode);

  mesh->computeInvLMM();

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

void loadKernels(mesh_t* mesh, occa::properties kernelInfo)
{
  if(platform->options.compareArgs("MOVING MESH", "TRUE")){
    std::string install_dir;
    install_dir.assign(getenv("NEKRS_INSTALL_DIR"));
    std::string oklpath = install_dir + "/okl/";
    {
        std::string filename = oklpath + "mesh/velocityBCHex3D.okl";
        mesh->velocityDirichletKernel =
          platform->device.buildKernel(filename,
                                   "velocityDirichletBCHex3D",
                                   kernelInfo);
        occa::properties meshKernelInfo = kernelInfo;
        meshKernelInfo["defines/" "p_cubNq"] = mesh->cubNq;
        meshKernelInfo["defines/" "p_cubNp"] = mesh->cubNp;

        filename = oklpath + "mesh/geometricFactorsHex3D.okl";
        mesh->geometricFactorsKernel =
          platform->device.buildKernel(filename,
                                   "geometricFactorsHex3D",
                                   meshKernelInfo);
        filename = oklpath + "mesh/surfaceGeometricFactorsHex3D.okl";
        mesh->surfaceGeometricFactorsKernel =
          platform->device.buildKernel(filename,
                                   "surfaceGeometricFactorsHex3D",
                                   meshKernelInfo);

        meshKernelInfo = kernelInfo;
        meshKernelInfo["defines/" "p_nAB"] = mesh->nAB;
        filename = oklpath + "core/nStagesSum.okl";
        mesh->nStagesSumVectorKernel =
          platform->device.buildKernel(filename,
                                   "nStagesSumVector",
                                   meshKernelInfo);
    }
  }
}
