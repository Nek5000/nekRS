#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void meshNekReaderHex3D(int N, mesh_t* mesh)
{
  
  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0) printf("loading mesh from nek ... "); fflush(stdout); 

  mesh->dim = 3; 
  mesh->Nverts = 8;
  mesh->Nfaces = 2 * mesh->dim;
  mesh->NfaceVertices = 4;
  mesh->Nelements = nekData.nelt;
  if(!mesh->cht) mesh->Nelements = nekData.nelv;

  const int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},
    {2,3,7,6},{3,0,4,7},{4,5,6,7}};
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices * mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices * mesh->Nfaces * sizeof(int));

  // pre-processor maps
  const int vtxMap[] = {0,1,3,2,4,5,7,6};
  const int faceMap[] = {1,2,3,4,0,5};

  // generate element vertex numbering
  mesh->Nnodes = nek::set_glo_num(2, mesh->cht);

  mesh->EToV
    = (hlong*) calloc(mesh->Nelements * mesh->Nverts, sizeof(hlong));
  for(int e = 0; e < mesh->Nelements; ++e)
    for(int j = 0; j < mesh->Nverts; j++)
      mesh->EToV[e * mesh->Nverts + j] = nekData.glo_num[e * mesh->Nverts + vtxMap[j]];

  // find number of boundary faces
  hlong NboundaryFaces = 0;
  int* bid = nekData.boundaryIDt;
  if(!mesh->cht) bid = nekData.boundaryID;
  for(int e = 0; e < mesh->Nelements; e++)
    for(int iface = 0; iface < mesh->Nfaces; iface++) {
      if(*bid) NboundaryFaces++;
      bid++;
    }
 
  MPI_Allreduce(MPI_IN_PLACE, &NboundaryFaces, 1, MPI_HLONG,
                MPI_SUM, platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) {
    int n = nekData.NboundaryIDt;
    if(!mesh->cht) n = nekData.NboundaryID;
    printf("NboundaryIDs: %d, NboundaryFaces: %lld ", n, NboundaryFaces);
  }

  // boundary face tags (face numbering is in pre-processor notation)
  mesh->EToB = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
  for(int i = 0; i < mesh->Nelements * mesh->Nfaces; i++) mesh->EToB[i] = -1;

  bid = nekData.boundaryIDt;
  if(!mesh->cht) bid = nekData.boundaryID;

  for(int e = 0; e < mesh->Nelements; e++) {
    for(int i = 0; i < mesh->Nfaces; i++) {
      const int ibc = bid[e * mesh->Nfaces + i];
      if(ibc > 0) mesh->EToB[e * mesh->Nfaces + faceMap[i]] = ibc; 
    }
  }

  // assign vertex coords
  mesh->elementInfo
    = (dlong*) calloc(mesh->Nelements, sizeof(dlong));

  double* VX = nekData.xc;
  double* VY = nekData.yc;
  double* VZ = nekData.zc;
  mesh->EX = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  for(int e = 0; e < mesh->Nelements; ++e) {
    for(int j = 0; j < mesh->Nverts; j++) {
      mesh->EX[e * mesh->Nverts + j] = VX[e * mesh->Nverts + j];
      mesh->EY[e * mesh->Nverts + j] = VY[e * mesh->Nverts + j];
      mesh->EZ[e * mesh->Nverts + j] = VZ[e * mesh->Nverts + j];
    }
    mesh->elementInfo[e] = 1; // solid
    if(e < nekData.nelv ) mesh->elementInfo[e] = 0;
  }

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}
