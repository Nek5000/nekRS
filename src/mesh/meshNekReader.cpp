#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "nrs.hpp"
#include "nekInterfaceAdapter.hpp"

void meshNekReaderHex3D(int N, mesh_t *mesh, int isMeshT){

  mesh->dim = nekData.ndim;
  if(mesh->dim != 3) {
    printf("ERROR: No support for %d-D meshes!\n", mesh->dim);
    exit(EXIT_FAILURE);
  }

  if(nekData.nx1 != N+1) {
    printf("ERROR: lx1=%d does not match N=%d\n!\n", nekData.nx1, N);
    exit(EXIT_FAILURE);
  }
 
  mesh->Nverts = 8;
  mesh->Nfaces = 2*mesh->dim;
  mesh->NfaceVertices = 4;
  mesh->Nelements = nekData.nelt;
  if(!isMeshT) mesh->Nelements = nekData.nelv;

  const int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},
                                  {2,3,7,6},{3,0,4,7},{4,5,6,7}};
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(int));

  const int vtxmap[8] = {0, 1, 3, 2, 4, 5, 7, 6};
   
  // build vertex numbering
  mesh->Nnodes = nek_set_glo_num(mesh->Nverts, isMeshT);

  mesh->EToV 
    = (hlong*) calloc(mesh->Nelements * mesh->Nverts, sizeof(hlong));
  for(int e=0; e<mesh->Nelements; ++e) {
    for(int j=0; j<mesh->Nverts; j++) {
      mesh->EToV[e*mesh->Nverts+j] = nekData.glo_num[e*mesh->Nverts+vtxmap[j]];
    }
  }

  // find number of boundary faces
  int nbc = 0;
  int *bid = nekData.boundaryIDt;
  if(!isMeshT) bid = nekData.boundaryID;  
  for(int e = 0; e < mesh->Nelements; e++) {
    for(int iface = 0; iface < mesh->Nfaces; iface++) {
      if(*bid) nbc++;
      bid++;
    }
  }
 
  int *recvCounts = (int *) calloc(mesh->size, sizeof(int));
  MPI_Allgather(&nbc, 1, MPI_INT, recvCounts, 1, MPI_INT, mesh->comm);
  int *displacement = (int *) calloc(mesh->size, sizeof(int));
  displacement[0] = 0;
  for(int i = 1; i < mesh->size; i++) 
    displacement[i] = displacement[i-1] + recvCounts[i-1];

  // build boundary info (for now every rank has all)
  mesh->NboundaryFaces = nbc;
  MPI_Allreduce(MPI_IN_PLACE, &mesh->NboundaryFaces, 1, MPI_HLONG, 
                MPI_SUM, mesh->comm);
  if(mesh->rank == 0) {
    int n = nekData.NboundaryIDt;
    if(!isMeshT) n = nekData.NboundaryID;
    printf("NboundaryIDs: %d\n", n);
    printf("NboundaryFaces: %d\n", mesh->NboundaryFaces);
  }

  int cnt = 0;
  bid = nekData.boundaryIDt;
  if(!isMeshT) bid = nekData.boundaryID;  
  int *eface1 = nekData.eface1;
  int *icface = nekData.icface;
  const double eps = 1e-6;
  const dlong Nfp = nekData.nx1 * nekData.nx1;
  mesh->boundaryInfo = (hlong*) calloc(mesh->NboundaryFaces*(mesh->NfaceVertices+1),sizeof(hlong));
     
  for(int e = 0; e < mesh->Nelements; e++) {
    for(int iface = 0; iface < mesh->Nfaces; iface++) {
      int ibc = *bid;
      if(ibc > 0) {
        hlong offset = displacement[mesh->rank]*(mesh->NfaceVertices+1);
        mesh->boundaryInfo[offset + cnt*(mesh->NfaceVertices+1)] = ibc;
        for(int j = 0; j < mesh->NfaceVertices; j++) {
          const int vertex = icface[j + mesh->NfaceVertices*(eface1[iface] - 1)] - 1;
          mesh->boundaryInfo[offset + cnt*(mesh->NfaceVertices + 1) + j + 1] =
            mesh->EToV[e*mesh->Nverts + vtxmap[vertex]]; 
        }
        cnt++;
      }
      bid++;
    }
  }

  // hack to avoid missing large-count in MPI
  MPI_Datatype bInfoType;
  MPI_Type_contiguous(mesh->NfaceVertices+1, MPI_HLONG, &bInfoType);
  MPI_Type_commit(&bInfoType);

  MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, mesh->boundaryInfo,
      (const int *)recvCounts, (const int *)displacement, bInfoType, mesh->comm);

  free(recvCounts);
  free(displacement);

  // assign vertex coords
  mesh->elementInfo
    = (hlong*) calloc(mesh->Nelements, sizeof(hlong));

  double *VX = nekData.xc;
  double *VY = nekData.yc;
  double *VZ = nekData.zc;
  mesh->EX = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nelements * mesh->Nverts, sizeof(dfloat));
  for(int e=0;e<mesh->Nelements;++e){
    for(int j=0; j<mesh->Nverts; j++) {
      mesh->EX[e*mesh->Nverts+j] = VX[e*mesh->Nverts+j];
      mesh->EY[e*mesh->Nverts+j] = VY[e*mesh->Nverts+j];
      mesh->EZ[e*mesh->Nverts+j] = VZ[e*mesh->Nverts+j];
    }
    mesh->elementInfo[e] = 1; // solid 
    if(e < nekData.nelv ) mesh->elementInfo[e] = 0;
  }
}
