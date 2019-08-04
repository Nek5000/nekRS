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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"
#include "libParanumal.hpp"
#include "nekInterfaceAdapter.hpp"

void meshNekReaderHex3D(int N, mesh_t *mesh){

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
  mesh->Nnodes = nekData.ngv;

  const int faceVertices[6][4] = {{0,1,2,3},{0,1,5,4},{1,2,6,5},
                                  {2,3,7,6},{3,0,4,7},{4,5,6,7}};
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices*mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices*mesh->Nfaces*sizeof(int));

  const int vtxmap[8] = {0, 1, 3, 2, 4, 5, 7, 6};
   
  // vertex numberbing
  mesh->EToV 
    = (hlong*) calloc(mesh->Nelements * mesh->Nverts, sizeof(hlong));
  hlong *glo_num = nekData.glo_num;
  for(int e=0; e<mesh->Nelements; ++e) {
    for(int j=0; j<mesh->Nverts; j++) {
      mesh->EToV[e*mesh->Nverts+j] = glo_num[e*mesh->Nverts+vtxmap[j]];
    }
  }

  // find number of boundary faces and check if supported
  int nbc = 0;
  char *cbc = nekData.cbc;
  for(int e = 0; e < mesh->Nelements; e++) {
    for(int iface = 0; iface < mesh->Nfaces; iface++) {
      if(!strncmp(cbc, "W  ", 3)) {
        nbc++;
      } else if(!strncmp(cbc, "v  ", 3) || !strncmp(cbc, "SYM", 3)) {
        nbc++;
      } else if(!strncmp(cbc, "o  ", 3) || !strncmp(cbc, "O  ", 3)) {
        nbc++;
      } else if(!strncmp(cbc, "E  ", 3) || !strncmp(cbc, "P  ",3 )) {
        ({;});
      } else {
        printf("ERROR: Unsupported nek bcType %s found!\n", cbc);
        exit(EXIT_FAILURE);
      }
      cbc += 3;
    }
  }
 
  int *recvCounts = (int *) calloc(mesh->size, sizeof(int));
  MPI_Allgather(&nbc, 1, MPI_INT, recvCounts, 1, MPI_INT, mesh->comm);
  int *displacement = (int *) calloc(mesh->size, sizeof(int));
  displacement[0] = 0;
  for(int i = 1; i < mesh->size; i++) {
    displacement[i] = displacement[i-1] + recvCounts[i-1];
  }

  // build boundary info (for now every rank has all)
  mesh->NboundaryFaces = nbc;
  MPI_Allreduce(MPI_IN_PLACE, &mesh->NboundaryFaces, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  if(mesh->rank == 0) printf("NboundaryFaces: %d\n", mesh->NboundaryFaces);

  int cnt = 0;
  cbc = nekData.cbc;
  int *eface1 = nekData.eface1;
  int *icface = nekData.icface;
  const double eps = 1e-6;
  const dlong Nfp = nekData.nx1 * nekData.nx1;
  mesh->boundaryInfo = (hlong*) calloc(mesh->NboundaryFaces*(mesh->NfaceVertices+1),sizeof(hlong));
     
  for(int e = 0; e < mesh->Nelements; e++) {
    for(int iface = 0; iface < mesh->Nfaces; iface++) {
      int ibc = -1;
      if(!strncmp(cbc, "W  ", 3)) {
        ibc = 1;
      } else if(!strncmp(cbc, "v  ", 3)) {
        ibc = 2;
      } else if(!strncmp(cbc, "o  ", 3) || !strncmp(cbc, "O  ", 3)) {
        ibc = 3;
      } else if(!strncmp(cbc, "SYM", 3)) {
        const int id = iface*Nfp + e*6*Nfp;
        const double nx = nekData.unx[id];
        const double ny = nekData.uny[id];
        const double nz = nekData.unz[id];
        if(abs(nx*nx-1) < eps && abs(ny*ny-0) < eps && abs(nz*nz-0) < eps) ibc=4; 
        if(abs(nx*nx-0) < eps && abs(ny*ny-1) < eps && abs(nz*nz-0) < eps) ibc=5; 
        if(abs(nx*nx-0) < eps && abs(ny*ny-0) < eps && abs(nz*nz-1) < eps) ibc=6; 
        if(ibc == -1) { 
          printf("ERROR: SYM boundary (n=%g,%g,%g) needs to be aligned!\n", nx, ny, nz);
          exit(EXIT_FAILURE);
        }
      }

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

      cbc+=3;
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
    = (hlong*) calloc(mesh->Nelements * mesh->Nverts, sizeof(hlong));

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
    mesh->elementInfo[e] = 1; // dummy value
  }
}
