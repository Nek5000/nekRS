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
#include  "mpi.h"

#include "mesh3D.h"

/*
   purpose: read gmsh tetrahedra mesh
 */
mesh3D* meshParallelReaderTet3D(char* fileName)
{
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE* fp = fopen(fileName, "r");

  char* status;

  mesh_t* mesh = new mesh_t();

  mesh->rank = rank;
  mesh->size = size;

  MPI_Comm_dup(MPI_COMM_WORLD, &mesh->comm);

  mesh->dim = 3;
  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;

  // vertices on each face
  int faceVertices[4][3] = {{0,1,2},{0,1,3},{1,2,3},{2,0,3}};
  mesh->NfaceVertices = 3;
  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices * mesh->Nfaces, sizeof(int));
  memcpy(mesh->faceVertices, faceVertices[0], 12 * sizeof(int));

  if(fp == NULL) {
    printf("meshReaderTet3D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat* VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat* VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat* VZ = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));

  /* load nodes */
  for(hlong n = 0; n < mesh->Nnodes; ++n) {
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           VX + n, VY + n, VZ + n);
  }

  /* look for section with Element node data */
  do{
    status = fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  hlong Nelements;
  status = fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &Nelements);

  /* find # of tets */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  hlong Ntets = 0, NboundaryFaces = 0;
  for(hlong n = 0; n < Nelements; ++n) {
    int elementType;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType == 4) ++Ntets; // tet code is 4
    if(elementType == 2) ++NboundaryFaces;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  hlong chunk = (hlong) Ntets / size;
  int remainder = (int) (Ntets - chunk * size);

  hlong NtetsLocal = chunk + (rank < remainder);

  /* where do these elements start ? */
  hlong start = rank * chunk + mymin(rank, remainder);
  hlong end = start + NtetsLocal - 1;

  /* allocate space for Element node index data */

  mesh->EToV
    = (hlong*) calloc(NtetsLocal * mesh->Nverts, sizeof(hlong));
  mesh->elementInfo
    = (hlong*) calloc(NtetsLocal,sizeof(hlong));

  /* scan through file looking for tetrahedra elements */
  hlong cnt = 0, bcnt = 0;
  Ntets = 0;

  mesh->boundaryInfo = (hlong*) calloc(NboundaryFaces * 4, sizeof(hlong));
  for(hlong n = 0; n < Nelements; ++n) {
    int elementType;
    hlong v1, v2, v3, v4;
    status = fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType == 2) { // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d" hlongFormat hlongFormat hlongFormat,
             mesh->boundaryInfo + bcnt * 4, &v1, &v2, &v3);
      mesh->boundaryInfo[bcnt * 4 + 1] = v1 - 1;
      mesh->boundaryInfo[bcnt * 4 + 2] = v2 - 1;
      mesh->boundaryInfo[bcnt * 4 + 3] = v3 - 1;
      ++bcnt;
    }

    if(elementType == 4) { // tet code is 4
      if(start <= Ntets && Ntets <= end) {
        sscanf(buf,
               "%*d%*d%*d " hlongFormat " %*d"
               hlongFormat hlongFormat hlongFormat hlongFormat,
               mesh->elementInfo + cnt,&v1, &v2, &v3, &v4);
        /* read vertex triplet for trianngle */
        mesh->EToV[cnt * mesh->Nverts + 0] = v1 - 1;
        mesh->EToV[cnt * mesh->Nverts + 1] = v2 - 1;
        mesh->EToV[cnt * mesh->Nverts + 2] = v3 - 1;
        mesh->EToV[cnt * mesh->Nverts + 3] = v4 - 1;
        ++cnt;
      }
      ++Ntets;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;

  /* record number of found tets */
  mesh->Nelements = (dlong) NtetsLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int n = 0; n < mesh->Nverts; ++n) {
      hlong vid = mesh->EToV[e * mesh->Nverts + n];
      mesh->EX[e * mesh->Nverts + n] = VX[vid];
      mesh->EY[e * mesh->Nverts + n] = VY[vid];
      mesh->EZ[e * mesh->Nverts + n] = VZ[vid];
    }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);

  return mesh;
}
