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
   purpose: read gmsh quadrilateral mesh
 */
mesh_t* meshParallelReaderQuad3D(char* fileName)
{
  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FILE* fp = fopen(fileName, "r");
  int n;

  //  mesh_t *mesh = (mesh_t*) calloc(1, sizeof(mesh_t));
  mesh_t* mesh = new mesh_t[1];

  mesh->rank = rank;
  mesh->size = size;

  MPI_Comm_dup(MPI_COMM_WORLD, &mesh->comm);

  mesh->dim = 3;
  mesh->Nverts = 4; // number of vertices per element
  mesh->Nfaces = 4;
  mesh->NfaceVertices = 2;

  int faceVertices[4][2] = {{0,1},{1,2},{2,3},{3,0}};

  mesh->faceVertices =
    (int*) calloc(mesh->NfaceVertices * mesh->Nfaces, sizeof(int));

  memcpy(mesh->faceVertices, faceVertices[0], mesh->NfaceVertices * mesh->Nfaces * sizeof(int));

  if(fp == NULL) {
    printf("meshReader2D: could not load file %s\n", fileName);
    exit(0);
  }

  char buf[BUFSIZ];
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Nodes"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, hlongFormat, &(mesh->Nnodes));

  /* allocate space for node coordinates */
  dfloat* VX = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat* VY = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));
  dfloat* VZ = (dfloat*) calloc(mesh->Nnodes, sizeof(dfloat));

  /* load nodes */
  for(n = 0; n < mesh->Nnodes; ++n) {
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d" dfloatFormat dfloatFormat dfloatFormat,
           VX + n, VY + n, VZ + n);
  }

  /* look for section with Element node data */
  do{
    fgets(buf, BUFSIZ, fp);
  }while(!strstr(buf, "$Elements"));

  /* read number of nodes in mesh */
  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d", &(mesh->Nelements));

  /* find # of quadrilaterals */
  fpos_t fpos;
  fgetpos(fp, &fpos);
  int Nquadrilaterals = 0;

  int NboundaryFaces = 0;
  for(n = 0; n < mesh->Nelements; ++n) {
    int elementType;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);
    if(elementType == 1) ++NboundaryFaces;
    if(elementType == 3) ++Nquadrilaterals;
  }
  // rewind to start of elements
  fsetpos(fp, &fpos);

  int chunk = Nquadrilaterals / size;
  int remainder = Nquadrilaterals - chunk * size;

  int NquadrilateralsLocal = chunk + (rank < remainder);

  /* where do these elements start ? */
  int start = rank * chunk + mymin(rank, remainder);
  int end = start + NquadrilateralsLocal - 1;

  /* allocate space for Element node index data */

  mesh->EToV
    = (hlong*) calloc(NquadrilateralsLocal * mesh->Nverts,
                      sizeof(hlong));

  mesh->elementInfo
    = (hlong*) calloc(NquadrilateralsLocal,sizeof(hlong));

  /* scan through file looking for quadrilateral elements */
  int cnt = 0, bcnt = 0;
  Nquadrilaterals = 0;

  mesh->boundaryInfo = (hlong*) calloc(NboundaryFaces * 3, sizeof(hlong));
  for(n = 0; n < mesh->Nelements; ++n) {
    int elementType;
    hlong v1, v2, v3, v4;
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d%d", &elementType);

    if(elementType == 1) { // boundary face
      sscanf(buf, "%*d%*d %*d" hlongFormat "%*d " hlongFormat hlongFormat,
             mesh->boundaryInfo + bcnt * 3, &v1, &v2);
      mesh->boundaryInfo[bcnt * 3 + 1] = v1 - 1;
      mesh->boundaryInfo[bcnt * 3 + 2] = v2 - 1;
      ++bcnt;
    }

    if(elementType == 3) { // quadrilateral
      if(start <= Nquadrilaterals && Nquadrilaterals <= end) {
        sscanf(buf, "%*d%*d%*d " hlongFormat " %*d" hlongFormat hlongFormat hlongFormat hlongFormat,
               mesh->elementInfo + cnt, &v1, &v2, &v3, &v4);

#if 0
        // check orientation
        dfloat xe1 = VX[v1 - 1], xe2 = VX[v2 - 1], xe4 = VX[v4 - 1];
        dfloat ye1 = VY[v1 - 1], ye2 = VY[v2 - 1], ye4 = VY[v4 - 1];
        dfloat J = 0.25 * ((xe2 - xe1) * (ye4 - ye1) - (xe4 - xe1) * (ye2 - ye1));
        if(J < 0) {
          int v4tmp = v4;
          v4 = v2;
          v2 = v4tmp;
          printf("unwarping element\n");
        }
#endif

        /* read vertex triplet for trianngle */
        mesh->EToV[cnt * mesh->Nverts + 0] = v1 - 1;
        mesh->EToV[cnt * mesh->Nverts + 1] = v2 - 1;
        mesh->EToV[cnt * mesh->Nverts + 2] = v3 - 1;
        mesh->EToV[cnt * mesh->Nverts + 3] = v4 - 1;
        ++cnt;
      }
      ++Nquadrilaterals;
    }
  }
  fclose(fp);

  /* record number of boundary faces found */
  mesh->NboundaryFaces = bcnt;

  /* record number of found quadrilaterals */
  mesh->Nelements = NquadrilateralsLocal;

  /* collect vertices for each element */
  mesh->EX = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EY = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  mesh->EZ = (dfloat*) calloc(mesh->Nverts * mesh->Nelements, sizeof(dfloat));
  for(int e = 0; e < mesh->Nelements; ++e) {
    for(n = 0; n < mesh->Nverts; ++n) {
      mesh->EX[e * mesh->Nverts + n] = VX[mesh->EToV[e * mesh->Nverts + n]];
      mesh->EY[e * mesh->Nverts + n] = VY[mesh->EToV[e * mesh->Nverts + n]];
      mesh->EZ[e * mesh->Nverts + n] = VZ[mesh->EToV[e * mesh->Nverts + n]];
#if 0
      printf("e %d v %d %g %g %g\n",
             e, n,
             mesh->EX[e * mesh->Nverts + n],
             mesh->EY[e * mesh->Nverts + n],
             mesh->EZ[e * mesh->Nverts + n]);
#endif
    }
  }

  /* release VX and VY (these are too big to keep) */
  free(VX);
  free(VY);
  free(VZ);

  return mesh;
}
