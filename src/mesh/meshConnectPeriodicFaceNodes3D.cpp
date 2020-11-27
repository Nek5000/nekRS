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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh3D.h"

int findBestPeriodicMatch(dfloat xper, dfloat yper, dfloat zper,
                          dfloat x1, dfloat y1, dfloat z1,
                          int Np2, int* nodeList, dfloat* x2, dfloat* y2, dfloat* z2, int* nP)
{
  int matchIndex;
  dfloat mindist2 = 1e9;
  int isFirst = 1;

  for(int n = 0; n < Np2; ++n) {
    /* next node */
    const int i2 = nodeList[n];
    for(int zp = 0; zp < 2; ++zp)
      for(int yp = 0; yp < 2; ++yp)
        for(int xp = 0; xp < 2; ++xp) {
          /* distance between target and next node */
          const dfloat dist2 =
            pow(fabs(x1 - x2[i2]) - xp * xper,2) +
            pow(fabs(y1 - y2[i2]) - yp * yper,2) +
            pow(fabs(z1 - z2[i2]) - zp * zper,2);

          /* if next node is closer to target update match */
          if(isFirst == 1 || dist2 < mindist2) {
            mindist2 = dist2;
            matchIndex = i2;
            *nP = n;
            isFirst = 0;
          }
        }
  }
  if(mindist2 > 1e-3) printf("arggh - bad match: x,y,z= %g,%g,%g => %g,%g,%g with mindist=%lg\n",
                             x1,y1,z1,  x2[matchIndex], y2[matchIndex],  z2[matchIndex], mindist2);

  return matchIndex;
}

// serial face-node to face-node connection
void meshConnectPeriodicFaceNodes3D(mesh3D* mesh, dfloat xper, dfloat yper, dfloat zper)
{
  /* volume indices of the interior and exterior face nodes for each element */
  mesh->vmapM = (dlong*) calloc(mesh->Nfp * mesh->Nfaces * mesh->Nelements, sizeof(dlong));
  mesh->vmapP = (dlong*) calloc(mesh->Nfp * mesh->Nfaces * mesh->Nelements, sizeof(dlong));
  mesh->mapP  = (dlong*) calloc(mesh->Nfp * mesh->Nfaces * mesh->Nelements, sizeof(dlong));

  /* assume elements already connected */
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int f = 0; f < mesh->Nfaces; ++f) {
      dlong eP = mesh->EToE[e * mesh->Nfaces + f];
      int fP = mesh->EToF[e * mesh->Nfaces + f];
      if(eP < 0 || fP < 0) { // fake connections for unconnected faces
        eP = e;
        fP = f;
      }
      /* for each node on this face find the neighbor node */
      for(int n = 0; n < mesh->Nfp; ++n) {
        dlong idM = mesh->faceNodes[f * mesh->Nfp + n] + e * mesh->Np;
        dfloat xM = mesh->x[idM];
        dfloat yM = mesh->y[idM];
        dfloat zM = mesh->z[idM];
        int nP;

        int idP = findBestPeriodicMatch(xper, yper, zper,
                                        xM, yM, zM,
                                        mesh->Nfp,
                                        mesh->faceNodes + fP * mesh->Nfp,
                                        mesh->x + eP * mesh->Np,
                                        mesh->y + eP * mesh->Np,
                                        mesh->z + eP * mesh->Np, &nP);

        dlong id = mesh->Nfaces * mesh->Nfp * e + f * mesh->Nfp + n;
        mesh->vmapM[id] = idM;
        mesh->vmapP[id] = idP + eP * mesh->Np;
        mesh->mapP[id] = eP * mesh->Nfaces * mesh->Nfp + fP * mesh->Nfp + nP;
      }
    }
}
