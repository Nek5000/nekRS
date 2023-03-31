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
#include <stdio.h>
#include <stdlib.h>
#include "mesh3D.h"

void interpolateFaceHex3D(int* faceNodes, dfloat* I, dfloat* x, int N, dfloat* Ix, int M)
{
  dfloat* Ix0 = (dfloat*) calloc(N * N, sizeof(dfloat));
  dfloat* Ix1 = (dfloat*) calloc(N * M, sizeof(dfloat));

  for(int j = 0; j < N; ++j)
    for(int i = 0; i < N; ++i)
      Ix0[j * N + i] = x[faceNodes[j * N + i]];

  for(int j = 0; j < N; ++j)
    for(int i = 0; i < M; ++i) {
      dfloat tmp = 0;
      for(int n = 0; n < N; ++n)
        tmp += I[i * N + n] * Ix0[j * N + n];
      Ix1[j * M + i] = tmp;
    }

  for(int j = 0; j < M; ++j)
    for(int i = 0; i < M; ++i) {
      dfloat tmp = 0;
      for(int n = 0; n < N; ++n)
        tmp += I[j * N + n] * Ix1[n * M + i];
      Ix[j * M + i] = tmp;
    }

  free(Ix0);
  free(Ix1);
}

/* compute outwards facing normals, surface Jacobian, and volume Jacobian for all face nodes */
void meshSurfaceGeometricFactorsHex3D(mesh_t *mesh)
{
  /* unified storage array for geometric factors */
  mesh->sgeo =
      (dfloat *)calloc((mesh->Nelements + mesh->totalHaloPairs) * mesh->Nsgeo * mesh->Nfp * mesh->Nfaces,
                       sizeof(dfloat));

  dfloat* xre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* xse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* xte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat* cubxre = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubxse = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubxte = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubyre = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubyse = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubyte = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubzre = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubzse = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));
  dfloat* cubzte = (dfloat*) calloc(mesh->cubNq * mesh->cubNq, sizeof(dfloat));

  for(dlong e = 0; e < mesh->Nelements + mesh->totalHaloPairs; ++e) { /* for each element */
    /* find vertex indices and physical coordinates */
    dlong id = e * mesh->Nverts;

    for(int n = 0; n < mesh->Np; ++n) {
      xre[n] = 0;
      xse[n] = 0;
      xte[n] = 0;
      yre[n] = 0;
      yse[n] = 0;
      yte[n] = 0;
      zre[n] = 0;
      zse[n] = 0;
      zte[n] = 0;
    }

    for(int k = 0; k < mesh->Nq; ++k)
      for(int j = 0; j < mesh->Nq; ++j)
        for(int i = 0; i < mesh->Nq; ++i) {
          int n = i + j * mesh->Nq + k * mesh->Nq * mesh->Nq;

          /* local node coordinates */
          dfloat rn = mesh->r[n];
          dfloat sn = mesh->s[n];
          dfloat tn = mesh->t[n];

          for(int m = 0; m < mesh->Nq; ++m) {
            int idr = e * mesh->Np + k * mesh->Nq * mesh->Nq + j * mesh->Nq + m;
            int ids = e * mesh->Np + k * mesh->Nq * mesh->Nq + m * mesh->Nq + i;
            int idt = e * mesh->Np + m * mesh->Nq * mesh->Nq + j * mesh->Nq + i;
            xre[n] += mesh->D[i * mesh->Nq + m] * mesh->x[idr];
            xse[n] += mesh->D[j * mesh->Nq + m] * mesh->x[ids];
            xte[n] += mesh->D[k * mesh->Nq + m] * mesh->x[idt];
            yre[n] += mesh->D[i * mesh->Nq + m] * mesh->y[idr];
            yse[n] += mesh->D[j * mesh->Nq + m] * mesh->y[ids];
            yte[n] += mesh->D[k * mesh->Nq + m] * mesh->y[idt];
            zre[n] += mesh->D[i * mesh->Nq + m] * mesh->z[idr];
            zse[n] += mesh->D[j * mesh->Nq + m] * mesh->z[ids];
            zte[n] += mesh->D[k * mesh->Nq + m] * mesh->z[idt];
          }
        }

    for(int f = 0; f < mesh->Nfaces; ++f) { // for each face
      for(int i = 0; i < mesh->Nfp; ++i) { // for each node on face
        /* volume index of face node */
        int n = mesh->faceNodes[f * mesh->Nfp + i];

        /* local node coordinates */
        dfloat rn = mesh->r[n];
        dfloat sn = mesh->s[n];
        dfloat tn = mesh->t[n];

        dfloat xr = xre[n], xs = xse[n], xt = xte[n];
        dfloat yr = yre[n], ys = yse[n], yt = yte[n];
        dfloat zr = zre[n], zs = zse[n], zt = zte[n];

        /* determinant of Jacobian matrix */
        dfloat J = xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) + zr * (xs * yt - ys * xt);

        dfloat rx =  (ys * zt - zs * yt) / J, ry = -(xs * zt - zs * xt) / J,
               rz =  (xs * yt - ys * xt) / J;
        dfloat sx = -(yr * zt - zr * yt) / J, sy =  (xr * zt - zr * xt) / J,
               sz = -(xr * yt - yr * xt) / J;
        dfloat tx =  (yr * zs - zr * ys) / J, ty = -(xr * zs - zr * xs) / J,
               tz =  (xr * ys - yr * xs) / J;

        /* face f normal and length */
        dfloat nx, ny, nz, d;
        switch(f) {
        case 0: nx = -tx;
          ny = -ty;
          nz = -tz;
          break;
        case 1: nx = -sx;
          ny = -sy;
          nz = -sz;
          break;
        case 2: nx = +rx;
          ny = +ry;
          nz = +rz;
          break;
        case 3: nx = +sx;
          ny = +sy;
          nz = +sz;
          break;
        case 4: nx = -rx;
          ny = -ry;
          nz = -rz;
          break;
        case 5: nx = +tx;
          ny = +ty;
          nz = +tz;
          break;
        }

        dfloat sJ = sqrt(nx * nx + ny * ny + nz * nz);
        nx /= sJ;
        ny /= sJ;
        nz /= sJ;
        sJ *= J;

        /* output index */
        dlong base = mesh->Nsgeo * (mesh->Nfaces * mesh->Nfp * e + mesh->Nfp * f + i);

        /* store normal, surface Jacobian, and reciprocal of volume Jacobian */
        mesh->sgeo[base + NXID] = nx;
        mesh->sgeo[base + NYID] = ny;
        mesh->sgeo[base + NZID] = nz;
        mesh->sgeo[base + SJID] = sJ;
        mesh->sgeo[base + IJID] = 1. / J;
        mesh->sgeo[base + WIJID] = 1. / (J * mesh->gllw[0]);
        mesh->sgeo[base + WSJID] = sJ * mesh->gllw[i % mesh->Nq] * mesh->gllw[i / mesh->Nq];

        const dfloat tol = 1e-04;
        dfloat vt1x = 0, vt1y = 0, vt1z = 0;
        dfloat vt2x = 0, vt2y = 0, vt2z = 0;
        if (std::abs(std::abs(nz) - 1.0) < tol) {
          vt1x = 1.0;
          vt1y = 0.0;
          vt1z = 0.0;
        }
        else {
          const dfloat mag = std::sqrt(nx * nx + ny * ny);
          vt1x = -ny / mag;
          vt1y = nx / mag;
          vt1z = 0.0;
        }

        mesh->sgeo[base + T1XID] = vt1x;
        mesh->sgeo[base + T1YID] = vt1y;
        mesh->sgeo[base + T1ZID] = vt1z;

        // vt2 = n \cross vt1
        vt2x = ny * vt1z - nz * vt1y;
        vt2y = nz * vt1x - nx * vt1z;
        vt2z = nx * vt1y - ny * vt1x;

        // normalize vt2
        const dfloat invMag = 1.0 / std::sqrt(vt2x * vt2x + vt2y * vt2y + vt2z * vt2z);
        vt2x *= invMag;
        vt2y *= invMag;
        vt2z *= invMag;

        mesh->sgeo[base + T2XID] = vt2x;
        mesh->sgeo[base + T2YID] = vt2y;
        mesh->sgeo[base + T2ZID] = vt2z;
      }
    }
  }

  free(xre);
  free(xse);
  free(xte);
  free(yre);
  free(yse);
  free(yte);
  free(zre);
  free(zse);
  free(zte);

  free(cubxre);
  free(cubxse);
  free(cubxte);
  free(cubyre);
  free(cubyse);
  free(cubyte);
  free(cubzre);
  free(cubzse);
  free(cubzte);
}
