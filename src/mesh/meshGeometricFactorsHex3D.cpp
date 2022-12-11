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
#include "mesh3D.h"
#include "platform.hpp"
#include "linAlg.hpp"

void interpolateHex3D(dfloat* I, dfloat* x, int N, dfloat* Ix, int M)
{
  dfloat* Ix1 = (dfloat*) calloc(N * N * M, sizeof(dfloat));
  dfloat* Ix2 = (dfloat*) calloc(N * M * M, sizeof(dfloat));

  for(int k = 0; k < N; ++k)
    for(int j = 0; j < N; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[i * N + n] * x[k * N * N + j * N + n];
        Ix1[k * N * M + j * M + i] = tmp;
      }

  for(int k = 0; k < N; ++k)
    for(int j = 0; j < M; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[j * N + n] * Ix1[k * N * M + n * M + i];
        Ix2[k * M * M + j * M + i] = tmp;
      }

  for(int k = 0; k < M; ++k)
    for(int j = 0; j < M; ++j)
      for(int i = 0; i < M; ++i) {
        dfloat tmp = 0;
        for(int n = 0; n < N; ++n)
          tmp += I[k * N + n] * Ix2[n * M * M + j * M + i];
        Ix[k * M * M + j * M + i] = tmp;
      }

  free(Ix1);
  free(Ix2);
}

void meshGeometricFactorsHex3D(mesh_t *mesh)
{
  double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("computing geometric factors ... "); fflush(stdout);

  /* note that we have volume geometric factors for each node */
  mesh->vgeo    = (dfloat*) calloc(mesh->Nelements * mesh->Nvgeo * mesh->Np, sizeof(dfloat));
  mesh->cubvgeo = (dfloat *)calloc(mesh->Nelements * mesh->Nvgeo * mesh->cubNp, sizeof(dfloat));
  mesh->ggeo    = (dfloat*) calloc(mesh->Nelements * mesh->Nggeo * mesh->Np,    sizeof(dfloat));

  dfloat minJ = 1e9, maxJ = -1e9, maxSkew = 0;

  dfloat* xre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* xse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* xte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* yte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zre = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zse = (dfloat*) calloc(mesh->Np, sizeof(dfloat));
  dfloat* zte = (dfloat*) calloc(mesh->Np, sizeof(dfloat));

  dfloat* cubxre = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubxse = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubxte = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubyre = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubyse = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubyte = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubzre = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubzse = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));
  dfloat* cubzte = (dfloat*) calloc(mesh->cubNp, sizeof(dfloat));

  mesh->volume = 0;

  int invalidJ = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e) {
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

    for (int k = 0; k < mesh->Nq; ++k) {
      for (int j = 0; j < mesh->Nq; ++j) {
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

          dfloat xr = xre[n], xs = xse[n], xt = xte[n];
          dfloat yr = yre[n], ys = yse[n], yt = yte[n];
          dfloat zr = zre[n], zs = zse[n], zt = zte[n];

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) + zr * (xs * yt - ys * xt);

          if (std::isnan(J) || std::isinf(J)) {
            invalidJ++;
          }

          dfloat hr = sqrt(xr * xr + yr * yr + zr * zr);
          dfloat hs = sqrt(xs * xs + ys * ys + zs * zs);
          dfloat ht = sqrt(xt * xt + yt * yt + zt * zt);
          minJ = mymin(J, minJ);
          maxJ = mymax(J, maxJ);
          maxSkew = mymax(maxSkew, hr / hs);
          maxSkew = mymax(maxSkew, hr / ht);
          maxSkew = mymax(maxSkew, hs / hr);
          maxSkew = mymax(maxSkew, hs / ht);
          maxSkew = mymax(maxSkew, ht / hr);
          maxSkew = mymax(maxSkew, ht / hs);

          //if(J<1e-12) printf("J = %g !!!!!!!!!!!!!\n", J);

          dfloat rx =  (ys * zt - zs * yt) / J, ry = -(xs * zt - zs * xt) / J,
                 rz =  (xs * yt - ys * xt) / J;
          dfloat sx = -(yr * zt - zr * yt) / J, sy =  (xr * zt - zr * xt) / J,
                 sz = -(xr * yt - yr * xt) / J;
          dfloat tx =  (yr * zs - zr * ys) / J, ty = -(xr * zs - zr * xs) / J,
                 tz =  (xr * ys - yr * xs) / J;

          dfloat JW = J * mesh->gllw[i] * mesh->gllw[j] * mesh->gllw[k];
          mesh->volume += JW;

          /* store geometric factors */
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * RXID] = rx;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * RYID] = ry;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * RZID] = rz;

          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * SXID] = sx;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * SYID] = sy;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * SZID] = sz;

          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * TXID] = tx;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * TYID] = ty;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * TZID] = tz;

          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * JID]  = J;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * JWID] = JW;
          mesh->vgeo[mesh->Nvgeo * mesh->Np * e + n + mesh->Np * IJWID] = 1. / JW;

          /* store second order geometric factors */
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G00ID] = JW *
                                                                          (rx * rx + ry * ry + rz *
                                                                           rz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G01ID] = JW *
                                                                          (rx * sx + ry * sy + rz *
                                                                           sz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G02ID] = JW *
                                                                          (rx * tx + ry * ty + rz *
                                                                           tz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G11ID] = JW *
                                                                          (sx * sx + sy * sy + sz *
                                                                           sz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G12ID] = JW *
                                                                          (sx * tx + sy * ty + sz *
                                                                           tz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * G22ID] = JW *
                                                                          (tx * tx + ty * ty + tz *
                                                                           tz);
          mesh->ggeo[mesh->Nggeo * mesh->Np * e + n + mesh->Np * GWJID] = JW;
        }
      }
    }

#if 1
    interpolateHex3D(mesh->cubInterp, xre, mesh->Nq, cubxre, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, xse, mesh->Nq, cubxse, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, xte, mesh->Nq, cubxte, mesh->cubNq);

    interpolateHex3D(mesh->cubInterp, yre, mesh->Nq, cubyre, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, yse, mesh->Nq, cubyse, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, yte, mesh->Nq, cubyte, mesh->cubNq);

    interpolateHex3D(mesh->cubInterp, zre, mesh->Nq, cubzre, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, zse, mesh->Nq, cubzse, mesh->cubNq);
    interpolateHex3D(mesh->cubInterp, zte, mesh->Nq, cubzte, mesh->cubNq);

    //geometric data for quadrature
    for (int k = 0; k < mesh->cubNq; ++k) {
      for (int j = 0; j < mesh->cubNq; ++j) {
        for(int i = 0; i < mesh->cubNq; ++i) {
          int n = k * mesh->cubNq * mesh->cubNq + j * mesh->cubNq + i;

          dfloat rn = mesh->cubr[i];
          dfloat sn = mesh->cubr[j];
          dfloat tn = mesh->cubr[k];

          /* Jacobian matrix */
          dfloat xr = cubxre[n], xs = cubxse[n], xt = cubxte[n];
          dfloat yr = cubyre[n], ys = cubyse[n], yt = cubyte[n];
          dfloat zr = cubzre[n], zs = cubzse[n], zt = cubzte[n];

          /* compute geometric factors for affine coordinate transform*/
          dfloat J = xr * (ys * zt - zs * yt) - yr * (xs * zt - zs * xt) + zr * (xs * yt - ys * xt);

          dfloat rx =  (ys * zt - zs * yt) / J, ry = -(xs * zt - zs * xt) / J,
                 rz =  (xs * yt - ys * xt) / J;
          dfloat sx = -(yr * zt - zr * yt) / J, sy =  (xr * zt - zr * xt) / J,
                 sz = -(xr * yt - yr * xt) / J;
          dfloat tx =  (yr * zs - zr * ys) / J, ty = -(xr * zs - zr * xs) / J,
                 tz =  (xr * ys - yr * xs) / J;

          dfloat JW = J * mesh->cubw[i] * mesh->cubw[j] * mesh->cubw[k];

          /* store geometric factors */
          dlong base = mesh->Nvgeo * mesh->cubNp * e + n;
          mesh->cubvgeo[base + mesh->cubNp * RXID] = rx;
          mesh->cubvgeo[base + mesh->cubNp * RYID] = ry;
          mesh->cubvgeo[base + mesh->cubNp * RZID] = rz;

          mesh->cubvgeo[base + mesh->cubNp * SXID] = sx;
          mesh->cubvgeo[base + mesh->cubNp * SYID] = sy;
          mesh->cubvgeo[base + mesh->cubNp * SZID] = sz;

          mesh->cubvgeo[base + mesh->cubNp * TXID] = tx;
          mesh->cubvgeo[base + mesh->cubNp * TYID] = ty;
          mesh->cubvgeo[base + mesh->cubNp * TZID] = tz;

          mesh->cubvgeo[base + mesh->cubNp * JID] = J;
          mesh->cubvgeo[base + mesh->cubNp * JWID] = JW;
          mesh->cubvgeo[base + mesh->cubNp * IJWID] = 1. / JW;
        }
      }
    }
#endif
  }

  MPI_Allreduce(MPI_IN_PLACE, &invalidJ, 1, MPI_INT, MPI_MAX, platform->comm.mpiComm);
  if (invalidJ) {
    if (platform->comm.mpiRank == 0) {
      std::cout << "Error: encountered nan or inf Jacobian!\n";
    }
    ABORT(1);
  }

  {
    dfloat globalMinJ = 0, globalMaxJ = 0, globalMaxSkew = 0;

    MPI_Allreduce(&minJ, &globalMinJ, 1, MPI_DFLOAT, MPI_MIN, platform->comm.mpiComm);
    MPI_Allreduce(&maxJ, &globalMaxJ, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);
    MPI_Allreduce(&maxSkew, &globalMaxSkew, 1, MPI_DFLOAT, MPI_MAX, platform->comm.mpiComm);

    if(platform->comm.mpiRank == 0)
      printf("J [%g,%g] ", globalMinJ, globalMaxJ);
      //printf("J [%g,%g] and max Skew = %g\n", globalMinJ, globalMaxJ, globalMaxSkew);

    if(globalMinJ < 0 || globalMaxJ < 0) {
      if (platform->comm.mpiRank == 0) printf("Jacobian < 0 !!! ");
      //EXIT_AND_FINALIZE(EXIT_FAILURE);
    }

    dfloat globalVolume;
    MPI_Allreduce(&mesh->volume, &globalVolume, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    mesh->volume = globalVolume;
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

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("done (%gs)\n", MPI_Wtime() - tStart); fflush(stdout);
}
