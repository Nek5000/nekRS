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

#include "elliptic.h"

void BuildLocalContinuousDiagHex3D (elliptic_t* elliptic,
                                    mesh_t* mesh,
                                    dlong eM,
                                    dfloat* B,
                                    dfloat* Br,
                                    dfloat* Bs,
                                    dfloat* Bt,
                                    dfloat* A);
void BuildLocalContinuousBlockDiagHex3D (elliptic_t* elliptic,
                                         mesh_t* mesh,
                                         dfloat* B,
                                         dfloat* Br,
                                         dfloat* Bs,
                                         dfloat* Bt,
                                         dfloat* A);

void ellipticUpdateJacobi(elliptic_t* elliptic)
{
  mesh_t* mesh       = elliptic->mesh;
  setupAide options  = elliptic->options;
  precon_t* precon   = elliptic->precon;

  const dfloat allNeumannScale = elliptic->allNeumannPenalty * elliptic->allNeumannScale *
                                 elliptic->allNeumannScale;
  const dlong Nlocal = mesh->Np * mesh->Nelements;

  elliptic->updateDiagonalKernel(mesh->Nelements,
                                 elliptic->Ntotal,
                                 elliptic->loffset,
                                 elliptic->allNeumann,
                                 allNeumannScale,
                                 elliptic->o_mapB,
                                 mesh->o_ggeo,
                                 mesh->o_Dmatrices,
                                 mesh->o_Smatrices,
                                 elliptic->o_lambda,
                                 precon->o_invDiagA);

  oogs::startFinish(precon->o_invDiagA, elliptic->Nfields, elliptic->Ntotal, ogsDfloat, ogsAdd, elliptic->oogs);

  const dfloat one = 1.0;
  if(elliptic->blockSolver)
    elliptic->scalarDivideManyKernel(Nlocal, elliptic->Ntotal, one, precon->o_invDiagA);
  else
    elliptic->scalarDivideKernel(Nlocal, one, precon->o_invDiagA);
}

void ellipticBuildJacobi(elliptic_t* elliptic, dfloat** invDiagA)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  MPI_Barrier(mesh->comm);
  const double tStart = MPI_Wtime();
  if(mesh->rank == 0) printf("building Jacobi ... ");
  fflush(stdout);

  // surface mass matrices MS = MM*LIFT
  dfloat* MS = (dfloat*) calloc(mesh->Nfaces * mesh->Nfp * mesh->Nfp,sizeof(dfloat));
  for (int f = 0; f < mesh->Nfaces; f++)
    for (int n = 0; n < mesh->Nfp; n++) {
      int fn = mesh->faceNodes[f * mesh->Nfp + n];

      for (int m = 0; m < mesh->Nfp; m++) {
        dfloat MSnm = 0;

        for (int i = 0; i < mesh->Np; i++)
          MSnm += mesh->MM[fn + i * mesh->Np] *
                  mesh->LIFT[i * mesh->Nfp * mesh->Nfaces + f * mesh->Nfp + m];

        MS[m + n * mesh->Nfp + f * mesh->Nfp * mesh->Nfp]  = MSnm;
      }
    }

  // build some monolithic basis arrays (for quads and hexes)
  dfloat* B  = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Br = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Bs = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));
  dfloat* Bt = (dfloat*) calloc(mesh->Np * mesh->Np, sizeof(dfloat));

  if (elliptic->elementType == QUADRILATERALS) {
    int mode = 0;
    for(int nj = 0; nj < mesh->N + 1; ++nj)
      for(int ni = 0; ni < mesh->N + 1; ++ni) {
        int node = 0;

        for(int j = 0; j < mesh->N + 1; ++j)
          for(int i = 0; i < mesh->N + 1; ++i) {
            if(nj == j && ni == i)
              B[mode * mesh->Np + node] = 1;
            if(nj == j)
              Br[mode * mesh->Np + node] = mesh->D[ni + mesh->Nq * i];
            if(ni == i)
              Bs[mode * mesh->Np + node] = mesh->D[nj + mesh->Nq * j];

            ++node;
          }
        ++mode;
      }
  }

  if (elliptic->elementType == HEXAHEDRA) {
    int mode = 0;
    for(int nk = 0; nk < mesh->N + 1; ++nk)
      for(int nj = 0; nj < mesh->N + 1; ++nj)
        for(int ni = 0; ni < mesh->N + 1; ++ni) {
          int node = 0;

          for(int k = 0; k < mesh->N + 1; ++k)
            for(int j = 0; j < mesh->N + 1; ++j)
              for(int i = 0; i < mesh->N + 1; ++i) {
                if(nk == k && nj == j && ni == i)
                  B[mode * mesh->Np + node] = 1;
                if(nj == j && nk == k)
                  Br[mode * mesh->Np + node] = mesh->D[ni + mesh->Nq * i];
                if(ni == i && nk == k)
                  Bs[mode * mesh->Np + node] = mesh->D[nj + mesh->Nq * j];
                if(ni == i && nj == j)
                  Bt[mode * mesh->Np + node] = mesh->D[nk + mesh->Nq * k];

                ++node;
              }

          ++mode;
        }
  }

  dlong diagNnum;
  if(elliptic->blockSolver)
    diagNnum = elliptic->Ntotal * elliptic->Nfields;
  else
    diagNnum = mesh->Np * mesh->Nelements;

  dfloat* diagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));

  switch(elliptic->elementType) {
  case HEXAHEDRA:
    if(elliptic->blockSolver) {
      BuildLocalContinuousBlockDiagHex3D(elliptic, mesh, B, Br, Bs, Bt, diagA);
      break;
    }else{
#pragma omp parallel for
      for(dlong eM = 0; eM < mesh->Nelements; ++eM)
        BuildLocalContinuousDiagHex3D(elliptic, mesh, eM, B, Br, Bs, Bt, diagA + eM * mesh->Np);
      break;
    }
  }

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    if(elliptic->blockSolver) {
      ogsGatherScatterMany(diagA,
                           elliptic->Nfields,
                           elliptic->Ntotal,
                           ogsDfloat,
                           ogsAdd,
                           elliptic->ogs);
      *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
      for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
        (*invDiagA)[n] = 1;

      for(int fld = 0; fld < elliptic->Nfields; fld++)
        for(int n = 0; n < mesh->Nelements * mesh->Np; n++)
          (*invDiagA)[n + fld * elliptic->Ntotal] = 1 / diagA[n + fld * elliptic->Ntotal];
    }else{
      ogsGatherScatter(diagA, ogsDfloat, ogsAdd, elliptic->ogs);
      *invDiagA = (dfloat*) calloc(diagNnum, sizeof(dfloat));
      for (dlong n = 0; n < mesh->Nelements * mesh->Np; n++)
        (*invDiagA)[n] = 1 / diagA[n];
    }
  }

  MPI_Barrier(mesh->comm);
  if(mesh->rank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);

  free(diagA);
  free(MS);
  free(B);
  free(Br);
  free(Bs);
  free(Bt);
}

void BuildLocalContinuousDiagHex3D(elliptic_t* elliptic,
                                   mesh_t* mesh,
                                   dlong eM,
                                   dfloat* B,
                                   dfloat* Br,
                                   dfloat* Bs,
                                   dfloat* Bt,
                                   dfloat* A)
{
  int var_coeff = elliptic->var_coeff;

  for (int nz = 0; nz < mesh->Nq; nz++)
    for (int ny = 0; ny < mesh->Nq; ny++)
      for (int nx = 0; nx < mesh->Nq; nx++) {
        int idn = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
        if (elliptic->mapB[idn + eM * mesh->Np] != 1) {
          A[idn] = 0;

          int id = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
          dlong base = eM * mesh->Np * mesh->Nggeo;

          dfloat lambda_0 =
            var_coeff ? elliptic->lambda[eM * mesh->Np + id + 0 * elliptic->Ntotal] : 1.0;
          dfloat lambda_1 =
            var_coeff ? elliptic->lambda[eM * mesh->Np + id + 1 *
                                         elliptic->Ntotal] : elliptic->lambda[0];

          dfloat Grs = mesh->ggeo[base + id + G01ID * mesh->Np];
          A[idn] += 2 * lambda_0 * Grs * mesh->D[nx + nx * mesh->Nq] * mesh->D[ny + ny * mesh->Nq];

          dfloat Grt = mesh->ggeo[base + id + G02ID * mesh->Np];
          A[idn] += 2 * lambda_0 * Grt * mesh->D[nx + nx * mesh->Nq] * mesh->D[nz + nz * mesh->Nq];

          dfloat Gst = mesh->ggeo[base + id + G12ID * mesh->Np];
          A[idn] += 2 * lambda_0 * Gst * mesh->D[ny + ny * mesh->Nq] * mesh->D[nz + nz * mesh->Nq];

          for (int k = 0; k < mesh->Nq; k++) {
            int iid = k + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
            dfloat Grr = mesh->ggeo[base + iid + G00ID * mesh->Np];
            lambda_0   =
              var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal] : 1.0;
            A[idn] += Grr * lambda_0 * mesh->D[nx + k * mesh->Nq] * mesh->D[nx + k * mesh->Nq];
          }

          for (int k = 0; k < mesh->Nq; k++) {
            int iid = nx + k * mesh->Nq + nz * mesh->Nq * mesh->Nq;
            dfloat Gss = mesh->ggeo[base + iid + G11ID * mesh->Np];
            lambda_0   =
              var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal] : 1.0;
            A[idn] += Gss * lambda_0 * mesh->D[ny + k * mesh->Nq] * mesh->D[ny + k * mesh->Nq];
          }

          for (int k = 0; k < mesh->Nq; k++) {
            int iid = nx + ny * mesh->Nq + k * mesh->Nq * mesh->Nq;
            dfloat Gtt = mesh->ggeo[base + iid + G22ID * mesh->Np];
            lambda_0   =
              var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal] : 1.0;
            A[idn] += Gtt * lambda_0 * mesh->D[nz + k * mesh->Nq] * mesh->D[nz + k * mesh->Nq];
          }

          dfloat JW = mesh->ggeo[base + id + GWJID * mesh->Np];
          A[idn] += JW * lambda_1;
        } else {
          A[idn] = 1; //just put a 1 so A is invertable
        }
      }

  //add the rank boost for the allNeumann Poisson problem
  if (elliptic->allNeumann)
    for(int n = 0; n < mesh->Np; ++n)
      if (elliptic->mapB[n + eM * mesh->Np] != 1) //dont fill rows for masked nodes
        A[n] += elliptic->allNeumannPenalty * elliptic->allNeumannScale * elliptic->allNeumannScale;
}

void BuildLocalContinuousBlockDiagHex3D(elliptic_t* elliptic,
                                        mesh_t* mesh,
                                        dfloat* B,
                                        dfloat* Br,
                                        dfloat* Bs,
                                        dfloat* Bt,
                                        dfloat* A)
{
  int var_coeff = elliptic->var_coeff;
#pragma omp parallel for
  for(dlong eM = 0; eM < mesh->Nelements; ++eM) {
    for(int fld = 0; fld < elliptic->Nfields; fld++) {
      const dlong offset = fld * elliptic->Ntotal;

      for (int nz = 0; nz < mesh->Nq; nz++)
        for (int ny = 0; ny < mesh->Nq; ny++)
          for (int nx = 0; nx < mesh->Nq; nx++) {
            int idn = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq + eM * mesh->Np + offset;
            if (elliptic->mapB[idn] != 1) {
              A[idn] = 0;

              int id     = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq + eM * mesh->Np;
              int gid    = nx + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
              dlong base = eM * mesh->Np * mesh->Nggeo;

              dfloat lambda_0 =
                var_coeff ? elliptic->lambda[id + 0 * elliptic->Ntotal + fld *
                                             elliptic->loffset] : 1.0;
              dfloat lambda_1 =
                var_coeff ? elliptic->lambda[id + 1 * elliptic->Ntotal + fld *
                                             elliptic->loffset] : elliptic->lambda[fld *
                                                                                   elliptic->loffset
                ];

              dfloat Grs = mesh->ggeo[base + gid + G01ID * mesh->Np];
              A[idn] += 2 * lambda_0 * Grs * mesh->D[nx + nx * mesh->Nq] *
                        mesh->D[ny + ny * mesh->Nq];

              dfloat Grt = mesh->ggeo[base + gid + G02ID * mesh->Np];
              A[idn] += 2 * lambda_0 * Grt * mesh->D[nx + nx * mesh->Nq] *
                        mesh->D[nz + nz * mesh->Nq];

              dfloat Gst = mesh->ggeo[base + gid + G12ID * mesh->Np];
              A[idn] += 2 * lambda_0 * Gst * mesh->D[ny + ny * mesh->Nq] *
                        mesh->D[nz + nz * mesh->Nq];

              for (int k = 0; k < mesh->Nq; k++) {
                int iid  = k + ny * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                dfloat Grr = mesh->ggeo[base + iid + G00ID * mesh->Np];
                lambda_0   =
                  var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal + fld *
                                               elliptic->loffset] : 1.0;
                A[idn] += Grr * lambda_0 * mesh->D[nx + k * mesh->Nq] * mesh->D[nx + k * mesh->Nq];
              }

              for (int k = 0; k < mesh->Nq; k++) {
                int iid = nx + k * mesh->Nq + nz * mesh->Nq * mesh->Nq;
                dfloat Gss = mesh->ggeo[base + iid + G11ID * mesh->Np];
                lambda_0   =
                  var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal + fld *
                                               elliptic->loffset] : 1.0;
                A[idn] += Gss * lambda_0 * mesh->D[ny + k * mesh->Nq] * mesh->D[ny + k * mesh->Nq];
              }

              for (int k = 0; k < mesh->Nq; k++) {
                int iid = nx + ny * mesh->Nq + k * mesh->Nq * mesh->Nq;
                dfloat Gtt = mesh->ggeo[base + iid + G22ID * mesh->Np];
                lambda_0   =
                  var_coeff ? elliptic->lambda[eM * mesh->Np + iid + 0 * elliptic->Ntotal + fld *
                                               elliptic->loffset] : 1.0;
                A[idn] += Gtt * lambda_0 * mesh->D[nz + k * mesh->Nq] * mesh->D[nz + k * mesh->Nq];
              }

              dfloat JW = mesh->ggeo[base + gid + GWJID * mesh->Np];
              A[idn] += JW * lambda_1;
            } else {
              A[idn] = 1; //just put a 1 so A is invertable
            }
          }
    }

    //add the rank boost for the allNeumann Poisson problem
    for(int fld = 0; fld < elliptic->Nfields; fld++) {
      const dlong offset = fld * elliptic->Ntotal;
      if (elliptic->allBlockNeumann[fld])
        for(int n = 0; n < mesh->Np; ++n)
          if (elliptic->mapB[n + eM * mesh->Np + offset] != 1) //dont fill rows for masked nodes
            A[n + eM * mesh->Np + offset] += elliptic->allNeumannPenalty *
                                             elliptic->allNeumannScale * elliptic->allNeumannScale;
    }
  }
}
