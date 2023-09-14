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
#include "timer.hpp"
#include "linAlg.hpp"

static dfloat tiny;

GmresData::GmresData(elliptic_t *elliptic)
    : nRestartVectors([&]() {
        int _nRestartVectors = 15;
        elliptic->options.getArgs("PGMRES RESTART", _nRestartVectors);
        return _nRestartVectors;
      }()),
      o_y(platform->device.malloc<dfloat>(nRestartVectors)),
      H((dfloat *)calloc((nRestartVectors + 1) * (nRestartVectors + 1), sizeof(dfloat))),
      sn((dfloat *)calloc(nRestartVectors, sizeof(dfloat))),
      cs((dfloat *)calloc(nRestartVectors, sizeof(dfloat))),
      s((dfloat *)calloc(nRestartVectors + 1, sizeof(dfloat)))
{
  tiny = 10*std::numeric_limits<dfloat>::min();

  int Nblock = (elliptic->mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
  const size_t N = nRestartVectors * Nblock;
  {
    h_scratch = platform->device.mallocHost<dfloat>(N);
    scratch = (dfloat *)h_scratch.ptr();

    h_y = platform->device.mallocHost<dfloat>(nRestartVectors);
    y = (dfloat *)h_y.ptr();
  }
  o_scratch = platform->device.malloc<dfloat>(N);
}

void initializeGmresData(elliptic_t *elliptic)
{
  GmresData *gmresData = new GmresData(elliptic);
  elliptic->gmresData = gmresData;
}

namespace {
void gmresUpdate(elliptic_t *elliptic, occa::memory o_x, int gmresUpdateSize)
{
  const int nRestartVectors = elliptic->gmresData->nRestartVectors;
  mesh_t *mesh = elliptic->mesh;
  dfloat *y = elliptic->gmresData->y;
  dfloat *H = elliptic->gmresData->H;
  dfloat *s = elliptic->gmresData->s;
  auto &o_V = elliptic->gmresData->o_V;
  auto &o_Z = elliptic->gmresData->o_Z;
  auto &o_y = elliptic->gmresData->o_y;
  auto &o_z = elliptic->o_z;
  auto &o_tmp = elliptic->o_p;

  for (int k = gmresUpdateSize - 1; k >= 0; --k) {
    y[k] = s[k];

    for (int m = k + 1; m < gmresUpdateSize; ++m) {
      y[k] -= H[k + m * (nRestartVectors + 1)] * y[m];
     }

    y[k] /= (H[k + k * (nRestartVectors + 1)] + tiny);
  }

  o_y.copyFrom(y, gmresUpdateSize);

  if (elliptic->options.compareArgs("SOLVER", "FLEXIBLE")) {
    elliptic->updatePGMRESSolutionKernel(mesh->Nlocal, elliptic->fieldOffset, gmresUpdateSize, o_y, o_Z, o_x);
  }
  else {
    platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_z);
    elliptic->updatePGMRESSolutionKernel(mesh->Nlocal, elliptic->fieldOffset, gmresUpdateSize, o_y, o_V, o_z);

    ellipticPreconditioner(elliptic, o_z, o_tmp);
    platform->linAlg->axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, 1.0, o_tmp, 1.0, o_x);

    double flopCount = 2 * gmresUpdateSize * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
    platform->flopCounter->add("gmresUpdate", flopCount);
  }
}
} // namespace

void free(elliptic_t *elliptic)
{ 
  elliptic->gmresData->o_V.free();
  elliptic->gmresData->o_Z.free();
}

// Ax=r
int pgmres(elliptic_t *elliptic,
           const dfloat tol,
           const int MAXIT,
           dfloat &rdotr,
           occa::memory &o_r,
           occa::memory &o_x)
{

  mesh_t *mesh = elliptic->mesh;
  linAlg_t &linAlg = *(platform->linAlg);

  const int nRestartVectors = elliptic->gmresData->nRestartVectors;
  const bool verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const bool serial = platform->device.mode() == "Serial" || platform->device.mode() == "OpenMP";
  const int flexible = elliptic->options.compareArgs("SOLVER", "FLEXIBLE");

  elliptic->gmresData->o_V = platform->o_memPool.reserve<dfloat>(static_cast<size_t>(elliptic->fieldOffset) * elliptic->Nfields * nRestartVectors);
  elliptic->gmresData->o_Z = platform->o_memPool.reserve<dfloat>(static_cast<size_t>(elliptic->fieldOffset) * elliptic->Nfields * ((flexible) ? nRestartVectors : 1));

  auto &o_w = elliptic->o_p;

  auto &o_Ax = elliptic->o_Ap;

  auto &o_V = elliptic->gmresData->o_V;
  auto &o_Z = elliptic->gmresData->o_Z;

  auto &o_y = elliptic->gmresData->o_y;
  auto &o_weight = elliptic->o_invDegree;

  auto &o_b = elliptic->o_z;
  o_b.copyFrom(o_r, elliptic->fieldOffset * elliptic->Nfields);

  auto y = elliptic->gmresData->y;
  auto H = elliptic->gmresData->H;
  auto sn = elliptic->gmresData->sn;
  auto cs = elliptic->gmresData->cs;
  auto s = elliptic->gmresData->s;

  const auto offset = elliptic->fieldOffset * elliptic->Nfields;
  const int Nblock = (mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;

  dfloat nr = rdotr / sqrt(elliptic->resNormFactor);
  dfloat error = rdotr;
  const dfloat TOL = tol;


  if (verbose && (platform->comm.mpiRank == 0)) {
    if (flexible)
      printf("PFGMRES ");
    else
      printf("PGMRES ");
    printf("%s: initial res norm %.15e WE NEED TO GET TO %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  int iter = 0;

  for (iter = 0; iter < MAXIT;) {

    s[0] = nr;

    // V(:,0) = r/nr
    linAlg.axpbyMany(mesh->Nlocal, 
                     elliptic->Nfields, 
                     elliptic->fieldOffset, 
                     1. / (nr + tiny), 
                     o_r, 
                     0.0, 
                     o_V);

    // Construct orthonormal basis via Gram-Schmidt
    for (int i = 0; i < nRestartVectors; ++i) {

      auto o_Mv = flexible ? o_Z + i*offset : o_Z;
      // z := M^{-1} V(:,i)
      ellipticPreconditioner(elliptic, o_V + i*offset, o_Mv);

      // w := A z
      ellipticOperator(elliptic, static_cast<const occa::memory>(o_Mv), o_w, dfloatString);

      // classical Gram-Schmidt
#if USE_WEIGHTED_INNER_PROD_MULTI_DEVICE
      linAlg.weightedInnerProdMulti(mesh->Nlocal,
                                    (i + 1),
                                    elliptic->Nfields,
                                    elliptic->fieldOffset,
                                    o_weight,
                                    o_V,
                                    o_w,
                                    platform->comm.mpiComm,
                                    o_y);
      o_y.copyTo(y, (i + 1) );
#else
      linAlg.weightedInnerProdMulti(mesh->Nlocal,
                                    (i + 1),
                                    elliptic->Nfields,
                                    elliptic->fieldOffset,
                                    o_weight,
                                    o_V,
                                    o_w,
                                    platform->comm.mpiComm,
                                    y);
      o_y.copyFrom(y, (i + 1));
#endif

      elliptic->gramSchmidtOrthogonalizationKernel(Nblock,
                                                   mesh->Nlocal,
                                                   elliptic->fieldOffset,
                                                   (i + 1),
                                                   o_weight,
                                                   o_y,
                                                   o_V,
                                                   o_w,
                                                   elliptic->gmresData->o_scratch);
      dfloat nw = 0.0;

      if (serial) {
        nw = *((dfloat *)elliptic->gmresData->o_scratch.ptr());
      }
      else {
        elliptic->gmresData->o_scratch.copyTo(elliptic->gmresData->scratch, Nblock);
        for (int k = 0; k < Nblock; ++k)
          nw += elliptic->gmresData->scratch[k];
      }
      MPI_Allreduce(MPI_IN_PLACE, &nw, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
      nw = sqrt(nw);

      {
        double flopCount = 5 * (i + 1) * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
        platform->flopCounter->add("gramSchmidt", flopCount);
      }

      //apply Givens rotations to new column
      // H(i+1,i) = ||w||_2
      H[i + 1 + i * (nRestartVectors + 1)] = nw;

      // V(:,i+1) = w/nw
      if (i < nRestartVectors - 1) {
        auto o_Vi = o_V + (i + 1)*offset;
        linAlg.axpbyMany(mesh->Nlocal,
                         elliptic->Nfields,
                         elliptic->fieldOffset,
                         1. / (nw + tiny),
                         o_w,
                         0,
                         o_Vi);
      }

      // apply Givens rotation
      for (int k = 0; k <= i; ++k)
        H[k + i * (nRestartVectors + 1)] = y[k];

      for (int k = 0; k < i; ++k) {
        const dfloat h1 = H[k + i * (nRestartVectors + 1)];
        const dfloat h2 = H[k + 1 + i * (nRestartVectors + 1)];

        H[k + i * (nRestartVectors + 1)] = cs[k] * h1 + sn[k] * h2;
        H[k + 1 + i * (nRestartVectors + 1)] = -sn[k] * h1 + cs[k] * h2;
      }

      // form i-th rotation matrix
      const dfloat h1 = H[i + i * (nRestartVectors + 1)];
      const dfloat h2 = H[i + 1 + i * (nRestartVectors + 1)];
      const dfloat hr = sqrt(h1 * h1 + h2 * h2) + tiny;
      cs[i] = h1 / hr;
      sn[i] = h2 / hr;

      H[i + i * (nRestartVectors + 1)] = cs[i] * h1 + sn[i] * h2;
      H[i + 1 + i * (nRestartVectors + 1)] = 0;

      // approximate residual norm
      s[i + 1] = -sn[i] * s[i];
      s[i] = cs[i] * s[i];

      iter++;
      error = fabs(s[i + 1]) * sqrt(elliptic->resNormFactor);
      rdotr = error;

      if (platform->comm.mpiRank == 0)
        nrsCheck(std::isnan(error),
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s\n",
                 "Detected invalid resiual norm while running linear solver!");

      if (verbose && (platform->comm.mpiRank == 0))
        printf("it %d r norm %.15e\n", iter, rdotr);

      if (error < TOL || iter == MAXIT) {
        // update approximation
        gmresUpdate(elliptic, o_x, i + 1);
        break;
      }
    }

    // exit if tolerance is reached
    if (error < TOL || iter == MAXIT)
      break;

    // update approximation
    gmresUpdate(elliptic, o_x, nRestartVectors);

    // nRestartVectors GMRES
    // compute A*x
    ellipticOperator(elliptic, o_x, o_Ax, dfloatString);

    elliptic->fusedResidualAndNormKernel(Nblock,
                                         mesh->Nlocal,
                                         elliptic->fieldOffset,
                                         elliptic->o_invDegree,
                                         o_b,
                                         o_Ax,
                                         o_r,
                                         elliptic->gmresData->o_scratch);

    if (serial) {
      nr = *((dfloat *)elliptic->gmresData->o_scratch.ptr());
    }
    else {
      nr = 0.0;
      elliptic->gmresData->o_scratch.copyTo(elliptic->gmresData->scratch, Nblock);
      for (dlong n = 0; n < Nblock; ++n)
        nr += elliptic->gmresData->scratch[n];
    }

    MPI_Allreduce(MPI_IN_PLACE, &nr, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    nr = sqrt(nr);

    {
      double flopCount = 4 * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
      platform->flopCounter->add("gmres evaluate residual and norm", flopCount);
    }

    error = nr * sqrt(elliptic->resNormFactor);
    rdotr = nr * sqrt(elliptic->resNormFactor);
    // exit if tolerance is reached
    if (error <= TOL) {
      free(elliptic);
      return iter;
    }
  }

  free(elliptic);
  return iter;
}
