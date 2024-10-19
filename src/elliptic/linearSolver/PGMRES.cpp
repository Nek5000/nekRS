#include "platform.hpp"
#include "elliptic.h"
#include "timer.hpp"
#include "linAlg.hpp"

namespace
{

int _pgmres(elliptic_t *elliptic,
            const dfloat tol,
            const int MAXIT,
            dfloat &rdotr,
            occa::memory &o_r,
            occa::memory &o_x)
{
  auto tiny = 10 * std::numeric_limits<dfloat>::min();

  auto mesh = elliptic->mesh;
  linAlg_t &linAlg = *(platform->linAlg);

  const auto offset = elliptic->fieldOffset * elliptic->Nfields;
  const int Nblock = (mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;

  const auto nRestartVectors = elliptic->gmresData->nRestartVectors;
  const auto verbose = platform->verbose;
  const auto serial = platform->serial;
  const auto flexible = elliptic->gmresData->flexible;

  auto &o_tmp = elliptic->gmresData->o_p;
  auto &o_w = elliptic->gmresData->o_Ap;
  auto &o_r0 = elliptic->gmresData->o_z;

  auto &o_V = elliptic->gmresData->o_V;
  auto &o_Z = elliptic->gmresData->o_Z;

  auto &o_y = elliptic->gmresData->o_y;

  auto y = static_cast<dfloat *>(elliptic->gmresData->_y.ptr());

  auto &H = elliptic->gmresData->H;
  auto &sn = elliptic->gmresData->sn;
  auto &cs = elliptic->gmresData->cs;
  auto &s = elliptic->gmresData->s;

  auto &o_weight = elliptic->o_residualWeight;

  o_r0.copyFrom(o_r, o_r.size());

  dfloat nr = rdotr;

  if (verbose && (platform->comm.mpiRank == 0)) {
    if (flexible) {
      printf("PFGMRES ");
    } else {
      printf("PGMRES ");
    }
    printf("%s: initial res norm %.15e target %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  auto update = [&](int gmresUpdateSize) {
    for (int k = gmresUpdateSize - 1; k >= 0; --k) {
      y[k] = s[k];

      for (int m = k + 1; m < gmresUpdateSize; ++m) {
        y[k] -= H[k + m * (nRestartVectors + 1)] * y[m];
      }

      y[k] /= (H[k + k * (nRestartVectors + 1)] + tiny);
    }

    o_y.copyFrom(y, gmresUpdateSize);

    if (flexible) {
      elliptic
          ->updatePGMRESSolutionKernel(mesh->Nlocal, elliptic->fieldOffset, gmresUpdateSize, o_y, o_Z, o_x);
    } else {
      platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_tmp);
      elliptic
          ->updatePGMRESSolutionKernel(mesh->Nlocal, elliptic->fieldOffset, gmresUpdateSize, o_y, o_V, o_tmp);

      auto &o_Mtmp = o_w;
      ellipticPreconditioner(elliptic, o_tmp, o_Mtmp);
      platform->linAlg
          ->axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, 1.0, o_Mtmp, 1.0, o_x);

      double flopCount = 2 * gmresUpdateSize * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
      platform->flopCounter->add("gmresUpdate", flopCount);
    }
  };

  int iter;
  for (iter = 0; iter < MAXIT;) {
    s[0] = nr;

    // normalize
    linAlg.axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, 1. / (nr + tiny), o_r, 0.0, o_V);

    for (int i = 0; i < nRestartVectors; ++i) {

      auto o_Mv = flexible ? o_Z + i * offset : o_Z;
      ellipticPreconditioner(elliptic, o_V + i * offset, o_Mv);

      ellipticOperator(elliptic, static_cast<const occa::memory>(o_Mv), o_w, dfloatString);

      // 1 pass classical Gram-Schmidt
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
      o_y.copyTo(y, (i + 1));
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
      } else {
        auto scratch = static_cast<dfloat *>(elliptic->gmresData->_scratch.ptr());
        elliptic->gmresData->o_scratch.copyTo(scratch, Nblock);
        for (int k = 0; k < Nblock; ++k) {
          nw += scratch[k];
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &nw, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
      nw = sqrt(nw);

      {
        double flopCount = 5 * (i + 1) * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
        platform->flopCounter->add("gramSchmidt", flopCount);
      }

      // apply Givens rotations to new column
      //  H(i+1,i) = ||w||_2
      H[i + 1 + i * (nRestartVectors + 1)] = nw;

      // nromalize
      if (i < nRestartVectors - 1) {
        auto o_Vi = o_V + (i + 1) * offset;
        linAlg.axpbyMany(mesh->Nlocal,
                         elliptic->Nfields,
                         elliptic->fieldOffset,
                         1. / (nw + tiny),
                         o_w,
                         0,
                         o_Vi);
      }

      // apply Givens rotation
      for (int k = 0; k <= i; ++k) {
        H[k + i * (nRestartVectors + 1)] = y[k];
      }

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
      rdotr = fabs(s[i + 1]);

      if (platform->comm.mpiRank == 0) {
        nekrsCheck(std::isnan(rdotr),
                   MPI_COMM_SELF,
                   EXIT_FAILURE,
                   "%s\n",
                   "Detected invalid resiual norm while running linear solver!");
      }

      if (verbose && (platform->comm.mpiRank == 0)) {
        printf("it %d r norm %.15e\n", iter, rdotr);
      }

      if (rdotr < tol || iter == MAXIT) {
        // update approximation
        update(i + 1);
        break;
      }
    }

    // exit if tolerance is reached
    if (rdotr < tol || iter == MAXIT) {
      break;
    }

    // update approximation
    update(nRestartVectors);

    ellipticOperator(elliptic, o_x, o_tmp, dfloatString);

    elliptic->fusedResidualAndNormKernel(Nblock,
                                         mesh->Nlocal,
                                         elliptic->fieldOffset,
                                         o_weight,
                                         o_r0,
                                         o_tmp,
                                         o_r,
                                         elliptic->gmresData->o_scratch);

    if (serial) {
      nr = *((dfloat *)elliptic->gmresData->o_scratch.ptr());
    } else {
      auto scratch = static_cast<dfloat *>(elliptic->gmresData->_scratch.ptr());
      elliptic->gmresData->o_scratch.copyTo(scratch, Nblock);
      nr = 0.0;
      for (dlong n = 0; n < Nblock; ++n) {
        nr += scratch[n];
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &nr, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    nr = sqrt(nr);

    {
      double flopCount = 4 * elliptic->Nfields * static_cast<double>(mesh->Nlocal);
      platform->flopCounter->add("gmres evaluate residual and norm", flopCount);
    }

    rdotr = nr;

    if (rdotr <= tol) {
      return iter;
    }
  }

  return iter;
}

} // namespace

void initializeGmresData(elliptic_t *elliptic)
{
  elliptic->gmresData = new GmresData(elliptic);
}

GmresData::GmresData(elliptic_t *elliptic)
{
  nRestartVectors = 15;
  if (elliptic->options.getArgs("PGMRES RESTART").empty()) {
    elliptic->options.setArgs("PGMRES RESTART", std::to_string(nRestartVectors));
  } else {
    elliptic->options.getArgs("PGMRES RESTART", nRestartVectors);
  }

  flexible = elliptic->options.compareArgs("SOLVER", "FLEXIBLE");

  H.resize((nRestartVectors + 1) * (nRestartVectors + 1));
  sn.resize(nRestartVectors);
  cs.resize(nRestartVectors);
  s.resize(nRestartVectors + 1);

  auto N = [&]() {
    auto Nblock = (elliptic->mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
    return nRestartVectors * Nblock;
  }();

  o_scratch = platform->device.malloc<dfloat>(N);
  _scratch = platform->device.mallocHost<dfloat>(o_scratch.size());

  o_y = platform->device.malloc<dfloat>(nRestartVectors);
  _y = platform->device.mallocHost<dfloat>(o_y.size());
}

int pgmres(elliptic_t *elliptic,
           const dfloat tol,
           const int MAXIT,
           dfloat &rdotr,
           occa::memory &o_r,
           occa::memory &o_x)
{
  const auto Nlocal = elliptic->Nfields * static_cast<size_t>(elliptic->fieldOffset);

  elliptic->gmresData->o_p = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  elliptic->gmresData->o_z = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  elliptic->gmresData->o_Ap = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);

  elliptic->gmresData->o_V =
      platform->deviceMemoryPool.reserve<dfloat>(Nlocal * elliptic->gmresData->nRestartVectors);
  elliptic->gmresData->o_Z = platform->deviceMemoryPool.reserve<dfloat>(
      Nlocal * ((elliptic->gmresData->flexible) ? elliptic->gmresData->nRestartVectors : 1));

  const int Niter = _pgmres(elliptic, tol, MAXIT, rdotr, o_r, o_x);

  elliptic->gmresData->o_p.free();
  elliptic->gmresData->o_z.free();
  elliptic->gmresData->o_Ap.free();
  elliptic->gmresData->o_V.free();
  elliptic->gmresData->o_Z.free();

  return Niter;
}
