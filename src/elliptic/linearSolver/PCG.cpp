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

#include <limits>
#include <array>
#include "elliptic.h"
#include "ellipticPrecon.h"
#include "timer.hpp"
#include "linAlg.hpp"

#if 0
#define DEBUG
#endif

namespace
{

constexpr auto tiny = 10 * std::numeric_limits<dfloat>::min();

occa::memory o_p;
occa::memory o_z;
occa::memory o_Ap;
occa::memory o_v;
occa::memory o_tmpReductions;
occa::memory h_tmpReductions;

dfloat update(elliptic_t *elliptic,
              const occa::memory &o_p,
              const occa::memory &o_Ap,
              const dfloat alpha,
              occa::memory &o_x,
              occa::memory &o_r)
{
  mesh_t *mesh = elliptic->mesh;

  const bool serial = platform->serial;

  // r <= r - alpha*A*p
  // dot(r,r)
  elliptic->updatePCGKernel(mesh->Nlocal,
                            elliptic->fieldOffset,
                            elliptic->o_residualWeight,
                            o_Ap,
                            alpha,
                            o_r,
                            o_tmpReductions);

  dfloat rdotr1 = 0;
#ifdef ELLIPTIC_ENABLE_TIMER
  platform->timer.tic("dotp");
#endif
  if (serial) {
    rdotr1 = *(o_tmpReductions.ptr<dfloat>());
  } else {
    auto tmp = h_tmpReductions.ptr<dfloat>();
    o_tmpReductions.copyTo(tmp);
    for (int n = 0; n < o_tmpReductions.size(); ++n) {
      rdotr1 += tmp[n];
    }
  }

  // x <= x + alpha*p
  platform->linAlg->axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, alpha, o_p, 1.0, o_x);

  MPI_Allreduce(MPI_IN_PLACE, &rdotr1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
#ifdef ELLIPTIC_ENABLE_TIMER
  platform->timer.toc("dotp");
#endif

  platform->flopCounter->add(elliptic->name + " ellipticUpdatePCG",
                             elliptic->Nfields * static_cast<double>(mesh->Nlocal) * 6 + mesh->Nlocal);

  return rdotr1;
}

void combinedPCGReductions(elliptic_t *elliptic,
                           const dlong preco,
                           const occa::memory &o_Minv,
                           const occa::memory &o_v,
                           const occa::memory &o_p,
                           const occa::memory &o_r,
                           std::array<dfloat, CombinedPCGId::nReduction> &reductions)
{
  constexpr auto nRed = CombinedPCGId::nReduction;
  auto mesh = elliptic->mesh;
  const bool serial = platform->serial;
  elliptic->combinedPCGPostMatVecKernel(mesh->Nlocal,
                                        elliptic->fieldOffset,
                                        preco,
                                        elliptic->o_residualWeight,
                                        elliptic->o_invDegree,
                                        o_Minv,
                                        o_v,
                                        o_p,
                                        o_r,
                                        o_tmpReductions);
  if (serial) {
    auto ptr = o_tmpReductions.ptr<dfloat>();
    std::copy(ptr, ptr + nRed, reductions.begin());
  } else {
    auto tmp = h_tmpReductions.ptr<dfloat>();
    o_tmpReductions.copyTo(tmp);
    std::fill(reductions.begin(), reductions.end(), 0.0);

    auto mesh = elliptic->mesh;
    const dlong Nlocal = mesh->Np * mesh->Nelements;
    const dlong Nblock = (Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int red = 0; red < nRed; ++red) {
      for (int n = 0; n < Nblock; ++n) {
        reductions[red] += tmp[n + Nblock * red];
      }
    }
  }

  // batch into single fused all-reduce
  MPI_Allreduce(MPI_IN_PLACE, reductions.data(), nRed, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);

  platform->flopCounter->add(elliptic->name + " ellipticCombinedPCGReductions",
                             elliptic->Nfields * static_cast<double>(mesh->Nlocal) * 3 * 7);
}

int standardPCG(elliptic_t *elliptic,
                const dfloat tol,
                const int MAXIT,
                dfloat &rdotr,
                occa::memory &o_r,
                occa::memory &o_x)
{

  mesh_t *mesh = elliptic->mesh;
  setupAide &options = elliptic->options;

  const int flexible = options.compareArgs("SOLVER", "FLEXIBLE");
  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");

  dfloat rdotz1;
  dfloat alpha;

  /*aux variables */
  auto &o_weight = elliptic->o_invDegree;
  platform->linAlg->fill(o_p.size(), 0.0, o_p);

  if (platform->comm.mpiRank == 0 && verbose) {
    if (flexible) {
      printf("PFCG ");
    } else {
      printf("PCG ");
    }
    printf("%s: initial res norm %.15e target %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  int iter = 0;
  do {
    iter++;
    const dfloat rdotz2 = rdotz1;
    if (!options.compareArgs("PRECONDITIONER", "NONE")) {
      ellipticPreconditioner(elliptic, o_r, o_z);

      rdotz1 = platform->linAlg->weightedInnerProdMany(mesh->Nlocal,
                                                       elliptic->Nfields,
                                                       elliptic->fieldOffset,
                                                       o_weight,
                                                       o_r,
                                                       o_z,
                                                       platform->comm.mpiComm);
    } else {
      rdotz1 = rdotr;
    }

    if (platform->comm.mpiRank == 0) {
      nekrsCheck(std::isnan(rdotz1),
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s\n",
                 "Detected invalid rdotz norm while running linear solver!");
    }

#ifdef DEBUG
    printf("norm rdotz1: %.15e\n", rdotz1);
#endif

    dfloat beta = 0;
    if (iter > 1) {
      beta = rdotz1 / rdotz2;
      if (flexible) {
        const dfloat zdotAp = platform->linAlg->weightedInnerProdMany(mesh->Nlocal,
                                                                      elliptic->Nfields,
                                                                      elliptic->fieldOffset,
                                                                      o_weight,
                                                                      o_z,
                                                                      o_Ap,
                                                                      platform->comm.mpiComm);
        beta = -alpha * zdotAp / rdotz2;
#ifdef DEBUG
        printf("norm zdotAp: %.15e\n", zdotAp);
#endif
      }
    }

#ifdef DEBUG
    printf("beta: %.15e\n", beta);
#endif

    platform->linAlg->axpbyMany(mesh->Nlocal, elliptic->Nfields, elliptic->fieldOffset, 1.0, o_z, beta, o_p);

    ellipticOperator(elliptic, o_p, o_Ap, dfloatString);
    const dfloat pAp = platform->linAlg->weightedInnerProdMany(mesh->Nlocal,
                                                               elliptic->Nfields,
                                                               elliptic->fieldOffset,
                                                               o_weight,
                                                               o_p,
                                                               o_Ap,
                                                               platform->comm.mpiComm);
    alpha = rdotz1 / (pAp + tiny);

#ifdef DEBUG
    printf("alpha: %.15e\n", alpha);
    printf("norm pAp: %.15e\n", pAp);
#endif

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr = std::sqrt(update(elliptic, o_p, o_Ap, alpha, o_x, o_r));
#ifdef DEBUG
    printf("rdotr: %.15e\n", rdotr);
#endif
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
  } while (rdotr > tol && iter < MAXIT);

  return iter;
}

// Algo 5 from https://arxiv.org/pdf/2205.08909.pdf
int combinedPCG(elliptic_t *elliptic,
                const dfloat tol,
                const int MAXIT,
                dfloat &rdotr,
                occa::memory &o_r,
                occa::memory &o_x)
{
  mesh_t *mesh = elliptic->mesh;
  setupAide &options = elliptic->options;
  auto &precon = elliptic->precon;

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int preco = !options.compareArgs("PRECONDITIONER", "NONE");


  dfloat betakm1 = 0;
  dfloat betakm2 = 0;
  dfloat alphak = 0;
  dfloat alphakm1 = 0;
  dfloat alphakm2 = 0;

  occa::memory o_null;

  /*aux variables */
  auto &o_Minv = (preco) ? precon->o_invDiagA : o_null;
  platform->linAlg->fill(o_p.size(), 0.0, o_p);
  platform->linAlg->fill(o_v.size(), 0.0, o_v);

  if (platform->comm.mpiRank == 0 && verbose) {
    printf("PCGC ");
    printf("%s: initial res norm %.15e target %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  constexpr int nRed = CombinedPCGId::nReduction;

  std::array<dfloat, nRed> reductions;

  int iter = 0;
  do {

    iter++;
    const dlong updateX = iter > 1 && iter % 2 == 1;

    elliptic->combinedPCGPreMatVecKernel(mesh->Nlocal,
                                         updateX,
                                         preco,
                                         elliptic->fieldOffset,
                                         alphakm1,
                                         alphakm2,
                                         betakm1,
                                         betakm2,
                                         (updateX) ? alphakm2 / betakm2 : static_cast<dfloat>(0),
                                         o_Minv,
                                         o_v,
                                         o_p,
                                         o_x,
                                         o_r);

    platform->flopCounter->add(elliptic->name + " ellipticCombinedPCGPreMatVecKernel",
                               elliptic->Nfields * static_cast<double>(mesh->Nlocal) * 0.5*(11 + 5));


    ellipticOperator(elliptic, o_p, o_v, dfloatString);

    combinedPCGReductions(elliptic, preco, o_Minv, o_v, o_p, o_r, reductions);

    const auto gammak = reductions[CombinedPCGId::gamma];
    const auto ak = reductions[CombinedPCGId::a];
    const auto bk = reductions[CombinedPCGId::b];
    const auto ck = reductions[CombinedPCGId::c];
    const auto dk = reductions[CombinedPCGId::d];
    const auto ek = reductions[CombinedPCGId::e];
    const auto fk = reductions[CombinedPCGId::f];

    alphak = dk / (ak + tiny);

    if (platform->comm.mpiRank == 0) {
      nekrsCheck(std::isnan(dk),
                 MPI_COMM_SELF,
                 EXIT_FAILURE,
                 "%s\n",
                 "Detected invalid rdotz norm while running linear solver!");
    }

#ifdef DEBUG
    printf("alpha: %.15e\n", alphak);
    printf("norm pAp: %.15e\n", ak); // ak = p^T A p
    printf("norm rdotz1: %.15e\n", dk);
#endif

    // r_{k+1}^T r_{k+1} = (r - alpha v)^T (r - alpha v)
    rdotr = gammak + alphak * (-2. * bk + alphak * ck);
    rdotr = std::sqrt(std::abs(rdotr));
#ifdef DEBUG
    printf("rdotr: %.15e\n", rdotr);
#endif
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

    // converged, update solution prior to exit
    if (rdotr <= tol) {
      const dlong singleVectorUpdate = iter % 2 == 1;
      if (platform->comm.mpiRank == 0) {
        nekrsCheck(!singleVectorUpdate && betakm1 == 0,
                   MPI_COMM_SELF,
                   EXIT_FAILURE,
                   "%s\n",
                   "Cannot update solution as beta == 0!");
      }
      const dfloat alphaInvBeta = (!singleVectorUpdate) ? alphakm1 / betakm1 : 0;

      elliptic->combinedPCGUpdateConvergedSolutionKernel(mesh->Nlocal,
                                                         singleVectorUpdate,
                                                         preco,
                                                         elliptic->fieldOffset,
                                                         alphak,
                                                         alphakm1,
                                                         betakm1,
                                                         alphaInvBeta,
                                                         o_Minv,
                                                         o_p,
                                                         o_r,
                                                         o_x);
    }

    betakm2 = betakm1;
    betakm1 = std::abs(1 + alphak * (-2. * ek + alphak * fk) / dk);

#ifdef DEBUG
    printf("beta: %.15e\n", betakm1);
#endif

    alphakm2 = alphakm1;
    alphakm1 = alphak;

  } while (rdotr > tol && iter < MAXIT);

  return iter;
}

} // namespace

int pcg(elliptic_t *elliptic,
        const dfloat tol,
        const int MAXIT,
        dfloat &rdotr,
        occa::memory &o_r,
        occa::memory &o_x)
{
  setupAide &options = elliptic->options;

  const auto Nlocal = (elliptic->Nfields > 1) ? elliptic->Nfields * static_cast<size_t>(elliptic->fieldOffset)
                                              : elliptic->mesh->Nlocal;

  o_p = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  o_z = (elliptic->options.compareArgs("PRECONDITIONER", "NONE"))
            ? o_r
            : platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  o_Ap = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  if (elliptic->options.compareArgs("SOLVER", "PCG+COMBINED")) {
    o_v = platform->deviceMemoryPool.reserve<dfloat>(Nlocal);
  }

  const auto Nblock = [&]()
  {
    auto mesh = elliptic->mesh;
    const dlong Nlocal = mesh->Np * mesh->Nelements;
    return (Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
  }();


  auto Nreductions = [&]() {
    int n = 1;
    if (options.compareArgs("SOLVER", "PCG+COMBINED")) n = CombinedPCGId::nReduction; 
    return n;
  }();

  h_tmpReductions = platform->memoryPool.reserve<dfloat>(Nreductions * Nblock);
  o_tmpReductions = platform->deviceMemoryPool.reserve<dfloat>(h_tmpReductions.size());

  const auto Niter = [&]() {
    if (elliptic->options.compareArgs("SOLVER", "PCG+COMBINED")) {
      return combinedPCG(elliptic, tol, MAXIT, rdotr, o_r, o_x);
    } else {
      return standardPCG(elliptic, tol, MAXIT, rdotr, o_r, o_x);
    }
  }();

  o_p.free();
  if (o_z != o_r) {
    o_z.free();
  }
  o_Ap.free();
  if (elliptic->options.compareArgs("SOLVER", "PCG+COMBINED")) {
    o_v.free();
  }

  o_tmpReductions.free();
  h_tmpReductions.free();

  return Niter;
}
