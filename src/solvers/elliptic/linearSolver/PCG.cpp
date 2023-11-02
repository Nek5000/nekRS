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

// #define DEBUG
#define TIMERS

static dfloat update(elliptic_t* elliptic,
                     const occa::memory &o_p, const occa::memory &o_Ap, const dfloat alpha,
                     occa::memory &o_x, occa::memory &o_r)
{
  mesh_t* mesh = elliptic->mesh;

  const bool serial = platform->serial;

  // r <= r - alpha*A*p
  // dot(r,r)
  elliptic->updatePCGKernel(mesh->Nlocal,
                            elliptic->fieldOffset,
                            elliptic->o_invDegree,
                            o_Ap,
                            alpha,
                            o_r,
                            elliptic->o_tmpHostScalars);

  dfloat rdotr1 = 0;
#ifdef ELLIPTIC_ENABLE_TIMER
    //platform->timer.tic("dotp",1);
#endif
  if(serial) {
    rdotr1 = *((dfloat *)elliptic->o_tmpHostScalars.ptr());
  } else {
    const dlong Nblock = (mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
    elliptic->o_tmpHostScalars.copyTo(elliptic->tmpHostScalars, Nblock);
    for(int n = 0; n < Nblock; ++n)
      rdotr1 += elliptic->tmpHostScalars[n];
  }

  // x <= x + alpha*p
  platform->linAlg->axpbyMany(
    mesh->Nlocal,
    elliptic->Nfields,
    elliptic->fieldOffset,
    alpha,
    o_p,
    1.0,
    o_x);

  MPI_Allreduce(MPI_IN_PLACE, &rdotr1, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
#ifdef ELLIPTIC_ENABLE_TIMER
    //platform->timer.toc("dotp");
#endif

  platform->flopCounter->add(elliptic->name + " ellipticUpdatePC",
                             elliptic->Nfields * static_cast<double>(mesh->Nlocal) * 6 + mesh->Nlocal);

  return rdotr1;
}

static void combinedPCGReductions(elliptic_t *elliptic,
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
#ifdef TIMERS
  platform->timer.tic("combinedPCGPostMatVec", 1);
#endif
  elliptic->combinedPCGPostMatVecKernel(mesh->Nlocal,
                                        elliptic->fieldOffset,
                                        preco,
                                        elliptic->o_invDegree,
                                        o_Minv,
                                        o_v,
                                        o_p,
                                        o_r,
                                        elliptic->o_tmpHostScalars);
#ifdef TIMERS
  platform->timer.toc("combinedPCGPostMatVec");
#endif
  if (serial) {
    auto ptr = elliptic->o_tmpHostScalars.ptr<dfloat>();
    std::copy(ptr, ptr + nRed, reductions.begin());
  } else {
    const dlong Nblock = (mesh->Nlocal + BLOCKSIZE - 1) / BLOCKSIZE;
    elliptic->o_tmpHostScalars.copyTo(elliptic->tmpHostScalars, nRed * Nblock);
    std::fill(reductions.begin(), reductions.end(), 0.0);
    for (int red = 0; red < nRed; ++red) {
      for (int n = 0; n < Nblock; ++n) {
        reductions[red] += elliptic->tmpHostScalars[n + Nblock * red];
      }
    }
  }

  // batch into single, large all-reduce
  MPI_Allreduce(MPI_IN_PLACE, reductions.data(), nRed, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
}

static int standardPCG(elliptic_t *elliptic,
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
  auto &o_p = elliptic->o_p;
  auto &o_z = (!options.compareArgs("PRECONDITIONER", "NONE")) ? elliptic->o_z : o_r;
  auto &o_Ap = elliptic->o_Ap;
  auto &o_weight = elliptic->o_invDegree;
  platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_p);

  if (platform->comm.mpiRank == 0 && verbose) {
    if (flexible) {
      printf("PFCG ");
    } else {
      printf("PCG ");
    }
    printf("%s: initial res norm %.15e WE NEED TO GET TO %e \n", elliptic->name.c_str(), rdotr, tol);
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
    alpha = rdotz1 / (pAp + 10 * std::numeric_limits<dfloat>::min());

#ifdef DEBUG
    printf("alpha: %.15e\n", alpha);
    printf("norm pAp: %.15e\n", pAp);
#endif

    //  x <= x + alpha*p
    //  r <= r - alpha*A*p
    //  dot(r,r)
    rdotr = sqrt(update(elliptic, o_p, o_Ap, alpha, o_x, o_r) * elliptic->resNormFactor);
#ifdef DEBUG
    printf("rdotr: %.15e\n", rdotr);
#endif
    if (platform->comm.mpiRank == 0) {
      nrsCheck(std::isnan(rdotr),
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
static int combinedPCG(elliptic_t *elliptic,
                       const dfloat tol,
                       const int MAXIT,
                       dfloat &rdotr,
                       occa::memory &o_r,
                       occa::memory &o_x)
{
  mesh_t* mesh = elliptic->mesh;
  setupAide& options = elliptic->options;
  auto &precon = elliptic->precon;

  const int verbose = platform->options.compareArgs("VERBOSE", "TRUE");
  const int preco = !options.compareArgs("PRECONDITIONER", "NONE");

  constexpr auto tiny = 10 * std::numeric_limits<dfloat>::min();

  dfloat betakm1 = 0;
  dfloat betakm2 = 0;
  dfloat alphak = 0;
  dfloat alphakm1 = 0;
  dfloat alphakm2 = 0;

  occa::memory o_null;

  /*aux variables */
  auto &o_v = elliptic->o_v;
  auto &o_p  = elliptic->o_p;
  auto &o_z = elliptic->o_z;
  auto &o_Minv = (preco) ? precon->o_invDiagA : o_null;
  auto &o_Ap = elliptic->o_Ap;
  auto &o_weight = elliptic->o_invDegree;
  platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_p);
  platform->linAlg->fill(elliptic->Nfields * elliptic->fieldOffset, 0.0, o_v);

  if(platform->comm.mpiRank == 0 && verbose) {
    printf("PCGC ");
    printf("%s: initial res norm %.15e WE NEED TO GET TO %e \n", elliptic->name.c_str(), rdotr, tol);
  }

  constexpr int nRed = CombinedPCGId::nReduction;

  std::array<dfloat, nRed> reductions;

  int iter = 0;
  do {

    iter++;
    const dlong updateX = iter > 1 && iter % 2 == 1;

#ifdef TIMERS
    platform->timer.tic("combinedPCGPreMatVec", 1);
#endif
    elliptic->combinedPCGPreMatVecKernel(mesh->Nlocal,
                                         updateX,
                                         preco,
                                         elliptic->fieldOffset,
                                         alphakm1,
                                         alphakm2,
                                         betakm1,
                                         betakm2,
                                         alphakm2 / betakm2,
                                         o_Minv,
                                         o_v,
                                         o_p,
                                         o_x,
                                         o_r);
#ifdef TIMERS
    platform->timer.toc("combinedPCGPreMatVec");
#endif

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

#ifdef DEBUG
    printf("alpha: %.15e\n", alphak);
    printf("norm pAp: %.15e\n", ak); // ak = p^T A p
    printf("norm rdotz1: %.15e\n", dk);
#endif

    rdotr = gammak - 2 * alphak * bk + alphak * alphak * ck;
    rdotr = sqrt(rdotr * elliptic->resNormFactor);
#ifdef DEBUG
    printf("rdotr: %.15e\n", rdotr);
#endif
    if (platform->comm.mpiRank == 0) {
      nrsCheck(std::isnan(rdotr),
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
      elliptic->combinedPCGUpdateConvergedSolutionKernel(mesh->Nlocal,
                                                         singleVectorUpdate,
                                                         preco,
                                                         elliptic->fieldOffset,
                                                         alphak,
                                                         alphakm1,
                                                         betakm1,
                                                         alphakm1 / betakm1,
                                                         o_Minv,
                                                         o_p,
                                                         o_r,
                                                         o_x);
    }

    betakm2 = betakm1;
    betakm1 = (dk - 2 * alphak * ek + alphak * alphak * fk) / (dk);

    alphakm2 = alphakm1;
    alphakm1 = alphak;
#ifdef DEBUG
    printf("beta: %.15e\n", betakm1);
#endif

  }
  while (rdotr > tol && iter < MAXIT);

  return iter;
}

int pcg(elliptic_t *elliptic,
        const dfloat tol,
        const int MAXIT,
        dfloat &rdotr,
        occa::memory &o_r,
        occa::memory &o_x)
{
  setupAide &options = elliptic->options;
  if (options.compareArgs("SOLVER", "PCG+COMBINED")) {
    return combinedPCG(elliptic, tol, MAXIT, rdotr, o_r, o_x);
  } else {
    return standardPCG(elliptic, tol, MAXIT, rdotr, o_r, o_x);
  }
}
