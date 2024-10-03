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

#ifndef ELLIPTIC_H
#define ELLIPTIC_H 1

#include "platform.hpp"
#include "mesh3D.h"

#include "ellipticApplyMask.hpp"
#include "ellipticBcTypes.h"

// #define ELLIPTIC_ENABLE_TIMER

class SolutionProjection;
class elliptic_t;
class precon_t;

struct GmresData {
  GmresData(elliptic_t *);
  int nRestartVectors;
  int flexible;

  occa::memory o_V;
  occa::memory o_Z;

  std::vector<dfloat> H;
  std::vector<dfloat> sn;
  std::vector<dfloat> cs;
  std::vector<dfloat> s;

  occa::memory o_scratch;
  occa::memory _scratch;

  occa::memory o_y;
  occa::memory _y;

  occa::memory o_p;
  occa::memory o_z;
  occa::memory o_Ap;
};

struct elliptic_t {
  int elementType = 12; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)
  int blockSolver = 0;
  int Nfields = 1;
  int stressForm = 0;
  int poisson = 0;

  const dlong loffset = 0; // same operator coeffs for all components

  bool mgLevel = false;

  std::string name = "unknown";

  int Niter;
  dfloat res00Norm, res0Norm, resNorm;

  dlong fieldOffset = -1;

  mesh_t *mesh = nullptr;

  precon_t *precon = nullptr;

  ogs_t *ogs;
  oogs_t *oogs;
  oogs_t *oogsAx;

  setupAide options;

  bool nullspace = 0;

  int *EToB = nullptr;

  // C0-FEM mask data
  dlong Nmasked;
  dlong NmaskedLocal;
  dlong NmaskedGlobal;

  occa::memory o_residualWeight;

  occa::memory o_maskIds;
  occa::memory o_maskIdsGlobal;
  occa::memory o_maskIdsLocal;

  occa::memory o_EToB;

  occa::memory o_invDegree;
  occa::memory o_interp;

  occa::memory o_rPfloat;
  occa::memory o_zPfloat;

  occa::kernel AxKernel;

  occa::kernel fusedCopyDfloatToPfloatKernel;

  // update for 1st Kind Chebyshev iteration
  occa::kernel updateChebyshevKernel;

  // fourth kind Chebyshev iteration
  occa::kernel updateFourthKindChebyshevKernel;

  occa::kernel updatePGMRESSolutionKernel;
  occa::kernel fusedResidualAndNormKernel;

  occa::kernel gramSchmidtOrthogonalizationKernel;

  occa::kernel updatePCGKernel;

  occa::kernel ellipticBlockBuildDiagonalKernel;
  occa::kernel ellipticBlockBuildDiagonalPfloatKernel;

  // specialized kernels needed for combined PCG iteration
  occa::kernel combinedPCGPreMatVecKernel;
  occa::kernel combinedPCGPostMatVecKernel;
  occa::kernel combinedPCGUpdateConvergedSolutionKernel;

  occa::memory o_lambda0;
  dfloat lambda0Avg;
  occa::memory o_lambda1;

  int nLevels = 0;
  int *levels = nullptr;

  SolutionProjection *solutionProjection = nullptr;
  GmresData *gmresData = nullptr;

  std::function<void(dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x)>
      applyZeroNormalMask;
  std::function<void(const occa::memory &o_r, occa::memory &o_z)> userPreconditioner;

  ~elliptic_t();
};

#include "ellipticSolutionProjection.hpp"
#include "MG/ellipticMultiGrid.h"

void ellipticMultiGridUpdateLambda(elliptic_t *elliptic);
void ellipticMultiGridSetup(elliptic_t *elliptic);

elliptic_t *ellipticBuildMultigridLevelFine(elliptic_t *elliptic);
elliptic_t *ellipticBuildMultigridLevel(elliptic_t *baseElliptic, int Nc, int Nf);

void ellipticPreconditioner(elliptic_t *elliptic, const occa::memory &o_r, occa::memory &o_z);
void ellipticPreconditionerSetup(elliptic_t *elliptic, ogs_t *ogs);
void ellipticBuildMultigridLevelKernels(elliptic_t *elliptic);

void ellipticSolve(elliptic_t *elliptic,
                   const occa::memory &o_lambda0,
                   const occa::memory &o_lambda1,
                   const occa::memory &o_r,
                   occa::memory o_x);

void ellipticSolveSetup(elliptic_t *elliptic, const occa::memory &o_lambda0, const occa::memory &o_lambda1);

// indices used in combinedPCG routines
struct CombinedPCGId {
  static constexpr int nReduction = 7;
  static constexpr int gamma = 0;
  static constexpr int a = 1;
  static constexpr int b = 2;
  static constexpr int c = 3;
  static constexpr int d = 4;
  static constexpr int e = 5;
  static constexpr int f = 6;
};

int pcg(elliptic_t *elliptic,
        const dfloat tol,
        const int MAXIT,
        dfloat &res,
        occa::memory &o_r,
        occa::memory &o_x);

void initializeGmresData(elliptic_t *);
int pgmres(elliptic_t *elliptic,
           const dfloat tol,
           const int MAXIT,
           dfloat &res,
           occa::memory &o_r,
           occa::memory &o_x);

void ellipticOperator(elliptic_t *elliptic,
                      const occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char *precision,
                      bool masked = true);

void ellipticAx(elliptic_t *elliptic,
                dlong NelementsList,
                const occa::memory &o_elementsList,
                const occa::memory &o_q,
                occa::memory &o_Aq,
                const char *precision);

void ellipticUpdateJacobi(elliptic_t *ellipticBase, occa::memory &o_invDiagA);
void ellipticUpdateJacobi(elliptic_t *elliptic);

dfloat ellipticUpdatePCG(elliptic_t *elliptic,
                         occa::memory &o_p,
                         occa::memory &o_Ap,
                         dfloat alpha,
                         occa::memory &o_x,
                         occa::memory &o_r);

void ellipticZeroMean(elliptic_t *elliptic, occa::memory &o_q);

void ellipticOgs(mesh_t *mesh,
                 dlong mNlocal,
                 int nFields,
                 dlong offset,
                 int *EToB,
                 dlong &Nmasked,
                 occa::memory &o_maskIds,
                 dlong &NmaskedLocal,
                 occa::memory &o_maskIdsLocal,
                 dlong &NmaskedGlobal,
                 occa::memory &o_maskIdsGlobal,
                 ogs_t **ogs);

void ellipticAllocateWorkspace(elliptic_t *elliptic);
void ellipticFreeWorkspace(elliptic_t *elliptic);

#endif
