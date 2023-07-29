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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nrssys.hpp"
#include "mesh3D.h"
#include "platform.hpp"

#include "ellipticApplyMask.hpp"

#define ELLIPTIC_ENABLE_TIMER

#define NO_OP 0

// lower id wins
#define DIRICHLET 1
#define ZERO_NORMAL 2
#define ZERO_TANGENTIAL 3
#define NEUMANN 4

class SolutionProjection;
class precon_t;
class elliptic_t;

struct GmresData{
  GmresData(elliptic_t*);
  int nRestartVectors;
  int flexible;
  occa::memory o_V;
  occa::memory o_Z;
  occa::memory o_y;
  occa::memory o_scratch;
  occa::memory h_scratch;
  occa::memory h_y;
  dfloat* y;
  dfloat* H;
  dfloat* sn;
  dfloat* cs;
  dfloat* s;
  dfloat* scratch;
};

struct elliptic_t
{
  static constexpr double targetTimeBenchmark {0.2};

  int elementType = 12;      // number of edges (3=tri, 4=quad, 6=tet, 12=hex)
  int blockSolver = 0;
  int Nfields = 1;
  int stressForm = 0;
  int poisson = 0; 
  dlong loffset = 0;

  bool mgLevel = false;

  std::string name;

  int Niter;
  dfloat res00Norm, res0Norm, resNorm;

  dlong fieldOffset; 

  mesh_t* mesh;

  precon_t *precon = nullptr;

  ogs_t* ogs;
  oogs_t* oogs;
  oogs_t* oogsAx;

  setupAide options;

  char* type;

  dfloat tau;

  bool allNeumann;

  int* EToB;

  // C0-FEM mask data
  dlong Nmasked;
  dlong NmaskedLocal;
  dlong NmaskedGlobal;

  occa::memory o_maskIds;
  occa::memory o_maskIdsGlobal;
  occa::memory o_maskIdsLocal;

  occa::memory o_EToB;

  occa::memory o_x;
  occa::memory o_x0;
  occa::memory o_r;
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_res;
  occa::memory o_Ap; // A*search direction
  occa::memory o_v;  // work array for combined PCG iteration
  occa::memory o_invDegree;
  occa::memory o_interp;

  occa::memory o_rPfloat;
  occa::memory o_zPfloat;

  occa::memory o_EXYZ; // element vertices for reconstructing geofacs (trilinear hexes only)

  occa::kernel AxKernel;

  occa::kernel fusedCopyDfloatToPfloatKernel;

  // update for 1st Kind Chebyshev iteration
  occa::kernel updateChebyshevKernel;

  // fourth kind Chebyshev iteration
  occa::kernel updateFourthKindChebyshevKernel;

  occa::kernel updatePGMRESSolutionKernel;
  occa::kernel fusedResidualAndNormKernel;

  occa::kernel gramSchmidtOrthogonalizationKernel;

  dfloat resNormFactor;

  // combined PCG update step
  dfloat *tmpHostScalars;
  occa::memory h_tmpHostScalars;
  occa::memory o_tmpHostScalars;
  occa::kernel updatePCGKernel;

  hlong NelementsGlobal;

  occa::kernel ellipticBlockBuildDiagonalKernel;
  occa::kernel ellipticBlockBuildDiagonalPfloatKernel;

  // specialized kernels needed for combined PCG iteration
  occa::kernel combinedPCGPreMatVecKernel;
  occa::kernel combinedPCGPostMatVecKernel;
  occa::kernel combinedPCGUpdateConvergedSolutionKernel;

  occa::memory o_lambda0;
  dfloat lambda0Avg;
  occa::memory o_lambda1;

  int nLevels;
  int* levels;

  SolutionProjection* solutionProjection;
  GmresData *gmresData;

  std::function<void(dlong Nelements, const occa::memory &o_elementList, occa::memory &o_x)> applyZeroNormalMask;
  std::function<void(const occa::memory & o_r, occa::memory & o_z)> userPreconditioner;

  ~elliptic_t();
};

#include "ellipticSolutionProjection.h"

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* elliptic);

void ellipticPreconditioner(elliptic_t* elliptic, const occa::memory& o_r, occa::memory& o_z);
void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs);
void ellipticBuildPreconditionerKernels(elliptic_t* elliptic);

void ellipticSolve(elliptic_t* elliptic, const occa::memory& o_r, occa::memory o_x);

void ellipticSolveSetup(elliptic_t* elliptic);

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

int pcg(elliptic_t* elliptic, const dfloat tol, const int MAXIT, dfloat &res, occa::memory &o_r, occa::memory &o_x);

void initializeGmresData(elliptic_t*);
int pgmres(elliptic_t* elliptic, const dfloat tol, const int MAXIT, dfloat &res, occa::memory &o_r, occa::memory &o_x);

void ellipticOperator(elliptic_t* elliptic,
                      const occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char* precision,
                      bool masked = true);

void ellipticAx(elliptic_t* elliptic,
                dlong NelementsList,
                const occa::memory &o_elementsList,
                const occa::memory &o_q,
                occa::memory &o_Aq,
                const char* precision);

void ellipticMultiGridUpdateLambda(elliptic_t* elliptic);
void ellipticUpdateJacobi(elliptic_t *ellipticBase, occa::memory &o_invDiagA);
void ellipticUpdateJacobi(elliptic_t* elliptic);

void ellipticMultiGridSetup(elliptic_t *elliptic, precon_t *precon);
elliptic_t* ellipticBuildMultigridLevel(elliptic_t* baseElliptic, int Nc, int Nf);

dfloat ellipticUpdatePCG(elliptic_t* elliptic, occa::memory &o_p, occa::memory &o_Ap, dfloat alpha,
                          occa::memory &o_x, occa::memory &o_r);

void ellipticZeroMean(elliptic_t* elliptic, occa::memory &o_q);

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

static void ellipticAllocateWorkspace(elliptic_t* elliptic)
{
  elliptic->o_p = platform->o_memPool.reserve<dfloat>(elliptic->Nfields * elliptic->fieldOffset);
  elliptic->o_z = platform->o_memPool.reserve<dfloat>(elliptic->Nfields * elliptic->fieldOffset);
  elliptic->o_Ap = platform->o_memPool.reserve<dfloat>(elliptic->Nfields * elliptic->fieldOffset);
  elliptic->o_x0 = platform->o_memPool.reserve<dfloat>(elliptic->Nfields * elliptic->fieldOffset); 
  elliptic->o_rPfloat = platform->o_memPool.reserve<pfloat>(elliptic->Nfields * elliptic->fieldOffset); 
  elliptic->o_zPfloat = platform->o_memPool.reserve<pfloat>(elliptic->Nfields * elliptic->fieldOffset);

  if (elliptic->options.compareArgs("SOLVER", "PCG+COMBINED")) {
    elliptic->o_v = platform->o_memPool.reserve<dfloat>(elliptic->Nfields * elliptic->fieldOffset);
  }
}

static void ellipticFreeWorkspace(elliptic_t* elliptic)
{
  elliptic->o_p.free();
  elliptic->o_z.free(); 
  elliptic->o_Ap.free(); 
  elliptic->o_x0.free(); 
  elliptic->o_rPfloat.free();
  elliptic->o_zPfloat.free();

  if (elliptic->options.compareArgs("SOLVER", "PCG+COMBINED")) {
    elliptic->o_v.free();
  }
}

 
#endif
