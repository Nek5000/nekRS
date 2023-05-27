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
  deviceVector_t o_V;
  deviceVector_t o_Z;
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
  static constexpr int NWorkspaceFields {6};

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

  static occa::memory o_wrk;

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
  occa::memory o_invDegree;
  occa::memory o_interp;

  occa::memory o_rPfloat;
  occa::memory o_zPfloat;

  occa::memory o_EXYZ; // element vertices for reconstructing geofacs (trilinear hexes only)

  occa::kernel AxKernel;
  occa::kernel AxPfloatKernel;

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
  dfloat* tmpNormr;
  occa::memory o_tmpNormr;
  occa::kernel updatePCGKernel;

  hlong NelementsGlobal;

  occa::kernel ellipticBlockBuildDiagonalKernel;
  occa::kernel ellipticBlockBuildDiagonalPfloatKernel;

  occa::memory o_lambda0;
  occa::memory o_lambda1;

  int nLevels;
  int* levels;

  SolutionProjection* solutionProjection;
  GmresData *gmresData;

  std::function<void(dlong Nelements, occa::memory &o_elementList, occa::memory &o_x)> applyZeroNormalMask;
  std::function<void(occa::memory & o_r, occa::memory & o_z)> userPreconditioner;

  ~elliptic_t();
};

#include "ellipticSolutionProjection.h"

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* elliptic);

void ellipticPreconditioner(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_z);
void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs);
void ellipticBuildPreconditionerKernels(elliptic_t* elliptic);

void ellipticSolve(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x);

void ellipticSolveSetup(elliptic_t* elliptic);

int pcg(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT, dfloat &res);

void initializeGmresData(elliptic_t*);
int pgmres(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x,
        const dfloat tol, const int MAXIT, dfloat &res);

void ellipticOperator(elliptic_t* elliptic,
                      occa::memory &o_q,
                      occa::memory &o_Aq,
                      const char* precision,
                      bool masked = true);

void ellipticAx(elliptic_t* elliptic,
                dlong NelementsList,
                occa::memory &o_elementsList,
                occa::memory &o_q,
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

static void ellipticUpdateWorkspace(elliptic_t* elliptic)
{
  const auto Nfields = 6; // first 6 slices are reserved as input to ellipticSolve
  const auto offsetBytesWrk = elliptic->fieldOffset * (Nfields * sizeof(dfloat));
  const auto offsetBytes = elliptic->fieldOffset * (elliptic->Nfields * sizeof(dfloat));

  const auto requiredBytes = offsetBytesWrk + elliptic_t::NWorkspaceFields * offsetBytes;
  nrsCheck(platform->o_mempool.o_ptr.size() < requiredBytes,
           MPI_COMM_SELF,
           EXIT_FAILURE,
           "platform mempool too small! (required %ld out of %ld bytes)\n", 
           requiredBytes, platform->o_mempool.o_ptr.size());

  elliptic_t::o_wrk = 
    platform->o_mempool.o_ptr.slice(offsetBytesWrk);

  elliptic->o_p = elliptic->o_wrk + 0 * offsetBytes;
  elliptic->o_z = elliptic->o_wrk + 1 * offsetBytes;
  elliptic->o_Ap = elliptic->o_wrk + 2 * offsetBytes;
  elliptic->o_x0 = elliptic->o_wrk + 3 * offsetBytes;
  elliptic->o_rPfloat = elliptic->o_wrk + 4 * offsetBytes;
  elliptic->o_zPfloat = elliptic->o_wrk + 5 * offsetBytes;
}
 
#endif
