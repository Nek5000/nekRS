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

#include "timer.hpp"
#include <functional>
#include "ellipticApplyMask.hpp"

#define ELLIPTIC_ENABLE_TIMER

#define NO_OP 0
#define DIRICHLET 1
#define NEUMANN 2
#define ZERO_NORMAL 3
#define ZERO_TANGENTIAL 4

class SolutionProjection;
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
  static constexpr int NScratchFields {6};
  static constexpr int minFDMBytesOverlap{2<<22};
  int dim = 1;
  int elementType = 12; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)
  int coeffField = 1;        // flag for variable coefficient (solver)
  int coeffFieldPreco = 1;   // flag for variable coefficient (preconditioner)
  int blockSolver = 0;
  int Nfields = 1;
  int stressForm = 0;
  int poisson; 

  bool mgLevel = false;

  std::string name;

  int Niter;
  dfloat res00Norm, res0Norm, resNorm;

  dlong fieldOffset; 

  mesh_t* mesh;

  void* precon = nullptr;

  ogs_t* ogs;
  oogs_t* oogs;
  oogs_t* oogsAx;

  setupAide options;

  char* type;

  dfloat tau;

  bool allNeumann;

  int* EToB;

  occa::memory o_wrk;

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
  occa::kernel axmyzManyPfloatKernel;

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
  occa::memory o_lambda;
  dlong loffset;
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

void ellipticStartHaloExchange(elliptic_t* elliptic,
                               occa::memory &o_q,
                               int Nentries,
                               dfloat* sendBuffer,
                               dfloat* recvBuffer);
void ellipticInterimHaloExchange(elliptic_t* elliptic,
                                 occa::memory &o_q,
                                 int Nentries,
                                 dfloat* sendBuffer,
                                 dfloat* recvBuffer);
void ellipticEndHaloExchange(elliptic_t* elliptic,
                             occa::memory &o_q,
                             int Nentries,
                             dfloat* recvBuffer);

//Linear solvers
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
void ellipticUpdateJacobi(elliptic_t* elliptic, occa::memory& o_invDiagA);

void ellipticMultiGridSetup(elliptic_t* elliptic, void* precon);
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
#endif
