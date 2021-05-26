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
#include "amgSolver/parAlmond/parAlmond.hpp"
#include "ellipticPrecon.h"
#include "platform.hpp"

#include "timer.hpp"

#define ELLIPTIC_ENABLE_TIMER

class ResidualProjection;
class elliptic_t;

struct GmresData{
  GmresData(elliptic_t*);
  int restart;
  int flexible;
  deviceVector_t o_V;
  deviceVector_t o_Z;
  occa::memory o_y;
  occa::memory o_scratch;
  occa::memory h_scratch;
  dfloat* y;
  dfloat* H;
  dfloat* sn;
  dfloat* cs;
  dfloat* s;
  dfloat* scratch;
};

struct elliptic_t
{
  int dim;
  int elementType; // number of edges (3=tri, 4=quad, 6=tet, 12=hex)
  int var_coeff;   // flag for variable coefficient
  int blockSolver, Nfields, stressForm; // flag for vector solver and number of fields

  string name;

  int Niter;
  dfloat res00Norm, res0Norm, resNorm;

  dlong Ntotal; // offset

  mesh_t* mesh;

  precon_t* precon;

  ogs_t* ogs;
  oogs_t* oogs;
  oogs_t* oogsAx;

  setupAide options;

  char* type;

  dfloat tau;

  int* BCType;
  int NBCType;

  int* allBlockNeumann;
  bool allNeumann;
  dfloat allNeumannPenalty;
  dfloat allNeumannScale;

  // HOST shadow copies
  dfloat* p, * z, * v, * Ap;
  dfloat* invDegree;

  int* EToB;

  dfloat* wrk;
  occa::memory o_wrk;

  //C0-FEM mask data
  int* mapB;      // boundary flag of face nodes
  dlong Nmasked;
  dlong* fNmasked;

  dlong* maskIds;
  hlong* maskedGlobalIds;

  occa::memory o_maskIds;
  occa::memory o_mapB;

  occa::stream defaultStream;
  occa::stream dataStream;

  occa::memory o_x;
  occa::memory o_x0;
  occa::memory o_r;
  occa::memory o_p; // search direction
  occa::memory o_z; // preconditioner solution
  occa::memory o_res;
  occa::memory o_Ap; // A*search direction
  occa::memory o_rtmp;
  occa::memory o_invDegree;
  occa::memory o_EToB;

  occa::memory o_EXYZ; // element vertices for reconstructing geofacs (trilinear hexes only)
  occa::memory o_gllzw; // GLL nodes and weights

  occa::kernel AxKernel;
  occa::kernel AxStressKernel;
  occa::kernel AxPfloatKernel;
  occa::kernel partialAxKernel;
  occa::kernel partialAxKernel2;
  occa::kernel partialAxPfloatKernel;
  occa::kernel partialCubatureAxKernel;

  occa::kernel rhsBCKernel;
  occa::kernel addBCKernel;
  occa::kernel scaledAddPfloatKernel;
  occa::kernel dotMultiplyPfloatKernel;
  occa::kernel copyDfloatToPfloatKernel;
  occa::kernel fusedCopyDfloatToPfloatKernel;
  occa::kernel copyPfloatToDPfloatKernel;
  
  // special kernels for single Chebyshev iteration
  occa::kernel updateSmoothedSolutionVecKernel;
  occa::kernel updateChebyshevSolutionVecKernel;

  // special kernel for two Chebyshev iterations
  occa::kernel updateIntermediateSolutionVecKernel;

  occa::kernel updatePGMRESSolutionKernel;
  occa::kernel fusedResidualAndNormKernel;

  occa::kernel gramSchmidtOrthogonalizationKernel;

  // SEMFEM kernels
  occa::kernel gatherKernel;
  occa::kernel scatterKernel;
  occa::memory o_dofMap;
  occa::memory o_SEMFEMBuffer1;
  occa::memory o_SEMFEMBuffer2;
  dlong numRowsSEMFEM;

  dfloat resNormFactor;

  // combined PCG update step
  dfloat* tmpNormr;
  occa::memory o_tmpNormr;
  occa::kernel updatePCGKernel;

  hlong NelementsGlobal;

  occa::kernel updateDiagonalKernel;
  occa::memory o_lambda;
  dfloat* lambda;
  dlong loffset;
  int nLevels;
  int* levels;

  ResidualProjection* residualProjection;
  GmresData* gmresData;
};

#include "ellipticMultiGrid.h"
#include "ellipticResidualProjection.h"

elliptic_t* ellipticBuildMultigridLevelFine(elliptic_t* elliptic);

void ellipticPreconditioner(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_z);
void ellipticPreconditionerSetup(elliptic_t* elliptic, ogs_t* ogs, occa::properties &kernelInfo);

void ellipticSEMFEMSetup(elliptic_t*);
void ellipticSEMFEMSolve(elliptic_t*, occa::memory&, occa::memory&);

void ellipticSolve(elliptic_t* elliptic, occa::memory &o_r, occa::memory &o_x);

void ellipticSolveSetup(elliptic_t* elliptic, occa::properties kernelInfo);
void ellipticBlockSolveSetup(elliptic_t* elliptic, occa::properties &kernelInfo);

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
                      const char* precision);

void ellipticAx(elliptic_t* elliptic,
                dlong NelementsList,
                occa::memory &o_elementsList,
                occa::memory &o_q,
                occa::memory &o_Aq,
                const char* precision);

void ellipticBuildContinuous(elliptic_t* elliptic, nonZero_t** A,
                             dlong* nnz, ogs_t** ogs, hlong* globalStarts);

void ellipticBuildContinuousGalerkinHex3D(elliptic_t* elliptic,
                                          elliptic_t* ellipticFine,
                                          dfloat lambda,
                                          nonZero_t** A,
                                          dlong* nnz,
                                          ogs_t** ogs,
                                          hlong* globalStarts);

void ellipticBuildJacobi(elliptic_t* elliptic, dfloat** invDiagA);
void ellipticUpdateJacobi(elliptic_t* elliptic);

void ellipticBuildLocalPatches(elliptic_t* elliptic, dfloat lambda, dfloat rateTolerance,
                               dlong* Npataches, dlong** patchesIndex, dfloat** patchesInvA);

void ellipticMultiGridSetup(elliptic_t* elliptic, precon_t* precon);
elliptic_t* ellipticBuildMultigridLevel(elliptic_t* baseElliptic, int Nc, int Nf);

dfloat ellipticUpdatePCG(elliptic_t* elliptic, occa::memory &o_p, occa::memory &o_Ap, dfloat alpha,
                          occa::memory &o_x, occa::memory &o_r);

occa::properties ellipticKernelInfo(mesh_t* mesh);

void ellipticZeroMean(elliptic_t* elliptic, occa::memory &o_q);

#endif
