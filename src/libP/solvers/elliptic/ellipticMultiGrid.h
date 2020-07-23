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

#ifndef ELLIPTIC_MGLEVEL_HPP
#define ELLIPTIC_MGLEVEL_HPP

typedef enum {RICHARDSON = 1,
              CHEBYSHEV = 2,
              SCHWARZ = 3} SmoothType;
typedef enum {JACOBI = 1,
              LOCALPATCH = 2,
              SCHWARZ_SMOOTH = 3} SmootherType;

class MGLevel : public parAlmond::multigridLevel
{
public:

  elliptic_t* elliptic;
  mesh_t* mesh;
  dfloat lambda;

  int degree;

  //coarsener
  dfloat* R;
  occa::memory o_R;
  int NpF;
  occa::memory o_invDegree;

  //smoothing params
  SmoothType stype;
  SmootherType smtypeUp;
  SmootherType smtypeDown;

  dfloat lambda1, lambda0;
  int ChebyshevIterations;

  static size_t smootherResidualBytes;
  static dfloat* smootherResidual;
  static occa::memory o_smootherResidual;
  static occa::memory o_smootherResidual2;
  static occa::memory o_smootherUpdate;
  occa::kernel preFDMKernel;
  occa::kernel fusedFDMKernel;
  occa::kernel postFDMKernel;
  occa::kernel collocateKernel;
  // Eigenvectors
  occa::memory o_Sx;
  occa::memory o_Sy;
  occa::memory o_Sz;

  occa::memory o_work1; // scratch space
  occa::memory o_work2; // scratch space
  occa::memory o_work3; // scratch space
  occa::memory o_wts; // weights to apply after operation

  // Eigenvalues
  occa::memory o_invL;

  //jacobi data
  occa::memory o_invDiagA;

  //local patch data
  occa::memory o_invAP, o_patchesIndex, o_invDegreeAP;
  //ogs_t* extendedOgs = nullptr;
  void* extendedOgs;
  void* ogs;
  void build(
    elliptic_t* pSolver);
  void generate_weights();

  setupAide options;

  //build a single level
  MGLevel(elliptic_t* ellipticBase, dfloat lambda_, int Nc,
          setupAide options_, parAlmond::KrylovType ktype_, MPI_Comm comm_
          );
  //build a level and connect it to the previous one
  MGLevel(elliptic_t* ellipticBase, //finest level
          mesh_t** meshLevels,
          elliptic_t* ellipticFine,          //previous level
          elliptic_t* ellipticCoarse,          //current level
          dfloat lambda_,
          int Nf, int Nc,
          setupAide options_,
          parAlmond::KrylovType ktype_,
          MPI_Comm comm_);

  void Ax(dfloat* x, dfloat* Ax) {}
  void Ax(occa::memory o_x, occa::memory o_Ax);

  void residual(dfloat* rhs, dfloat* x, dfloat* res) {}
  void residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res);

  void coarsen(dfloat* x, dfloat* Cx) {}
  void coarsen(occa::memory o_x, occa::memory o_Cx);

  void prolongate(dfloat* x, dfloat* Px) {}
  void prolongate(occa::memory o_x, occa::memory o_Px);

  //smoother ops
  void smooth(dfloat* rhs, dfloat* x, bool x_is_zero) {}
  void smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero);

  void smoother(occa::memory o_x, occa::memory o_Sx, bool xIsZero);

  void smoothRichardson(occa::memory &o_r, occa::memory &o_x, bool xIsZero);
  void smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero);
  void smoothSchwarz (occa::memory &o_r, occa::memory &o_x, bool xIsZero);

  void smootherLocalPatch(occa::memory &o_r, occa::memory &o_Sr);
  void smootherJacobi    (occa::memory &o_r, occa::memory &o_Sr);

  void Report();

  void setupSmoother(elliptic_t* base);
  dfloat maxEigSmoothAx();

  void buildCoarsenerTriTet(mesh_t** meshLevels, int Nf, int Nc);
  void buildCoarsenerQuadHex(mesh_t** meshLevels, int Nf, int Nc);
};

void MGLevelAllocateStorage(MGLevel* level, int k, parAlmond::CycleType ctype);

#endif
