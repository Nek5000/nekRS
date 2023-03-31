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

#include "elliptic.h"
#include "MG/MGSolver.hpp"
#include <vector>

enum class SmootherType
{
  CHEBYSHEV,
  FOURTH_CHEBYSHEV,
  OPT_FOURTH_CHEBYSHEV,
  ASM,
  RAS,
  JACOBI,
};
enum class ChebyshevSmootherType
{  
  NONE,
  JACOBI,
  ASM,
  RAS,
};

std::vector<pfloat> optimalCoeffs(int ChebyshevDegree);

class pMGLevel : public MGSolver_t::multigridLevel
{
public:
  static constexpr hlong Narnoldi {10};

  elliptic_t* elliptic;
  mesh_t* mesh;

  int degree;

  //coarsener
  dfloat* R;
  occa::memory o_R;
  int NpF;
  occa::memory o_invDegreeFine;

  //smoothing params
  SmootherType smootherType;
  ChebyshevSmootherType chebySmootherType = ChebyshevSmootherType::NONE;

  dfloat lambda1, lambda0;
  dfloat maxEig;

  int DownLegChebyshevDegree;
  int UpLegChebyshevDegree;

  inline static pfloat* smootherResidual;
  inline static occa::memory o_smootherResidual;
  inline static occa::memory o_smootherResidual2;
  inline static occa::memory o_smootherUpdate;
  occa::kernel preFDMKernel;
  occa::kernel fusedFDMKernel;
  occa::kernel postFDMKernel;
  // Eigenvectors
  occa::memory o_Sx;
  occa::memory o_Sy;
  occa::memory o_Sz;

  occa::memory o_work1; // scratch space
  occa::memory o_work2; // scratch space
  occa::memory o_wts; // weights to apply after operation

  // Eigenvalues
  occa::memory o_invL;

  //jacobi data
  occa::memory o_invDiagA;

  //local patch data
  occa::memory o_invAP, o_patchesIndex, o_invDegreeAP;
  void* ogs;
  void* ogsOverlap;

  void* ogsExt;
  void* ogsExtOverlap;

  void build(
    elliptic_t* pSolver);
  void generate_weights();

  setupAide options;

  bool isCoarse;

  std::vector<pfloat> UpLegBetas;
  std::vector<pfloat> DownLegBetas;

  //build a single level
  pMGLevel(elliptic_t* ellipticBase, int Nc,
          setupAide options_, MPI_Comm comm_,
          bool _isCoarse = false
          );
  //build a level and connect it to the previous one
  pMGLevel(elliptic_t* ellipticBase, //finest level
          mesh_t** meshLevels,
          elliptic_t* ellipticFine,          //previous level
          elliptic_t* ellipticCoarse,          //current level
          int Nf, int Nc,
          setupAide options_,
          MPI_Comm comm_,
          bool _isCoarse = false
          );

  void Ax(dfloat* /*x*/, dfloat* /*Ax*/) {}
  void Ax(occa::memory o_x, occa::memory o_Ax) final;

  void residual(dfloat* /*rhs*/, dfloat* /*x*/, dfloat* /*res*/) {}
  void residual(occa::memory o_rhs, occa::memory o_x, occa::memory o_res) final;

  void coarsen(dfloat* /*x*/, dfloat* /*Cx*/) {}
  void coarsen(occa::memory o_x, occa::memory o_Cx) final;

  void prolongate(dfloat* /*x*/, dfloat* /*Px*/) {}
  void prolongate(occa::memory o_x, occa::memory o_Px) final;

  //smoother ops
  void smooth(dfloat* /*rhs*/, dfloat* /*x*/, bool /*x_is_zero*/) {}
  void smooth(occa::memory o_rhs, occa::memory o_x, bool x_is_zero) final;

  void smoother(occa::memory o_x, occa::memory o_Sx, bool xIsZero);

  void smoothChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero);
  void smoothFourthKindChebyshev (occa::memory &o_r, occa::memory &o_x, bool xIsZero);
  void smoothSchwarz (occa::memory &o_r, occa::memory &o_x, bool xIsZero);
  void smoothJacobi (occa::memory &o_r, occa::memory &o_x, bool xIsZero);

  void smootherJacobi    (occa::memory &o_r, occa::memory &o_Sr);

  void Report() final;

  void setupSmoother(elliptic_t* base);
  dfloat maxEigSmoothAx();

  void buildCoarsenerQuadHex(mesh_t **meshLevels, int Nf, int Nc);
};

#endif
