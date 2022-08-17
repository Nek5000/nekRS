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

#include "elliptic.h"
#include "platform.hpp"
#include "linAlg.hpp"

size_t MGLevel::smootherResidualBytes;
pfloat* MGLevel::smootherResidual;
occa::memory MGLevel::o_smootherResidual;
occa::memory MGLevel::o_smootherResidual2;
occa::memory MGLevel::o_smootherUpdate;

//build a single level
MGLevel::MGLevel(elliptic_t* ellipticBase, int Nc,
                 setupAide options_, MPI_Comm comm_, bool _isCoarse) :
  multigridLevel(ellipticBase->mesh->Nelements * ellipticBase->mesh->Np,
                 (ellipticBase->mesh->Nelements + ellipticBase->mesh->totalHaloPairs) * ellipticBase->mesh->Np,
                 comm_)
{
  isCoarse = _isCoarse;
  
  elliptic = ellipticBase;
  mesh = elliptic->mesh;
  options = options_;
  degree = Nc;

  if(!isCoarse || options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
    this->setupSmoother(ellipticBase);

}

//build a level and connect it to the previous one
MGLevel::MGLevel(elliptic_t* ellipticBase, //finest level
                 mesh_t** meshLevels,
                 elliptic_t* ellipticFine, //previous level
                 elliptic_t* ellipticCoarse, //current level
                 int Nf, int Nc,
                 setupAide options_,
                 MPI_Comm comm_,
                 bool _isCoarse
                 )
  :
  multigridLevel(ellipticCoarse->mesh->Nelements * ellipticCoarse->mesh->Np,
                 (ellipticCoarse->mesh->Nelements + ellipticCoarse->mesh->totalHaloPairs) * ellipticCoarse->mesh->Np,
                 comm_)
{
  
  isCoarse = _isCoarse;
  elliptic = ellipticCoarse;
  mesh = elliptic->mesh;
  options = options_;
  degree = Nc;

  NpF = ellipticFine->mesh->Np;
  o_invDegreeFine = ellipticFine->o_invDegree;

  /* build coarsening and prologation operators to connect levels */
  this->buildCoarsenerQuadHex(meshLevels, Nf, Nc);

  if(!isCoarse || options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE"))
    this->setupSmoother(ellipticBase);
}

void MGLevel::setupSmoother(elliptic_t* ellipticBase)
{

  dfloat minMultiplier;
  options.getArgs("MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR", minMultiplier);

  dfloat maxMultiplier;
  options.getArgs("MULTIGRID CHEBYSHEV MAX EIGENVALUE BOUND FACTOR", maxMultiplier);

  if (options.compareArgs("MULTIGRID SMOOTHER","ASM") ||
      options.compareArgs("MULTIGRID SMOOTHER","RAS")) {
    stype = SmootherType::SCHWARZ;
    smtypeUp = SecondarySmootherType::JACOBI;
    smtypeDown = SecondarySmootherType::JACOBI;
    build(ellipticBase);

    if(options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      smtypeUp = SecondarySmootherType::SCHWARZ;
      smtypeDown = SecondarySmootherType::SCHWARZ;
      stype = SmootherType::CHEBYSHEV;

      if(isCoarse) {
        ChebyshevIterations = 8;
        options.getArgs("COARSE MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations);
      } else {
        ChebyshevIterations = 2;
        options.getArgs("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations);
      }

      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();
      lambda1 = maxMultiplier * rho;
      lambda0 = minMultiplier * rho;
    }
    if(options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI") ||
       options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI")) {
      o_invDiagA = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat));
      ellipticUpdateJacobi(elliptic,o_invDiagA);
      if(options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI"))
        smtypeUp = SecondarySmootherType::JACOBI;
      if(options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI"))
        smtypeDown = SecondarySmootherType::JACOBI;
    }
  } else if (options.compareArgs("MULTIGRID SMOOTHER","DAMPEDJACOBI")) { //default to damped jacobi
    stype = SmootherType::JACOBI;
    smtypeUp = SecondarySmootherType::JACOBI;
    smtypeDown = SecondarySmootherType::JACOBI;

    o_invDiagA = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat));
    ellipticUpdateJacobi(elliptic,o_invDiagA);

    if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      stype = SmootherType::CHEBYSHEV;

      if (!options.getArgs("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations))
        ChebyshevIterations = 2; //default to degree 2

      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();

      lambda1 = maxMultiplier * rho;
      lambda0 = minMultiplier * rho;
    }
  }
}

void MGLevel::Report()
{
  platform_t * platform = platform_t::getInstance();
  hlong hNrows = (hlong) Nrows;

  dlong minNrows = 0, maxNrows = 0;
  hlong totalNrows = 0;
  dfloat avgNrows;

  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, platform->comm.mpiComm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);
  avgNrows = (dfloat) totalNrows / platform->comm.mpiCommSize;

  if (Nrows == 0) Nrows = maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, platform->comm.mpiComm);

  char smootherString[BUFSIZ];
  if (!isCoarse || options.compareArgs("MULTIGRID COARSE SOLVE", "FALSE")) {
    if (stype == SmootherType::CHEBYSHEV && smtypeDown == SecondarySmootherType::JACOBI)
      strcpy(smootherString, "Chebyshev+Jacobi ");
    else if (stype == SmootherType::SCHWARZ)
      strcpy(smootherString, "Schwarz          ");
    else if (stype == SmootherType::JACOBI)
      strcpy(smootherString, "Jacobi           ");
    else if (stype == SmootherType::CHEBYSHEV && smtypeDown == SecondarySmootherType::SCHWARZ)
      strcpy(smootherString, "Chebyshev+Schwarz");
    else
      strcpy(smootherString, "???");
  }

  if (platform->comm.mpiRank == 0) {
    if(isCoarse && options.compareArgs("MULTIGRID COARSE SOLVE","TRUE")) {
      std::string solver;
      if(options.getArgs("COARSE SOLVER", solver)) strcpy(smootherString, solver.c_str());
      printf(     "|    AMG     |   Matrix        | %s\n", smootherString);
      printf("     |            |     Degree %2d   |\n", degree);
    } else {
      printf(     "|    pMG     |   Matrix-free   | %s\n", smootherString);
      printf("     |            |     Degree %2d   |\n", degree);
    }
  }
}

void MGLevel::buildCoarsenerQuadHex(mesh_t** meshLevels, int Nf, int Nc)
{

  const int Nfq = Nf + 1;
  const int Ncq = Nc + 1;
  dfloat *cToFInterp = (dfloat *)calloc(Nfq * Ncq, sizeof(dfloat));
  InterpolationMatrix1D(Nc, Ncq, meshLevels[Nc]->r, Nfq, meshLevels[Nf]->r, cToFInterp);

  pfloat *R = (pfloat *)calloc(Nfq * Ncq, sizeof(pfloat));
  // transpose
  for (int i = 0; i < Ncq; i++) {
    for (int j = 0; j < Nfq; j++) {
      R[i * Nfq + j] = cToFInterp[j * Ncq + i];
    }
  }

  o_R = platform->device.malloc(Nfq * Ncq * sizeof(pfloat), R);

  free(R);
  free(cToFInterp);
}

static void eig(const int Nrows, double* A, double* WR, double* WI)
{
  int NB  = 256;
  char JOBVL  = 'V';
  char JOBVR  = 'V';
  int N = Nrows;
  int LDA = Nrows;
  int LWORK  = (NB + 2) * N;

  double* WORK  = new double[LWORK];
  double* VL  = new double[Nrows * Nrows];
  double* VR  = new double[Nrows * Nrows];

  auto invalid = 0;
  for(int i = 0; i < Nrows * Nrows; i++) {
    if(std::isnan(A[i]) || std::isinf(A[i])) invalid++;
  }
  if(invalid) {
    if(platform->comm.mpiRank == 0) printf("invalid matrix entries!\n");
    ABORT(1);
  }

  int INFO = -999;
  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
          VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);

  if(INFO != 0) {
    if(platform->comm.mpiRank == 0) printf("failed");
    ABORT(INFO);
  }

  delete [] VL;
  delete [] VR;
  delete [] WORK;
}

dfloat MGLevel::maxEigSmoothAx()
{
  MPI_Barrier(platform->comm.mpiComm);
  const double tStart = MPI_Wtime();
  if(platform->comm.mpiRank == 0)  printf("estimating maxEigenvalue ... "); fflush(stdout);
     
  const dlong N = Nrows;
  const dlong M = Ncols;

  hlong Nlocal = (hlong) Nrows;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, platform->comm.mpiComm);

  occa::memory o_invDegree = platform->device.malloc(Nlocal*sizeof(dfloat), elliptic->ogs->invDegree);

  int k;
  if(Ntotal > 10) 
    k = 10;
  else
    k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double* H = (double*) calloc(k * k,sizeof(double));

  // allocate memory for basis
  dfloat* Vx = (dfloat*) calloc(M, sizeof(dfloat));
  occa::memory* o_V = new occa::memory[k + 1];

  size_t offset = 0;
  const size_t vectorSize = ((M * sizeof(dfloat))/ALIGN_SIZE + 1) * ALIGN_SIZE ;

  for(int i = 0; i <= k; i++) {
    if(offset + vectorSize < platform->o_mempool.o_ptr.size()) {
      o_V[i] = platform->o_mempool.o_ptr.slice(offset, vectorSize);
      offset += vectorSize;
    } else {
      o_V[i]  = platform->device.malloc(vectorSize);
    }
  }

  occa::memory o_Vx;
  if(offset + vectorSize < platform->o_mempool.o_ptr.size()) {
    o_Vx = platform->o_mempool.o_ptr.slice(offset, vectorSize);
    offset += vectorSize;
  } else {
    o_Vx  = platform->device.malloc(vectorSize);
  }

  occa::memory o_AVx;
  if(offset + vectorSize < platform->o_mempool.o_ptr.size()) {
    o_AVx = platform->o_mempool.o_ptr.slice(offset, vectorSize);
    offset += vectorSize;
  } else {
    o_AVx  = platform->device.malloc(vectorSize);
  }

  occa::memory o_AVxPfloat = platform->device.malloc(M, sizeof(pfloat));
  occa::memory o_VxPfloat = platform->device.malloc(M, sizeof(pfloat));

  // generate a random vector for initial basis vector
  for (dlong i = 0; i < N; i++) Vx[i] = (dfloat) drand48();

  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ogsGatherScatter(Vx, ogsDfloat, ogsAdd, mesh->ogs);

    if(elliptic->Nmasked > 0){
      dlong* maskIds = (dlong*) calloc(elliptic->Nmasked, sizeof(dlong));
      elliptic->o_maskIds.copyTo(maskIds, elliptic->Nmasked * sizeof(dlong));
      for (dlong i = 0; i < elliptic->Nmasked; i++) Vx[maskIds[i]] = 0.;
      free(maskIds);
    }
  }

  o_Vx.copyFrom(Vx, M*sizeof(dfloat));
  dfloat norm_vo = platform->linAlg->weightedInnerProdMany(
    Nlocal,
    elliptic->Nfields,
    elliptic->Ntotal,
    o_invDegree,
    o_Vx,
    o_Vx,
    platform->comm.mpiComm
  );
  norm_vo = sqrt(norm_vo);

  platform->linAlg->axpbyMany(
    Nlocal,
    elliptic->Nfields,
    elliptic->Ntotal,
    1. / norm_vo,
    o_Vx,
    0.0,
    o_V[0]
  );

  for(int j = 0; j < k; j++) {
    // v[j+1] = invD*(A*v[j])
    platform->copyDfloatToPfloatKernel(M, o_V[j], o_VxPfloat);
    ellipticOperator(elliptic, o_VxPfloat, o_AVxPfloat, pfloatString);
    this->smoother(o_AVxPfloat, o_VxPfloat, true);
    platform->copyPfloatToDfloatKernel(M, o_VxPfloat, o_V[j + 1]);

    // modified Gram-Schmidth
    for(int i = 0; i <= j; i++) {
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = platform->linAlg->weightedInnerProdMany(
        Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        o_invDegree,
        o_V[i],
        o_V[j+1],
        platform->comm.mpiComm
      );

      // v[j+1] = v[j+1] - hij*v[i]
      platform->linAlg->axpbyMany(
        Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        -hij,
        o_V[i],
        1.0,
        o_V[j+1]
      );

      H[i + j * k] = (double) hij;
    }

    if(j + 1 < k) {
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj = platform->linAlg->weightedInnerProdMany(
        Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        o_invDegree,
        o_V[j+1],
        o_V[j+1],
        platform->comm.mpiComm
      );
      norm_vj = sqrt(norm_vj);
      platform->linAlg->scaleMany(
        Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        1 / norm_vj,
        o_V[j+1]
      );

      H[j + 1 + j * k] = (double) norm_vj;
    }
  }

  double* WR = (double*) calloc(k,sizeof(double));
  double* WI = (double*) calloc(k,sizeof(double));

  eig(k, H, WR, WI);

  double rho = 0.;

  for(int i = 0; i < k; i++) {
    double rho_i  = sqrt(WR[i] * WR[i] + WI[i] * WI[i]);

    if(rho < rho_i)
      rho = rho_i;
  }

  free(H);
  free(WR);
  free(WI);

  free(Vx);
  o_Vx.free();
  o_AVx.free();
  o_AVxPfloat.free();
  o_VxPfloat.free();
  o_invDegree.free();

  for(int i = 0; i <= k; i++) o_V[i].free();
  delete[] o_V;

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("%g done (%gs)\n", rho, MPI_Wtime() - tStart); fflush(stdout);

  return rho;
}
