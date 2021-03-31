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
MGLevel::MGLevel(elliptic_t* ellipticBase, dfloat lambda_, int Nc,
                 setupAide options_, parAlmond::KrylovType ktype_, MPI_Comm comm_) :
  multigridLevel(ellipticBase->mesh->Nelements * ellipticBase->mesh->Np,
                 (ellipticBase->mesh->Nelements + ellipticBase->mesh->totalHaloPairs) * ellipticBase->mesh->Np,
                 ktype_,
                 comm_)
{
  
  elliptic = ellipticBase;
  mesh = elliptic->mesh;
  options = options_;
  lambda = lambda_;
  degree = Nc;
  weighted = false;

  //use weighted inner products
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    weighted = true;
    o_weight = elliptic->o_invDegree;
    weight   = elliptic->invDegree;
  }

  this->setupSmoother(ellipticBase);

  o_xPfloat = platform->device.malloc(Nrows ,  sizeof(pfloat));
  o_rhsPfloat = platform->device.malloc(Nrows ,  sizeof(pfloat));
}

//build a level and connect it to the previous one
MGLevel::MGLevel(elliptic_t* ellipticBase, //finest level
                 mesh_t** meshLevels,
                 elliptic_t* ellipticFine, //previous level
                 elliptic_t* ellipticCoarse, //current level
                 dfloat lambda_,
                 int Nf, int Nc,
                 setupAide options_,
                 parAlmond::KrylovType ktype_,
                 MPI_Comm comm_)
  :
  multigridLevel(ellipticCoarse->mesh->Nelements * ellipticCoarse->mesh->Np,
                 (ellipticCoarse->mesh->Nelements + ellipticCoarse->mesh->totalHaloPairs) * ellipticCoarse->mesh->Np,
                 ktype_,
                 comm_)
{
  
  elliptic = ellipticCoarse;
  mesh = elliptic->mesh;
  options = options_;
  lambda = lambda_;
  degree = Nc;
  weighted = false;

  //use weighted inner products
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    weighted = true;
    o_weight = elliptic->o_invDegree;
    weight   = elliptic->invDegree;

    NpF = ellipticFine->mesh->Np;
    o_invDegree = ellipticFine->ogs->o_invDegree;
  }

  this->setupSmoother(ellipticBase);

  /* build coarsening and prologation operators to connect levels */
  this->buildCoarsenerQuadHex(meshLevels, Nf, Nc);

  o_xPfloat = platform->device.malloc(Nrows ,  sizeof(pfloat));
  o_rhsPfloat = platform->device.malloc(Nrows ,  sizeof(pfloat));
}

void MGLevel::setupSmoother(elliptic_t* ellipticBase)
{
  
  if (degree == 1) return; // solved by coarse grid solver

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
      if (!options.getArgs("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations))
        ChebyshevIterations = 2;   //default to degree 2
      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();
      lambda1 = 1.1 * rho;
      lambda0 = rho / 10.;
    }
    if(options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI") ||
       options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI")) {
      dfloat* invDiagA;
      std::vector<pfloat> casted_invDiagA(mesh->Np * mesh->Nelements, 0.0);
      ellipticBuildJacobi(elliptic,&invDiagA);
      for(dlong i = 0; i < mesh->Np * mesh->Nelements; ++i)
        casted_invDiagA[i] = static_cast<pfloat>(invDiagA[i]);
      o_invDiagA = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat), casted_invDiagA.data());
      if(options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI"))
        smtypeUp = SecondarySmootherType::JACOBI;
      if(options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI"))
        smtypeDown = SecondarySmootherType::JACOBI;
    }
  } else if (options.compareArgs("MULTIGRID SMOOTHER","DAMPEDJACOBI")) { //default to damped jacobi
    smtypeUp = SecondarySmootherType::JACOBI;
    smtypeDown = SecondarySmootherType::JACOBI;
    dfloat* invDiagA;
    ellipticBuildJacobi(elliptic,&invDiagA);
    std::vector<pfloat> casted_invDiagA(mesh->Np * mesh->Nelements, 0.0);
    for(dlong i = 0; i < mesh->Np * mesh->Nelements; ++i)
      casted_invDiagA[i] = static_cast<pfloat>(invDiagA[i]);

    o_invDiagA = platform->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat), casted_invDiagA.data());

    if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      stype = SmootherType::CHEBYSHEV;

      if (!options.getArgs("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations))
        ChebyshevIterations = 2; //default to degree 2

      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();

      lambda1 = 1.1 * rho;
      lambda0 = rho / 10.;
    }else {

      std::string invalidSmootherName;
      options.getArgs("MULTIGRID SMOOTHER", invalidSmootherName);
      if(platform->comm.mpiRank == 0) printf("Smoother %s is not supported!\n", invalidSmootherName.c_str());
      ABORT(EXIT_FAILURE);
    }
    free(invDiagA);
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
  if (degree != 1) {
    if (stype == SmootherType::CHEBYSHEV && smtypeDown == SecondarySmootherType::JACOBI)
      strcpy(smootherString, "Chebyshev+Jacobi ");
    else if (stype == SmootherType::SCHWARZ)
      strcpy(smootherString, "Schwarz          ");
    else if (stype == SmootherType::CHEBYSHEV && smtypeDown == SecondarySmootherType::SCHWARZ)
      strcpy(smootherString, "Chebyshev+Schwarz");
  }

  if (platform->comm.mpiRank == 0) {
    if(degree == 1) {
      strcpy(smootherString, "BoomerAMG        ");
      printf(     "|    AMG     |   Matrix        | %s |\n", smootherString);
      printf("     |            |     Degree %2d   |                   |\n", degree);
    } else {
      printf(     "|    pMG     |   Matrix-free   | %s |\n", smootherString);
      printf("     |            |     Degree %2d   |                   |\n", degree);
    }
  }
}

void MGLevel::buildCoarsenerQuadHex(mesh_t** meshLevels, int Nf, int Nc)
{
  
  int NqFine   = Nf + 1;
  int NqCoarse = Nc + 1;
  dfloat* P    = (dfloat*) calloc(NqFine * NqCoarse,sizeof(dfloat));
  dfloat* Ptmp = (dfloat*) calloc(NqFine * NqCoarse,sizeof(dfloat));

  //initialize P as identity
  for (int i = 0; i < NqCoarse; i++) P[i * NqCoarse + i] = 1.0;

  for (int n = Nc; n < Nf; n++) {
    int Nqp1 = n + 2;
    int Nq   = n + 1;

    //copy P
    for (int i = 0; i < Nq * NqCoarse; i++) Ptmp[i] = P[i];

    //Multiply by the raise op
    for (int i = 0; i < Nqp1; i++)
      for (int j = 0; j < NqCoarse; j++) {
        P[i * NqCoarse + j] = 0.;
        for (int k = 0; k < Nq; k++)
          P[i * NqCoarse + j] += meshLevels[n]->interpRaise[i * Nq + k] * Ptmp[k * NqCoarse + j];
      }
  }

  //the coarsen matrix is P^T
  R = (dfloat*) calloc(NqFine * NqCoarse,sizeof(dfloat));
  for (int i = 0; i < NqCoarse; i++)
    for (int j = 0; j < NqFine; j++)
      R[i * NqFine + j] = P[j * NqCoarse + i];
  o_R = platform->device.malloc(NqFine * NqCoarse * sizeof(dfloat), R);

  free(P);
  free(Ptmp);
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

  int INFO = -999;

  dgeev_ (&JOBVL, &JOBVR, &N, A, &LDA, WR, WI,
          VL, &LDA, VR, &LDA, WORK, &LWORK, &INFO);

  assert(INFO == 0);

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
  //  occa::memory *o_V = (occa::memory *) calloc(k+1, sizeof(occa::memory));
  occa::memory* o_V = new occa::memory[k + 1];

  occa::memory o_Vx  = platform->device.malloc(M * sizeof(dfloat),Vx);
  occa::memory o_AVx = platform->device.malloc(M * sizeof(dfloat),Vx);
  occa::memory o_AVxPfloat = platform->device.malloc(M ,  sizeof(pfloat));
  occa::memory o_VxPfloat = platform->device.malloc(M ,  sizeof(pfloat));

  for(int i = 0; i <= k; i++)
    o_V[i] = platform->device.malloc(M * sizeof(dfloat),Vx);

  // generate a random vector for initial basis vector
  for (dlong i = 0; i < N; i++) Vx[i] = (dfloat) drand48();

  //gather-scatter
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ogsGatherScatter(Vx, ogsDfloat, ogsAdd, mesh->ogs);

    for (dlong i = 0; i < elliptic->Nmasked; i++) Vx[elliptic->maskIds[i]] = 0.;
  }

  o_Vx.copyFrom(Vx); //copy to device
  dfloat norm_vo = platform->linAlg->weightedInnerProdMany(
    Nlocal,
    elliptic->Nfields,
    elliptic->Ntotal,
    elliptic->o_invDegree,
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
    //this->Ax(o_V[j],o_AVx);
    ellipticOperator(elliptic,o_V[j],o_AVx,dfloatString);
    elliptic->copyDfloatToPfloatKernel(M, o_AVxPfloat, o_AVx);
    this->smoother(o_AVxPfloat, o_VxPfloat, true);
    elliptic->copyPfloatToDPfloatKernel(M, o_VxPfloat, o_V[j + 1]);

    // modified Gram-Schmidth
    for(int i = 0; i <= j; i++) {
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij = platform->linAlg->weightedInnerProdMany(
        Nlocal,
        elliptic->Nfields,
        elliptic->Ntotal,
        elliptic->o_invDegree,
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
        elliptic->o_invDegree,
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

  // free memory
  free(H);
  free(WR);
  free(WI);

  free(Vx);
  o_Vx.free();
  o_AVx.free();
  o_AVxPfloat.free();
  o_VxPfloat.free();
  for(int i = 0; i <= k; i++) o_V[i].free();
  //free((void*)o_V);
  delete[] o_V;

  MPI_Barrier(platform->comm.mpiComm);
  if(platform->comm.mpiRank == 0)  printf("%g done (%gs)\n", rho, MPI_Wtime() - tStart); fflush(stdout);

  return rho;
}
