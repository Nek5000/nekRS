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

  o_xPfloat = mesh->device.malloc(Nrows*sizeof(pfloat));
  o_rhsPfloat = mesh->device.malloc(Nrows*sizeof(pfloat));

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
  if (elliptic->elementType == TRIANGLES || elliptic->elementType == TETRAHEDRA)
    this->buildCoarsenerTriTet(meshLevels, Nf, Nc);
  else
    this->buildCoarsenerQuadHex(meshLevels, Nf, Nc);

  o_xPfloat = mesh->device.malloc(Nrows*sizeof(pfloat));
  o_rhsPfloat = mesh->device.malloc(Nrows*sizeof(pfloat));
}

void MGLevel::setupSmoother(elliptic_t* ellipticBase)
{
  if (degree == 1) return; // solved by coarse grid solver

  if (options.compareArgs("MULTIGRID SMOOTHER","ASM") ||
             options.compareArgs("MULTIGRID SMOOTHER","RAS")) {
    stype = SCHWARZ;
    smtypeUp = JACOBI;
    smtypeDown = JACOBI;
    build(ellipticBase);

    if(options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      smtypeUp = SCHWARZ_SMOOTH;
      smtypeDown = SCHWARZ_SMOOTH;
      stype = CHEBYSHEV;
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
      std::vector<pfloat> casted_invDiagA(mesh->Np*mesh->Nelements, 0.0);
      ellipticBuildJacobi(elliptic,&invDiagA);
      for(dlong i = 0 ; i < mesh->Np*mesh->Nelements; ++i){
        casted_invDiagA[i] = static_cast<pfloat>(invDiagA[i]);
      }
      o_invDiagA = mesh->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat), casted_invDiagA.data());
      if(options.compareArgs("MULTIGRID UPWARD SMOOTHER","JACOBI"))
        smtypeUp = JACOBI;
      if(options.compareArgs("MULTIGRID DOWNWARD SMOOTHER","JACOBI"))
        smtypeDown = JACOBI;
    }
  } else if (options.compareArgs("MULTIGRID SMOOTHER","DAMPEDJACOBI")) { //default to damped jacobi
    smtypeUp = JACOBI;
    smtypeDown = JACOBI;
    dfloat* invDiagA;
    ellipticBuildJacobi(elliptic,&invDiagA);
    std::vector<pfloat> casted_invDiagA(mesh->Np*mesh->Nelements, 0.0);
    for(dlong i = 0 ; i < mesh->Np*mesh->Nelements; ++i){
      casted_invDiagA[i] = static_cast<pfloat>(invDiagA[i]);
    }

    o_invDiagA = mesh->device.malloc(mesh->Np * mesh->Nelements * sizeof(pfloat), casted_invDiagA.data());

    if (options.compareArgs("MULTIGRID SMOOTHER","CHEBYSHEV")) {
      stype = CHEBYSHEV;

      if (!options.getArgs("MULTIGRID CHEBYSHEV DEGREE", ChebyshevIterations))
        ChebyshevIterations = 2; //default to degree 2

      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();

      lambda1 = 1.1 * rho;
      lambda0 = rho / 10.;
    }else {
      stype = RICHARDSON;

      //estimate the max eigenvalue of S*A
      dfloat rho = this->maxEigSmoothAx();

      //set the stabilty weight (jacobi-type interation)
      lambda0 = (4. / 3.) / rho;

      for (dlong n = 0; n < mesh->Np * mesh->Nelements; n++)
        casted_invDiagA[n] *= static_cast<pfloat>(lambda0);

      //update diagonal with weight
      o_invDiagA.copyFrom(casted_invDiagA.data());
    }
    free(invDiagA);
  }
}

void MGLevel::Report()
{
  hlong hNrows = (hlong) Nrows;

  dlong minNrows = 0, maxNrows = 0;
  hlong totalNrows = 0;
  dfloat avgNrows;

  MPI_Allreduce(&Nrows, &maxNrows, 1, MPI_DLONG, MPI_MAX, mesh->comm);
  MPI_Allreduce(&hNrows, &totalNrows, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  avgNrows = (dfloat) totalNrows / mesh->size;

  if (Nrows == 0) Nrows = maxNrows; //set this so it's ignored for the global min
  MPI_Allreduce(&Nrows, &minNrows, 1, MPI_DLONG, MPI_MIN, mesh->comm);

  char smootherString[BUFSIZ];
  if (degree != 1){
    if (stype == RICHARDSON && smtypeDown == JACOBI)
      strcpy(smootherString, "Damped Jacobi    ");
    else if (stype == CHEBYSHEV && smtypeDown == JACOBI)
      strcpy(smootherString, "Chebyshev+Jacobi ");
    else if (stype == SCHWARZ)
      strcpy(smootherString, "Schwarz          ");
    else if (stype == CHEBYSHEV && smtypeDown == SCHWARZ_SMOOTH)
      strcpy(smootherString, "Chebyshev+Schwarz");
  }

  if (mesh->rank == 0) {
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

void MGLevel::buildCoarsenerTriTet(mesh_t** meshLevels, int Nf, int Nc)
{
  int NpFine   = meshLevels[Nf]->Np;
  int NpCoarse = meshLevels[Nc]->Np;
  dfloat* P    = (dfloat*) calloc(NpFine * NpCoarse,sizeof(dfloat));
  dfloat* Ptmp = (dfloat*) calloc(NpFine * NpCoarse,sizeof(dfloat));

  //initialize P as identity (which it is for SPARSE)
  for (int i = 0; i < NpCoarse; i++) P[i * NpCoarse + i] = 1.0;

  for (int n = Nc; n < Nf; n++) {
    int Npp1 = meshLevels[n + 1]->Np;
    int Np   = meshLevels[n  ]->Np;

    //copy P
    for (int i = 0; i < Np * NpCoarse; i++) Ptmp[i] = P[i];

    //Multiply by the raise op
    for (int i = 0; i < Npp1; i++)
      for (int j = 0; j < NpCoarse; j++) {
        P[i * NpCoarse + j] = 0.;
        for (int k = 0; k < Np; k++)
          P[i * NpCoarse + j] += meshLevels[n]->interpRaise[i * Np + k] * Ptmp[k * NpCoarse + j];
      }
  }

  if (elliptic->options.compareArgs("BASIS","BERN")) {
    dfloat* BBP = (dfloat*) calloc(NpFine * NpCoarse,sizeof(dfloat));
    for (int j = 0; j < NpFine; j++)
      for (int i = 0; i < NpCoarse; i++)
        for (int k = 0; k < NpCoarse; k++)
          for (int l = 0; l < NpFine; l++)
            BBP[i + j * NpCoarse] += meshLevels[Nf]->invVB[l + j * NpFine] * P[k + l * NpCoarse] *
                                     meshLevels[Nc]->VB[i + k * NpCoarse];
    for (int j = 0; j < NpFine; j++)
      for (int i = 0; i < NpCoarse; i++)
        P[i + j * NpCoarse] = BBP[i + j * NpCoarse];
    free(BBP);
  }

  //the coarsen matrix is P^T
  R = (dfloat*) calloc(NpFine * NpCoarse,sizeof(dfloat));
  for (int i = 0; i < NpCoarse; i++)
    for (int j = 0; j < NpFine; j++)
      R[i * NpFine + j] = P[j * NpCoarse + i];
  o_R = elliptic->mesh->device.malloc(NpFine * NpCoarse * sizeof(dfloat), R);

  free(P);
  free(Ptmp);
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
  o_R = elliptic->mesh->device.malloc(NqFine * NqCoarse * sizeof(dfloat), R);

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
  const dlong N = Nrows;
  const dlong M = Ncols;

  int k = 10;

  hlong Nlocal = (hlong) Nrows;
  hlong Ntotal = 0;
  MPI_Allreduce(&Nlocal, &Ntotal, 1, MPI_HLONG, MPI_SUM, mesh->comm);
  if(k > Ntotal) k = (int) Ntotal;

  // do an arnoldi

  // allocate memory for Hessenberg matrix
  double* H = (double*) calloc(k * k,sizeof(double));

  // allocate memory for basis
  dfloat* Vx = (dfloat*) calloc(M, sizeof(dfloat));
  //  occa::memory *o_V = (occa::memory *) calloc(k+1, sizeof(occa::memory));
  occa::memory* o_V = new occa::memory[k + 1];

  occa::memory o_Vx  = mesh->device.malloc(M * sizeof(dfloat),Vx);
  occa::memory o_AVx = mesh->device.malloc(M * sizeof(dfloat),Vx);
  occa::memory o_AVxPfloat = mesh->device.malloc(M * sizeof(pfloat));
  occa::memory o_VxPfloat = mesh->device.malloc(M * sizeof(pfloat));

  for(int i = 0; i <= k; i++)
    o_V[i] = mesh->device.malloc(M * sizeof(dfloat),Vx);

  // generate a random vector for initial basis vector
  for (dlong i = 0; i < N; i++) Vx[i] = (dfloat) drand48();

  //gather-scatter
  if (options.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ogsGatherScatter(Vx, ogsDfloat, ogsAdd, mesh->ogs);

    for (dlong i = 0; i < elliptic->Nmasked; i++) Vx[elliptic->maskIds[i]] = 0.;
  }

  o_Vx.copyFrom(Vx); //copy to device
  dfloat norm_vo = ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_Vx, o_Vx);
  norm_vo = sqrt(norm_vo);

  ellipticScaledAdd(elliptic, 1. / norm_vo, o_Vx, 0., o_V[0]);

  for(int j = 0; j < k; j++) {
    // v[j+1] = invD*(A*v[j])
    //this->Ax(o_V[j],o_AVx);
    ellipticOperator(elliptic,o_V[j],o_AVx,dfloatString);
    elliptic->copyDfloatToPfloatKernel(M, o_AVxPfloat, o_AVx);
    this->smoother(o_AVxPfloat, o_VxPfloat, true);
    elliptic->copyPfloatToDPfloatKernel(M, o_VxPfloat, o_V[j+1]);

    // modified Gram-Schmidth
    for(int i = 0; i <= j; i++) {
      // H(i,j) = v[i]'*A*v[j]
      dfloat hij =
        ellipticWeightedInnerProduct(elliptic, elliptic->o_invDegree, o_V[i], o_V[j + 1]);

      // v[j+1] = v[j+1] - hij*v[i]
      ellipticScaledAdd(elliptic, -hij, o_V[i], 1., o_V[j + 1]);

      H[i + j * k] = (double) hij;
    }

    if(j + 1 < k) {
      // v[j+1] = v[j+1]/||v[j+1]||
      dfloat norm_vj = ellipticWeightedInnerProduct(elliptic,
                                                    elliptic->o_invDegree,
                                                    o_V[j + 1],
                                                    o_V[j + 1]);
      norm_vj = sqrt(norm_vj);
      ellipticScaledAdd(elliptic, 1 / norm_vj, o_V[j + 1], 0., o_V[j + 1]);

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
  for(int i = 0; i <= k; i++) o_V[i].free();
  //free((void*)o_V);
  delete[] o_V;

  // if((mesh->rank==0)&&(options.compareArgs("VERBOSE","TRUE"))) printf("weight = %g \n", rho);

  return rho;
}
