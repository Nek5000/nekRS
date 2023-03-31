#include "nrs.hpp"
#include "filter.hpp"

namespace {
dfloat filterFactorial(int n)
{
  if(n == 0)
    return 1;
  else
    return n * filterFactorial(n - 1);
}
// low Pass
void filterFunctionRelaxation1D(int Nmodes, int Nc, dfloat* A)
{
  // Set all diagonal to 1
  for(int n = 0; n < Nmodes; n++)
    A[n * Nmodes + n] = 1.0;

  int k0 = Nmodes - Nc;
  for (int k = k0; k < Nmodes; k++) {
    dfloat amp = ((k + 1.0 - k0) * (k + 1.0 - k0)) / (Nc * Nc);
    A[k + Nmodes * k] = 1.0 - amp;
  }
}

// jacobi polynomials at [-1,1] for GLL
dfloat filterJacobiP(dfloat a, dfloat alpha, dfloat beta, int N)
{
  dfloat ax = a;

  dfloat* P = (dfloat*) calloc((N + 1), sizeof(dfloat));

  // Zero order
  dfloat gamma0 =
    pow(2,
        (alpha + beta + 1)) / (alpha + beta + 1) * filterFactorial(alpha) * filterFactorial(beta) /
    filterFactorial(alpha + beta);
  dfloat p0     = 1.0 / sqrt(gamma0);

  if (N == 0) {
    free(P);
    return p0;
  }
  P[0] = p0;

  // first order
  dfloat gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0;
  dfloat p1     = ((alpha + beta + 2) * ax / 2 + (alpha - beta) / 2) / sqrt(gamma1);
  if (N == 1) {
    free(P);
    return p1;
  }

  P[1] = p1;

  /// Repeat value in recurrence.
  dfloat aold = 2 / (2 + alpha + beta) * sqrt((alpha + 1.) * (beta + 1.) / (alpha + beta + 3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for(int i = 1; i <= N - 1; ++i) {
    dfloat h1 = 2. * i + alpha + beta;
    dfloat anew = 2. / (h1 + 2.) *
                  sqrt(
      (i + 1.) * (i + 1. + alpha + beta) * (i + 1 + alpha) * (i + 1 + beta) / (h1 + 1) / (h1 + 3));
    dfloat bnew = -(alpha * alpha - beta * beta) / h1 / (h1 + 2);
    P[i + 1] = 1. / anew * ( -aold * P[i - 1] + (ax - bnew) * P[i]);
    aold = anew;
  }

  dfloat pN = P[N];
  free(P);
  return pN;
}

void filterVandermonde1D(int N, int Np, dfloat* r, dfloat* V)
{
  int sk = 0;
  for(int i = 0; i <= N; i++) {
    for(int n = 0; n < Np; n++)
      V[n * Np + sk] = filterJacobiP(r[n],0,0,i);
    sk++;
  }
}


}

double* filterSetup(mesh_t* mesh, const dlong filterNc)
{
  nrsCheck(filterNc < 1, platform->comm.mpiComm, EXIT_FAILURE,
           "filterNc must be at least 1, but is set to %d\n", filterNc);

  // Construct Filter Function
  int Nmodes = mesh->N + 1; // N+1, 1D GLL points

  // Vandermonde matrix
  dfloat* V = (dfloat*) calloc(Nmodes * Nmodes, sizeof(dfloat));
  // Filter matrix, diagonal
  dfloat* A = (dfloat*) calloc(Nmodes * Nmodes, sizeof(dfloat));

  // Construct Filter Function
  filterFunctionRelaxation1D(Nmodes, filterNc, A);

  // Construct Vandermonde Matrix
  filterVandermonde1D(mesh->N, Nmodes, mesh->r, V);

  // Invert the Vandermonde
  int INFO;
  int N        = Nmodes;
  int LWORK    = N * N;
  double* WORK = (double*) calloc(LWORK, sizeof(double));
  int* IPIV  = (int*)   calloc(Nmodes + 1,sizeof(int));
  double* iV   = (double*) calloc(Nmodes * Nmodes, sizeof(double));

  for(int n = 0; n < (Nmodes + 1); n++)
    IPIV[n] = 1;
  for(int n = 0; n < Nmodes * Nmodes; ++n)
    iV[n] = V[n];

  dgetrf_(&N,&N,(double*)iV,&N,IPIV,&INFO);
  nrsCheck(INFO, MPI_COMM_SELF, EXIT_FAILURE, "dgetrf failed", "");

  dgetri_(&N,(double*)iV,&N,IPIV,(double*)WORK,&LWORK,&INFO);
  nrsCheck(INFO, MPI_COMM_SELF, EXIT_FAILURE, " dgetrifailed", "");

  // V*A*V^-1 in row major
  char TRANSA = 'T';
  char TRANSB = 'T';
  double ALPHA = 1.0, BETA = 0.0;
  int MD = Nmodes;
  int ND = Nmodes;
  int KD = Nmodes;
  int LDA = Nmodes;
  int LDB = Nmodes;

  double* C = (double*) calloc(Nmodes * Nmodes, sizeof(double));
  int LDC   = Nmodes;
  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, A, &LDA, iV, &LDB, &BETA, C, &LDC);

  TRANSA = 'T';
  TRANSB = 'N';
  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, V, &LDA, C, &LDB, &BETA, A, &LDC);

  //occa::memory o_filterMT =  platform->device.malloc(Nmodes * Nmodes * sizeof(dfloat), A); // copy Tranpose

  free(C);
  free(V);
  free(iV);
  free(IPIV);
  free(WORK);

  return A;
}
