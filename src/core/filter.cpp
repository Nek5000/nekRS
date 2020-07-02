#include "nrs.hpp"
#include "ins.h"

void filterFunctionRelaxation1D(int Nmodes, int Nc, dfloat *A); 

void filterVandermonde1D(int N, int Np, dfloat *r, dfloat *V);

dfloat filterSimplex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k);
dfloat filterJacobiP(dfloat a, dfloat alpha, dfloat beta, int N);

dfloat filterFactorial(int n); 

void filterSetup(ins_t* ins){

  mesh_t *mesh = ins->mesh; 

  // First construct filter function
  ins->filterS = 10.0; // filter Weight...  
  ins->options.getArgs("HPFRT STRENGTH", ins->filterS); 
  ins->options.getArgs("HPFRT MODES", ins->filterNc); 
  ins->filterS = -1.0*fabs(ins->filterS);

  // Construct Filter Function
  int Nmodes = mesh->N+1; // N+1, 1D GLL points 
  
  // Vandermonde matrix
  dfloat *V = (dfloat *) calloc(Nmodes*Nmodes, sizeof(dfloat));
  // Filter matrix, diagonal
  dfloat *A = (dfloat *) calloc(Nmodes*Nmodes, sizeof(dfloat));

  // Construct Filter Function
  filterFunctionRelaxation1D(Nmodes, ins->filterNc, A); 

  // Construct Vandermonde Matrix
  filterVandermonde1D(mesh->N, Nmodes, mesh->r, V); 


  // Invert the Vandermonde 
  int INFO;
  int N        = Nmodes; 
  int LWORK    = N*N;
  double *WORK = (double *) calloc(LWORK, sizeof(double));
  int   *IPIV  = (int *)   calloc(Nmodes+1,sizeof(int)); 
  double *iV   = (double*) calloc(Nmodes*Nmodes, sizeof(double));

  for(int n=0; n<(Nmodes+1);n++)
    IPIV[n] = 1; 
  for(int n=0;n<Nmodes*Nmodes;++n)
    iV[n] = V[n];

  dgetrf_(&N,&N,(double*)iV,&N,IPIV,&INFO);
  dgetri_(&N,(double*)iV,&N,IPIV,(double*)WORK,&LWORK,&INFO);

  if(INFO){
    printf("DGE_TRI/TRF error: %d \n", INFO);
    ABORT(EXIT_FAILURE);
  }
 

  // Some ugly operations to perform V*A*V^-1 in row major, watchout the transposes AK. !!!
  char TRANSA = 'T'; 
  char TRANSB = 'T'; 
  double ALPHA= 1.0, BETA= 0.0;
  int MD = Nmodes;
  int ND = Nmodes;
  int KD = Nmodes;
  int LDA = Nmodes;
  int LDB = Nmodes;

  double *C = (double *) calloc(Nmodes*Nmodes, sizeof(double));
  int LDC   = Nmodes;
  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, A, &LDA, iV, &LDB, &BETA, C, &LDC);
  
  TRANSA = 'T'; 
  TRANSB = 'N'; 
  dgemm_(&TRANSA, &TRANSB, &MD, &ND, &KD, &ALPHA, V, &LDA, C, &LDB, &BETA, A, &LDC);
  
  // store filter matrix, row major again !!!
  ins->filterM = (dfloat *) calloc(Nmodes*Nmodes, sizeof(dfloat)); 
  for(int c=0; c<Nmodes; c++){
    for(int r=0; r<Nmodes; r++){
      ins->filterM[c+ r*Nmodes] = A[r + c*Nmodes];
    }
  }

  // Copy To Device
  ins->o_filterMT =  mesh->device.malloc(Nmodes*Nmodes*sizeof(dfloat), A); // copy Tranpose

  if(mesh->rank==0)
    printf("High pass filter relaxation: chi = %.4f using %d mode(s)\n", 
           fabs(ins->filterS), ins->filterNc);

  free(A); 
  free(C); 
  free(V); 
  free(iV); 
  free(IPIV);
  free(WORK);
}


// low Pass
void filterFunctionRelaxation1D(int Nmodes, int Nc, dfloat *A){
  
  // Set all diagonal to 1 
  for(int n=0; n<Nmodes; n++)
    A[n*Nmodes + n] = 1.0; 

  int k0 = Nmodes - Nc;
  for (int k=k0; k<Nmodes; k++){
    dfloat amp = ((k+1.0 - k0)*(k+1.0 - k0))/(Nc*Nc);
    A[k + Nmodes*k] = 1.0 - amp; 
  }
}

void filterVandermonde1D(int N, int Np, dfloat *r, dfloat *V){
  int sk=0;
  for(int i=0; i<=N; i++){
    for(int n=0; n<Np; n++){
      V[n*Np + sk] = filterJacobiP(r[n],0,0,i);
    }
    sk++;
  }
}

dfloat filterSimplex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k){
  // 
  dfloat p1 = filterJacobiP(a,0,0,i); 
  dfloat p2 = filterJacobiP(b,2*i+1,0,j); 
  dfloat p3 = filterJacobiP(c,2*(i+j)+2,0,k);
  dfloat P  = 2.0*sqrt(2.0)*p1*p2*pow(1-b,i)*p3*pow(1-c,i+j);
  return P; 
}


// jacobi polynomials at [-1,1] for GLL
dfloat filterJacobiP(dfloat a, dfloat alpha, dfloat beta, int N){

  dfloat ax = a; 

  dfloat *P = (dfloat *) calloc((N+1), sizeof(dfloat));

  // Zero order
  dfloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*filterFactorial(alpha)*filterFactorial(beta)/filterFactorial(alpha+beta);
  dfloat p0     = 1.0/sqrt(gamma0);

  if (N==0){ free(P); return p0;}
  P[0] = p0; 

  // first order
  dfloat gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
  dfloat p1     = ((alpha+beta+2)*ax/2 + (alpha-beta)/2)/sqrt(gamma1);
  if (N==1){free(P); return p1;} 

  P[1] = p1;

  /// Repeat value in recurrence.
  dfloat aold = 2/(2+alpha+beta)*sqrt((alpha+1.)*(beta+1.)/(alpha+beta+3.));
  /// Forward recurrence using the symmetry of the recurrence.
  for(int i=1;i<=N-1;++i){
    dfloat h1 = 2.*i+alpha+beta;
    dfloat anew = 2./(h1+2.)*sqrt( (i+1.)*(i+1.+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
    dfloat bnew = -(alpha*alpha-beta*beta)/h1/(h1+2);
    P[i+1] = 1./anew*( -aold*P[i-1] + (ax-bnew)*P[i]);
    aold =anew;
  }
  
  dfloat pN = P[N]; 
  free(P);
  return pN;

}


dfloat filterFactorial(int n){

  if(n==0)
    return 1;
  else
    return n*filterFactorial(n-1);
}
