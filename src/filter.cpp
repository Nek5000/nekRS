/*
  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the SoftwarfilterV is
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

#include "ins.h"

// AK. Keep them here to make filterSetup self contained
void filterFunctionRelaxation1D(int Nmodes, int Nc, dfloat *A); 
void filterFunctionRelaxationTri2D(int N, int Nc, dfloat *A); 
void filterFunctionRelaxationTet3D(int N, int Nc, dfloat *A); 

void filterVandermonde1D(int N, int Np, dfloat *r, dfloat *V);
void filterVandermondeTri2D(int N, int Np, dfloat *r, dfloat *s, dfloat *V);
void filterVandermondeTet3D(int N, int Np, dfloat *r, dfloat *s, dfloat *t, dfloat *V);

dfloat filterSimplex2D(dfloat a, dfloat b, int i, int j);
dfloat filterSimplex3D(dfloat a, dfloat b, dfloat c, int i, int j, int k);
dfloat filterJacobiP(dfloat a, dfloat alpha, dfloat beta, int N);

dfloat filterFactorial(int n); 

void filterSetup(ins_t* ins){

  // Setting Filter Here.....
  mesh_t *mesh = ins->mesh; 

  // First construct filter function
  ins->filterS = 10.0; // filter Weight...  
  dfloat filterC = 0.90; // Nc/N i.e. percentage of modes that is not touched 
  ins->options.getArgs("FILTER STRENGTH", ins->filterS); 
  ins->options.getArgs("FILTER CUTOFF RATIO", filterC); 

  // Construc Filter Function
  int Nmodes = 1; 
  // Constructing Dense Operators
  if(ins->elementType==TRIANGLES || ins->elementType==TETRAHEDRA){
    Nmodes = mesh->Np; 
    ins->filterNc = (mesh->N - mymax( (int)round(mesh->N*(1.0-filterC)) -1, 0)) -1;
  }
  else if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA){
    Nmodes = mesh->N+1; // N+1, 1D GLL points 
    // Approximate cutoff order from percentage 
    ins->filterNc = (Nmodes - mymax( (int)round(Nmodes*(1.0-filterC)) -1, 0)) -1;
  }

  
  // Vandermonde matrix
  dfloat *V = (dfloat *) calloc(Nmodes*Nmodes, sizeof(dfloat));
  // Filter matrix, diagonal
  dfloat *A = (dfloat *) calloc(Nmodes*Nmodes, sizeof(dfloat));

  // Construct Filter Function
  if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA)
  filterFunctionRelaxation1D(Nmodes, ins->filterNc, A); 
  else if(ins->elementType==TRIANGLES)
  filterFunctionRelaxationTri2D(mesh->N, ins->filterNc, A); 
  else if(ins->elementType==TETRAHEDRA)
  filterFunctionRelaxationTet3D(mesh->N, ins->filterNc, A); 

  // Construct Vandermonde Matrix
  if(ins->elementType==TRIANGLES)
    filterVandermondeTri2D(mesh->N, Nmodes, mesh->r, mesh->s, V); 
  else if(ins->elementType==TETRAHEDRA)
    filterVandermondeTet3D(mesh->N, Nmodes, mesh->r, mesh->s, mesh->t, V); 
  else if(ins->elementType==QUADRILATERALS || ins->elementType==HEXAHEDRA)
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
    exit(EXIT_FAILURE);
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
  printf("Filter is Activated: filter Strength: %.4f and cutoff Order=%d \n", ins->filterS, ins->filterNc);

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

  // use quadratic damping as Nek, exponential could be better...
  for (int k=Nc; k<Nmodes; k++){
    dfloat amp = ((k-Nc+1.0)*(k-Nc+1.0))/((Nmodes-Nc)*(Nmodes-Nc));
    A[k + Nmodes*k] = 1.0 - amp; 
  }
}

// low Pass
void filterFunctionRelaxationTri2D(int N, int Nc, dfloat *A){
  
  int Np = (N+1)*(N+2)/2; 
  // Set all diagonal to 1 
  for(int n=0; n<Np; n++)
    A[n*Np + n] = 1.0; 
  
#if 0
  // Quadratic, could be too agressive, need to be tested
  int sk = 0; 
  for(int i=0; i<=N; i++){
    for(int j=0; j<=(N-i); j++){
      if((i+j)>=Nc){
        int k = i+j; 
        A[sk*Np + sk] = 1.0- ((k-Nc+0.0)*(k-Nc+0.0))/( (N-Nc)*(N-Nc) ); 
      }
      sk++; 
      printf("sk:%d\n", sk);
    } 
  }
#else
 // exponential
  dfloat eps = 1.e-16; // hard coded now, AK.
  // dfloat eps = 2.220446049250313e-16;
  dfloat sp = 32.0; // hard coded now, AK.
  dfloat alpha = -log(eps); 

  int sk = 0; 
  for(int i=0; i<=N; i++){
    for(int j=0; j<=(N-i); j++){
      if((i+j)>=Nc){
        int k = i+j; 
        A[sk*Np + sk] = exp(-alpha*pow( (double)(k - Nc)/(N-Nc),sp));
      }
      sk++; 
    } 
  }
#endif
}

// low Pass
void filterFunctionRelaxationTet3D(int N, int Nc, dfloat *A){
  
  int Np = (N+1)*(N+2)*(N+3)/6; 
  // Set all diagonal to 1 
  for(int n=0; n<Np; n++)
    A[n*Np + n] = 1.0; 
  
#if 0
  // Quadratic, could be too agressive, need to be tested
  int sk = 0; 
  for(int i=0; i<=N; i++){
    for(int j=0; j<=(N-i); j++){
      for(int k=0; k<=(N-i-j); k++){
      if((i+j+k)>=Nc){
        int k = i+j+k; 
        A[sk*Np + sk] = 1.0- ((k-Nc+0.0)*(k-Nc+0.0))/( (N-Nc)*(N-Nc) ); 
      }
      sk++; 
      } 
    }
  }
#else
  // exponential
  dfloat eps = 1.e-16; // hard coded now, AK.
  // dfloat eps = 2.220446049250313e-16;
  dfloat sp = 16.0; // hard coded now, AK.
  dfloat alpha = -log(eps); 

  int sk = 0; 
   for(int i=0; i<=N; i++){
    for(int j=0; j<=N-i; j++){
      for(int k=0; k<=(N-i-j); k++){
      if((i+j+k)>=Nc){
        A[sk*Np + sk] = exp(-alpha*pow( (double)((i+j+k) - Nc)/(N-Nc),sp));
      }
      sk++; 
    } 
  }
  }
#endif
}

void filterVandermondeTet3D(int N, int Np, dfloat *r, dfloat *s, dfloat *t, dfloat *V){

  // First convert to ab coordinates
  dfloat *a = (dfloat *) calloc(Np, sizeof(dfloat));
  dfloat *b = (dfloat *) calloc(Np, sizeof(dfloat));
  dfloat *c = (dfloat *) calloc(Np, sizeof(dfloat));

  dfloat TOL = 1e-8; 
  // Duffy Transform
  for(int n=0; n<Np; n++){
    if((fabs(s[n]+t[n])>TOL)){
      a[n] = 2*(1+r[n])/(-s[n]-t[n])-1.0;
    }else{
      a[n] = -1.0; 
    }
    if((fabs(t[n]-1.0)>TOL)){
      b[n] = 2.0*(1.0+s[n])/(1.0-t[n])-1.0;
    }else{
      b[n] = -1.0; 
    }
    //
    c[n] = t[n]; 
  }
  
  int sk=0;
  for(int i=0; i<=N; i++){
    for(int j=0; j<=N-i; j++){
      for(int k=0; k<=(N-i-j); k++){
      	for(int n=0; n<Np; n++){
      	  V[n*Np + sk] = filterSimplex3D(a[n], b[n], c[n], i, j, k);
      	}
	     sk++;
      }
    }
  }

  free(a);
  free(b); 
  free(c); 
}



void filterVandermondeTri2D(int N, int Np, dfloat *r, dfloat *s, dfloat *V){

  // First convert to ab coordinates
  dfloat *a = (dfloat *) calloc(Np, sizeof(dfloat));
  dfloat *b = (dfloat *) calloc(Np, sizeof(dfloat));
  for(int n=0; n<Np; n++){

    if(fabs(s[n]-1.0)>1e-8)
      a[n] = 2.0*(1.+r[n])/(1.0-s[n])-1.0;
    else
      a[n] = -1.0; 

    b[n] = s[n];

  }
  
  int sk=0;

  for(int i=0; i<=N; i++){
    for(int j=0; j<=N-i; j++){
      for(int n=0; n<Np; n++){
        V[n*Np + sk] = filterSimplex2D(a[n], b[n], i, j);
      }
      sk++;
    }
  }

  free(a);
  free(b); 
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


dfloat filterSimplex2D(dfloat a, dfloat b, int i, int j){
  // 
  dfloat p1 = filterJacobiP(a,0,0,i);
  dfloat p2 = filterJacobiP(b,2*i+1,0,j);
  dfloat P = sqrt(2.0)*p1*p2*pow(1-b,i);
  return P; 
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
