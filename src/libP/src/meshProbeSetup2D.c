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

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mesh2D.h"

void meshProbeSetup2D(mesh2D *mesh, dfloat *pX, dfloat *pY){
 

 dfloat mindist = -1.e-12; // Set minum distance
  //

  mesh->probeN = 0; 


  #if 1
 
  dfloat *probeR           = (dfloat *) calloc(10*mesh->probeNTotal,sizeof(dfloat));
  dfloat *probeS           = (dfloat *) calloc(10*mesh->probeNTotal,sizeof(dfloat));
  mesh->probeElementIds    = (int *)   calloc(10*mesh->probeNTotal,sizeof(int));
  mesh->probeIds           = (int *)   calloc(10*mesh->probeNTotal,sizeof(int));

 
  double *A        = (double *) calloc((mesh->dim+1)*mesh->Nverts,sizeof(double)); 
  int   *IPIV     = (int *)   calloc(mesh->Nverts,sizeof(int)); 
  int   *IPIV2    = (int *)   calloc((mesh->Np+1),sizeof(int)); 
 
  dfloat *b       = (dfloat *) calloc((mesh->dim+1),sizeof(dfloat));
  double *q       = (double *) calloc(mesh->Nverts,sizeof(double));

  for(int n=0; n<mesh->Nverts;n++)
    IPIV[n] = 1; 

  for(int n=0; n<(mesh->Np+1);n++)
    IPIV2[n] = 1; 


  int N    = (mesh->dim+1); // A->Nrwos
  int NRHS =  1; // B->Ncolumns
  int LDA  = N; 
  int LDB  = (mesh->dim + 1); // B->Nrows
  int INFO;

  // for(int n=0; n<mesh->mesh->probeNTotal; n++){
  // // Coordinates of probe
  // printf("Probe %d  pX: %g pY:%g \n",n, pX[n], pY[n]);
  // }

// Assumes each probe is in one element, may change later 
  for(int n=0; n<mesh->probeNTotal; n++){
    // Coordinates of probe
    b[0] = 1.0; 
    b[1] = pX[n]; 
    b[2] = pY[n]; 

    for (int e=0;e<mesh->Nelements;e++) {
    // Create A[1 vx vy]
      for (int v=0;v<mesh->Nverts;v++) {
        dfloat vx = mesh->EX[e*mesh->Nverts+v];
        dfloat vy = mesh->EY[e*mesh->Nverts+v];
        //
        A[v*mesh->Nverts + 0] = 1.0;
        A[v*mesh->Nverts + 1] = vx;
        A[v*mesh->Nverts + 2] = vy;
      } 
     
     for(int l=0; l<(mesh->dim+1); l++)
     q[l] = b[l]; 


    dgesv_(&N,&NRHS,(double*) A,&LDA,IPIV, (double*) q,&LDB,&INFO);

    if(INFO)
      printf("DGSEV error: %d \n", INFO);

    // Check all non-negative barycentric coordinates 
    // Assumes a probe can be represented by single element!!!!

      dfloat qmin = q[0]; 
      for(int i =1; i<(mesh->dim+1); i++)
       qmin = mymin(qmin, q[i]);
        
      
      // Catch the element
      if(qmin>mindist){

        // Increase number of probes        
        mesh->probeIds[mesh->probeN] = n+1;
        // hold element ids
        mesh->probeElementIds[mesh->probeN] = e; 
        // hold local r,s coordinates
        dfloat l1 =  q[2]; 
        dfloat l2 =  q[0]; 
        dfloat l3 =  q[1]; 

        probeR[mesh->probeN] = 2.*l3 -1.; 
        probeS[mesh->probeN] = 2.*l1 -1.;

        printf("element: %d probe Id: %d qmin:%.5e R: %.5e S:%.5e\n", e, n, qmin, probeR[mesh->probeN],probeS[mesh->probeN]); 
        //
        mesh->probeN++;

      }
    }



  }




#else

mesh->probeN = 0; 
 
  dfloat *probeR           = (dfloat *) calloc(mesh->probeNTotal,sizeof(dfloat));
  dfloat *probeS           = (dfloat *) calloc(mesh->probeNTotal,sizeof(dfloat));
  mesh->probeElementIds    = (int *)   calloc(mesh->probeNTotal,sizeof(int));

 
  double *A        = (double *) calloc((mesh->dim+1)*mesh->Nverts,sizeof(double)); 
  int   *IPIV     = (int *)   calloc(mesh->Nverts,sizeof(int)); 
  int   *IPIV2    = (int *)   calloc((mesh->Np+1),sizeof(int)); 
 
  dfloat *b       = (dfloat *) calloc((mesh->dim+1)*mesh->probeNTotal,sizeof(dfloat));
  double *q       = (double *) calloc(mesh->Nverts*mesh->probeNTotal,sizeof(double));

  for(int n=0; n<mesh->Nverts;n++)
    IPIV[n] = 1; 

  for(int n=0; n<(mesh->Np+1);n++)
    IPIV2[n] = 1; 


  int N    = (mesh->dim+1); // A->Nrwos
  int NRHS =  mesh->probeNTotal; // B->Ncolumns
  int LDA  = N; 
  int LDB  = (mesh->dim + 1); // B->Nrows
  int INFO;



  // Assumes each probe is in one element, may change later 
  for(int n=0; n<mesh->probeNTotal; n++){
    // Coordinates of probe
    b[n*mesh->probeNTotal + 0] = 1.0; 
    b[n*mesh->probeNTotal + 1] = pX[n]; 
    b[n*mesh->probeNTotal + 2] = pY[n]; 
  }
  
  
  //
  for (int e=0;e<mesh->Nelements;e++) {
    // Create A[1 vx vy]
    for (int v=0;v<mesh->Nverts;v++) {
      dfloat vx = mesh->EX[e*mesh->Nverts+v];
      dfloat vy = mesh->EY[e*mesh->Nverts+v];
      //
      A[v*mesh->Nverts + 0] = 1.0;
      A[v*mesh->Nverts + 1] = vx;
      A[v*mesh->Nverts + 2] = vy;
    } 


    for(int l=0; l<mesh->probeNTotal*(mesh->dim+1); l++)
      q[l] = b[l]; 
    
  
    // Find barycentric coordinates
    // q = A^-1*b

    dgesv_(&N,&NRHS,(double*) A,&LDA,IPIV, (double*) q,&LDB,&INFO);

    if(INFO)
      printf("DGSEV error: %d \n", INFO);

    // Check all non-negative barycentric coordinates 
    // Assumes a probe can be represented by single element!!!!
    for(int n=0; n<mesh->probeNTotal; n++){

      dfloat qmin = q[n*mesh->probeNTotal + 0]; 
      for(int i =1; i<(mesh->dim+1); i++)
	     qmin = mymin(qmin, q[n*mesh->probeNTotal + i]);
        
      // Catch the element
      if(qmin>mindist){
      	mesh->probeN++;
      	// hold element ids
      	mesh->probeElementIds[n] = e; 
      	// hold local r,s coordinates
      	dfloat l1 =  q[n*mesh->probeNTotal + 2]; 
      	dfloat l2 =  q[n*mesh->probeNTotal + 0]; 
      	dfloat l3 =  q[n*mesh->probeNTotal + 1]; 

      	probeR[n] = 2.*l3 -1.; 
      	probeS[n] = 2.*l1 -1.;

      	printf("element: %d probe %d qmin:%.5e R: %.5e S:%.5e\n", e, n, qmin, probeR[n],probeS[n]); 

      }
    }
  }
 
#endif


 printf("probe Number: %d \n", mesh->probeN); 



  if(mesh->probeN){
    //Reallocate ProbeIds and Element Ids, Now take cares of  cares 
    mesh->probeIds        = (int *)   realloc(mesh->probeIds,        mesh->probeN*sizeof(int));
    mesh->probeElementIds = (int *)   realloc(mesh->probeElementIds, mesh->probeN*sizeof(int));
    probeR                = (dfloat *) realloc(probeR, mesh->probeN*sizeof(dfloat));
    probeS                = (dfloat *) realloc(probeS, mesh->probeN*sizeof(dfloat));

    // Compute Vandermonde Matrix and invert  it
    dfloat *V = (dfloat *) calloc(mesh->Np* (mesh->N+1)*(mesh->N+2)/2, sizeof(dfloat));
    meshVandermonde2D(mesh->N, mesh->Np, mesh->r, mesh->s, V);

    double *dV = (double*) calloc(mesh->Np* (mesh->N+1)*(mesh->N+2)/2, sizeof(double));
    for(int n=0;n<mesh->Np*(mesh->N+1)*(mesh->N+2)/2;++n)
      dV[n] = V[n];
    
    //
    N    = mesh->Np; 
    int LWORK = mesh->Np*mesh->Np;
    double *WORK = (double *) calloc(LWORK, sizeof(double));

    dgetrf_(&N,&N,(double*)dV,&N,IPIV2,&INFO);
    dgetri_(&N,(double*)dV,&N,IPIV2,(double*)WORK,&LWORK,&INFO);
    
    if(INFO)
      printf("DGE_TRI/TRF error: %d \n", INFO);

  
    // Compute Vandermonde matrix of probes
    dfloat *Vprobe = (dfloat *) calloc(mesh->probeN*mesh->Np,sizeof(dfloat));
    meshVandermonde2D(mesh->N, mesh->probeN, probeR, probeS, Vprobe);
  

    mesh->probeI = (dfloat *) calloc(mesh->probeN*mesh->Np, sizeof(dfloat));

    for(int r=0; r<mesh->probeN; r++){
      for(int c=0; c<mesh->Np; c++){
      	dfloat s = 0;
      	for(int i=0; i<mesh->Np; i++){
      	  s += Vprobe[r*mesh->Np+i]*dV[i*mesh->Np + c];
      	}
      	mesh->probeI[r*mesh->Np + c] = s;
      }

    }

    free(V);
    free(dV);
    free(WORK);

  }

  free(IPIV);
  free(IPIV2);
  free(probeR);
  free(probeS);
  free(A);
  free(b);
  free(q);

}


void meshVandermonde2D(int N, int Npoints, dfloat *r, dfloat *s, dfloat *V){

  // First convert to ab coordinates
  dfloat *a = (dfloat *) calloc(Npoints, sizeof(dfloat));
  dfloat *b = (dfloat *) calloc(Npoints, sizeof(dfloat));
  for(int n=0; n<Npoints; n++){

    if(fabs(s[n]-1.0)>1e-8)
      a[n] = 2.0*(1.+r[n])/(1.0-s[n])-1.0;
    else
      a[n] = -1.0; 

    b[n] = s[n];

  }
  
  int sk=0;

  int Np = (N+1)*(N+2)/2; 

  for(int i=0; i<=N; i++){
    for(int j=0; j<=N-i; j++){
      for(int n=0; n<Npoints; n++){
        V[n*Np + sk] = meshSimplex2D(a[n], b[n], i, j);
      }
      sk++;
    }
  }


  free(a);
  free(b); 

}



dfloat meshSimplex2D(dfloat a, dfloat b, int i, int j){
  // 
  dfloat p1 = meshJacobiP(a,0,0,i);
  dfloat p2 = meshJacobiP(b,2*i+1,0,j);
  dfloat P = sqrt(2.0)*p1*p2*pow(1-b,i);

  return P; 
}


dfloat meshJacobiP(dfloat a, dfloat alpha, dfloat beta, int N){

  dfloat ax = a; 

  dfloat *P = (dfloat *) calloc((N+1), sizeof(dfloat));

  // Zero order
  dfloat gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*meshFactorial(alpha)*meshFactorial(beta)/meshFactorial(alpha+beta);
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


dfloat meshFactorial(int n){

  if(n==0)
    return 1;
  else
    return n*meshFactorial(n-1);
}
