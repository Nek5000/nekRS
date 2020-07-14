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

int parallelCompareRowColumn(const void *a, const void *b);

void ellipticBuildIpdgTri2D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);
void ellipticBuildIpdgTri3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);
void ellipticBuildIpdgQuad2D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);
void ellipticBuildIpdgQuad3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);
void ellipticBuildIpdgTet3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);
void ellipticBuildIpdgHex3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts);



void ellipticBuildIpdg(elliptic_t *elliptic, int basisNp, dfloat *basis,
                       dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  switch(elliptic->elementType){
  case TRIANGLES:{
    if(elliptic->dim==2)
      ellipticBuildIpdgTri2D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
    else
      ellipticBuildIpdgTri3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts);

    break;
  }
  case QUADRILATERALS:{
    if(elliptic->dim==2)
      ellipticBuildIpdgQuad2D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
    else
      ellipticBuildIpdgQuad3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); 
    break;
  }  
  case TETRAHEDRA:
    ellipticBuildIpdgTet3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); break;
  case HEXAHEDRA:
    ellipticBuildIpdgHex3D(elliptic, basisNp, basis,lambda, A, nnzA, globalStarts); break;
  }

}

void ellipticBuildIpdgTri2D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  if(!basis) { // default to degree N Lagrange basis
    basisNp = Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // number of degrees of freedom on this rank
  hlong Nnum = Np*Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
    rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = basisNp*basisNp*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(Nfaces*Nfp*Nfp,sizeof(dfloat));
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<Nfp;n++) {
      int fn = mesh->faceNodes[f*Nfp+n];

      for (int m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<Np;i++){
          MSnm += mesh->MM[fn+i*Np]*mesh->LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }


  // reset non-zero counter
  dlong nnz = 0;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  dfloat *SM = (dfloat*) calloc(Np*Np,sizeof(dfloat));
  dfloat *SP = (dfloat*) calloc(Np*Np,sizeof(dfloat));

  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // loop over all elements
  for(dlong eM=0;eM<Nelements;++eM){

    dlong vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        SM[n*Np+m]  = J*lambda*mesh->MM[n*Np+m];
        SM[n*Np+m] += J*drdx*drdx*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdx*dsdx*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

        SM[n*Np+m] += J*drdy*drdy*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdy*dsdy*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdy*dsdy*mesh->Sss[n*Np+m];
      }
    }

    for (int fM=0;fM<Nfaces;fM++) {

      for (int n=0;n<Np*Np;n++) SP[n] =0;

      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = elliptic->tau*hinv;

      dlong eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      dlong vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // reset eP
      eP = mesh->EToE[eM*Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP  = (int) (mesh->vmapP[idM]%Np);

          dfloat MSfnm = sJ*MSf[n*Nfp+m];

          SM[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          SP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];
            int iP = (int) (mesh->vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = drdx*mesh->Dr[iM*Np+m] + dsdx*mesh->Ds[iM*Np+m];
            dfloat DyMim = drdy*mesh->Dr[iM*Np+m] + dsdy*mesh->Ds[iM*Np+m];

            dfloat DxPim = drdxP*mesh->Dr[iP*Np+m] + dsdxP*mesh->Ds[iP*Np+m];
            dfloat DyPim = drdyP*mesh->Dr[iP*Np+m] + dsdyP*mesh->Ds[iP*Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            SM[nM*Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;

            SP[nM*Np+m] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            SP[nM*Np+m] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = (int) (mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];

            SM[n*Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;

            SP[n*Np+mP] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            SP[n*Np+mP] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
          }
        }
      }

      // store non-zeros for off diagonal block
      for(int j=0;j<basisNp;++j){
        for(int i=0;i<basisNp;++i){
          dfloat val = 0;
          for(int n=0;n<Np;++n){
            for(int m=0;m<Np;++m){
              val += basis[n*Np+j]*SP[n*Np+m]*basis[m*Np+i];
            }
          }

          if(fabs(val)>tol){
            (*A)[nnz].row = globalIds[eM*Np + j];
            (*A)[nnz].col = globalIds[eP*Np + i];
            (*A)[nnz].val = val;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
    // store non-zeros for diagonal block
    for(int j=0;j<basisNp;++j){
      for(int i=0;i<basisNp;++i){
        dfloat val = 0;
        for(int n=0;n<Np;++n){
          for(int m=0;m<Np;++m){
            val += basis[n*Np+j]*SM[n*Np+m]*basis[m*Np+i];
          }
        }

        if(fabs(val)>tol){
          (*A)[nnz].row = globalIds[eM*Np + j];
          (*A)[nnz].col = globalIds[eM*Np + i];
          (*A)[nnz].val = val;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  //printf("nnz = %d\n", nnz);

  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

#if 0
  dfloat* Ap = (dfloat *) calloc(Np*Np*Nelements*Nelements,sizeof(dfloat));
  for (int n=0;n<nnz;n++) {
    int row = (*A)[n].row;
    int col = (*A)[n].col;

    Ap[col+row*Np*Nelements] = (*A)[n].val;
  }

  for (int i=0;i<Np*Nelements;i++) {
    for (int j =0;j<Nelements*Np;j++) {
      printf("%4.2f \t", Ap[j+i*Np*Nelements]);
    }
    printf("\n");
  }
#endif

  free(globalIds);

  free(SM); free(SP);
  free(MS);
}

void ellipticBuildIpdgTri3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  if(!basis) { // default to degree N Lagrange basis
    basisNp = Np;
    basis = (dfloat*) calloc(basisNp*basisNp, sizeof(dfloat));
    for(int n=0;n<basisNp;++n){
      basis[n+n*basisNp] = 1;
    }
  }

  // number of degrees of freedom on this rank
  hlong Nnum = Np*Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
    rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = basisNp*basisNp*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(Nfaces*Nfp*Nfp,sizeof(dfloat));
  for (int f=0;f<Nfaces;f++) {
    for (int n=0;n<Nfp;n++) {
      int fn = mesh->faceNodes[f*Nfp+n];

      for (int m=0;m<Nfp;m++) {
        dfloat MSnm = 0;

        for (int i=0;i<Np;i++){
          MSnm += mesh->MM[fn+i*Np]*mesh->LIFT[i*Nfp*Nfaces+f*Nfp+m];
        }
        MS[m+n*Nfp + f*Nfp*Nfp]  = MSnm;
      }
    }
  }


  // reset non-zero counter
  dlong nnz = 0;

  *A = (nonZero_t*) calloc(nnzLocalBound, sizeof(nonZero_t));

  dfloat *SM = (dfloat*) calloc(Np*Np,sizeof(dfloat));
  dfloat *SP = (dfloat*) calloc(Np*Np,sizeof(dfloat));

  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // loop over all elements
  for(dlong eM=0;eM<Nelements;++eM){

    dlong vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat drdz = mesh->vgeo[vbase+RZID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat dsdz = mesh->vgeo[vbase+SZID];
    dfloat J = mesh->vgeo[vbase+JID];

    /* start with stiffness matrix  */
    for(int n=0;n<Np;++n){
      for(int m=0;m<Np;++m){
        SM[n*Np+m]  = J*lambda*mesh->MM[n*Np+m];
        SM[n*Np+m] += J*drdx*drdx*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdx*dsdx*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdx*drdx*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdx*dsdx*mesh->Sss[n*Np+m];

        SM[n*Np+m] += J*drdy*drdy*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdy*dsdy*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdy*drdy*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdy*dsdy*mesh->Sss[n*Np+m];

        SM[n*Np+m] += J*drdz*drdz*mesh->Srr[n*Np+m];
        SM[n*Np+m] += J*drdz*dsdz*mesh->Srs[n*Np+m];
        SM[n*Np+m] += J*dsdz*drdz*mesh->Ssr[n*Np+m];
        SM[n*Np+m] += J*dsdz*dsdz*mesh->Sss[n*Np+m];

      }
    }

    for (int fM=0;fM<Nfaces;fM++) {

      for (int n=0;n<Np*Np;n++) SP[n] =0;

      // load surface geofactors for this face
      dlong sid = mesh->Nsgeo*(eM*Nfaces+fM);
      dfloat nx = mesh->sgeo[sid+NXID];
      dfloat ny = mesh->sgeo[sid+NYID];
      dfloat nz = mesh->sgeo[sid+NZID];
      dfloat sJ = mesh->sgeo[sid+SJID];
      dfloat hinv = mesh->sgeo[sid+IHID];
      dfloat penalty = elliptic->tau*hinv;

      dlong eP = mesh->EToE[eM*Nfaces+fM];
      if (eP < 0) eP = eM;

      dlong vbaseP = eP*mesh->Nvgeo;
      dfloat drdxP = mesh->vgeo[vbaseP+RXID];
      dfloat drdyP = mesh->vgeo[vbaseP+RYID];
      dfloat drdzP = mesh->vgeo[vbaseP+RZID];
      dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
      dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
      dfloat dsdzP = mesh->vgeo[vbaseP+SZID];

      int bcD = 0, bcN =0;
      int bc = mesh->EToB[fM+Nfaces*eM]; //raw boundary flag
      int bcType = 0;

      if(bc>0) bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)

      // this needs to be double checked (and the code where these are used)
      if(bcType==1){ // Dirichlet
        bcD = 1;
        bcN = 0;
      } else if(bcType==2){ // Neumann
        bcD = 0;
        bcN = 1;
      }

      // reset eP
      eP = mesh->EToE[eM*Nfaces+fM];

      // mass matrix for this face
      dfloat *MSf = MS+fM*Nfp*Nfp;

      // penalty term just involves face nodes
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Nfp;++m){
          dlong idM = eM*Nfp*Nfaces+fM*Nfp+m;
          int nM = mesh->faceNodes[fM*Nfp+n];
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP  = (int) (mesh->vmapP[idM]%Np);

          dfloat MSfnm = sJ*MSf[n*Nfp+m];

          SM[nM*Np+mM] +=  0.5*(1.-bcN)*(1.+bcD)*penalty*MSfnm;
          SP[nM*Np+mP] += -0.5*(1.-bcN)*(1.-bcD)*penalty*MSfnm;
        }
      }

      // now add differential surface terms
      for(int n=0;n<Nfp;++n){
        for(int m=0;m<Np;++m){
          int nM = mesh->faceNodes[fM*Nfp+n];

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];
            int iP = (int) (mesh->vmapP[i + fM*Nfp+eM*Nfp*Nfaces]%Np);

            dfloat MSfni = sJ*MSf[n*Nfp+i]; // surface Jacobian built in

            dfloat DxMim = drdx*mesh->Dr[iM*Np+m] + dsdx*mesh->Ds[iM*Np+m];
            dfloat DyMim = drdy*mesh->Dr[iM*Np+m] + dsdy*mesh->Ds[iM*Np+m];
	    dfloat DzMim = drdz*mesh->Dr[iM*Np+m] + dsdz*mesh->Ds[iM*Np+m];
	    
            dfloat DxPim = drdxP*mesh->Dr[iP*Np+m] + dsdxP*mesh->Ds[iP*Np+m];
            dfloat DyPim = drdyP*mesh->Dr[iP*Np+m] + dsdyP*mesh->Ds[iP*Np+m];
	    dfloat DzPim = drdzP*mesh->Dr[iP*Np+m] + dsdzP*mesh->Ds[iP*Np+m];

            // OP11 = OP11 + 0.5*( - mmE*Dn1)
            SM[nM*Np+m] += -0.5*nx*(1+bcD)*(1-bcN)*MSfni*DxMim;
            SM[nM*Np+m] += -0.5*ny*(1+bcD)*(1-bcN)*MSfni*DyMim;
	    SM[nM*Np+m] += -0.5*nz*(1+bcD)*(1-bcN)*MSfni*DzMim;

            SP[nM*Np+m] += -0.5*nx*(1-bcD)*(1-bcN)*MSfni*DxPim;
            SP[nM*Np+m] += -0.5*ny*(1-bcD)*(1-bcN)*MSfni*DyPim;
	    SP[nM*Np+m] += -0.5*nz*(1-bcD)*(1-bcN)*MSfni*DzPim;
          }
        }
      }

      for(int n=0;n<Np;++n){
        for(int m=0;m<Nfp;++m){
          int mM = mesh->faceNodes[fM*Nfp+m];
          int mP = (int) (mesh->vmapP[m + fM*Nfp+eM*Nfp*Nfaces]%Np);

          for(int i=0;i<Nfp;++i){
            int iM = mesh->faceNodes[fM*Nfp+i];

            dfloat MSfim = sJ*MSf[i*Nfp+m];

            dfloat DxMin = drdx*mesh->Dr[iM*Np+n] + dsdx*mesh->Ds[iM*Np+n];
            dfloat DyMin = drdy*mesh->Dr[iM*Np+n] + dsdy*mesh->Ds[iM*Np+n];
	    dfloat DzMin = drdz*mesh->Dr[iM*Np+n] + dsdz*mesh->Ds[iM*Np+n];

            SM[n*Np+mM] +=  -0.5*nx*(1+bcD)*(1-bcN)*DxMin*MSfim;
            SM[n*Np+mM] +=  -0.5*ny*(1+bcD)*(1-bcN)*DyMin*MSfim;
	    SM[n*Np+mM] +=  -0.5*nz*(1+bcD)*(1-bcN)*DzMin*MSfim;

            SP[n*Np+mP] +=  +0.5*nx*(1-bcD)*(1-bcN)*DxMin*MSfim;
            SP[n*Np+mP] +=  +0.5*ny*(1-bcD)*(1-bcN)*DyMin*MSfim;
	    SP[n*Np+mP] +=  +0.5*nz*(1-bcD)*(1-bcN)*DzMin*MSfim;
          }
        }
      }

      // store non-zeros for off diagonal block
      for(int j=0;j<basisNp;++j){
        for(int i=0;i<basisNp;++i){
          dfloat val = 0;
          for(int n=0;n<Np;++n){
            for(int m=0;m<Np;++m){
              val += basis[n*Np+j]*SP[n*Np+m]*basis[m*Np+i];
            }
          }

          if(fabs(val)>tol){
            (*A)[nnz].row = globalIds[eM*Np + j];
            (*A)[nnz].col = globalIds[eP*Np + i];
            (*A)[nnz].val = val;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
    // store non-zeros for diagonal block
    for(int j=0;j<basisNp;++j){
      for(int i=0;i<basisNp;++i){
        dfloat val = 0;
        for(int n=0;n<Np;++n){
          for(int m=0;m<Np;++m){
            val += basis[n*Np+j]*SM[n*Np+m]*basis[m*Np+i];
          }
        }

        if(fabs(val)>tol){
          (*A)[nnz].row = globalIds[eM*Np + j];
          (*A)[nnz].col = globalIds[eM*Np + i];
          (*A)[nnz].val = val;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  //printf("nnz = %d\n", nnz);

  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

  free(globalIds);

  free(SM); free(SP);
  free(MS);
}



void ellipticBuildIpdgQuad2D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
    rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj=0;nj<mesh->N+1;++nj){
    for(int ni=0;ni<mesh->N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh->N+1;++j){
        for(int i=0;i<mesh->N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh->Np+node] = 1;
          if(nj==j)
            Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
          if(ni==i)
            Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
          
          ++node;
        }
      }
      
      ++mode;
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));
  
  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // reset non-zero counter
  dlong nnz = 0;
      
  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){
    
    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          dlong base = eM*mesh->Np*mesh->Nvgeo + i;
          dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
          dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
          dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
          dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
          dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
          
          int idn = n*mesh->Np+i;
          int idm = m*mesh->Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn];
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn];
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm];
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm];
          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];

            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP      = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh->vmapP[idM]/mesh->Np;
            dlong baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
            
            // grab surface geometric factors
            dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];
            
            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM; 
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM;
            dfloat lmM = B[idmM];
            
            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP;
            dfloat lmP = B[idmP];
            
            dfloat penalty = elliptic->tau*hinv;     

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP < 0) {
              int qSgn, gradqSgn;
              int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
              int bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              Anm += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              Anm += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              Anm += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-) 
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
          }
          if(fabs(AnmP)>tol){
            // remote info
            dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            (*A)[nnz].row = globalIds[eM*mesh->Np + n];
            (*A)[nnz].col = globalIds[eP*mesh->Np + m];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          } 
        }
        if(fabs(Anm)>tol){
          // local block
          (*A)[nnz].row = globalIds[eM*mesh->Np+n];
          (*A)[nnz].col = globalIds[eM*mesh->Np+m];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

  free(globalIds);
  free(B);  free(Br); free(Bs); 
}


void ellipticBuildIpdgQuad3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG, rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nj=0;nj<mesh->N+1;++nj){
    for(int ni=0;ni<mesh->N+1;++ni){

      int node = 0;

      for(int j=0;j<mesh->N+1;++j){
        for(int i=0;i<mesh->N+1;++i){

          if(nj==j && ni==i)
            B[mode*mesh->Np+node] = 1;
          if(nj==j)
            Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
          if(ni==i)
            Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
          
          ++node;
        }
      }
      
      ++mode;
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));
  
  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // reset non-zero counter
  dlong nnz = 0;
      
  // loop over all elements
  for(dlong eM=0;eM<mesh->Nelements;++eM){
    
    /* build Dx,Dy (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          dlong base = eM*mesh->Np*mesh->Nvgeo + i;
          dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
          dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
          dfloat drdz = mesh->vgeo[base+mesh->Np*RZID];
          dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
          dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
          dfloat dsdz = mesh->vgeo[base+mesh->Np*SZID];
          dfloat dtdx = mesh->vgeo[base+mesh->Np*TXID];
          dfloat dtdy = mesh->vgeo[base+mesh->Np*TYID];
          dfloat dtdz = mesh->vgeo[base+mesh->Np*TZID];
          dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
          
          int idn = n*mesh->Np+i;
          int idm = m*mesh->Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn]; + dtdx;
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn]; + dtdy;
          dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn]; + dtdz;
          					     	    
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm]; + dtdx;
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm]; + dtdy;
          dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm]; + dtdz;

	  Anm += JW*(dlndx*dlmdx+dlndy*dlmdy + dlndz*dlmdz);
	  Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
            
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
            dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];

            dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
            dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
            dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];


            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP      = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh->vmapP[idM]/mesh->Np;
            dlong baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat drdzP = mesh->vgeo[baseP+mesh->Np*RZID];
            
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
            dfloat dsdzP = mesh->vgeo[baseP+mesh->Np*SZID];

            dfloat dtdxP = mesh->vgeo[baseP+mesh->Np*TXID];
            dfloat dtdyP = mesh->vgeo[baseP+mesh->Np*TYID];
            dfloat dtdzP = mesh->vgeo[baseP+mesh->Np*TZID];
            
            // grab surface geometric factors
            dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat nz = mesh->sgeo[base+NZID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];
            
            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM; 
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM]; + dtdxM;
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM]; + dtdyM;
            dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM]; + dtdzM;

            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM + nz*dlndzM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM]; + dtdxM;
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM]; + dtdyM;
            dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM]; + dtdzM;
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM + nz*dlmdzM;
            dfloat lmM = B[idmM];
            
            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP]; + dtdxP;
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP]; + dtdyP;
            dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP]; + dtdzP;
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;
            dfloat lmP = B[idmP];
            
            dfloat penalty = elliptic->tau*hinv;     

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
            AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
            AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
          }
          if(fabs(AnmP)>tol){
            // remote info
            dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            (*A)[nnz].row = globalIds[eM*mesh->Np + n];
            (*A)[nnz].col = globalIds[eP*mesh->Np + m];
            (*A)[nnz].val = AnmP;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          } 
        }
	
        if(fabs(Anm)>tol){
          // local block
          (*A)[nnz].row = globalIds[eM*mesh->Np+n];
          (*A)[nnz].col = globalIds[eM*mesh->Np+m];
          (*A)[nnz].val = Anm;
          (*A)[nnz].ownerRank = rankM;
          ++nnz;
        }
      }
    }
  }

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

#if 0
  {
    FILE *fp = fopen("DGS.dat", "w");
    for(int n=0;n<nnz;++n){
      fprintf(fp, "%d %d %17.15lf\n",
	      (*A)[n].row+1,
	      (*A)[n].col+1,
	      (*A)[n].val);
    }
    fclose(fp);
  }
#endif
  
  free(globalIds);
  free(B);  free(Br); free(Bs); 
}






void ellipticBuildIpdgTet3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  // number of degrees of freedom on this rank
  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
    for(int r=0;r<mesh->size;++r)
      globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  dlong Nelements = mesh->Nelements;
  MPI_Allgather(&(mesh->Nelements), 1, MPI_DLONG,
                      rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<mesh->Nelements;e++) {
    for (int n=0;n<mesh->Np;n++) {
      globalIds[e*mesh->Np +n] = n + (e + rankStarts[rankM])*mesh->Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(mesh->Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, mesh->Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + mesh->Nelements*mesh->Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = mesh->Np*mesh->Np*(1+mesh->Nfaces)*mesh->Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // surface mass matrices MS = MM*LIFT
  dfloat *MS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Nfp;m++) {
        dfloat MSnm = 0;
        for (int i=0;i<mesh->Np;i++)
          MSnm += mesh->MM[n+i*mesh->Np]*mesh->LIFT[i*mesh->Nfp*mesh->Nfaces+f*mesh->Nfp+m];

        MS[m+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]  = MSnm;
      }
    }
  }

  // DrT*MS, DsT*MS, DtT*MS
  dfloat *DrTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DsTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  dfloat *DtTMS = (dfloat *) calloc(mesh->Nfaces*mesh->Np*mesh->Nfp,sizeof(dfloat));
  for (int f=0;f<mesh->Nfaces;f++) {
    for (int n=0;n<mesh->Np;n++) {
      for (int i=0;i<mesh->Nfp;i++) {
        DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        DtTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np] = 0.;
        for (int m=0;m<mesh->Np;m++) {
          DrTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dr[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DsTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Ds[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
          DtTMS[i+n*mesh->Nfp + f*mesh->Nfp*mesh->Np]
            += mesh->Dt[n+m*mesh->Np]*MS[i+m*mesh->Nfp+f*mesh->Nfp*mesh->Np];
        }
      }
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));

  // reset non-zero counter
  dlong nnz = 0;

  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // loop over all elements
  #pragma omp parallel
{

  dfloat *BM = (dfloat *) calloc(mesh->Np*mesh->Np,sizeof(dfloat));

  dfloat *qmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *qmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmM = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));
  dfloat *ndotgradqmP = (dfloat *) calloc(mesh->Nfp,sizeof(dfloat));

  #pragma omp for
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    dlong gbase = eM*mesh->Nggeo;
    dfloat Grr = mesh->ggeo[gbase+G00ID];
    dfloat Grs = mesh->ggeo[gbase+G01ID];
    dfloat Grt = mesh->ggeo[gbase+G02ID];
    dfloat Gss = mesh->ggeo[gbase+G11ID];
    dfloat Gst = mesh->ggeo[gbase+G12ID];
    dfloat Gtt = mesh->ggeo[gbase+G22ID];
    dfloat J   = mesh->ggeo[gbase+GWJID];

    /* start with stiffness matrix  */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){
        BM[m+n*mesh->Np]  = J*lambda*mesh->MM[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grr*mesh->Srr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Srs[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grt*mesh->Srt[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grs*mesh->Ssr[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gss*mesh->Sss[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gst*mesh->Sst[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Grt*mesh->Str[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gst*mesh->Sts[m+n*mesh->Np];
        BM[m+n*mesh->Np] += Gtt*mesh->Stt[m+n*mesh->Np];
      }
    }

    dlong vbase = eM*mesh->Nvgeo;
    dfloat drdx = mesh->vgeo[vbase+RXID];
    dfloat drdy = mesh->vgeo[vbase+RYID];
    dfloat drdz = mesh->vgeo[vbase+RZID];
    dfloat dsdx = mesh->vgeo[vbase+SXID];
    dfloat dsdy = mesh->vgeo[vbase+SYID];
    dfloat dsdz = mesh->vgeo[vbase+SZID];
    dfloat dtdx = mesh->vgeo[vbase+TXID];
    dfloat dtdy = mesh->vgeo[vbase+TYID];
    dfloat dtdz = mesh->vgeo[vbase+TZID];

    for (int m=0;m<mesh->Np;m++) {
      for (int fM=0;fM<mesh->Nfaces;fM++) {
        // load surface geofactors for this face
        dlong sid = mesh->Nsgeo*(eM*mesh->Nfaces+fM);
        dfloat nx = mesh->sgeo[sid+NXID];
        dfloat ny = mesh->sgeo[sid+NYID];
        dfloat nz = mesh->sgeo[sid+NZID];
        dfloat sJ = mesh->sgeo[sid+SJID];
        dfloat hinv = mesh->sgeo[sid+IHID];

        dlong eP = mesh->EToE[eM*mesh->Nfaces+fM];
        if (eP < 0) eP = eM;
        dlong vbaseP = eP*mesh->Nvgeo;
        dfloat drdxP = mesh->vgeo[vbaseP+RXID];
        dfloat drdyP = mesh->vgeo[vbaseP+RYID];
        dfloat drdzP = mesh->vgeo[vbaseP+RZID];
        dfloat dsdxP = mesh->vgeo[vbaseP+SXID];
        dfloat dsdyP = mesh->vgeo[vbaseP+SYID];
        dfloat dsdzP = mesh->vgeo[vbaseP+SZID];
        dfloat dtdxP = mesh->vgeo[vbaseP+TXID];
        dfloat dtdyP = mesh->vgeo[vbaseP+TYID];
        dfloat dtdzP = mesh->vgeo[vbaseP+TZID];

        // extract trace nodes
        for (int i=0;i<mesh->Nfp;i++) {
          // double check vol geometric factors are in halo storage of vgeo
          int idM    = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
          int vidM   = mesh->faceNodes[i+fM*mesh->Nfp];
          int vidP   = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo

          qmM[i] =0;
          if (vidM == m) qmM[i] =1;
          qmP[i] =0;
          if (vidP == m) qmP[i] =1;

          ndotgradqmM[i] = (nx*drdx+ny*drdy+nz*drdz)*mesh->Dr[m+vidM*mesh->Np]
                          +(nx*dsdx+ny*dsdy+nz*dsdz)*mesh->Ds[m+vidM*mesh->Np]
                          +(nx*dtdx+ny*dtdy+nz*dtdz)*mesh->Dt[m+vidM*mesh->Np];
          ndotgradqmP[i] = (nx*drdxP+ny*drdyP+nz*drdzP)*mesh->Dr[m+vidP*mesh->Np]
                          +(nx*dsdxP+ny*dsdyP+nz*dsdzP)*mesh->Ds[m+vidP*mesh->Np]
                          +(nx*dtdxP+ny*dtdyP+nz*dtdzP)*mesh->Dt[m+vidP*mesh->Np];
        }

        dfloat penalty = elliptic->tau*hinv;
        eP = mesh->EToE[eM*mesh->Nfaces+fM];
        for (int n=0;n<mesh->Np;n++) {
          for (int i=0;i<mesh->Nfp;i++) {
            BM[m+n*mesh->Np] += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
            BM[m+n*mesh->Np] += -0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                -0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                -0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
            BM[m+n*mesh->Np] += +0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
          }

          dfloat AnmP = 0;
          if (eP < 0) {
            int qSgn, gradqSgn;
            int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
            int bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)
            if(bcType==1){ // Dirichlet
              qSgn     = -1;
              gradqSgn =  1;
            } else if (bcType==2){ // Neumann
              qSgn     =  1;
              gradqSgn = -1;
            } else { // Neumann for now
              qSgn     =  1;
              gradqSgn = -1;
            }

            for (int i=0;i<mesh->Nfp;i++) {
              BM[m+n*mesh->Np] += -0.5*gradqSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmM[i];
              BM[m+n*mesh->Np] += +0.5*qSgn*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i]
                                  +0.5*qSgn*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmM[i];
              BM[m+n*mesh->Np] += -0.5*qSgn*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmM[i];
            }
          } else {
            for (int i=0;i<mesh->Nfp;i++) {
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*ndotgradqmP[i];
              AnmP += +0.5*sJ*(nx*drdx+ny*drdy+nz*drdz)*DrTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dsdx+ny*dsdy+nz*dsdz)*DsTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i]
                      +0.5*sJ*(nx*dtdx+ny*dtdy+nz*dtdz)*DtTMS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*qmP[i];
              AnmP += -0.5*sJ*MS[i+n*mesh->Nfp+fM*mesh->Nfp*mesh->Np]*penalty*qmP[i];
            }
          }

          if(fabs(AnmP)>tol){
            #pragma omp critical
            {
              // remote info
              (*A)[nnz].row = globalIds[eM*mesh->Np+n];
              (*A)[nnz].col = globalIds[eP*mesh->Np+m];
              (*A)[nnz].val = AnmP;
              (*A)[nnz].ownerRank = rankM;
              ++nnz;
            }
          }
        }
      }
    }

    for (int n=0;n<mesh->Np;n++) {
      for (int m=0;m<mesh->Np;m++) {
        dfloat Anm = BM[m+n*mesh->Np];

        if(fabs(Anm)>tol){
          #pragma omp critical
          {
            (*A)[nnz].row = globalIds[eM*mesh->Np+n];
            (*A)[nnz].col = globalIds[eM*mesh->Np+m];
            (*A)[nnz].val = Anm;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
  }
  
  free(BM);
  free(qmM); free(qmP);
  free(ndotgradqmM); free(ndotgradqmP);
}
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  // free up unused storage
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;
  
  if(rankM==0) printf("done.\n");
  
  free(globalIds);

  free(MS);
  free(DrTMS); free(DsTMS); free(DtTMS);
}

void ellipticBuildIpdgHex3D(elliptic_t *elliptic, int basisNp, dfloat *basis,
                            dfloat lambda, nonZero_t **A, dlong *nnzA, hlong *globalStarts){

  mesh_t *mesh = elliptic->mesh;
  setupAide options = elliptic->options;

  int rankM = mesh->rank;
  
  int Np = mesh->Np;
  int Nfp = mesh->Nfp;
  int Nfaces = mesh->Nfaces;
  dlong Nelements = mesh->Nelements;

  hlong Nnum = mesh->Np*mesh->Nelements;

  // create a global numbering system
  hlong *globalIds = (hlong *) calloc((Nelements+mesh->totalHaloPairs)*Np,sizeof(hlong));

  // every degree of freedom has its own global id
  MPI_Allgather(&Nnum, 1, MPI_HLONG, globalStarts+1, 1, MPI_HLONG, mesh->comm);
  for(int r=0;r<mesh->size;++r)
    globalStarts[r+1] = globalStarts[r]+globalStarts[r+1];

  /* so find number of elements on each rank */
  dlong *rankNelements = (dlong*) calloc(mesh->size, sizeof(dlong));
  hlong *rankStarts = (hlong*) calloc(mesh->size+1, sizeof(hlong));
  MPI_Allgather(&Nelements, 1, MPI_DLONG,
    rankNelements, 1, MPI_DLONG, mesh->comm);
  //find offsets
  for(int r=0;r<mesh->size;++r){
    rankStarts[r+1] = rankStarts[r]+rankNelements[r];
  }
  //use the offsets to set a global id
  for (dlong e =0;e<Nelements;e++) {
    for (int n=0;n<Np;n++) {
      globalIds[e*Np +n] = n + (e + rankStarts[rankM])*Np;
    }
  }

  /* do a halo exchange of global node numbers */
  if (mesh->totalHaloPairs) {
    hlong *idSendBuffer = (hlong *) calloc(Np*mesh->totalHaloPairs,sizeof(hlong));
    meshHaloExchange(mesh, Np*sizeof(hlong), globalIds, idSendBuffer, globalIds + Nelements*Np);
    free(idSendBuffer);
  }

  dlong nnzLocalBound = Np*Np*(1+Nfaces)*Nelements;

  // drop tolerance for entries in sparse storage
  dfloat tol = 1e-8;

  // build some monolithic basis arrays (use Dr,Ds,Dt and insert MM instead of weights for tet version)
  dfloat *B  = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Br = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bs = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));
  dfloat *Bt = (dfloat*) calloc(mesh->Np*mesh->Np, sizeof(dfloat));

  int mode = 0;
  for(int nk=0;nk<mesh->N+1;++nk){
    for(int nj=0;nj<mesh->N+1;++nj){
      for(int ni=0;ni<mesh->N+1;++ni){

        int node = 0;

        for(int k=0;k<mesh->N+1;++k){
          for(int j=0;j<mesh->N+1;++j){
            for(int i=0;i<mesh->N+1;++i){

              if(nk==k && nj==j && ni==i)
                B[mode*mesh->Np+node] = 1;
              if(nj==j && nk==k)
                Br[mode*mesh->Np+node] = mesh->D[ni+mesh->Nq*i]; 
              if(ni==i && nk==k)
                Bs[mode*mesh->Np+node] = mesh->D[nj+mesh->Nq*j]; 
              if(ni==i && nj==j)
                Bt[mode*mesh->Np+node] = mesh->D[nk+mesh->Nq*k]; 
              
              ++node;
            }
          }
        }
        
        ++mode;
      }
    }
  }

  *A = (nonZero_t*) calloc(nnzLocalBound,sizeof(nonZero_t));
  
  if(rankM==0) printf("Building full IPDG matrix...");fflush(stdout);

  // reset non-zero counter
  dlong nnz = 0;
  
  // loop over all elements
  //#pragma omp parallel for
  for(dlong eM=0;eM<mesh->Nelements;++eM){

    /* build Dx,Dy,Dz (forget the TP for the moment) */
    for(int n=0;n<mesh->Np;++n){
      for(int m=0;m<mesh->Np;++m){ // m will be the sub-block index for negative and positive trace
        dfloat Anm = 0;

        // (grad phi_n, grad phi_m)_{D^e}
        for(int i=0;i<mesh->Np;++i){
          dlong base = eM*mesh->Np*mesh->Nvgeo + i;
          dfloat drdx = mesh->vgeo[base+mesh->Np*RXID];
          dfloat drdy = mesh->vgeo[base+mesh->Np*RYID];
          dfloat drdz = mesh->vgeo[base+mesh->Np*RZID];
          dfloat dsdx = mesh->vgeo[base+mesh->Np*SXID];
          dfloat dsdy = mesh->vgeo[base+mesh->Np*SYID];
          dfloat dsdz = mesh->vgeo[base+mesh->Np*SZID];
          dfloat dtdx = mesh->vgeo[base+mesh->Np*TXID];
          dfloat dtdy = mesh->vgeo[base+mesh->Np*TYID];
          dfloat dtdz = mesh->vgeo[base+mesh->Np*TZID];
          dfloat JW   = mesh->vgeo[base+mesh->Np*JWID];
          
          int idn = n*mesh->Np+i;
          int idm = m*mesh->Np+i;
          dfloat dlndx = drdx*Br[idn] + dsdx*Bs[idn] + dtdx*Bt[idn];
          dfloat dlndy = drdy*Br[idn] + dsdy*Bs[idn] + dtdy*Bt[idn];
          dfloat dlndz = drdz*Br[idn] + dsdz*Bs[idn] + dtdz*Bt[idn];    
          dfloat dlmdx = drdx*Br[idm] + dsdx*Bs[idm] + dtdx*Bt[idm];
          dfloat dlmdy = drdy*Br[idm] + dsdy*Bs[idm] + dtdy*Bt[idm];
          dfloat dlmdz = drdz*Br[idm] + dsdz*Bs[idm] + dtdz*Bt[idm];
          Anm += JW*(dlndx*dlmdx+dlndy*dlmdy+dlndz*dlmdz);
          Anm += lambda*JW*B[idn]*B[idm];
        }

        // loop over all faces in this element
        for(int fM=0;fM<mesh->Nfaces;++fM){
          // accumulate flux terms for negative and positive traces
          dfloat AnmP = 0;
          for(int i=0;i<mesh->Nfp;++i){
            int vidM = mesh->faceNodes[i+fM*mesh->Nfp];

            // grab vol geofacs at surface nodes
            dlong baseM = eM*mesh->Np*mesh->Nvgeo + vidM;
            dfloat drdxM = mesh->vgeo[baseM+mesh->Np*RXID];
            dfloat drdyM = mesh->vgeo[baseM+mesh->Np*RYID];
            dfloat drdzM = mesh->vgeo[baseM+mesh->Np*RZID];
            dfloat dsdxM = mesh->vgeo[baseM+mesh->Np*SXID];
            dfloat dsdyM = mesh->vgeo[baseM+mesh->Np*SYID];
            dfloat dsdzM = mesh->vgeo[baseM+mesh->Np*SZID];
            dfloat dtdxM = mesh->vgeo[baseM+mesh->Np*TXID];
            dfloat dtdyM = mesh->vgeo[baseM+mesh->Np*TYID];
            dfloat dtdzM = mesh->vgeo[baseM+mesh->Np*TZID];

            // double check vol geometric factors are in halo storage of vgeo
            dlong idM     = eM*mesh->Nfp*mesh->Nfaces+fM*mesh->Nfp+i;
            int vidP    = (int) (mesh->vmapP[idM]%mesh->Np); // only use this to identify location of positive trace vgeo
            dlong localEP = mesh->vmapP[idM]/mesh->Np;
            dlong baseP   = localEP*mesh->Np*mesh->Nvgeo + vidP; // use local offset for vgeo in halo
            dfloat drdxP = mesh->vgeo[baseP+mesh->Np*RXID];
            dfloat drdyP = mesh->vgeo[baseP+mesh->Np*RYID];
            dfloat drdzP = mesh->vgeo[baseP+mesh->Np*RZID];
            dfloat dsdxP = mesh->vgeo[baseP+mesh->Np*SXID];
            dfloat dsdyP = mesh->vgeo[baseP+mesh->Np*SYID];
            dfloat dsdzP = mesh->vgeo[baseP+mesh->Np*SZID];
            dfloat dtdxP = mesh->vgeo[baseP+mesh->Np*TXID];
            dfloat dtdyP = mesh->vgeo[baseP+mesh->Np*TYID];
            dfloat dtdzP = mesh->vgeo[baseP+mesh->Np*TZID];
            
            // grab surface geometric factors
            dlong base = mesh->Nsgeo*(eM*mesh->Nfp*mesh->Nfaces + fM*mesh->Nfp + i);
            dfloat nx = mesh->sgeo[base+NXID];
            dfloat ny = mesh->sgeo[base+NYID];
            dfloat nz = mesh->sgeo[base+NZID];
            dfloat wsJ = mesh->sgeo[base+WSJID];
            dfloat hinv = mesh->sgeo[base+IHID];
            
            // form negative trace terms in IPDG
            int idnM = n*mesh->Np+vidM; 
            int idmM = m*mesh->Np+vidM;
            int idmP = m*mesh->Np+vidP;

            dfloat dlndxM = drdxM*Br[idnM] + dsdxM*Bs[idnM] + dtdxM*Bt[idnM];
            dfloat dlndyM = drdyM*Br[idnM] + dsdyM*Bs[idnM] + dtdyM*Bt[idnM];
            dfloat dlndzM = drdzM*Br[idnM] + dsdzM*Bs[idnM] + dtdzM*Bt[idnM];
            dfloat ndotgradlnM = nx*dlndxM+ny*dlndyM+nz*dlndzM;
            dfloat lnM = B[idnM];

            dfloat dlmdxM = drdxM*Br[idmM] + dsdxM*Bs[idmM] + dtdxM*Bt[idmM];
            dfloat dlmdyM = drdyM*Br[idmM] + dsdyM*Bs[idmM] + dtdyM*Bt[idmM];
            dfloat dlmdzM = drdzM*Br[idmM] + dsdzM*Bs[idmM] + dtdzM*Bt[idmM];
            dfloat ndotgradlmM = nx*dlmdxM+ny*dlmdyM+nz*dlmdzM;
            dfloat lmM = B[idmM];
            
            dfloat dlmdxP = drdxP*Br[idmP] + dsdxP*Bs[idmP] + dtdxP*Bt[idmP];
            dfloat dlmdyP = drdyP*Br[idmP] + dsdyP*Bs[idmP] + dtdyP*Bt[idmP];
            dfloat dlmdzP = drdzP*Br[idmP] + dsdzP*Bs[idmP] + dtdzP*Bt[idmP];
            dfloat ndotgradlmP = nx*dlmdxP+ny*dlmdyP+nz*dlmdzP;
            dfloat lmP = B[idmP];
            
            dfloat penalty = elliptic->tau*hinv;     

            Anm += -0.5*wsJ*lnM*ndotgradlmM;  // -(ln^-, N.grad lm^-)
            Anm += -0.5*wsJ*ndotgradlnM*lmM;  // -(N.grad ln^-, lm^-)
            Anm += +0.5*wsJ*penalty*lnM*lmM; // +((tau/h)*ln^-,lm^-)

            dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
            if (eP<0) {
              int qSgn, gradqSgn;
              int bc = mesh->EToB[fM+mesh->Nfaces*eM]; //raw boundary flag
              int bcType = elliptic->BCType[bc];          //find its type (Dirichlet/Neumann)
              if(bcType==1){ // Dirichlet
                qSgn     = -1;
                gradqSgn =  1;
              } else if (bcType==2){ // Neumann
                qSgn     =  1;
                gradqSgn = -1;
              } else { // Neumann for now
                qSgn     =  1;
                gradqSgn = -1;
              }

              Anm += -0.5*gradqSgn*wsJ*lnM*ndotgradlmM;  // -(ln^-, -N.grad lm^-)
              Anm += +0.5*qSgn*wsJ*ndotgradlnM*lmM;  // +(N.grad ln^-, lm^-)
              Anm += -0.5*qSgn*wsJ*penalty*lnM*lmM; // -((tau/h)*ln^-,lm^-) 
            } else {
              AnmP += -0.5*wsJ*lnM*ndotgradlmP;  // -(ln^-, N.grad lm^+)
              AnmP += +0.5*wsJ*ndotgradlnM*lmP;  // +(N.grad ln^-, lm^+)
              AnmP += -0.5*wsJ*penalty*lnM*lmP; // -((tau/h)*ln^-,lm^+)
            }
          }
          if(fabs(AnmP)>tol){
            //#pragma omp critical
            {
              // remote info
              dlong eP    = mesh->EToE[eM*mesh->Nfaces+fM];
              (*A)[nnz].row = globalIds[eM*mesh->Np + n];
              (*A)[nnz].col = globalIds[eP*mesh->Np + m];
              (*A)[nnz].val = AnmP;
              (*A)[nnz].ownerRank = rankM;
              ++nnz;
            }
          } 
        }
        if(fabs(Anm)>tol){
          //#pragma omp critical
          {
            // local block
            (*A)[nnz].row = globalIds[eM*mesh->Np+n];
            (*A)[nnz].col = globalIds[eM*mesh->Np+m];
            (*A)[nnz].val = Anm;
            (*A)[nnz].ownerRank = rankM;
            ++nnz;
          }
        }
      }
    }
  }

  // sort received non-zero entries by row block (may need to switch compareRowColumn tests)
  qsort((*A), nnz, sizeof(nonZero_t), parallelCompareRowColumn);
  //*A = (nonZero_t*) realloc(*A, nnz*sizeof(nonZero_t));
  *nnzA = nnz;

  if(rankM==0) printf("done.\n");

  free(globalIds);
  free(B);  free(Br); free(Bs); free(Bt);
}
