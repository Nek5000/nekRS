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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mesh2D.h"


// serial face-mode to face-mode connection
void meshConnectFaceModes2D(mesh2D *mesh, int *faceModes, dfloat *V){

  /* volume indices of the interior and exterior face modes for each element */
  mesh->mmapM = (dlong*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(dlong));
  mesh->mmapP = (dlong*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(dlong));
  mesh->mmapS = (int*) calloc(mesh->Nfp*mesh->Nfaces*mesh->Nelements, sizeof(int));

  dfloat *VM = (dfloat *) calloc(mesh->Np,sizeof(dfloat));
  dfloat *VP = (dfloat *) calloc(mesh->Np,sizeof(dfloat));

  /* assume elements already connected */
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      dlong eP = mesh->EToE[e*mesh->Nfaces+f];
      int fP = mesh->EToF[e*mesh->Nfaces+f];
      if(eP<0 || fP<0){ // fake connections for unconnected faces
        eP = e;
        fP = f;
      }
      
      /* for each mode on this face find the neighbor mode */
      for(int n=0;n<mesh->Nfp;++n){
        int m = faceModes[n+f*mesh->Nfp]; //get face mode number

        for (int i=0;i<mesh->Nfp;i++) {
          int k = mesh->faceNodes[i+f*mesh->Nfp];
          VM[i] = V[m+k*mesh->Np]; //evaluate mode at WB nodes on face  
        }
        
        dfloat mindist = 1E9;
        int s; 
        int mMatch;
        for (int nP=0;nP<mesh->Nfp;nP++) {
          //test the modes on face fP
          int mP = faceModes[nP+fP*mesh->Nfp];

          for (int i=0;i<mesh->Nfp;i++) {
            //get neighbouring node
            dlong id = i+f*mesh->Nfp+e*mesh->Nfp*mesh->Nfaces;
            int k = mesh->vmapP[id]%mesh->Np;

            VP[i] = V[mP+k*mesh->Np]; //evaluate mode at WB nodes on face     
          }

          dfloat dist1=0, dist2=0;
          for (int i=0;i<mesh->Nfp;i++){
            dist1 += pow(VM[i]-VP[i],2);
            dist2 += pow(VM[i]+VP[i],2);
          }
          dist1 = sqrt(dist1); 
          dist2 = sqrt(dist2);

          /* if next node is closer to target update match */
          if(dist1<mindist){
            mindist = dist1;
            mMatch = mP;
            s = 1;
          }
          if(dist2<mindist){
            mindist = dist2;
            mMatch = mP;
            s = -1;
          }
        }
        if(mindist>1e-3) printf("arggh - bad match: e="dlongFormat",f=%d, mode=%d\n", e,f, m);

        dlong id  = mesh->Nfaces*mesh->Nfp*e + f*mesh->Nfp + n;
        dlong idM = faceModes[f*mesh->Nfp+n] + e*mesh->Np;
        dlong idP = mMatch + eP*mesh->Np;

        mesh->mmapM[id] = idM;
        mesh->mmapP[id] = idP;
        mesh->mmapS[id] = s;
      }
    }
  }
}