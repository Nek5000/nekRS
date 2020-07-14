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

#include <stdlib.h>
#include <stdio.h>
#include "mesh.h"

// structure used to encode vertices that make 
// each face, the element/face indices, and
// the neighbor element/face indices (if any)
typedef struct{

  dlong element;
  int face;

  int NfaceVertices;
  
  hlong v[4]; // max number of face vertices

  int bctype;

}boundaryFace_t;

// comparison function that orders vertices 
// based on their combined vertex indices
int compareBoundaryFaces(const void *a, 
                         const void *b){

  boundaryFace_t *fa = (boundaryFace_t*) a;
  boundaryFace_t *fb = (boundaryFace_t*) b;

  for(int n=0;n<fa->NfaceVertices;++n){
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;

}


/* routine to find EToB (Element To Boundary)*/
void meshConnectBoundary(mesh_t *mesh){

  /* count number of boundary faces (i.e. not yet connected) */
  hlong bcnt = 0;
  for(dlong e=0;e<mesh->Nelements;++e)
    for(int f=0;f<mesh->Nfaces;++f)
      if(mesh->EToE[e*mesh->Nfaces+f]==-1) // || mesh->EToE[e*mesh->Nfaces+f]==e)
        ++bcnt;

#if 0
  printf("Nbf = %d\n", mesh->NboundaryFaces);
  printf("Nfv = %d\n", mesh->NfaceVertices);
  printf("bcnt = %d\n", bcnt);
  printf("Nelements = %d\n", mesh->Nelements);
#endif
  
  /* build list of boundary faces */
  boundaryFace_t *boundaryFaces = (boundaryFace_t*) calloc(bcnt+mesh->NboundaryFaces,
                                                           sizeof(boundaryFace_t));

  bcnt = 0; // reset counter
  for(dlong e=0;e<mesh->Nelements;++e){    
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToE[e*mesh->Nfaces+f]==-1) { 
        
        for(int n=0;n<mesh->NfaceVertices;++n){
          dlong vid = e*mesh->Nverts + mesh->faceVertices[f*mesh->NfaceVertices+n];
          boundaryFaces[bcnt].v[n] = mesh->EToV[vid];
        }
      
        mysort(boundaryFaces[bcnt].v,mesh->NfaceVertices, "descending");

        boundaryFaces[bcnt].NfaceVertices = mesh->NfaceVertices;
        boundaryFaces[bcnt].element = e;
        boundaryFaces[bcnt].face = f;
        boundaryFaces[bcnt].bctype = -1;
        ++bcnt;
      }
    }
  }
  
  /* add boundary info */
  for(hlong b=0;b<mesh->NboundaryFaces;++b){
    
    for(int n=0;n<mesh->NfaceVertices;++n)
      boundaryFaces[bcnt].v[n] = mesh->boundaryInfo[b*(mesh->NfaceVertices+1)+n+1];
    
    mysort(boundaryFaces[bcnt].v,mesh->NfaceVertices, "descending");

    boundaryFaces[bcnt].NfaceVertices = mesh->NfaceVertices;
    boundaryFaces[bcnt].element = -1;
    boundaryFaces[bcnt].face = -1;
    boundaryFaces[bcnt].bctype = mesh->boundaryInfo[b*(mesh->NfaceVertices+1)];

    ++bcnt;
  }

#if 0
  for(int b=0;b<bcnt;++b){
    printf("%d: e=%d, f=%d, bc=%d, v=",
           b,
           boundaryFaces[b].element,
           boundaryFaces[b].face,
           boundaryFaces[b].bctype);
    for(int n=0;n<mesh->NfaceVertices;++n)
      printf("%d ", boundaryFaces[b].v[n]);
    printf("\n");
  }
#endif
  
  /* sort boundaryFaces by their vertex number pairs */
  qsort(boundaryFaces, bcnt, sizeof(boundaryFace_t), compareBoundaryFaces);

  /* scan through sorted face lists looking for element-boundary matches */
  mesh->EToB = (int*) calloc(mesh->Nelements*mesh->Nfaces, sizeof(int));
  for(dlong n=0;n<mesh->Nelements*mesh->Nfaces;++n) mesh->EToB[n] = -1;

  for(hlong cnt=0;cnt<bcnt-1;++cnt){
    if(!compareBoundaryFaces(boundaryFaces+cnt, boundaryFaces+cnt+1)){
      dlong e = mymax(boundaryFaces[cnt].element, boundaryFaces[cnt+1].element);
      int f   = mymax(boundaryFaces[cnt].face,    boundaryFaces[cnt+1].face);

      mesh->EToB[e*mesh->Nfaces+f] =
        mymax(boundaryFaces[cnt].bctype, boundaryFaces[cnt+1].bctype);
    }
  }

#if 0
  int cnt = 0;
  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      printf("EToE(%d,%d) = %d \n", e,f, mesh->EToE[cnt]);
      ++cnt;
    }
  }
#endif

  free(boundaryFaces);
}

