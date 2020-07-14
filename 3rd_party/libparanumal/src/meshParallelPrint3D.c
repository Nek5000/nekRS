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

#include "mesh3D.h"

void meshParallelPrint3D(mesh3D *mesh){

  int rank, size;
  rank = mesh->rank;
  size = mesh->size;

  printf("rank %d: Nelements=" dlongFormat " Nnodes=" hlongFormat "\n", 
         rank, mesh->Nelements, mesh->Nnodes);
  
#if 0
  printf("EToV:\n");
  for(int e=0;e<mesh->Nelements;++e){
    printf("%d %d %d\n", 
           mesh->EToV[e*mesh->Nverts+0],
           mesh->EToV[e*mesh->Nverts+1],
           mesh->EToV[e*mesh->Nverts+2]);
  }
#endif

  dlong *otherNelements = (dlong*) calloc(size, sizeof(dlong));
  MPI_Allgather(&(mesh->Nelements), 1, MPI_DLONG,
                    otherNelements, 1, MPI_DLONG, 
                    mesh->comm);
  
  hlong *elementStarts = (hlong*) calloc(size, sizeof(hlong));
  for(int r=1;r<size;++r){
    elementStarts[r] = elementStarts[r-1]+otherNelements[r-1];
  }

  for(int r1=0;r1<size;++r1){
    MPI_Barrier(mesh->comm);
    if(rank==r1){
      fflush(stdout);
      if(r1==0)
        printf("EToE:\n");
      for(dlong e1=0;e1<mesh->Nelements;++e1){
        dlong id = e1*mesh->Nfaces;
        for(int f1=0;f1<mesh->Nfaces;++f1){
          hlong e2 = (hlong) mesh->EToE[id+f1];
          int f2 = mesh->EToF[id+f1];
          int r2 = mesh->EToP[id+f1];
          if(e2==-1 || f2==-1) 
            printf("(" hlongFormat " " "%d" ")=>X (" hlongFormat "," "%d" ")\n", 
                   e1+elementStarts[r1], f1, e2, f2);
          else{
            
            if(r2!=-1)
              e2 += elementStarts[r2];
            else
              e2 += elementStarts[r1];
            
            
            printf("(" hlongFormat " " "%d" ")=>(" hlongFormat " " "%d" ")\n", 
                   e1+elementStarts[r1], f1, e2, f2);
          }
        }
        fflush(stdout);
      }
    }
    MPI_Barrier(mesh->comm);
  }
  free(otherNelements);
  free(elementStarts);
}
