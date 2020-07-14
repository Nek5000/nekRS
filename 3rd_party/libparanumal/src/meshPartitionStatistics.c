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
#include "mesh.h"

void meshPartitionStatistics(mesh_t *mesh){

  /* get MPI rank and size */
  int rank, size;
  rank = mesh->rank;
  size = mesh->size;
  
  /* now gather statistics on connectivity between processes */
  int *comms = (int*) calloc(size, sizeof(int));
  int Ncomms = 0;

  /* count elements with neighbors on each other rank ranks */
  for(dlong e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      if(mesh->EToP[e*mesh->Nfaces+f]!=-1){
        ++comms[mesh->EToP[e*mesh->Nfaces+f]];
        ++Ncomms;
      }
    }
  }

  int Nmessages = 0;
  for(int r=0;r<size;++r)
    if(comms[r]>0)
      ++Nmessages;

  for(int r=0;r<size;++r){
    MPI_Barrier(mesh->comm);
    if(r==rank){
      fflush(stdout);
      printf("r: %02d [", rank);
      for(int s=0;s<size;++s){
        printf(" %04d", comms[s]);
      }
      printf("] (Nelements=" dlongFormat ", Nmessages=%d, Ncomms=%d)\n", mesh->Nelements,Nmessages, Ncomms);
      fflush(stdout);
    }
  }
  
  free(comms);
}
