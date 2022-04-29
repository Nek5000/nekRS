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
#include <string.h>

#include "nrssys.hpp"

void mergeLists(size_t sz,
                int N1, char* v1,
                int N2, char* v2,
                char* v3,
                int (* compare)(const void*, const void*),
                void (* match)(void*, void*))
{
  int n1 = 0, n2 = 0, n3 = 0;

  // merge two lists from v1 and v2
  for(n3 = 0; n3 < N1 + N2; ++n3) {
    if(n1 < N1 && n2 < N2) {
      int c = compare(v1 + n1 * sz,v2 + n2 * sz);
      if(c == -1) {
        memcpy(v3 + n3 * sz, v1 + n1 * sz, sz);
        ++n1;
      }else {
        memcpy(v3 + n3 * sz, v2 + n2 * sz, sz);
        ++n2;
      }
    }else if(n1 < N1) {
      memcpy(v3 + n3 * sz, v1 + n1 * sz, sz);
      ++n1;
    }else if(n2 < N2) {
      memcpy(v3 + n3 * sz, v2 + n2 * sz, sz);
      ++n2;
    }
  }

  // scan for matches
  for(n3 = 0; n3 < N1 + N2 - 1; ++n3)
    if(!compare(v3 + n3 * sz,v3 + (n3 + 1) * sz))
      match(v3 + n3 * sz, v3 + (n3 + 1) * sz);

  /* copy result back to v1, v2 */
  memcpy(v1, v3,       N1 * sz);
  memcpy(v2, v3 + sz * N1, N2 * sz);
}

// assumes N is even and the same on all ranks
void parallelSort(int size, int rank, MPI_Comm comm,
                  int N, void* vv, size_t sz,
                  int (* compare)(const void*, const void*),
                  void (* match)(void*, void*)
                  )
{
  /* cast void * to char * */
  char* v = (char*) vv;

  /* sort faces by their vertex number pairs */
  qsort(v, N, sz, compare);

  /* now do progressive merges */
  int NA = N / 2, NB = N / 2, NC = N / 2;

  MPI_Request recvA, recvC;
  MPI_Request sendA, sendC;
  MPI_Status status;
  int tag = 999;

  /* temporary buffer for incoming data */
  void* A = (void*) calloc(NA, sz);
  void* B = v;
  void* C = v + NB * sz;

  /* temporary space for merge sort */
  void* tmp = (void*) calloc(N, sz);

  /* max and min elements out of place hop one process at each step */
  for(int step = 0; step < size - 1; ++step) {
    /* send C, receive A */
    if(rank < size - 1)
      MPI_Isend(C, NC * sz, MPI_CHAR,  rank + 1, tag, comm, &sendC);
    if(rank > 0)
      MPI_Irecv(A, NA * sz, MPI_CHAR,  rank - 1, tag, comm, &recvA);

    if(rank < size - 1)
      MPI_Wait(&sendC, &status);
    if(rank > 0)
      MPI_Wait(&recvA, &status);

    /* merge sort A & B */
    if(rank > 0)
      mergeLists(sz, NA, (char*)A, NB, (char*)B, (char*)tmp, compare, match);

    /* send A, receive C */
    if(rank > 0)
      MPI_Isend(A, NA * sz, MPI_CHAR, rank - 1, tag, comm, &sendA);
    if(rank < size - 1)
      MPI_Irecv(C, NC * sz, MPI_CHAR, rank + 1, tag, comm, &recvC);

    if(rank > 0)
      MPI_Wait(&sendA, &status);
    if(rank < size - 1)
      MPI_Wait(&recvC, &status);

    /* merge sort B & C */
    mergeLists(sz, NB, (char*)B, NC, (char*)C, (char*)tmp, compare, match);
  }

  free(tmp);
  free(A);
}
