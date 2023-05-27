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
#include "platform.hpp"

struct parallelFace_t
{
  hlong v[4]; // vertices on face
  dlong element, elementN;
  int NfaceVertices;
  int face, rank;    // face info
  int faceN, rankN; // N for neighbor face info
};

// comparison function that orders vertices
// based on their combined vertex indices
int parallelCompareVertices(const void* a,
                            const void* b)
{
  parallelFace_t* fa = (parallelFace_t*) a;
  parallelFace_t* fb = (parallelFace_t*) b;

  for(int n = 0; n < fa->NfaceVertices; ++n) {
    if(fa->v[n] < fb->v[n]) return -1;
    if(fa->v[n] > fb->v[n]) return +1;
  }

  return 0;
}

/* comparison function that orders element/face
   based on their indexes */
int parallelCompareFaces(const void* a,
                         const void* b)
{
  parallelFace_t* fa = (parallelFace_t*) a;
  parallelFace_t* fb = (parallelFace_t*) b;

  if(fa->rank < fb->rank) return -1;
  if(fa->rank > fb->rank) return +1;

  if(fa->element < fb->element) return -1;
  if(fa->element > fb->element) return +1;

  if(fa->face < fb->face) return -1;
  if(fa->face > fb->face) return +1;

  return 0;
}

// mesh is the local partition
void meshParallelConnect(mesh_t* mesh)
{
  
  int rank, size;
  rank = platform->comm.mpiRank;
  size = platform->comm.mpiCommSize;

  //MPI_Barrier(platform->comm.mpiComm);
  //const double tStart = MPI_Wtime();
  //if(platform->comm.mpiRank == 0) printf("Building parallel face connectivity ... ");

  // serial connectivity on each process
  meshConnect(mesh);

  // count # of elements to send to each rank based on
  // minimum {vertex id % size}
  int* Nsend = (int*) calloc(size, sizeof(int));
  int* Nrecv = (int*) calloc(size, sizeof(int));
  int* sendOffsets = (int*) calloc(size, sizeof(int));
  int* recvOffsets = (int*) calloc(size, sizeof(int));

  // WARNING: In some corner cases, the number of faces to send may overrun int storage
  int allNsend = 0;
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int f = 0; f < mesh->Nfaces; ++f)
      if(mesh->EToE[e * mesh->Nfaces + f] == -1) {
        // find rank of destination for sorting based on max(face vertices)%size
        hlong maxv = 0;
        for(int n = 0; n < mesh->NfaceVertices; ++n) {
          int nid = mesh->faceVertices[f * mesh->NfaceVertices + n];
          hlong id = mesh->EToV[e * mesh->Nverts + nid];
          maxv = std::max(maxv, id);
        }
        int destRank = (int) (maxv % size);

        // increment send size for
        ++Nsend[destRank];
        ++allNsend;
      }

  // find send offsets
  for(int r = 1; r < size; ++r)
    sendOffsets[r] = sendOffsets[r - 1] + Nsend[r - 1];

  // reset counters
  for(int r = 0; r < size; ++r)
    Nsend[r] = 0;

  // buffer for outgoing data
  parallelFace_t* sendFaces = (parallelFace_t*) calloc(allNsend, sizeof(parallelFace_t));

  // Make the MPI_PARALLELFACE_T data type
  MPI_Datatype MPI_PARALLELFACE_T;
  MPI_Datatype dtype[8] = {MPI_HLONG, MPI_DLONG, MPI_DLONG, MPI_INT,
                           MPI_INT, MPI_INT, MPI_INT, MPI_INT};
  int blength[8] = {4, 1, 1, 1, 1, 1, 1, 1};
  MPI_Aint addr[8], displ[8];
  MPI_Get_address ( &(sendFaces[0]              ), addr + 0);
  MPI_Get_address ( &(sendFaces[0].element      ), addr + 1);
  MPI_Get_address ( &(sendFaces[0].elementN     ), addr + 2);
  MPI_Get_address ( &(sendFaces[0].NfaceVertices), addr + 3);
  MPI_Get_address ( &(sendFaces[0].face         ), addr + 4);
  MPI_Get_address ( &(sendFaces[0].rank         ), addr + 5);
  MPI_Get_address ( &(sendFaces[0].faceN        ), addr + 6);
  MPI_Get_address ( &(sendFaces[0].rankN        ), addr + 7);
  displ[0] = 0;
  displ[1] = addr[1] - addr[0];
  displ[2] = addr[2] - addr[0];
  displ[3] = addr[3] - addr[0];
  displ[4] = addr[4] - addr[0];
  displ[5] = addr[5] - addr[0];
  displ[6] = addr[6] - addr[0];
  displ[7] = addr[7] - addr[0];
  MPI_Type_create_struct (8, blength, displ, dtype, &MPI_PARALLELFACE_T);
  MPI_Type_commit (&MPI_PARALLELFACE_T);

  // pack face data
  for(dlong e = 0; e < mesh->Nelements; ++e)
    for(int f = 0; f < mesh->Nfaces; ++f)
      if(mesh->EToE[e * mesh->Nfaces + f] == -1) {
        // find rank of destination for sorting based on max(face vertices)%size
        hlong maxv = 0;
        for(int n = 0; n < mesh->NfaceVertices; ++n) {
          int nid = mesh->faceVertices[f * mesh->NfaceVertices + n];
          hlong id = mesh->EToV[e * mesh->Nverts + nid];
          maxv = std::max(maxv, id);
        }
        int destRank = (int) (maxv % size);

        // populate face to send out staged in segment of sendFaces array
        int id = sendOffsets[destRank] + Nsend[destRank];

        sendFaces[id].element = e;
        sendFaces[id].face = f;
        for(int n = 0; n < mesh->NfaceVertices; ++n) {
          int nid = mesh->faceVertices[f * mesh->NfaceVertices + n];
          sendFaces[id].v[n] = mesh->EToV[e * mesh->Nverts + nid];
        }

        mysort(sendFaces[id].v,mesh->NfaceVertices, "descending");

        sendFaces[id].NfaceVertices = mesh->NfaceVertices;
        sendFaces[id].rank = rank;

        sendFaces[id].elementN = -1;
        sendFaces[id].faceN = -1;
        sendFaces[id].rankN = -1;

        ++Nsend[destRank];
      }

  // exchange byte counts
  MPI_Alltoall(Nsend, 1, MPI_INT,
               Nrecv, 1, MPI_INT,
               platform->comm.mpiComm);

  // count incoming faces
  int allNrecv = 0;
  for(int r = 0; r < size; ++r)
    allNrecv += Nrecv[r];

  // find offsets for recv data
  for(int r = 1; r < size; ++r)
    recvOffsets[r] = recvOffsets[r - 1] + Nrecv[r - 1]; // byte offsets

  // buffer for incoming face data
  parallelFace_t* recvFaces = (parallelFace_t*) calloc(allNrecv, sizeof(parallelFace_t));

  // exchange parallel faces
  MPI_Alltoallv(sendFaces, Nsend, sendOffsets, MPI_PARALLELFACE_T,
                recvFaces, Nrecv, recvOffsets, MPI_PARALLELFACE_T,
                platform->comm.mpiComm);

  // local sort allNrecv received faces
  qsort(recvFaces, allNrecv, sizeof(parallelFace_t), parallelCompareVertices);

  // find matches
  for(int n = 0; n < allNrecv - 1; ++n)
    // since vertices are ordered we just look for pairs
    if(!parallelCompareVertices(recvFaces + n, recvFaces + n + 1)) {
      recvFaces[n].elementN = recvFaces[n + 1].element;
      recvFaces[n].faceN = recvFaces[n + 1].face;
      recvFaces[n].rankN = recvFaces[n + 1].rank;

      recvFaces[n + 1].elementN = recvFaces[n].element;
      recvFaces[n + 1].faceN = recvFaces[n].face;
      recvFaces[n + 1].rankN = recvFaces[n].rank;
    }

  // sort back to original ordering
  qsort(recvFaces, allNrecv, sizeof(parallelFace_t), parallelCompareFaces);

  // send faces back from whence they came
  MPI_Alltoallv(recvFaces, Nrecv, recvOffsets, MPI_PARALLELFACE_T,
                sendFaces, Nsend, sendOffsets, MPI_PARALLELFACE_T,
                platform->comm.mpiComm);

  // extract connectivity info
  mesh->EToP = (int*) calloc(mesh->Nelements * mesh->Nfaces, sizeof(int));
  for(dlong cnt = 0; cnt < mesh->Nelements * mesh->Nfaces; ++cnt)
    mesh->EToP[cnt] = -1;

  for(int cnt = 0; cnt < allNsend; ++cnt) {
    dlong e = sendFaces[cnt].element;
    dlong eN = sendFaces[cnt].elementN;
    int f = sendFaces[cnt].face;
    int fN = sendFaces[cnt].faceN;
    int rN = sendFaces[cnt].rankN;

    if(e >= 0 && f >= 0 && eN >= 0 && fN >= 0) {
      mesh->EToE[e * mesh->Nfaces + f] = eN;
      mesh->EToF[e * mesh->Nfaces + f] = fN;
      mesh->EToP[e * mesh->Nfaces + f] = rN;
    }
  }

  MPI_Barrier(platform->comm.mpiComm);
  MPI_Type_free(&MPI_PARALLELFACE_T);
  free(sendFaces);
  free(recvFaces);

  //MPI_Barrier(platform->comm.mpiComm);
  //if(platform->comm.mpiRank == 0) printf("done (%gs)\n", MPI_Wtime() - tStart);
}
