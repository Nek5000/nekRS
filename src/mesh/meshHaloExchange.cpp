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

#include "mesh.h"
#include "platform.hpp"
#include <cstring>

// send data from partition boundary elements
// and receive data to ghost elements
void meshHaloExchange(mesh_t* mesh,
                      size_t Nbytes, // message size per element
                      void* sourceBuffer,
                      void* sendBuffer, // temporary buffer
                      void* recvBuffer)
{
  
  // MPI info
  int rank, size;
  rank = platform->comm.mpiRank;
  size = platform->comm.mpiCommSize;

  // count outgoing and incoming meshes
  int tag = 999;

  // copy data from outgoing elements into temporary send buffer
  for(int i = 0; i < mesh->totalHaloPairs; ++i) {
    // outgoing element
    int e = mesh->haloElementList[i];
    // copy element e data to sendBuffer
    memcpy(((char*)sendBuffer) + i * Nbytes, ((char*)sourceBuffer) + e * Nbytes, Nbytes);
  }

  // initiate immediate send  and receives to each other process as needed
  int offset = 0, message = 0;
  for(int r = 0; r < size; ++r)
    if(r != rank) {
      size_t count = mesh->NhaloPairs[r] * Nbytes;
      if(count) {
        //	printf("rank %d sending %d bytes to rank %d\n", rank, count, r);
        MPI_Irecv(((char*)recvBuffer) + offset, count, MPI_CHAR, r, tag,
                  platform->comm.mpiComm, (MPI_Request*)mesh->haloRecvRequests + message);

        MPI_Isend(((char*)sendBuffer) + offset, count, MPI_CHAR, r, tag,
                  platform->comm.mpiComm, (MPI_Request*)mesh->haloSendRequests + message);
        offset += count;
        ++message;
      }
    }

  //  printf("mesh->NhaloMessages = %d\n", mesh->NhaloMessages);

  // Wait for all sent messages to have left and received messages to have arrived
  MPI_Status* sendStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
  MPI_Status* recvStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));

  MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloRecvRequests, recvStatus);
  MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloSendRequests, sendStatus);

  free(recvStatus);
  free(sendStatus);
}

// start halo exchange (for q)
void meshHaloExchangeStart(mesh_t* mesh,
                           size_t Nbytes, // message size per element
                           void* sendBuffer, // temporary buffer
                           void* recvBuffer)
{
  
  if(mesh->totalHaloPairs > 0) {
    // MPI info
    //    int rank, size;
    //    MPI_Comm_rank(platform->comm.mpiComm, &rank);
    //    MPI_Comm_size(platform->comm.mpiComm, &size);

    int rank = platform->comm.mpiRank;
    int size = platform->comm.mpiCommSize;

    // count outgoing and incoming meshes
    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int offset = 0, message = 0;
    for(int r = 0; r < size; ++r)
      if(r != rank) {
        size_t count = mesh->NhaloPairs[r] * Nbytes;
        if(count) {
          MPI_Irecv(((char*)recvBuffer) + offset, count, MPI_CHAR, r, tag,
                    platform->comm.mpiComm, ((MPI_Request*)mesh->haloRecvRequests) + message);

          MPI_Isend(((char*)sendBuffer) + offset, count, MPI_CHAR, r, tag,
                    platform->comm.mpiComm, ((MPI_Request*)mesh->haloSendRequests) + message);
          offset += count;
          ++message;
        }
      }
  }
}

void meshHaloExchangeFinish(mesh_t* mesh)
{
  if(mesh->totalHaloPairs > 0) {
    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Status* sendStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));
    MPI_Status* recvStatus = (MPI_Status*) calloc(mesh->NhaloMessages, sizeof(MPI_Status));

    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloRecvRequests, recvStatus);
    MPI_Waitall(mesh->NhaloMessages, (MPI_Request*)mesh->haloSendRequests, sendStatus);

    free(recvStatus);
    free(sendStatus);
  }
}

// start halo exchange (for q)
void meshHaloExchange(mesh_t* mesh,
                      size_t Nbytes, // message size per element
                      void* sendBuffer, // temporary buffer
                      void* recvBuffer)
{
  
  if(mesh->totalHaloPairs > 0) {
    int rank = platform->comm.mpiRank;
    int size = platform->comm.mpiCommSize;

    // count outgoing and incoming meshes
    int tag = 999;

    int Nmessages = 0;

    for(int r = 0; r < size; ++r)
      if(r != rank) {
        size_t count = mesh->NhaloPairs[r] * Nbytes;
        if(count)
          ++Nmessages;
      }

    MPI_Request* sendRequests = (MPI_Request*) calloc(Nmessages, sizeof(MPI_Request));
    MPI_Request* recvRequests = (MPI_Request*) calloc(Nmessages, sizeof(MPI_Request));

    // initiate immediate send  and receives to each other process as needed
    int offset = 0;
    Nmessages = 0;
    for(int r = 0; r < size; ++r)
      if(r != rank) {
        size_t count = mesh->NhaloPairs[r] * Nbytes;
        if(count) {
          MPI_Irecv(((char*)recvBuffer) + offset,
                    count,
                    MPI_CHAR,
                    r,
                    tag,
                    platform->comm.mpiComm,
                    (recvRequests) + Nmessages);
          MPI_Isend(((char*)sendBuffer) + offset,
                    count,
                    MPI_CHAR,
                    r,
                    tag,
                    platform->comm.mpiComm,
                    (sendRequests) + Nmessages);
          offset += count;
          ++Nmessages;
        }
      }

    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Status* sendStatus = (MPI_Status*) calloc(Nmessages, sizeof(MPI_Status));
    MPI_Status* recvStatus = (MPI_Status*) calloc(Nmessages, sizeof(MPI_Status));

    MPI_Waitall(Nmessages, recvRequests, recvStatus);
    MPI_Waitall(Nmessages, sendRequests, sendStatus);

    free(recvStatus);
    free(sendStatus);

    free(recvRequests);
    free(sendRequests);
  }
}
