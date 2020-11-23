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

typedef struct
{
  hlong vertex;
  hlong element;
  hlong rank;
  hlong rankN;    // neighbor rank
  hlong sortTag;
}vertex_t;

// generic comparator
int compareSortTag(const void* a,
                   const void* b)
{
  vertex_t* va = (vertex_t*) a;
  vertex_t* vb = (vertex_t*) b;

  if(va->sortTag < vb->sortTag) return -1;
  if(va->sortTag > vb->sortTag) return +1;

  return 0;
}

// use this to sort list of elements to send to each neighbor rank
int compareRankNElement(const void* a,
                        const void* b)
{
  vertex_t* va = (vertex_t*) a;
  vertex_t* vb = (vertex_t*) b;

  if(va->rankN < vb->rankN) return -1;
  if(va->rankN > vb->rankN) return +1;

  if(va->element < vb->element) return -1;
  if(va->element > vb->element) return +1;

  return 0;
}

// start one ring exchange (for q)
void ellipticOneRingExchangeStart(MPI_Comm &comm,
                                  size_t Nbytes, // message size per element
                                  hlong NoneRingSendTotal,
                                  int* NoneRingSend,
                                  void* sendBuffer, // temporary buffer
                                  MPI_Request* sendRequests,
                                  int* NsendMessages,
                                  hlong NoneRingRecvTotal,
                                  int* NoneRingRecv,
                                  void* recvBuffer,
                                  MPI_Request* recvRequests,
                                  int* NrecvMessages)
{
  // WATCH OUT - LOOPING OVER ALL RANKS BAD
  if(NoneRingRecvTotal + NoneRingSendTotal > 0) {
    // MPI info
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // count outgoing and incoming meshes
    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int sendOffset = 0, recvOffset = 0, sendMessage = 0, recvMessage = 0;
    for(int r = 0; r < size; ++r)
      if(r != rank) {
        size_t recvCount = NoneRingRecv[r] * Nbytes;
        if(recvCount) {
          MPI_Irecv(((char*)recvBuffer) + recvOffset, recvCount, MPI_CHAR, r, tag,
                    comm, recvRequests + recvMessage);
          recvOffset += recvCount;
          ++recvMessage;
        }

        size_t sendCount = NoneRingSend[r] * Nbytes;
        if(sendCount) {
          MPI_Isend(((char*)sendBuffer) + sendOffset, sendCount, MPI_CHAR, r, tag,
                    comm, sendRequests + sendMessage);

          sendOffset += sendCount;
          ++sendMessage;
        }
      }

    *NsendMessages = sendMessage;
    *NrecvMessages = recvMessage;
  }
}

void ellipticOneRingExchangeFinish(MPI_Comm &comm,
                                   int NsendMessages,
                                   MPI_Request* sendRequests,
                                   int NrecvMessages,
                                   MPI_Request* recvRequests)
{
  if(NrecvMessages) {
    // Wait for all sent messages to have left and received messages to have arrived
    MPI_Status* recvStatus = (MPI_Status*) calloc(NrecvMessages, sizeof(MPI_Status));

    MPI_Waitall(NrecvMessages, recvRequests, recvStatus);

    free(recvStatus);
  }

  if(NsendMessages) {
    MPI_Status* sendStatus = (MPI_Status*) calloc(NsendMessages, sizeof(MPI_Status));

    MPI_Waitall(NsendMessages, sendRequests, sendStatus);

    free(sendStatus);
  }
}

void ellipticOneRingExchange(MPI_Comm &comm,
                             hlong Nelements,
                             size_t Nbytes, // message size per element
                             void* q,
                             hlong NoneRingSendTotal,
                             hlong* oneRingSendList,
                             hlong* NoneRingSend,
                             void* sendBuffer, // temporary buffer
                             MPI_Request* sendRequests,
                             hlong NoneRingRecvTotal,
                             hlong* NoneRingRecv,
                             MPI_Request* recvRequests,
                             void* qOneRing)
{
  // do oneRing extract
  for(hlong n = 0; n < NoneRingSendTotal; ++n) {
    hlong e = oneRingSendList[n];
    memcpy((char*)sendBuffer + n * Nbytes, (char*)q + e * Nbytes, Nbytes);
  }

  void* recvBuffer = (char*)qOneRing + Nelements * Nbytes; // fix later

  int NsendMessages, NrecvMessages;

  ellipticOneRingExchangeStart(comm,
                               Nbytes,
                               NoneRingSendTotal,
                               NoneRingSend,
                               sendBuffer,
                               sendRequests,
                               &NsendMessages,
                               NoneRingRecvTotal,
                               NoneRingRecv,
                               recvBuffer,
                               recvRequests,
                               &NrecvMessages);

  // copy from q to qOneRing while data in transit
  memcpy(qOneRing, q, Nbytes * Nelements);

  ellipticOneRingExchangeFinish(comm,
                                NsendMessages, sendRequests,
                                NrecvMessages, recvRequests);
}

// occa memory version
void ellipticOneRingExchange(elliptic_t* elliptic,
                             elliptic_t* elliptic1,
                             size_t Nbytes, // message size per element
                             occa::memory &o_q,
                             occa::memory &o_qOneRing)
{
  // extract from original mesh
  mesh_t* mesh = elliptic->mesh;
  mesh_t* mesh1 = elliptic1->mesh;

  precon_t* precon = elliptic->precon;

  hlong NelementSend = precon->NoneRingSendTotal;
  hlong NelementRecv = precon->NoneRingRecvTotal;

  if(NelementSend + NelementRecv > 0) {
    mesh->device.finish();
    mesh1->device.finish();
    mesh->haloExtractKernel(NelementSend, mesh->Np,
                            precon->o_oneRingSendList,
                            o_q, precon->o_oneRingSendBuffer);

    precon->o_oneRingSendBuffer.copyTo(precon->oneRingSendBuffer, Nbytes * NelementSend, 0);

    dfloat* sendBuffer = (dfloat*) precon->oneRingSendBuffer;
    dfloat* recvBuffer = (dfloat*) precon->oneRingRecvBuffer;

    MPI_Request* sendRequests = precon->oneRingSendRequests;
    MPI_Request* recvRequests = precon->oneRingRecvRequests;

    // do exchange via MPI

    // count outgoing and incoming meshes
    int tag = 999;

    // initiate immediate send  and receives to each other process as needed
    int sendOffset = 0, recvOffset = 0, sendMessage = 0, recvMessage = 0;
    for(int r = 0; r < mesh->size; ++r)
      if(r != mesh->rank) {
        size_t recvCount = precon->NoneRingRecv[r] * Nbytes;
        if(recvCount) {
          MPI_Irecv(((char*)recvBuffer) + recvOffset,
                    recvCount, MPI_CHAR, r, tag, mesh->comm, recvRequests + recvMessage);

          recvOffset += recvCount;
          ++recvMessage;
        }

        size_t sendCount = precon->NoneRingSend[r] * Nbytes;
        if(sendCount) {
          MPI_Isend(((char*)sendBuffer) + sendOffset,
                    sendCount, MPI_CHAR, r, tag, mesh->comm, sendRequests + sendMessage);

          sendOffset += sendCount;
          ++sendMessage;
        }
      }

    MPI_Status* sendStatus = (MPI_Status*) calloc(sendMessage, sizeof(MPI_Status));
    MPI_Status* recvStatus = (MPI_Status*) calloc(recvMessage, sizeof(MPI_Status));

    MPI_Waitall(recvMessage, recvRequests, recvStatus);
    MPI_Waitall(sendMessage, sendRequests, sendStatus);

    free(recvStatus);
    free(sendStatus);

    mesh1->device.finish();

    // copy incoming to end of o_qOneRing
    o_qOneRing.copyFrom(recvBuffer, Nbytes * NelementRecv, Nbytes * mesh->Nelements); // offset into end of oneRing
  }

  mesh1->device.finish();

  // copy core
  o_qOneRing.copyFrom(o_q, Nbytes * mesh->Nelements, 0);
}

// build one ring including MPI exchange information

void ellipticBuildOneRing(elliptic_t* elliptic, dfloat lambda, occa::properties &kernelInfo)
{
  mesh_t* mesh = elliptic->mesh;

  //  occa::properties kernelInfo = ellipticKernelInfo(mesh);

  vertex_t* vertexSendList = (vertex_t*) calloc(mesh->Nelements * mesh->Nverts, sizeof(vertex_t));

  hlong* vertexSendCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  hlong* vertexRecvCounts = (hlong*) calloc(mesh->size, sizeof(hlong));

  hlong cnt = 0;
  for(hlong e = 0; e < mesh->Nelements; ++e)
    for(int v = 0; v < mesh->Nverts; ++v) {
      vertexSendList[cnt].vertex = mesh->EToV[e * mesh->Nverts + v];
      vertexSendList[cnt].element = e;
      vertexSendList[cnt].rank = mesh->rank;
      vertexSendList[cnt].rankN = mesh->rank;

      vertexSendList[cnt].sortTag = vertexSendList[cnt].vertex % mesh->size;
      ++vertexSendCounts[vertexSendList[cnt].sortTag];
      ++cnt;
    }

  // sort based on sortTag (=vertex%size)
  qsort(vertexSendList, cnt, sizeof(vertex_t), compareSortTag);

  // send sortTagCounts (hackety)
  MPI_Alltoall(vertexSendCounts, 1, MPI_HLONG,
               vertexRecvCounts, 1, MPI_HLONG,
               mesh->comm);

  // exchange vertices
  hlong* vertexSendDispls = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  hlong* vertexRecvDispls = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  hlong NvertexSend = 0;
  hlong NvertexRecv = 0;
  for(int r = 0; r < mesh->size; ++r) {
    NvertexSend += vertexSendCounts[r];
    NvertexRecv += vertexRecvCounts[r];

    vertexSendCounts[r] *= sizeof(vertex_t); // hack-hack-hack
    vertexRecvCounts[r] *= sizeof(vertex_t); // hack-hack-hack

    vertexSendDispls[r + 1] = vertexSendDispls[r] + vertexSendCounts[r];
    vertexRecvDispls[r + 1] = vertexRecvDispls[r] + vertexRecvCounts[r];
  }

  // hack-hack-hack
  vertex_t* vertexRecvList = (vertex_t*) calloc(NvertexRecv, sizeof(vertex_t));

  MPI_Alltoallv(vertexSendList, vertexSendCounts, vertexSendDispls, MPI_CHAR,
                vertexRecvList, vertexRecvCounts, vertexRecvDispls, MPI_CHAR,
                mesh->comm);

  for(int v = 0; v < NvertexRecv; ++v)
    vertexRecvList[v].sortTag = vertexRecvList[v].vertex;

  // sort received vertex based on sortTag (=vertex number)
  qsort(vertexRecvList, NvertexRecv, sizeof(vertex_t), compareSortTag);

  // count number of unique received vertices
  hlong NvertexUniqueRecv = (NvertexRecv > 0) ? 1:0;
  for(hlong n = 1; n < NvertexRecv; ++n)
    if(compareSortTag(vertexRecvList + n, vertexRecvList + n - 1) != 0) // new vertex
      ++NvertexUniqueRecv;

  // find offset of the start of each new unique vertex  in sorted list
  hlong* vertexUniqueRecvOffsets = (hlong*) calloc(NvertexUniqueRecv + 1, sizeof(hlong));

  cnt = 1;
  vertexUniqueRecvOffsets[0] = 0;
  for(hlong n = 1; n < NvertexRecv; ++n)
    if(compareSortTag(vertexRecvList + n, vertexRecvList + n - 1) != 0) { // new vertex
      vertexUniqueRecvOffsets[cnt] = n;
      ++cnt;
    }
  vertexUniqueRecvOffsets[cnt] = NvertexRecv; // cap at end

  // now count how many vertices to send to each rank
  hlong* vertexOneRingSendCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  hlong Ntotal = 0;
  for(hlong n = 0; n < NvertexUniqueRecv; ++n) {
    hlong start = vertexUniqueRecvOffsets[n];
    hlong end   = vertexUniqueRecvOffsets[n + 1];

    int NuniqueRecvMultiplicity = end - start;
    for(hlong m = start; m < end; ++m) {
      vertexOneRingSendCounts[vertexRecvList[m].rank]
        += NuniqueRecvMultiplicity; // watch out for this
      Ntotal += NuniqueRecvMultiplicity;
    }
  }

  vertex_t* vertexOneRingSendList = (vertex_t*) calloc(Ntotal, sizeof(vertex_t));
  cnt = 0;
  for(hlong n = 0; n < NvertexUniqueRecv; ++n) {
    hlong start = vertexUniqueRecvOffsets[n];
    hlong end   = vertexUniqueRecvOffsets[n + 1];

    for(hlong v1 = start; v1 < end; ++v1) // vertex v1 to be sent back with list of conns
      for(hlong v2 = start; v2 < end; ++v2) {
        vertexOneRingSendList[cnt] = vertexRecvList[v1];
        vertexOneRingSendList[cnt].rankN    = vertexRecvList[v2].rank;

        vertexOneRingSendList[cnt].sortTag  = vertexRecvList[v1].rank;
        ++cnt;
      }
  }

  hlong NvertexOneRingSend = cnt;

  // sort OneRing send list based on sort rank (=destination tag)
  qsort(vertexOneRingSendList, NvertexOneRingSend, sizeof(vertex_t), compareSortTag);   // check qsort counts

  // now figure out how many oneRing vertices to expect
  hlong* vertexOneRingRecvCounts = (hlong*) calloc(mesh->size, sizeof(hlong));
  MPI_Alltoall(vertexOneRingSendCounts, 1, MPI_HLONG,
               vertexOneRingRecvCounts, 1, MPI_HLONG,
               mesh->comm);

  // find displacements for
  hlong* vertexOneRingSendDispls = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  hlong* vertexOneRingRecvDispls = (hlong*) calloc(mesh->size + 1, sizeof(hlong));
  hlong NvertexOneRingRecv = 0;

  for(int r = 0; r < mesh->size; ++r) {
    NvertexOneRingRecv += vertexOneRingRecvCounts[r];
    vertexOneRingSendCounts[r] *= sizeof(vertex_t);
    vertexOneRingRecvCounts[r] *= sizeof(vertex_t);
    vertexOneRingSendDispls[r + 1] = vertexOneRingSendDispls[r] + vertexOneRingSendCounts[r];
    vertexOneRingRecvDispls[r + 1] = vertexOneRingRecvDispls[r] + vertexOneRingRecvCounts[r];
  }

  vertex_t* vertexOneRingRecvList =
    (vertex_t*) calloc(NvertexOneRingRecv, sizeof(vertex_t)); // hack-hack-hack

  // send element lists to the relevant ranks
  MPI_Alltoallv(vertexOneRingSendList, vertexOneRingSendCounts, vertexOneRingSendDispls, MPI_CHAR,
                vertexOneRingRecvList, vertexOneRingRecvCounts, vertexOneRingRecvDispls, MPI_CHAR,
                mesh->comm);

  // finally we now have a list of all elements that we need to send to form the 1-ring (to rule them all)
  vertex_t* vertexOneRingOut  = (vertex_t*) calloc(NvertexOneRingRecv, sizeof(vertex_t));
  memcpy(vertexOneRingOut,  vertexOneRingRecvList, NvertexOneRingRecv * sizeof(vertex_t));

  // sort the list by "neighbor rank then element"
  qsort(vertexOneRingOut,  NvertexOneRingRecv, sizeof(vertex_t), compareRankNElement);

  // remove elements connected to this rank from oneRing list
  cnt = 0;
  for(hlong v = 0; v < NvertexOneRingRecv; ++v)
    if(vertexOneRingOut[v].rankN != mesh->rank) // only connect connections with off rank elements
      vertexOneRingOut[cnt++] = vertexOneRingOut[v];

  hlong NvertexOneRingOut = cnt;

  // remove duplicate connections from oneRingInOut list
  if(NvertexOneRingOut) {
    cnt = 1; // assumes at least one oneRing element
    for(hlong v = 1; v < NvertexOneRingOut; ++v)
      if(!(vertexOneRingOut[v].element == vertexOneRingOut[cnt - 1].element
           && vertexOneRingOut[v].rank == vertexOneRingOut[cnt - 1].rank
           && vertexOneRingOut[v].rankN == vertexOneRingOut[cnt - 1].rankN
           ))
        vertexOneRingOut[cnt++] = vertexOneRingOut[v];
    NvertexOneRingOut = cnt;
  }

  printf("NvertexOneRingOut = %d, Nelements = %d\n", NvertexOneRingOut, mesh->Nelements);

  // next: put new stuff in elliptic
  //-1. count how many elements send to each rankN
  // 0. send count to each rankN
  // 1. populate NoneRingExchanges[0:size),
  // 4. adapt halo exchange to oneRingExchange
  // 5. oneRingExchange: globalNumbers for gs stuff
  // 3. set up the gs info  using exchange globalNumbers [ need to understand how to populate from the local elements on each rank to the oneRing ]
  // 6. oneRingExchange: geofacs (ggeo)
  // 7. build local continuous numbering and local global continuous numbering (see meshParallelConnectNodes)
  // 8. o_qOneRing
  // 9. how to precondition patch problem ?

  hlong NoneRingSendTotal = NvertexOneRingOut; // should rename things above
  hlong* oneRingSendList = (hlong*) calloc(NoneRingSendTotal + 1, sizeof(hlong));
  hlong* NoneRingSend = (hlong*) calloc(mesh->size, sizeof(hlong));
  hlong* NoneRingRecv = (hlong*) calloc(mesh->size, sizeof(hlong));

  for(hlong e = 0; e < NoneRingSendTotal; ++e) {
    vertex_t v = vertexOneRingOut[e];
    oneRingSendList[e] = v.element;
    ++NoneRingSend[v.rankN];
  }

  MPI_Alltoall(NoneRingSend, 1, MPI_HLONG,
               NoneRingRecv, 1, MPI_HLONG, mesh->comm);

  hlong NoneRingRecvTotal = 0;
  for(int r = 0; r < mesh->size; ++r)
    NoneRingRecvTotal += NoneRingRecv[r];

  int maxNbytes = mesh->Np * sizeof(dfloat); // fingers crossed.
  char* sendBuffer = (char*) calloc(maxNbytes * NoneRingSendTotal, sizeof(char));
  char* recvBuffer = (char*) calloc(maxNbytes * NoneRingRecvTotal, sizeof(char));

  MPI_Request* sendRequests = (MPI_Request*) calloc(mesh->size, sizeof(MPI_Request));
  MPI_Request* recvRequests = (MPI_Request*) calloc(mesh->size, sizeof(MPI_Request));

  //  mesh_t *mesh1 = new mesh_t[1];

  mesh_t* mesh1 = (mesh_t*) calloc(1, sizeof(mesh_t)); // check

  // single process communicator for mesh1
  MPI_Comm_split(mesh->comm, mesh->rank, mesh->rank, &mesh1->comm);

  MPI_Comm_rank(mesh1->comm, &mesh1->rank);
  MPI_Comm_size(mesh1->comm, &mesh1->size);

  mesh1->dim = mesh->dim;
  mesh1->Nverts = mesh->Nverts;
  mesh1->Nfaces = mesh->Nfaces;
  mesh1->NfaceVertices = mesh->NfaceVertices;
  mesh1->Nnodes = mesh->Nnodes;
  mesh1->Nfields = mesh->Nfields;

  mesh1->N   = mesh->N;

  mesh1->faceVertices =
    (int*) calloc(mesh1->NfaceVertices * mesh1->Nfaces, sizeof(int));

  memcpy(mesh1->faceVertices, mesh->faceVertices, mesh->NfaceVertices * mesh->Nfaces * sizeof(int));

  mesh1->Nelements = mesh->Nelements + NoneRingRecvTotal;

  mesh1->EX = (dfloat*) calloc(mesh1->Nelements * mesh1->Nverts, sizeof(dfloat));
  mesh1->EY = (dfloat*) calloc(mesh1->Nelements * mesh1->Nverts, sizeof(dfloat));
  mesh1->EZ = (dfloat*) calloc(mesh1->Nelements * mesh1->Nverts, sizeof(dfloat));
  ellipticOneRingExchange(mesh->comm,
                          mesh->Nelements, mesh->Nverts * sizeof(dfloat), mesh->EX,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->EX);

  ellipticOneRingExchange(mesh->comm,
                          mesh->Nelements, mesh->Nverts * sizeof(dfloat), mesh->EY,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->EY);

  ellipticOneRingExchange(mesh->comm,
                          mesh->Nelements, mesh->Nverts * sizeof(dfloat), mesh->EZ,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->EZ);

  mesh1->NboundaryFaces = mesh->NboundaryFaces;
  mesh1->boundaryInfo =
    (hlong*) calloc(mesh1->NboundaryFaces * (mesh1->NfaceVertices + 1), sizeof(hlong));
  memcpy(mesh1->boundaryInfo,
         mesh->boundaryInfo,
         mesh1->NboundaryFaces * (mesh1->NfaceVertices + 1) * sizeof(hlong));

  mesh1->EToV = (hlong*) calloc(mesh1->Nelements * mesh1->Nverts, sizeof(hlong));
  ellipticOneRingExchange(mesh->comm,
                          mesh->Nelements, mesh->Nverts * sizeof(hlong), mesh->EToV,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->EToV);

  meshParallelConnect(mesh1);

  meshConnectBoundary(mesh1);

  // correct bcs (replaces unconnected faces with Dirichlet)
  for(hlong e = 0; e < mesh1->Nelements; ++e)
    for(int f = 0; f < mesh1->Nfaces; ++f) {
      hlong id = e * mesh1->Nfaces + f;
      if(mesh1->EToE[id] == -1 &&
         mesh1->EToB[id] == -1) {
        mesh1->EToB[id] = 1; // hack to 1 assume Dirichlet
        mesh1->EToE[id] = e; // hack to 1 assume Dirichlet
      }
    }

  meshLoadReferenceNodesHex3D(mesh1, mesh1->N);

  mesh1->x = (dfloat*) calloc(mesh1->Nelements * mesh1->Np, sizeof(dfloat));
  mesh1->y = (dfloat*) calloc(mesh1->Nelements * mesh1->Np, sizeof(dfloat));
  mesh1->z = (dfloat*) calloc(mesh1->Nelements * mesh1->Np, sizeof(dfloat));

  ellipticOneRingExchange(mesh->comm, mesh->Nelements, mesh->Np * sizeof(dfloat), mesh->x,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->x);

  ellipticOneRingExchange(mesh->comm, mesh->Nelements, mesh->Np * sizeof(dfloat), mesh->y,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->y);

  ellipticOneRingExchange(mesh->comm,  mesh->Nelements, mesh->Np * sizeof(dfloat), mesh->z,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, mesh1->z);

  // this is all vanilla HEX --->
  meshGeometricFactorsHex3D(mesh1);

  meshHaloSetup(mesh1); // nada

  meshConnectFaceNodes3D(mesh1);

  meshSurfaceGeometricFactorsHex3D(mesh1);

  meshParallelConnectNodes(mesh1); // data
  // <------

  setupAide options1 = elliptic->options; // check this

  // manually specify preconditioner for oneRing grid
  //  options1.setArgs(string("PRECONDITIONER"),     string("MULTIGRID"));
  // options1.setArgs(string("PRECONDITIONER"),     string("JACOBI"));/
  //  options1.setArgs(string("MAXIMUM ITERATIONS"), string("1"));
  //  options1.setArgs(string("FIXED ITERATION COUNT"), string("1"));

  //options1.setArgs(string("PRECONDITIONER"), string("SEMFEM"));
  //  options1.setArgs(string("POLYNOMIAL DEGREE"),  string("1"));
  //options1.setArgs("PRECONDITIONER", "JACOBI");

  //  occa::properties kernelInfo1 = kernelInfo;

  options1.setArgs(string("KRYLOV SOLVER"),        string("PCG+FLEXIBLE"));
  options1.setArgs(string("PRECONDITIONER"),       string("MULTIGRID"));
  //  options1.setArgs(string("PRECONDITIONER"),       string("JACOBI"));
  options1.setArgs(string("MAXIMUM ITERATIONS"),   string("2"));
  options1.setArgs(string("FIXED ITERATION COUNT"),string("TRUE"));
  options1.setArgs(string("PARALMOND CYCLE"),string("KCYCLE"));

  mesh1->device = mesh->device; // check this
#if 1
  mesh1->defaultStream = mesh->defaultStream;
  mesh1->dataStream = mesh->dataStream;
  mesh1->computeStream = mesh->computeStream;
  mesh1->device.setStream(mesh->defaultStream);
#endif

  meshOccaPopulateDevice3D(mesh1, options1, kernelInfo);

  // set up
  elliptic_t* elliptic1 = ellipticSetup(mesh1, lambda, kernelInfo, options1);

  cout << "options1: " << elliptic1->options << endl;
  //  cout << "options: " << elliptic->options << endl;

  dfloat* ggeoNoJW = (dfloat*) calloc(mesh1->Np * mesh1->Nelements * 6,sizeof(dfloat));
  for(int e = 0; e < mesh1->Nelements; ++e)
    for(int n = 0; n < mesh1->Np; ++n) {
      ggeoNoJW[e * mesh1->Np * 6 + n + 0 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G00ID * mesh1->Np];
      ggeoNoJW[e * mesh1->Np * 6 + n + 1 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G01ID * mesh1->Np];
      ggeoNoJW[e * mesh1->Np * 6 + n + 2 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G02ID * mesh1->Np];
      ggeoNoJW[e * mesh1->Np * 6 + n + 3 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G11ID * mesh1->Np];
      ggeoNoJW[e * mesh1->Np * 6 + n + 4 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G12ID * mesh1->Np];
      ggeoNoJW[e * mesh1->Np * 6 + n + 5 *
               mesh1->Np] = mesh1->ggeo[e * mesh1->Np * mesh1->Nggeo + n + G22ID * mesh1->Np];
    }

  elliptic1->o_ggeoNoJW = mesh1->device.malloc(mesh1->Np * mesh1->Nelements * 6 * sizeof(dfloat),
                                               ggeoNoJW);

  // store the extended patch ring
  elliptic->precon->ellipticOneRing = (void*) elliptic1;

  elliptic->precon->NoneRingSendTotal = NoneRingSendTotal;
  elliptic->precon->NoneRingRecvTotal = NoneRingRecvTotal;

  elliptic->precon->oneRingSendList = oneRingSendList;

  elliptic->precon->NoneRingSend = NoneRingSend;
  elliptic->precon->NoneRingRecv = NoneRingRecv;

  elliptic->precon->oneRingSendBuffer = sendBuffer;
  elliptic->precon->oneRingRecvBuffer = recvBuffer;

  elliptic->precon->oneRingSendRequests = sendRequests;
  elliptic->precon->oneRingRecvRequests = recvRequests;

#if 0

  ellipticOneRingDiagnostics(elliptic, elliptic1, lambda);

#endif

  // build gs op to gather all contributions
  hlong* globalNums = (hlong*) calloc(mesh1->Nelements * mesh1->Np, sizeof(hlong));

  ellipticOneRingExchange(mesh->comm,
                          mesh->Nelements, mesh->Np * sizeof(hlong), mesh->globalIds,
                          NoneRingSendTotal, oneRingSendList, NoneRingSend, sendBuffer,
                          sendRequests,
                          NoneRingRecvTotal, NoneRingRecv, recvRequests, globalNums);

  elliptic->precon->oasOgs = ogsSetup(mesh1->Nelements * mesh1->Np,
                                      globalNums,
                                      mesh->comm,
                                      1,
                                      mesh->device);

  elliptic->precon->o_oneRingSendList =
    mesh->device.malloc(elliptic->precon->NoneRingSendTotal * sizeof(hlong),
                        elliptic->precon->oneRingSendList);

  elliptic->precon->o_oneRingSendBuffer =
    mesh->device.malloc(elliptic->precon->NoneRingSendTotal * mesh->Np * sizeof(dfloat));

  elliptic->precon->o_oneRingRecvBuffer =
    mesh->device.malloc(elliptic->precon->NoneRingRecvTotal * mesh->Np * sizeof(dfloat));

  free(vertexSendList);
  free(vertexSendCounts);
  free(vertexRecvCounts);
  free(vertexSendDispls);
  free(vertexRecvDispls);
  free(vertexRecvList);
  free(vertexUniqueRecvOffsets);
  free(vertexOneRingSendCounts);
  free(vertexOneRingSendList);
  free(vertexOneRingRecvCounts);
  free(vertexOneRingSendDispls);
  free(vertexOneRingRecvDispls);
  free(vertexOneRingRecvList);
}
