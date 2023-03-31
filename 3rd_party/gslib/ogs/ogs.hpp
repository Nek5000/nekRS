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

/*
  OCCA Gather/Scatter Library

  The code
  
  	int N;
    long long int id[N];    // the long long int and int types are defined in "types.h"
    int   haloFlag[N];    
    ...
    struct ogs_t *ogs = ogsSetup(N, id, &comm, verbose);
    
  defines a partition of the set of (processor, local index) pairs,
    (p,i) \in S_j  iff   abs(id[i]) == j  on processor p
  That is, all (p,i) pairs are grouped together (in group S_j) that have the
    same id (=j).
  S_0 is treated specially --- it is ignored completely
    (i.e., when id[i] == 0, local index i does not participate in any
    gather/scatter operation
  
  When a partition of pairs, is shared between MPI processes, one specific 
  member of the parition is chosen in an arbitrary but consistent way. This member
  is considered the 'owning' member and is used to define the nonsymmetric gather 
  and scatter operations.

  When "ogs" is no longer needed, free it with
  
    ogsFree(ogs);
  
  A basic gatherScatter operation is, e.g.,
  
    occa::memory o_v; 
    ...
    ogsGatherScatter(o_v, ogsDouble, ogsAdd, ogs);
    
  This gs call has the effect,
  
    o_v[i] <--  \sum_{ (p,j) \in S_{id[i]} } o_v_(p) [j]
    
  where o_v_(p) [j] means o_v[j] on proc p. In other words, every o_v[i] is replaced
  by the sum of all o_v[j]'s with the same id, given by id[i]. This accomplishes
  "direct stiffness summation" corresponding to the action of QQ^T, where
  "Q" is a boolean matrix that copies from a global vector (indexed by id)
  to the local vectors indexed by (p,i) pairs.
  
  Summation on doubles is not the only operation and datatype supported. Support 
  includes the operations
    ogsAdd, ogsMul, ogsMax, ogsMin
  and datatypes
    ogsDfloat, ogsDouble, ogsFloat, ogsInt, ogsLong, ogsDlong, ogsHlong.
   
  For the nonsymmetric behavior, the operations are
  
    ogsGather (o_Gv, o_v, gsDouble, gsAdd, ogs);
    ogsScatter(o_Sv, o_v, gsDouble, gsAdd, ogs);

  A version for vectors (contiguously packed) is, e.g.,
  
    occa::memory o_v[k];
    ogsGatherScatterVec(o_v,k, ogsDouble,ogsAdd, transpose, ogs);
  
  which is like "gs" operating on the datatype double[k],
  with summation here being vector summation. Number of messages sent
  is independent of k.
  
  For combining the communication for "gs" on multiple arrays:
      
    occa::memory o_v1, o_v2, ..., o_vk;
    
    ogsGatherScatterMany(o_v, k, stride, ogsDouble, op, ogs);
  
  when the arrays o_v1, o_v2, ..., o_vk are packed in o_v as 

    o_v1 = o_v + 0*stride, o_v2 = o_v + 1*stride, ...

  This call is equivalent to
  
    ogsGatherScatter(o_v1, gsDouble, op, ogs);
    ogsGatherScatter(o_v2, gsDouble, op, ogs);
    ...
    ogsGatherScatter(o_vk, gsDouble, op, ogs);
    
  except that all communication is done together.

*/  

#ifndef OGS_HPP
#define OGS_HPP 1

#include <functional>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <occa.hpp>
#include <string>

#include "mpi.h"

#define ogsFloat  "float"
#define ogsDouble "double"
#define ogsDfloat dfloatString
#define ogsPfloat pfloatString
#define ogsInt  "int"
#define ogsLong "long long int"
#define ogsDlong intString
#define ogsHlong hlongString

#define ogsAdd "add"
#define ogsMul "mul"
#define ogsMax "max"
#define ogsMin "min"

//#define OGS_ENABLE_TIMER

typedef std::function<occa::kernel(const std::string &fileName,
                                   const std::string &kernelName,
                                   const occa::properties &props)> ogsBuildKernel_t;

// OCCA+gslib gather scatter
typedef struct {
  
  MPI_Comm comm;
  occa::device device;

  int         N;
  int         Ngather;        //  total number of gather nodes
  int         Nlocal;         //  number of local nodes
  int         NlocalGather;   //  number of local gathered nodes 
  int         Nhalo;          //  number of halo nodes
  int         NhaloGather;    //  number of gathered nodes on halo
  hlong       NhaloGatherGlobal;
  int         NownedHalo;     //  number of owned halo nodes
  int         NrowBlocks;

  int         *localGatherOffsets;
  int         *localGatherIds;
  occa::memory o_localGatherOffsets;  
  occa::memory o_localGatherIds;      
  occa::memory o_blockRowStarts;

  int         *haloGatherOffsets;
  int         *haloGatherIds;
  occa::memory o_haloGatherOffsets;
  occa::memory o_haloGatherIds;    

  void         *hostGsh;          // gslib gather 
  void         *haloGshSym;       // gslib gather 
  void         *haloGshNonSym;    // gslib gather 
  
  //degree vectors
  dfloat *invDegree, *gatherInvDegree;
  occa::memory o_invDegree;
  occa::memory o_gatherInvDegree;

}ogs_t;


ogs_t *ogsSetup(int N, long long int *ids, MPI_Comm &comm, 
                int verbose, occa::device device);

void ogsFree(ogs_t* ogs);

// Host array versions
void ogsGatherScatter    (void  *v, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call
void ogsGatherScatterVec (void  *v, const int k, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call
void ogsGatherScatterMany(void  *v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call

void ogsGather    (void  *gv, void  *v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherVec (void  *gv, void  *v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherMany(void  *gv, void  *v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);

void ogsScatter    (void  *sv, void  *v, const char *type, const char *op, ogs_t *ogs);
void ogsScatterVec (void  *sv, void  *v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsScatterMany(void  *sv, void  *v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);

void ogsHostTic(MPI_Comm comm, int ifSync);
void ogsHostToc();
double ogsTime(bool reportHostTime);
void ogsResetTime();

// Synchronous device buffer versions
void ogsGatherScatter    (occa::memory  o_v, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call
void ogsGatherScatterVec (occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call
void ogsGatherScatterMany(occa::memory  o_v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs); //wrapper for gslib call

void ogsGather    (occa::memory  o_gv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherVec (occa::memory  o_gv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherMany(occa::memory  o_gv, occa::memory  o_v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);

void ogsScatter    (occa::memory  o_sv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsScatterVec (occa::memory  o_sv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsScatterMany(occa::memory  o_sv, occa::memory  o_v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);

// Asynchronous device buffer versions
void ogsGatherScatterStart     (occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherScatterFinish    (occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherScatterVecStart  (occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherScatterVecFinish (occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherScatterManyStart (occa::memory  o_v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);
void ogsGatherScatterManyFinish(occa::memory  o_v, const int k, const int stride, const char *type, const char *op, ogs_t *ogs);

void ogsGatherStart     (occa::memory  o_Gv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherFinish    (occa::memory  o_Gv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsGatherVecStart  (occa::memory  o_Gv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherVecFinish (occa::memory  o_Gv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsGatherManyStart (occa::memory  o_Gv, occa::memory  o_v, const int k, const int gstride, const int stride, const char *type, const char *op, ogs_t *ogs);
void ogsGatherManyFinish(occa::memory  o_Gv, occa::memory  o_v, const int k, const int gstride, const int stride, const char *type, const char *op, ogs_t *ogs);

void ogsScatterStart     (occa::memory  o_Sv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsScatterFinish    (occa::memory  o_Sv, occa::memory  o_v, const char *type, const char *op, ogs_t *ogs);
void ogsScatterVecStart  (occa::memory  o_Sv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsScatterVecFinish (occa::memory  o_Sv, occa::memory  o_v, const int k, const char *type, const char *op, ogs_t *ogs);
void ogsScatterManyStart (occa::memory  o_Sv, occa::memory  o_v, const int k, const int sstride, const int stride, const char *type, const char *op, ogs_t *ogs);
void ogsScatterManyFinish(occa::memory  o_Sv, occa::memory  o_v, const int k, const int sstride, const int stride, const char *type, const char *op, ogs_t *ogs);

void *ogsHostMallocPinned(occa::device &device, size_t size, void *source, occa::memory &mem, occa::memory &h_mem);

#define USE_OOGS

enum oogs_mode { OOGS_LOCAL, OOGS_DEFAULT, OOGS_HOSTMPI, OOGS_DEVICEMPI, OOGS_AUTO };
enum oogs_modeExchange { OOGS_EX_PW, OOGS_EX_NBC };

typedef struct {
  int* sendcounts;
  int* senddispls;
  int* recvcounts;
  int* recvdispls;
  MPI_Comm comm;
} nbc_t;

typedef struct {

  ogs_t *ogs;

  MPI_Comm comm;
  int rank;

  occa::memory h_buffSend, h_buffRecv;
  unsigned char *bufSend, *bufRecv;

  occa::memory o_bufSend, o_bufRecv;

  occa::memory o_scatterOffsets, o_gatherOffsets;
  occa::memory o_scatterIds, o_gatherIds;

  int earlyPrepostRecv;

  oogs_mode mode;
  oogs_modeExchange modeExchange;

  nbc_t nbc;

} oogs_t;

namespace oogs{

void start(occa::memory &o_v, const int k, const int stride, const char *type, const char *op, oogs_t *h);
void finish(occa::memory &o_v, const int k, const int stride, const char *type, const char *op, oogs_t *h);
void startFinish(void *v, const int k, const int stride, const char *type, const char *op, oogs_t *h);
void startFinish(occa::memory &o_v, const int k, const int stride, const char *type, const char *op, oogs_t *h);
void compile(const occa::device& device, ogsBuildKernel_t buildKernel, std::string mode, MPI_Comm comm, bool verbose = false);
oogs_t *setup(ogs_t *ogs, int nVec, int stride, const char *type, std::function<void()> callback, oogs_mode gsMode);
oogs_t *setup(int N, long long int *ids, const int k, const int stride, const char *type, MPI_Comm &comm,
              int verbose, occa::device device, std::function<void()> callback, oogs_mode mode);
void gpu_mpi(int val);
void overlap(int val);
int gpu_mpi();
void destroy(oogs_t *h);

}

#endif
