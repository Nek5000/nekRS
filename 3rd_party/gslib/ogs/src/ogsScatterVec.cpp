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

#include "ogstypes.h"
#include "ogs.hpp"
#include "ogsKernels.hpp"
#include "ogsInterface.h"

#include "scatterVec.tpp"

void ogsScatterVec_op(void *gv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs);

void ogsScatterVec(occa::memory o_sv, 
               occa::memory o_v, 
               const int k,
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  ogsScatterVecStart (o_sv, o_v, k, type, op, ogs);
  ogsScatterVecFinish(o_sv, o_v, k, type, op, ogs);
}

void ogsScatterVecStart(occa::memory o_sv, 
                    occa::memory o_v, 
                    const int k,
                    const char *type, 
                    const char *op, 
                    ogs_t *ogs){
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    if (ogs::o_haloBuf.size() < ogs->NhaloGather*Nbytes*k) {
      if (ogs::o_haloBuf.size()) ogs::o_haloBuf.free();
      //      ogs::o_haloBuf = ogs->device.mappedAlloc(ogs->NhaloGather*Nbytes*k);
      //      ogs::haloBuf = ogs::o_haloBuf.getMappedPointer();
      ogs::haloBuf = ogsHostMallocPinned(ogs->device, ogs->NhaloGather*Nbytes*k, NULL, ogs::o_haloBuf, ogs::h_haloBuf);
    }
  }

  if(ogs->NlocalGather) {
    occaScatterVec(ogs->NlocalGather, k, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v, o_sv);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);

    if (ogs->NownedHalo)
      o_v.copyTo(ogs::haloBuf, ogs->NownedHalo*Nbytes*k, 
                              ogs->NlocalGather*Nbytes*k, "async: true");

    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsScatterVecFinish(occa::memory o_sv, 
                     occa::memory o_v, 
                     const int k,
                     const char *type, 
                     const char *op, 
                     ogs_t *ogs){
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    // MPI based scatter using gslib
    ogsHostScatterVec(ogs::haloBuf, k, type, op, ogs->haloGshNonSym);

    // copy totally scattered halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");

    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);

    occaScatterVec(ogs->NhaloGather, k, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_sv);
  }
}

void ogsScatterVec(void *sv, 
               void *v, 
               const int k,
               const char *type, 
               const char *op, 
               ogs_t *ogs){
  
  size_t Nbytes;
  if (!strcmp(type, "float")) 
    Nbytes = sizeof(float);
  else if (!strcmp(type, "double")) 
    Nbytes = sizeof(double);
  else if (!strcmp(type, "int")) 
    Nbytes = sizeof(int);
  else if (!strcmp(type, "long long int")) 
    Nbytes = sizeof(long long int);

  if (ogs->NhaloGather) {
    if (ogs::hostBufSize < ogs->NhaloGather*Nbytes*k) {
      if (ogs::hostBufSize) free(ogs::hostBuf);
      ogs::hostBuf = (void *) malloc(ogs->NhaloGather*Nbytes*k);
    }
  }

  ogsScatterVec_op(sv, v, k, Nbytes, type, ogs);
}

void ogsScatterVec_op(void *sv, void *v, const int k, const size_t Nbytes, const char *type, ogs_t *ogs){

  if (!strcmp(type, "float")) 
    scatterVec<float>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (float*)v, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatterVec<double>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (double*)v, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatterVec<int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (int*)v, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatterVec<long long int>(ogs->NlocalGather, k, ogs->localGatherOffsets,
                      ogs->localGatherIds, (long long int*)v, (long long int*)sv);

  if (ogs->NhaloGather) {
    if (ogs->NownedHalo)
      memcpy(ogs::hostBuf, (char*) v+ogs->NlocalGather*Nbytes*k, ogs->NownedHalo*Nbytes*k);

    // MPI based scatterVec using gslib
    ogsHostScatterVec(ogs::hostBuf, k, type, ogsAdd, ogs->haloGshNonSym);
  }

  if (!strcmp(type, "float")) 
    scatterVec<float>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (float*)ogs::hostBuf, (float*)sv);
  else if (!strcmp(type, "double")) 
    scatterVec<double>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (double*)ogs::hostBuf, (double*)sv);
  else if (!strcmp(type, "int")) 
    scatterVec<int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (int*)ogs::hostBuf, (int*)sv);
  else if (!strcmp(type, "long long int")) 
    scatterVec<long long int>(ogs->NhaloGather, k, ogs->haloGatherOffsets,
                      ogs->haloGatherIds, (long long int*)ogs::hostBuf, (long long int*)sv);
}

void occaScatterVec(const  dlong Nscatter,
                const int Nentries,
                occa::memory o_scatterStarts,
                occa::memory o_scatterIds,
                const char* type,
                const char* op,
                occa::memory  o_v,
                occa::memory  o_sv) {
  
  if      (!strcmp(type, "float")) 
    ogs::scatterVecKernel_float(Nscatter, Nentries, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "double")) 
    ogs::scatterVecKernel_double(Nscatter, Nentries, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "int")) 
    ogs::scatterVecKernel_int(Nscatter, Nentries, o_scatterStarts, o_scatterIds, o_v, o_sv);
  else if (!strcmp(type, "long long int")) 
    ogs::scatterVecKernel_long(Nscatter, Nentries, o_scatterStarts, o_scatterIds, o_v, o_sv);
}
