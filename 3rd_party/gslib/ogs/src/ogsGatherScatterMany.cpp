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

void ogsGatherScatterMany(void *v, 
                      const int k,
                      const dlong stride,
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

  void* V[k];
  for (int i=0;i<k;i++) V[i] = (char*)v + i*stride*Nbytes;

  ogsHostGatherScatterMany(V, k, type, op, ogs->hostGsh);
}

void ogsGatherScatterMany(occa::memory o_v, 
                      const int k,
                      const dlong stride,
                      const char *type, 
                      const char *op, 
                      ogs_t *ogs){
  ogsGatherScatterManyStart (o_v, k, stride, type, op, ogs);
  ogsGatherScatterManyFinish(o_v, k, stride, type, op, ogs);
}

void ogsGatherScatterManyStart(occa::memory o_v, 
                          const int k,
                          const dlong stride,
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

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    occaGatherMany(ogs->NhaloGather, k, stride, ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherScatterManyFinish(occa::memory o_v, 
                          const int k,
                          const dlong stride,
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

  if(ogs->NlocalGather) {
    occaGatherScatterMany(ogs->NlocalGather, k, stride, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    void* H[k];
    for (int i=0;i<k;i++) H[i] = (char*)ogs::haloBuf + i*ogs->NhaloGather*Nbytes;

    // MPI based gather scatter using libgs
    ogsHostGatherScatterMany(H, k, type, op, ogs->haloGshSym);

    // copy totally gather halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");

    // do scatter back to local nodes
    occaScatterMany(ogs->NhaloGather, k, ogs->NhaloGather, stride, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_v);
    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}

void occaGatherScatterMany(const  dlong Ngather,
                const int Nentries,
                const dlong stride,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v) {
  
  if      ((!strcmp(type, "float"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterManyKernel_floatAdd(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterManyKernel_floatMul(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterManyKernel_floatMin(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterManyKernel_floatMax(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterManyKernel_doubleAdd(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterManyKernel_doubleMul(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterManyKernel_doubleMin(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterManyKernel_doubleMax(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterManyKernel_intAdd(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterManyKernel_intMul(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterManyKernel_intMin(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterManyKernel_intMax(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterManyKernel_longAdd(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterManyKernel_longMul(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterManyKernel_longMin(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterManyKernel_longMax(Ngather, Nentries, stride, o_gatherStarts, o_gatherIds, o_v);
}
