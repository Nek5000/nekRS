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

void ogsGatherScatterVec(void *v, 
                      const int k,
                      const char *type, 
                      const char *op, 
                      ogs_t *ogs){
  ogsHostGatherScatterVec(v, k, type, op, ogs->hostGsh);
}

void ogsGatherScatterVec(occa::memory o_v, 
                      const int k,
                      const char *type, 
                      const char *op, 
                      ogs_t *ogs){
  ogsGatherScatterVecStart (o_v, k, type, op, ogs);
  ogsGatherScatterVecFinish(o_v, k, type, op, ogs);
}

void ogsGatherScatterVecStart(occa::memory o_v, 
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

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    occaGatherVec(ogs->NhaloGather, k, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, o_v, ogs::o_haloBuf);
    
    ogs->device.finish();
    ogs->device.setStream(ogs::dataStream);
    ogs::o_haloBuf.copyTo(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");
    ogs->device.setStream(ogs::defaultStream);
  }
}


void ogsGatherScatterVecFinish(occa::memory o_v, 
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

  if(ogs->NlocalGather) {
    occaGatherScatterVec(ogs->NlocalGather, k, ogs->o_localGatherOffsets, ogs->o_localGatherIds, type, op, o_v);
  }

  if (ogs->NhaloGather) {
    ogs->device.setStream(ogs::dataStream);
    ogs->device.finish();

    // MPI based gather scatter using libgs
    ogsHostGatherScatterVec(ogs::haloBuf, k, type, op, ogs->haloGshSym);

    // copy totally gather halo data back from HOST to DEVICE
    ogs::o_haloBuf.copyFrom(ogs::haloBuf, ogs->NhaloGather*Nbytes*k, 0, "async: true");

    // do scatter back to local nodes
    occaScatterVec(ogs->NhaloGather, k, ogs->o_haloGatherOffsets, ogs->o_haloGatherIds, type, op, ogs::o_haloBuf, o_v);
    ogs->device.finish();
    ogs->device.setStream(ogs::defaultStream);
  }
}

void occaGatherScatterVec(const  dlong Ngather,
                const int Nentries,
                occa::memory o_gatherStarts,
                occa::memory o_gatherIds,
                const char* type,
                const char* op,
                occa::memory  o_v) {
  
  if      ((!strcmp(type, "float"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterVecKernel_floatAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterVecKernel_floatMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterVecKernel_floatMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "float"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterVecKernel_floatMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterVecKernel_doubleAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterVecKernel_doubleMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterVecKernel_doubleMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "double"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterVecKernel_doubleMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterVecKernel_intAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterVecKernel_intMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterVecKernel_intMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "int"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterVecKernel_intMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "add"))) 
    ogs::gatherScatterVecKernel_longAdd(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "mul"))) 
    ogs::gatherScatterVecKernel_longMul(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "min"))) 
    ogs::gatherScatterVecKernel_longMin(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
  else if ((!strcmp(type, "long long int"))&&(!strcmp(op, "max"))) 
    ogs::gatherScatterVecKernel_longMax(Ngather, Nentries, o_gatherStarts, o_gatherIds, o_v);
}
