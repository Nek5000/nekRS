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
#include <stddef.h>

#include "mesh.h"

// ok to use o_v = o_gsv
void meshParallelGatherScatter(mesh_t *mesh,
                               ogs_t *ogs, 
                               occa::memory &o_v){
  int one = 1;
  dlong dOne = 1;

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, one, dOne, o_v, ogs->o_haloGatherTmp);
    ogs->o_haloGatherTmp.copyTo(ogs->haloGatherTmp);
  }

  if(ogs->NnonHaloGather) {
    mesh->gatherScatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, one, dOne, o_v);
  }

  if (ogs->NhaloGather) {
    // MPI based gather scatter using libgs
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

    // copy totally gather halo data back from HOST to DEVICE
    mesh->device.setStream(mesh->dataStream);
    ogs->o_haloGatherTmp.copyFrom(ogs->haloGatherTmp,"async: true");

    // wait for async copy
    mesh->device.finish();

    // do scatter back to local nodes
    mesh->scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, one, dOne, ogs->o_haloGatherTmp, o_v);
    mesh->device.finish();
    mesh->device.setStream(mesh->defaultStream);
  }
}

void meshParallelGather(mesh_t *mesh,
                        ogs_t *ogs, 
                        occa::memory &o_v,
                        occa::memory &o_gv){

  int one = 1;
  dlong dOne = 1;

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    mesh->gatherKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, one, dOne, o_v, ogs->o_haloGatherTmp);
    ogs->o_haloGatherTmp.copyTo(ogs->haloGatherTmp);
  }

  if(ogs->NnonHaloGather) {
    mesh->gatherKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, one, dOne, o_v, o_gv);
  }

  if (ogs->NhaloGather) {
    // MPI based gather scatter using libgs
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

    mesh->device.setStream(mesh->dataStream);
    
    // copy gathered halo data back from HOST to DEVICE
    ogs->o_haloGatherTmp.copyFrom(ogs->haloGatherTmp,"async: true");
    
    // wait for async copy
    mesh->device.finish();

    // insert totally gathered halo node data - need this kernel 
    mesh->putKernel(ogs->NhaloGather, ogs->o_haloGatherTmp, ogs->o_ownedHaloGatherIds, o_gv); 

    mesh->device.finish();
    mesh->device.setStream(mesh->defaultStream);
  }
}

void meshParallelScatter(mesh_t *mesh,
                        ogs_t *ogs, 
                        occa::memory &o_v,
                        occa::memory &o_sv){

  int one = 1;
  dlong dOne = 1;

  // gather halo nodes on device
  if (ogs->NhaloGather) {
    mesh->getKernel(ogs->NhaloGather, o_v, ogs->o_ownedHaloGatherIds, ogs->o_haloGatherTmp);
    ogs->o_haloGatherTmp.copyTo(ogs->haloGatherTmp);
  }

  if(ogs->NnonHaloGather) {
    mesh->scatterKernel(ogs->NnonHaloGather, ogs->o_nonHaloGatherOffsets, ogs->o_nonHaloGatherLocalIds, one, dOne, o_v, o_sv);
  }

  if (ogs->NhaloGather) {
    // MPI based gather scatter using libgs
    gsParallelGatherScatter(ogs->haloGsh, ogs->haloGatherTmp, dfloatString, "add");

    // copy totally gather halo data back from HOST to DEVICE
    mesh->device.setStream(mesh->dataStream);
    
    // copy owned and gathered halo data back from HOST to DEVICE
    ogs->o_haloGatherTmp.copyFrom(ogs->haloGatherTmp,"async: true");
    
    // wait for async copy
    mesh->device.finish();

    // insert totally gathered halo node data - need this kernel 
    mesh->scatterKernel(ogs->NhaloGather, ogs->o_haloGatherOffsets, ogs->o_haloGatherLocalIds, one, dOne, ogs->o_haloGatherTmp, o_sv);

    mesh->device.finish();
    mesh->device.setStream(mesh->defaultStream);
  }
}