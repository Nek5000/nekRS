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
#include <unistd.h>

#include "omp.h"
#include "mesh.h"

occa::device occaDeviceConfig(setupAide &options, MPI_Comm comm)
{
  // OCCA build stuff
  char deviceConfig[BUFSIZ];
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int device_id = 0;

  if(options.compareArgs("DEVICE NUMBER", "LOCAL-RANK")) {
    long int hostId = gethostid();

    long int* hostIds = (long int*) calloc(size,sizeof(long int));
    MPI_Allgather(&hostId,1,MPI_LONG,hostIds,1,MPI_LONG,comm);

    int totalDevices = 0;
    for (int r = 0; r < rank; r++)
      if (hostIds[r] == hostId) device_id++;
    for (int r = 0; r < size; r++)
      if (hostIds[r] == hostId) totalDevices++;
  } else {
    options.getArgs("DEVICE NUMBER",device_id);
  }

  occa::properties deviceProps;

  if(options.compareArgs("THREAD MODEL", "CUDA")) {
    sprintf(deviceConfig, "mode: 'CUDA', device_id: %d", device_id);
  }else if(options.compareArgs("THREAD MODEL", "HIP"))  {
    sprintf(deviceConfig, "mode: 'HIP', device_id: %d",device_id);
  }else if(options.compareArgs("THREAD MODEL", "OPENCL"))  {
    int plat;
    options.getArgs("PLATFORM NUMBER", plat);
    sprintf(deviceConfig, "mode: 'OpenCL', device_id: %d, platform_id: %d", device_id, plat);
  }else if(options.compareArgs("THREAD MODEL", "OPENMP"))  {
    sprintf(deviceConfig, "mode: 'OpenMP' ");
  }else  {
    sprintf(deviceConfig, "mode: 'Serial', memory: { use_host_pointer: true }");
    sprintf(deviceConfig, "mode: 'Serial'");
    options.setArgs("THREAD MODEL", "SERIAL");
  }

  if(rank == 0) printf("Initializing device\n");
  occa::device device;
  device.setup((std::string)deviceConfig);

  if (device.mode() == "Serial")
    options.setArgs("THREAD MODEL", "SERIAL");

  if(rank == 0)
    std::cout << "active occa mode: " << device.mode() << "\n\n";

#ifdef USE_OCCA_MEM_BYTE_ALIGN
  occa::env::OCCA_MEM_BYTE_ALIGN = USE_OCCA_MEM_BYTE_ALIGN;
#endif

  int Nthreads = 1;
  omp_set_num_threads(Nthreads);
  //if(rank==0) printf("Number of OMP threads: %d\n", omp_get_num_threads());

  return device;
}


void* occaHostMallocPinned(occa::device &device,
                           size_t size,
                           void* source,
                           occa::memory &mem,
                           occa::memory &h_mem)
{
  occa::properties props;
  props["mapped"] = true;

  if(source != NULL)
    mem =  device.malloc(size, source);
  else
    mem =  device.malloc(size);

  h_mem =  device.malloc(size, props);

  void* ptr = h_mem.ptr(props);

  return ptr;
}
