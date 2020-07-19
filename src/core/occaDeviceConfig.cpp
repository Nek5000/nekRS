#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "omp.h"
#include "mpi.h"
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
//    sprintf(deviceConfig, "mode: 'Serial', memory: { use_host_pointer: true }");
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

  occa::initTimer(device);

  return device;
}
