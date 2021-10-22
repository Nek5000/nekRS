#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
#include <unistd.h>
#include "mpi.h"

#include "nrssys.hpp"
#include "setupAide.hpp"
#include "platform.hpp"
#include "configReader.hpp"

namespace {

occa::kernel fdmKernel;

occa::memory o_Sx;
occa::memory o_Sy;
occa::memory o_Sz;
occa::memory o_invL;
occa::memory o_u;
occa::memory o_Su;

int Np; 
int Nelements; 

double run(int Ntests)
{
  platform->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);
  const double start = MPI_Wtime();

  for(int test = 0; test < Ntests; ++test) {
    fdmKernel(Nelements, o_Su, o_Sx, o_Sy, o_Sz, o_invL, o_u);
  }

  platform->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);
  return (MPI_Wtime() - start) / Ntests;
} 

void* (*randAlloc)(int);

void* rand32Alloc(int N)
{
  float* v = (float*) malloc(N * sizeof(float));

  for(int n = 0; n < N; ++n)
    v[n] = drand48();

  return v;
}

void* rand64Alloc(int N)
{
  double* v = (double*) malloc(N * sizeof(double));

  for(int n = 0; n < N; ++n)
    v[n] = drand48();

  return v;
}

} // namespace

int main(int argc, char** argv)
{
  int rank = 0, size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  configRead(MPI_COMM_WORLD);
  std::string installDir(getenv("NEKRS_HOME"));
  setupAide options; 

  int err = 0;
  int cmdCheck = 0;

  int N;
  int okl = 1;
  int Ntests = -1;
  size_t wordSize = 8;

  while(1) {
    static struct option long_options[] =
    {
      {"p-order", required_argument, 0, 'p'},
      {"elements", required_argument, 0, 'e'},
      {"backend", required_argument, 0, 'b'},
      {"arch", required_argument, 0, 'a'},
      {"fp32", no_argument, 0, 'f'},
      {"help", required_argument, 0, 'h'},
      {"iterations", required_argument, 0, 'i'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    int c = getopt_long (argc, argv, "", long_options, &option_index);

    if (c == -1)
      break;

    switch(c) {
    case 'p':
      N = atoi(optarg); 
      cmdCheck++; 
      break;
    case 'e':
      Nelements = atoi(optarg);
      cmdCheck++;
      break;
    case 'b':
      options.setArgs("THREAD MODEL", std::string(optarg));
      cmdCheck++;
      break;
    case 'f':
      wordSize = 4;;
      break;
    case 'i':
      Ntests = atoi(optarg);
      break;
    case 'h':
      err = 1;
      break;
    default:
      err = 1;
    }
  }

  if(err || cmdCheck != 3) {
    if(rank == 0)
      printf("Usage: ./nekrs-fdm  --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>\n"
             "                    [--fp32] [--iterations <n>]\n"); 
    exit(1); 
  }

  Nelements = std::max(1, Nelements/size);
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;

  platform = platform_t::getInstance(options, MPI_COMM_WORLD, MPI_COMM_WORLD); 
  const int Nthreads =  omp_get_max_threads();

  // build+load kernel
  occa::properties props = platform->kernelInfo + meshKernelProperties(N);
  if(wordSize == 4) props["defines/pfloat"] = "float";
  else props["defines/pfloat"] = "dfloat";

  props["defines/p_Nq_e"] = Nq;
  props["defines/p_Np_e"] = Np;
  props["defines/p_overlap"] = 0;

  // always benchmark ASM
  props["defines/p_restrict"] = 0;

  std::string kernelName = "fusedFDM";
  const std::string ext = (platform->device.mode() == "Serial") ? ".c" : ".okl";
  const std::string fileName = 
    installDir + "/okl/elliptic/ellipticSchwarzSolverHex3D" + ext;

  fdmKernel = platform->device.buildKernel(fileName, props, true);

  // populate arrays
  randAlloc = &rand64Alloc; 
  if(wordSize == 4) randAlloc = &rand32Alloc;

  void *Sx   = randAlloc(Nelements * Nq * Nq);
  void *Sy   = randAlloc(Nelements * Nq * Nq);
  void *Sz   = randAlloc(Nelements * Nq * Nq);
  void *invL = randAlloc(Nelements * Np);
  void *Su   = randAlloc(Nelements * Np);
  void *u    = randAlloc(Nelements * Np);

  o_Sx = platform->device.malloc(Nelements * Nq * Nq * wordSize, Sx);
  free(Sx);
  o_Sy = platform->device.malloc(Nelements * Nq * Nq * wordSize, Sy);
  free(Sy);
  o_Sz = platform->device.malloc(Nelements * Nq * Nq * wordSize, Sz);
  free(Sz);
  o_invL = platform->device.malloc(Nelements * Np * wordSize, invL);
  free(invL);
  o_Su = platform->device.malloc(Nelements * Np * wordSize, Su);
  free(Su);
  o_u = platform->device.malloc(Nelements * Np * wordSize, u);
  free(u);

  // warm-up
  double elapsed = run(10);
  const int elapsedTarget = 10;
  if(Ntests < 0) Ntests = elapsedTarget/elapsed;

  // ***** 
  elapsed = run(Ntests);
  // ***** 
 
  // print statistics
  const dfloat GDOFPerSecond = (size * Nelements * (N* N * N) / elapsed) / 1.e9;

  size_t bytesPerElem = (3 * Np + 3 * Nq * Nq) * wordSize;
  const double bw = (size * Nelements * bytesPerElem / elapsed) / 1.e9;

  double flopsPerElem = 12 * Nq * Np + Np;
  const double gflops = (size * flopsPerElem * Nelements / elapsed) / 1.e9;

  if(rank == 0)
    std::cout << "MPItasks=" << size
              << " OMPthreads=" << Nthreads
              << " NRepetitions=" << Ntests
              << " N=" << N
              << " Nelements=" << size * Nelements
              << " elapsed time=" << elapsed
              << " wordSize=" << 8*wordSize
              << " GDOF/s=" << GDOFPerSecond
              << " GB/s=" << bw
              << " GFLOPS/s=" << gflops
              << "\n";

  MPI_Finalize();
  exit(0);
}
