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

occa::kernel subcyclingKernel;
occa::memory o_elementList;
occa::memory o_cubD;
occa::memory o_cubInterpT;
occa::memory o_BdivW;
occa::memory o_conv;
occa::memory o_Ud;
occa::memory o_NU;
occa::memory o_invLMM;

int Nelements;
int Np; 
int cubNp;
dlong fieldOffset;
dlong cubatureOffset;
bool dealias;

double run(int Ntests)
{
  platform->device.finish();
  MPI_Barrier(MPI_COMM_WORLD);
  const double start = MPI_Wtime();

  const dfloat c0 = 0.1;
  const dfloat c1 = 0.2;
  const dfloat c2 = 0.3;

  for(int test = 0; test < Ntests; ++test) {
    if(!dealias) {
      subcyclingKernel(Nelements, o_elementList, o_cubD, fieldOffset,
        0, o_invLMM, o_BdivW, c0, c1, c2, o_conv, o_Ud, o_NU);
    } else {
      subcyclingKernel(Nelements, o_elementList, o_cubD, o_cubInterpT, fieldOffset,
        cubatureOffset, 0, o_invLMM, o_BdivW, c0, c1, c2, o_conv, o_Ud, o_NU);
    }
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

  std::string threadModel;
  int N;
  int cubN = -1;
  int okl = 1;
  int Ntests = -1;
  int nEXT = 2;
  dealias = true;
  static constexpr size_t wordSize = 8;

  while(1) {
    static struct option long_options[] =
    {
      {"p-order", required_argument, 0, 'p'},
      {"ext-order", required_argument, 0, 'x'},
      {"c-order", required_argument, 0, 'c'},
      {"no-cubature", no_argument, 0, 'd'},
      {"elements", required_argument, 0, 'e'},
      {"backend", required_argument, 0, 'b'},
      {"arch", required_argument, 0, 'a'},
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
    case 'c':
      cubN = atoi(optarg); 
      break;
    case 'd':
      dealias = false;
      break;
    case 'x':
      nEXT = atoi(optarg); 
      if(nEXT <= 0 || nEXT > 3){
        if(rank == 0){
          printf("Error, 0 < nEXT <= 3!\n");
        }
        exit(1);
      }
      break;
    case 'e':
      Nelements = atoi(optarg);
      cmdCheck++;
      break;
    case 'b':
      options.setArgs("THREAD MODEL", std::string(optarg));
      cmdCheck++;
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
      printf("Usage: ./nekrs-bench-advsub  --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>\n"
             "                    [--c-order <n>] [--no-cubature] [--ext-order <n>] [--iterations <n>]\n"); 
    exit(1); 
  }

  if(cubN < 0) {
    if(dealias) cubN = round((3./2) * (N+1) - 1) - 1;
    else cubN = N;
  }
  if(cubN < N){
    if(rank == 0)
      printf("Error: cubature order (%d) must be larger than or equal to the quadrature order (%d)!\n",
        cubN,
        N);
    exit(1);
  }
  Nelements = std::max(1, Nelements/size);
  const int Nq = N + 1;
  Np = Nq * Nq * Nq;
  const int cubNq = cubN + 1;
  cubNp = cubNq * cubNq * cubNq;
  fieldOffset = Np * Nelements;
  const int pageW = ALIGN_SIZE / sizeof(dfloat);
  if (fieldOffset % pageW) fieldOffset = (fieldOffset / pageW + 1) * pageW;
  cubatureOffset = std::max(fieldOffset, Nelements * cubNp);

  platform = platform_t::getInstance(options, MPI_COMM_WORLD, MPI_COMM_WORLD); 
  const int Nthreads =  omp_get_max_threads();

  // build+load kernel
  occa::properties props = platform->kernelInfo + meshKernelProperties(N);
  static constexpr int nFields = 3;
  props["defines/p_cubNq"] = cubNq;
  props["defines/p_cubNp"] = cubNp;
  props["defines/p_nEXT"] = nEXT;
  props["defines/p_NVfields"] = nFields;
  props["defines/p_MovingMesh"] = 0;

  std::string kernelName;
  if(dealias){
    kernelName = "subCycleStrongCubatureVolumeHex3D";
  } else {
    kernelName = "subCycleStrongVolumeHex3D";
  }

  const std::string ext = (platform->device.mode() == "Serial") ? ".c" : ".okl";
  std::string fileName = 
    installDir + "/okl/nrs/" + kernelName + ext;
  
  // currently lacking a native implementation of the non-dealiased kernel
  if(!dealias) fileName = installDir + "/okl/nrs/subCycleHex3D.okl";

  subcyclingKernel = platform->device.buildKernel(fileName, props, true);

  // populate arrays

  dlong* elementList = (dlong*) calloc(Nelements, sizeof(dlong));
  for(int e = 0; e < Nelements; ++e){
    elementList[e] = e;
  }
  o_elementList = platform->device.malloc(Nelements * sizeof(dlong), elementList);
  free(elementList);


  randAlloc = &rand64Alloc; 

  void *invLMM   = randAlloc(Nelements * Np);
  void *cubD  = randAlloc(cubNq * cubNq);
  void *NU  = randAlloc(nFields * fieldOffset);
  void *conv  = randAlloc(nFields * cubatureOffset * nEXT);
  void *cubInterpT  = randAlloc(Nq * cubNq);
  void *Ud  = randAlloc(nFields * fieldOffset);

  o_invLMM = platform->device.malloc(Nelements * Np * wordSize, invLMM);
  free(invLMM);
  o_cubD = platform->device.malloc(cubNq * cubNq * wordSize, cubD);
  free(cubD);
  o_NU = platform->device.malloc(nFields * fieldOffset * wordSize, NU);
  free(NU);
  o_conv = platform->device.malloc(nFields * cubatureOffset * nEXT * wordSize, conv);
  free(conv);
  o_cubInterpT = platform->device.malloc(Nq * cubNq * wordSize, cubInterpT);
  free(cubInterpT);
  o_Ud = platform->device.malloc(nFields * fieldOffset * wordSize, Ud);
  free(Ud);

  // warm-up
  double elapsed = run(10);
  const int elapsedTarget = 10;
  if(Ntests < 0) Ntests = elapsedTarget/elapsed;

  // ***** 
  elapsed = run(Ntests);
  // ***** 
 
  // print statistics
  const dfloat GDOFPerSecond = nFields * (size * Nelements * (N * N * N) / elapsed) / 1.e9;

  size_t bytesPerElem = 2 * nFields * Np; // Ud, NU
  bytesPerElem += Np; // inv mass matrix
  bytesPerElem += nFields * cubNp * nEXT; // U(r,s,t)

  size_t otherBytes = cubNq * cubNq; // D
  if(cubNq > Nq){
    otherBytes += Nq * cubNq; // interpolator
  }
  otherBytes   *= wordSize;
  bytesPerElem *= wordSize;
  const double bw = (size * (Nelements * bytesPerElem + otherBytes) / elapsed) / 1.e9;

  double flopCount = 0.0; // per elem basis
  if(cubNq > Nq){
    flopCount += 6. * cubNp * nEXT; // extrapolate U(r,s,t) to current time
    flopCount += 18. * cubNp * cubNq; // apply Dcub
    flopCount += 9. * Np; // compute NU
    flopCount += 12. * Nq * (cubNp + cubNq * cubNq * Nq + cubNq * Nq * Nq); // interpolation
  } else {
    flopCount = Nq * Nq * Nq * (18. * Nq + 6. * nEXT + 24.);
  }
  const double gflops = (size * flopCount * Nelements / elapsed) / 1.e9;

  if(rank == 0)
    std::cout << "MPItasks=" << size
              << " OMPthreads=" << Nthreads
              << " NRepetitions=" << Ntests
              << " N=" << N
              << " cubN=" << cubN
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
