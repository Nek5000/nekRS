#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "omp.h"
#include <unistd.h>
#include "mpi.h"

#include "nrssys.hpp"
#include "setupAide.hpp"
#include "platform.hpp"
#include "configReader.hpp"

#include "benchmarkAdvsub.hpp"

namespace {

int Nfields = 3;
int Nelements;
int Np; 
int cubNp;
dlong fieldOffset;
dlong cubatureOffset;
bool dealias;

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
      {"block-dim", required_argument, 0, 'n'},
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
    case 'n':
      Nfields = atoi(optarg); 
      break;
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
             "                    [--block-dim <n>] [--c-order <n>] [--no-cubature] [--ext-order <n>] [--iterations <n>]\n"); 
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

  if (Nfields != 1 && Nfields != 3){
      printf("Error: Nfields (%d) must be 1 or 3!\n",
        Nfields);
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
  if (cubatureOffset % pageW)
    cubatureOffset = (cubatureOffset / pageW + 1) * pageW;

  platform = platform_t::getInstance(options, MPI_COMM_WORLD, MPI_COMM_WORLD); 
  platform->options.setArgs("BUILD ONLY", "FALSE");
  const int Nthreads =  omp_get_max_threads();

  bool isScalar = Nfields == 1;

  if(Ntests != -1){
    benchmarkAdvsub(Nfields, Nelements, Nq, cubNq, nEXT, dealias, isScalar, 2, Ntests, true);
  } else {
    benchmarkAdvsub(Nfields, Nelements, Nq, cubNq, nEXT, dealias, isScalar, 2, 10.0, true);
  }

  MPI_Finalize();
  exit(0);
}
