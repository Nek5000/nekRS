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

#include "benchmarkAx.hpp"

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
  int Nelements;
  int Ng = -1;
  int Ndim = 1;
  int okl = 1;
  int BKmode = 0;
  int Ntests = -1;
  size_t wordSize = 8;
  int computeGeom = 0;

  while(1) {
    static struct option long_options[] =
    {
      {"p-order", required_argument, 0, 'p'},
      {"g-order", required_argument, 0, 'g'},
      {"computeGeom", no_argument, 0, 'c'},
      {"block-dim", required_argument, 0, 'd'},
      {"elements", required_argument, 0, 'e'},
      {"backend", required_argument, 0, 'b'},
      {"arch", required_argument, 0, 'a'},
      {"bk-mode", no_argument, 0, 'm'},
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
    case 'g':
      Ng = atoi(optarg); 
      break;
    case 'c':
      computeGeom = 1; 
      break;
    case 'd':
      Ndim = atoi(optarg);
      break;
    case 'e':
      Nelements = atoi(optarg);
      cmdCheck++;
      break;
    case 'b':
      options.setArgs("THREAD MODEL", std::string(optarg));
      cmdCheck++;
      break;
    case 'm':
      BKmode = 1;
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
      printf("Usage: ./nekrs-axhelm  --p-order <n> --elements <n> --backend <CPU|CUDA|HIP|OPENCL>\n"
             "                    [--block-dim <n>]\n"
             "                    [--g-order <n>] [--computeGeom]\n"
             "                    [--bk-mode] [--fp32] [--iterations <n>]\n"); 
    exit(1); 
  }

  if(Ng < 0) Ng = N; 
  Nelements = std::max(1, Nelements/size);
  constexpr int p_Nggeo {7};
  const int Nq = N + 1;
  const int Np = Nq * Nq * Nq;
  const int Nq_g = Ng + 1;
  const int Np_g = Nq_g * Nq_g * Nq_g; 

  // BKmode <-> both constant coeff AND poisson
  bool poisson = false;
  bool constCoeff = false;
  if(BKmode){
    poisson = true;
    constCoeff = true;
  }

  platform = platform_t::getInstance(options, MPI_COMM_WORLD, MPI_COMM_WORLD);
  platform->options.setArgs("BUILD ONLY", "FALSE");
  const int verbosity = 2;
  if (Ntests != -1) {
    benchmarkAx(Nelements,
                Nq,
                Ng,
                poisson,
                constCoeff,
                computeGeom,
                wordSize,
                Ndim,
                false, // no stress formulation
                verbosity,
                Ntests,
                true,
                "");
  }
  else {
    const double targetTime = 10.0;
    benchmarkAx(Nelements,
                Nq,
                Ng,
                poisson,
                constCoeff,
                computeGeom,
                wordSize,
                Ndim,
                false, // no stress formulation
                verbosity,
                targetTime,
                true,
                "");
  }
  MPI_Finalize();
  exit(0);
}
