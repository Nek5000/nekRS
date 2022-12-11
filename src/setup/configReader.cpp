#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

#include <stdlib.h>

#include "inipp.hpp"
#include "nrs.hpp"

#define UPPER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::toupper)); \
}
#define LOWER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::tolower)); \
}

void configRead(MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  
  char* nekrs_home = getenv("NEKRS_HOME");
  if (nekrs_home == nullptr) {
    if (rank == 0)
      std::cout << "\nERROR: The environment variable NEKRS_HOME is not defined!\n";
    EXIT_AND_FINALIZE(1);
  }
  std::string installDir{nekrs_home};
  std::string configFile = installDir + "/nekrs.conf";

  const char* ptr = realpath(configFile.c_str(), NULL);
  if (!ptr) {
    if (rank == 0) 
      std::cout << "\nERROR: Cannot find " << configFile << "!\n";
    EXIT_AND_FINALIZE(1);
  }

  char* rbuf;
  long fsize;
  if(rank == 0) {
    FILE* f = fopen(configFile.c_str(), "rb");
    fseek(f, 0, SEEK_END);
    fsize = ftell(f);
    fseek(f, 0, SEEK_SET);
    rbuf = new char[fsize];
    fread(rbuf, 1, fsize, f);
    fclose(f);
  }
  MPI_Bcast(&fsize, sizeof(fsize), MPI_BYTE, 0, comm);
  if(rank != 0) rbuf = new char[fsize];
  MPI_Bcast(rbuf, fsize, MPI_CHAR, 0, comm);
  std::stringstream is;
  is.write(rbuf, fsize);

  inipp::Ini ini;
  ini.parse(is, false);
  ini.interpolate();

  std::string buf;
  ini.extract("general", "cxx", buf);
  setenv("NEKRS_CXX", buf.c_str(), 1);

  ini.extract("general", "cxxflags", buf);
  setenv("NEKRS_CXXFLAGS", buf.c_str(), 1);

  ini.extract("general", "cc", buf);
  setenv("NEKRS_CC", buf.c_str(), 1);

  ini.extract("general", "cflags", buf);
  setenv("NEKRS_CFLAGS", buf.c_str(), 1);

  ini.extract("general", "fc", buf);
  setenv("NEKRS_FC", buf.c_str(), 1);

  ini.extract("general", "fflags", buf);
  setenv("NEKRS_FFLAGS", buf.c_str(), 1);

  ini.extract("general", "nek5000_pplist", buf);
  setenv("NEKRS_NEK5000_PPLIST", buf.c_str(), 1);

  ini.extract("general", "nekrs_mpi_underlying_compiler", buf);
  if(!getenv("NEKRS_MPI_UNDERLYING_COMPILER")) setenv("NEKRS_MPI_UNDERLYING_COMPILER", buf.c_str(), 1);

  ini.extract("general", "nekrs_udf_includes", buf);
  if(getenv("NEKRS_UDF_INCLUDES")) buf = std::string(getenv("NEKRS_UDF_INCLUDES")) + " " + buf;
  if(!buf.empty()) setenv("NEKRS_UDF_INCLUDES", buf.c_str(), 1);

  ini.extract("general", "nekrs_udf_libs", buf);
  if(getenv("NEKRS_UDF_LDFLAGS")) buf = std::string(getenv("NEKRS_UDF_LDFLAGS")) + " " + buf;
  if(!buf.empty()) setenv("NEKRS_UDF_LDFLAGS", buf.c_str(), 1);

  ini.extract("general", "occa_cxx", buf);
  if(!getenv("OCCA_CXX")) setenv("OCCA_CXX", buf.c_str(), 1);

  ini.extract("general", "occa_cxxflags", buf);
  if(!getenv("OCCA_CXXFLAGS")) setenv("OCCA_CXXFLAGS", buf.c_str(), 1);

  ini.extract("general", "occa_cuda_compiler_flags", buf);
  if(!getenv("OCCA_CUDA_COMPILER_FLAGS")) setenv("OCCA_CUDA_COMPILER_FLAGS", buf.c_str(), 1);

  ini.extract("general", "occa_hip_compiler_flags", buf);
  if(!getenv("OCCA_HIP_COMPILER_FLAGS")) setenv("OCCA_HIP_COMPILER_FLAGS", buf.c_str(), 1);

  ini.extract("general", "occa_opencl_compiler_flags", buf);
  if(!getenv("OCCA_OPENCL_COMPILER_FLAGS")) setenv("OCCA_OPENCL_COMPILER_FLAGS", buf.c_str(), 1);

  ini.extract("general", "nekrs_gpu_mpi", buf);
  if(!getenv("NEKRS_GPU_MPI")) setenv("NEKRS_GPU_MPI", buf.c_str(), 1);

  buf = installDir;
  setenv("NEKRS_INSTALL_DIR", buf.c_str(), 1);

  ini.extract("general", "occa_mode_default", buf);
  if(!getenv("NEKRS_OCCA_MODE_DEFAULT")) setenv("NEKRS_OCCA_MODE_DEFAULT", buf.c_str(), 1);

  setenv("OCCA_MEM_BYTE_ALIGN", "1024", 1);
}
