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

#define UPPER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::toupper)); }
#define LOWER(a)  { transform(a.begin(), a.end(), a.begin(), std::ptr_fun<int, int>(std::tolower)); }

void configRead(MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  
  string install_dir;
  install_dir.assign(getenv("NEKRS_HOME"));
  string configFile = install_dir + "/nekrs.conf";

  const char *ptr = realpath(configFile.c_str(), NULL);
  if (!ptr) {
     if (rank == 0) cout << "\nERROR: Cannot find " << configFile << "!\n";
     ABORT(1);
  }

  char *rbuf;
  long fsize; 
  if(rank == 0) {
    FILE *f = fopen(configFile.c_str(), "rb");
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
  stringstream is;
  is.write(rbuf, fsize);

  inipp::Ini<char> ini;
  ini.parse(is, false);
  ini.interpolate();

  string buf;
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

  buf = install_dir + "/nek5000";
  setenv("NEKRS_NEK5000_DIR", buf.c_str(), 1);

  ini.extract("general", "nek5000_pplist", buf);
  setenv("NEKRS_NEK5000_PPLIST", buf.c_str(), 1);

  buf = install_dir + "/libparanumal";
  setenv("NEKRS_LIBP_DIR", buf.c_str(), 1);

  ini.extract("general", "libp_defines", buf);
  setenv("NEKRS_LIBP_DEFINES", buf.c_str(), 1);

  buf = install_dir + "/udf";
  setenv("NEKRS_UDF_DIR", buf.c_str(), 1);

  buf = install_dir + "/nekInterface";
  setenv("NEKRS_NEKINTERFACE_DIR", buf.c_str(), 1);

  ini.extract("general", "occa_cxx", buf);
  if(!getenv("OCCA_CXX")) setenv("OCCA_CXX", buf.c_str(), 1);

  ini.extract("general", "occa_cxxflags", buf);
  if(!getenv("OCCA_CXXFLAGS")) setenv("OCCA_CXXFLAGS", buf.c_str(), 1);

  buf = install_dir + "/occa";
  setenv("OCCA_DIR", buf.c_str(), 1);

  buf = install_dir;
  setenv("NEKRS_INSTALL_DIR", buf.c_str(), 1);
}
