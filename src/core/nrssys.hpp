#if !defined(nekrs_nrssys_hpp_)
#define nekrs_nrssys_hpp_

//float data type
#if 0
using dfloat = float;
#define DFLOAT_SINGLE
#define MPI_DFLOAT MPI_FLOAT
#define dfloatFormat "%f"
#define dfloatString "float"
#else
using dfloat = double;
#define DFLOAT_DOUBLE
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

//smoother float data type
#if 1
using pfloat = float;
#define MPI_PFLOAT MPI_FLOAT
#define pfloatFormat "%f"
#define pfloatString "float"
#else
using pfloat = double;
#define MPI_PFLOAT MPI_DOUBLE
#define pfloatFormat "%lf"
#define pfloatString "double"
#endif

//host index data type
#if 0
using hlong = int;
#define MPI_HLONG MPI_INT
#define hlongFormat "%d"
#define hlongString "int"
#else
using hlong = long long int;
#define MPI_HLONG MPI_LONG_LONG_INT
#define hlongFormat "%lld"
#define hlongString "long long int"
#endif

//device index data type
#if 1
using dlong = int;
#define MPI_DLONG MPI_INT
#define dlongFormat "%d"
#define dlongString "int"
#else
using dlong = long long int;
#define MPI_DLONG MPI_LONG_LONG_INT;
#define dlongFormat "%lld"
#define dlongString "long long int"
#endif

// Workaround for https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1

#include <mpi.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <fstream>
#include <limits>
#include <array>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>

#include "occa.hpp"
#include "ogs.hpp"
#include "setupAide.hpp"
#include "timer.hpp"

#define nrsCheck(_nrsCheckCond, _nrsCheckComm, _nrsCheckExitCode, _nrsCheckMessage, ...) \
  do { \
    int _nrsCheckErr = 0; \
    if(_nrsCheckCond) _nrsCheckErr = 1; \
    if(_nrsCheckComm != MPI_COMM_SELF) MPI_Allreduce(MPI_IN_PLACE, &_nrsCheckErr, 1, MPI_INT, MPI_SUM, _nrsCheckComm); \
    if(_nrsCheckErr) { \
      int rank = 0; \
      MPI_Comm_rank(_nrsCheckComm, &rank); \
      if(rank == 0) { \
        fprintf(stderr, "Error in %s: ", __func__);\
        fprintf(stderr, _nrsCheckMessage, __VA_ARGS__); \
      } \
      fflush(stderr); \
      fflush(stdout); \
      MPI_Barrier(_nrsCheckComm); \
      MPI_Abort(MPI_COMM_WORLD, _nrsCheckExitCode); \
    } \
  } while (0)

#define nrsAbort(...) \
  do { \
    nrsCheck(true, __VA_ARGS__); \
  } while (0)


struct platform_t;
extern platform_t* platform;

namespace {

constexpr unsigned int BLOCKSIZE = 256;
constexpr unsigned int ALIGN_SIZE = 1024;
constexpr unsigned int NSCALAR_MAX = 100;
occa::memory o_NULL;

template<typename T>
unsigned int alignStride(unsigned int stride)
{
  const auto pageW = ALIGN_SIZE / sizeof(T);
  if (stride % pageW)
    stride = (stride / pageW + 1) * pageW;

  return stride;
}

const std::string scalarDigitStr(int i)
{
  const int scalarWidth = std::to_string(NSCALAR_MAX - 1).length();
  std::stringstream ss;
  ss << std::setfill('0') << std::setw(scalarWidth) << i;
  return ss.str();
};

// std::to_string might be not accurate enough
std::string to_string_f(double a)
{
  std::stringstream s;
  constexpr auto maxPrecision{std::numeric_limits<double>::digits10 + 1};
  s << std::setprecision(maxPrecision) << std::scientific << a;
  return s.str();
}

std::vector<std::string> serializeString(const std::string sin, char dlim)
{
  std::vector<std::string> slist;
  std::string s(sin);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
  std::stringstream ss;
  ss.str(s);
  while (ss.good()) {
    std::string substr;
    std::getline(ss, substr, dlim);
    if(!substr.empty()) slist.push_back(substr);
  }
  return slist;
}

void lowerCase(std::string& s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
}

void lowerCase(std::vector<std::string>& stringVec)
{
  for(auto && s : stringVec) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
  }
}

void upperCase(std::vector<std::string>& stringVec)
{
  for(auto && s : stringVec) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::toupper(c); });
  }
}

void upperCase(std::string& s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
}

}

#endif
