#if !defined(nekrs_nrssys_hpp_)
#define nekrs_nrssys_hpp_

#define BLOCKSIZE 256
#define ALIGN_SIZE 4096

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
#include "omp.h"
#include <limits>
#include <string>
#include "occa.hpp"
#include "ogs.hpp"
#include "setupAide.hpp"

static occa::memory o_NULL;

struct platform_t;
extern platform_t* platform;

bool useNodeLocalCache();

#define EXIT_AND_FINALIZE(a)  { fflush(stdout); MPI_Finalize(); exit(a); }
#define ABORT(a)  { fflush(stdout); MPI_Abort(MPI_COMM_WORLD, a); }

#endif
