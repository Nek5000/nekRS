#if !defined(nekrs_nrssys_hpp_)
#define nekrs_nrssys_hpp_

#define BLOCKSIZE 256

//float data type
#if 0
#define DFLOAT_SINGLE
#define dfloat float
#define MPI_DFLOAT MPI_FLOAT
#define dfloatFormat "%f"
#define dfloatString "float"
#else
#define DFLOAT_DOUBLE
#define dfloat double
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

//smoother float data type
#if 1
#define pfloat float
#define MPI_PFLOAT MPI_FLOAT
#define pfloatFormat "%f"
#define pfloatString "float"
#else
#define pfloat double
#define MPI_PFLOAT MPI_DOUBLE
#define pfloatFormat "%lf"
#define pfloatString "double"
#endif

//host index data type
#if 0
#define hlong int
#define MPI_HLONG MPI_INT
#define hlongFormat "%d"
#define hlongString "int"
#else
#define hlong long long int
#define MPI_HLONG MPI_LONG_LONG_INT
#define hlongFormat "%lld"
#define hlongString "long long int"
#endif

//device index data type
#if 1
#define dlong int
#define MPI_DLONG MPI_INT
#define dlongFormat "%d"
#define dlongString "int"
#else
#define dlong long long int
#define MPI_DLONG MPI_LONG_LONG_INT
#define dlongFormat "%lld"
#define dlongString "long long int"
#endif

// Workaround for https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1

#include <mpi.h>
#include "occa.hpp"
#include "ogs.hpp"
#include "setupAide.hpp"

#define NEKRS_VERSION "20"
#define NEKRS_SUBVERSION "1"

#define EXIT(a)  { fflush(stdout); MPI_Finalize(); exit(a); }
#define ABORT(a) { fflush(stdout); MPI_Abort(MPI_COMM_WORLD,a); }

#endif
