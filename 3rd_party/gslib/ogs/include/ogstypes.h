#if !defined(ogstypes_h)
#define ogstypes_h

#define BLOCKSIZE 256

// float data type
#if 0
typedef float dfloat;
#define DFLOAT_SINGLE
#define MPI_DFLOAT MPI_FLOAT
#define dfloatFormat "%f"
#define dfloatString "float"
#else
typedef double dfloat;
#define DFLOAT_DOUBLE
#define MPI_DFLOAT MPI_DOUBLE
#define dfloatFormat "%lf"
#define dfloatString "double"
#endif

// host index data type
#if 0
typedef int hlong;
#define MPI_HLONG MPI_INT
#define hlongFormat "%d"
#define hlongString "int"
#else
typedef long long int hlong;
#define MPI_HLONG MPI_LONG_LONG_INT
#define hlongFormat "%lld"
#define hlongString "long long int"
#endif

// device index data type
#if 1
typedef int dlong;
#define MPI_DLONG MPI_INT
#define dlongFormat "%d"
#define dlongString "int"
#else
typedef long long int dlong;
#define MPI_DLONG MPI_LONG_LONG_INT;
#define dlongFormat "%lld"
#define dlongString "long long int"
#endif

#endif
