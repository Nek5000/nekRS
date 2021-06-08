#if !defined(ogstypes_h)
#define ogstypes_h

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

#endif
