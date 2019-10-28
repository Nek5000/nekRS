#ifndef _GENMAP_TYPES_H_
#define _GENMAP_TYPES_H_

#include <genmap-gslib.h>

typedef long long GenmapLong;
typedef unsigned long long GenmapULong;
#define GenmapLongFormat "%lld"
#define GenmapULongFormat "%llu"

#define GENMAP_LONG MPI_LONG_LONG
#define GENMAP_UNSIGNED_LONG MPI_UNSIGNED_LONG_LONG

#define genmap_gs_long gs_long_long

typedef int GenmapInt;
typedef unsigned int GenmapUInt;
#define GenmapIntFormat "%d"
#define GenmapUIntFormat "%u"

#define GENMAP_INT MPI_INT
#define GENMAP_UNSIGNED_INT MPI_UNSIGNED_INT

#define genmap_gs_int gs_int

typedef double GenmapScalar;
#define genmap_gs_scalar gs_double
#define GenmapScalarFormat "%lf"

#define GENMAP_SCALAR MPI_DOUBLE

#endif
