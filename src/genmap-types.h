#ifndef _GENMAP_TYPES_H_
#define _GENMAP_TYPES_H_

#include <genmap-gslib.h>
//
// Genmap types
//
#if defined(GENMAP_LONG_LONG)
  typedef long long GenmapLong;
  typedef unsigned long long GenmapULong;
  #define GenmapLongFormat "%lld"
  #define GenmapULongFormat "%llu"
  
  #if defined(GENMAP_MPI)
  #define GENMAP_LONG MPI_LONG_LONG
  #define GENMAP_UNSIGNED_LONG MPI_UNSIGNED_LONG_LONG
  #else
  #define GENMAP_LONG 0
  #define GENMAP_UNSIGNED_LONG 0
  #endif
  
  #define genmap_gs_long gs_long_long
#else
  typedef long GenmapLong;
  typedef unsigned long GenmapULong;
  #define GenmapLongFormat "%ld"
  #define GenmapULongFormat "%lu"
  
  #if defined(GENMAP_MPI)
  #define GENMAP_LONG MPI_LONG
  #define GENMAP_UNSIGNED_LONG MPI_UNSIGNED_LONG
  #else
  #define GENMAP_LONG 0
  #define GENMAP_UNSIGNED_LONG 0
  #endif

  #define genmap_gs_long gs_long
#endif

typedef int GenmapInt;
typedef unsigned int GenmapUInt;
#define GenmapIntFormat "%d"
#define GenmapUIntFormat "%u"

#if defined(GENMAP_MPI)
#define GENMAP_INT MPI_INT
#define GENMAP_UNSIGNED_INT MPI_UNSIGNED_INT
#else
#define GENMAP_INT 0
#define GENMAP_UNSIGNED_INT 0
#endif

#define genmap_gs_int gs_int

typedef double GenmapScalar;
#define genmap_gs_scalar gs_double
#define GenmapScalarFormat "%lf"

#if defined(GENMAP_MPI)
#define GENMAP_SCALAR MPI_DOUBLE
#else
#define GENMAP_SCALAR 0
#endif

#endif
