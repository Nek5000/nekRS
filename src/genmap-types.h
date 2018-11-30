#ifndef _GENMAP_TYPES_H_
#define _GENMAP_TYPES_H_

#include <genmap-gslib.h>
//
// Genmap types
//
typedef int GenmapInt32;
typedef unsigned int GenmapUInt32;
typedef long GenmapInt64;
typedef unsigned long GenmapUInt64;
typedef long long GenmapLongLong;
typedef unsigned long long GenmapULongLong;

#if defined(GLOBAL_LONG_LONG)
  typedef GenmapLongLong GenmapLong;
  typedef GenmapULongLong GenmapULong;
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
#elif defined(GLOBAL_LONG)
  typedef GenmapInt64 GenmapLong;
  typedef GenmapUInt64 GenmapULong;
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
#else
  typedef GenmapInt32 GenmapLong;
  typedef GenmapUInt32 GenmapULong;
  #define GenmapLongFormat "%d"
  #define GenmapULongFormat "%u"
  
  #if defined(GENMAP_MPI)
  #define GENMAP_LONG MPI_INT
  #define GENMAP_UNSIGNED_LONG MPI_UNSIGNED_INT
  #else
  #define GENMAP_LONG 0
  #define GENMAP_UNSIGNED_LONG 0
  #endif
  
  #define genmap_gs_long gs_int
#endif

#if defined(USE_LONG_LONG)
  typedef GenmapLongLong GenmapInt;
  typedef GenmapULongLong GenmapUInt;
  #define GenmapIntFormat "%lld"
  #define GenmapUIntFormat "%llu"

  #if defined(GENMAP_MPI)
  #define GENMAP_INT MPI_LONG_LONG
  #define GENMAP_UNSIGNED_INT MPI_UNSIGNED_LONG_LONG
  #endif

  #define genmap_gs_int gs_long_long
#elif defined(USE_LONG)
  typedef GenmapInt64 GenmapInt;
  typedef GenmapUInt64 GenmapUInt;
  #define GenmapIntFormat "%ld"
  #define GenmapUIntFormat "%lu"

  #if defined(GENMAP_MPI)
  #define GENMAP_INT MPI_LONG
  #define GENMAP_UNSIGNED_INT MPI_UNSIGNED_LONG
  #else
  #define GENMAP_INT 0
  #define GENMAP_UNSIGNED_INT 0
  #endif

  #define genmap_gs_int gs_long
#else
  typedef GenmapInt32 GenmapInt;
  typedef GenmapUInt32 GenmapUInt;
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
#endif

typedef double GenmapScalar;
#define genmap_gs_scalar gs_double
#define GenmapScalarFormat "%lf"

#if defined(GENMAP_MPI)
#define GENMAP_SCALAR MPI_DOUBLE
#else
#define GENMAP_SCALAR 0
#endif

#endif
