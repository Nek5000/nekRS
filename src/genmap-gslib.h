#ifndef _GENMAP_GSLIB_H_
#define _GENMAP_GSLIB_H_

// Data type sint/uint
// (defualt) int
// #define USE_LONG long
// #define USE_LONG_LONG long long

// Data type slong/ulong
// (default) int
// #define GLOBAL_LONG long
// #define GLOBAL_LONG_LONG long long
// #define GLOBAL_LONG_LONG

#include "gslib.h"

#if !defined(MPI)
#error "gslib needs to be compiled with MPI"
#endif

#if !defined(GLOBAL_LONG_LONG)
#error "gslib needs to be compiled with GLOBAL_LONG_LONG"
#endif

#define TYPE_INT    0
#define TYPE_LONG   1
#define TYPE_FLOAT  2
#define TYPE_DOUBLE 3

#endif
