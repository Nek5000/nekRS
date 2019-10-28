#ifndef _GENMAP_GSLIB_H_
#define _GENMAP_GSLIB_H_

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
