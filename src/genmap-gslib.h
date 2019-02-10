#ifndef _GENMAP_GSLIB_H_
#define _GENMAP_GSLIB_H_

#if defined(GENMAP_UNDERSCORE)
#  define UNDERSCORE
#endif

#if defined(GENMAP_MPI)
#  define MPI
#endif

#define PREFIX gslib_
#define FPREFIX fgslib_

// Data type sint/uint
// (defualt) int
// #define USE_LONG long
// #define USE_LONG_LONG long long

// Data type slong/ulong
// (default) int
// #define GLOBAL_LONG long
// #define GLOBAL_LONG_LONG long long
// #define GLOBAL_LONG_LONG
#define GLOBAL_LONG_LONG

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "c99.h"
#include "types.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "gs_defs.h"
#include "comm.h"
#include "gs.h"
#include "sort.h"
#include "sarray_sort.h"
#include "crystal.h"
#include "sarray_transfer.h"

#define TYPE_INT    0
#define TYPE_LONG   1
#define TYPE_FLOAT  2
#define TYPE_DOUBLE 3

#endif
