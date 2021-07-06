
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "c99.h"
#include "types.h"
#include "name.h"
#include "fail.h"
#include "mem.h"
#include "obbox.h"
#include "poly.h"
#include "gs_defs.h"
#include "sort.h"
#include "comm.h"
#include "crystal.h"
#include "sarray_transfer.h"
#include "sarray_sort.h"
#include "findpts_el.h"
#include "findpts_local.h"
#include "findpts.h"

#include "ogs_findpts.h"

#define CODE_INTERNAL 0
#define CODE_BORDER 1
#define CODE_NOT_FOUND 2

static slong lfloor(double x) { return floor(x); }
static slong lceil (double x) { return ceil (x); }

static ulong hash_index_aux(double low, double fac, ulong n, double x)
{
  const slong i = lfloor((x-low)*fac);
  return i<0 ? 0 : (n-1<(ulong)i ? n-1 : (ulong)i);
}

#define D 2
#define WHEN_3D(a)
#include "ogs_findpts_imp.h"
#undef WHEN_3D
#undef D

#define D 3
#define WHEN_3D(a) a
#include "ogs_findpts_imp.h"
#undef WHEN_3D
#undef D
