#define UNDERSCORE 1
#define USE_NAIVE_BLAS 
#define NO_NEX_EXITT 1
#define GLOBAL_LONG_LONG 1

#define MPI 1

#ifdef __cplusplus
extern "C" {
#endif
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "gs_defs.h"
#include "comm.h"
#include "gs.h"
#include "crystal.h"
#include "sort.h"
#include "sarray_sort.h"
#include "sarray_transfer.h"
#include "tensor.h"
#include "poly.h"
#include "lob_bnd.h"
#include "obbox.h"
#include "findpts.h"
#include "findpts_el.h"
#include "findpts_local.h"
#ifdef __cplusplus
}
#endif
