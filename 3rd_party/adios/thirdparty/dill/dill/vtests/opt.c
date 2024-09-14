#include "../config.h"
#include "stdio.h"
#ifdef HAVE_MALLOC_H
#include "malloc.h"
#endif
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include "unistd.h"
#endif

#include "dill.h"
#ifdef USE_MMAP_CODE_SEG
#include "sys/mman.h"
#endif

static int verbose = 0;

void a () {
    dill_stream s = dill_create_stream();
    dill_reg a;
    dill_exec_ctx ec;
    dill_exec_handle handle;
    int (*ip)(dill_exec_ctx, int);
    int result;
    dill_start_proc(s, "a_gen", DILL_I, "%EC%i");
    
    a = dill_getreg(s, DILL_I);
    dill_seti(s, a, 0);
    dill_addii(s, a, a, 5);
    dill_addii(s, a, a, 15);
    dill_addi(s, a, a, dill_vparam(s, 1));
    dill_addii(s, a, a, 20);
    dill_reti(s, a);
    handle = dill_finalize(s);
    ip = (int(*)(dill_exec_ctx, int))dill_get_fp(handle);
    
    if (verbose) dill_dump(s);
     
    ec = dill_get_exec_context(s);
    result = (*ip)(ec, 1);
    if (result != 41) {
	printf("A failed.  Expected 41, got %d\n", result);
    }
    dill_free_exec_context(ec);
    dill_free_stream(s);
    dill_free_handle(handle);
}

void b () {
    dill_stream s = dill_create_stream();
    dill_reg a;
    dill_exec_ctx ec;
    dill_exec_handle handle;
    int (*ip)(dill_exec_ctx, int);
    int result;
    dill_start_proc(s, "a_gen", DILL_I, "%EC%i");
    
    a = dill_getreg(s, DILL_I);
    dill_seti(s, a, 0);
    dill_addii(s, a, a, 5);
    dill_addii(s, a, a, 15);
    dill_addii(s, a, a, 20);
    dill_addii(s, a, a, 1);
    dill_reti(s, a);
    handle = dill_finalize(s);
    ip = (int(*)(dill_exec_ctx, int))dill_get_fp(handle);
    
    if (verbose) dill_dump(s);
     
    ec = dill_get_exec_context(s);
    result = (*ip)(ec, 1);
    if (result != 41) {
	printf("B failed.  Expected 41, got %d\n", result);
    }
    dill_free_exec_context(ec);
    dill_free_stream(s);
    dill_free_handle(handle);
}

int 
main(int argc, char **argv)
{
    if (argc > 1) verbose++;
    if (verbose) printf("########## A\n");
    a();
    if (verbose) printf("########## B\n");
    b();
    if (verbose) printf("########## end\n");
    return 0;
}
