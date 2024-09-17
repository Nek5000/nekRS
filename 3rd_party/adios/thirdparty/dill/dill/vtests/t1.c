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
    char *target;
    dill_reg a,b,p3,d,e,f;
    dill_exec_ctx ec;
    dill_exec_handle handle;
    int (*ip)(dill_exec_ctx, int , int );

     dill_start_proc(s, "a_gen", DILL_I, "%EC%i%i");

     a = dill_vparam(s, 1);
     b = dill_vparam(s, 2);
     p3 = dill_getreg(s, DILL_I);
     d = dill_getreg(s, DILL_I);
     e = dill_getreg(s, DILL_I);
     f = dill_getreg(s, DILL_I);

     dill_addii(s, p3, a, 5);
     dill_addi(s, d, a, b);
     dill_addi(s, e, d, p3);
     dill_movi(s, f, e);

     p3 = dill_getreg(s, DILL_I);
     d = dill_getreg(s, DILL_I);
     e = dill_getreg(s, DILL_I);
     dill_addii(s, p3, a, 5);
     dill_addi(s, d, a, b);
     dill_addi(s, e, d, p3);

     dill_addi(s, f, f, e);
     dill_reti(s, f);
     handle = dill_finalize(s);
     ip = (int(*)(dill_exec_ctx, int, int))dill_get_fp(handle);

#ifdef USE_MMAP_CODE_SEG
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif
     {
	 int size = dill_code_size(s);
	 static long ps = -1;
	 if (ps == -1) {
	     ps = (getpagesize ());
	 }
	 if (ps > size) size = ps;
	 target = (void*)mmap(0, size, 
			      PROT_EXEC | PROT_READ | PROT_WRITE, 
			      MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
     }
     if (target == (void*)-1) perror("mmap");
#else
     target = (void*)malloc(dill_code_size(s));
#endif
     if (verbose) dill_dump(s);

     ec = dill_get_exec_context(s);
     printf("**18=%d\n", (*ip)(ec, 1, 2));
     dill_free_exec_context(ec);
     dill_free_stream(s);
     dill_free_handle(handle);
}

int 
main(int argc, char **argv)
{
    if (argc > 1) verbose++;
    printf("########## A\n");
    a();
    printf("########## end\n");
    return 0;
}
