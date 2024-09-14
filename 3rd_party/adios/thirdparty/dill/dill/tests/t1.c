#include <assert.h>
#include <stdio.h>
#include "dill.h"
#include <stdlib.h>

#ifdef TESTING
unsigned int x86_64_disassemble(unsigned char *bytes, unsigned int max, int offset, char *output);

int test (int val)
{
    return val;
}
#endif

int main(int argc, char **argv) 
{ 
    dill_stream s = dill_create_raw_stream();
    int (*func)();
    int verbose = 0;

    if (argc > 1) verbose++;

    {
	int result;
	dill_exec_handle h;
	dill_start_simple_proc(s, "foo", DILL_I);
	dill_retii(s, 5);
	h = dill_finalize(s);
	func = (int (*)())dill_get_fp(h);
	if (verbose) dill_dump(s);
	fprintf(stderr, "T1, func is %p\n", (void*) func);
	result = func();
	dill_free_handle(h);

	dill_free_stream(s);
	if (result != 5) {
	    printf("Test 1 failed, got %d instead of 5\n", result);
	    exit(1);
	}
    }
#ifdef TESTING
    {
	dill_stream s = dill_create_raw_stream();
	dill_reg param0;
	dill_exec_handle h;
	int result;
	int (*proc)(int);
	int s2i = 0xdeadbeef;
	dill_start_proc(s, "param1_i", DILL_I, "%i");
	param0 = dill_param_reg(s, 0);
	dill_reti(s, param0);
	h = dill_finalize(s);
	proc = (int (*)(int)) dill_get_fp(h);
	result = proc(s2i);
	dill_dump(s);
	if (result != s2i) {
	    printf("test for 1 arguments of type \"i\" failed, expected %x, got %x\n",
		   s2i, result);
	    dill_dump(s);
	    printf("\n*************\n\n");
	}
    }
    {
	printf("\nDisassembly of Test *************\n\n");
	unsigned char *tmp = (unsigned char*) &test;
	int i = 0;
	for (i = 0; i < 40; ) {
	    char out[128];
	    int ret = x86_64_disassemble(tmp + i, sizeof(out), 0, out);
	    printf("%s\n", out);
	    i += ret;
	}
    }
#endif
    return 0;
}
