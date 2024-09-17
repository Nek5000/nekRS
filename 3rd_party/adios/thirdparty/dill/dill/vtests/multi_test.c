#include <stdio.h>
#include "dill.h"
#include <stdlib.h>
#include <string.h>

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

int main(int argc, char **argv) 
{ 
    int verbose = 0;
    int i;
    for (i=1; i < argc; i++) {
	if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	}
    }

    {
	dill_stream s = dill_create_stream();
	int (*func1)();
	int (*func2)();
	int (*func3)();
	dill_exec_handle handle1, handle2, handle3;
	dill_start_simple_proc(s, "foo", DILL_I);
	dill_retii(s, 1);
	handle1 = dill_finalize(s);
	dill_free_handle(handle1);
	handle1 = dill_get_handle(s);
	func1 = (int(*)())dill_get_fp(handle1);

	dill_start_simple_proc(s, "foo", DILL_I);
	dill_retii(s, 2);
	handle2 = dill_finalize(s);
	dill_free_handle(handle2);
	handle2 = dill_get_handle(s);
	func2 = (int(*)())dill_get_fp(handle2);


	dill_start_simple_proc(s, "foo", DILL_I);
	dill_retii(s, 3);
	handle3 = dill_finalize(s);
	dill_free_handle(handle3);
	handle3 = dill_get_handle(s);
	func3 = (int(*)())dill_get_fp(handle3);


	dill_free_stream(s);
	
	assert(func1() == 1);
	assert(func2() == 2);
	assert(func3() == 3);

	dill_free_handle(handle2);

	assert(func1() == 1);
	assert(func3() == 3);

	dill_ref_handle(handle1);

	assert(func1() == 1);

	dill_free_handle(handle1);

	func1 = (int(*)())dill_get_fp(handle1);
	assert(func1() == 1);

	dill_free_handle(handle1);

	assert(func3() == 3);

	dill_free_handle(handle3);
    }
    return 0;
}
