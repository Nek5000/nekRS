#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "dill.h"

int main(int argc, char **argv) 
{ 
    dill_reg r;
    int i;
    dill_stream s = dill_create_raw_stream();
    dill_exec_handle h;
    void (*func)();
    int verbose = 0, no_float = 0;

    for (i=1; i < argc; i++) {
	if (strcmp(argv[i], "-no_float") == 0) {
	    no_float++;
	} else if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	}
    }

    dill_start_simple_proc(s, "foo", DILL_V);
    dill_raw_getreg(s, &r, DILL_I, DILL_TEMP);
	
    dill_push_init(s);
    if (!dill_do_reverse_vararg_push(s)) {
	dill_push_argpi(s, "Hello: %d %d %d %d\n");
	dill_seti(s, r, 10);
	dill_push_argi(s, r);
	dill_seti(s, r,20);
	dill_push_argi(s, r);
	dill_seti(s, r,30);
	dill_push_argi(s, r);
	dill_seti(s, r,40);
	dill_push_argi(s, r);
    } else {
	dill_seti(s, r,40);
	dill_push_argi(s, r);
	dill_seti(s, r,30);
	dill_push_argi(s, r);
	dill_seti(s, r,20);
	dill_push_argi(s, r);
	dill_seti(s, r, 10);
	dill_push_argi(s, r);
	dill_push_argpi(s, "Hello: %d %d %d %d\n");
    }
    dill_callv(s, (void*)printf, "printf");
    h = dill_finalize(s);
    func = (void (*)())dill_get_fp(h);
    if (verbose) dill_dump(s);
    func();
    dill_free_handle(h);

    dill_start_simple_proc(s, "foo", DILL_V);
    dill_push_init(s);
    if (!dill_do_reverse_vararg_push(s)) {
	dill_push_argpi(s, "Hello: %d %d %d %d\n");
	dill_push_argii(s, 10);
	dill_push_argii(s, 20);
	dill_push_argii(s, 30);
	dill_push_argii(s, 40);
    } else {
	dill_push_argii(s, 40);
	dill_push_argii(s, 30);
	dill_push_argii(s, 20);
	dill_push_argii(s, 10);
	dill_push_argpi(s, "Hello: %d %d %d %d\n");
    }
    dill_callv(s, (void*)printf, "printf");
    h = dill_finalize(s);
    func = (void (*)()) dill_get_fp(h);
    if (verbose) dill_dump(s);
    func();
    dill_free_handle(h);

    dill_start_simple_proc(s, "foo", DILL_V);
    dill_push_init(s);
    if (!dill_do_reverse_vararg_push(s)) {
	dill_push_argpi(s, "Hello: %d %d %d %d %d %d %d %d %d %d\n");
	dill_push_argii(s, 10);
	dill_push_argii(s, 20);
	dill_push_argii(s, 30);
	dill_push_argii(s, 40);
	dill_push_argii(s, 50);
	dill_push_argii(s, 60);
	dill_push_argii(s, 70);
	dill_push_argii(s, 80);
	dill_push_argii(s, 90);
	dill_push_argii(s, 100);
    } else {
	dill_push_argii(s, 100);
	dill_push_argii(s, 90);
	dill_push_argii(s, 80);
	dill_push_argii(s, 70);
	dill_push_argii(s, 60);
	dill_push_argii(s, 50);
	dill_push_argii(s, 40);
	dill_push_argii(s, 30);
	dill_push_argii(s, 20);
	dill_push_argii(s, 10);
	dill_push_argpi(s, "Hello: %d %d %d %d %d %d %d %d %d %d\n");
    }
    dill_callv(s, (void*)printf, "printf");
    h = dill_finalize(s);
    func = (void (*)())dill_get_fp(h);
    if (verbose) dill_dump(s);
    func();
    dill_free_handle(h);

    if (!no_float) {
	dill_start_simple_proc(s, "foo", DILL_V);
	dill_push_varidiac_init(s, 1);
	if (!dill_do_reverse_vararg_push(s)) {
	    dill_push_argpi(s, "Hello: %e %e %e %e %e %e %e %e %e %e\n");
	    dill_push_argdi(s, 10.0);
	    dill_push_argdi(s, 20.0);
	    dill_push_argdi(s, 30.0);
	    dill_push_argdi(s, 40.0);
	    dill_push_argdi(s, 50.0);
	    dill_push_argdi(s, 60.0);
	    dill_push_argdi(s, 70.0);
	    dill_push_argdi(s, 80.0);
	    dill_push_argdi(s, 90.0);
	    dill_push_argdi(s, 100.0);
	} else {
	    dill_push_argdi(s, 100.0);
	    dill_push_argdi(s, 90.0);
	    dill_push_argdi(s, 80.0);
	    dill_push_argdi(s, 70.0);
	    dill_push_argdi(s, 60.0);
	    dill_push_argdi(s, 50.0);
	    dill_push_argdi(s, 40.0);
	    dill_push_argdi(s, 30.0);
	    dill_push_argdi(s, 20.0);
	    dill_push_argdi(s, 10.0);
	    dill_push_argpi(s, "Hello: %e %e %e %e %e %e %e %e %e %e\n");
	}
	dill_callv(s, (void*)printf, "printf");
	h = dill_finalize(s);
	func = (void (*)())dill_get_fp(h);
	if (verbose) dill_dump(s);
	func();
	dill_free_handle(h);
    } else {
	/* just do the printf so the output is the same */
	printf("Hello: %e %e %e %e %e %e %e %e %e %e\n", 10.0, 20.0, 30.0, 
	       40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0);
    }
    dill_free_stream(s);
    return 0;
}
