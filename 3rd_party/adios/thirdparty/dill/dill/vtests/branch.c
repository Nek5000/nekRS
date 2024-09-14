#include "stdio.h"
#include "dill.h"

static int verbose = 0;

int p() {
    static int a = 0;
    return a++;
}

int 
main(int argc, char **argv)
{
    if (argc > 1) verbose++;
    dill_stream s = dill_create_stream();
    dill_exec_handle handle;
    int (*ip)();
    int failure = 0;

    dill_start_simple_proc(s, "a_gen", DILL_I);
    int end_label, check_success;
    int (*proc)() = &p;
    dill_reg reg, ret_reg = dill_getreg(s, DILL_I);
    reg = dill_scalli(s, proc, "P", "");
    dill_movi(s, ret_reg, reg);
    check_success = dill_alloc_label(s, "check_success");
    dill_bneii(s, ret_reg, (int)1, check_success);
    dill_retii(s, -1);
    dill_mark_label(s, check_success);
    end_label = dill_alloc_label(s, "end_label");
    dill_bneii(s, ret_reg, (int)2, end_label);
    dill_retii(s, -2);
    dill_mark_label(s, end_label);
    dill_retii(s, -3);
    handle = dill_finalize(s);
    ip = (int(*)())dill_get_fp(handle);
    
    if (verbose) dill_dump(s);
    
    if (ip() != -3) {
	if (verbose) {
	    printf("Failed at point 1\n");
	}
	failure = 1;
    }
    if (ip() != -1)  {
	if (verbose) {
	    printf("Failed at point 2\n");
	}
	failure = 2;
    }
    if (ip() != -2) {
	if (verbose) {
	    printf("Failed at point 3\n");
	}
	failure = 3;
    }
    if (ip() != -3) {
	if (verbose) {
	    printf("Failed at point 4\n");
	}
	failure = 4;
    }
    dill_free_stream(s);
    dill_free_handle(handle);
    return failure;
}
