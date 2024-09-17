#include "config.h"
#include "cod.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef NO_EMULATION
#define EC_param0
#define EC_param1
#else
#define EC_param0 ec
#define EC_param1 ec,
#endif

static int verbose = 0;

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

int
main(int argc, char **argv)
{
    int test_to_run = -1;
    while (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
	} else if (sscanf(argv[1], "%d", &test_to_run) == 1) {
	    /* got a test_to_run */
	} else {
	    printf("Unknown argument %s\n", argv[1]);
	}
	argc--; argv++;
    }

    /* test of ECL static state saving and restoring */

    if ((test_to_run == 1) || (test_to_run == -1)) {
	/* test dynamic array fields in static vars */
	char code_string[] = "\
{\n\
    static int run_number = 0;\n\
    run_number = run_number + 1;\n\
    return run_number;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)(cod_exec_context);
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	cod_subroutine_declaration("int proc(cod_exec_context ec)", context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    cod_exec_context ec, ec2;
	    cod_code gen_code2;
	    void *state;
	    int state_size = 0;
	    ec = cod_create_exec_context(gen_code);
	    func = (long(*)(cod_exec_context)) (intptr_t) gen_code->func;
	    result = func(ec);
	    if (verbose) printf("Run number after first run is %ld, expected 1\nGrabbing state...", result);
	    assert(result == 1);
	    
	    state = cod_extract_state(ec, &state_size);
	    result = func(ec);
	    if (verbose) printf("Run number after second run is %ld, expected 2\n", result);
	    if (result != 2) {
	      printf("First point, Expected result to be 2, found %ld\n", result);
	    }
	    ec2 = cod_create_exec_context(gen_code);

	    gen_code2 = cod_code_gen(code_string, context);
	    if (verbose) printf("Installing state in second context\n");
	    assert(cod_install_state(ec, state, state_size) == 1);

	    result = func(ec);
	    if (verbose) printf("Run number using state after first run is %ld, expected 2\n", result);
	    if (result != 2) {
	      printf("Second point, Expected result to be 2, found %ld\n", result);
	    }
	    cod_code_free(gen_code);
	    cod_code_free(gen_code2);
	    cod_exec_context_free(ec);
	    cod_exec_context_free(ec2);
	    free(state);
	}
	cod_free_parse_context(context);
    }
    if ((test_to_run == 2) || (test_to_run == -1)) {
	/* test dynamic array fields in static vars */
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    int size;\n\
	    int array[size];\n\
	} test_struct;\n\
\n\
    int i;\n\
    static int first = 1;\n\
    static test_struct t1;\n\
    if (first) {\n\
      first = 0;\n\
      t1.size = 5;\n\
      t1.array[0] = 6;\n\
      t1.array[1] = 5;\n\
      t1.array[2] = 4;\n\
      t1.array[3] = 3;\n\
      t1.array[4] = 2;\n\
      return 0;\n\
    } else {\n\
      int sum = 0;\n\
      int i;\n\
      for (i = 0; i < t1.size ; i++) {\n\
	sum = sum + t1.array[i];\n\
      }\n\
      return sum;\n\
    }\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)(cod_exec_context);
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	cod_subroutine_declaration("int proc(cod_exec_context ec)", context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    cod_exec_context ec, ec2;
	    cod_code gen_code2;
	    void *state;
	    int state_size = 0;
	    ec = cod_create_exec_context(gen_code);
	    func = (long(*)(cod_exec_context)) (intptr_t) gen_code->func;
	    result = func(ec);
	    assert(result == 0);
	    
	    state = cod_extract_state(ec, &state_size);
	    ec2 = cod_create_exec_context(gen_code);

	    gen_code2 = cod_code_gen(code_string, context);
	    assert(cod_install_state(ec, state, state_size) == 1);

	    result = func(ec);
	    assert(result == 20);
	    free(state);
	    cod_code_free(gen_code);
	    cod_code_free(gen_code2);
	    cod_exec_context_free(ec);
	    cod_exec_context_free(ec2);
	}
	cod_free_parse_context(context);
    }

    return 0;
}

