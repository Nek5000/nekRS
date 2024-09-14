/*
 *   cod - T3
 *     
 *       static variables?  Simple stuff with a cod_exec_context.  Static array.
 */

#include "config.h"
#include "cod.h"
#include <string.h>
#include <stdlib.h>
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>

static double testd(){return 1.0;}
static int testi(){return 4;}

#ifdef NO_EMULATION
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
if (output_file) cod_set_error_func(x, error_func);
#define EC_param0
#define EC_param1
#define EC_param0_decl
#define EC_param1_decl
#else
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
cod_add_param("ec", "cod_exec_context", 0, x);\
if (output_file) cod_set_error_func(x, error_func);
#define EC_param0 ec
#define EC_param1 ec,
#define EC_param0_decl cod_exec_context ec
#define EC_param1_decl cod_exec_context ec,
#endif

static int verbose = 0;
static FILE *output_file;

static void
error_func(void *client_data, char *string)
{
    fprintf(output_file, "%s", string);
}

int
main(int argc, char **argv)
{
    int test_to_run = -1;

    while (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[1], "-o") == 0) {
	    sscanf(argv[2], "%d", &test_to_run);
	    argc--; argv++;
	} else if (strcmp(argv[1], "-output") == 0) {
	    output_file = fopen(argv[2], "w");
	    if (!output_file) {
		printf("Couldn't open output file \"%s\"\n", argv[2]);
		exit(1);
	    }
	    argc--; argv++;
	}
	argc--; argv++;
    }
    if ((test_to_run == 1) || (test_to_run == -1)) {
	/* test the basics */
	char code_string[] = "\
{\n\
    static int j = 4;\n\
    static long k = 10;\n\
    static short l = 3;\n\
    static int m = 0;\n\
 static float thresh = 1.5;\n\
    static int base_len;\n\
    static int epoch;      \n\
    static float * base_csym;\n\
    static int have_data = 0;\n\
\n\
    int i;\n\
    int condition;\n\
\n\
    j = j + 1;\n\
    k = k + 2;\n\
    l = l + 3;\n\
\n\
    return j + k + l + m;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	cod_code gen_code;
	long ret;
	long (*func)(EC_param0_decl);

	if (verbose) printf("Running test 1 (-o 1)\n");
	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	ret = func(EC_param0);
	assert(ret == 23);
	ret = func(EC_param0);
	assert(ret == 29);
	ret = func(EC_param0);
	assert(ret == 35);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 2) || (test_to_run == -1)) {

	/* test the ability to have a parameter */
	char code_string[] = "\
{\n\
    static int j = 4;\n\
    static long k = 10;\n\
    static short l = 3;\n\
\n\
    return l * (j + k + i);\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	cod_code gen_code;
    	long (*func)(EC_param1_decl int);

	if (verbose) printf("Running test 2 (-o 2)\n");
#ifdef NO_EMULATION
	cod_subroutine_declaration("int proc(int i)", context);
#else
	cod_subroutine_declaration("int proc(cod_exec_context ec, int i)", context);
#endif
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param1_decl int)) (intptr_t) gen_code->func;
        assert(func(EC_param1 15) == 87);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 3) || (test_to_run == -1)) {
	/* structured types */
	char code_string[] = "\
{\n\
    return input.l * (input.j + input.k + input.i);\n\
}";

	typedef struct test {
	    int i;
	    int j;
	    long k;
	    short l;
	} test_struct, *test_struct_p;

	FMField struct_fields[] = {
	    {"i", "integer", sizeof(int), FMOffset(test_struct_p, i)},
	    {"j", "integer", sizeof(int), FMOffset(test_struct_p, j)},
	    {"k", "integer", sizeof(long), FMOffset(test_struct_p, k)},
	    {"l", "integer", sizeof(short), FMOffset(test_struct_p, l)},
	    {(void*)0, (void*)0, 0, 0}};

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	test_struct str;
	cod_code gen_code;
	long (*func)(EC_param1_decl test_struct_p);

	if (verbose) printf("Running test 3 (-o 3)\n");
	cod_add_simple_struct_type("struct_type", struct_fields, context);
#ifdef NO_EMULATION
	cod_subroutine_declaration("int proc(struct_type *input)", context);
#else
	cod_subroutine_declaration("int proc(cod_exec_context ec, struct_type *input)", context);
#endif

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param1_decl test_struct_p)) (intptr_t) gen_code->func;

	str.i = 15;
	str.j = 4;
	str.k = 10;
	str.l = 3;
	assert(func(EC_param1 &str) == 87);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if ((test_to_run == 4) || (test_to_run == -1)) {
	/* structured types */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
char code_string[] = {"\
  {\n\
    static int localSum = 0;\n\
    localSum = localSum + input.i;\n\
    input.i = localSum;\n\
    return 1;\n\
  }\n\
"};
	typedef struct test {
	    int i;
	    int j;
	    long k;
	    short l;
	} test_struct, *test_struct_p;

	FMField struct_fields[] = {
	    {"i", "integer", sizeof(int), FMOffset(test_struct_p, i)},
	    {"j", "integer", sizeof(int), FMOffset(test_struct_p, j)},
	    {"k", "integer", sizeof(long), FMOffset(test_struct_p, k)},
	    {"l", "integer", sizeof(short), FMOffset(test_struct_p, l)},
	    {(void*)0, (void*)0, 0, 0}};

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	test_struct str;
	cod_code gen_code;
	long (*func)(EC_param1_decl test_struct_p);

	if (verbose) printf("Running test 4 (-o 4)\n");
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	cod_add_simple_struct_type("struct_type", struct_fields, context);
#ifdef NO_EMULATION
	cod_subroutine_declaration("int proc(struct_type *input)", context);
#else
	cod_subroutine_declaration("int proc(cod_exec_context ec, struct_type *input)", context);
#endif

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param1_decl test_struct_p)) (intptr_t) gen_code->func;

	str.i = 15;
	str.j = 4;
	str.k = 10;
	str.l = 3;
	func(EC_param1 &str);
	assert(func(EC_param1 &str) == 1);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if ((test_to_run == 5) || (test_to_run == -1)) {
	static char code[] = "{\
		    int i;\
		    int j;\
		    double sum = 0.0;\
		    double average = 0.0;\
		    for(i = 0; i<37; i= i+1) {\
		        for(j = 0; j<253; j=j+1) {\
			sum = sum + input.levels[j][i];\
		        }\
		    }\
		    average = sum / (37 * 253);\
		    return average;\
		}";

	static FMField input_field_list[] =
	{
	    {"levels", "float[253][37]", sizeof(double), 0},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	int i, j;
	double levels[253][37], result;
	cod_code gen_code;
	double (*func)(EC_param1_decl double*);


	if (verbose) printf("Running test 5 (-o 5)\n");
	cod_add_simple_struct_type("input_type", input_field_list, context);
#ifdef NO_EMULATION
	cod_subroutine_declaration("double proc(input_type *input)", context);
#else
	cod_subroutine_declaration("double proc(cod_exec_context ec, input_type *input)", context);
#endif

	for(i=0; i< 253; i++) {
	    for (j=0; j< 37; j++) {
	        levels[i][j] = i + 1000*j;
	    }
	}

	gen_code = cod_code_gen(code, context);
	ec = cod_create_exec_context(gen_code);
	func = (double (*)(EC_param1_decl double*))(intptr_t) gen_code->func;
	result = (func)(EC_param1 &levels[0][0]);
	if (result != 18126.00) {
	    printf("Expected %g, got %g\n", 18126.0, result);
	}
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if ((test_to_run == 6) || (test_to_run == -1)) {
	static char extern_string[] = "int printf(string format, ...);\
					double testd();";
	static cod_extern_entry externs[] = 
	{
	    {"testd", (void*)(intptr_t)testd},
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
				   double b = testd();\n\
				   return b;\
		}";

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	double result;
	cod_code gen_code;
	double (*func)(EC_param0_decl);

	if (verbose) printf("Running test 6 (-o 6)\n");
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

#ifdef NO_EMULATION
	cod_subroutine_declaration("double proc()", context);
#else
	cod_subroutine_declaration("double proc(cod_exec_context ec)", context);
#endif
	gen_code = cod_code_gen(code, context);
	ec = cod_create_exec_context(gen_code);
	func = (double (*)(EC_param0_decl))(intptr_t) gen_code->func;
	result = (func)(EC_param0);
	if (result != 1.0) {
	    printf("Expected %g, got %g\n", 1.0, result);
	}
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 7) || (test_to_run == -1)) {
	/* test static arrays */
	char code_string[] = "\
{\n\
    static int n[2*2];\n\
    if (n[0] + n[1] + n[2] + n[3] == 0) {\n\
        /* first time */\n\
        n[0] = 4;\n\
        n[1] = 10;\n\
        n[2] = 3;\n\
        n[3] = 0;\n\
    }\n\
    n[0] = n[0] + 1;\n\
    n[1] = n[1] + 2;\n\
    n[2] = n[2] + 3;\n\
    return n[0] + n[1] + n[2] + n[3];\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	cod_code gen_code;
	long (*func)(EC_param0_decl);

	if (verbose) printf("Running test 7 (-o 7)\n");
	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	assert(func(EC_param0) == 23);
	assert(func(EC_param0) == 29);
	assert(func(EC_param0) == 35);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 8) || (test_to_run == -1)) {
	static char extern_string[] = "int printf(string format, ...);\
					int testi();";
	static cod_extern_entry externs[] = 
	{
	    {"testi", (void*)(intptr_t)testi},
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
				   static int count = 0;\n\
				   return count % testi();\n\
		}";

	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;
	int result;
	cod_code gen_code;
	int (*func)(EC_param0_decl);

	if (verbose) printf("Running test 8 (-o 8)\n");
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

#ifdef NO_EMULATION
	cod_subroutine_declaration("int proc()", context);
#else
	cod_subroutine_declaration("int proc(cod_exec_context ec)", context);
#endif
	gen_code = cod_code_gen(code, context);
	ec = cod_create_exec_context(gen_code);
	func = (int (*)(EC_param0_decl))(intptr_t) gen_code->func;
	result = (func)(EC_param0);
	if (result != 0) {
	    printf("Expected %d, got %d\n", 0, result);
	}
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    return 0;
}
