#include "config.h"
#include "data_funcs.h"
#include "cod.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

static int verbose = 0;
#ifdef NO_EMULATION
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();
#define EC_param0
#define EC_param1
#else
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
cod_add_param("ec", "cod_exec_context", 0, x);
#define EC_param0 ec
#define EC_param1 ec,
#define EC_param0_decl cod_exec_context
#define EC_param1_decl cod_exec_context,
#endif

extern void
write_buffer(char *filename, FMStructDescList desc, void *data, 
             int test_num);
extern char *read_buffer(FMContext c, char *read_file, int test_num);

static int count(cod_exec_context ec, long queue) {return queue;}
static int discard(cod_exec_context ec, long queue, long index) {return queue + index;}

int
main(int argc, char**argv)
{
    int test_num = 0;
    int run_only = -1;
    while (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[1], "-o") == 0) {
	    sscanf(argv[2], "%d", &run_only);
	    argc--; argv++;
	}
	argc--; argv++;
    }
    if ((run_only == -1) || (run_only == test_num)) {
	/* test the basics */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	int j = 0;\n\
	int k = 1;\n\
        int or_count = 0;\n\
	int and_count = 0;\n\
	(i || or_count++);\n\
	(k || (or_count+=4));\n\
	(i && (and_count+=8));\n\
	(k && (and_count+=16));\n\
	return or_count + and_count;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 17);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 1  -  goto */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	i = 5;\n\
top:\n\
	i *= 2;\n\
	if (i < 20) {\n\
	   goto top;\n\
	}\n\
	return i;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 20);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 2  -  break */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	for(i=5; i < 40; i++) {\n\
	    if (i > 10) break;\n\
	}\n\
	return i;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 11);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 3  -  break */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	int j = 0;\n\
	int count = 0;\n\
        for (j=0; j< 4; j++) {\n\
	    for(i=0; i < 10; i++) {\n\
		count++;\n\
	        if (i > 3) break;\n\
	    }\n\
	    count +=100;\n\
	}\n\
	return count;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 420);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 4  -  continue */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	int j = 0;\n\
	int count = 0;\n\
        for (j=0; j< 4; j++) {\n\
	    for(i=0; i < 10; i++) {\n\
	        if (i > 4) continue;\n\
		count++;\n\
	    }\n\
	    count +=100;\n\
	}\n\
	return count;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 420);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 5  -  mixed decls and statements */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 5;\n\
	int j = 0;\n\
	int count = 0;\n\
        for (j=0; j< 4; j++) {\n\
	    for(i=0; i < 10; i++) {\n\
	        if (i > 4) continue;\n\
		count++;\n\
	    }\n\
	    count +=100;\n\
	}\n\
	return count;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 420);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 6  -  subroutine closures */
	static char extern_string[] = "int printf(string format, ...);\n"
        "int EVdiscard(cod_exec_context ec, cod_closure_context type, int index);\n"
        "int EVcount(cod_exec_context ec, cod_closure_context type);\n";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"EVdiscard", (void*)(intptr_t)discard},
	    {"EVcount", (void*)(intptr_t)count},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i, j;\n\
	i = EVdiscard(5);\n\
	j = EVcount();\n\
	return i + j;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	cod_set_closure("EVdiscard", (void*) 7, context);
	cod_set_closure("EVcount", (void*) 9, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 5 + 7 + 9);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 7  -  do while */
	static char extern_string[] = "int printf(string format, ...);";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	int i = 0;\n\
	do {\n\
	   i += 7;\n\
	} while (i < 40);\n\
	return i;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 42);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    return 0;
}
