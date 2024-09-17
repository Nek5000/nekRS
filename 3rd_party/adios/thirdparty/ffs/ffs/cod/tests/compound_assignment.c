#include "config.h"
#include "cod.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

#ifdef NO_EMULATION
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
if (output_file) cod_set_error_func(x, error_func);
#define EC_param0
#define EC_param1
#else
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
cod_add_param("ec", "cod_exec_context", 0, x);\
if (output_file) cod_set_error_func(x, error_func);
#define EC_param0 ec
#define EC_param1 ec,
#define EC_param0_decl cod_exec_context
#define EC_param1_decl cod_exec_context,
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
	/* test +=, -=, *= */
	char code_string[] = "\
{\n\
    static int j = 4;\n\
    static long k = 10;\n\
    static int l = 3;\n\
\n\
    j += 1;\n\
    k -= 2;\n\
    l *= 2;\n\
\n\
    return j + k + l;\n\
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

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	ret = func(EC_param0);
	assert(ret == 19);
	ret = func(EC_param0);
	assert(ret == 24);
	ret = func(EC_param0);
	assert(ret == 35);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 2) || (test_to_run == -1)) {
	/* test /=, %= */
	char code_string[] = "\
{\n\
    static int m = 1024;\n\
    static int n = 14;\n\
\n\
    m /= 2;\n\
    n %= 5;\n\
    n=n*3;\n\
\n\
    return m + n ;\n\
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

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	ret = func(EC_param0);
	assert(ret == 524);
	ret = func(EC_param0);
	assert(ret == 262);
	ret = func(EC_param0);
	assert(ret == 131);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if ((test_to_run == 3) || (test_to_run == -1)) {
	/* test <<=, >>= */
	char code_string[] = "\
{\n\
    static long o = 0x456789;\n\
    static long p = 0x456789;\n\
\n\
    o <<= 2;\n\
    p >>= 2;\n\
\n\
    return o + p;\n\
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

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	ret = func(EC_param0);
	assert(ret == 0x1159e24 + 0x1159e2);
	ret = func(EC_param0);
	assert(ret == 0x4567890 + 0x45678);
	ret = func(EC_param0);
	assert(ret == 0x1159e240 + 0x1159e);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if ((test_to_run == 4) || (test_to_run == -1)) {
	/* test &=, ^=, |= */
	char code_string[] = "\
{\n\
    static long q = 0x6789;\n\
    static long r = 0x6789;\n\
    static long s = 0x6789;\n\
\n\
printf(\"coming in  q = %x, r = %x, s = %x\\n\", q, r, s);\n\
    q &= 0xff;\n\
    r ^= 0xff;\n\
    s |= 0xff;\n\
\n\
    s = s + 3;\n\
    q = q + 3;\n\
printf(\"after q = %x, r = %x, s = %x\\n\", q, r, s);\n\
    return q + r + s;\n\
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

	GEN_PARSE_CONTEXT(context);
	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	ret = func(EC_param0);
	assert(ret == 0x8c + 0x6776 + 0x6802);
	ret = func(EC_param0);
	assert(ret == 0x8f + 0x6789 + 0x6902);
	ret = func(EC_param0);
	assert(ret == 0x92 + 0x6776 + 0x6a02);
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }


    return 0;
}

