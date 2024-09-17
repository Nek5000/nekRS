/*
 *   cod - T2
 *     
 *       Printf working?   Float arrays?  dynamic arrays?
 */
#include "config.h"
#include "data_funcs.h"
#include "cod.h"
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include "malloc.h"
#endif
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>
#include <string.h>

int
main(int argc, char **argv)
{
    int test_num = 0;
    int run_only = -1;
    char *read_file = NULL;
    char *write_file = NULL;
    int verbose = 0;
    while (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[1], "-w") == 0) {
	    if (argc <= 1) {
		printf("Need argument to \"-w\"\n");
	    } else {
		write_file = strdup(argv[2]);
	    }
	    argc--; argv++;
	} else if (strcmp(argv[1], "-r") == 0) {
	    if (argc <= 1) {
		printf("Need argument to \"-r\"\n");
	    } else {
		read_file = strdup(argv[2]);
	    }
	    argc--; argv++;
	} else if (strcmp(argv[1], "-o") == 0) {
	    sscanf(argv[2], "%d", &run_only);
	    argc--; argv++;
	}
	argc--; argv++;
    }
    if ((run_only == -1) || (run_only == test_num)) {
	/* 0 */
	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	/* test external call */
	static char code[] = "{\
			printf(\"values are is %d, %g, %s\\n\", i, d, s);\
		}";

	cod_parse_context context = new_cod_parse_context();

	cod_code gen_code;
	void (*func)(int, double, char*);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	cod_subroutine_declaration("void proc(int i, double d, string s)", 
				   context);
	gen_code = cod_code_gen(code, context);
	func = (void (*)(int, double, char*))(intptr_t)gen_code->func;
	printf("Expect -> \"values are is %d, %g, %s\"\n", 5, 3.14159, "hello!");
	(func)(5, (double)3.14159, "hello!");
	
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* 1 */
	typedef struct test {
	    int count;
	    double *vals;
	} test_struct, *test_struct_p;

	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
		    int i;\n\
		    double sum = 0.0;\n\
		    for(i = 0; i<input.count; i= i+1) {\n\
			sum = sum + input.vals[i];\n\
		    }\n\
		    return sum;\n\
		}";

	static FMField input_field_list[] =
	{
	    {"count", "integer", sizeof(int), 
	     FMOffset(test_struct_p, count)},
	    {"vals", "float[count]", sizeof(double), 
	     FMOffset(test_struct_p, vals)},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	int i;
	test_struct tmp;
	test_struct *param = &tmp;
	cod_code gen_code;
	double (*func)(test_struct_p);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	if (read_file) {
	    FMContext c = create_local_FMcontext();
	    char *buf = read_buffer(c, read_file, test_num);
	    param = (test_struct*)buf;
	    cod_add_encoded_param("input", buf, 0, c, context);
	    cod_set_return_type("double", context);
	} else {
	    cod_add_simple_struct_type("input_type", input_field_list, context);
	    cod_subroutine_declaration("double proc(input_type *input)", context);
	}
	tmp.count = 10;
	tmp.vals = (double*) malloc(tmp.count * sizeof(double));
	for(i=0; i< tmp.count; i++) {
	    tmp.vals[i] = i + 0.1;
	}

	if (write_file) {
	    FMStructDescRec formats[] = {
		{"struct", input_field_list, sizeof(tmp), NULL},
		{NULL, NULL, 0, NULL}};
	    write_buffer(write_file, &formats[0], &tmp, test_num);
	}
	gen_code = cod_code_gen(code, context);
	func = (double (*)(test_struct_p))(intptr_t) gen_code->func;
	assert((func)(param) == 46.00);
	free(tmp.vals);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* 2 */
	typedef struct test {
	    int count;
	    int *vals;
	} test_struct, *test_struct_p;

	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
		    int i;\n\
		    double sum = 0.0;\n\
		    for(i = 0; i<input.count; i= i+1) {\n\
			sum = sum + input.vals[i];\n\
		    }\n\
		    return sum;\n\
/* comment */\n\
		}";

	static FMField input_field_list[] =
	{
	    {"count", "integer", sizeof(int), 
	     FMOffset(test_struct_p, count)},
	    {"vals", "integer[count]", sizeof(int), 
	     FMOffset(test_struct_p, vals)},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	int i;
	test_struct tmp;
	test_struct *param = &tmp;
	cod_code gen_code;
	int (*func)(test_struct_p);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	if (read_file) {
	    FMContext c = create_local_FMcontext();
	    char *buf = read_buffer(c, read_file, test_num);
	    param = (test_struct*)buf;
	    cod_add_encoded_param("input", buf, 0, c, context);
	    cod_set_return_type("int", context);
	} else {
	    cod_add_simple_struct_type("input_type", input_field_list, context);
	    cod_subroutine_declaration("int proc(input_type *input)", context);
	}
	tmp.count = 10;
	tmp.vals = (int*) malloc(tmp.count * sizeof(int));
	for(i=0; i< tmp.count; i++) {
	    tmp.vals[i] = i + 10;
	}

	if (write_file) {
	    FMStructDescRec formats[] = 
		{{"struct", input_field_list, sizeof(tmp), NULL},
		 {NULL, NULL, 0, NULL}};
	    write_buffer(write_file, &formats[0], &tmp, test_num);
	}
	gen_code = cod_code_gen(code, context);
	func = (int (*)(test_struct_p))(intptr_t) gen_code->func;
	assert((func)(param) == 145);
	free(tmp.vals);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* 3 */
	typedef struct test {
	    double *vals;
	} test_struct, *test_struct_p;

	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
		    int i;\n\
		    double sum = 0.0;\n\
		    int count = 10;\n\
		    for(i = 0; i<count; i= i+1) {\n\
			sum = sum + input.vals[i];\n\
		    }\n\
		    return sum;\n\
		}";

	static FMField input_field_list[] =
	{
	    {"vals", "*float[10]", sizeof(double), 
	     FMOffset(test_struct_p, vals)},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	int i, count = 10;
	test_struct tmp;
	test_struct *param = &tmp;
	cod_code gen_code;
	double (*func)(test_struct_p);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	if (read_file) {
	    FMContext c = create_local_FMcontext();
	    char *buf = read_buffer(c, read_file, test_num);
	    param = (test_struct*)buf;
	    cod_add_encoded_param("input", buf, 0, c, context);
	    cod_set_return_type("double", context);
	} else {
	    cod_add_simple_struct_type("input_type", input_field_list, context);
	    cod_subroutine_declaration("double proc(input_type *input)", context);
	}
	count = 10;
	tmp.vals = (double*) malloc(count * sizeof(double));
	for(i=0; i< count; i++) {
	    tmp.vals[i] = i + 0.1;
	}

	if (write_file) {
	    FMStructDescRec formats[] = {
		{"struct", input_field_list, sizeof(tmp), NULL},
		{NULL, NULL, 0, NULL}};
	    write_buffer(write_file, &formats[0], &tmp, test_num);
	}
	gen_code = cod_code_gen(code, context);
	func = (double (*)(test_struct_p))(intptr_t) gen_code->func;
	assert((func)(param) == 46.00);
	free(tmp.vals);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* 4 */
	typedef struct test {
	    int count;
	    double *vals;
	} test_struct, *test_struct_p;

	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
		    int i;\n\
		    double sum = 0.0;\n\
		    for(i = 0; i<input.count; i= i+1) {\n\
		        double *ptr = input.vals + i;\n			 \
			sum = sum + (*ptr);\n\
		    }\n\
		    return sum;\n\
		}";

	static FMField input_field_list[] =
	{
	    {"count", "integer", sizeof(int), 
	     FMOffset(test_struct_p, count)},
	    {"vals", "float[count]", sizeof(double), 
	     FMOffset(test_struct_p, vals)},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	int i;
	test_struct tmp;
	test_struct *param = &tmp;
	cod_code gen_code;
	double (*func)(test_struct_p);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	if (read_file) {
	    FMContext c = create_local_FMcontext();
	    char *buf = read_buffer(c, read_file, test_num);
	    param = (test_struct*)buf;
	    cod_add_encoded_param("input", buf, 0, c, context);
	    cod_set_return_type("double", context);
	} else {
	    cod_add_simple_struct_type("input_type", input_field_list, context);
	    cod_subroutine_declaration("double proc(input_type *input)", context);
	}
	tmp.count = 10;
	tmp.vals = (double*) malloc(tmp.count * sizeof(double));
	for(i=0; i< tmp.count; i++) {
	    tmp.vals[i] = i + 0.1;
	}

	if (write_file) {
	    FMStructDescRec formats[] = {
		{"struct", input_field_list, sizeof(tmp), NULL},
		{NULL, NULL, 0, NULL}};
	    write_buffer(write_file, &formats[0], &tmp, test_num);
	}
	gen_code = cod_code_gen(code, context);
	func = (double (*)(test_struct_p))(intptr_t) gen_code->func;
	assert((func)(param) == 46.00);
	free(tmp.vals);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* 5 */
	/* 
	 *  This test tries to see if we're correctly dereferencing static pointers.  Earlier this would overwrite the 'testing' variable, indexing through
	 *  static memory instead of the malloc'd block.
	 */

	typedef struct test {
	    int count;
	    double *vals;
	} test_struct, *test_struct_p;

	static char extern_string[] = "int printf(string format, ...); void*malloc(int size);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {(void*)0, (void*)0}
	};
	static char code[] = "{\
		    static float *ptr;\n\
		    static float testing = 46;\n\
ptr = malloc(3*sizeof(float));\n\
ptr[0] = 5;\n\
ptr[1] = 6;\n\
ptr[2] = 7;\n\
return testing;\n\
		}";

	static FMField input_field_list[] =
	{
	    {"count", "integer", sizeof(int), 
	     FMOffset(test_struct_p, count)},
	    {"vals", "float[count]", sizeof(double), 
	     FMOffset(test_struct_p, vals)},
	    {(void*)0, (void*)0, 0, 0}
	};

	cod_parse_context context = new_cod_parse_context();
	int i;
	test_struct tmp;
	test_struct *param = &tmp;
	cod_code gen_code;
	double (*func)(test_struct_p);

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	if (read_file) {
	    FMContext c = create_local_FMcontext();
	    char *buf = read_buffer(c, read_file, test_num);
	    param = (test_struct*)buf;
	    cod_add_encoded_param("input", buf, 0, c, context);
	    cod_set_return_type("double", context);
	} else {
	    cod_add_simple_struct_type("input_type", input_field_list, context);
	    cod_subroutine_declaration("double proc(input_type *input)", context);
	}
	tmp.count = 10;
	tmp.vals = (double*) malloc(tmp.count * sizeof(double));
	for(i=0; i< tmp.count; i++) {
	    tmp.vals[i] = i + 0.1;
	}

	if (write_file) {
	    FMStructDescRec formats[] = {
		{"struct", input_field_list, sizeof(tmp), NULL},
		{NULL, NULL, 0, NULL}};
	    write_buffer(write_file, &formats[0], &tmp, test_num);
	}
	gen_code = cod_code_gen(code, context);
	func = (double (*)(test_struct_p))(intptr_t) gen_code->func;
	assert((func)(param) == 46.00);
	(void)func;
	free(tmp.vals);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    return 0;
}
