/*
 * t7 - not actually a regression test.
 *   subtest 0 - this code was created to try to verify the handling of different types of 
 *       static and dynamic arrays 
 *       Added a bit of code that dumps out types that are passed into routines like submit.
 *   subtest 1 - this code really has to be run in verbose mode to see its effects.  Created to debug
 *       "const int blah = 5;" sorts of things.
 *
 *     no success or failure is associated with the execution of this file.
 */
#include "cod.h"
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include "malloc.h"
#endif
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>
#include <string.h>

static int verbose = 0;
#ifdef NO_EMULATION
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();
#define EC_param0
#define EC_param1
#define EC_param0_decl
#define EC_param1_decl
#else
#define GEN_PARSE_CONTEXT(x) \
x = new_cod_parse_context();\
cod_add_param("ec", "cod_exec_context", 0, x);
#define EC_param0 ec
#define EC_param1 ec,
#define EC_param0_decl cod_exec_context
#define EC_param1_decl cod_exec_context,
#endif

typedef struct _multi_array_rec {
    long	ifield;
    double	double_array[2][2][2][2];
    int		(*int_array)[2];
    int		(*int_array2)[4];
    int		(*int_array3)[4][4];
} multi_array_rec, *multi_array_rec_ptr;

FMField multi_array_flds[] = {
    {"ifield", "integer", sizeof(int), FMOffset(multi_array_rec_ptr, ifield)},
    {"double_array", "float[2][2][2][2]", sizeof(double),
     FMOffset(multi_array_rec_ptr, double_array)},
    {"int_array", "integer[ifield][2]", sizeof(int),
     FMOffset(multi_array_rec_ptr, int_array)},
    {"int_array2", "integer[2][ifield]", sizeof(int),
     FMOffset(multi_array_rec_ptr, int_array2)},
    {"int_array3", "integer[ifield][ifield][ifield]", sizeof(int),
    FMOffset(multi_array_rec_ptr, int_array3)},
    {(char *) 0, (char *) 0, 0, 0}
};


typedef struct _complex_rec {
    double r;
    double i;
} complex, *complex_ptr;

typedef struct _nested_rec {
    complex item;
} nested, *nested_ptr;

static FMField nested_field_list[] =
{
    {"item", "complex", sizeof(complex), FMOffset(nested_ptr, item)},
    {NULL, NULL, 0, 0}
};

static FMField complex_field_list[] =
{
    {"r", "double", sizeof(double), FMOffset(complex_ptr, r)},
    {"i", "double", sizeof(double), FMOffset(complex_ptr, i)},
    {NULL, NULL, 0, 0}
};

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;

static FMField simple_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {"short_field", "integer",
     sizeof(short), FMOffset(simple_rec_ptr, short_field)},
    {"long_field", "integer",
     sizeof(long), FMOffset(simple_rec_ptr, long_field)},
    {"nested_field", "nested",
     sizeof(nested), FMOffset(simple_rec_ptr, nested_field)},
    {"double_field", "float",
     sizeof(double), FMOffset(simple_rec_ptr, double_field)},
    {"char_field", "char",
     sizeof(char), FMOffset(simple_rec_ptr, char_field)},
    {"scan_sum", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, scan_sum)},
    {NULL, NULL, 0, 0}
};

int
submit(cod_exec_context ec, int port, void *data, void *type_info)
{
    FMStructDescList formats = (FMStructDescList) type_info;
    printf("In submit, ec is %p, port is %d, data is %p, type_info is %p\n",
	   ec, port, data, type_info);
    while (formats[0].format_name != NULL) {
	FMFieldList tmp = formats[0].field_list;
	printf("Format \"%s\" - \n", formats[0].format_name);
	while(tmp[0].field_name != NULL) {
	    printf("{%s, %s, %d, %d},\n", tmp[0].field_name,
		   tmp[0].field_type, tmp[0].field_size, tmp[0].field_offset);
	    tmp++;
	}
	formats++;
    }
    return 1;
}

int
main(int argc, char **argv) 
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
	static char extern_string[] = "int printf(string format, ...);\n\
	int submit(cod_exec_context ec, int port, void *d, cod_type_spec dt);";
	static cod_extern_entry externs[] =
	    {
		{"printf", (void*)(intptr_t)printf},
		{"submit", (void*)(intptr_t)submit},
		{(void*)0, (void*)0}
	    };
	static char code[] = "{\n\
    submit(1, input);\n		 \
     }";
	
/*

*/
    int i, j, k, l;
    cod_code gen_code;
    cod_exec_context ec;
    void (*func)(cod_exec_context, void*);

    multi_array_rec multi_array;
   
    cod_parse_context context = new_cod_parse_context();

    multi_array.ifield = 4;
    for (i = 0; i < 2; i++) {
	for (j = 0; j < 2; j++) {
	    for (k = 0; k < 2; k++) {
		for (l = 0; l < 2; l++) {
		    multi_array.double_array[i][j][k][l] = 
			1000*i + 100*j + 10*k +l;
		}
	    }
	}
    }
    multi_array.int_array = malloc(2*4*sizeof(int));
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 2; j++) {
	    multi_array.int_array[i][j] = 1000*i + 100*j;
	}
    }
    multi_array.int_array2 = malloc(4*2*sizeof(int));
    for (i = 0; i < 2; i++) {
	for (j = 0; j < 4; j++) {
	    printf("element [%d][%d] is at addr %p\n", i, j, &multi_array.int_array2[i][j]);
	    multi_array.int_array2[i][j] = 1000*i + 100*j;
	}
    }
    multi_array.int_array3 = malloc(4*4*4*sizeof(int));
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    for (k = 0; k < 4; k++) {
		multi_array.int_array3[i][j][k] = 1000*i + 100*j + 10*k;
	    }
	}
    }


    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);

    cod_add_simple_struct_type("multi_array", multi_array_flds, context);
    cod_add_simple_struct_type("nested", nested_field_list, context);
    cod_add_simple_struct_type("complex", complex_field_list, context);
    cod_add_simple_struct_type("simple", simple_field_list, context);
    cod_subroutine_declaration("void proc(cod_exec_context ec, simple *input)", context);
   
    gen_code = cod_code_gen(code, context);
    func = (void (*)(cod_exec_context, void*))(intptr_t)gen_code->func;

    cod_dump(gen_code);
    ec = cod_create_exec_context(gen_code);
    printf("Main ec is %p\n", ec);
    func(ec, &multi_array);
    cod_exec_context_free(ec);
    cod_code_free(gen_code);
    cod_free_parse_context(context);
/*    if ((data.num_points != 1) || (data.image_data[0].num_points != 1)) 
      return 1;*/
    return 0;
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {
	/* test 1 */
	char code_string[] = "\
{\n\
    const int j = 4;\n				\
    const long k = 10;\n				\
    const short l = 3;\n					\
\n\
    return j + j;\n\
}";

	cod_parse_context context;
	cod_exec_context ec;
	cod_code gen_code;
	long (*func)(EC_param0_decl);
	long result;

	GEN_PARSE_CONTEXT(context);
	gen_code = cod_code_gen(code_string, context);
	ec = cod_create_exec_context(gen_code);
	func = (long(*)(EC_param0_decl)) (intptr_t) gen_code->func;
	if (verbose) cod_dump(gen_code);
	result = func(EC_param0);
	assert(result == 8);
	cod_code_free(gen_code);
	cod_exec_context_free(ec);
	cod_free_parse_context(context);
    }
    test_num++; /* 2 */
}
