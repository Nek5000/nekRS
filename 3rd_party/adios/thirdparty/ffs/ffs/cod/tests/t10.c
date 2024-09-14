#include "config.h"
#include "cod.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#define assert(e)  \
    ((void) ((e) ? 0 : __assert (#e, __FILE__, __LINE__)))
#define __assert(e, file, line) \
    ((void)printf ("%s:%u: failed assertion `%s'\n", file, line, e), abort())

int
main(int argc, char **argv)
{
    int test_to_run = -1;
    if (argc > 1) {
	sscanf(argv[1], "%d", &test_to_run);
    }

    if ((test_to_run == 1) || (test_to_run == -1)) {
	/* test the basics */
	static char extern_string[] = "int printf(string format, ...);\n\
void *malloc(int size);\n\
void free(void *pointer);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    int ti,tj;\n\
	    long tk;\n\
	    short tl;\n\
	} test_struct;\n\
\n\
    int j = 4;\n\
    long k = 10;\n\
    short l = 3;\n\
    test_struct *a[2];\n\
    test_struct t;\n\
    t.tj = j;\n\
    t.tk = k;\n\
    t.tl = l;\n\
    a[0] = &t;\n\
    a[0].ti = j;\n\
    return a[0].tl * (t.tj + t.tk);\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == 42);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    if ((test_to_run == 2) || (test_to_run == -1)) {
	/* test the basics */
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    void *next;\n\
	    int i;\n\
	} test_struct;\n\
\n\
    int sum = 0;\n\
    test_struct *t1;\n\
    test_struct *t2;\n\
    t1 = malloc(sizeof(test_struct));\n\
    t1.i = 1; \n\
    t1.next = (void*)malloc(sizeof(test_struct));\n\
    t2 = (test_struct *) t1.next;\n\
    t2.i = 2;\n\
    t2.next = (void*)malloc(sizeof(test_struct));\n\
    t2 = (test_struct *) t2.next;\n\
    t2.i = 3;\n\
    t2.next = (void*)malloc(sizeof(test_struct));\n\
    t2 = (test_struct *) t2.next;\n\
    t2.i = 4;\n\
    t2.next = (void*)malloc(sizeof(test_struct));\n\
    t2 = (test_struct *) t2.next;\n\
    t2.i = 5;\n\
    t2.next = (void*)malloc(sizeof(test_struct));\n\
    t2 = (test_struct *) t2.next;\n\
    t2.i = 6;\n\
    t2.next = (void*)0;\n\
    while ( t1 != (void*)0 ) {\n\
	sum = sum + t1.i;\n\
	t2 = (test_struct *)t1.next;\n\
	free(t1); t1 = t2;\n\
    }\n\
    return sum;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n\
void *malloc(int size);\n\
void free(void *pointer);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 2\n");
	} else {
	    func = (long(*)()) (intptr_t) gen_code->func;
	    result = func();
	    assert(result == 21);
	    cod_code_free(gen_code);
	}
	cod_free_parse_context(context);
    }

    if ((test_to_run == 3) || (test_to_run == -1)) {
	/* test the array fields */
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    int size;\n\
	    int array[5];\n\
	} test_struct;\n\
\n\
    int i;\n\
    int sum = 0;\n\
    test_struct t1;\n\
    t1.size = 5;\n\
    t1.array[0] = 6;\n\
    t1.array[1] = 5;\n\
    t1.array[2] = 4;\n\
    t1.array[3] = 3;\n\
    t1.array[4] = 2;\n\
    for (i = 0; i < t1.size ; i++) {\n\
	sum = sum + t1.array[i];\n\
    }\n\
    return sum;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n\
void *malloc(int size);\n\
void free(void *pointer);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    func = (long(*)()) (intptr_t) gen_code->func;
	    result = func();
	    assert(result == 20);
	    cod_code_free(gen_code);
	}
	cod_free_parse_context(context);
    }
    if ((test_to_run == 4) || (test_to_run == -1)) {
	/* test dynamic array fields */
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    int size;\n\
	    int array[size];\n\
	} test_struct;\n\
\n\
    int i;\n\
    int sum = 0;\n\
    test_struct t1;\n\
    t1.size = 5;\n\
    t1.array[0] = 6;\n\
    t1.array[1] = 5;\n\
    t1.array[2] = 4;\n\
    t1.array[3] = 3;\n\
    t1.array[4] = 2;\n\
    for (i = 0; i < t1.size ; i++) {\n\
	sum = sum + t1.array[i];\n\
    }\n\
    return sum;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n\
void *malloc(int size);\n\
void free(void *pointer);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    func = (long(*)()) (intptr_t) gen_code->func;
	    result = func();
	    assert(result == 20);
	    cod_code_free(gen_code);
	}
	cod_free_parse_context(context);
    }
    if ((test_to_run == 5) || (test_to_run == -1)) {
	/* test dynamic array fields in static vars */
	char code_string[] = "\
{\n\
	typedef struct test {\n\
	    int size;\n\
	    int array[size];\n\
	} test_struct;\n\
\n\
    int i;\n\
    static int sum = 0;\n\
    static test_struct t1;\n\
    t1.size = 5;\n\
    t1.array[0] = 6;\n\
    t1.array[1] = 5;\n\
    t1.array[2] = 4;\n\
    t1.array[3] = 3;\n\
    t1.array[4] = 2;\n\
    for (i = 0; i < t1.size ; i++) {\n\
	sum = sum + t1.array[i];\n\
    }\n\
    return sum;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    func = (long(*)()) (intptr_t) gen_code->func;
	    result = func();
	    assert(result == 20);
	    cod_code_free(gen_code);
	}
	cod_free_parse_context(context);
    }
    if ((test_to_run == 6) || (test_to_run == -1)) {
typedef struct _dyn_arrays {
  int dim1;
  int dim2;
  int dim3;
  int *array1;
  int *array2;
  int *array3;
  int *array4;
  int *array5;
  int *array6;
  int *array7;
  int *array8;
} dyn_arrays,*dyn_arrays_t;

static  FMField dyn_arrays_field_list[] =
{
    {"dim1", "integer", sizeof(int), FMOffset(dyn_arrays_t, dim1)},
    {"dim2", "integer", sizeof(int), FMOffset(dyn_arrays_t, dim2)},
    {"dim3", "integer", sizeof(int), FMOffset(dyn_arrays_t, dim3)},
    {"array1", "*(float[7][5][3])", sizeof(int), FMOffset(dyn_arrays_t, array1)},
    {"array2", "float[dim1][5][3]", sizeof(int), FMOffset(dyn_arrays_t, array2)},
    {"array3", "float[7][dim2][3]", sizeof(int), FMOffset(dyn_arrays_t, array3)},
    {"array4", "float[dim2][dim2][3]", sizeof(int), FMOffset(dyn_arrays_t, array4)},
    {"array5", "float[7][5][dim3]", sizeof(int), FMOffset(dyn_arrays_t, array5)},
    {"array6", "float[dim1][5][dim3]", sizeof(int), FMOffset(dyn_arrays_t, array6)},
    {"array7", "float[7][dim2][dim3]", sizeof(int), FMOffset(dyn_arrays_t, array7)},
    {"array8", "float[dim2][dim2][dim3]", sizeof(int), FMOffset(dyn_arrays_t, array8)},
    {NULL, NULL, 0, 0}
};
	/* test dynamic array fields in static vars */
	char code_string[] = "\
{\n\
	int ret = 0;\n\
	int *p = &(input.array1[1][2][2]);\n		\
	if (p != &(input.array2[1][2][2])) printf(\"DIED 2\\n\"); else ret=ret+1;\n\
	if (p != &(input.array3[1][2][2])) printf(\"DIED 3\\n\"); else ret=ret+1;\n\
	if (p != &(input.array4[1][2][2])) printf(\"DIED 4\\n\"); else ret=ret+1;\n\
	if (p != &(input.array5[1][2][2])) printf(\"DIED 5\\n\"); else ret=ret+1;\n\
	if (p != &(input.array6[1][2][2])) printf(\"DIED 6\\n\"); else ret=ret+1;\n\
	if (p != &(input.array7[1][2][2])) printf(\"DIED 7\\n\"); else ret=ret+1;\n\
	if (p != &(input.array8[1][2][2])) printf(\"DIED 8\\n\"); else ret=ret+1;\n\
\n\
	return ret;\n\
}\n";

	static char extern_string[] = "int printf(string format, ...);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)(dyn_arrays *);
	long result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	cod_add_simple_struct_type("dyn_arrays", dyn_arrays_field_list, context);
	cod_subroutine_declaration("int proc(dyn_arrays *input)", context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 3\n");
	} else {
	    dyn_arrays  input;
	    func = (long(*)(dyn_arrays*)) (intptr_t) gen_code->func;
	    input.dim1 = 7;
	    input.dim2 = 5;
	    input.dim3 = 3;
	    input.array1 = malloc(input.dim1 * input.dim2 * input.dim3 * sizeof(input.array1[0]));
	    input.array2 = input.array1;
	    input.array3 = input.array1;
	    input.array4 = input.array1;
	    input.array5 = input.array1;
	    input.array6 = input.array1;
	    input.array7 = input.array1;
	    input.array8 = input.array1;
	    result = func(&input);
	    assert(result == 7);
	    cod_code_free(gen_code);
	    free(input.array1);
	}
	cod_free_parse_context(context);
    }
    return 0;
}

