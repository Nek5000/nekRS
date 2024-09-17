/* 
 * t8 - test the handling of pointers and basic pointer arithmetic.
 */
#include "config.h"
#include "cod.h"
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

static int verbose = 0;

int
main(int argc, char** argv)
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
	/* test pointer basics */
	char code_string[] = "\
{\n\
    int j = 2;\n\
    int * k;\n\
    k = &j;\n\
    *k = 4;\n\
    return (*k + 3);\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == 7);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 1 */
	/* test cast of complex type (&j) to basic type (intptr_t) */
	char code_string[] = "\
{\n\
    int j = 2;\n\
    long addr;\n\
    int * k;\n\
    addr = (long)&j;\n\
    k = (unsigned long)addr;\n\
    return (*k + 3 + j);\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	if (result != 7) {
	    printf("Expected 7, got %ld in t8, subtest 2\n", result);
	}
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 2 */
	/*
	 * Testing of pointer assignments and conversions:
	 * DR_P -> DR_L/DR_UL
	 * DR_C/DR_UC/DR_S/DR_US/DR_I/DR_U/DR_L/DR_UL -> DR_P
	 */
	char code_string[] = "\
{\n\
    char c = 1;\n\
    unsigned char uc = 2;\n\
    short s = 3;\n\
    unsigned short us = 4;\n\
    int i = 5;\n\
    unsigned int ui = 6;\n\
    long l = 7;\n\
    unsigned long ul = 8;\n\
    \n\
    int * p;\n\
    p  = l;\n\
    p  = ul;\n\
    \n\
    c  = p;\n\
    uc = p;\n\
    s  = p;\n\
    us = p;\n\
    i  = p;\n\
    ui = p;\n\
    l  = p;\n\
    ul = p;\n\
    return ui;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result==8);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 3 */
	/*
	 * Test op_inc pointer arithmetic on pointers
	 */
	char code_string[] = "\
{\n\
    int * p = (long)0;\n\
    p++;\n\
    return p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("int *proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == sizeof(int));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 4 */
	/*
	 * Test op_dec pointer arithmetic on pointers
	 */
	char code_string[] = "\
{\n\
    int * p = (long)0;\n\
    p--;\n\
    return p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("int *proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == -(long) sizeof(int));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 5 */
	/*
	 * Test op_inc pointer arithmetic on pointers to pointers
	 */
	char code_string[] = "\
{\n\
    double ** p = (long)0;\n\
    p++;\n\
    return p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == sizeof(double *));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 6 */
	/*
	 * Test op_plus pointer arithmetic on pointers
	 */
	char code_string[] = "\
{\n\
    int i;\n\
    double * p = (long)0;\n\
    double ** dp = &p;\n\
    p = (*dp) + 2;\n\
    return p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == 2*sizeof(double));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 7 */
	/*
	 * Test op_plus pointer arithmetic on pointers to pointers
	 */
	char code_string[] = "\
{\n\
    int i;\n\
    double ** dp = (long)0;\n\
    dp = dp + 2;\n\
    return dp;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == 2*sizeof(double *));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 8 */
	/*
	 * Test op_minus pointer arithmetic between pointers
	 */
	char code_string[] = "\
{\n\
     long r;\n\
    double * dp1 = (long)24;\n\
    double * dp2 = (long)8;\n\
    r = dp1 - dp2;\n\
    return r;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
        assert(result == 2);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 9 */
	/*
	 * Test op_minus pointer arithmetic between pointers
	 */
	char code_string[] = "\
{\n\
     long r;\n\
    double * dp1 = (long)24;\n\
    double * dp2 = (long)8;\n\
    r = dp2 - dp1;\n\
    return r;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("long proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
        assert(result == -2);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 10 */
	/*
	 * Test op_minus pointer arithmetic between pointers to pointers
	 */
	char code_string[] = "\
{\n\
     long r;\n\
    double ** dp1 = (long)24;\n\
    double ** dp2 = (long)8;\n\
    r = dp1 - dp2;\n\
    return r;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("long proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	if (result != 16/sizeof(void*)) { 
	    printf(" op minus result is %ld\n", result);
	    exit(1);
	}
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 11 */
	/*
	 * Test op_minus pointer arithmetic between pointers to pointers
	 */
	char code_string[] = "\
{\n\
     long r;\n\
    double ** dp1 = (long)24;\n\
    double ** dp2 = (long)8;\n\
    r = dp2 - dp1;\n\
    return r;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("long proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
        if (result != -(16/(long)sizeof(void*))) {
	    printf(" 2nd op minus result is %ld\n", result);
	    exit(1);
	}
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 12 */
	/*
	 * Test op_minus pointer arithmetic on a pointer and an integral
	 */
	char code_string[] = "\
{\n\
    int i;\n\
    double * p = (long)0;\n\
    double ** dp = &p;\n\
    p = (*dp) - 2;\n\
    return p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("int *proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == -2*(long)sizeof(double));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 13 */
	/*
	 * Test op_minus pointer arithmetic on a pointer to pointer and an integral
	 */
	char code_string[] = "\
{\n\
    int i;\n\
    double ** dp = (long)0;\n\
    dp = dp - 2;\n\
    return dp;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	long (*func)();
	long result;

	cod_subroutine_declaration("int *proc()", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (long(*)()) (intptr_t) gen_code->func;
	result = func();
	assert(result == -2*(long)sizeof(double *));
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 14 */
	/*
	 * Test pointer casting
	 */
	char code_string[] = "\
{\n\
    int i;\n\
    unsigned char *p;\n\
    p = param;\n\
    return (long) p;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	unsigned char *p = (unsigned char *) "hiya";
	intptr_t (*func)(unsigned char *);
	intptr_t result;

	cod_subroutine_declaration("long test(unsigned char * param)", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (intptr_t(*)(unsigned char *)) (intptr_t) gen_code->func;
	result = func(p);
	assert(result == (intptr_t) p);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 15 */
	/*
	 * Test postincrement to a pointer dereference
	 */
	char code_string[] = "\
{\n\
    return (*param)++;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	int i = 8;
	int (*func)(int *);
	int result;

	cod_subroutine_declaration("int test(int* param)", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (int(*)(int *)) (intptr_t) gen_code->func;
	result = func(&i);
	if (result != 8) printf("result was %d, not 8\n", result);
	if (i != 9) printf("i was %d, not 9\n", i);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    test_num++;
    if ((run_only == -1) || (run_only == test_num)) {  /* 16 */
	/*
	 * Test preincrement to a pointer dereference
	 */
	char code_string[] = "\
{\n\
    return ++(*param);\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	int i = 8;
	int (*func)(int *);
	int result;

	cod_subroutine_declaration("int test(int* param)", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (int(*)(int *)) (intptr_t) gen_code->func;
	result = func(&i);
	if (result != 9) printf("result was %d, not 8\n", result);
	if (i != 9) printf("i was %d, not 9\n", i);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    {
	/*
	 * Test preincrement to a parameter pointer
	 */
	char code_string[] = "\
{\n\
    return ++param;\n\
}";

	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	int i = 8;
	int *(*func)(int *);
	int *result;

	cod_subroutine_declaration("int test(int* param)", context);
	gen_code = cod_code_gen(code_string, context);
	if(!gen_code) {
	  printf("Code generation failed!\n");
	  return -1;
	} 
	func = (int*(*)(int *)) (intptr_t) gen_code->func;
	result = func(&i);
	if (result != (&i + 1)) printf("result was %p, not %p\n", result,
		&i + 1);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }

    return 0;
}
