#include "config.h"
#include "cod.h"
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include "malloc.h"
#endif

#include <string.h>
#include <stdio.h>
#include <stdint.h>
#ifndef HAVE_ATL_H
int
main()
{return 0;}
#else
#include "atl.h"

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

static int verbose = 0;

int
main(int argc, char **argv)
{
    char *test_only = NULL;
    int i = 1;
    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[i], "-test_only") == 0){
	    test_only = argv[++i];
	} else {
	    printf("Unknown argument %s\n", argv[i]);
	}
    }
    if (!test_only || (strcmp(test_only, "paramattr") == 0)) {
	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	/* test external call */
	static char code[] = "{\
			if (attr_set(l, \"test_value\")) {\n\
				return attr_ivalue(l, \"test_value\");\n\
			}\n\
			return 0;\n\
		}";

	cod_parse_context context = new_cod_parse_context();
	cod_parse_context context2;

	cod_code gen_code;
	int (*func)(attr_list l);
	attr_list l;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	context2 = cod_copy_context(context);

	l = create_attr_list();
	set_attr(l, attr_atom_from_string("test_value"), Attr_Int4, (attr_value)15);

	cod_subroutine_declaration("int proc(attr_list l)", context);
	gen_code = cod_code_gen(code, context);
	func = (int (*)(attr_list))(intptr_t)gen_code->func;

	if ((func)(l) != 15) {
	    printf("Function didn't return 15\n");
	}
	
	cod_code_free(gen_code);
	cod_free_parse_context(context);

	cod_subroutine_declaration("int proc(attr_list l)", context2);
	gen_code = cod_code_gen(code, context2);
	func = (int (*)(attr_list))(intptr_t)gen_code->func;

	if ((func)(l) != 15) {
	    printf("Function didn't return 15\n");
	}
	free_attr_list(l);
	cod_code_free(gen_code);
	cod_free_parse_context(context2);
    }
    if (!test_only || (strcmp(test_only, "globalattr") == 0)) {
	static char extern_string[] = "int printf(string format, ...);\n\
				       attr_list l;";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"l", (void*) 0},
	    {(void*)0, (void*)0}
	};
	/* test external call */
	static char code[] = "{\
			if (attr_set(l, \"test_value\")) {\n\
				return attr_ivalue(l, \"test_value\");\n\
			}\n\
			return 0;\n\
		}";

	cod_parse_context context = new_cod_parse_context();
	cod_parse_context context2;

	cod_code gen_code;
	int (*func)(attr_list l);
	attr_list l = create_attr_list();
	
	externs[1].extern_value = &l;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	context2 = cod_copy_context(context);

	set_attr(l, attr_atom_from_string("test_value"), Attr_Int4, (attr_value)15);

	cod_subroutine_declaration("int proc()", context);
	gen_code = cod_code_gen(code, context);
	func = (int (*)(attr_list))(intptr_t)gen_code->func;

	if ((func)(l) != 15) {
	    printf("Function didn't return 15\n");
	}
	
	free_attr_list(l);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
	cod_free_parse_context(context2);

    }
    if (!test_only || (strcmp(test_only, "attrarray") == 0)) {
	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	/* test external call */
	static char code[] = "{\
			if (attr_set(l[i], \"test_value\")) {\n\
				return attr_ivalue(l[i], \"test_value\");\n\
			}\n\
			return 0;\n\
		}";

	cod_parse_context context = new_cod_parse_context();

	cod_code gen_code;
	int (*func)(attr_list *l, int i);
	attr_list l[2];

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	cod_subroutine_declaration("int proc(attr_list l[], int i)", context);
	gen_code = cod_code_gen(code, context);
	func = (int (*)(attr_list*, int))(intptr_t)gen_code->func;

	l[0] = create_attr_list();
	set_attr(l[0], attr_atom_from_string("test_value"), Attr_Int4, (attr_value)15);
	l[1] = create_attr_list();
	set_attr(l[1], attr_atom_from_string("test_value"), Attr_Int4, (attr_value)25);

	if ((func)(&l[0], 0) != 15) {
	    printf("Function didn't return 15\n");
	}
	
	if ((func)(&l[0], 1) != 25) {
	    printf("Function didn't return 25\n");
	}
	
	free_attr_list(l[0]);
	free_attr_list(l[1]);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    if (!test_only || (strcmp(test_only, "timer") == 0)) {
	static char extern_string[] = "int printf(string format, ...);";

	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {(void*)0, (void*)0}
	};
	/* test external call */
	static char code[] = "{\
			chr_time timer;\n\
			int i;\n\
			int tmp = 0;\n\
			chr_timer_start(&timer);\n\
			for(i=0; i < limit; i++) {\n\
				tmp++;\n\
			}\n\
			chr_timer_stop(&timer);\n\
			return chr_time_to_millisecs(&timer);\n\
		}";

	cod_parse_context context = new_cod_parse_context();

	cod_code gen_code;
	double (*func)(int);
	double big, little;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);

	cod_subroutine_declaration("double proc(int limit)", context);
	gen_code = cod_code_gen(code, context);
	func = (double (*)(int))(intptr_t)gen_code->func;

	big =  (func)(1000000);
	little = (func)(10000);

	if (verbose) {
	    printf("200K autoincrements took %g milliseconds\n", big);
	    printf("2000 autoincrements took %g milliseconds\n", little);
	}

	if (((big / 100.0) < (.5 *little)) || 
	    ((big / 100.0) > (2.0 *little))) {
	    printf("chr_timer_test values inconsistent:  big, little out of range, %g, %g\n",
		   big, little);
	}
	
	cod_code_free(gen_code);
	cod_free_parse_context(context);
    }
    return 0;
}
#endif
