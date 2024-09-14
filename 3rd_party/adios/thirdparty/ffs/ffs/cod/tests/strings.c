#include "config.h"
#include "cod.h"
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

int
main(int argc, char **argv)
{
    int test_to_run = -1;
    if (argc > 1) {
	sscanf(argv[1], "%d", &test_to_run);
    }

    if ((test_to_run == 1) || (test_to_run == -1)) {
	char code_string[] = "\
{\n\
	string a = strstr(input, \"soda\");\n\
	return a;\n\
}";

	static char extern_string[] = "int printf(string format, ...);\n\
void *malloc(int size);\n\
char *strstr(char* s1, char *s2);\n\
void free(void *pointer);\n";
	static cod_extern_entry externs[] = 
	{
	    {"printf", (void*)(intptr_t)printf},
	    {"malloc", (void*)(intptr_t)malloc},
	    {"strstr", (void*)(intptr_t)strstr},
	    {"free", (void*)(intptr_t)free},
	    {(void*)0, (void*)0}
	};
	cod_parse_context context = new_cod_parse_context();
	cod_code gen_code;
	char *(*func)(char*);
	char*result;

	cod_assoc_externs(context, externs);
	cod_parse_for_context(extern_string, context);
	cod_subroutine_declaration("char *proc(char *input)", context);
	gen_code = cod_code_gen(code_string, context);
	if (gen_code == NULL) {
	    printf("Code generation failed for test 2\n");
	} else {
	    func = (char*(*)(char*)) (intptr_t) gen_code->func;
	    result = func("I'd really like a cream soda please!");
	    assert(strcmp("soda please!", result) == 0);
	    cod_code_free(gen_code);
	}
	cod_free_parse_context(context);
    }

    return 0;
}

