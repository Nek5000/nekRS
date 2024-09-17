#include "config.h"
#include "cod.h"
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>
char code_string[] = "{return sizeof(int);}";
char code_string2[] = "{return sizeof(int*);}";
char code_string3[] = "{   int j;      return sizeof j;}";

int
main()
{
    char *code_blocks[] = {code_string, code_string2, code_string3, NULL};
    int results[] = {sizeof(int), sizeof(int*), sizeof(int)};
    int test = 0;
	
    while (code_blocks[test] != NULL) {
	int ret;
	cod_parse_context context = new_cod_parse_context();
	cod_exec_context ec;

	cod_code gen_code;
	int (*func)(cod_exec_context);

	cod_add_param("ec", "cod_exec_context", 0, context);
	gen_code = cod_code_gen(code_blocks[test], context);
	ec = cod_create_exec_context(gen_code);
	func = (int(*)(cod_exec_context))gen_code->func;
	ret = (func)(ec);
	if (ret != results[test]) {
	    printf("bad test %d, ret was %d\n", test, ret);
	}
	cod_exec_context_free(ec);
	cod_code_free(gen_code);
	cod_free_parse_context(context);
	test++;
    }
    return 0;
}
