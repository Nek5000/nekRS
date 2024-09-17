#include "config.h"
#include "cod.h"
#include <stdio.h>
#include <stdint.h>

#define IMGDEBUG

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

#if defined(IMGDEBUG)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

int fd;
#endif
static char mixImage[] = "\n\
{\n\
  int col;\n\
  int row;\n\
  int indx;\n\
  long longfactor;\n\
\n\
  longfactor = (long)(fade * 65536);\n\
  indx = 0;\n\
  for (row = 0; row < 100; row = row + 1) {\n\
    for (col = 0; col < 100; col = col + 1) {\n\
\n\
      output.out[indx] = input1.in[indx] + \n\
	(((input2.in[indx] - input1.in[indx]) \n\
	  * longfactor) >> 16);\n\
\n\
      output.out[indx+1] = input1.in[indx+1] + \n\
	(((input2.in[indx+1] - input1.in[indx+1]) \n\
	  * longfactor) >> 16);\n\
\n\
      output.out[indx+2] = input1.in[indx+2] + \n\
	(((input2.in[indx+2] - input1.in[indx+2]) \n\
	  * longfactor) >> 16);\n\
      \n\
      indx = indx + 3;\n\
    }\n\
  }\n\
\n\
  return 0;\n\
}";

static char extern_string[] = "int printf(string format, ...);";
static cod_extern_entry externs[] = 
{
    {"printf", (void*)(long)printf},
    {(void*)0, (void*)0}
};

/*
  int mixImage(unsigned char in1[100*100*3],
	       unsigned char in2[100*100*3],
	       unsigned char out[100*100*3],
	       float fade)
 */


int main()
{
    static FMField input1_field_list[] = {
	{"in", "unsigned integer[30000]", sizeof(unsigned char), 0},
	{(void*)0, (void*)0, 0, 0}
    };
  
    static FMField input2_field_list[] = {
	{"in", "unsigned integer[30000]", sizeof(unsigned char), 0},
	{(void*)0, (void*)0, 0, 0}
    };

    static FMField output_field_list[] = {
	{"out", "unsigned integer[30000]", sizeof(unsigned char), 0},
	{(void*)0, (void*)0, 0, 0}
    };
    int i, index;

    unsigned char out[100*100*3];
    unsigned char img1a[100*100*3];
    unsigned char img2a[100*100*3];

    cod_parse_context context = new_cod_parse_context();
    cod_code gen_code;
    long (*func)(unsigned char *, unsigned char *, unsigned char *, float);
    
    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);

    cod_add_simple_struct_type("input1_type", input1_field_list, context);
    cod_add_simple_struct_type("input2_type", input2_field_list, context);
    cod_add_simple_struct_type("output_type", output_field_list, context);
    cod_subroutine_declaration("void proc(input1_type *input1, input2_type *input2, output_type output, float fade)", context);


    gen_code = cod_code_gen(mixImage, context);
    func = (long (*)(unsigned char *, unsigned char *, unsigned char *, float)) (long) gen_code->func;
    
    if(!func) return 1;

    index = 0;
    for(i = 0; i < 100; i++) {
	int j;
	for(j = 0; j < 100; j++) {
	    int k;
	    for(k = 0; k < 3; k++) {
		img1a[index] = 0xfa;
		img2a[index] = 0x1a;
		index++;
	    }
	}
    }
    (void) func(img1a, img2a, out, 0.5);

    index = 0;
    for(i = 0; i < 100; i++) {
	int j;
	for(j = 0; j < 100; j++) {
	    int k;
	    for(k = 0; k < 3; k++) {
		if (out[index] != 0x8a) {
		    printf("Out element [%d,%d,%d] is incorrect -> 0x%x\n", 
			   i, j, k, out[index]);
		}
		index++;
	    }
	}
    }
    
    cod_code_free(gen_code);
    cod_free_parse_context(context);
    return 0;
}

