#include "config.h"
#include "cod.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))

/*
 * This includes two MASSIVE global unsigned char arrays containing
 * the original image (scaled down Carianna Nebulae image) and a
 * correctly greyscaled copy of it (for comparison purposes).
 *
 * The format of the PPM images is as follows:
 * First a signature "P6\n" for color PPMs, "P5\n" for grayscale PPMs.
 * Then "X Y\n" -- two ASCII representations of X and Y size in pixels.
 * Then "DEPTH\n" -- ASCII representation of the intensity maxvalue.
 * Then a series of (binary!) RGB/Gray values. P6 pixel is thus 3 bytes,
 * and a grayscale pixel is thus 1 byte.
 *
 * E.g., this is how the header needs to be converted in our case:
 *
 * P6\n             P5\n
 * 57 57\n    =>    57 57\n
 * 255\n            255\n
 *
 * The rest of the image gets shortened to 1/3 of its size due to
 * compressing R+G+B bytes into one gray value byte.
 */
#include "gray-images.c"

/*
 * Buffer for the output image.
 * 10000 ought to be enough for all eventualities.
 * Input image is 9760 bytes long.
 * Grayscale result should be 3262 bytes long. 
 */
unsigned char output[10000];

/*
 * Global state for the subroutine to use.
 * [0]     = pixel fraction,
 * [1],[2] = up to two bytes of the pixel's RGB info. in order.
 */
unsigned char global[3];

/*
 * This should be as long or longer that the PPM header.
 */
#define HEADER       16

void doOneTest(void (*f)(), int step)
{
	long result;
        long offset;

	unsigned char * p, * q;

	long (*func)(unsigned char *,
		     unsigned char *,
		     unsigned long,
		     unsigned long,
		     unsigned char *);

	/*
	 * Convert the pointer we're passed to correct type
	 */
	func = (long (*)(unsigned char *, unsigned char *,
			 unsigned long,	  unsigned long,
			 unsigned char *)) (intptr_t) f;
		
	/* pre-clear the output buffer */
	for(p=output; p < output + sizeof(output); ) *p++ = 0x0;


	/* run a first chuck with the header and a pixel */
	global[0] = 0;
        result =  func(image, output, 0, HEADER, global);
        offset = HEADER;

	/* run the rest in chunks of some size */
        for(; offset + step  <= sizeof(image);) {
       	    result += func(&image[offset], &output[result],
			   offset, step, global);
	    offset += step;
        }

	/* run any leftover smaller than the chunksize */
        if(offset < sizeof(image)) {
	    result += func(&image[offset], &output[result],
		           offset, sizeof(image)-offset, global);
	}

	/* first, the easy size check */
	if (result != sizeof(gray)) { printf("result is %ld, not %zu\n", result,
					     sizeof(gray));}
	assert(result == sizeof(gray));

        /* then compare the whole image w/ expected gray */
        for(p=output, q=gray; q < (gray + sizeof(gray)); ) {
	    if(*p++ != *q++) break;
	}
        assert(q == (gray + sizeof(gray)));

	/* If we return then all is good */
}

static int verbose = 0;

int
main(int argc, char **argv)
{
/*
 * In this code header processing works, as well as the
 * main processing cycle (works fine for size == N * 3)
 */
    char code_string[] = "\
{\n\
   unsigned char * p;\n\
   unsigned char * r;\n\
\n\
   /* init pointers */\n\
   p = (unsigned char *)src;\n\
   r = (unsigned char *)dst;\n\
\n\
   if(offset == 0) {\n\
     /* process signature */\n\
     if(*p++ != 0x50) return -1;\n\
     else *r++ = 0x50;\n\
\n\
     if(*p++ != 0x36) return -1;\n\
     else *r++ = 0x35;\n\
\n\
     for(;*p != 0xa;) *r++ = *p++;\n\
     *r++ = *p++;\n\
\n\
     /* x/y size */\n\
     for(;*p != 0xa;) *r++ = *p++;\n\
     *r++ = *p++;\n\
\n\
     /* colordepth */\n\
     for(;*p != 0xa;) *r++ = *p++;\n\
     *r++ = *p++;\n\
   }\n\
\n\
   /* SO FAR SO GOOD */\n\
\n\
   /* deal with previous fragmented pixel */\n\
   if(*global == 2) {\n\
     *r++ = (*(global+1) >> 2) + (5 * (*(global+2) >> 3)) + (*p >> 3);\n\
     p++;\n\
   }\n\
   if(*global == 1) {\n\
     *r++ = (*(global+1) >> 2) + (5 * (*p >> 3)) + (*(p+1) >> 3);\n\
     p = p + 2;\n\
   }\n\
\n\
   /* main loop running through all the following full pixels */\n\
   for(;p <= (src + size - 3);) {\n\
     *r++ = (*p >> 2) + (5 * (*(p+1) >> 3)) + (*(p+2) >> 3);\n\
     p = p + 3;\n\
   }\n\
\n\
   /* deal with leftover fragment pixels */\n\
   *global = 0;\n\
   if((p - src) == size - 2) {\n\
     *global = 2;\n\
     *(global + 1) = *p;\n\
     *(global + 2) = *(p+1);\n\
     p = p + 2;\n\
   }\n\
   if((p - src) == size - 1) {\n\
     *global = 1;\n\
     *(global + 1) = *p;\n\
     p = p + 2;\n\
   }\n\
\n\
   return (long)(r - dst);\n\
}";

    static char extern_string[] = "int printf(string format, ...);";
    static cod_extern_entry externs[] = 
	{
            {"printf", (void*)(intptr_t)printf},
            {(void*)0, (void*)0}
        };
    
    cod_parse_context context = new_cod_parse_context();
    cod_code gen_code;

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);
    
    cod_subroutine_declaration("long proc("
			       "unsigned char * src,"
			       "unsigned char * dst,"
			       "unsigned long offset,"
			       "unsigned long size,"
			       "unsigned char * global)", context);
    gen_code = cod_code_gen(code_string, context);
    if(!gen_code) {
	printf("Code generation failed!\n");
	return -1;
    } 


    /*
     * Run a tests with different fragment sizes
     */

    /*
     * This bunch will succeed because they all use fragment size
     * which is a multiple of 3 bytes (the pixel size) and do not
     * exercise the code which has problems...
     */
    
    if (argc > 1) verbose++;
    if (verbose) printf("running test with fragment size 3\n");
    doOneTest(gen_code->func, 3);
    if (verbose) printf("running test with fragment size 6\n");
    doOneTest(gen_code->func, 6);
    if (verbose) printf("running test with fragment size 18\n");
    doOneTest(gen_code->func, 18);
    if (verbose) printf("running test with fragment size 60\n");
    doOneTest(gen_code->func, 60);
    if (verbose) printf("running test with fragment size 180\n");
    doOneTest(gen_code->func, 180);
    if (verbose) printf("running test with fragment size 1002\n");
    doOneTest(gen_code->func, 1002);
    if (verbose) printf("running test with fragment size 2004\n");
    doOneTest(gen_code->func, 2004);
    if (verbose) printf("running test with fragment size 3006\n");
    doOneTest(gen_code->func, 3006);

    /*
     * These exercise the fragmented pixel code
     */

    if (verbose) printf("running test with fragment size 1000\n");
    doOneTest(gen_code->func, 1000);    /* 1000 % 3 = 1 */
    if (verbose) printf("running test with fragment size 2000\n");
    doOneTest(gen_code->func, 2000);    /* 2000 % 3 = 2 */
    if (verbose) printf("running test with fragment size 4\n");
    doOneTest(gen_code->func, 4);       /* 4 % 3 = 1    */
    if (verbose) printf("running test with fragment size 5\n");
    doOneTest(gen_code->func, 5);       /* 5 % 3 = 2    */
    if (verbose) printf("running test with fragment size 7\n");
    doOneTest(gen_code->func, 7);       /* 7 % 3 = 1    */
    if (verbose) printf("running test with fragment size 8\n");
    doOneTest(gen_code->func, 8);       /* 8 % 3 = 2    */

    cod_code_free(gen_code);
    cod_free_parse_context(context);
    
    return 0;
}
