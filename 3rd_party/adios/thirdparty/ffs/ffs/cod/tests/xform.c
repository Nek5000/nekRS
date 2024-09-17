#include "config.h"
#include "ecl.h"
#define assert(EX) ((EX) ? (void)0 : (fprintf(stderr, "\"%s\" failed, file %s, line %d\n", #EX, __FILE__, __LINE__), exit(1)))
#include <stdio.h>
#include <stdint.h>

#define IMGDEBUG

#if defined(IMGDEBUG)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "img1.h"
#include "img2.h"

int fd;
#endif

static char extern_string[] = "int printf(string format, ...);";
static ecl_extern_entry externs[] = 
{
    {"printf", (void*)(long)printf},
    {(void*)0, (void*)0}
};

/*
  int catImage(char in1[100*100*3],
               char in2[100*100*3],
               char out[200*100*3])
*/
char catImage[] = "\n\
{\n\
  int col;\n\
  int row;\n\
  int indx1;\n\
  int indx2;\n\
  int outIndx;\n\
\n\
  indx1 = 0;\n\
  indx2 = 0;\n\
  outIndx = 0;\n\
\n\
  for (row = 0; row < 100; row = row + 1) {\n\
    for (col = 0; col < 100; col = col + 1) {\n\
      output.out[outIndx]   = input1.in[indx1];   /* red */\n\
      output.out[outIndx+1] = input1.in[indx1+1]; /* green */\n\
      output.out[outIndx+2] = input1.in[indx1+2]; /* blue */\n\
      outIndx = outIndx + 3;\n\
      indx1 = indx1 + 3;\n\
    }\n\
    for (col = 0; col < 100; col = col + 1) {\n\
      output.out[outIndx]   = input2.in[indx2];   /* red */\n\
      output.out[outIndx+1] = input2.in[indx2+1]; /* green */\n\
      output.out[outIndx+2] = input2.in[indx2+2]; /* blue */\n\
      outIndx = outIndx + 3;\n\
      indx2 = indx2 + 3;\n\
    }\n\
  }\n\
\n\
  return 0;\n\
}";


/*
  int cropImage(char in [100*100*3],
                char out[50*50*3],
                int inRow, int inCol)
 */
char cropImage[] = "\n\
{\n\
\n\
  int row;\n\
  int col;\n\
  int outIndx;\n\
  int inIndx;\n\
\n\
  for (row = 0; row < 50; row = row + 1) {\n\
    for (col = 0; col < 50; col = col + 1) {\n\
      outIndx = (row * 50 * 3) + (col * 3);\n\
      inIndx = ((inRow + row) * 100 * 3) + (inCol + col) * 3;\n\
      output.out[outIndx]   = input.in[inIndx];\n\
      output.out[outIndx+1] = input.in[inIndx+1];\n\
      output.out[outIndx+2] = input.in[inIndx+2];\n\
    }\n\
  }\n\
  \n\
  return 0;\n\
}";


/*
  int edgeImage(char in [100*100*3],
                char out[100*100*3])
 */
char edgeImage[] = "\n\
{\n\
  int col;\n\
  int row;\n\
  int indx;\n\
  int outIndx;\n\
  int row0Indx;\n\
  int row1Indx;\n\
  int row2Indx;\n\
  int tmprowIndx;\n\
\n\
  double sum1;\n\
  double sum2;\n\
  double sum;\n\
  double abssum1;\n\
  double abssum2;\n\
\n\
  /* Get the first two rows. */\n\
  row0Indx = 0;\n\
  row1Indx = 100;\n\
\n\
  /* Write out the first row, all zeros. */\n\
  for (col = 0; col < 100; col = col + 1)\n\
    output.out[col] = 0;\n\
\n\
  /* Now the rest of the image -- read in the next row, and write\n\
  ** write out the current row.\n\
  */\n\
  for (row = 1; row < (100 - 1); row = row + 1) {\n\
    indx = (row - 1) * 100;\n\
    outIndx = (row - 1) * 100;\n\
    row2Indx = indx;\n\
\n\
    for ( col = 1; col < 100 - 1; col = col + 1 ) {\n\
      sum1 = (double) input.in[row0Indx+col+1] - \n\
	(double) input.in[row0Indx+col-1] +\n\
	2.0 * ( (double) input.in[row1Indx+col+1] - \n\
		(double) input.in[row1Indx+col-1] ) +\n\
	(double) input.in[row2Indx+col+1] - \n\
	(double) input.in[row2Indx+col-1];\n\
\n\
      sum2 = ( (double) input.in[row2Indx+col-1] + \n\
	       2.0 * (double) input.in[row2Indx+col] +\n\
	       (double) input.in[row2Indx+col+1] ) -\n\
	( (double) input.in[row0Indx+col-1] + \n\
	  2.0 * (double) input.in[row0Indx+col] +\n\
	  (double) input.in[row0Indx+col+1] );\n\
      \n\
      /* \n\
       * The following 8 lines of code are just to avoid using sqrt since\n\
       * we cannot 'move' it into E-code.\n\
       * Now, \"why 0.2 of the smallest one?\" you would ask: because I say\n\
       * so.\n\
       *\n\
       * sum = sqrt( sum1 * sum1 + sum2 * sum2 ); \n\
       *\n\
       */\n\
\n\
      if (sum1 < 0.0) abssum1 = -sum1;\n\
      else abssum1 = sum1;\n\
      if (sum2 < 0.0) abssum2 = -sum2;\n\
      else abssum2 = sum2;\n\
\n\
      if (abssum1 == 0.0) sum = abssum2;\n\
      if (abssum2 == 0.0) sum = abssum1;\n\
\n\
      if (abssum1 > abssum2) sum = abssum1 + 0.2 * abssum2;\n\
      else sum = 0.2 * abssum1 + abssum2;\n\
\n\
      sum = sum / 1.8;		/* arbitrary scaling factor */\n\
      if ( sum > 255 ) sum = 255;\n\
      output.out[outIndx + col] = sum;\n\
    }\n\
    \n\
    /* Rotate rows. */\n\
    tmprowIndx = row0Indx;\n\
    row0Indx = row1Indx;\n\
    row1Indx = row2Indx;\n\
    row2Indx = tmprowIndx;\n\
  }\n\
  \n\
  /* And write the last row, zeros again. */\n\
  outIndx = (100 - 2) * 100;\n\
  for (col = 0; col < 100; col = col + 1)\n\
    output.out[outIndx + col] = 0;\n\
  \n\
  return 0;\n\
}";

/*
  int enlargeImage(char in [100*100*3],
       	           char out[200*200*3],
		   int  n = 2)
 */
char enlargeImage[] = "\n\
{\n\
  int col;\n\
  int subcol;\n\
  int row;\n\
  int subrow;\n\
\n\
  int rowIndx;\n\
  int outIndx;\n\
  int colIndx;\n\
  int outColIndx;\n\
  int outRowIndx;\n\
\n\
  outRowIndx = 0;\n\
  for (row = 0; row < 100; row = row + 1) {\n\
    rowIndx = row * 100 * 3;\n\
    for (subrow = 0; subrow < n; subrow = subrow + 1) {\n\
    col = 0;\n\
    colIndx = 0;\n\
      for (; col < 100;) {\n\
	outIndx = outRowIndx * 200 * 3;\n\
\n\
        subcol = 0;\n\
        outColIndx = (col * n * 3);\n\
	for (; subcol < n;) {\n\
	  output.out[outIndx+outColIndx+0] = input.in[rowIndx+colIndx+0];\n\
	  output.out[outIndx+outColIndx+1] = input.in[rowIndx+colIndx+1];\n\
	  output.out[outIndx+outColIndx+2] = input.in[rowIndx+colIndx+2];\n\
	  subcol = subcol + 1;\n\
	  outColIndx = outColIndx + 3;\n\
	}\n\
        col = col + 1;\n\
	colIndx = colIndx + 3;\n\
      }\n\
      outRowIndx = outRowIndx + 1;\n\
    }\n\
  }\n\
\n\
  return 0;\n\
}";

/*
  int greyImage(unsigned char in [100*100*3],
	        unsigned char out[100*100])
 */
char greyImage[] = "\n\
{\n\
  int col;\n\
  int row;\n\
  int indx;\n\
  int outIndx;\n\
\n\
  indx = 0;\n\
  outIndx = 0;\n\
  for (row = 0; row < 100; row = row + 1) {\n\
    for ( col = 0; col < 100; col = col + 1) {\n\
      output.out[outIndx] =  \n\
	(unsigned char) ( (0.299 * input.in[indx]   +\n\
			   0.587 * input.in[indx+1] +\n\
			   0.114 * input.in[indx+2]) + \n\
			  0.5 );\n\
      indx = indx + 3;\n\
      outIndx = outIndx + 1;\n\
    }\n\
  }\n\
  return 0;\n\
}";

/*
  int intGreyImage(unsigned char in [100*100*3],
	           unsigned char out[100*100])
 */
char intGreyImage[] = "\n\
{\n\
  int col;\n\
  int row;\n\
  int indx;\n\
  int outIndx;\n\
\n\
  indx = 0;\n\
  outIndx = 0;\n\
  for (row = 0; row < 100; row = row + 1) {\n\
    for ( col = 0; col < 100; col = col + 1) {\n\
      output.out[outIndx] =  \n\
	(unsigned char) (1 * (input.in[indx]   >> 2) +\n\
			 5 * (input.in[indx+1] >> 3) +\n\
			 1 * (input.in[indx+2] >> 3));\n\
      indx = indx + 3;\n\
      outIndx = outIndx + 1;\n\
    }\n\
  }\n\
  return 0;\n\
}";


/*
  int mixImage(unsigned char in1[100*100*3],
	       unsigned char in2[100*100*3],
	       unsigned char out[100*100*3],
	       float fade)
 */
char mixImage[] = "\n\
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

/*
  int reduceImage(unsigned char in [100*100*3],
		  unsigned char out[50*50*3],
		  int  n = 2)
 */
char reduceImage[] = "\n\
{\n\
  int row;\n\
  int inRow;\n\
  int col;\n\
  int inCol;\n\
  int suprow;\n\
  int supcol;\n\
  int inIndx;\n\
  int outIndx;\n\
  int sumR;\n\
  int sumG;\n\
  int sumB;\n\
  int n2 = n * n;\n\
\n\
  for ( row = 0; row < 50; row = row + 1) {\n\
    for (col = 0; col < 50; col = col + 1) {\n\
      outIndx = (row * 50 * 3) + col * 3;\n\
      sumR = 0; sumG = 0; sumB = 0;\n\
      for (suprow = 0; suprow < n; suprow = suprow + 1) {\n\
	inRow = row * n + suprow;\n\
	for (supcol = 0; supcol < n; supcol = supcol + 1) {\n\
	  inCol = col * n + supcol;\n\
	  inIndx = (inRow * 100 * 3) + inCol * 3;\n\
	  sumR = sumR + input.in[inIndx];\n\
	  sumG = sumG + input.in[inIndx+1];\n\
	  sumB = sumB + input.in[inIndx+2];\n\
	}\n\
      }\n\
      output.out[outIndx+0] = sumR / n2;\n\
      output.out[outIndx+1] = sumG / n2;\n\
      output.out[outIndx+2] = sumB / n2;\n\
    }\n\
  }\n\
  \n\
  return 0;\n\
}";

#define rdtsc(x)  __asm__ __volatile__("rdtsc\n\t" : "=A" (x));

void test_catImage()
{
/*
  int catImage(char in1[100*100*3],
               char in2[100*100*3],
               char out[200*100*3])
*/
  unsigned long long t1, t2;

  static IOField input1_field_list[] = {
    {"in", "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField input2_field_list[] = {
    {"in", "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  static IOField output_field_list[] = {
    {"out", "char[60000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  char out[200*100*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *, char *), result;

  ecl_add_struct_type("input1_type", input1_field_list, context);
  ecl_add_struct_type("input2_type", input2_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input1_type input1, input2_type input2, output_type output)", context);


  gen_code = ecl_code_gen(catImage, context);
  func = (long (*)(char *, char *, char *)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, img2, out);

  rdtsc(t1);
  result = func(img1, img2, out);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P6\n200 100\n255\n", 15) != 15) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 200*100*3) != 200*100*3) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_enlargeImage()
{
/*
  int enlargeImage(char in [100*100*3],
	           char out[200*200*3],
		   int n)
 */

  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out",  "char[120000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  char out[200*200*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *, int), result;

  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output, int n)", context);

  gen_code = ecl_code_gen(enlargeImage, context);
  func = (long (*)(char *, char *, int)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, out, 2);

  rdtsc(t1);
  result = func(img1, out, 2);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P6\n200 200\n255\n", 15) != 15) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 200*200*3) != 200*200*3) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_reduceImage()
{
/*
  int reduceImage(unsigned char in [100*100*3],
		  unsigned char out[50*50*3],
		  int  n = 2)
 */
  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out",  "unsigned integer[7500]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  unsigned char out[50*50*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *, int), result;

  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output, int n)", context);


  gen_code = ecl_code_gen(reduceImage, context);
  func = (long (*)(char *, char *, int)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, out, 2);

  rdtsc(t1);
  result = func(img1, out, 2);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P6\n50 50\n255\n", 13) != 13) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 50*50*3) != 50*50*3) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_intGreyImage()
{
/*
  int intGreyImage(unsigned char in [100*100*3],
	           unsigned char out[100*100])
 */

  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out", "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  unsigned char out[100*100];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *), result;


  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output)", context);

  gen_code = ecl_code_gen(intGreyImage, context);
  func = (long (*)(char *, char *)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, out);

  rdtsc(t1);
  result = func(img1, out);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P5\n100 100\n255\n", 15) != 15) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 100*100) != 100*100) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_greyImage()
{
/*
  int greyImage(unsigned char in [100*100*3],
	        unsigned char out[100*100])
 */

  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out", "unsigned integer[10000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  unsigned char out[100*100];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *), result;


  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output)", context);

  gen_code = ecl_code_gen(greyImage, context);
  func = (long (*)(char *, char *)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, out);

  rdtsc(t1);
  result = func(img1, out);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P5\n100 100\n255\n", 15) != 15) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 100*100) != 100*100) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_edgeImage()
{
/*
  int edgeImage(char in [100*100*3],
                char out[100*100*3])
 */
  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out", "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  char in [100*100*3];
  char out[100*100*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *), result;


  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output)", context);

  gen_code = ecl_code_gen(edgeImage, context);
  func = (long (*)(char *, char *)) (long) gen_code->func;

  if(!func) return;

  result = func(in, out);

  rdtsc(t1);
  result = func(in, out);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);
}

void test_cropImage()
{
/*
  int cropImage(char in [100*100*3],
                char out[50*50*3],
                int inRow, int inCol)
 */
  unsigned long long t1, t2;

  static IOField input_field_list[] = {
    {"in",  "char[30000]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField output_field_list[] = {
    {"out",  "char[7500]", sizeof(char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  char out[50*50*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *, int, int), result;

  ecl_add_struct_type("input_type", input_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input_type input, output_type output, int inRow, int inCol)", context);


  gen_code = ecl_code_gen(cropImage, context);
  func = (long (*)(char *, char *, int, int)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, out, 25, 25);

  rdtsc(t1);
  result = func(img1, out, 25, 25);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P6\n50 50\n255\n", 13) != 13) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 50*50*3) != 50*50*3) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

void test_mixImage()
{
/*
  int mixImage(unsigned char in1[100*100*3],
	       unsigned char in2[100*100*3],
	       unsigned char out[100*100*3],
	       float fade)
 */
  unsigned long long t1, t2;

  static IOField input1_field_list[] = {
    {"in", "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };
  
  static IOField input2_field_list[] = {
    {"in", "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  static IOField output_field_list[] = {
    {"out", "unsigned integer[30000]", sizeof(unsigned char), 0},
    {(void*)0, (void*)0, 0, 0}
  };

  unsigned char out[100*100*3];

  ecl_parse_context context = new_ecl_parse_context();
  ecl_code gen_code;
  long (*func)(char *, char *, char *, int), result;

  ecl_assoc_externs(context, externs);
  ecl_parse_for_context(extern_string, context);

  ecl_add_struct_type("input1_type", input1_field_list, context);
  ecl_add_struct_type("input2_type", input2_field_list, context);
  ecl_add_struct_type("output_type", output_field_list, context);
  ecl_subroutine_declaration("int proc(input1_type input1, input2_type input2, output_type output, float fade)", context);

  gen_code = ecl_code_gen(mixImage, context);
  func = (long (*)(char *, char *, char *, int)) (long) gen_code->func;

  if(!func) return;

  result = func(img1, img2, out, 1.0);

  rdtsc(t1);
  result = func(img1, img2, out, 1.0);
  rdtsc(t2);

  t2 -= t1;
  printf("%0lld\n", t2);

  ecl_code_free(gen_code);
  ecl_free_parse_context(context);

#if defined(IMGDEBUG)
  unlink("out.pnm");
  fd = open("out.pnm", O_CREAT|O_TRUNC|O_RDWR, 0644);
  if(fd < 0) {
    perror("open");
    exit(-1);
  }

  if(write(fd, "P6\n100 100\n255\n", 15) != 15) {
    perror("write");
    exit(-1);
  }

  if(write(fd, out, 100*100*3) != 100*100*3) {
    perror("write");
    exit(-1);
  }

  close(fd);
#endif
}

int main(int argc, char **argv)
{
  if(argc != 2 || !*argv[1] || argv[1][1] || *argv[1]<'1' || *argv[1]>'7') {
    fprintf(stderr, "Usage: imgxform <test-number>\n");
    fprintf(stderr, "\t(1) enlargeImage()\n"
	    "\t(2) reduceImage()\n"
	    "\t(3) intGreyImage()\n"
	    "\t(4) greyImage()\n"
	    "\t(5) catImage()\n"
	    "\t(6) cropImage()\n"
	    "\t(7) mixImage()\n");
    return -1;
  }
  
  switch(*argv[1]) {
  case '1':
    test_enlargeImage();
    break;
  case '2':
    test_reduceImage();
    break;
  case '3':
    test_intGreyImage();
    break;
  case '4':
    test_greyImage();
    break;
  case '5':
    test_catImage();
    break;
  case '6':
    test_cropImage();
    break;
  case '7':
    test_mixImage();
    break;
  default:
    printf("No such choice!\n");
    return -1;
  }
  return 0;
}
