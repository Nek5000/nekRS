#include <stdio.h>
#include <stdlib.h>

#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 1418)
#endif
#ifdef _MSC_VER
#include <stdlib.h>
#define srand48(s) srand(s)
#define drand48() (rand()*(1./RAND_MAX))
#define lrand48() ((long long)rand() << 32 | rand())
#endif
typedef struct _complex_rec {
    double r;
    double i;
} complex, *complex_ptr;

typedef struct _nested_rec {
    complex item;
} nested, *nested_ptr;

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;


int 
filter(void* input, void* output, void *data)
{
    (void) output; (void) data;
    return ((simple_rec_ptr)input)->long_field %2;
}

typedef struct _output_rec {
    int random_field;
    int integer_field;
    short sum1_field;
    long sum2_field;
    nested nested_field;
    double double_field;
    char char_field;
    int sum3_field;
} output_rec, *output_rec_ptr;

extern int
transform(simple_rec_ptr input, output_rec_ptr output)
{
    output->random_field = (lrand48() % 10);
    output->integer_field = input->integer_field;
    output->sum1_field = input->short_field;
    output->sum2_field = input->long_field - output->random_field;
    output->nested_field.item.i = input->nested_field.item.i;
    output->nested_field.item.r = input->nested_field.item.r;
    output->double_field = input->double_field;
    output->char_field = input->char_field;
    output->sum3_field = input->scan_sum;
    return input->long_field % 2;
}

