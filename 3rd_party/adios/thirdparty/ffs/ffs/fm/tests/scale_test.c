
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include <string.h>
#include <time.h>
#include "fm.h"
#include "ffs.h"

#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define sleep(x) Sleep(1000*x)
#else
extern int sleep();
#endif

#include "test_funcs.h"

char *gen_name(int i)
{
    char tmp_name[128];
    snprintf(tmp_name, sizeof(tmp_name), "SST_Variable_FieldName that's really really long because I can't imagine why %d", i);
    return strdup(tmp_name);
}

struct base_elem {
  int64_t elem_count1;
  int64_t elem_count2;
  int64_t *array1;
  int64_t *array2;
};

#ifdef _MSC_VER
// this is not important on windows
#define clock_gettime(x,y)
#endif

int
main(int argc, char **argv)
{

    FMStructDescRec str_list[5];
    struct timespec start = {0,0}, stop = {0,0};
    
    FMContext context;
    int field_count = 20000;
    field_count = ((field_count >> 2 ) << 2); // ensure field count is divisible by 4;
    FMFieldList list = malloc(sizeof(struct _FMField) * (field_count + 1));
    int cur_count = 0;
    while (cur_count < field_count) {
        /* do 4 at a time */
        char tmp[128];
        char *n1 = gen_name(cur_count);
        char *n2 = gen_name(cur_count + 1);
        char *n3 = gen_name(cur_count + 2);
        char *n4 = gen_name(cur_count + 3);
        list[cur_count].field_name = n1;
        list[cur_count].field_type = strdup("integer");
        list[cur_count].field_size = 8;
        list[cur_count].field_offset = cur_count * 8;
        list[cur_count+1].field_name = n2;
        list[cur_count+1].field_type = strdup("integer");
        list[cur_count+1].field_size = 8;
        list[cur_count+1].field_offset = (cur_count+1) * 8;
        list[cur_count+2].field_name = n3;
        snprintf(tmp, sizeof(tmp), "integer[%s]", n1);
        list[cur_count+2].field_type = strdup(tmp);
        list[cur_count+2].field_size = 8;
        list[cur_count+2].field_offset = (cur_count+2) * 8;
        list[cur_count+3].field_name = n4;
        snprintf(tmp, sizeof(tmp), "integer[%s]", n2);
        list[cur_count+3].field_type = strdup(tmp);
        list[cur_count+3].field_size = 8;
        list[cur_count+3].field_offset = (cur_count+3) * 8;
        cur_count +=4;
    }
    list[cur_count].field_name = list[cur_count].field_type = NULL;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    str_list[0].format_name = "first format";
    str_list[0].field_list = list;
    str_list[0].struct_size = sizeof(first_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    FFSContext fc = create_FFSContext();
    context = FMContext_from_FFS(fc);

    FMFormat format = register_data_format(context, str_list);

    clock_gettime(CLOCK_MONOTONIC, &stop);
    double duration = (stop.tv_sec + 1.0e-9*stop.tv_nsec) - (start.tv_sec + 1.0e-9*start.tv_nsec);
    printf("Registration took %g seconds\n", duration);
    
    char * buf= malloc(sizeof(struct base_elem) * (field_count / 4));
    int i;
    for (i = 0; i < field_count ; i+= 4) {
      struct base_elem *tmp = (struct base_elem *) (buf + i * sizeof(int64_t));
      tmp->elem_count1 = 3;
      tmp->elem_count2 = 3;
      tmp->array1 = malloc(3*sizeof(tmp->array1[0]));
      tmp->array2 = malloc(3*sizeof(tmp->array2[0]));
    }
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    FFSBuffer b = create_FFSBuffer();
    
    size_t buf_size;
    char *encode = FFSencode(b, format, buf, &buf_size);
    FFSTypeHandle th = FFSTypeHandle_from_encode(fc, encode);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    duration = (stop.tv_sec + 1.0e-9*stop.tv_nsec) - (start.tv_sec + 1.0e-9*start.tv_nsec);
    printf("encode took %g seconds\n", duration);
    clock_gettime(CLOCK_MONOTONIC, &start);
    establish_conversion(fc, th, str_list);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    duration = (stop.tv_sec + 1.0e-9*stop.tv_nsec) - (start.tv_sec + 1.0e-9*start.tv_nsec);
    printf("establish took %g seconds\n", duration);

}


