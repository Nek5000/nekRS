#include "../config.h"
#include <assert.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <stdio.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include "fm.h"
#include "cod.h"
#include "ffs.h"


    /* 
       NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE

    The APIs exercised in this test are not necessarily for you.  They are
    designed for a specific ADIOS scenario where the information to be
    transfered becomes available piece-wise.  In particular, the atomic
    variables (which include size values for dynamic arrays) are available
    first.  Later, the actual data for those arrays becomes available
    perhaps one-by-one, *AND* the semantics are such that that array data
    must be copied at the time it becomes available.  Also, the data
    structure in question is relatively simple.  I.E. it's a base structure
    with dynamic arrays hanging off of it, but those arrays are of simple
    types (have no pointers).  To optimize this, we fill in the atomic
    values in the base data structure and use FFSencode_no_leaf_copy().
    This results in a contiguous encoded block with spots reserved for the
    array data (but no memcpy's to fill that in have been performed).
    Later, we get that array data, memcpy the data into the appropriate
    place, and our encode is complete.
    */


typedef struct Struct_array1 {
	int x;
	int y;
	int *array1;
} Struct_array1, *Struct_array1_ptr;

FMField struct1_flds[] = {
    {"x", "integer", sizeof(int), FMOffset(Struct_array1_ptr, x)},
    {"y", "integer", sizeof(int), FMOffset(Struct_array1_ptr, y)},
    {"array1", "integer[x][y]", sizeof(int), FMOffset(Struct_array1_ptr, array1)},
    {(char *) 0, (char *) 0, 0, 0}
};
    
typedef struct Struct_array2 {
	int x;
	double *array2;
} Struct_array2, *Struct_array2_ptr;

FMField struct2_flds[] = {
    {"x", "integer", sizeof(int), FMOffset(Struct_array2_ptr, x)},
    {"array2", "double[x]", sizeof(int), FMOffset(Struct_array2_ptr, array2)},
    {(char *) 0, (char *) 0, 0, 0}
};

typedef struct data {
    int x;
    int y;
    Struct_array1 s1;
    Struct_array2 s2;
} data_struct, *data_struct_ptr;

FMField data_flds[] = {
    {"x", "integer", sizeof(int), FMOffset(data_struct_ptr, x)},
    {"y", "integer", sizeof(int), FMOffset(data_struct_ptr, y)},
    {"s1", "Struct_array1", sizeof(Struct_array1), FMOffset(data_struct_ptr, s1)},
    {"s2", "Struct_array2", sizeof(Struct_array2), FMOffset(data_struct_ptr, s2)},
    {(char *) 0, (char *) 0, 0, 0}
};
    
FMStructDescRec data_struct_list[] = {
    {"Data", data_flds, sizeof(data_struct), NULL},
    {"Struct_array1", struct1_flds, sizeof(struct Struct_array1), NULL},
    {"Struct_array2", struct2_flds, sizeof(struct Struct_array2), NULL},
    {NULL, NULL, 0, NULL}
};
    
int
main(int argc, char **argv)
{
    data_struct data;
    void *encode_buffer, *base_of_encoded_array1, *base_of_encoded_array2;
    FMContext fmcontext;
    cod_parse_context cod_context = new_cod_parse_context();
    cod_code gen_code1;
    cod_code gen_code2;
    void* (*func1)(void*);
    void* (*func2)(void*);
    FMFormat format;
    FFSBuffer buf;
    size_t buf_size;
    int test_array1[10][5];
    double test_array2[10];
    int x, y;
    typedef int (*pointer_to_2darray)[10][5];
    /* setup a test record to play with */
    data.s1.x = data.s2.x = data.x = 10;
    data.s1.y = data.y = 5;
    data.s1.array1 = (void*)(intptr_t)(0xdeadbeef);    /* should get a seg fault if we dereference this */
    data.s2.array2 = (void*)(intptr_t)(0xD15EA5E);    /* should get a seg fault if we dereference this */


    /* register a format so we can do more stuff */
    fmcontext = create_local_FMcontext();
    format = FMregister_data_format(fmcontext, data_struct_list);

    /* encode the data to a buffer to operate upon later */
    encode_buffer = malloc(4096); /* should be plenty */
    buf = create_fixed_FFSBuffer(encode_buffer, 4096);
    FFSencode_no_leaf_copy(buf, format, &data, &buf_size);

    /* now see if we can operate on the encoded data */
    cod_add_encoded_param("encode", encode_buffer, 0, fmcontext, cod_context);
    cod_set_return_type("Struct_array1", cod_context);
    
    gen_code1 = cod_code_gen("{return encode.s1.array1;}", cod_context);
    func1 = (void*(*)(void*)) (intptr_t) gen_code1->func;
    base_of_encoded_array1 = func1(encode_buffer);
    
    gen_code2 = cod_code_gen("{return encode.s2.array2;}", cod_context);
    func2 = (void*(*)(void*)) (intptr_t) gen_code2->func;
    base_of_encoded_array2 = func2(encode_buffer);

    for (x=0; x < 10; x++) {
	for(y=0; y < 5; y++) {
	    test_array1[x][y] = x*100 + y;
	}
	test_array2[x] = 10*x;
    }
    memcpy(base_of_encoded_array1, &test_array1[0][0], sizeof(test_array1));
    memcpy(base_of_encoded_array2, &test_array2[0], sizeof(test_array2));
    
    /* what follows from here is simply testing to see the stuff above went well */
    data_struct *actual_struct;
    FFSContext FFScontext = create_FFSContext_FM(fmcontext);
    int ret = 0;
    (void) FFSset_fixed_target(FFScontext, data_struct_list);
    (void) FFSdecode_in_place(FFScontext, encode_buffer, (void**)&actual_struct);
    
    if ((*((pointer_to_2darray)actual_struct->s1.array1))[0][0] != 000) {
	printf("Struct.s1.array1[0][0] is bad!\n"); ret++;
    }
    if ((*((pointer_to_2darray)actual_struct->s1.array1))[2][3] != 203) {
	printf("Struct.s1.array1[2][3] is bad!\n"); ret++;
    }
    if ((*((pointer_to_2darray)actual_struct->s1.array1))[9][4] != 904) {
	printf("Struct.s1.array1[9][4] is bad!\n"); ret++;
    }
    if (actual_struct->s2.array2[0] != 0) {
	printf("Struct.s2.array2[0] is bad!  (%g)\n", actual_struct->s2.array2[0]); ret++;
    }
    if (actual_struct->s2.array2[5] != 50) {
	printf("Struct.s2.array2[5] is bad! (%g)\n", actual_struct->s2.array2[5]); ret++;
    }
    if (actual_struct->s2.array2[9] != 90) {
	printf("Struct.s2.array2[9] is bad!  (%g)\n", actual_struct->s2.array2[9]); ret++;
    }
    free_FFSBuffer(buf);
    free(encode_buffer);
    free_FMcontext(fmcontext);
    free_FFSContext(FFScontext);
    cod_code_free(gen_code1);
    cod_code_free(gen_code2);
    cod_free_parse_context(cod_context);
    
    return ret;
}


