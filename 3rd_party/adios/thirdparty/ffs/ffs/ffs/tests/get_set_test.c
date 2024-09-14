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
#include "ffs.h"

#include "test_funcs.h"

    /* 
       NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE NOTE

    The APIs exercised in this test are not recommended for general use.
    They are designed to work between homogeneous and very tightly coupled
    applications where fixed data types are known to both the encoder and
    decoder.  Use at your own risk.
    */

int
main(int argc, char **argv)
{
    third_rec test_rec;
    void *data_ptr, *tmp_ptr, *encode_buffer, *base_data;
    char *tmp_string;
    FMContext context;
    FMFormat third_format;
    FFSBuffer buf;
    size_t buf_size;
    
    /* setup a test record to play with */
    test_rec.integer_field = -15;
    test_rec.long_field = -2;
    test_rec.string = "penny arcade";
    test_rec.string2 = "xkcd";
    test_rec.double_field = 3.14159;
    
    data_ptr = &test_rec;

    /* stop using test_rec and see if we can access fields using the data_ptr and the field list */
    tmp_ptr = get_FMfieldAddr_by_name(field_list3, "integer field", data_ptr);
    if (*((int*)tmp_ptr) != -15) {
	printf("Failed integer field\n");
	exit(1);
    }

    tmp_ptr = get_FMfieldAddr_by_name(field_list3, "long field", data_ptr);
    if (*((long*)tmp_ptr) != -2) {
	printf("Failed long field\n");
	exit(1);
    }

    /* OK, now get pointers (here, strings) and test */
    tmp_string = get_FMPtrField_by_name(field_list3, "string field", data_ptr, 0);
    if (strcmp(tmp_string, "penny arcade") != 0) {
	printf("failed string field\n");
	exit(1);
    }

    tmp_string = get_FMPtrField_by_name(field_list3, "string field2", data_ptr, 0);
    if (strcmp(tmp_string, "xkcd") != 0) {
	printf("failed string field2\n");
	exit(1);
    }
    
    /* OK, now SET pointers (here, strings) and test */
    if (set_FMPtrField_by_name(field_list3, "string field", data_ptr, strdup("Ack! Thbbft!")) != 1) exit(1);

    if (set_FMPtrField_by_name(field_list3, "string field2", data_ptr, strdup("svelte buoyant waterfowl")) != 1) exit(1);
	   
    if (strcmp(test_rec.string, "Ack! Thbbft!") != 0) {
	printf("Failed Ack! Thbbft!\n");
	exit(1);
    }
    if (strcmp(test_rec.string2, "svelte buoyant waterfowl") != 0) {
	printf("Failed svelte bouyant waterfowl\n");
	exit(1);
    }

    /* register a format so we can do more stuff */
    context = create_local_FMcontext();
    third_format = FMregister_simple_format(context, "third", field_list3, sizeof(test_rec));

    /* encode the data to a buffer to operate upon later */
    encode_buffer = malloc(4096); /* should be plenty */
    buf = create_fixed_FFSBuffer(encode_buffer, 4096);
    FFSencode(buf, third_format, data_ptr, &buf_size);

    /* free the pointer part of the original data */
    FMfree_var_rec_elements(third_format, data_ptr);

    /* now see if we can operate on the encoded data */
    base_data = FMheader_skip(context, encode_buffer);
    
    /* stop using test_rec and see if we can access fields using the data_ptr and the field list */
    tmp_ptr = get_FMfieldAddr_by_name(field_list3, "integer field", base_data);
    if (*((int*)tmp_ptr) != -15) {
	printf("failed second integer field\n");
	exit(1);
    }

    tmp_ptr = get_FMfieldAddr_by_name(field_list3, "long field", base_data);
    if (*((long*)tmp_ptr) != -2) {
	printf("failed second long field\n");
	exit(1);
    }

    /* OK, now get pointers (here, strings) and test */
    tmp_string = get_FMPtrField_by_name(field_list3, "string field", base_data, 1);
    if (strcmp(tmp_string, "Ack! Thbbft!") != 0) {
	printf("Failed second Ack!\n");
	exit(1);
    }

    tmp_string = get_FMPtrField_by_name(field_list3, "string field2", base_data, 1);
    if (strcmp(tmp_string, "svelte buoyant waterfowl") != 0) {
	printf("Failed second svelte\n");
	exit(1);
    }
    
    free_FFSBuffer(buf);
    free(encode_buffer);
    free_FMcontext(context);
    
    return 0;
}


