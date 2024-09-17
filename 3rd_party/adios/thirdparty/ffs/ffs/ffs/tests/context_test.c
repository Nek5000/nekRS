
#include "config.h"
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_NETINET_IN_H
#include <netinet/in.h>
#endif
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <string.h>
#include <assert.h>
#include "ffs.h"

#include "test_funcs.h"
#include <sys/types.h>
#include <sys/stat.h>

#if defined(_MSC_VER)
#include "winsock.h"   // ntohl
#endif

extern void dump_FMFormat(FMFormat ioformat);

static void test_receive(char *buffer, size_t buf_size, int finished,
			       int test_level);
static void test_all_receive(char *buffer, size_t buf_size, int finished);
static void write_buffer(FMFormat format, char *buf, size_t size);
static void read_test_only();

static int write_output = 0;
static char *output_file = NULL;
static char *read_file = NULL;
static int fail = 0;
static char *test_only = NULL;

static FFSContext rcv_context = NULL;
static int verbose = 0;
static FMContext loaded_FMcontext = NULL;

int
main(int argc, char **argv)
{
    FMContext src_context;
    FFSBuffer encode_buffer;
    FMFormat first_rec_ioformat, second_rec_ioformat, third_rec_ioformat;
    FMFormat fourth_rec_ioformat, later_ioformat, nested_ioformat;
    FMFormat embedded_rec_ioformat, fifth_rec_ioformat, triangle_ioformat;
    FMStructDescRec str_list[5];
    char *xfer_buffer;
    size_t buf_size;
    first_rec rec1;
    first_rec array1[10];
    second_rec rec2;
    third_rec rec3;
    fourth_rec rec4;
    later_rec rec5;
    later_rec2 rec6;
    nested_rec rec7;
    fifth_rec emb_array;
    sixth_rec var_array;
    ninth_rec var_var;
    string_array_rec str_array;
    int i, j;
    FMFormat sixth_rec_ioformat, ninth_rec_ioformat, string_array_ioformat;
    FMFormat derive_ioformat, multi_array_ioformat, add_action_ioformat;

    init_written_data();

    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-w") == 0) {
	    output_file = argv[++i];
	    write_output++;
	} else if (strcmp(argv[i], "-r") == 0) {
	    read_file = argv[++i];
	} else if (strcmp(argv[i], "-v") == 0) {
	    verbose++;
	} else if (strcmp(argv[i], "-t") == 0) {
	    test_only = argv[++i];
	} else {
	    printf("Unknown argument %s\n", argv[i]);
	    printf("Usage:\n\t-w\t  write test output\n");
	    printf("\t-r file\t  rest/process test data in <file>\n");
	    printf("\t-t test\t  test only format <test>\n");
	    exit(1);
	}
    }

    if (read_file) {
	read_test_only();
	free_written_data();
	if (rcv_context != NULL) {
	    free_FFSContext(rcv_context);
	    rcv_context = NULL;
	}
	if (loaded_FMcontext) {
	    free_FMcontext(loaded_FMcontext);
	}
	if (fail) exit(1);
	exit(0);
    }
    src_context = create_FMcontext();
    encode_buffer = create_FFSBuffer();
    str_list[0].format_name = "first format";
    str_list[0].field_list = field_list;
    str_list[0].struct_size = sizeof(first_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    first_rec_ioformat = register_data_format(src_context, str_list);
    if (!first_rec_ioformat) {
	printf("Format server failure\n");
	exit(1);
    }

    str_list[0].format_name = "string format";
    str_list[0].field_list = field_list2;
    str_list[0].struct_size = sizeof(second_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    second_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "two string format";
    str_list[0].field_list = field_list3;
    str_list[0].struct_size = sizeof(third_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    third_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "internal array format";
    str_list[0].field_list = field_list4;
    str_list[0].struct_size = sizeof(fourth_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    fourth_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "embedded";
    str_list[0].field_list = embedded_field_list;
    str_list[0].struct_size = sizeof(embedded_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    embedded_rec_ioformat = register_data_format(src_context, str_list);
    (void) embedded_rec_ioformat;

    str_list[0].format_name = "structured array format";
    str_list[0].field_list = field_list5;
    str_list[0].struct_size = sizeof(fifth_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "embedded";
    str_list[1].field_list = embedded_field_list;
    str_list[1].struct_size = sizeof(embedded_rec);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = NULL;
    fifth_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "variant array format";
    str_list[0].field_list = field_list6;
    str_list[0].struct_size = sizeof(sixth_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "string format";
    str_list[1].field_list = field_list2;
    str_list[1].struct_size = sizeof(second_rec);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = NULL;
    sixth_rec_ioformat = register_data_format(src_context, str_list);

    memset((char *) &rec1, 0, sizeof(rec1));
    rec1.integer_field = 14;
    rec1.double_field = 2.717;
    rec1.char_field = 'A';
    xfer_buffer = FFSencode(encode_buffer, first_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(first_rec_ioformat, xfer_buffer, buf_size);
    memset((char *) &emb_array, 0, sizeof(emb_array));
    emb_array.earray[0].dfield = 4.0;
    emb_array.earray[0].ifield = 4;
    emb_array.earray[0].string = (char *) malloc(10);
    memset(emb_array.earray[0].string, 0, 10);
    sprintf(emb_array.earray[0].string, "string%d", emb_array.earray[0].ifield * 5);
    emb_array.earray[1].dfield = 3.0;
    emb_array.earray[1].ifield = 3;
    emb_array.earray[1].string = (char *) malloc(10);
    memset(emb_array.earray[1].string, 0, 10);
    sprintf(emb_array.earray[1].string, "string%d", emb_array.earray[1].ifield * 5);
    emb_array.earray[2].dfield = 2.0;
    emb_array.earray[2].ifield = 2;
    emb_array.earray[2].string = (char *) malloc(10);
    memset(emb_array.earray[2].string, 0, 10);
    sprintf(emb_array.earray[2].string, "string%d", emb_array.earray[2].ifield * 5);
    emb_array.earray[3].dfield = 1.0;
    emb_array.earray[3].ifield = 1;
    emb_array.earray[3].string = (char *) malloc(10);
    memset(emb_array.earray[3].string, 0, 10);
    sprintf(emb_array.earray[3].string, "string%d", emb_array.earray[3].ifield * 5);
    xfer_buffer = FFSencode(encode_buffer, fifth_rec_ioformat, &emb_array, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(fifth_rec_ioformat, xfer_buffer, buf_size);
    free(emb_array.earray[0].string);
    free(emb_array.earray[1].string);
    free(emb_array.earray[2].string);
    free(emb_array.earray[3].string);

    memset((char *) &rec2, 0, sizeof(rec2));
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "testing";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = FFSencode(encode_buffer, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(second_rec_ioformat, xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = FFSencode(encode_buffer, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(second_rec_ioformat, xfer_buffer, buf_size);
    rec1.integer_field = 17;
    rec1.double_field *= 3.0;
    rec1.char_field = 'B';
    xfer_buffer = FFSencode(encode_buffer, first_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(first_rec_ioformat, xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = FFSencode(encode_buffer, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(second_rec_ioformat, xfer_buffer, buf_size);
    rec1.integer_field *= 2;
    rec1.double_field *= 2.717;
    rec1.char_field = 'C';
/*    write_comment_IOfile(iofile, "this is another comment in the file");*/
    xfer_buffer = FFSencode(encode_buffer, first_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(first_rec_ioformat, xfer_buffer, buf_size);
    memset((char *) &rec3, 0, sizeof(rec3));
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3.ulong_field = 0xf7e589ceec9dd130;
#else
    rec3.ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Red_Stripe;
    xfer_buffer = FFSencode(encode_buffer, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(third_rec_ioformat, xfer_buffer, buf_size);

    str_list[0].format_name = "later format";
    str_list[0].field_list = later_field_list;
    str_list[0].struct_size = sizeof(later_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    later_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "nested format";
    str_list[0].field_list = nested_field_list;
    str_list[0].struct_size = sizeof(nested_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "string format";
    str_list[1].field_list = field_list2;
    str_list[1].struct_size = sizeof(second_rec);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = NULL;
    nested_ioformat = register_data_format(src_context, str_list);

    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Paulaner;
    xfer_buffer = FFSencode(encode_buffer, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(third_rec_ioformat, xfer_buffer, buf_size);
    memset((char *) &rec7, 0, sizeof(rec7));
    rec7.integer_field = 47;
    rec7.nested_rec.integer_field = 14;
    rec7.nested_rec.short_field = 27;
    rec7.nested_rec.long_field = 987234;
    rec7.nested_rec.string = "Another string";
    rec7.nested_rec.double_field = 2.717;
    rec7.nested_rec.char_field = 'A';
    rec7.string = "Yet another string";
    xfer_buffer = FFSencode(encode_buffer, nested_ioformat, &rec7, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(nested_ioformat, xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    xfer_buffer = FFSencode(encode_buffer, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(third_rec_ioformat, xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Red_Stripe;
    xfer_buffer = FFSencode(encode_buffer, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(third_rec_ioformat, xfer_buffer, buf_size);
    memset((char *) &rec5, 0, sizeof(rec5));
    rec5.integer_field = 9872346;
    rec5.string = "ABCD";
    rec5.double_field = 3.14159265358797323;
    xfer_buffer = FFSencode(encode_buffer, later_ioformat, &rec5, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(later_ioformat, xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    xfer_buffer = FFSencode(encode_buffer, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(third_rec_ioformat, xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "the end";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = FFSencode(encode_buffer, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(second_rec_ioformat, xfer_buffer, buf_size);
    for (i = 0; i < 10; i++) {
	memset((char *) &array1[i], 0, sizeof(array1[i]));
	array1[i].integer_field = 2 * i * i;
	array1[i].double_field = 2.717 * (i * i);
	array1[i].char_field = 'D' + i;
    }
    str_list[0].format_name = "later format";
    str_list[0].field_list = later_field_list2;
    str_list[0].struct_size = sizeof(later_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    later_ioformat = register_data_format(src_context, str_list);

/*    if (!write_array_IOfile(iofile, first_rec_ioformat, &array1[0],
			    10, sizeof(array1[0])))
	IOperror(iofile, "write failed\n");*/

    memset((char *) &rec4, 0, sizeof(rec4));
    for (i = 0; i < ARRAY_SIZE; i++) {
	rec4.int_array[i] = 297 + i;
    }
    rec4.double_array[0][0] = 1.0;
    rec4.double_array[0][1] = 2.0;
    rec4.double_array[1][0] = 3.0;
    rec4.double_array[1][1] = 4.0;
    rec4.ifield = -rec4.int_array[ARRAY_SIZE - 1];
    xfer_buffer = FFSencode(encode_buffer, fourth_rec_ioformat, &rec4, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(fourth_rec_ioformat, xfer_buffer, buf_size);
    memset((char *) &rec6, 0, sizeof(rec6));
    rec6.integer_field = 23462346;
    rec6.string = "Efghij";
    rec6.double_field = 3.14159265358797323 * 2.0;
    xfer_buffer = FFSencode(encode_buffer, later_ioformat, &rec6, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(later_ioformat, xfer_buffer, buf_size);
    rec6.integer_field = 2346987;
    rec6.string = "Klmn";
    rec6.double_field = 3.14159265358797323 * 3.0;
    xfer_buffer = FFSencode(encode_buffer, later_ioformat, &rec6, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(later_ioformat, xfer_buffer, buf_size);

    for (i = 1; i < 20; i += 5) {
        memset((char *) &var_array, 0, sizeof(var_array));
	var_array.string = malloc(10);
	memset(var_array.string, 0, 10);
	sprintf(var_array.string, "variant%d", i);
	var_array.icount = 2 * i;
	var_array.var_int_array = malloc(sizeof(((sixth_rec_ptr) 0)->var_int_array[0]) * var_array.icount);
	var_array.var_double_array = malloc(sizeof(double) * var_array.icount);
	var_array.var_string_array = malloc(sizeof(second_rec) * var_array.icount);
        memset((char *) var_array.var_string_array, 0, 
	       sizeof(second_rec) * var_array.icount);
	for (j = 0; j < var_array.icount; j++) {
	    var_array.var_int_array[j] = 297 + j;
	    var_array.var_double_array[j] = 2.717 * j;
	    var_array.var_string_array[j].integer_field = 345 * j;
	    var_array.var_string_array[j].short_field = j;
	    var_array.var_string_array[j].long_field = 785 * j;
	    var_array.var_string_array[j].string = malloc(15);
	    memset(var_array.var_string_array[j].string, 0, 15);
	    sprintf(var_array.var_string_array[j].string,
		    "substring%d", j);
	    var_array.var_string_array[j].double_field = 3.1415 * j;
	    var_array.var_string_array[j].char_field = 'a' + 2 * j;
	}

	xfer_buffer = FFSencode(encode_buffer, sixth_rec_ioformat, 
					      &var_array, &buf_size);
	test_all_receive(xfer_buffer, buf_size, 0);
	write_buffer(sixth_rec_ioformat, xfer_buffer, buf_size);
	for (j = 0; j < var_array.icount; j++) {
	    free(var_array.var_string_array[j].string);
	}
	free(var_array.var_string_array);
	free(var_array.string);
	free(var_array.var_int_array);
	free(var_array.var_double_array);
    }
    str_list[0].format_name = "EventV";
    str_list[0].field_list = field_list9;
    str_list[0].struct_size = sizeof(ninth_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "EventVecElem";
    str_list[1].field_list = event_vec_elem_fields;
    str_list[1].struct_size = sizeof(struct _io_encode_vec);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = NULL;
    ninth_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "string_array";
    str_list[0].field_list = string_array_field_list;
    str_list[0].struct_size = sizeof(string_array_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    string_array_ioformat = register_data_format(src_context, str_list);

    for (i = 1; i < 10; i += 2) {
        memset((char *) &var_var, 0, sizeof(var_var));
	var_var.vec_length = i;
	var_var.eventv = malloc(sizeof(((ninth_rec_ptr) 0)->eventv[0]) * var_var.vec_length);
        memset((char *) var_var.eventv, 0, sizeof(((ninth_rec_ptr) 0)->eventv[0]) * var_var.vec_length);

        memset((char *) &str_array, 0, sizeof(str_array));
	str_array.array_len = i;
	str_array.array = malloc(sizeof(char*) * i);

	for (j = 0; j < var_var.vec_length; j++) {
	    int k;
	    var_var.eventv[j].iov_len = j + i;
	    var_var.eventv[j].iov_base = malloc(j + i);
	    str_array.array[j] = malloc(i + j + 2);
	    for (k=0; k<j+i; k++) {
		((char*)var_var.eventv[j].iov_base)[k] = 'A' + k + i/5;
		str_array.array[j][k] = 'a' + k + i;
	    }
	    str_array.array[j][k] = 0;
	}
	if ((i %4) == 1) {
	    str_array.base_string = NULL;
	} else {
	    str_array.base_string = strdup("Whoa there!");
	}
	xfer_buffer = FFSencode(encode_buffer, ninth_rec_ioformat, 
					      &var_var, &buf_size);
	test_all_receive(xfer_buffer, buf_size, 0);
	write_buffer(ninth_rec_ioformat, xfer_buffer, buf_size);

	xfer_buffer = FFSencode(encode_buffer, 
					      string_array_ioformat, 
					      &str_array, &buf_size);
	test_all_receive(xfer_buffer, buf_size, 0);
	write_buffer(string_array_ioformat, xfer_buffer, buf_size);

	for (j = 0; j < var_var.vec_length; j++) {
	    free(var_var.eventv[j].iov_base);
	    free(str_array.array[j]);
	}
	if (str_array.base_string) free(str_array.base_string);
	free(var_var.eventv);
	free(str_array.array);
    }
    str_list[0].format_name = "Channel Derive";
    str_list[0].field_list = derive_msg_field_list;
    str_list[0].struct_size = sizeof(DeriveMsg);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "IOfield_list";
    str_list[1].field_list = field_list_flds;
    str_list[1].struct_size = sizeof(FMField);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = "DEFormatList";
    str_list[2].field_list = format_list_field_list;
    str_list[2].struct_size = sizeof(format_list_element);
    str_list[2].opt_info = NULL;
    str_list[3].format_name = "channel_ID";
    str_list[3].field_list = channel_id_flds;
    str_list[3].struct_size = sizeof(channel_ID_struct);
    str_list[3].opt_info = NULL;

    str_list[4].format_name = NULL;
    derive_ioformat = register_data_format(src_context, str_list);
    xfer_buffer = FFSencode(encode_buffer, derive_ioformat,
					  &derive, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(derive_ioformat, xfer_buffer, buf_size);

    str_list[0].format_name = "multi_array";
    str_list[0].field_list = multi_array_flds;
    str_list[0].struct_size = sizeof(multi_array);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    multi_array_ioformat = register_data_format(src_context, str_list);
    xfer_buffer = FFSencode(encode_buffer, multi_array_ioformat,
					  &multi_array, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(multi_array_ioformat, xfer_buffer, buf_size);

    str_list[0].format_name = "triangle_param";
    str_list[0].field_list = triangle_field;
    str_list[0].struct_size = sizeof(triangle);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "compressed_mesh_param";
    str_list[1].field_list = compressed_mesh;
    str_list[1].struct_size = sizeof(compressed_mesh_param);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = NULL;
    triangle_ioformat = register_data_format(src_context, str_list);

    xfer_buffer = FFSencode(encode_buffer, triangle_ioformat,
					  &triangle, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(triangle_ioformat, xfer_buffer, buf_size);
    

    str_list[0].format_name = "add_action";
    str_list[0].field_list = add_field_list;
    str_list[0].struct_size = sizeof(add_action_record);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = "XMLFormatList";
    str_list[1].field_list = xml_format_list_flds;
    str_list[1].struct_size = sizeof(msg_format_list_element);
    str_list[1].opt_info = NULL;
    str_list[2].format_name = "IOfield_list";
    str_list[2].field_list = field_list_flds;
    str_list[2].struct_size = sizeof(FMField);
    str_list[2].opt_info = NULL;
    str_list[3].format_name = NULL;
    add_action_ioformat = register_data_format(src_context, str_list);

    xfer_buffer = FFSencode(encode_buffer, add_action_ioformat,
					  &add_action_record, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(add_action_ioformat, xfer_buffer, buf_size);
    

    str_list[0].format_name = "node";
    str_list[0].field_list = node_field_list;
    str_list[0].struct_size = sizeof(struct node);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
/*    FMFormat node_ioformat;
    node_ioformat = register_data_format(src_context, str_list);

    for (i = 0; i < sizeof(nodes)/sizeof(nodes[0]); i++) {
	nodes[i].node_num = i;
	nodes[i].link1 = nodes[i].link2 = NULL;
    }

    for (i=0; i <  sizeof(nodes)/sizeof(nodes[0]) - 1; i++) {
	nodes[i].link1 = &nodes[i+1];
    }
    nodes[0].link2 = &nodes[sizeof(nodes)/sizeof(nodes[0])-1];
    v.node_count = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);

    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(node_ioformat, xfer_buffer, buf_size);

    nodes[0].link2 = NULL;
    nodes[sizeof(nodes)/sizeof(nodes[0]) - 1].link1 = &nodes[2];
    v.node_count = 0;
    nodes[0].node_num = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);
    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(node_ioformat, xfer_buffer, buf_size);

    for (i=0; i <  sizeof(nodes)/sizeof(nodes[0]) - 1; i++) {
	nodes[i].link1 = nodes[i].link2 = NULL;
    }
    nodes[0].link1 = &nodes[1];
    nodes[0].link2 = &nodes[2];
    nodes[1].link1 = &nodes[3];
    nodes[1].link2 = &nodes[4];
    nodes[2].link1 = &nodes[5];
    nodes[2].link2 = &nodes[6];
    nodes[3].link1 = &nodes[7];
    nodes[3].link2 = &nodes[8];
    nodes[4].link1 = &nodes[9];

    v.node_count = 0;
    nodes[0].node_num = 0;
    nodes[0].node_num = calc_signature(&nodes[0], &v);
    xfer_buffer = FFSencode(encode_buffer, node_ioformat,
					  &nodes[0], &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(node_ioformat, xfer_buffer, buf_size);
*/    
    free_FMcontext(src_context);
    free_FFSBuffer(encode_buffer);
    src_context = NULL;
    test_all_receive(NULL, 0, 1);
    write_buffer(first_rec_ioformat, NULL, 0);
    free_written_data();
    if (rcv_context != NULL) {
	free_FFSContext(rcv_context);
	rcv_context = NULL;
    }
    if (loaded_FMcontext) {
	free_FMcontext(loaded_FMcontext);
    }
    if (fail) exit(1);
    return 0;
}

/* NT needs O_BINARY, but it doesn't exist elsewhere */
#ifndef O_BINARY
#define O_BINARY 0
#endif

static char *
get_buffer(size_t *size_p)
{
    static int file_fd = 0;
    static char *buffer = NULL;
    char *tmp_buffer;
    static int last_size = -1;
    int to_read;
    int tmp_size;
    unsigned int beef = 0xdeadbeef;
    int indicator;

    if (read_file == NULL) exit(1);

    if (loaded_FMcontext == NULL) {
	loaded_FMcontext = create_local_FMcontext();
    }
    if (file_fd == 0) {
	file_fd = open(read_file, O_RDONLY|O_BINARY, 0777);
	buffer = malloc(1);
    }
    if (last_size != -1) {
	if (memcmp(buffer+last_size, &beef, 4) != 0) {
	    printf("memory overwrite error\n");
	}
    }
    if(read(file_fd, &indicator, 4) != 4) exit(1);

    indicator = ntohl(indicator);
    if ((indicator >> 24) == 0x2) {
	/* got a format coming in */
	int format_rep_size, format_id_size;
	char *format_id, *format_rep;
	if (read(file_fd, &format_rep_size, 4) != 4) exit(1);
	format_rep_size = ntohl(format_rep_size);
	format_id_size = indicator & 0xff;
	format_id = malloc(format_id_size);
	format_rep = malloc(format_rep_size);
	if (read(file_fd, format_id, format_id_size) != format_id_size) exit(1);
	if (read(file_fd, format_rep, format_rep_size) != format_rep_size) exit(1);
	(void) load_external_format_FMcontext(loaded_FMcontext, format_id,
					      format_id_size, format_rep);
	if(read(file_fd, &indicator, 4) != 4) exit(1);
	indicator = ntohl(indicator);
	free(format_id);
    }
    if ((indicator >> 24) != 0x3) printf("BAD!\n");
    to_read = indicator & 0xffffff;
    last_size = to_read;
    buffer = realloc(buffer, to_read+4);
    tmp_buffer = buffer;
    while((tmp_size = read(file_fd, tmp_buffer, to_read)) != to_read) {
	if (tmp_size == 0) {
	    free(buffer);
	    return NULL;
	} else if (tmp_size == -1) {
	    perror("Read failure");
	    free(buffer);
	    return NULL;
	}
	to_read -= tmp_size;
	tmp_buffer += tmp_size;
    } 
    memcpy(buffer+last_size, &beef, 4);
    if (last_size == 0) {
	free(buffer);
	close(file_fd);
	file_fd = 0;
	return NULL;
    } else {
	*size_p = last_size;
	return buffer;
    }
}

static void
read_test_only()
{
    char *input;
    size_t size;
    while ((input = get_buffer(&size)) != NULL) {
	test_all_receive(input, size, 0);
    }
    test_all_receive(NULL, 0, 1);
    free(input);
}

static FFSTypeHandle first_rec_ioformat, second_rec_ioformat, third_rec_ioformat;
static FFSTypeHandle fourth_rec_ioformat, later_rec_ioformat, nested_rec_ioformat;
static FFSTypeHandle fifth_rec_ioformat, sixth_rec_ioformat;
static FFSTypeHandle ninth_rec_ioformat, string_array_ioformat, derive_ioformat;
static FFSTypeHandle multi_array_ioformat, triangle_ioformat, add_action_ioformat;
static FFSTypeHandle node_ioformat;

static void
set_targets(FFSContext context)
{
    if ((test_only == NULL) || (strcmp(test_only, "first_rec") == 0))
	first_rec_ioformat = FFSset_fixed_target(context, first_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "second_rec") == 0))
	second_rec_ioformat = FFSset_fixed_target(context, string_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "third_rec") == 0))
	third_rec_ioformat = FFSset_fixed_target(context, two_string_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "fourth_rec") == 0))
	fourth_rec_ioformat = FFSset_fixed_target(context, fourth_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "fifth_rec") == 0))
	fifth_rec_ioformat = FFSset_fixed_target(context, structured_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "sixth_rec") == 0))
	sixth_rec_ioformat = FFSset_fixed_target(context, variant_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "later_rec") == 0))
	later_rec_ioformat = FFSset_fixed_target(context, later_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "nested_rec") == 0))
	nested_rec_ioformat = FFSset_fixed_target(context, nested_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "ninth_rec") == 0))
	ninth_rec_ioformat = FFSset_fixed_target(context, ninth_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "string_array") == 0))
	string_array_ioformat = FFSset_fixed_target(context, string_array_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "derive") == 0))
	derive_ioformat = FFSset_fixed_target(context, derive_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "multi_array") == 0))
	multi_array_ioformat = FFSset_fixed_target(context, multi_array_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "triangle_param") == 0))
	triangle_ioformat = FFSset_fixed_target(context, triangle_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "add_action") == 0))
	add_action_ioformat = FFSset_fixed_target(context, add_action_format_list);
    if ((test_only == NULL) || (strcmp(test_only, "node") == 0))
	node_ioformat = FFSset_fixed_target(context, node_format_list);
}

size_t base_size_func(FFSContext context, char *src, size_t rec_len,
		   size_t native_struct_size)
{
    return native_struct_size;
}

size_t total_size_func(FFSContext context, char *src, size_t rec_len, 
		    size_t native_struct_size)
{
    return FFS_est_decode_length(context, src, rec_len);
}

static int 
decode_in_place(FFSContext context, char *src, size_t src_len, void *dest)
{
    if (FFSdecode_in_place_possible(FFSTypeHandle_from_encode(context, src))) {
	int ret, header_len;
	char *real_dest;
	ret = FFSdecode_in_place(context, src, (void**)&real_dest);
	header_len = (int)(intptr_t)(real_dest - src);
	memcpy(dest, real_dest, src_len - header_len);
	return ret;
    } else {
	return FFSdecode_to_buffer(context, src, dest);
    }
}

static int
decode_IOcontext_wrapper(FFSContext context, char *src, size_t src_len, void *dest)
{
    return FFSdecode(context, src, dest);
}

static int
decode_to_buffer_IOcontext_wrapper(FFSContext context, char *src, size_t src_len,
				   void *dest)
{
    return FFSdecode_to_buffer(context, src, dest);
}

typedef size_t (*size_func_t)(FFSContext context, char *src, size_t buf_size, 
				   size_t nativ_struct);

typedef int (*decode_func_t)(FFSContext context, char *src, size_t src_len, 
				   void *dest);

size_func_t size_funcs[] = {base_size_func, total_size_func, total_size_func};
decode_func_t decode_funcs[] = {decode_IOcontext_wrapper, 
				decode_to_buffer_IOcontext_wrapper,
				decode_in_place};


#define NUM_TESTS 3

static void
test_all_receive(char *buffer, size_t buf_size, int finished)
{
    int test_type = 0;
    char *tmp_buffer = malloc(buf_size);
    for ( test_type = 0; test_type < NUM_TESTS; test_type++) {
	memcpy(tmp_buffer, buffer, buf_size);
	test_receive(tmp_buffer, buf_size, finished, test_type);
    }
    free(tmp_buffer);
}
	
static void*
get_mem(size_t size)
{
    char *buffer;
    unsigned int beef = 0xdeadbeef;

    buffer = malloc(size+4);
    memcpy(buffer+size, &beef, 4);
    return buffer;
}

static void
check_mem(size_t size, char *buffer)
{
    unsigned int beef = 0xdeadbeef;
    if (memcmp(buffer+size, &beef, 4) != 0) {
	printf("memory overwrite error\n");
    }
}

    
static void
test_receive(char *buffer, size_t buf_size, int finished, int test_level)
{
    static FFSContext c = NULL;
/*    static int comment_count[NUM_TESTS] = {0,0,0};*/
    static int first_rec_count[NUM_TESTS] = {0,0,0};
    static int second_rec_count[NUM_TESTS] = {0,0,0};
    static int third_rec_count[NUM_TESTS] = {0,0,0};
    static int fourth_rec_count[NUM_TESTS] = {0,0,0};
    static int fifth_rec_count[NUM_TESTS] = {0,0,0};
    static int sixth_rec_count[NUM_TESTS] = {0,0,0};
    static int nested_rec_count[NUM_TESTS] = {0,0,0};
    static int later_rec_count[NUM_TESTS] = {0,0,0};
    static int ninth_rec_count[NUM_TESTS] = {0,0,0};
    static int string_array_count[NUM_TESTS] = {0,0,0};

    static int unknown_rec_count[NUM_TESTS] = {0,0,0};
    size_func_t size_func = size_funcs[test_level];
    decode_func_t decode_func = decode_funcs[test_level];
    if (c == NULL) {
	if (loaded_FMcontext != NULL) {
	    c = create_FFSContext_FM(loaded_FMcontext);
	} else {
	    c = create_FFSContext();
	}
	rcv_context = c;
	set_targets(rcv_context);
    }
    if (!finished) {
/*	char *comment;*/
	FFSTypeHandle buffer_format = FFS_target_from_encode(rcv_context, buffer);

	if ((buffer_format == NULL) && (test_only == NULL)) {
	    printf("No format!\n");
	    exit(1);
	}
	if (verbose && (test_only == NULL)) {
	    printf("Testing format %s, test level %d\n",
		   FFSTypeHandle_name(buffer_format), test_level);
	}
	if (((test_only == NULL) || (strcmp(test_only, "first_rec") == 0)) &&
	    (buffer_format == first_rec_ioformat)) {
	    first_rec read_data[10];
/*	    if (get_IOcontext_record_count(iocontext, buffer) == 1) {*/
		memset(&read_data[0], 0, sizeof(first_rec));
		if (!decode_func(rcv_context, buffer, buf_size, &read_data[0]))
		    printf("decode failed, first data\n");
		
		if (!first_rec_eq(&read_data[0], &rec1_array[first_rec_count[test_level]++])) {
		    printf("Rec1 failure\n");
		    fail++;
		}
/*	    } else {
		int count = next_IOrecord_count(iofile);
		if (count == 10) {
		    memset(&read_data[0], 0, sizeof(first_rec) * 10);
		    if (read_array_IOfile(iofile, &read_data[0], 10, sizeof(first_rec)) != 10)
			printf("decode failed, first array\n");
		    if (memcmp(&read_data[0], &rec1_array[first_rec_count[test_level]],
			       sizeof(first_rec) * 10) != 0) {
			printf("Rec1 failure\n");
			fail++;
		    }
		    first_rec_count[test_level] += 10;
		} else {
		    printf("Rec1 failure\n");
		    fail++;
		}
	    }		*/
	} else if (((test_only == NULL) || (strcmp(test_only, "second_rec") == 0)) &&
		   (buffer_format == second_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(second_rec));
	    second_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, second data failed\n");
	    if (!second_rec_eq(read_data, &rec2_array[second_rec_count[test_level]++])) {
		printf("Rec2 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "third_rec") == 0)) &&
		   (buffer_format == third_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(third_rec));
	    third_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, third data failed\n");
	    if (!third_rec_eq(read_data, &rec3_array[third_rec_count[test_level]++])) {
		printf("Rec3 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "fourth_rec") == 0)) &&
		   (buffer_format == fourth_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(fourth_rec));
	    fourth_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, fourth data failed\n");
	    if (!fourth_rec_eq(read_data, &rec4)) {
		printf("Rec4 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	    fourth_rec_count[test_level]++;
	} else if (((test_only == NULL) || (strcmp(test_only, "fifth_rec") == 0)) &&
		   (buffer_format == fifth_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(fifth_rec));
	    fifth_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, fifth data failed\n");
	    if (!fifth_rec_eq(read_data, &rec5)) {
		printf("Rec5 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	    fifth_rec_count[test_level]++;
	} else if (((test_only == NULL) || (strcmp(test_only, "sixth_rec") == 0)) &&
		   (buffer_format == sixth_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(sixth_rec));
	    sixth_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, variant format\n");
	    if (!sixth_rec_eq(read_data, &rec6_array[sixth_rec_count[test_level]++])) {
		printf("sixth_rec failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "nested_rec") == 0)) &&
		   (buffer_format == nested_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(nested_rec));
	    nested_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, variant format\n");
	    if (!nested_rec_eq(read_data, &rec7_array[nested_rec_count[test_level]++])) {
		printf("Rec7 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "later_rec") == 0)) &&
		   (buffer_format == later_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(later_rec));
	    later_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, variant format\n");
	    if (!later_rec_eq(read_data, &rec8_array[later_rec_count[test_level]++])) {
		printf("Rec8 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "ninth_rec") == 0)) &&
		   (buffer_format == ninth_rec_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, sizeof(ninth_rec));
	    ninth_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, variant format\n");
	    if (!ninth_rec_eq(read_data, &rec9_array[ninth_rec_count[test_level]++])) {
		printf("Rec9 failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "string_array") == 0)) &&
		   (buffer_format == string_array_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(string_array_rec));
	    string_array_rec *sread_data = get_mem(size);
	    memset(sread_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, sread_data))
		printf("decode failed, string array format\n");
	    if (!string_array_eq(sread_data, 
				 &string_array_array[string_array_count[test_level]++])) {
		printf("string array failure\n");
		fail++;
	    }
	    check_mem(size, (char*)sread_data);
	    free(sread_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "derive") == 0)) &&
		   (buffer_format == derive_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(DeriveMsg));
	    DeriveMsgPtr read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    if (!derive_eq(read_data, &derive)) {
		printf("derive msg failure, decode func %d\n", test_level);
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "multi_array") == 0)) &&
		   (buffer_format == multi_array_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(multi_array));
	    multi_array_rec *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    if (!multi_array_eq(read_data, &multi_array)) {
		printf("multi array failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "triangle_param") == 0)) &&
		   (buffer_format == triangle_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(multi_array));
	    triangle_param *read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    if (!triangle_param_eq(read_data, &triangle)) {
		printf("triangle param failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "add_action") == 0)) &&
		   (buffer_format == add_action_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(add_rec));
	    add_rec_ptr read_data = get_mem(size);
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    if (!add_rec_eq(read_data, &add_action_record)) {
		printf("add_action failure\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (((test_only == NULL) || (strcmp(test_only, "node") == 0)) &&
		   (buffer_format == node_ioformat)) {
	    size_t size = size_func(rcv_context, buffer, buf_size, 
				 sizeof(add_rec));
	    node_ptr read_data = get_mem(size);
	    struct visit_table v;
	    int signature;
	    memset(read_data, 0, size);
	    if (!decode_func(rcv_context, buffer, buf_size, read_data))
		printf("decode failed, decode msg format\n");
	    signature = read_data->node_num;
	    read_data->node_num = 0;
	    v.node_count = 0;
	    if (signature != calc_signature(read_data, &v)) {
		printf("Node sig not right\n");
		fail++;
	    }
	    check_mem(size, (char*)read_data);
	    free(read_data);
	} else if (test_only == NULL) {
	    printf("discarding a record\n");
	    unknown_rec_count[test_level]++;
	}
    } else {
	/* finished */
	if (test_level == 0) {
	    free_FFSContext(rcv_context);
	    rcv_context = NULL;
	}
	if (test_only != NULL) return;
	if (first_rec_count[test_level] != 3) {
	    printf("Missed first\n");
	    fail++;
	}
	if (second_rec_count[test_level] != sizeof(rec2_array) / sizeof(rec2_array[0])) {
	    printf("Missed second\n");
	    fail++;
	}
	if (third_rec_count[test_level] != sizeof(rec3_array) / sizeof(rec3_array[0])) {
	    printf("Missed third\n");
	    fail++;
	}
	if (fourth_rec_count[test_level] != 1) {
	    printf("Missed fourth\n");
	    fail++;
	}
	if (fifth_rec_count[test_level] != 1) {
	    printf("Missed fifth\n");
	    fail++;
	}
	if (sixth_rec_count[test_level] != sizeof(rec6_array) / sizeof(rec6_array[0])) {
	    printf("Missed sixth\n");
	    fail++;
	}
	if (nested_rec_count[test_level] != sizeof(rec7_array) / sizeof(rec7_array[0])) {
	    printf("Missed nested\n");
	    fail++;
	}
	if (later_rec_count[test_level] != 3) {
	    printf("Missed later\n");
	    fail++;
	}
	if (unknown_rec_count[test_level] != 0) {
	    printf("Got unknown\n");
	    fail++;
	}
    }
}

static FMFormat seen_formats[100];
static int seen_count = 0;

static void
write_buffer(FMFormat format, char *buf, size_t size)
{
    static int file_fd = 0;
    int i;
    int written = 0;
    int indicator;
    if (output_file == NULL) return;

    if (file_fd == 0) {
	file_fd = open(output_file, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY, 0777);
    }
    for (i=0; i < seen_count; i++) {
	if (format == seen_formats[i]) {
	    written++;
	}
    }
    if (!written) {
	struct {
	    int indicator;
	    int format_len;
	} format_header;
	struct FFSEncodeVec vec[4];
	char *server_id;
	int id_len;
	char *server_rep;
	int rep_len;
	if (verbose) dump_FMFormat(format);
	server_id = get_server_ID_FMformat(format, &id_len);
	server_rep = get_server_rep_FMformat(format, &rep_len);

	/*
	 * next_data indicator is a 2 4-byte chunks in network byte order.
	 * In the first chunk, 
	 *    the top byte is 0x2, middle 2 are unused and
	 *    the bottom byte is the size of the format id;
	 * The second chunk holds the length of the format rep;
	 */
	
	format_header.indicator = htonl((id_len & 0xff) | 0x2 << 24);
	format_header.format_len = htonl(rep_len);
	
	vec[0].iov_len = 8;
	vec[0].iov_base = &format_header;
	vec[1].iov_len = id_len;
	vec[1].iov_base = server_id;
	vec[2].iov_len = rep_len;
	vec[2].iov_base = server_rep;
	vec[3].iov_len = 0;
	vec[3].iov_base = NULL;
#ifndef _MSC_VER
	if (writev(file_fd, (const struct iovec *) vec, 3) == -1) {
	    printf("Writev failed\n");
	    exit(1);
	}
#else
	write(file_fd, vec[0].iov_base, (int)vec[0].iov_len);
	write(file_fd, vec[1].iov_base, (int)vec[1].iov_len);
	write(file_fd, vec[2].iov_base, (int)vec[2].iov_len);
#endif

	seen_formats[seen_count++] = format;
    }

    /*
     * next_data indicator is a 4-byte chunk in network byte order.
     * The top byte is 0x3.  The bottom 3 bytes are the size of the data.
     */
    indicator = htonl((size & 0xffffff) | 0x3 << 24);

    if (write(file_fd, &indicator, 4) != 4) {printf("Write 4 failed\n");exit(1);}
    if (write(file_fd, buf, (unsigned int)size) != size) {printf("Write size failed\n");exit(1);}
}
