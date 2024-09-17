
#include "config.h"
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
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
#include "unix_defs.h"

#include "test_funcs.h"

static void test_all_receive ARGS((char *buffer, int buf_size, int finished));
static void write_buffer ARGS((char *buf, int size));
static void read_test_only();

static int write_output = 0;
static char *output_file = NULL;
static char *read_file = NULL;


int
main(argc, argv)
int argc;
char **argv;
{
    IOContext src_context = create_IOcontext(NULL);
    IOFormat first_rec_ioformat, newer_rec_ioformat;
    IOFormat second_rec_ioformat, third_rec_ioformat;
    IOFormat fourth_rec_ioformat, later_ioformat, nested_ioformat;
    IOFormat embedded_rec_ioformat, fifth_rec_ioformat;
    char *xfer_buffer;
    IOOptInfo opt_info[2];
    int buf_size;
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
    int i, j;
    IOFormat sixth_rec_ioformat, ninth_rec_ioformat;
    opt_info[0].info_type = COMPAT_OPT_INFO;
    opt_info[1].info_type = 0;

    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "-w") == 0) {
	    output_file = argv[++i];
	    write_output++;
	} else if (strcmp(argv[i], "-r") == 0) {
	    read_file = argv[++i];
	}
    }

    if (read_file) {
	read_test_only();
	exit(0);
    }
    opt_info[0].info_type = XML_OPT_INFO;
    opt_info[0].info_len = strlen(first_xml) +1;
    opt_info[0].info_block = first_xml;
    first_rec_ioformat = register_opt_format("first format", field_list,
					     opt_info, src_context);
    {
	char * xform_code = 
"{\n\
	old.integer_field = new.ganzzahl;\n\
	old.double_field = new.doppltezahl;\n\
	old.char_field = new.buchstabe;\n\
}\n";
	opt_info[0].info_type = COMPAT_OPT_INFO;
	opt_info[0].info_block = create_compat_info(first_rec_ioformat,
						    xform_code, 
						    &opt_info[0].info_len);
	newer_rec_ioformat = register_opt_format("newer format", 
						 newer_field_list,
						 opt_info, src_context);
    }

    opt_info[0].info_type = XML_OPT_INFO;
    opt_info[0].info_len = strlen(string_xml) +1;
    opt_info[0].info_block = string_xml;
    second_rec_ioformat = register_opt_format("string format", field_list2,
					      opt_info, src_context);
    opt_info[0].info_len = strlen(third_xml) +1;
    opt_info[0].info_block = third_xml;
    third_rec_ioformat = register_opt_format("two string format", field_list3,
					     opt_info, src_context);
    opt_info[0].info_len = strlen(fourth_xml) +1;
    opt_info[0].info_block = fourth_xml;
    fourth_rec_ioformat = register_opt_format("internal array format",
					      field_list4, opt_info,
					      src_context);
    opt_info[0].info_len = strlen(embedded_xml) +1;
    opt_info[0].info_block = embedded_xml;
    embedded_rec_ioformat = register_opt_format("embedded",
						embedded_field_list, opt_info,
						src_context);
    (void)embedded_rec_ioformat;
    opt_info[0].info_len = strlen(struct_array_xml) +1;
    opt_info[0].info_block = struct_array_xml;
    fifth_rec_ioformat = register_opt_format("structured array format",
					     field_list5, opt_info,
					     src_context);
    opt_info[0].info_len = strlen(var_array_xml) +1;
    opt_info[0].info_block = var_array_xml;
    sixth_rec_ioformat = register_opt_format("variant array format",
					     field_list6, opt_info,
					     src_context);
/*    write_comment_IOfile(src_context, "this is a comment in the file");*/
    memset((char *) &rec1, 0, sizeof(rec1));
    rec1.integer_field = 14;
    rec1.double_field = 2.717;
    rec1.char_field = 'A';
    xfer_buffer = encode_IOcontext_buffer(src_context, newer_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
#ifdef NDEF
    write_buffer(xfer_buffer, buf_size);
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
    xfer_buffer = encode_IOcontext_buffer(src_context, fifth_rec_ioformat, &emb_array, &buf_size);
    free(emb_array.earray[0].string);
    free(emb_array.earray[1].string);
    free(emb_array.earray[2].string);
    free(emb_array.earray[3].string);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    memset((char *) &rec2, 0, sizeof(rec2));
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "testing";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = encode_IOcontext_buffer(src_context, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = encode_IOcontext_buffer(src_context, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec1.integer_field = 17;
    rec1.double_field *= 3.0;
    rec1.char_field = 'B';
    xfer_buffer = encode_IOcontext_buffer(src_context, first_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = encode_IOcontext_buffer(src_context, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec1.integer_field *= 2;
    rec1.double_field *= 2.717;
    rec1.char_field = 'C';
/*    write_comment_IOfile(iofile, "this is another comment in the file");*/
    xfer_buffer = encode_IOcontext_buffer(src_context, first_rec_ioformat, &rec1, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
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
    xfer_buffer = encode_IOcontext_buffer(src_context, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    opt_info[0].info_len = strlen(later_xml) +1;
    opt_info[0].info_block = later_xml;
    later_ioformat = register_opt_format("later format",later_field_list,
					 opt_info, src_context);
    opt_info[0].info_len = strlen(nested_xml) +1;
    opt_info[0].info_block = nested_xml;
    nested_ioformat = register_opt_format("nested format", nested_field_list,
					  opt_info, src_context);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Paulaner;
    xfer_buffer = encode_IOcontext_buffer(src_context, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    memset((char *) &rec7, 0, sizeof(rec7));
    rec7.integer_field = 47;
    rec7.nested_rec.integer_field = 14;
    rec7.nested_rec.short_field = 27;
    rec7.nested_rec.long_field = 987234;
    rec7.nested_rec.string = "Another string";
    rec7.nested_rec.double_field = 2.717;
    rec7.nested_rec.char_field = 'A';
    rec7.string = "Yet another string";
    xfer_buffer = encode_IOcontext_buffer(src_context, nested_ioformat, &rec7, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    xfer_buffer = encode_IOcontext_buffer(src_context, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Red_Stripe;
    xfer_buffer = encode_IOcontext_buffer(src_context, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    memset((char *) &rec5, 0, sizeof(rec5));
    rec5.integer_field = 9872346;
    rec5.string = "ABCD";
    rec5.double_field = 3.14159265358797323;
    xfer_buffer = encode_IOcontext_buffer(src_context, later_ioformat, &rec5, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    xfer_buffer = encode_IOcontext_buffer(src_context, third_rec_ioformat, &rec3, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "the end";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    xfer_buffer = encode_IOcontext_buffer(src_context, second_rec_ioformat, &rec2, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    for (i = 0; i < 10; i++) {
	memset((char *) &array1[i], 0, sizeof(array1[i]));
	array1[i].integer_field = 2 * i * i;
	array1[i].double_field = 2.717 * (i * i);
	array1[i].char_field = 'D' + i;
    }
    opt_info[0].info_len = strlen(later2_xml) +1;
    opt_info[0].info_block = later2_xml;
    later_ioformat = register_opt_format("later format", later_field_list2,
					 opt_info, src_context);
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
    xfer_buffer = encode_IOcontext_buffer(src_context, fourth_rec_ioformat, &rec4, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    memset((char *) &rec6, 0, sizeof(rec6));
    rec6.integer_field = 23462346;
    rec6.string = "Efghij";
    rec6.double_field = 3.14159265358797323 * 2.0;
    xfer_buffer = encode_IOcontext_buffer(src_context, later_ioformat, &rec6, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);
    rec6.integer_field = 2346987;
    rec6.string = "Klmn";
    rec6.double_field = 3.14159265358797323 * 3.0;
    xfer_buffer = encode_IOcontext_buffer(src_context, later_ioformat, &rec6, &buf_size);
    test_all_receive(xfer_buffer, buf_size, 0);
    write_buffer(xfer_buffer, buf_size);

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

	xfer_buffer = encode_IOcontext_buffer(src_context, sixth_rec_ioformat, 
					      &var_array, &buf_size);
	test_all_receive(xfer_buffer, buf_size, 0);
	write_buffer(xfer_buffer, buf_size);
	for (j = 0; j < var_array.icount; j++) {
	    free(var_array.var_string_array[j].string);
	}
	free(var_array.var_string_array);
	free(var_array.string);
	free(var_array.var_int_array);
	free(var_array.var_double_array);
    }
    opt_info[0].info_len = strlen(event_xml) +1;
    opt_info[0].info_block = event_xml;
    register_opt_format("EventVecElem", event_vec_elem_fields,
			opt_info, src_context);

    opt_info[0].info_len = strlen(event_vec_xml) +1;
    opt_info[0].info_block = event_vec_xml;
    ninth_rec_ioformat = register_opt_format("EventV", field_list9,
					     opt_info, src_context);
    for (i = 1; i < 10; i += 2) {
        memset((char *) &var_var, 0, sizeof(var_var));
	var_var.vec_length = i;
	var_var.eventv = malloc(sizeof(((ninth_rec_ptr) 0)->eventv[0]) * var_var.vec_length);

	for (j = 0; j < var_var.vec_length; j++) {
	    int k;
	    var_var.eventv[j].iov_len = j + i;
	    var_var.eventv[j].iov_base = malloc(j + i);
	    for (k=0; k<j+i; k++) {
		((char*)var_var.eventv[j].iov_base)[k] = 'A' + k + i/5;
	    }
	}
	xfer_buffer = encode_IOcontext_buffer(src_context, ninth_rec_ioformat, 
					      &var_var, &buf_size);
	test_all_receive(xfer_buffer, buf_size, 0);
	write_buffer(xfer_buffer, buf_size);
	for (j = 0; j < var_var.vec_length; j++) {
	    free(var_var.eventv[j].iov_base);
	}
	free(var_var.eventv);
    }
#endif
    free_IOcontext(src_context);
    write_buffer(NULL, 0);
    return 0;
}

/* NT needs O_BINARY, but it doesn't exist elsewhere */
#ifndef O_BINARY
#define O_BINARY 0
#endif

static char *
get_buffer(size_p)
int *size_p;
{
    static int file_fd = 0;
    static char *buffer = NULL;
    char *tmp_buffer;
    static int last_size = -1;
    unsigned short ssize;
    int to_read;
    int tmp_size;
    unsigned char csize;
    unsigned int beef = 0xdeadbeef;

    if (read_file == NULL) exit(1);

    if (file_fd == 0) {
	file_fd = open(read_file, O_RDONLY|O_BINARY, 0777);
	buffer = malloc(1);
    }
    if (last_size != -1) {
	if (memcmp(buffer+last_size, &beef, 4) != 0) {
	    printf("memory overwrite error\n");
	}
    }
    read(file_fd, &csize, 1);	/* low byte of 2-byte size */
    ssize = csize;
    read(file_fd, &csize, 1);	/* high byte of 2-byte size */
    ssize += ((csize << 8) & 0xff00);
    to_read = ssize;
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
    last_size = ssize;
    memcpy(buffer+last_size, &beef, 4);
    if (ssize == 0) {
	free(buffer);
	close(file_fd);
	file_fd = 0;
	return NULL;
    } else {
	*size_p = ssize;
	return buffer;
    }
}

static void
read_test_only()
{
    char *input;
    int size;
    while ((input = get_buffer(&size)) != NULL) {
	test_all_receive(input, size, 0);
    }
}


static void
test_all_receive(buffer, buf_size, finished)
char *buffer;
int buf_size;
int finished;
{
    IOContext iocontext = create_IOcontext();
    char *str;
    IOFormat format = get_format_IOcontext(iocontext, buffer);
    IOcompat_formats compats = get_compat_formats(format);
    printf("\nCurrent record format is -> \n");
    dump_IOFormat(format);
    if (compats != NULL) {
	int i = 0;
	printf("\n Record has backwards compatability information!\n\n");
	while (compats[i].prior_format != NULL) {
	    printf("   BC format/code group %d\n", i);
	    printf("\tTransform code is \n");
	    printf("\t##### %s #####\n", compats[i].xform_code);
	    printf("\tTransform creates data in format -> \n");
	    dump_IOFormat(compats[i].prior_format);
	    i++;
	}
    }
    printf("\n current data in XML is\n");
    str = IOencoded_to_XML_string(iocontext, buffer);
    printf("%s\n==================================================\n\n", str);
    free(str);
    free_IOcontext(iocontext);
}
	
static void
write_buffer(buf, size)
char *buf;
int size;
{
    static int file_fd = 0;
    unsigned short ssize;
    unsigned char csize;
    if (output_file == NULL) return;

    if (file_fd == 0) {
	file_fd = open(output_file, O_WRONLY|O_CREAT|O_TRUNC|O_BINARY, 0777);
    }
    ssize = size;
    csize = ssize & 0xff;
    write(file_fd, &csize, 1);	/* low byte of 2-byte size */
    csize = ((ssize >> 8) & 0xff);
    write(file_fd, &csize, 1);	/* high byte of 2-byte size */
    write(file_fd, buf, size);
}
