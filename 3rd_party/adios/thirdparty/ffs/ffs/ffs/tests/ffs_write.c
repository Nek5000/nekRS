
#include "config.h"
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdio.h>
#include <string.h>
#include "ffs.h"

#include "test_funcs.h"
#define UPDATE_ATTR(x)     if (x) {\
	int tmp;\
	get_int_attr(x, iteration_atom, &tmp);	\
	tmp++;\
	set_int_attr(x, iteration_atom, tmp);\
    }\
	

int
main(int argc, char **argv)
{
    FMContext src_context;
    FFSFile ffsfile;
    FMFormat first_rec_ioformat, second_rec_ioformat, third_rec_ioformat;
    FMFormat fourth_rec_ioformat, later_ioformat, nested_ioformat;
    FMFormat embedded_rec_ioformat, fifth_rec_ioformat, sixth_rec_ioformat;
    FMFormat ninth_rec_ioformat;
    FMFormat string_array_ioformat;
    first_rec rec1;
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
    FMStructDescRec str_list[5];
    FMOptInfo opt_info[2], opt_info2[2];
    int verbose = 0, indexed = 0, attributes = 0, reopen = 0;
    char *output_file = "test_output";

    atom_t level_atom = -1, iteration_atom = -1;
    level_atom = attr_atom_from_string("level");
    iteration_atom = attr_atom_from_string("iteration");
    attr_list a1 = NULL;
    attr_list a2 = NULL;

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
	    /* argument */
	    if (argv[i][1] == 'v') {
		verbose++;
	    } else if (argv[i][1] == 'i') {
		indexed++;
	    } else if (argv[i][1] == 'a') {
		attributes++;
		indexed++;
	    } else if (argv[i][1] == 'r') {
		reopen++;
	    } else {
		printf("Unknown argument \"%s\"\n", argv[i]);
	    }
	} else {
	    output_file = argv[i];
	}
    }

    if (attributes) {
	a1 = create_attr_list();
	a2 = create_attr_list();
	set_int_attr(a1, level_atom, 5);
	set_int_attr(a1, iteration_atom, 1);
	set_int_attr(a2, iteration_atom, 4);
    }
    if (!indexed) {
	ffsfile = open_FFSfile(output_file, "w");
    } else {
	ffsfile = open_FFSfile(output_file, "wi");
    }

    
    src_context = FMContext_of_file(ffsfile);

    opt_info[0].info_type = 0x584D4C20;   /* XML */
    opt_info[1].info_type = 0;
    opt_info[1].info_len = 0;
    opt_info[1].info_block = NULL;
    str_list[0].format_name = "first format";
    str_list[0].field_list = field_list;
    str_list[0].struct_size = sizeof(first_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;

    opt_info[0].info_len = (int)strlen(first_xml) +1;
    opt_info[0].info_block = first_xml;
    first_rec_ioformat = register_data_format(src_context, str_list);
    str_list[0].format_name = "string format";
    str_list[0].field_list = field_list2;
    str_list[0].struct_size = sizeof(second_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(string_xml) +1;
    opt_info[0].info_block = string_xml;
    second_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "two string format";
    str_list[0].field_list = field_list3;
    str_list[0].struct_size = sizeof(third_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(third_xml) +1;
    opt_info[0].info_block = third_xml;
    third_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "internal array format";
    str_list[0].field_list = field_list4;
    str_list[0].struct_size = sizeof(fourth_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(fourth_xml) +1;
    opt_info[0].info_block = fourth_xml;
    fourth_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "embedded";
    str_list[0].field_list = embedded_field_list;
    str_list[0].struct_size = sizeof(embedded_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(embedded_xml) +1;
    opt_info[0].info_block = embedded_xml;
    embedded_rec_ioformat = register_data_format(src_context, str_list);


    fifth_rec_ioformat = register_data_format(src_context, structured_format_list);

    str_list[0].format_name = "variant array format";
    str_list[0].field_list = field_list6;
    str_list[0].struct_size = sizeof(sixth_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = "string format";
    str_list[1].field_list = field_list2;
    str_list[1].struct_size = sizeof(second_rec);
    str_list[1].opt_info = &opt_info2[0];
    str_list[2].format_name = NULL;
    opt_info[0].info_len = (int)strlen(var_array_xml) +1;
    opt_info[0].info_block = var_array_xml;
    opt_info2[0].info_type = 0x584D4C20;   /* XML */
    opt_info2[0].info_len = (int)strlen(string_xml) +1;
    opt_info2[0].info_block = string_xml;
    opt_info2[1].info_type = 0;
    opt_info2[1].info_len = 0;
    opt_info2[1].info_block = NULL;
    sixth_rec_ioformat = register_data_format(src_context, str_list);

    memset((char *) &rec1, 0, sizeof(rec1));
    rec1.integer_field = 14;
    rec1.double_field = 2.717;
    rec1.char_field = 'A';
    if (!write_FFSfile_attrs(ffsfile, first_rec_ioformat, &rec1, a1))
	printf("write failed\n");
    UPDATE_ATTR(a1); /* new 2 */
    memset((char *) &emb_array, 0, sizeof(emb_array));
    memset((char *) &var_array, 0, sizeof(var_array));
    emb_array.earray[0].dfield = 4.0;
    emb_array.earray[0].ifield = 4;
    emb_array.earray[0].string = (char *) malloc(10);
    sprintf(emb_array.earray[0].string, "string%d", emb_array.earray[0].ifield * 5);
    emb_array.earray[1].dfield = 3.0;
    emb_array.earray[1].ifield = 3;
    emb_array.earray[1].string = (char *) malloc(10);
    sprintf(emb_array.earray[1].string, "string%d", emb_array.earray[1].ifield * 5);
    emb_array.earray[2].dfield = 2.0;
    emb_array.earray[2].ifield = 2;
    emb_array.earray[2].string = (char *) malloc(10);
    sprintf(emb_array.earray[2].string, "string%d", emb_array.earray[2].ifield * 5);
    emb_array.earray[3].dfield = 1.0;
    emb_array.earray[3].ifield = 1;
    emb_array.earray[3].string = (char *) malloc(10);
    sprintf(emb_array.earray[3].string, "string%d", emb_array.earray[3].ifield * 5);
    if (attributes) {
	if (!write_FFSfile_attrs(ffsfile, fifth_rec_ioformat, &emb_array, a1))
	    printf("write failed\n");
	UPDATE_ATTR(a1);  /* new 3 */
    } else {
	if (!write_FFSfile(ffsfile, fifth_rec_ioformat, &emb_array))
	    printf("write failed\n");
    }
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
    if (!write_FFSfile_attrs(ffsfile, second_rec_ioformat, &rec2, a2))  /* attr = 4 */
	printf("write failed\n");
    UPDATE_ATTR(a2); 
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    if (!write_FFSfile_attrs(ffsfile, second_rec_ioformat, &rec2, a1))
	printf("write failed\n");
    UPDATE_ATTR(a1); /* new 4 */
    rec1.integer_field = 17;
    rec1.double_field *= 3.0;
    rec1.char_field = 'B';
    write_FFSfile_attrs(ffsfile, first_rec_ioformat, &rec1, a1);
    UPDATE_ATTR(a1); /* new 5 */
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = NULL;
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    if (!write_FFSfile_attrs(ffsfile, second_rec_ioformat, &rec2, a1))
	printf("write failed\n");
    UPDATE_ATTR(a1); /* new 6 */
    rec1.integer_field *= 2;
    rec1.double_field *= 2.717;
    rec1.char_field = 'C';
    write_comment_FFSfile(ffsfile, "this is another comment in the file");
    write_FFSfile_attrs(ffsfile, first_rec_ioformat, &rec1, a1);
    UPDATE_ATTR(a1); /* new 7 */
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
    if (!write_FFSfile_attrs(ffsfile, third_rec_ioformat, &rec3, a1))
	printf("write failed\n");

    if (reopen) {
	close_FFSfile(ffsfile);
	free_FFSfile(ffsfile);
	if (!indexed) {
	    ffsfile = open_FFSfile(output_file, "a");
	} else {
	    ffsfile = open_FFSfile(output_file, "ai");
	}
	src_context = FMContext_of_file(ffsfile);
	opt_info[0].info_type = 0x584D4C20;   /* XML */
	opt_info[1].info_type = 0;
	opt_info[1].info_len = 0;
	opt_info[1].info_block = NULL;
	str_list[0].format_name = "first format";
	str_list[0].field_list = field_list;
	str_list[0].struct_size = sizeof(first_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = NULL;

	opt_info[0].info_len = (int)strlen(first_xml) +1;
	opt_info[0].info_block = first_xml;
	first_rec_ioformat = register_data_format(src_context, str_list);


	str_list[0].format_name = "string format";
	str_list[0].field_list = field_list2;
	str_list[0].struct_size = sizeof(second_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = NULL;
	opt_info[0].info_len = (int)strlen(string_xml) +1;
	opt_info[0].info_block = string_xml;
	second_rec_ioformat = register_data_format(src_context, str_list);
	
	str_list[0].format_name = "two string format";
	str_list[0].field_list = field_list3;
	str_list[0].struct_size = sizeof(third_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = NULL;
	opt_info[0].info_len = (int)strlen(third_xml) +1;
	opt_info[0].info_block = third_xml;
	third_rec_ioformat = register_data_format(src_context, str_list);
	
	str_list[0].format_name = "internal array format";
	str_list[0].field_list = field_list4;
	str_list[0].struct_size = sizeof(fourth_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = NULL;
	opt_info[0].info_len = (int)strlen(fourth_xml) +1;
	opt_info[0].info_block = fourth_xml;
	fourth_rec_ioformat = register_data_format(src_context, str_list);

	str_list[0].format_name = "embedded";
	str_list[0].field_list = embedded_field_list;
	str_list[0].struct_size = sizeof(embedded_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = NULL;
	opt_info[0].info_len = (int)strlen(embedded_xml) + 1;
	opt_info[0].info_block = embedded_xml;
	embedded_rec_ioformat = register_data_format(src_context, str_list);
        (void)embedded_rec_ioformat;
	
	fifth_rec_ioformat = register_data_format(src_context, structured_format_list);
	
	str_list[0].format_name = "variant array format";
	str_list[0].field_list = field_list6;
	str_list[0].struct_size = sizeof(sixth_rec);
	str_list[0].opt_info = &opt_info[0];
	str_list[1].format_name = "string format";
	str_list[1].field_list = field_list2;
	str_list[1].struct_size = sizeof(second_rec);
	str_list[1].opt_info = &opt_info2[0];
	str_list[2].format_name = NULL;
	opt_info[0].info_len = (int)strlen(var_array_xml) +1;
	opt_info[0].info_block = var_array_xml;
	opt_info2[0].info_type = 0x584D4C20;   /* XML */
	opt_info2[0].info_len = (int)strlen(string_xml) +1;
	opt_info2[0].info_block = string_xml;
	opt_info2[1].info_type = 0;
	opt_info2[1].info_len = 0;
	opt_info2[1].info_block = NULL;
	sixth_rec_ioformat = register_data_format(src_context, str_list);
    }
    
    UPDATE_ATTR(a1); /* new 8 */
    str_list[0].format_name = "later format";
    str_list[0].field_list = later_field_list;
    str_list[0].struct_size = sizeof(later_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(later_xml) +1;
    opt_info[0].info_block = later_xml;
    later_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "nested format";
    str_list[0].field_list = nested_field_list;
    str_list[0].struct_size = sizeof(nested_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = "string format";
    str_list[1].field_list = field_list2;
    str_list[1].struct_size = sizeof(second_rec);
    str_list[1].opt_info = &opt_info2[0];
    str_list[2].format_name = NULL;
    opt_info[0].info_len = (int)strlen(nested_xml) +1;
    opt_info[0].info_block = nested_xml;
    nested_ioformat = register_data_format(src_context, str_list);

    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Paulaner;
    if (!write_FFSfile_attrs(ffsfile, third_rec_ioformat, &rec3, a1))
	printf("write failed\n");
    UPDATE_ATTR(a1); /* new 9 */
    memset((char *) &rec7, 0, sizeof(rec7));
    rec7.integer_field = 47;
    rec7.nested_rec.integer_field = 14;
    rec7.nested_rec.short_field = 27;
    rec7.nested_rec.long_field = 987234;
    rec7.nested_rec.string = "Another string";
    rec7.nested_rec.double_field = 2.717;
    rec7.nested_rec.char_field = 'A';
    rec7.string = "Yet another string";
    if (!write_FFSfile(ffsfile, nested_ioformat, &rec7))
	printf("write failed\n");
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    if (!write_FFSfile(ffsfile, third_rec_ioformat, &rec3))
	printf("write failed\n");
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.double_field = 2.717;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    rec3.enum_field = Red_Stripe;
    if (!write_FFSfile(ffsfile, third_rec_ioformat, &rec3))
	printf("write failed\n");
    memset((char *) &rec5, 0, sizeof(rec5));
    rec5.integer_field = 9872346;
    rec5.string = "ABCD";
    rec5.double_field = 3.14159265358797323;
    if (!write_FFSfile(ffsfile, later_ioformat, &rec5))
	printf("write failed\n");

    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.double_field = 2.717;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    rec3.enum_field = Pilsner;
    if (!write_FFSfile_attrs(ffsfile, third_rec_ioformat, &rec3, a1))
	printf("write failed\n");

    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "the end";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    if (!write_FFSfile_attrs(ffsfile, second_rec_ioformat, &rec2, a1))
	printf("writev failed\n");

    str_list[0].format_name = "later format";
    str_list[0].field_list = later_field_list2;
    str_list[0].struct_size = sizeof(later_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = (int)strlen(later2_xml) +1;
    opt_info[0].info_block = later2_xml;

    later_ioformat = register_data_format(src_context, str_list);

    memset((char *) &rec4, 0, sizeof(rec4));
    for (i = 0; i < ARRAY_SIZE; i++) {
	rec4.int_array[i] = 297 + i;
    }
    rec4.double_array[0][0] = 1.0;
    rec4.double_array[0][1] = 2.0;
    rec4.double_array[1][0] = 3.0;
    rec4.double_array[1][1] = 4.0;
    rec4.ifield = -rec4.int_array[ARRAY_SIZE - 1];
    if (!write_FFSfile(ffsfile, fourth_rec_ioformat, &rec4))
	printf("writev failed\n");

    memset((char *) &rec6, 0, sizeof(rec6));
    rec6.integer_field = 23462346;
    rec6.string = "Efghij";
    rec6.double_field = 3.14159265358797323 * 2.0;
    if (!write_FFSfile(ffsfile, later_ioformat, &rec6))
	printf("writev failed\n");

    rec6.integer_field = 2346987;
    rec6.string = "Klmn";
    rec6.double_field = 3.14159265358797323 * 3.0;
    
    if (!write_FFSfile(ffsfile, later_ioformat, &rec6))
	printf("write failed\n");

    for (i = 1; i < 20; i += 5) {
        memset((char *) &var_array, 0, sizeof(var_array));
	var_array.string = malloc(10);
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
	    var_array.var_string_array[j].string = malloc(20);
	    sprintf(var_array.var_string_array[j].string,
		    "substring%d", j);
	    var_array.var_string_array[j].double_field = 3.1415 * j;
	    var_array.var_string_array[j].char_field = 'a' + 2 * j;
	}
	if (!write_FFSfile(ffsfile, sixth_rec_ioformat, &var_array))
	    printf("write failed");
	free(var_array.string);
	free(var_array.var_int_array);
	free(var_array.var_double_array);
	for (j = 0; j < var_array.icount; j++) {
	    free(var_array.var_string_array[j].string);
	}
	free(var_array.var_string_array);
    }

    str_list[0].format_name = "EventV";
    str_list[0].field_list = field_list9;
    str_list[0].struct_size = sizeof(ninth_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = "EventVecElem";
    str_list[1].field_list = event_vec_elem_fields;
    str_list[1].struct_size = sizeof(struct _io_encode_vec);
    str_list[1].opt_info = &opt_info2[0];
    str_list[2].format_name = NULL;
    opt_info[0].info_len = (int)strlen(event_xml) +1;
    opt_info[0].info_block = event_xml;

    opt_info2[0].info_len = (int)strlen(event_vec_xml) +1;
    opt_info2[0].info_block = event_vec_xml;
    ninth_rec_ioformat = register_data_format(src_context, str_list);

    str_list[0].format_name = "string_array";
    str_list[0].field_list = string_array_field_list;
    str_list[0].struct_size = sizeof(string_array_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    string_array_ioformat = register_data_format(src_context, str_list);

    for (i = 1; i < 10; i += 2) {
        memset((char *) &var_var, 0, sizeof(var_var));
	memset((char *) &str_array, 0, sizeof(str_array));
	var_var.vec_length = i;
	var_var.eventv = malloc(sizeof(((ninth_rec_ptr) 0)->eventv[0]) * var_var.vec_length);

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
	if ((i % 4) == 1) {
	    str_array.base_string = NULL;
	} else {
	    str_array.base_string = strdup("Whoa there!");
	}
	if (!write_FFSfile(ffsfile, ninth_rec_ioformat, &var_var))
	    printf("write failed");
	if (!write_FFSfile(ffsfile, string_array_ioformat, &str_array))
	    printf("write failed");
	if (str_array.base_string) free(str_array.base_string);
	for (j = 0; j < var_var.vec_length; j++) {
	    free(var_var.eventv[j].iov_base);
	    free(str_array.array[j]);
	}
	free(var_var.eventv);
	free(str_array.array);
    }

    close_FFSfile(ffsfile);
    free_FFSfile(ffsfile);
    free_attr_list(a1);
    free_attr_list(a2);
    return 0;
}
