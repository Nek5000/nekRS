#include "config.h"
#include <assert.h>
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
#include <string.h>
#include "ffs.h"
#include "unix_defs.h"

#include "test_funcs.h"

static int verbose = 0;
static void init_def_written_data();

int
size_func_sizeof(iofile, size)
IOFile iofile;
int size;
{
    return size;
}

int
size_func_next_size(iofile, size)
IOFile iofile;
int size;
{
    return next_IOrecord_length(iofile);
}

int
read_func_no_buffer(iofile, data, size)
IOFile iofile;
void *data;
int size;
{
    return read_IOfile(iofile, data);
}

int
read_func_buffer(iofile, data, size)
IOFile iofile;
void *data;
int size;
{
    return read_to_buffer_IOfile(iofile, data, size);
}


void
do_test(input_file, size_func, read_func)
char *input_file;
int (*size_func) (IOFile, int);
int (*read_func) ();
{
    IOFile iofile;
    IOFormat first_rec_ioformat = NULL, second_rec_ioformat = NULL;
    IOFormat third_rec_ioformat = NULL, fourth_rec_ioformat = NULL;
    IOFormat later_rec_ioformat = NULL, nested_rec_ioformat = NULL;
    IOFormat embedded_rec_ioformat = NULL, fifth_rec_ioformat = NULL;
    IOFormat sixth_rec_ioformat = NULL, ninth_rec_ioformat = NULL;
    IOFormat string_array_ioformat = NULL;
    int finished = 0;
    int comment_count = 0;
    int first_rec_count = 0;
    int second_rec_count = 0;
    int third_rec_count = 0;
    int fourth_rec_count = 0;
    int fifth_rec_count = 0;
    int sixth_rec_count = 0;
    int nested_rec_count = 0;
    int later_rec_count = 0;
    int ninth_rec_count = 0;
    int string_array_rec_count = 0;

    int unknown_rec_count = 0;
    iofile = open_IOfile(input_file, "r");

    init_def_written_data();

    while (!finished) {
	char *comment;
	IOFormat new_format;
	char *format_name;

	switch (next_IOrecord_type(iofile)) {
	case IOcomment:
	    comment = read_comment_IOfile(iofile);
	    if (strcmp(comment, comment_array[comment_count++]) != 0) {
		printf("From IOfile -> \"%s\"\n", comment);
		exit(1);
	    }
	    break;
	case IOformat:
	    new_format = read_format_IOfile(iofile);
	    assert(new_format != NULL);
	    format_name = name_of_IOformat(new_format);
	    if (strcmp(format_name, "first format") == 0) {
		first_rec_ioformat = new_format;
		set_IOconversion(iofile, "first format", field_list,
				 sizeof(first_rec));
	    } else if (strcmp(format_name, "string format") == 0) {
		second_rec_ioformat = new_format;
		set_IOconversion(iofile, "string format", field_list2,
				 sizeof(second_rec));
	    } else if (strcmp(format_name, "two string format") == 0) {
		third_rec_ioformat = new_format;
		set_IOconversion(iofile, "two string format", field_list3,
				 sizeof(third_rec));
	    } else if (strcmp(format_name, "internal array format") == 0) {
		fourth_rec_ioformat = new_format;
		set_IOconversion(iofile, "internal array format", field_list4,
				 sizeof(fourth_rec));
	    } else if (strcmp(format_name, "embedded") == 0) {
		embedded_rec_ioformat = new_format;
		set_IOconversion(iofile, "embedded", embedded_field_list,
				 sizeof(embedded_rec));
	    } else if (strcmp(format_name, "structured array format") == 0) {
		fifth_rec_ioformat = new_format;
		set_IOconversion(iofile, "structured array format", field_list5,
				 sizeof(fifth_rec));
	    } else if (strcmp(format_name, "variant array format") == 0) {
		sixth_rec_ioformat = new_format;
		set_IOconversion(iofile, "variant array format", field_list6,
				 sizeof(sixth_rec));
	    } else if (strcmp(format_name, "later format") == 0) {
		later_rec_ioformat = new_format;
		set_IOconversion(iofile, "later format",
				 later_field_list,
				 sizeof(later_rec));
	    } else if (strcmp(format_name, "nested format") == 0) {
		nested_rec_ioformat = new_format;
		set_IOconversion(iofile, "nested format", nested_field_list,
				 sizeof(nested_rec));
	    } else if (strcmp(format_name, "EventVecElem") == 0) {
		set_IOconversion(iofile, "EventVecElem", event_vec_elem_fields,
				 sizeof(struct _io_encode_vec));
	    } else if (strcmp(format_name, "EventV") == 0) {
		ninth_rec_ioformat = new_format;
		set_IOconversion(iofile, "EventV", field_list9,
				 sizeof(ninth_rec));
	    } else if (strcmp(format_name, "string_array") == 0) {
		string_array_ioformat = new_format;
		set_IOconversion(iofile, "string_array", 
				 string_array_field_list,
				 sizeof(string_array_rec));
	    } else {
		printf("Got unexpected format %s\n", format_name);
		exit(1);
	    }
	    break;
	case IOdata:
	    if (verbose) 
		printf("Record of format %s\n",
		       name_of_IOformat(next_IOrecord_format(iofile)));
	    if (next_IOrecord_format(iofile) == first_rec_ioformat) {
		first_rec read_data[10];
		if (next_IOrecord_count(iofile) == 1) {
		    memset(&read_data[0], 0, sizeof(first_rec));
		    if (!read_IOfile(iofile, &read_data[0]))
			IOperror(iofile, "read first data");
		    if (memcmp(&read_data[0], &rec1_array[first_rec_count++],
			       sizeof(first_rec)) != 0) {
			if (read_data[0].integer_field != 
			    rec1_array[first_rec_count-1].integer_field) {
			    printf("read integer = %d (0x%x), exp integer = %d (0x%x)\n",
				   read_data[0].integer_field, 
				   read_data[0].integer_field, 
				   rec1_array[first_rec_count-1].integer_field,
				   rec1_array[first_rec_count-1].integer_field);
			}
			if (read_data[0].double_field != 
			    rec1_array[first_rec_count-1].double_field) {
			    printf("read double = %g, exp double = %g\n",
				   read_data[0].double_field, 
				   rec1_array[first_rec_count-1].double_field);
			}
			if (read_data[0].char_field != 
			    rec1_array[first_rec_count-1].char_field) {
			    printf("read char = %c, exp char = %c\n",
				   read_data[0].char_field, 
				   rec1_array[first_rec_count-1].char_field);
			}
			printf("Rec1 failure\n");
			exit(1);
		    }
		} else {
		    int count = next_IOrecord_count(iofile);
		    if (count == 10) {
			memset(&read_data[0], 0, sizeof(first_rec) * 10);
			if (read_array_IOfile(iofile, &read_data[0], 10, sizeof(first_rec)) != 10)
			    IOperror(iofile, "read first array");
			if (memcmp(&read_data[0], &rec1_array[first_rec_count],
				   sizeof(first_rec) * 10) != 0) {
			    printf("Rec1 array failure\n");
			    exit(1);
			}
			first_rec_count += 10;
		    } else {
			printf("Rec1 failure\n");
			exit(1);
		    }
		}
	    } else if (next_IOrecord_format(iofile) == second_rec_ioformat) {
		int size = size_func(iofile, sizeof(second_rec));
		second_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read second data failed");
		if (!second_rec_eq(read_data, &rec2_array[second_rec_count++])) {
		    printf("Rec2 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == third_rec_ioformat) {
		int size = size_func(iofile, sizeof(third_rec));
		third_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read third data failed");
		if (!third_rec_eq(read_data, &rec3_array[third_rec_count++])) {
		    printf("Rec3 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == fourth_rec_ioformat) {
		int size = size_func(iofile, sizeof(fourth_rec));
		fourth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read fourth data failed");
		if (!fourth_rec_eq(read_data, &rec4)) {
		    printf("Rec4 failure\n");
		    exit(1);
		}
		free(read_data);
		fourth_rec_count++;
	    } else if (next_IOrecord_format(iofile) == embedded_rec_ioformat) {

		printf("Emb Rec failure\n");
		exit(1);
	    } else if (next_IOrecord_format(iofile) == fifth_rec_ioformat) {
		int size = size_func(iofile, sizeof(fifth_rec));
		fifth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read fifth data failed");
		if (!fifth_rec_eq(read_data, &rec5)) {
		    printf("Rec5 failure\n");
		    exit(1);
		}
		free(read_data);
		fifth_rec_count++;
	    } else if (next_IOrecord_format(iofile) == sixth_rec_ioformat) {
		int size = size_func(iofile, sizeof(sixth_rec));
		sixth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read variant format");
		if (!sixth_rec_eq(read_data, &rec6_array[sixth_rec_count++])) {
		    printf("Rec6 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == nested_rec_ioformat) {
		int size = size_func(iofile, sizeof(nested_rec));
		nested_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read variant format");
		if (!nested_rec_eq(read_data, &rec7_array[nested_rec_count++])) {
		    printf("Rec7 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == later_rec_ioformat) {
		int size = size_func(iofile, sizeof(later_rec));
		later_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read variant format");
		if (!later_rec_eq(read_data, &rec8_array[later_rec_count++])) {
		    printf("Rec8 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == ninth_rec_ioformat) {
		int size = size_func(iofile, sizeof(ninth_rec));
		ninth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read variant format");
		if (!ninth_rec_eq(read_data, &rec9_array[ninth_rec_count++])) {
		    printf("Rec9 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (next_IOrecord_format(iofile) == string_array_ioformat) {
		int size = size_func(iofile, sizeof(string_array_rec));
		string_array_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    IOperror(iofile, "read variant format");
		if (!string_array_eq(read_data, &string_array_array[string_array_rec_count++])) {
		    printf("string_array failure\n");
		    exit(1);
		}
		free(read_data);

	    } else {
		printf("discarding a record\n");
		read_raw_IOfile(iofile, NULL, 0, NULL);
		unknown_rec_count++;
	    }
	    break;
	case IOerror:
	case IOend:
	    finished++;
	    break;
	}
    }

    close_IOfile(iofile);
    free_IOfile(iofile);
    free_written_data();
    if (first_rec_count != sizeof(rec1_array) / sizeof(rec1_array[0])) {
	printf("Missed first\n");
	exit(1);
    }
    if (second_rec_count != sizeof(rec2_array) / sizeof(rec2_array[0])) {
	printf("Missed second\n");
	exit(1);
    }
    if (third_rec_count != sizeof(rec3_array) / sizeof(rec3_array[0])) {
	printf("Missed third\n");
	exit(1);
    }
    if (fourth_rec_count != 1) {
	printf("Missed fourth\n");
	exit(1);
    }
    if (fifth_rec_count != 1) {
	printf("Missed fifth\n");
	exit(1);
    }
    if (sixth_rec_count != sizeof(rec6_array) / sizeof(rec6_array[0])) {
	printf("Missed sixth\n");
	exit(1);
    }
    if (nested_rec_count != sizeof(rec7_array) / sizeof(rec7_array[0])) {
	printf("Missed sixth\n");
	exit(1);
    }
    if (later_rec_count != 3) {
	printf("Missed later\n");
	exit(1);
    }
    if (unknown_rec_count != 0) {
	printf("Got unknown\n");
	exit(1);
    }
    if (ninth_rec_count != 5) {
	printf("missed rec9\n");
	exit(1);
    }
    if ((string_array_rec_count != 5) && (string_array_rec_count != 0)){
	printf("missed string_array\n");
	exit(1);
    }
}


int
main(argc, argv)
int argc;
char **argv;
{
    int i;
    char *input_file = "def_test_output";

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
	    /* argument */
	    if (argv[i][1] == 'v') {
		verbose++;
	    } else {
		printf("Unknown argument \"%s\"\n", argv[i]);
	    }
	} else {
	    input_file = argv[i];
	}
    }
    do_test(input_file, size_func_sizeof, read_func_no_buffer);
    do_test(input_file, size_func_next_size, read_func_buffer);
    return 0;
}

static void
init_def_written_data()
{
    int i, j, k, l, index;
    memset((char *) &rec1_array, 0, sizeof(rec1_array));
    rec1_array[0].integer_field = 14;
    rec1_array[0].double_field = 2.717;
    rec1_array[0].char_field = 'A';
    rec1_array[1].integer_field = 14;
    rec1_array[1].double_field = rec1_array[0].double_field * 3.0;
    rec1_array[1].char_field = 'A';
    rec1_array[2].integer_field = 14;
    rec1_array[2].double_field = rec1_array[1].double_field * 2.717;
    rec1_array[2].char_field = 'A';
    for (i = 0; i < 10; i++) {
	rec1_array[i + 3].integer_field = 14;
	rec1_array[i + 3].double_field = 2.717 * (i * i);
	rec1_array[i + 3].char_field = 'A';
    }

    memset((char *) &rec2_array[0], 0, sizeof(rec2_array));
    rec2_array[0].integer_field = -14;
    rec2_array[0].short_field = 27;
    rec2_array[0].long_field = -9876129;
    rec2_array[0].string = "testing";
    rec2_array[0].double_field = 2.717;
    rec2_array[0].char_field = 'A';

    rec2_array[1].integer_field = -14;
    rec2_array[1].short_field = 27;
    rec2_array[1].long_field = -9876129;
    rec2_array[1].string = NULL;
    rec2_array[1].double_field = 2.717;
    rec2_array[1].char_field = 'A';

    rec2_array[2].integer_field = -14;
    rec2_array[2].short_field = 27;
    rec2_array[2].long_field = -9876129;
    rec2_array[2].string = NULL;
    rec2_array[2].double_field = 2.717;
    rec2_array[2].char_field = 'A';

    rec2_array[3].integer_field = -14;
    rec2_array[3].short_field = 27;
    rec2_array[3].long_field = -9876129;
    rec2_array[3].string = "the end";
    rec2_array[3].double_field = 2.717;
    rec2_array[3].char_field = 'A';

    memset((char *) &rec3_array[0], 0, sizeof(rec3_array));
    rec3_array[0].integer_field = 14;
    rec3_array[0].long_field = 987234;
    rec3_array[0].uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3_array[0].ulong_field = 0xf7e589ceec9dd130;
#else
    rec3_array[0].ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3_array[0].string = "testing";
    rec3_array[0].double_field = 2.717;
    rec3_array[0].string2 = "jambalaya";
    rec3_array[0].char_field = 'A';
    rec3_array[0].enum_field = Paulaner;

    rec3_array[1].integer_field = 14;
    rec3_array[1].long_field = 987234;
    rec3_array[1].uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3_array[1].ulong_field = 0xf7e589ceec9dd130;
#else
    rec3_array[1].ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3_array[1].string = NULL;
    rec3_array[1].double_field = 2.717;
    rec3_array[1].string2 = "jambalaya";
    rec3_array[1].char_field = 'A';
    rec3_array[1].enum_field = Paulaner;

    rec3_array[2].integer_field = 14;
    rec3_array[2].long_field = 987234;
    rec3_array[2].uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3_array[2].ulong_field = 0xf7e589ceec9dd130;
#else
    rec3_array[2].ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3_array[2].string = "testing";
    rec3_array[2].double_field = 2.717;
    rec3_array[2].string2 = NULL;
    rec3_array[2].char_field = 'A';
    rec3_array[2].enum_field = Paulaner;

    rec3_array[3].integer_field = 14;
    rec3_array[3].long_field = 987234;
    rec3_array[3].uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3_array[3].ulong_field = 0xf7e589ceec9dd130;
#else
    rec3_array[3].ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3_array[3].string = NULL;
    rec3_array[3].double_field = 2.717;
    rec3_array[3].string2 = NULL;
    rec3_array[3].char_field = 'A';
    rec3_array[3].enum_field = Paulaner;

    rec3_array[4].integer_field = 14;
    rec3_array[4].long_field = 987234;
    rec3_array[4].uint_field = 0xf7e589ce;	/* = 4159015374 */
#if SIZEOF_LONG==64
    rec3_array[4].ulong_field = 0xf7e589ceec9dd130;
#else
    rec3_array[4].ulong_field = 0xec9dd130;	/* = 3969765680 */
#endif
    rec3_array[4].string = "testing";
    rec3_array[4].double_field = 2.717;
    rec3_array[4].string2 = "jambalaya";
    rec3_array[4].char_field = 'A';
    rec3_array[4].enum_field = Paulaner;

    memset((char *) &rec4, 0, sizeof(rec4));
    for (i = 0; i < ARRAY_SIZE; i++) {
	rec4.int_array[i] = 297 + i;
    }
    rec4.double_array[0][0] = 1.0;
    rec4.double_array[0][1] = 2.0;
    rec4.double_array[1][0] = 3.0;
    rec4.double_array[1][1] = 4.0;
    rec4.ifield = -rec4.int_array[ARRAY_SIZE - 1];

    memset((char *) &rec5, 0, sizeof(rec5));
    rec5.earray[0].dfield = 4.0;
    rec5.earray[0].ifield = 4;
    rec5.earray[0].string = "string20";
    rec5.earray[1].dfield = 3.0;
    rec5.earray[1].ifield = 3;
    rec5.earray[1].string = "string15";
    rec5.earray[2].dfield = 2.0;
    rec5.earray[2].ifield = 2;
    rec5.earray[2].string = "string10";
    rec5.earray[3].dfield = 1.0;
    rec5.earray[3].ifield = 1;
    rec5.earray[3].string = "string5";

    k = 0;
    for (i = 1; i < 20; i += 5) {
	rec6_array[k].string = malloc(10);
	memset(rec6_array[k].string, 0, 10);
	sprintf(rec6_array[k].string, "variant%d", i);
	rec6_array[k].icount = 2 * i;
	rec6_array[k].var_int_array = malloc(sizeof(((sixth_rec_ptr) 0)->var_int_array[0]) * rec6_array[k].icount);
	rec6_array[k].var_double_array = malloc(sizeof(double) * rec6_array[k].icount);
	rec6_array[k].var_string_array = malloc(sizeof(second_rec) * rec6_array[k].icount);
	for (j = 0; j < rec6_array[k].icount; j++) {
	    rec6_array[k].var_int_array[j] = 297 + j;
	    rec6_array[k].var_double_array[j] = 2.717 * j;
	    rec6_array[k].var_string_array[j].integer_field = -14;
	    rec6_array[k].var_string_array[j].short_field = j;
	    rec6_array[k].var_string_array[j].long_field = -9876129;
	    rec6_array[k].var_string_array[j].string = malloc(15);
	    memset(rec6_array[k].var_string_array[j].string, 0, 15);
	    sprintf(rec6_array[k].var_string_array[j].string,
		    "substring%d", j);
	    rec6_array[k].var_string_array[j].double_field = 2.717;
	    rec6_array[k].var_string_array[j].char_field = 'a' + 2 * j;
	}
	k++;
    }
    rec7_array[0].integer_field = 47;
    rec7_array[0].nested_rec.integer_field = -14;
    rec7_array[0].nested_rec.short_field = 27;
    rec7_array[0].nested_rec.long_field = -9876129;
    rec7_array[0].nested_rec.string = "Another string";
    rec7_array[0].nested_rec.double_field = 2.717;
    rec7_array[0].nested_rec.char_field = 'A';
    rec7_array[0].string = "Yet another string";

    memset((char *) &rec8_array[0], 0, sizeof(rec8_array));
    rec8_array[0].integer_field = 9872346;
    rec8_array[0].string = "ABCD";
    rec8_array[0].double_field = 3.14159265358797323;

    rec8_array[1].integer_field = 23462346;
    rec8_array[1].string = "Efghij";
    rec8_array[1].double_field = 3.14159265358797323 * 2.0;

    rec8_array[2].integer_field = 2346987;
    rec8_array[2].string = "Klmn";
    rec8_array[2].double_field = 3.14159265358797323 * 3.0;

    index = 0;
    for (i = 1; i < 10; i += 2) {
        memset((char *) &rec9_array[index], 0, sizeof(ninth_rec));
	memset((char *) &string_array_array[index], 0, 
	       sizeof(string_array_rec));
	rec9_array[index].vec_length = i;
	rec9_array[index].eventv = malloc(sizeof(((ninth_rec_ptr) 0)->eventv[0]) * rec9_array[index].vec_length);

	string_array_array[index].array_len = i;
	string_array_array[index].array = malloc(sizeof(char*) * i);

	for (j = 0; j < rec9_array[index].vec_length; j++) {
	    int k;
	    rec9_array[index].eventv[j].iov_len = j + i;
	    rec9_array[index].eventv[j].iov_base = malloc(j + i);
	    string_array_array[index].array[j] = malloc(i + j + 2);

	    for (k=0; k<j+i; k++) {
		((char*)rec9_array[index].eventv[j].iov_base)[k] = 'A' + k + i/5;
		string_array_array[index].array[j][k] = 'a' + k + i;
	    }
	    string_array_array[index].array[j][k] = 0;
	}
	if ((i %4) == 1) {
	    string_array_array[index].base_string = NULL;
	} else {
	    string_array_array[index].base_string = strdup("Whoa there!");
	}
	index++;
    }
    derive.chan_str = "b013050800000000";
    derive.cond = 2; 
    derive.client_channel_id.len = 12; 
    derive.client_channel_id.channel = malloc(12);
    for(i=0; i < 12; i++) derive.client_channel_id.channel[i] = 2 * 1 + 1;
    derive.client_contact_str = "AAIAAEFQSUPEBc+CUFBJQ3yUAAA=";
    derive.filter = "{ return 1;}"; 
    derive.field_list_len = 3; 
    derive.field_list = malloc(3*sizeof(derive.field_list[0]));
    derive.field_list[0].field_name = "num_points"; 
    derive.field_list[0].field_type = "integer"; 
    derive.field_list[0].field_size = 4;  
    derive.field_list[0].field_offset = 0; 
    derive.field_list[1].field_name = "image_data";  
    derive.field_list[1].field_type = "PolygonPoints[num_points]"; 
    derive.field_list[1].field_size = 8;  
    derive.field_list[1].field_offset = 4; 
    derive.field_list[2].field_name = ""; 
    derive.field_list[2].field_type = ""; 
    derive.field_list[2].field_size = 0; 
    derive.field_list[2].field_offset = 0;
    derive.format_list_len = 2; 
    derive.format_list = malloc(2*sizeof(derive.format_list[0]));
    derive.format_list[0].format_name = "PipelinedPoint"; 
    derive.format_list[0].field_list_len = 7; 
    derive.format_list[0].field_list = 
      malloc(7*sizeof(derive.format_list[0].field_list[0]));
    derive.format_list[0].field_list[0].field_name = "x"; 
    derive.format_list[0].field_list[0].field_type = "integer"; 
    derive.format_list[0].field_list[0].field_size = 2; 
    derive.format_list[0].field_list[0].field_offset = 0;
    derive.format_list[0].field_list[1].field_name = "y"; 
    derive.format_list[0].field_list[1].field_type = "integer"; 
    derive.format_list[0].field_list[1].field_size = 2; 
    derive.format_list[0].field_list[1].field_offset = 2;
    derive.format_list[0].field_list[2].field_name = "z"; 
    derive.format_list[0].field_list[2].field_type = "integer"; 
    derive.format_list[0].field_list[2].field_size = 2; 
    derive.format_list[0].field_list[2].field_offset = 4; 
    derive.format_list[0].field_list[3].field_name = "r"; 
    derive.format_list[0].field_list[3].field_type = "integer"; 
    derive.format_list[0].field_list[3].field_size = 2; 
    derive.format_list[0].field_list[3].field_offset = 6;
    derive.format_list[0].field_list[4].field_name = "g"; 
    derive.format_list[0].field_list[4].field_type = "integer"; 
    derive.format_list[0].field_list[4].field_size = 2; 
    derive.format_list[0].field_list[4].field_offset = 8;
    derive.format_list[0].field_list[5].field_name = "b"; 
    derive.format_list[0].field_list[5].field_type = "integer"; 
    derive.format_list[0].field_list[5].field_size = 2; 
    derive.format_list[0].field_list[5].field_offset = 10;
    derive.format_list[0].field_list[6].field_name = ""; 
    derive.format_list[0].field_list[6].field_type = ""; 
    derive.format_list[0].field_list[6].field_size = 0; 
    derive.format_list[0].field_list[6].field_offset = 0; 
    derive.format_list[1].format_name = "PolygonPoints"; 
    derive.format_list[1].field_list_len = 3; 
    derive.format_list[1].field_list = 
      malloc(3*sizeof(derive.format_list[1].field_list[0]));
    derive.format_list[1].field_list[0].field_name = "num_points"; 
    derive.format_list[1].field_list[0].field_type = "integer"; 
    derive.format_list[1].field_list[0].field_size = 4; 
    derive.format_list[1].field_list[0].field_offset = 0;
    derive.format_list[1].field_list[1].field_name = "polygon_points"; 
    derive.format_list[1].field_list[1].field_type = "PipelinedPoint[num_points]"; 
    derive.format_list[1].field_list[1].field_size = 12; 
    derive.format_list[1].field_list[1].field_offset = 4;
    derive.format_list[1].field_list[2].field_name = ""; 
    derive.format_list[1].field_list[2].field_type = ""; 
    derive.format_list[1].field_list[2].field_size = 0; 
    derive.format_list[1].field_list[2].field_offset = 0;
    derive.init_data_block = NULL; 
    derive.init_data_len = 0;

    multi_array.ifield = 4;
    for (i = 0; i < 2; i++) {
	for (j = 0; j < 2; j++) {
	    for (k = 0; k < 2; k++) {
		for (l = 0; l < 2; l++) {
		    multi_array.double_array[i][j][k][l] = 
			1000*i + 100*j + 10*k +l;
		}
	    }
	}
    }
    multi_array.int_array = malloc(2*4*sizeof(int));
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 2; j++) {
	    multi_array.int_array[i][j] = 1000*i + 100*j;
	}
    }
    multi_array.int_array2 = malloc(4*2*sizeof(int));
    for (i = 0; i < 2; i++) {
	for (j = 0; j < 4; j++) {
	    multi_array.int_array2[i][j] = 1000*i + 100*j;
	}
    }
    multi_array.int_array3 = malloc(4*4*4*sizeof(int));
    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    for (k = 0; k < 4; k++) {
		multi_array.int_array3[i][j][k] = 1000*i + 100*j + 10*k;
	    }
	}
    }
}

