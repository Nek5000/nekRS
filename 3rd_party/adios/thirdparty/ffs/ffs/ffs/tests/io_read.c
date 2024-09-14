
#include "config.h"
#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <string.h>
#include "ffs.h"
#include "unix_defs.h"

typedef struct _first_rec {
    int integer_field;
    double double_field;
    char char_field;
} first_rec, *first_rec_ptr;

typedef enum {
    Red_Stripe, Paulaner, Pilsner
} enum_type;

static IOField field_list[] =
{
    {"integer field", "integer",
     sizeof(int), IOOffset(first_rec_ptr, integer_field)},
    {"double field", "float",
     sizeof(double), IOOffset(first_rec_ptr, double_field)},
    {"char field", "char",
     sizeof(char), IOOffset(first_rec_ptr, char_field)},
    {NULL, NULL, 0, 0}
};

#define ARRAY_SIZE 14

typedef struct _fourth_rec {
    long ifield;
    int int_array[ARRAY_SIZE];
#if SIZEOF_LONG_DOUBLE != 0 && SIZEOF_LONG_DOUBLE != SIZEOF_DOUBLE
    long double double_array[2][2];
#else
    double double_array[2][2];
#endif
} fourth_rec, *fourth_rec_ptr;

static IOField field_list4[] =
{
    {"ifield", "integer",
     sizeof(long), IOOffset(fourth_rec_ptr, ifield)},
    {"int_array", IOArrayDecl(integer, ARRAY_SIZE),
     sizeof(int), IOOffset(fourth_rec_ptr, int_array[0])},
    {"double field", "float[2][2]",
     sizeof(((fourth_rec_ptr) 0)->double_array[0][0]),
     IOOffset(fourth_rec_ptr, double_array[0][0])},
    {NULL, NULL, 0, 0}
};

int
main(argc, argv)
int argc;
char **argv;
{
    IOFile iofile = NULL;
    IOFormat two_string_format, first_ioformat, fourth_ioformat;
    IOFieldPtr long_field_descr, ulong_field_descr, uint_field_descr,
     string2_descr, integer_field_descr, double_field_descr, string_field_descr,
     char_field_descr, enum_field_descr;
    char *buffer;
    int buffer_size;
    int check_only = 0;

    switch (argc) {
    case 1:
	iofile = open_IOfile("test_output", "r");
	break;
    case 2:
	if (strcmp(argv[1], "-check") == 0) {
	    iofile = open_IOfile("test_output", "r");
	    check_only++;
	} else {
	    iofile = open_IOfile(argv[1], "r");
	}
	break;
    case 3:
	if (strcmp(argv[1], "-check") == 0) {
	    check_only++;
	} else {
	    printf("Unknown arg \"%s\"\n", argv[1]);
	    exit(1);
	}
	iofile = open_IOfile(argv[2], "r");
	break;
    }

    if (iofile == NULL) {
	printf("Open failed\n");
	exit(1);
    }
    if (!check_only)
	dump_IOFile(iofile);

    /* open two formats (Ignore the third) */
    two_string_format = get_IOformat_by_name(iofile, "two string format");
    assert(two_string_format != NULL);
    first_ioformat = get_IOformat_by_name(iofile, "first format");
    assert(first_ioformat != NULL);

    set_IOconversion(iofile, "first format", field_list, sizeof(first_rec));

    fourth_ioformat = get_IOformat_by_name(iofile, "internal array format");
    assert(fourth_ioformat != NULL);

    set_IOconversion(iofile, "internal array format", field_list4,
		     sizeof(fourth_rec));

    /* 
     * for the second, we'll get each field individually. 
     * (Necessary if the formats differ)  Get IOFieldPtr values for
     * each of the fields. 
     */
    long_field_descr = get_IOfieldPtr(iofile, "two string format",
				      "long field");
    assert(long_field_descr != NULL);

    uint_field_descr = get_IOfieldPtr(iofile, "two string format",
				      "uint field");
    assert(uint_field_descr != NULL);

    ulong_field_descr = get_IOfieldPtr(iofile, "two string format",
				       "ulong field");
    assert(ulong_field_descr != NULL);

    string2_descr = get_IOfieldPtr(iofile, "two string format",
				   "string field2");
    assert(string2_descr != NULL);

    integer_field_descr = get_IOfieldPtr(iofile, "two string format",
					 "integer field");
    assert(integer_field_descr != NULL);

    double_field_descr = get_IOfieldPtr(iofile, "two string format",
					"double field");
    assert(double_field_descr != NULL);

    string_field_descr = get_IOfieldPtr(iofile, "two string format",
					"string field");
    assert(string_field_descr != NULL);

    char_field_descr = get_IOfieldPtr(iofile, "two string format",
				      "char field");
    assert(char_field_descr != NULL);

    enum_field_descr = get_IOfieldPtr(iofile, "two string format",
				      "enum field");
    assert(enum_field_descr != NULL);


    /* start reading the file. allocate buffer space if you know
     * the max record size before hand, use that here.  variant records
     * (with strings) have variant size though, so you'll have to do
     * some dynamic stuff or set a max string length somehow 
     */
    buffer = malloc(2048);
    buffer_size = 2048;

    while (next_IOrecord_format(iofile)) {
	if (buffer_size < next_raw_IOrecord_length(iofile)) {
	    buffer_size = next_raw_IOrecord_length(iofile);
	    buffer = realloc(buffer, buffer_size);
	}
	if (next_IOrecord_format(iofile) == first_ioformat) {
	    if (next_IOrecord_count(iofile) == 1) {
		/* read data directly to mem structure */
		first_rec read_data;
		if (!read_IOfile(iofile, &read_data))
		    fprintf(stderr, "read failed for format1\n");
		printf("first format rec had %d, %g, %c\n",
		       read_data.integer_field, read_data.double_field,
		       read_data.char_field);
	    } else {
		int i, count = 7;
		first_rec *read_data = (first_rec *) malloc(sizeof(first_rec) *
							    count);
		if (next_IOrecord_count(iofile) < count) {
		    count = next_IOrecord_count(iofile);
		}
		if (read_array_IOfile(iofile, read_data,
				    count, sizeof(first_rec)) != count) {
		    printf("array read failed\n");
		}
		printf("read %d elements of first_rec:\n", count);
		for (i = 0; i < count; i++) {
		    printf("\t rec %d had %d, %g, %c\n", i,
			   read_data[i].integer_field,
			   read_data[i].double_field,
			   read_data[i].char_field);
		}
	    }
	} else if (next_IOrecord_format(iofile) == two_string_format) {
	    char *string1, *string2;
	    int int_var;
	    long long_var;
	    unsigned int uint_var;
	    unsigned long ulong_var;
	    enum_type enum_var;
	    double double_var;
	    char char_var;
	    if (!read_raw_IOfile(iofile, buffer, buffer_size, NULL))
		fprintf(stderr, "read failed for format2\n");
	    string1 = get_IOstring(string_field_descr, buffer);
	    string2 = get_IOstring(string2_descr, buffer);
	    if (string1 == NULL)
		string1 = "(null)";
	    if (string2 == NULL)
		string2 = "(null)";
	    int_var = get_IOint(integer_field_descr, buffer);
	    double_var = get_IOdouble(double_field_descr, buffer);
	    long_var = get_IOlong(long_field_descr, buffer);
	    ulong_var = get_IOulong(ulong_field_descr, buffer);
	    uint_var = get_IOuint(uint_field_descr, buffer);
	    char_var = get_IOchar(char_field_descr, buffer);
	    enum_var = (enum_type) get_IOenum(enum_field_descr, buffer);
	    printf("int=%d, long=%ld, uint=%u, ulong=%lu, string1=%s, double=%g, string2=%s, char=%c, enum=",
	     int_var, long_var, uint_var, ulong_var, string1, double_var,
		   string2, char_var);
	    switch (enum_var) {
	    case Red_Stripe:
		printf("Red_Stripe\n");
		break;
	    case Paulaner:
		printf("Paulaner\n");
		break;
	    case Pilsner:
		printf("Pilsner\n");
		break;
	    default:
		printf("Unknown_ enum_%d\n", (int) enum_var);
		break;
	    }
	} else if (next_IOrecord_format(iofile) == fourth_ioformat) {
	    /* read data directly to mem structure */
	    int i;
	    fourth_rec rec4;
	    if (!read_IOfile(iofile, &rec4))
		fprintf(stderr, "read failed for format1\n");
	    printf("fourth format rec had %ld, in_array = (%d",
		   rec4.ifield, rec4.int_array[0]);
	    for (i = 1; i < ARRAY_SIZE; i++) {
		printf(", %d", rec4.int_array[i]);
	    }
	    printf(")\n\t\tdouble array is ((%g, %g), (%g, %g))\n",
		   (double) rec4.double_array[0][0],
		   (double) rec4.double_array[0][1],
		   (double) rec4.double_array[1][0],
		   (double) rec4.double_array[1][1]);
	} else {
	    /* read and discard */
	    printf("discarding a record\n");
	    read_IOfile(iofile, NULL);
	}
    }
    close_IOfile(iofile);
    free_IOfile(iofile);
    return 0;
}
