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
#include "ffs.h"
#include "unix_defs.h"

typedef struct _first_rec {
    int integer_field;
    double double_field;
    char char_field;
} first_rec, *first_rec_ptr;

typedef struct _direct_first_rec {
    short integer_field;
    double double_field;
    char char_field;
} direct_first_rec, *direct_first_rec_ptr;

typedef struct _first_rec2 {
    double double_field;
    int integer_field;
    char char_field;
} first_rec2, *first_rec2_ptr;

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

static IOField discard_field_list[] =
{
    {"integer field", "integer",
     sizeof(int), IOOffset(first_rec_ptr, integer_field)},
    {"char field", "char",
     sizeof(char), IOOffset(first_rec_ptr, char_field)},
    {NULL, NULL, 0, 0}
};

static IOField direct_list[] =
{
    {"integer field", "integer",
     sizeof(short), IOOffset(direct_first_rec_ptr, integer_field)},
    {"char field", "char",
     sizeof(char), IOOffset(direct_first_rec_ptr, char_field)},
    {NULL, NULL, 0, 0}
};

typedef struct _float_first_rec {
    int integer_field;
    float double_field;
    char char_field;
} float_first_rec, *float_first_rec_ptr;

static IOField float_field_list[] =
{
    {"integer field", "integer",
     sizeof(int), IOOffset(first_rec_ptr, integer_field)},
    {"double field", "float",
     sizeof(float), IOOffset(first_rec_ptr, double_field)},
    {NULL, NULL, 0, 0}
};

static IOField field2_list[] =
{
    {"double field", "float",
     sizeof(double), IOOffset(first_rec2_ptr, double_field)},
    {"integer field", "integer",
     sizeof(int), IOOffset(first_rec2_ptr, integer_field)},
    {NULL, NULL, 0, 0}
};

typedef enum {
    Red_Stripe, Paulaner, Pilsner
} enum_type;

typedef struct _third_rec {
    int integer_field;
    long long_field;
    unsigned int uint_field;
    unsigned long ulong_field;
    char *string;
    double double_field;
    char *string2;
    char char_field;
    enum_type enum_field;
} third_rec, *third_rec_ptr;

static IOField field_list3[] =
{
    {"integer field", "integer",
     sizeof(int), IOOffset(third_rec_ptr, integer_field)},
    {"long field", "integer",
     sizeof(long), IOOffset(third_rec_ptr, long_field)},
    {"uint field", "unsigned integer",
     sizeof(int), IOOffset(third_rec_ptr, uint_field)},
    {"ulong field", "unsigned integer",
     sizeof(long), IOOffset(third_rec_ptr, ulong_field)},
    {"string field", "string",
     sizeof(char *), IOOffset(third_rec_ptr, string)},
    {"double field", "float",
     sizeof(double), IOOffset(third_rec_ptr, double_field)},
    {"string field2", "string",
     sizeof(char *), IOOffset(third_rec_ptr, string2)},
    {"char field", "char",
     sizeof(char), IOOffset(third_rec_ptr, char_field)},
    {"enum field", "enumeration",
     sizeof(enum_type), IOOffset(third_rec_ptr, enum_field)},
    {NULL, NULL, 0, 0}
};

int
main(argc, argv)
int argc;
char **argv;
{
    IOFile iofile;
    IOFormat first_ioformat, third_ioformat;
    typedef enum {
	straight = 0, discard = 1, float_conv = 2, field2 = 3, direct = 4, last = 5
    } conv_types;
    conv_types conv_index = straight;

    if (argc != 2) {
	iofile = open_IOfile("test_output", "r");
    } else {
	iofile = open_IOfile(argv[1], "r");
    }

    dump_IOFile(iofile);

    first_ioformat = get_IOformat_by_name(iofile, "first format");
    assert(first_ioformat != NULL);

    third_ioformat = get_IOformat_by_name(iofile, "two string format");
    assert(third_ioformat != NULL);

    set_IOconversion(iofile, "two string format", field_list3,
		     sizeof(third_rec));

    /* start reading the file.  */
    while (next_IOrecord_format(iofile)) {
	if (next_IOrecord_format(iofile) == first_ioformat) {
	    switch (conv_index) {
	    case straight:
		{
		    first_rec read_data;
		    set_IOconversion(iofile, "first format",
				     field_list, sizeof(first_rec));

		    if (!read_IOfile(iofile, &read_data)) {
			fprintf(stderr, "read failed for format1\n");
		    }
		    printf("first format rec had %d, %g, %c\n",
			 read_data.integer_field, read_data.double_field,
			   read_data.char_field);
		    break;
		}
	    case discard:
		{
		    first_rec read_data;
		    set_IOconversion(iofile, "first format",
				  discard_field_list, sizeof(first_rec));

		    if (!read_IOfile(iofile, &read_data)) {
			fprintf(stderr, "read failed for format1\n");
		    }
		    printf("first format rec had %d, %g, %c\n",
			 read_data.integer_field, read_data.double_field,
			   read_data.char_field);
		    break;
		}
	    case float_conv:
		{
		    float_first_rec read_data;
		    set_IOconversion(iofile, "first format",
				     float_field_list,
				     sizeof(float_first_rec));

		    if (!read_IOfile(iofile, &read_data)) {
			fprintf(stderr, "read failed for float_conv\n");
		    }
		    printf("float rec had %d, %g\n",
			read_data.integer_field, read_data.double_field);
		    break;
		}
	    case field2:
		{
		    first_rec2 read_data;
		    set_IOconversion(iofile, "first format",
				     field2_list, sizeof(first_rec2));

		    if (!read_IOfile(iofile, &read_data)) {
			fprintf(stderr, "read failed for float2_conv\n");
		    }
		    printf("field2 format rec had %d, %g\n",
			   read_data.integer_field,
			   read_data.double_field);
		    break;
		}
	    case direct:
		{
		    direct_first_rec read_data;
		    set_IOconversion(iofile, "first format",
				     direct_list, sizeof(first_rec));

		    if (!read_IOfile(iofile, &read_data)) {
			fprintf(stderr, "read failed for direct\n");
		    }
		    printf("direct rec had %d, %g, %c\n",
			 read_data.integer_field, read_data.double_field,
			   read_data.char_field);
		    break;
		}
	    default:
		assert(0);
	    }
	    conv_index = (conv_types) ((((int) conv_index) + 1) % last);
	} else if (next_IOrecord_format(iofile) == third_ioformat) {
	    third_rec read_data;
	    if (!read_IOfile(iofile, &read_data)) {
		fprintf(stderr, "read failed for third rec\n");
	    }
	    if (read_data.string == NULL)
		read_data.string = "(null)";
	    if (read_data.string2 == NULL)
		read_data.string2 = "(null)";
	    printf("third rec had %d, %ld, %u, %lu, %s, %g, %s, %c, %d\n",
		   read_data.integer_field, read_data.long_field,
		   read_data.uint_field, read_data.ulong_field,
		   read_data.string, read_data.double_field,
		   read_data.string2, read_data.char_field,
		   (int) read_data.enum_field);
	} else {
	    /* read and discard */
	    read_raw_IOfile(iofile, NULL, 0, NULL);
	}
    }

    return 0;
}
