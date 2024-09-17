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
#include <stdio.h>
#include <string.h>
#include "ffs.h"

#include "test_funcs.h"

static int verbose = 0;
static char *test_only = NULL;

size_t
size_func_sizeof(FFSFile iofile, size_t size)
{
    return size;
}

size_t
size_func_next_size(FFSFile iofile, size_t size)
{
    return FFSfile_next_decode_length(iofile);
}

int
read_func_no_buffer(FFSFile iofile, void *data, size_t size)
{
    return FFSread(iofile, data);
}

int
read_func_buffer(FFSFile iofile, void *data, size_t size)
{
    FFSBuffer b = create_fixed_FFSBuffer(data, size);
    int ret = FFSread_to_buffer(iofile, b, NULL);
    free(b);
    return ret;
}


static FFSTypeHandle first_rec_ioformat, second_rec_ioformat, third_rec_ioformat;
static FFSTypeHandle fourth_rec_ioformat, later_rec_ioformat, nested_rec_ioformat;
static FFSTypeHandle fifth_rec_ioformat, sixth_rec_ioformat;
static FFSTypeHandle ninth_rec_ioformat, string_array_ioformat, derive_ioformat;
static FFSTypeHandle multi_array_ioformat, triangle_ioformat, add_action_ioformat;
static FFSTypeHandle node_ioformat, embedded_rec_ioformat;

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
    if ((test_only == NULL) || (strcmp(test_only, "derive") == 0)) {
	derive_ioformat = FFSset_fixed_target(context, derive_format_list);
	(void) derive_ioformat;
    }
    if ((test_only == NULL) || (strcmp(test_only, "multi_array") == 0)) {
	multi_array_ioformat = FFSset_fixed_target(context, multi_array_format_list);
	(void) multi_array_ioformat;
    }
    if ((test_only == NULL) || (strcmp(test_only, "triangle_param") == 0)) {
	triangle_ioformat = FFSset_fixed_target(context, triangle_format_list);
	(void) triangle_ioformat;
    }
    if ((test_only == NULL) || (strcmp(test_only, "add_action") == 0)) {
	add_action_ioformat = FFSset_fixed_target(context, add_action_format_list);
	(void) add_action_ioformat;
    }
    if ((test_only == NULL) || (strcmp(test_only, "node") == 0)) {
	node_ioformat = FFSset_fixed_target(context, node_format_list);
	(void) node_ioformat;
    }
    if ((test_only == NULL) || (strcmp(test_only, "embedded") == 0))
	embedded_rec_ioformat = FFSset_fixed_target(context, embedded_format_list);
}

void
do_test(char *input_file, size_t (*size_func) (FFSFile, size_t), int (*read_func) (FFSFile iofile, void *data, size_t size))
{
    FFSFile iofile;
    int finished = 0;
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
    iofile = open_FFSfile(input_file, "r");

    init_written_data();

    set_targets(FFSContext_of_file(iofile));

    while (!finished) {
	char *comment;
	FFSRecordType next = FFSnext_record_type(iofile);

	switch (next) {
	case FFScomment:
	    comment = FFSread_comment(iofile);
	    (void) comment;
	    break;
	case FFSformat:
	    /* ignore format, get to the next data record */
	    (void) FFSnext_type_handle(iofile);
	    break;
	case FFSdata:
	    if (verbose) 
		printf("Record of format %s\n",
		       FFSTypeHandle_name(FFSnext_type_handle(iofile)));
	    if (FFSnext_type_handle(iofile) == first_rec_ioformat) {
		first_rec read_data[10];
		memset(&read_data[0], 0, sizeof(first_rec));
		if (!FFSread(iofile, &read_data[0]))
		    printf("read first data");
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
	    } else if (FFSnext_type_handle(iofile) == second_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(second_rec));
		second_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read second data failed\n");
		if (!second_rec_eq(read_data, &rec2_array[second_rec_count++])) {
		    printf("Rec2 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == third_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(third_rec));
		third_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read third data failed\n");
		if (!third_rec_eq(read_data, &rec3_array[third_rec_count++])) {
		    printf("Rec3 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == fourth_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(fourth_rec));
		fourth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read fourth data failed\n");
		if (!fourth_rec_eq(read_data, &rec4)) {
		    printf("Rec4 failure\n");
		    exit(1);
		}
		free(read_data);
		fourth_rec_count++;
	    } else if (FFSnext_type_handle(iofile) == embedded_rec_ioformat) {

		printf("Emb Rec failure\n");
		FFSread(iofile, NULL);
	    } else if (FFSnext_type_handle(iofile) == fifth_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(fifth_rec));
		fifth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read fifth data failed\n");
		if (!fifth_rec_eq(read_data, &rec5)) {
		    printf("Rec5 failure\n");
		    exit(1);
		}
		free(read_data);
		fifth_rec_count++;
	    } else if (FFSnext_type_handle(iofile) == sixth_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(sixth_rec));
		sixth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read variant format failed\n");
		if (!sixth_rec_eq(read_data, &rec6_array[sixth_rec_count++])) {
		    printf("Rec6 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == nested_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(nested_rec));
		nested_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read variant format");
		if (!nested_rec_eq(read_data, &rec7_array[nested_rec_count++])) {
		    printf("Rec7 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == later_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(later_rec));
		later_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read variant format");
		if (!later_rec_eq(read_data, &rec8_array[later_rec_count++])) {
		    printf("Rec8 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == ninth_rec_ioformat) {
		size_t size = size_func(iofile, sizeof(ninth_rec));
		ninth_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read variant format");
		if (!ninth_rec_eq(read_data, &rec9_array[ninth_rec_count++])) {
		    printf("Rec9 failure\n");
		    exit(1);
		}
		free(read_data);
	    } else if (FFSnext_type_handle(iofile) == string_array_ioformat) {
		size_t size = size_func(iofile, sizeof(string_array_rec));
		string_array_rec *read_data = malloc(size);
		memset(read_data, 0, size);
		if (!read_func(iofile, read_data, size))
		    printf("read variant format");
		if (!string_array_eq(read_data, &string_array_array[string_array_rec_count++])) {
		    printf("string_array failure\n");
		    exit(1);
		}
		free(read_data);

	    } else {
		if (test_only == NULL) {
		    printf("discarding a record\n");
		}
		FFSread(iofile, NULL);
		unknown_rec_count++;
	    }
	    break;
	case FFSerror:
	case FFSend:
	default:
	    finished++;
	    break;
	}
    }

    close_FFSfile(iofile);
    free_FFSfile(iofile);
    free_written_data();
    if (first_rec_count != 3) {
	printf("Missed first, expected 3, got %d\n", first_rec_count);
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
	printf("Missed fifth, expected 1, got %d\n", fifth_rec_count);
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
main(int argc, char **argv)
{
    int i;
    char *input_file = "test_output";

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
	    /* argument */
	    if (argv[i][1] == 'v') {
		verbose++;
	    } else if (argv[i][1] == 't') {
		test_only = argv[++i];
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
