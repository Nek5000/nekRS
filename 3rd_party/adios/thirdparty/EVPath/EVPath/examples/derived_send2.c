#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "evpath.h"

typedef struct _simple_rec {
    int integer_field;
} simple_rec, *simple_rec_ptr;

static FMField simple_field_list[] =
{
    {"integer_field", "integer", sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {NULL, NULL, 0, 0}
};
static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {NULL, NULL}
};

typedef struct _second_rec {
    double data_field;
    char data_type;
} second_rec, *second_rec_ptr;

static FMField second_field_list[] =
{
    {"data_field", "float", sizeof(double), FMOffset(second_rec_ptr, data_field)},
    {"data_type", "char", sizeof(char), FMOffset(second_rec_ptr, data_type)},
    {NULL, NULL, 0, 0}
};
static FMStructDescRec second_format_list[] =
{
    {"second", second_field_list, sizeof(second_rec), NULL},
    {NULL, NULL}
};

/* this file is evpath/examples/derived_send2.c */
int main(int argc, char **argv)
{
    CManager cm;
    simple_rec data;
    second_rec data2;
    EVstone split_stone;
    EVaction split_action;
    EVsource source, source2;
    int i;

    cm = CManager_create();
    CMlisten(cm);

    split_stone = EValloc_stone(cm);
    split_action = EVassoc_split_action(cm, split_stone, NULL);

/* this file is evpath/examples/derived_send2.c */
    for (i = 1; i < argc; i++) {
	char string_list[20480];
	attr_list contact_list;
	char **filter_specs = NULL, *filter_spec, *next;
	EVstone remote_stone, output_stone;
        if (sscanf(argv[i], "%d:%s", &remote_stone, &string_list[0]) != 2) {
	    printf("Bad argument \"%s\"\n", argv[i]);
	    exit(0);
	}
	filter_spec = strchr(string_list, ':');
	if (filter_spec != NULL) {	/* if there is a filter spec */
	    int filter_count = 0;
	    *filter_spec = 0;           /* terminate the contact list */
	    filter_spec++;		/* advance pointer to string start */
	    filter_specs = malloc(sizeof(filter_specs[0]) * 2);
	    while (filter_spec != NULL) {
		next = strchr(filter_spec, ':');
		if (next != NULL) {
		    *next = 0;
		    next++;
		}
		atl_base64_decode((unsigned char *)filter_spec, NULL);  /* decode in place */
		filter_specs = realloc(filter_specs, sizeof(filter_specs[0]) * (filter_count + 2));
		filter_specs[filter_count++] = filter_spec;
		filter_spec = next;
	    }
	    filter_specs[filter_count] = NULL;
	}

	/* regardless of filtering or not, we'll need an output stone */
	output_stone = EValloc_stone(cm);
	contact_list = attr_list_from_string(string_list);
	EVassoc_bridge_action(cm, output_stone, contact_list, remote_stone);

	if (filter_specs == NULL) {
	    EVaction_add_split_target(cm, split_stone, split_action, output_stone);
	} else {
	    int i = 0;
	    EVstone filter_stone = EValloc_stone(cm);
	    while (filter_specs[i] != NULL) {
		EVaction filter_action = EVassoc_immediate_action(cm, filter_stone, filter_specs[i], NULL);
		EVaction_set_output(cm, filter_stone, filter_action, 0, output_stone);
		i++;
	    }
	    EVaction_add_split_target(cm, split_stone, split_action, filter_stone);
	}
    }

    source = EVcreate_submit_handle(cm, split_stone, simple_format_list);
    source2 = EVcreate_submit_handle(cm, split_stone, second_format_list);
    data.integer_field = 318;
    data2.data_field = 18.8;
    data2.data_type = 'A';
    for (i=0; i < 20; i++) {
	if ((i % 2) == 0) {
	    EVsubmit(source, &data, NULL);
	    data.integer_field++;
	} else {
	    EVsubmit(source2, &data2, NULL);
	    data2.data_field += 1.0;
	    data2.data_type++;
	}
    }
    return 1;
}
