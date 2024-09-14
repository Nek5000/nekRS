#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "evpath.h"

typedef struct _simple_rec {
    int integer_field;
    char *str;
} simple_rec, *simple_rec_ptr;

static FMField simple_field_list[] =
{
    {"integer_field", "integer", sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {"str", "string", sizeof(char*), FMOffset(simple_rec_ptr, str)},
    {NULL, NULL, 0, 0}
};
static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {NULL, NULL}
};

/* this file is evpath/examples/derived_send.c */
int main(int argc, char **argv)
{
    CManager cm;
    simple_rec data;
    EVstone split_stone;
    EVaction split_action;
    EVsource source;
    int i;

    cm = CManager_create();
    CMlisten(cm);

    split_stone = EValloc_stone(cm);
    split_action = EVassoc_split_action(cm, split_stone, NULL);

/* this file is evpath/examples/derived_send.c */
    for (i = 1; i < argc; i++) {
	char string_list[2048];
	attr_list contact_list;
	char *filter_spec;
	EVstone remote_stone, output_stone;
        if (sscanf(argv[i], "%d:%s", &remote_stone, &string_list[0]) != 2) {
	    printf("Bad argument \"%s\"\n", argv[i]);
	    exit(0);
	}
	filter_spec = strchr(string_list, ':');
	if (filter_spec != NULL) {	/* if there is a filter spec */
	    *filter_spec = 0;           /* terminate the contact list */
	    filter_spec++;		/* advance pointer to string start */
	    atl_base64_decode((unsigned char *)filter_spec, NULL);  /* decode in place */
	    printf("String list is %s\n", string_list);
	}

	/* regardless of filtering or not, we'll need an output stone */
	output_stone = EValloc_stone(cm);
	contact_list = attr_list_from_string(string_list);
	printf("This is the contact list   --------");
	dump_attr_list(contact_list);
	EVassoc_bridge_action(cm, output_stone, contact_list, remote_stone);

	if (filter_spec == NULL) {
	    EVaction_add_split_target(cm, split_stone, split_action, output_stone);
	} else {
	    EVstone filter_stone = EValloc_stone(cm);
	    EVaction filter_action = EVassoc_immediate_action(cm, filter_stone, filter_spec, NULL);
	    EVaction_set_output(cm, filter_stone, filter_action, 0, output_stone);
	    EVaction_add_split_target(cm, split_stone, split_action, filter_stone);
	}
    }

    source = EVcreate_submit_handle(cm, split_stone, simple_format_list);
    data.integer_field = 318;
    data.str = "kraut";
    for (i=0; i < 10; i++) {
	EVsubmit(source, &data, NULL);
	data.integer_field++;
    }
    return 0;
}
