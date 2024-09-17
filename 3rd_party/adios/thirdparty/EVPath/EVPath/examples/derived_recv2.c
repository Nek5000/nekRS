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

static int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    printf("I got %d\n", event->integer_field);
}

/* this file is evpath/examples/derived_recv2.c */
int main(int argc, char **argv)
{
    CManager cm;
    EVstone stone;
    char *string_list, *filter_spec, *encoded_filter_spec;

    cm = CManager_create();
    CMlisten(cm);

    stone = EValloc_stone(cm);
    EVassoc_terminal_action(cm, stone, simple_format_list, simple_handler, NULL);
    string_list = attr_list_to_string(CMget_contact_list(cm));
    filter_spec = create_filter_action_spec(simple_format_list, 
					    "{ return input.integer_field % 2;}");
    encoded_filter_spec = atl_base64_encode(filter_spec, strlen(filter_spec) + 1);
    printf("Contact list \"%d:%s:%s\"\n", stone, string_list, encoded_filter_spec);
    free(filter_spec);
    free(encoded_filter_spec);
    CMsleep(cm, 600);
}
