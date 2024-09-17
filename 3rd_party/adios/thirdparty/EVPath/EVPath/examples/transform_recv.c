#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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

/* this file is evpath/examples/transform_recv.c */
typedef struct _output_rec {
    int integer_field;
    double average;
    char *str;
} output_rec, *output_rec_ptr;

static FMField output_field_list[] =
{
    {"integer_field", "integer", sizeof(int), FMOffset(output_rec_ptr, integer_field)},
    {"average", "double", sizeof(double), FMOffset(output_rec_ptr, average)},
    {"str", "string", sizeof(char*), FMOffset(output_rec_ptr, str)},
    {NULL, NULL, 0, 0}
};
static FMStructDescRec output_format_list[] =
{
    {"simple2", output_field_list, sizeof(output_rec), NULL},
    {NULL, NULL}
};

static int
output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    output_rec_ptr event = vevent;
    printf("I got %d, average is now %g\n", event->integer_field, event->average);
    printf("Base event is %p event, string is %p\n", event, event->str);
    printf("real string is %s\n", event->str);
    return 1;
}

int main(int argc, char **argv)
{
    CManager cm;
    EVstone stone;
    char *string_list, *trans_spec, *encoded_trans_spec;
    char *trans_func = "\
{\
    static double sum = 0.0;\
    static int count = 0;\
    sum = sum + input.integer_field;\
    count++;\
    output.integer_field = input.integer_field;\
    output.average = sum / count;\
    output.str = input.str;\
    return (count % 5) == 0;  /* pass filter every fifth*/ \
}";

    cm = CManager_create();
    CMlisten(cm);

    stone = EValloc_stone(cm);
    EVassoc_terminal_action(cm, stone, output_format_list, output_handler, NULL);
    string_list = attr_list_to_string(CMget_contact_list(cm));
    trans_spec = create_transform_action_spec(simple_format_list, output_format_list,
					    trans_func);
    printf("trns spec is %s\n", trans_spec);
    encoded_trans_spec = atl_base64_encode(trans_spec, (unsigned int) strlen(trans_spec) + 1);
    printf("Contact list \"%d:%s:%s\"\n", stone, string_list, encoded_trans_spec);
    free(trans_spec);
    free(encoded_trans_spec);
    CMsleep(cm, 600);
    return 0;
}
