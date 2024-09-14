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

typedef struct _output_rec {
    int integer_field;
    double average;
} output_rec, *output_rec_ptr;

static FMField output_field_list[] =
{
    {"integer_field", "integer", sizeof(int), FMOffset(output_rec_ptr, integer_field)},
    {"average", "double", sizeof(double), FMOffset(output_rec_ptr, average)},
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
    return 0;
}

int main(int argc, char **argv)
{
    CManager cm;
    EVstone stone;
    char *string_list, *trans_spec, *encoded_trans_spec, *trans_spec2, *encoded_trans_spec2;
/* this file is evpath/examples/transform_recv3.c */
    char *trans_func = "\
{\
    double sum = 0.0;\
    int count = 0;\
    if (attr_set(stone_attrs, \"sum\")) {\n\
	sum = attr_dvalue(stone_attrs, \"sum\");\n\
	count = attr_ivalue(stone_attrs, \"count\");\n\
    }\n\
    sum = sum + input.integer_field;\
    count++;\
    output.integer_field = input.integer_field;\
    output.average = sum / count;\
    set_double_attr(stone_attrs, \"sum\", sum);\n\
    set_int_attr(stone_attrs, \"count\", count);\n\
    return (count % 5) == 0;  /* pass filter every fifth*/ \
}";

    char *trans_func2 = "\
{\
    double sum = 0.0;\
    int count = 0;\
    if (attr_set(stone_attrs, \"sum\")) {\n\
	sum = attr_dvalue(stone_attrs, \"sum\");\n\
	count = attr_ivalue(stone_attrs, \"count\");\n\
    }\n\
    sum = sum + input.data_field;\
    count++;\
    output.integer_field = input.data_field;\
    output.average = sum / count;\
    set_double_attr(stone_attrs, \"sum\", sum);\n\
    set_int_attr(stone_attrs, \"count\", count);\n\
    return (count % 5) == 0;  /* pass filter every fifth*/ \
}";

    cm = CManager_create();
    CMlisten(cm);

    stone = EValloc_stone(cm);
    EVassoc_terminal_action(cm, stone, output_format_list, output_handler, NULL);
    string_list = attr_list_to_string(CMget_contact_list(cm));
    trans_spec = create_transform_action_spec(simple_format_list, output_format_list,
					    trans_func);
    trans_spec2 = create_transform_action_spec(second_format_list, output_format_list,
					    trans_func2);
    encoded_trans_spec = atl_base64_encode(trans_spec, (unsigned int) strlen(trans_spec) + 1);
    encoded_trans_spec2 = atl_base64_encode(trans_spec2, (unsigned int) strlen(trans_spec2) + 1);
    printf("Contact list \"%d:%s:%s:%s\"\n", stone, string_list, encoded_trans_spec,
	   encoded_trans_spec2);
    free(trans_spec);
    free(trans_spec2);
    free(encoded_trans_spec);
    free(encoded_trans_spec2);
    CMsleep(cm, 600);
    return 0;
}

