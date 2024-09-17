#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

/* this file is evpath/examples/net_send.c */
int main(int argc, char **argv)
{
    CManager cm;
    simple_rec data;
    EVstone stone;
    EVsource source;
    char string_list[2048];
    attr_list contact_list;
    EVstone remote_stone;

    if (sscanf(argv[1], "%d:%s", &remote_stone, &string_list[0]) != 2) {
	printf("Bad arguments \"%s\"\n", argv[1]);
	exit(0);
    }

    cm = CManager_create();
    CMlisten(cm);

    stone = EValloc_stone(cm);
    contact_list = attr_list_from_string(string_list);
    EVassoc_bridge_action(cm, stone, contact_list, remote_stone);

    source = EVcreate_submit_handle(cm, stone, simple_format_list);
    data.integer_field = 318;
    EVsubmit(source, &data, NULL);
    return 0;
}
