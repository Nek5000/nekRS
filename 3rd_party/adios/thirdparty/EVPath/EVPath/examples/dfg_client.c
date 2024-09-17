#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "evpath.h"
#include "ev_dfg.h"

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
    return 1;
}

/* this file is evpath/examples/dfg_client.c */
int main(int argc, char **argv)
{
    CManager cm;
    EVsource source_handle;
    EVclient test_client;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, -1, simple_format_list);
    source_capabilities = EVclient_register_source("event source", source_handle);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

    /* We're node argv[1] in the set of process clients, contact list is argv[2] */
    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);

    if (EVclient_ready_wait(test_client) != 1) {
	/* initialization failed! */
	exit(1);
    }

    
    if (EVclient_source_active(source_handle)) {
	simple_rec rec;
	rec.integer_field = 318;
	/* submit would be quietly ignored if source is not active */
	EVsubmit(source_handle, &rec, NULL);
    }

    CMrun_network(cm);
    return 0;
}
