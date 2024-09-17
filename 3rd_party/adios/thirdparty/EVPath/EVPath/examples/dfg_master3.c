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

EVclient test_client;

static int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    printf("I got %d\n", event->integer_field);
    EVclient_shutdown(test_client, event->integer_field == 318);
    return 1;
}

/* this file is evpath/examples/dfg_master3.c */
int main(int argc, char **argv)
{
/*! [Changed nodes array] */
    char *nodes[] = {"a", "b", "c", NULL};
/*! [Changed nodes array] */
    CManager cm;
    EVmaster test_master;
    EVdfg test_dfg;
    char *str_contact;
    EVdfg_stone src, mid, sink;
    EVsource source_handle;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, -1, simple_format_list);
    source_capabilities = EVclient_register_source("event source", source_handle);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

/*
**  DFG CREATION
*/
    test_dfg = EVdfg_create(test_master);

/*! [Changed DFG Creation] */
    src = EVdfg_create_source_stone(test_dfg, "event source");
    EVdfg_assign_node(src, "a");

    mid = EVdfg_create_stone(test_dfg, NULL);
    EVdfg_assign_node(mid, "b");
    EVdfg_link_port(src, 0, mid);

    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_assign_node(sink, "c");
    EVdfg_link_port(mid, 0, sink);

    EVdfg_realize(test_dfg);
/*! [Changed DFG Creation] */

    /* We're node "a" in the DFG */
    test_client = EVclient_assoc_local(cm, "a", test_master, source_capabilities, sink_capabilities);

    printf("Contact list is \"%s\"\n", str_contact);
    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    
    if (EVclient_source_active(source_handle)) {
	simple_rec rec;
	rec.integer_field = 318;
	/* submit would be quietly ignored if source is not active */
	EVsubmit(source_handle, &rec, NULL);
    }

    if (EVclient_active_sink_count(test_client) > 0) {
	/* if there are active sinks, the handler will call EVclient_shutdown() */
    } else {
	if (EVclient_source_active(source_handle)) {
	    /* we had a source and have already submitted, indicate success */
	    EVclient_shutdown(test_client, 0 /* success */);
	} else {
	    /* we had neither a source or sink, ready to shutdown, no opinion */
	    EVclient_ready_for_shutdown(test_client);
	}
    }

    return(EVclient_wait_for_shutdown(test_client));
}
