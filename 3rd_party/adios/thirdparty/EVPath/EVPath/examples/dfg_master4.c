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

void
join_handler_func(EVmaster master, char *identifier, void*cur_unused1, void*cur_unused2)
{
    EVdfg_stone src, mid, sink;
    static int i = 0;
    EVdfg dfg;
    char *expected_nodes[] = {"a", "b", "c"};

    /* Nth client joining is named for that element in expected_nodes[] */
    EVmaster_assign_canonical_name(master, identifier, expected_nodes[i]);
    i++;

    if (i < sizeof(expected_nodes)/sizeof(expected_nodes[0])) {
	/* not everyone is here yet, so return */
	return;
    }

    /* actually create the DFG and assign stones to nodes */
    dfg = EVdfg_create(master);

    /* EVdfg stone creation and assignment occurs here */
    src = EVdfg_create_source_stone(dfg, "event source");
    EVdfg_assign_node(src, "a");

    mid = EVdfg_create_stone(dfg, NULL);
    EVdfg_assign_node(mid, "b");
    EVdfg_link_port(src, 0, mid);

    sink = EVdfg_create_sink_stone(dfg, "simple_handler");
    EVdfg_assign_node(sink, "c");
    EVdfg_link_port(mid, 0, sink);

    EVdfg_realize(dfg);
}

/* this file is evpath/examples/dfg_master4.c */
int main(int argc, char **argv)
{
    CManager cm;
    char *str_contact;
    EVsource source_handle;
    EVmaster test_master;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_node_join_handler (test_master, join_handler_func);

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
