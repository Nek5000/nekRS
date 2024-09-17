/*
 *  Test the ability of a stone to have sink handlers for multiple two different types of incoming data.
 */

#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

typedef struct _rec_a {
    int a_field;
} rec_a, *rec_a_ptr;

static FMField a_field_list[] =
{
    {"a_field", "integer",
     sizeof(int), FMOffset(rec_a_ptr, a_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec a_format_list[] =
{
    {"a_rec", a_field_list, sizeof(rec_a), NULL},
    {NULL, NULL, 0, NULL}
};

static int status;
static EVclient test_client;

static int message_count = 0;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    checksum_simple_record(event, attrs, quiet);
    message_count++;
    if (message_count == 2) {
	EVclient_shutdown(test_client, 0);
    }
    return 0;
}

static
int
simple_handler_2(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    rec_a_ptr event = vevent;
    (void)cm;
    if (event->a_field % 2 != 1) {
	printf("Received record should be odd, got %d\n", event->a_field);
    }
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	a_field = %d\n", event->a_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    message_count++;
    if (message_count == 2) {
	EVclient_shutdown(test_client, 0);
    }
    return 0;
}

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone src, src2, sink;
    EVsource source_handle, source_handle2;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_handle2 = EVcreate_submit_handle(cm, DFG_SOURCE, a_format_list);
    (void)EVclient_register_source("master_source", source_handle);
    source_capabilities = EVclient_register_source("a_source", source_handle2);
    (void)EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler_2", a_format_list,
				(EVSimpleHandlerFunc) simple_handler_2, NULL);

/*
**  MASTER and DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);

    test_dfg = EVdfg_create(test_master);
    src = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(src, "b");
    src2 = EVdfg_create_source_stone(test_dfg, "a_source");
    EVdfg_assign_node(src2, "b");
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_add_sink_action(sink, "simple_handler_2");
    EVdfg_assign_node(sink, "a");
    EVdfg_link_dest(src, sink);
    EVdfg_link_dest(src2, sink);

    EVdfg_realize(test_dfg);

/* We're node 0 in the process group */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

/* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    free(str_contact);
    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    
    if (EVclient_source_active(source_handle)) {
	simple_rec rec;
	generate_simple_record(&rec);
	/* submit would be quietly ignored if source is not active */
	printf("Submitting simple\n");
	EVsubmit(source_handle, &rec, NULL);

	rec_a a;
	a.a_field = 3;
	/* submit would be quietly ignored if source is not active */
	printf("Submitting a\n");
	EVsubmit(source_handle2, &a, NULL);
    }

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);

    EVfree_source(source_handle);
    EVfree_source(source_handle2);
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource src, src2;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  multi_sink  <nodename> <mastercontact>\n");
	exit(1);
    }

    src = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    src2 = EVcreate_submit_handle(cm, DFG_SOURCE, a_format_list);
    (void)EVclient_register_source("master_source", src);
    source_capabilities = EVclient_register_source("a_source", src2);
    (void)EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler_2", a_format_list,
				(EVSimpleHandlerFunc) simple_handler_2, NULL);
    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    if (EVclient_source_active(src)) {
	simple_rec rec;
	generate_simple_record(&rec);
	/* submit would be quietly ignored if source is not active */
	EVsubmit(src, &rec, NULL);

	rec_a a;
	a.a_field = 3;
	/* submit would be quietly ignored if source is not active */
	EVsubmit(src2, &a, NULL);
    }

    EVfree_source(src);
    EVfree_source(src2);
    return EVclient_wait_for_shutdown(test_client);
}
