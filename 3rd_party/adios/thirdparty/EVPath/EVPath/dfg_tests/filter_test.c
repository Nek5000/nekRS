#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;
static int repeat_count = 10;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    static int count = 0;
    simple_rec_ptr event = vevent;
    (void)cm;
    (void) client_data;
    checksum_simple_record(event, attrs, quiet);
    count++;
    if (count == repeat_count) 
	EVclient_shutdown(test_client, 0);
    return 0;
}


static FMField filter_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec filter_format_list[] =
{
    {"simple", filter_field_list, sizeof(simple_rec), NULL},
    {NULL, NULL}
};

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", "c", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone src, filter, sink;
    EVsource source_handle;
    char *filter_action_spec;
    EVdfg test_dfg;
    EVmaster test_master;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_capabilities = EVclient_register_source("master_source", source_handle);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

/*
**  MASTER AND DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);

    test_dfg = EVdfg_create(test_master);
    src = EVdfg_create_source_stone(test_dfg, "master_source");

    filter_action_spec = create_filter_action_spec(filter_format_list, "{int ret = input.long_field % 2;return ret;}\0\0");
    filter = EVdfg_create_stone(test_dfg, filter_action_spec);
    free(filter_action_spec);
    EVdfg_link_dest(src, filter);
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_link_dest(filter, sink);

    if ((argc != 1) || ((argc == 1) && (strcmp(argv[0], "3") == 0))) {
	EVmaster_register_node_list(test_master, &nodes[0]);
	EVdfg_assign_node(src, "a");
	EVdfg_assign_node(filter, "b");
	EVdfg_assign_node(sink, "c");
    } else if (strcmp(argv[0], "2a") == 0) {
	nodes[2] = NULL;
	EVmaster_register_node_list(test_master, &nodes[0]);
	EVdfg_assign_node(src, "a");
	EVdfg_assign_node(filter, "a");
	EVdfg_assign_node(sink, "b");
    } else if (strcmp(argv[0], "2b") == 0) {
	nodes[2] = NULL;
	EVmaster_register_node_list(test_master, &nodes[0]);
	EVdfg_assign_node(src, "a");
	EVdfg_assign_node(filter, "b");
	EVdfg_assign_node(sink, "b");
    } else if (strcmp(argv[0], "1") == 0) {
	nodes[1] = NULL;
	EVmaster_register_node_list(test_master, &nodes[0]);
	EVdfg_assign_node(src, "a");
	EVdfg_assign_node(filter, "a");
	EVdfg_assign_node(sink, "a");
    }
	

    EVdfg_realize(test_dfg);

/* We're node 0 in the DFG */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

/* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    
    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    if (EVclient_source_active(source_handle)) {
	int count = repeat_count;
	while (count != 0) {
	    simple_rec rec;
	    generate_simple_record(&rec);
	    EVsubmit(source_handle, &rec, NULL);
	    if ((rec.long_field%2 == 1) && (count != -1)) {
		count--;
	    }
	}
    }
    status = EVclient_wait_for_shutdown(test_client);
    free(str_contact);
    EVfree_source(source_handle);
    wait_for_children(nodes);

    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource src;
    int i;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  evtest  <nodename> <mastercontact>\n");
	exit(1);
    }

    src = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_capabilities = EVclient_register_source("master_source", src);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);
    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);
    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    if (EVclient_source_active(src)) {
	for (i=0; i < 20 ; i++) {
	    simple_rec rec;
	    generate_simple_record(&rec);
	    EVsubmit(src, &rec, NULL);
	}
    }
    EVfree_source(src);
    return EVclient_wait_for_shutdown(test_client);
}
