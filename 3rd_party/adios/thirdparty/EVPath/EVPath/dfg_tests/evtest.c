#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    checksum_simple_record(event, attrs, quiet);
    EVclient_shutdown(test_client, 0);
    return 0;
}


extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone src, sink;
    EVsource source_handle;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sources source_capabilities;
    EVclient_sinks sink_capabilities;

    (void)argc; (void)argv;
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
**  DFG CREATION
*/
    
    test_master = EVmaster_create(cm);
    EVmaster_register_node_list(test_master, &nodes[0]);

/* We're node 0 in the DFG */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

    str_contact = EVmaster_get_contact_list(test_master);


    test_dfg = EVdfg_create(test_master);
    src = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(src, "b");
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_assign_node(sink, "a");
    EVdfg_link_dest(src, sink);

    EVdfg_realize(test_dfg);


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
	EVsubmit(source_handle, &rec, NULL);
    }

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);

    EVfree_source(source_handle);
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource src;
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
	simple_rec rec;
	generate_simple_record(&rec);
	/* submit would be quietly ignored if source is not active */
	EVsubmit(src, &rec, NULL);
    }
    EVfree_source(src);
    status = EVclient_wait_for_shutdown(test_client);

    CManager_close(cm);
    return status;
}
