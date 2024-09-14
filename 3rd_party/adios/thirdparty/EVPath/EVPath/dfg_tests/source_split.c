/*
 *  Test the ability of a source stone to have multiple outputs .
 */

#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;

static int message_count[2] = {0, 0};

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    EVstone parent_stone = EVexecuting_stone(cm);
    attr_list sink_attrs = EVextract_attr_list(cm, parent_stone);
    atom_t DFG_SINK_ID;
    int sink_index = 0;
    if (quiet <= 0) {
	if (!sink_attrs) {
	    printf("Sink Stone Attributes is NULL\n");
	} else {
	    printf("Sink Stone Attributes are: ");
	    dump_attr_list(sink_attrs);
	}
    }
    DFG_SINK_ID = attr_atom_from_string("dfg_sink_id");
    if (get_int_attr(sink_attrs, DFG_SINK_ID, &sink_index) == 0) {
	printf("Attr not found\n");
	exit(1);
    }
    if (quiet <= 0) {
	printf("In the handler for stone %d\n", sink_index);
    }
    checksum_simple_record(event, attrs, quiet);
    message_count[sink_index]++;
    if ((message_count[0] == 1) && (message_count[1] == 1)) {
	EVclient_shutdown(test_client, 0);
    }
    return 0;
}

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone src, sink, sink2;
    EVsource source_handle;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;
    atom_t DFG_SINK_ID;
    DFG_SINK_ID = attr_atom_from_string("dfg_sink_id");
    attr_list attrs1, attrs2;

    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    (void)EVclient_register_source("master_source", source_handle);
    source_capabilities = EVclient_register_source("a_source", source_handle);
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", simple_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

/*
**  MASTER and DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);

    test_dfg = EVdfg_create(test_master);
    src = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(src, "a");
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    sink2 = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    attrs1 = create_attr_list();
    attrs2 = create_attr_list();
    add_int_attr(attrs1, DFG_SINK_ID, 0);
    add_int_attr(attrs2, DFG_SINK_ID, 1);
    EVdfg_set_attr_list(sink, attrs1);
    EVdfg_set_attr_list(sink2, attrs2);
    EVdfg_assign_node(sink, "a");
    EVdfg_assign_node(sink2, "a");
    EVdfg_link_dest(src, sink);
    EVdfg_link_dest(src, sink2);
    free_attr_list(attrs1);
    free_attr_list(attrs2);
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
	printf("Child usage:  multi_sink  <nodename> <mastercontact>\n");
	exit(1);
    }

    src = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_capabilities = EVclient_register_source("a_source", src);
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
    return EVclient_wait_for_shutdown(test_client);
}
