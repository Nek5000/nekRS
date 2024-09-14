#include <stdio.h>
#include <stdlib.h>
#include "string.h"

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;

static char *filter_func = "{\n\
int hop_count;\n\
hop_count = attr_ivalue(event_attrs, \"hop_count_atom\");\n\
hop_count++;\n\
set_int_attr(event_attrs, \"hop_count_atom\", hop_count);\n\
return 1;\n\
}\0\0";

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    int hop_count;
    atom_t hop_count_atom = -1;
    if (hop_count_atom == -1) {
	hop_count_atom = attr_atom_from_string("hop_count_atom");
    }
    get_int_attr(attrs, hop_count_atom, &hop_count);
    checksum_simple_record(event, attrs, quiet);
    EVclient_shutdown(test_client, 0);
    if (!quiet) printf("\nreceived had %d hops\n", hop_count);
    return 0;
}


extern int
be_test_master(int argc, char **argv)
{
    char **nodes;
    CManager cm;
    attr_list contact_list;
    char *str_contact;
    EVdfg_stone src, last, tmp, sink;
    EVsource source_handle;
    int node_count = 3;
    int i;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    if (argc == 1) {
	sscanf(argv[0], "%d", &node_count);
    }
    nodes = malloc(sizeof(nodes[0]) * (node_count+1));
    for (i=0; i < node_count; i++) {
	nodes[i] = malloc(25);
	sprintf(nodes[i], "N%d", i);
    }
    nodes[i] = NULL;
    cm = CManager_create();
    CMlisten(cm);
    contact_list = CMget_contact_list(cm);
    str_contact = attr_list_to_string(contact_list);
    free_attr_list(contact_list);

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
    EVmaster_register_node_list(test_master, &nodes[0]);
    test_dfg = EVdfg_create(test_master);
    src = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(src, nodes[0]);
    char *filter;
    filter = create_filter_action_spec(NULL, filter_func);

    last = src;

    for (i=1; i < node_count -1; i++) {
        tmp = EVdfg_create_stone(test_dfg, filter);
	EVdfg_link_dest(last, tmp);
	EVdfg_assign_node(tmp, nodes[i]);
	last = tmp;
    }
    free(filter);
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_link_dest(last, sink);
    EVdfg_assign_node(sink, nodes[node_count-1]);

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
	simple_rec rec;
	atom_t hop_count_atom;
	attr_list attrs = create_attr_list();
	hop_count_atom = attr_atom_from_string("hop_count_atom");
	add_int_attr(attrs, hop_count_atom, 1);
	generate_simple_record(&rec);
	/* submit would be quietly ignored if source is not active */
	EVsubmit(source_handle, &rec, attrs);
	free_attr_list(attrs);
    }

    EVfree_source(source_handle);
    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);
    free(str_contact);
    for (i=0; i < node_count; i++) {
	free(nodes[i]);
    }
    free(nodes);
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
    return EVclient_wait_for_shutdown(test_client);
}
