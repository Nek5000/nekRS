#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "config.h"
#include "ev_dfg.h"
#include "test_support.h"

#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#endif

static int status;
static EVclient test_client;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    static int handle_count = 0;
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    handle_count++;
    if (!quiet) {
	printf("Got %d handles, waiting on %d\n", handle_count, event->integer_field);
    }
    if (event->integer_field == handle_count) {
	EVclient_shutdown(test_client, 0);
    }
    return 0;
}


static char *router_function = "\
{\n\
    static int count = 0;\n\
    return (count++) % EVmax_output();\n\
}\0\0";

extern int
be_test_master(int argc, char **argv)
{
    char **nodes;
    CManager cm;
    char *str_contact;
    EVdfg_stone source, router;
    EVsource source_handle;
    int out_count, node_count;
    int ndig = 5;
    int i;
    int repeat_count = 40;
    char *router_action;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    srand48(time(NULL));
    out_count = lrand48() % 4 + 2;
    if (argc == 1) {
	sscanf(argv[0], "%d", &out_count);
    }
    if (!quiet) {
	printf("Running with out_count = %d\n", out_count);
    }
    repeat_count = ((int)(repeat_count / out_count) + 1) * out_count;

    node_count = out_count + 1;

    nodes = malloc(sizeof(nodes[0]) * (node_count+1));
    for (i=0; i < node_count; i++) {
	nodes[i] = malloc(ndig+10);
	sprintf(nodes[i], "N%d", i+1);
    }
    nodes[node_count] = NULL;
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
    EVmaster_register_node_list(test_master, &nodes[0]);
    test_dfg = EVdfg_create(test_master);

    router_action = create_router_action_spec(simple_format_list, router_function);

    source = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(source, nodes[0]);

    router = EVdfg_create_stone(test_dfg, router_action);
    free(router_action);
    EVdfg_assign_node(router, nodes[0]);

    EVdfg_link_dest(source, router);

    for (i=1; i < node_count; i++) {
	EVdfg_stone terminal = EVdfg_create_sink_stone(test_dfg,"simple_handler");
	EVdfg_link_port(router, i-1, terminal);
	EVdfg_assign_node(terminal, nodes[i]);
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
	int i;
	simple_rec rec;
	generate_simple_record(&rec);
	for (i=0 ; i < repeat_count; i++) {
	    /* encode shutdown in fields */
	    rec.integer_field = repeat_count / out_count;
	    EVsubmit(source_handle, &rec, NULL);
	}
    }

    status = EVclient_wait_for_shutdown(test_client);
    free(str_contact);
    EVfree_source(source_handle);
    wait_for_children(nodes);

    CManager_close(cm);
    for (i=0; i < node_count; i++) {
	free(nodes[i]);
    }
    free(nodes);
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
	/* submit will be quietly ignored if source is not active */
	EVsubmit(src, &rec, NULL);
    }
    EVfree_source(src);
    status = EVclient_wait_for_shutdown(test_client);
    CManager_close(cm);
    return status;
}
