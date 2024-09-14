/*
 *  Test the ability of a CoD handler to reference other stones by name, see notes below for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

typedef struct _result {
    int tally;
} result, *result_ptr;

static FMField tally_list[] =
{
    {"tally", "integer",
     sizeof(int), FMOffset(result_ptr, tally)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec result_format_list[] =
{
    {"result", tally_list, sizeof(result), NULL},
    {NULL, NULL, 0, NULL}
};

static int status;
static EVclient test_client;

static
int
event_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    EVstone this_stone;
    attr_list stone_attrs;
    simple_rec_ptr event = vevent;
    int message_count = 0;
    static atom_t MSG_COUNT_ATOM = -1;
    (void)cm;
    (void)client_data;
    if (quiet <= 0) printf("In handler for stone %d\n", EVexecuting_stone(cm));
    if (MSG_COUNT_ATOM == -1) {
	MSG_COUNT_ATOM = attr_atom_from_string("MSG_COUNT");
    }
    checksum_simple_record(event, attrs, quiet);
    /* get handle to the terminal stone */
    this_stone = EVexecuting_stone(cm);
    stone_attrs = EVextract_attr_list(cm, this_stone);
    /* get the message count attribute */
    if (!get_int_attr(stone_attrs, MSG_COUNT_ATOM, &message_count)) {
	message_count = 0;
    }
    /* increment the message count */
    message_count++;
    /* set the message count attribute */
    set_int_attr(stone_attrs, MSG_COUNT_ATOM, message_count);
    return 0;
}

static
int
result_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    result_ptr event = vevent;
    (void)cm;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	tally = %d\n", event->tally);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (event->tally > 10) {
	static int shutdown_done = 0;
	if (!shutdown_done) {
	    EVclient_force_shutdown(test_client, 0);
	    shutdown_done++;
	}
    }
    return 0;
}

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone S1, T1, A1, T2;
    EVsource source_handle;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;
    char *A1_action_spec;
    attr_list T1_name_attrs;
    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_capabilities = EVclient_register_source("master_source", source_handle);
    (void)EVclient_register_sink_handler(cm, "event_handler", simple_format_list,
				(EVSimpleHandlerFunc) event_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "result_handler", result_format_list,
				(EVSimpleHandlerFunc) result_handler, NULL);

/*
**  MASTER and DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);

    /* create:
       - one source (S1) to generate a series of events.
       - one terminal stone (T1) to catch those events and update it's own
         attribute list.  This stone is "named", so it's attribute list is
         accessible from other stones.
       - one auto-stone (A1, a transform stone with automatic event submission
         at a particular frequency) which "monitors" the status of the
         terminal stone by watching its attributes.  When it notes a
         termination condition, it generates a "result" event.
       - a terminal stone (T2) that will catch the final event and report the result.

       In this arrangement, the initial source (S1) is linked to the first
       terminal stone (T1) and the autostone (A1) is linked to the final
       terminal stone (T2).  The actual monitoring occurs via CoD-based
       access from A1 to the attribute lists of T1.  These stones *MUST* be
       co-located for this to occur.  (Stone names are local to a node, and
       direct attribute access is only supported between stones on a single
       node.)
    */

char *COD_monitor = "{\n\
    attr_list tmp = EVget_stone_attrs(\"my_data_stone\");\n	\
    int message_count = attr_ivalue(tmp, \"MSG_COUNT\");\n	\
    output.tally = message_count;\n\
\n\
    return 1;\n\
}";


    test_dfg = EVdfg_create(test_master);
    S1 = EVdfg_create_source_stone(test_dfg, "master_source");
    EVdfg_assign_node(S1, "b");
    T1 = EVdfg_create_sink_stone(test_dfg, "event_handler");
    T1_name_attrs = create_attr_list();
    set_string_attr(T1_name_attrs, attr_atom_from_string("EVP_STONE_NAME"), strdup("my_data_stone"));
    EVdfg_set_attr_list(T1, T1_name_attrs);  /* this one will cause reconfig */
    free_attr_list(T1_name_attrs);
    T2 = EVdfg_create_sink_stone(test_dfg, "result_handler");

    A1_action_spec = create_transform_action_spec(NULL, result_format_list, COD_monitor);
    A1 = EVdfg_create_stone(test_dfg, A1_action_spec);
    free(A1_action_spec);
    EVdfg_assign_node(S1, "b");
    EVdfg_assign_node(T1, "a");
    EVdfg_assign_node(A1, "a");
    EVdfg_assign_node(T2, "b");
    EVdfg_link_dest(S1, T1);
    EVdfg_link_dest(A1, T2);
    EVdfg_enable_auto_stone(A1, 1, 0);  /* 1 sec intervals for monitoring */

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
	while (!EVclient_test_for_shutdown(test_client)) {
	    simple_rec rec;
	    generate_simple_record(&rec);
	    /* submit would be quietly ignored if source is not active */
	    EVsubmit(source_handle, &rec, NULL);
	    CMusleep(cm, 999999); /* wait for .1 sec */
	}
	/* exit if shutdown is imminent */
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
    source_capabilities = EVclient_register_source("master_source", src);
    (void)EVclient_register_sink_handler(cm, "event_handler", simple_format_list,
				(EVSimpleHandlerFunc) event_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "result_handler", result_format_list,
				(EVSimpleHandlerFunc) result_handler, NULL);
    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);

    if (EVclient_source_active(src)) {
	while (!EVclient_test_for_shutdown(test_client)) {
	    simple_rec rec;
	    generate_simple_record(&rec);
	    /* submit would be quietly ignored if source is not active */
	    EVsubmit(src, &rec, NULL);
	    CMusleep(cm, 999999); /* wait for .1 sec */
	}
	/* exit if shutdown is imminent */
    }

    EVfree_source(src);
    return EVclient_wait_for_shutdown(test_client);
}
