/*
 *   Test queue manipulation techniques in EVPath queue stones.  In
 *   particular, this test uses messages with anticipated types, as well as
 *   testing manipulation of messages whose types were unanticipated (WRT
 *   the queue stone handler) and which remain anonymous.  The idea is that
 *   we'll dump a bunch of messages into a queue, then send 'command'
 *   messages to cause certain manipulations of the queue, releasing some of
 *   those queued messages.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "ev_dfg.h"
#include "test_support.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#else
#include <sys/wait.h>
#endif

static int status;
static EVclient test_client;


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

typedef struct _rec_b {
    int b_field;
} rec_b, *rec_b_ptr;

static FMField b_field_list[] =
{
    {"b_field", "integer",
     sizeof(int), FMOffset(rec_b_ptr, b_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec b_format_list[] =
{
    {"b_rec", b_field_list, sizeof(rec_b), NULL},
    {NULL, NULL, 0, NULL}
};

typedef struct _rec_c {
    int c_field;
} rec_c, *rec_c_ptr;

static FMField c_field_list[] =
{
    {"c_field", "integer",
     sizeof(int), FMOffset(rec_c_ptr, c_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec c_format_list[] =
{
    {"c_rec", c_field_list, sizeof(rec_c), NULL},
    {NULL, NULL, 0, NULL}
};

typedef struct _rec_command {
    int command_field;
} rec_command, *rec_command_ptr;

static FMField command_field_list[] =
{
    {"command_field", "integer",
     sizeof(int), FMOffset(rec_command_ptr, command_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec command_format_list[] =
{
    {"command", command_field_list, sizeof(rec_command), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec anonymous_format_list[] =
{
    {"anonymous", NULL, 0, NULL},
    {NULL, NULL, 0, NULL}
};

/*
 *  'b' and 'c' are unknown to the queue stone
 */
static FMStructDescList queue_list[] = {a_format_list, command_format_list, anonymous_format_list, NULL};

#define EVENT_COUNT 20
#define COMMAND_COUNT 2

static int sequence_num = 0;

static
void
generate_a_record(rec_a_ptr event)
{
    event->a_field = ((int) lrand48() % 50) * 2 + 10000*++sequence_num;
}

static
void
generate_b_record(rec_b_ptr event)
{
    /* always odd */
    event->b_field = ((int) lrand48() % 50) * 2 + 1 * 2 + 10000*++sequence_num;
}

static
void
generate_c_record(rec_c_ptr event)
{
    /* always odd */
    event->c_field = ((int) lrand48() % 50) * 2 + 1 * 2 + 10000*++sequence_num;
}

static int expected_sequences[20]  = {10,1,6,7,20,18,16,15,14,12, 19,17,13,11,9,2,3,4,5,8};
static int received_sequences[EVENT_COUNT];

static int message_count = 0;

static void check_termination(int sequence)
{
    int i;
    int result = 0;
    received_sequences[message_count++] = sequence/10000;
    if (message_count < EVENT_COUNT) return;
    for (i=0; i < EVENT_COUNT; i++) {
	if (expected_sequences[i] != received_sequences[i]) {
	    printf("In slot %d, expected event %d, got %d\n", i, expected_sequences[i], received_sequences[i]);
	    result=1;
	}
    }
    EVclient_shutdown(test_client, result);
}

static
int
a_output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    rec_a_ptr event = vevent;
    (void)cm;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	a_field = %d\n", event->a_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    check_termination(event->a_field);
    return 0;
}

static
int
b_output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    rec_b_ptr event = vevent;
    (void)cm;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	b_field = %d\n", event->b_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    check_termination(event->b_field);
    return 0;
}

static
int
c_output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    rec_c_ptr event = vevent;
    (void)cm;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	c_field = %d\n", event->c_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    check_termination(event->c_field);
    return 0;
}

static char *trans = "{\n\
    int found = 0;\n\
    a_rec *a;\n\
    printf(\"==== Event count is a_rec = %d, b_rec = %d, command = %d, anon = %d\\n\", EVcount(0), EVcount(1), EVcount_command(), EVcount_anonymous());\n\
    if (EVcount_command() > 0) {\n\
	command *c = EVdata_command(0);\n\
	if (c->command_field == 1) {\n\
	    EVdiscard_and_submit_a_rec(0, EVcount_a_rec()-1);\n	\
	    while(EVcount_a_rec()) {\n\
	        printf(\"==== Submit and discard item 0(oldest in queue), count is %d\\n\", EVcount_a_rec());\n\
	        EVdiscard_and_submit_a_rec(0, 0);\n\
	    }\n\
	} else {\n\
	    int half = EVcount_anonymous() / 2;\n\
	    printf(\"==== There are %d anonymous events, discard half (%d) from top down \\n\", EVcount_anonymous(), half);\n\
	    while (half--) {\n\
	        printf(\"==== Submit and discard anonymous  %d, count is %d\\n\", EVcount_anonymous()-1,EVcount_anonymous());\n\
		EVdiscard_and_submit_anonymous(0, EVcount_anonymous()-1);\n\
	    }\n\
	    half = EVcount_full() / 2;\n					\
	    printf(\"==== There are %d total events, discard half (%d) from top (recent) down (skipping command)\\n\", EVcount_full(), half);\n\
	    while (half--) {\n\
	        printf(\"==== Submit and discard full  %d, count is %d\\n\", EVcount_full()-1,EVcount_full());\n\
		EVdiscard_and_submit_full(0, EVcount_full()-2);\n\
	    }\n\
	    while(EVcount_a_rec()) {\n\
	        printf(\"==== Submit and discard a_rec 0, count is %d\\n\", EVcount_a_rec());\n\
	        EVdiscard_and_submit_a_rec(0, 0);\n\
	    }\n\
	    while(EVcount_anonymous()) {\n\
	        printf(\"==== Submit and discard anonymous item 0, count is %d\\n\", EVcount_anonymous());\n\
	        EVdiscard_and_submit_anonymous(0, 0);\n\
	    }\n\
	}\n\
	EVdiscard_command(0);\n\
    }\n\
}\0\0";

static void
data_free(void *event_data, void *client_data)
{
    (void) client_data;
    free(event_data);
}

static int
raw_handler(CManager cm, void *vevent, int len, void *client_data,
	    attr_list attrs)
{
    static FFSContext c = NULL;
    FFSTypeHandle f;
    simple_rec incoming;
    (void)len;
    if (c == NULL) {
	c = create_FFSContext();
    }
    
    f = FFSTypeHandle_from_encode(c, vevent);
    if (!FFShas_conversion(f)) {
	establish_conversion(c, f, simple_format_list);
    }
    FFSdecode_to_buffer(c, vevent, &incoming);
    return c_output_handler(cm, (void*) &incoming, client_data, attrs);
}

extern void
generate_event_sequence(EVsource a_handle, EVsource b_handle, EVsource c_handle, EVsource command_handle)
{
/*
  generate event sequence:  
0: a=1
1: c=2
2: c=3
3: b=4
4: b=5
5: a=6
6: a=7
7: b=8
8: c=9
9: a=10
10: command = 1
11: a=11
12: c=12
13: a=13
14: b=14
15: b=15
16: c=16
17: a=17
18: b=18
19: a=19
20: c=20
21: command = 2

on command = 1 release top a, then all a from 0
0: a=10
1: a=1
2: a=6
3: a=7
release half (6) of anonymous from top
4: c=20    
5: b=18
6: c=16
7: b=15
8: b=14
9: c=12
There are 11 total events, discard half (5) from top (recent) down (skipping command)
10: a=19
11: a=17
12: a=13
13: a=11
14: c=9
submit any remaining a from bottom (oldest)
submit any remaining anon from bottom (oldest)
15: c=2
16: c=3
17: b=4
18: b=5
19: b=8
*/
    /* we know the sources are here, just do submits */
    int i;
    for (i=0; i < EVENT_COUNT + COMMAND_COUNT; i++) {
	switch (i) {
	case 0:
	case 5:
	case 6:
	case 9:
	case 11:
	case 13:
	case 17:
	case 19:{
	    rec_a_ptr a = malloc(sizeof(*a));
	    a->a_field = 1;
	    generate_a_record(a);
	    
	    if (quiet <=0) {printf("submitting a -> %d\n", a->a_field);}
	    EVsubmit_general(a_handle, a, (EVFreeFunction)free, NULL);
	    break;
	}
	case 3:
	case 4:
	case 7: 
	case 14:
	case 15:
	case 18: {
	    rec_b_ptr b = malloc(sizeof(*b));
	    b->b_field = 2;
	    generate_b_record(b);
	    
	    if (quiet <=0) {printf("submitting b -> %d\n", b->b_field);}
	    EVsubmit_general(b_handle, b, (EVFreeFunction)free, NULL);

	    break;
	}
	case 1:
	case 2:
	case 8:	
	case 12:
	case 16:
	case 20: {
	    rec_c_ptr c = malloc(sizeof(*c));
	    c->c_field = 3;
	    generate_c_record(c);
	    
	    if (quiet <=0) {printf("submitting c -> %d\n", c->c_field);}
	    EVsubmit_general(c_handle, c, (EVFreeFunction)free, NULL);

	    break;
	}
	case 21: {
	    rec_command_ptr com = malloc(sizeof(*com));
	    com->command_field = 2;
	    EVsubmit(command_handle, com, NULL);
	    break;
	}
	case 10: {
	    rec_command_ptr com = malloc(sizeof(*com));
	    com->command_field = 1;
	    EVsubmit(command_handle, com, NULL);
	    break;
	}
	}
    }
}

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", "c", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone srca, srcb, srcc, srccommand, multiq, sink;
    EVsource a_handle, b_handle, c_handle, command_handle;
    char * q_action_spec;
    int node_count = 3;
    char *src_node = "a";
    char *queue_node = "b";
    char *terminal_node = "c";
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    (void)argc; (void)argv;
    cm = CManager_create();
    if (argc == 1) {
	if (!sscanf(argv[0], "%d", &node_count)) printf("arg \"%s\" not understood\n", argv[0]);
	if ((node_count < 1) || (node_count > 3)) {
	    printf("bad node count %d\n", node_count);
	    node_count = 3;
	}
	nodes[node_count] = NULL;  /* slice off last few nodes */
	switch (node_count) {
	case 1:
	    src_node = nodes[0];
	    queue_node = nodes[0];
	    terminal_node = nodes[0];
	    break;
	case 2:
	    src_node = nodes[1];
	    queue_node = nodes[0];
	    terminal_node = nodes[1];
	    break;
	case 3:
	    src_node = nodes[0];
	    queue_node = nodes[1];
	    terminal_node = nodes[2];
	    break;
	} 
    }
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    a_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, a_format_list,
					   data_free, NULL);
    b_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, b_format_list,
					   data_free, NULL);
    c_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, c_format_list,
					   data_free, NULL);
    command_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, command_format_list,
					   data_free, NULL);
    (void) EVclient_register_source("a_source", a_handle);
    (void) EVclient_register_source("b_source", b_handle);
    (void) EVclient_register_source("c_source", c_handle);
    source_capabilities = EVclient_register_source("command_source", command_handle);
    (void) EVclient_register_sink_handler(cm, "a_output_handler", a_format_list,
						       (EVSimpleHandlerFunc) a_output_handler, NULL);
    (void) EVclient_register_sink_handler(cm, "b_output_handler", b_format_list,
				(EVSimpleHandlerFunc) b_output_handler, NULL);
    (void) EVclient_register_sink_handler(cm, "c_output_handler", c_format_list,
				(EVSimpleHandlerFunc) c_output_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "raw_output_handler", NULL,
				(EVSimpleHandlerFunc) raw_handler, NULL);
/*
**  MASTER AND DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);
    test_dfg = EVdfg_create(test_master);

    srca = EVdfg_create_source_stone(test_dfg, "a_source");
    srcb = EVdfg_create_source_stone(test_dfg, "b_source");
    srcc = EVdfg_create_source_stone(test_dfg, "c_source");
    srccommand = EVdfg_create_source_stone(test_dfg, "command_source");
    sink = EVdfg_create_sink_stone(test_dfg, "a_output_handler");
    EVdfg_add_sink_action(sink, "b_output_handler");
    EVdfg_add_sink_action(sink, "c_output_handler");
    q_action_spec = create_multityped_action_spec(queue_list, trans);
    multiq = EVdfg_create_stone(test_dfg, q_action_spec);
    free(q_action_spec);
    EVdfg_link_dest(srca, multiq);
    EVdfg_link_dest(srcb, multiq);
    EVdfg_link_dest(srcc, multiq);
    EVdfg_link_dest(srccommand, multiq);
    EVdfg_link_port(multiq, 0, sink);

    EVdfg_assign_node(srca, src_node);
    EVdfg_assign_node(srcb, src_node);
    EVdfg_assign_node(srcc, src_node);
    EVdfg_assign_node(srccommand, src_node);
    EVdfg_assign_node(multiq, queue_node);
    EVdfg_assign_node(sink, terminal_node);

    EVdfg_realize(test_dfg);

/* We're node 0 in the DFG */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

/* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    if (EVclient_source_active(a_handle)) {
	generate_event_sequence(a_handle, b_handle, c_handle, command_handle);
    }
    CMsleep(cm, 1);

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);
    free(str_contact);
    EVfree_source(command_handle);
    EVfree_source(a_handle);
    EVfree_source(b_handle);
    EVfree_source(c_handle);
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource a_handle, b_handle, c_handle, command_handle;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  evtest  <nodename> <mastercontact>\n");
	exit(1);
    }

    (void) EVclient_register_sink_handler(cm, "a_output_handler", a_format_list,
				(EVSimpleHandlerFunc) a_output_handler, NULL);
    (void) EVclient_register_sink_handler(cm, "b_output_handler", b_format_list,
				(EVSimpleHandlerFunc) b_output_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "c_output_handler", c_format_list,
				(EVSimpleHandlerFunc) c_output_handler, NULL);
    a_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, a_format_list,
					   data_free, NULL);
    b_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, b_format_list,
					   data_free, NULL);
    c_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, c_format_list,
					   data_free, NULL);
    command_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, command_format_list,
					   data_free, NULL);
    (void) EVclient_register_source("a_source", a_handle);
    (void) EVclient_register_source("b_source", b_handle);
    (void) EVclient_register_source("c_source", c_handle);
    source_capabilities = EVclient_register_source("command_source", command_handle);

    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);

    if (EVclient_source_active(a_handle)) {
	generate_event_sequence(a_handle, b_handle, c_handle, command_handle);
    }
    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    /* sources are with the master, don't do submits here */

    return EVclient_wait_for_shutdown(test_client);
}
