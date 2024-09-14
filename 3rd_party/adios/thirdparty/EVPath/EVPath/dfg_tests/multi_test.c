#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "config.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;


typedef struct _rec_a {
    int a_field;
} rec_a, *rec_a_ptr;

typedef struct _rec_b {
    int b_field;
} rec_b, *rec_b_ptr;

typedef struct _rec_c {
    int c_field;
} rec_c, *rec_c_ptr;

static FMField a_field_list[] =
{
    {"a_field", "integer",
     sizeof(int), FMOffset(rec_a_ptr, a_field)},
    {NULL, NULL, 0, 0}
};

static FMField b_field_list[] =
{
    {"b_field", "integer",
     sizeof(int), FMOffset(rec_b_ptr, b_field)},
    {NULL, NULL, 0, 0}
};

static FMField c_field_list[] =
{
    {"c_field", "integer",
     sizeof(int), FMOffset(rec_c_ptr, c_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec a_format_list[] =
{
    {"a_rec", a_field_list, sizeof(rec_a), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec b_format_list[] =
{
    {"b_rec", b_field_list, sizeof(rec_b), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec c_format_list[] =
{
    {"c_rec", c_field_list, sizeof(rec_c), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescList queue_list[] = {a_format_list, b_format_list, c_format_list, NULL};

static int repeat_count = 100;

static
void
generate_a_record(rec_a_ptr event)
{
    /* always even */
    event->a_field = ((int) rand() % 50) * 2;
}

static
void
generate_b_record(rec_b_ptr event)
{
    /* always odd */
    event->b_field = ((int) rand() % 50) * 2 + 1;
}

static
int
output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    static int message_count = 0;
    rec_c_ptr event = vevent;
    (void)cm;
    if (event->c_field % 2 != 1) {
	printf("Received record should be odd, got %d\n", event->c_field);
    }
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	c_field = %d\n", event->c_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
    message_count++;
    if (message_count == repeat_count/2) {
	EVclient_shutdown(test_client, 0);
    }
    return 0;
}

static char *trans = "{\n\
    int found = 0;\n\
    a_rec *a;\n\
    b_rec *b;\n\
    c_rec c;\n\
    if (EVpresent(a_rec_ID, 0)) {\n\
        a = EVdata_a_rec(0); ++found;\n\
    }\n\
    if (EVpresent(b_rec_ID, 0)) {\n\
        b = EVdata_b_rec(0); ++found;\n\
    }\n\
    if (found == 2) {\n\
        c.c_field = a.a_field + b.b_field;\n\
        if (!EVpresent_b_rec(0))\n\
            printf(\"??? <1> not present (1)\\n\");\n\
        EVdiscard_a_rec(0);\n\
        if (!EVpresent_b_rec(0))\n\
            printf(\"??? <2> not present (1)\\n\");\n\
        EVdiscard_b_rec(0);\n\
        EVsubmit(0, c);\n\
    }\n\
}\0\0";

static void
data_free(void *event_data, void *client_data)
{
    (void) client_data;
    free(event_data);
}



extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", "c", "d", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone srca, srcb, multiq, sink;
    EVsource a_handle, b_handle;
    char * q_action_spec;
    int count, i;
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

    a_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, a_format_list,
					   data_free, NULL);
    b_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, b_format_list,
					   data_free, NULL);
    (void) EVclient_register_source("a_source", a_handle);
    source_capabilities = EVclient_register_source("b_source", b_handle);
    sink_capabilities = EVclient_register_sink_handler(cm, "c_output_handler", c_format_list,
				(EVSimpleHandlerFunc) output_handler, NULL);

/*
**  DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);
    test_dfg = EVdfg_create(test_master);

    srca = EVdfg_create_source_stone(test_dfg, "a_source");
    srcb = EVdfg_create_source_stone(test_dfg, "b_source");
    sink = EVdfg_create_sink_stone(test_dfg, "c_output_handler");
    q_action_spec = create_multityped_action_spec(queue_list, trans);
    multiq = EVdfg_create_stone(test_dfg, q_action_spec);
    free(q_action_spec);
    EVdfg_link_dest(srca, multiq);
    EVdfg_link_dest(srcb, multiq);
    EVdfg_link_port(multiq, 0, sink);

    EVdfg_assign_node(srca, "a");
    EVdfg_assign_node(srcb, "b");
    EVdfg_assign_node(multiq, "c");
    EVdfg_assign_node(sink, "d");

    EVdfg_realize(test_dfg);

/* We're node 0 in the DFG */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

/* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    
    count = repeat_count;
    for (i=0; i < count/2 ; i++) {
	if (EVclient_source_active(a_handle)) {
	    rec_a_ptr a = malloc(sizeof(*a));
	    generate_a_record(a);
	    if (quiet <=0) {printf("submitting a -> %d\n", a->a_field);}
	    EVsubmit(a_handle, a, NULL);
	}
	if (EVclient_source_active(b_handle)) {
	    rec_b_ptr b = malloc(sizeof(*b));
	    generate_b_record(b);
	    if (quiet <=0) {printf("submitting b -> %d\n", b->b_field);}
	    EVsubmit(b_handle, b, NULL);
	}
    }
    CMsleep(cm, 1);

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);

    EVfree_source(b_handle);
    EVfree_source(a_handle);
    free(str_contact);
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource a_handle, b_handle;
    int count, i;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  evtest  <nodename> <mastercontact>\n");
	exit(1);
    }

    sink_capabilities = EVclient_register_sink_handler(cm, "c_output_handler", c_format_list,
				(EVSimpleHandlerFunc) output_handler, NULL);
    a_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, a_format_list,
					   data_free, NULL);
    b_handle = EVcreate_submit_handle_free(cm, DFG_SOURCE, b_format_list,
					   data_free, NULL);
    (void) EVclient_register_source("a_source", a_handle);
    source_capabilities = EVclient_register_source("b_source", b_handle);

    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);

    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    count = repeat_count;

    for (i=0; i < count/2 ; i++) {
	if (EVclient_source_active(a_handle)) {
	    rec_a_ptr a = malloc(sizeof(*a));
	    generate_a_record(a);
	    if (quiet <=0) {printf("submitting a -> %d\n", a->a_field);}
	    EVsubmit(a_handle, a, NULL);
	}
	if (EVclient_source_active(b_handle)) {
	    rec_b_ptr b = malloc(sizeof(*b));
	    generate_b_record(b);
	    if (quiet <=0) {printf("submitting b -> %d\n", b->b_field);}
	    EVsubmit(b_handle, b, NULL);
	}
    }
    status = EVclient_wait_for_shutdown(test_client);
    EVfree_source(b_handle);
    EVfree_source(a_handle);
    return status;
}
