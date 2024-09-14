#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "config.h"
#include "ev_dfg.h"
#include "test_support.h"

CManager cm;

static int status;
static EVdfg dfg;
static int base=2;
static int nbase = 2;
static int level_count = 3;

static EVdfg_stone *tmp;
static EVdfg_stone *last;

const int reconfig_node_count = 3;
const int num_terminals = 2;

char *str_contact;
char **reconfig_list = NULL;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    (void)client_data;
    static int count;
    checksum_simple_record(event, attrs, quiet);
    if (++count == 300) {
		EVclient_shutdown(dfg, 0);
    }
    printf("\nStatic configuration working, ready for reconfiguration! All the best!\n");
    fflush(stdout);
    return 0;
}


/* ****** Algorithm ******
 
 - The initial tree will consist of two nodes, root and one child
 
 - First child freezes the dataflow and deletes the link between the root and it's only child
 - First child then forks n number of children which will now attempt to join the dfg
 - First child does not call realize, dataflow is frozen
 
 - Every node joining in creates stones and links them appropriately
 - According to the number of the node joining in, it will either be an intermediate node or a leaf node
 - Leaf nodes will have terminal actions associated with them and intermidate nodes will have split actions (as of now)
 - When appropriate number of nodes have joined in, the graph is realized
 - The above three points can be implemented by basically unrolling the loops in the static tree_test code
 
 - Timing measurements can be made by starting the timer in the node_register_handler only when the 2nd node for reconfiguration joins in
 - The timer will be stopped when the graph has been realized, mostly inside node_register_handler
 
 ****** ********* ****** */

static void
join_handler(EVdfg dfg, char *identifier, void* available_sources, void *available_sinks) 
{
	
    static int joined_node_count;
    static int i = 2;
    static int j;
    static int k;
	
    char *thandle;
    char *canon_name = malloc(20 * sizeof(char));
	
    if (strcmp(identifier, "origin") == 0) {
	EVdfg_stone terminal = NULL;
	
	canon_name = strdup("origin");
	EVmaster_assign_canonical_name(dfg, identifier, canon_name);
	
	EVclient_register_sink_handler(cm, "thandler", simple_format_list, (EVSimpleHandlerFunc) simple_handler);
	terminal = EVdfg_create_sink_stone(dfg, "thandler");
	
	EVdfg_link_port(last[0], 0, terminal);
	EVdfg_assign_node(last[0], "origin");
	EVdfg_assign_node(terminal, "origin");
	
	++joined_node_count;
	EVdfg_realize(dfg);
    } else {
	if (strcmp(identifier, "forker") == 0) {
	    //      EVdfg_reconfig_delete_link(dfg, 0, 1);
	    
	    //      canon_name = malloc(20 * sizeof(char));
	    ++joined_node_count;
	    
	    //      delayed_fork_children(cm, &reconfig_list[0], str_contact, 5);
	} else {
	    if (i < level_count) {
		//        if (j < nbase) {
		tmp[j] = EVdfg_create_stone(dfg, NULL);
		
		sprintf(canon_name, "client%d", joined_node_count);
		EVmaster_assign_canonical_name(dfg, identifier, canon_name);
		EVdfg_assign_node(tmp[j], canon_name);
		
		EVdfg_reconfig_link_port(last[j / base], j % base, tmp[j], NULL);
		
		++joined_node_count;
		++j;
		//        }
		if (j >= nbase) {
		    for (j = 0; j < nbase; ++j) {
			last[j] = tmp[j];
		    }
		    nbase *= base;
		    j = 0;
		    ++i;
		}
	    } else {
		thandle = strdup("thandler");
		//	sprintf(thandle, "thandler%d", joined_node_count);
		EVclient_register_sink_handler(cm, thandle, simple_format_list, (EVSimpleHandlerFunc) simple_handler);
		tmp[k] = EVdfg_create_sink_stone(dfg, thandle);
		
		sprintf(canon_name, "terminal_node%d", joined_node_count);
		EVmaster_assign_canonical_name(dfg, identifier, canon_name);
		EVdfg_assign_node(tmp[k], canon_name);
		
		EVdfg_reconfig_link_port(last[k / base], k % base, tmp[k], NULL);
		
		++k;
		++joined_node_count;
		if (k >= nbase) {
		    EVdfg_realize(dfg);
		}
	    }
	}
    }
}


extern int
be_test_master(int argc, char **argv)
{
    attr_list contact_list;
    EVsource source_handle;
    int last_row_size;
    int i;
    char **nodes;
    
    nodes = malloc(2 * sizeof(nodes[0]));
	
    nodes[0] = strdup("origin");
    nodes[1] = NULL;
	
    last_row_size = (int) pow(base, level_count - 1);
	
    tmp = malloc(last_row_size * sizeof(EVdfg_stone));
    last = malloc(last_row_size * sizeof(EVdfg_stone));
    
    cm = CManager_create();
    CMlisten(cm);
    contact_list = CMget_contact_list(cm);
    str_contact = attr_list_to_string(contact_list);
    
    printf("\nMaster's contact = %s\n", str_contact);
    fflush(stdout);
	
    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE ,simple_format_list);
    source_capabilities = EVclient_register_source("master_source", source_handle);
	
    dfg = EVdfg_create(cm);
    EVmaster_node_join_handler(dfg, join_handler);;
	
    last[0] = EVdfg_create_source_stone(dfg, "master_source");
	
    /* pprabhu: creating list of reconfiguration node names */
    reconfig_list = malloc(sizeof(reconfig_list[0]) * (reconfig_node_count + 2));
    for (i=0; i < (reconfig_node_count + 1); i++) {
        reconfig_list[i] = strdup("client");
    }
    reconfig_list[reconfig_node_count + 1] = NULL;
    
    EVdfg_join_dfg(dfg, "origin", str_contact);
    
    if (EVclient_ready_wait(dfg) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }
    
	
    if (EVclient_active_sink_count(dfg) == 0) {
	EVclient_ready_for_shutdown(dfg);
    }
	
    if (EVclient_source_active(source_handle)) {
	simple_rec rec;
	for (i = 0; i < 300; ++i) {
	    CMsleep(cm, 1);
	    generate_simple_record(&rec);
	    /* submit would be quietly ignored if source is not active */
	    EVsubmit(source_handle, &rec, NULL);
	}
    }
    
    status = EVclient_wait_for_shutdown(dfg);
    
    wait_for_children(nodes);
    
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource src;
    char *chandle;
	
    int i;
	
    cm = CManager_create();
	
    if (argc != 3) {
        printf("Child usage:  evtest  <nodename> <mastercontact>\n");
        exit(1);
    }
	
    dfg = EVdfg_create(cm);
	
	
    reconfig_list = malloc(sizeof(reconfig_list[0]) * (reconfig_node_count + 2));
    for (i=0; i < (reconfig_node_count + 1); i++) {
        reconfig_list[i] = strdup("client");
    }
    reconfig_list[reconfig_node_count + 1] = NULL;
	
    src = EVcreate_submit_handle(cm, DFG_SOURCE, simple_format_list);
    source_capabilities = EVclient_register_source("master_source", src);
    chandle = strdup("thandler");
	
    sink_capabilities = EVclient_register_sink_handler(cm,chandle, simple_format_list,
                                (EVSimpleHandlerFunc) simple_handler);
    EVdfg_join_dfg(dfg, argv[1], argv[2]);
	
    if (strcmp(argv[1], "forker") == 0) {
		test_fork_children(&reconfig_list[0], argv[2]);
    }
	
    EVclient_ready_wait(dfg);
    if (EVclient_active_sink_count(dfg) == 0) {
        EVclient_ready_for_shutdown(dfg);
    }
	
    if (EVclient_source_active(src)) {
        simple_rec rec;
        generate_simple_record(&rec);
        /* submit will be quietly ignored if source is not active */
        EVsubmit(src, &rec, NULL);
    }
    return EVclient_wait_for_shutdown(dfg);
}
