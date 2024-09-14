#include "config.h"
#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>

#ifdef HAVE_COD_H
#include "cod.h"
#else
#define cod_get_client_data(x,y) NULL
typedef void *cod_exec_context;
#endif
#include "evpath.h"
#include "cm_internal.h"
#include "ev_dfg.h"
#include "revpath.h"
#include "ev_dfg_internal.h"
#include "revp_internal.h"
#undef NDEBUG
#include <assert.h>

/*
 *  Some notes about the (current) operation of EVdfg.
 *  
 *  EVdfg operation is based on CM-level messages.  It currently supports
 *  only one DFG per process/CM (I.E. there is no DFG identifier passed in
 *  any operation, so the DFG is implied by the contacted CM).  DFGs pass
 *  through several states, including : Joining (while waiting for
 *  participating nodes), Starting (between last node join and actual
 *  operation), Running (normal operation), Reconfiguring (master changing
 *  DFG structure), and Shutdown (all nodes have voted for shutdown and
 *  we're killing the system).  (State is maintained in the master, and is
 *  not reliable in the client nodes.)
 *  
 *  Generally, most operation triggers have the same sort of structure.
 *  I.E. they check to see if they are the master.  If they are, the perform
 *  the operation directly, and if not they send a message to the master.
 *  The message handler then does the subroutine call that performs the
 *  operation directly.  So, the message handler and the wrapper subroutine
 *  are structured similarly for most operations.  However, some are handled
 *  differently.  In particular, notification of client departure
 *  (connection/node failure) and client join can't be handled except when
 *  we're in the Running state.  (Well, Join is handled in Joining too.)
 *  Also, voluntary reconfiguration, in addition to only being handled in
 *  the Running state, *must* be queued for later handling if it happens
 *  locally (because it's triggered by a call inside a CoD event handler,
 *  and reconfiguring while running inside a handler seems bad).  So, those
 *  messages must be queued for later handling.
 *
 *  Currently, we handle each of these messages on transition from another
 *  state into Running, or, in the case of voluntary reconfiguration if
 *  triggered locally, we handle it inside a CM delayed task (which should
 *  assure that we're at least at a message handling point).
 */


char *str_state[] = {"DFG_Joining", "DFG_Starting", "DFG_Running", "DFG_Reconfiguring", "DFG_Shutting_Down"};
static char *master_msg_str[] = {"DFGnode_join", "DFGdeploy_ack", "DFGshutdown_contrib", "DFGconn_shutdown", 
			  "DFGflush_reconfig", NULL};
char *stone_condition_str[EVstone_condition_last] = {"Undeployed", "Deployed", "Frozen", "Lost"};
char *ACT_string[] = {"ACT_no_op", "ACT_create", "ACT_add_action", "ACT_set_auto_period", "ACT_link_port", "ACT_link_dest", 
		      "ACT_unlink_port", "ACT_unlink_dest", "ACT_set_attrs", "ACT_destroy", "ACT_freeze", "ACT_unfreeze",
		      "ACT_assign_node", "ACT_create_bridge"};


static void handle_conn_shutdown(EVmaster master, EVmaster_msg_ptr msg);
static void handle_node_join(EVmaster master, EVmaster_msg_ptr msg);
static void handle_flush_reconfig(EVmaster master, EVmaster_msg_ptr);
static void handle_deploy_ack(EVmaster master, EVmaster_msg_ptr);
static void handle_deploy_ack_wrapper(EVmaster master, EVmaster_msg_ptr);
static void handle_shutdown_contrib(EVmaster master, EVmaster_msg_ptr);

static void
queue_master_msg(EVmaster master, void*vmsg, EVmaster_msg_type msg_type, CMConnection conn, int copy);
static void free_master_msg(EVmaster_msg *msg);
static void free_dfg_state(EVdfg_configuration state);

static void free_attrs_msg(EVflush_attrs_reconfig_ptr msg);

static FMField EVleaf_element_flds[] = {
    {"name", "string", sizeof(char*), FMOffset(leaf_element*, name)},
    {"FMtype", "string", sizeof(char*), FMOffset(leaf_element*, FMtype)},
    {NULL, NULL, 0, 0}
};

static FMField EVnode_join_msg_flds[] = {
    {"node_name", "string", sizeof(char*), FMOffset(EVnode_join_ptr, node_name)},
    {"contact_string", "string", sizeof(char*), FMOffset(EVnode_join_ptr, contact_string)},
    {"source_count", "integer", sizeof(int), FMOffset(EVnode_join_ptr, source_count)},
    {"sink_count", "integer", sizeof(int), FMOffset(EVnode_join_ptr, sink_count)},
    {"sources", "source_element[source_count]", sizeof(leaf_element), FMOffset(EVnode_join_ptr, sources)},
    {"sinks", "sink_element[sink_count]", sizeof(leaf_element), FMOffset(EVnode_join_ptr, sinks)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_node_join_format_list[] = {
    {"EVdfg_node_join", EVnode_join_msg_flds, sizeof(EVnode_join_msg), NULL},
    {"sink_element", EVleaf_element_flds, sizeof(leaf_element), NULL},
    {"source_element", EVleaf_element_flds, sizeof(leaf_element), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVready_msg_flds[] = {
    {"node_id", "integer", sizeof(int), FMOffset(EVready_ptr, node_id)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_ready_format_list[] = {
    {"EVdfg_ready", EVready_msg_flds, sizeof(EVready_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVdeploy_ack_msg_flds[] = {
    {"node_id", "string", sizeof(char*), FMOffset(EVdeploy_ack_ptr, node_id)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_deploy_ack_format_list[] = {
    {"EVdfg_deploy_ack", EVdeploy_ack_msg_flds, sizeof(EVdeploy_ack_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVshutdown_msg_flds[] = {
    {"value", "integer", sizeof(int), FMOffset(EVshutdown_ptr, value)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVclient_shutdown_format_list[] = {
    {"EVclient_shutdown", EVshutdown_msg_flds, sizeof(EVshutdown_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVshutdown_contribution_msg_flds[] = {
    {"value", "integer", sizeof(int), FMOffset(EVshutdown_contribution_ptr, value)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVclient_shutdown_contribution_format_list[] = {
    {"EVclient_shutdown_contribution", EVshutdown_contribution_msg_flds, sizeof(EVshutdown_contribution_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVconn_shutdown_msg_flds[] = {
    {"stone", "integer", sizeof(int), FMOffset(EVconn_shutdown_ptr, stone)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_conn_shutdown_format_list[] = {
    {"EVdfg_conn_shutdown", EVconn_shutdown_msg_flds, sizeof(EVconn_shutdown_msg), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVattr_stone_flds[] = {
    {"stone", "integer", sizeof(long), FMOffset(EVattr_stone_ptr, stone)},
    {"attr_str", "string", sizeof(char*), FMOffset(EVattr_stone_ptr, attr_str)},
    {NULL, NULL, 0, 0}
};

static FMField EVflush_attrs_reconfig_msg_flds[] = {
    {"reconfig", "integer", sizeof(int), FMOffset(EVflush_attrs_reconfig_ptr, reconfig)},
    {"count", "integer", sizeof(long), FMOffset(EVflush_attrs_reconfig_ptr, count)},
    {"attr_stone_list", "attr_stone_element[count]", sizeof(EVattr_stone_struct), FMOffset(EVflush_attrs_reconfig_ptr, attr_stone_list)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_flush_attrs_reconfig_format_list[] = {
    {"EVflush_attrs_reconfig", EVflush_attrs_reconfig_msg_flds, sizeof(EVflush_attrs_reconfig_msg), NULL},
    {"attr_stone_element", EVattr_stone_flds, sizeof(EVattr_stone_struct), NULL},
    {NULL, NULL, 0, NULL}
};

static FMField EVdfg_stone_flds[] = {
    {"global_stone_id", "integer", sizeof(int),
     FMOffset(deploy_msg_stone, global_stone_id)},
    {"attrs", "string", sizeof(char*),
     FMOffset(deploy_msg_stone, attrs)},
    {"period_secs", "integer", sizeof(int),
     FMOffset(deploy_msg_stone, period_secs)},
    {"period_usecs", "integer", sizeof(int),
     FMOffset(deploy_msg_stone, period_usecs)},
    {"out_count", "integer", sizeof(int),
     FMOffset(deploy_msg_stone, out_count)},
    {"out_links", "integer[out_count]", sizeof(int),
     FMOffset(deploy_msg_stone, out_links)},
    {"action", "string", sizeof(char*),
     FMOffset(deploy_msg_stone, action)},
    {"extra_actions", "integer", sizeof(int),
     FMOffset(deploy_msg_stone, extra_actions)},
    {"xactions", "string[extra_actions]", sizeof(char*),
     FMOffset(deploy_msg_stone, xactions)},
    {NULL, NULL, 0, 0}
};

static FMField EVdfg_deploy_msg_flds[] = {
    {"canonical_name", "string", sizeof(char*),
     FMOffset(EVdfg_deploy_ptr, canonical_name)},
    {"stone_count", "integer", sizeof(int),
     FMOffset(EVdfg_deploy_ptr, stone_count)},
    {"stone_list", "EVdfg_deploy_stone[stone_count]", sizeof(struct _EVdfg_msg_stone), FMOffset(EVdfg_deploy_ptr, stone_list)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec EVdfg_deploy_format_list[] = {
    {"EVdfg_deploy", EVdfg_deploy_msg_flds, sizeof(EVdfg_deploy_msg), NULL},
    {"EVdfg_deploy_stone", EVdfg_stone_flds, sizeof(struct _EVdfg_msg_stone), NULL},
    {NULL, NULL, 0, NULL}
};


/* msg action model
 *
 For each state/ for each master msg one of these possibilities:
	H - handle - dequeue and call handler (may change state, start over )
	U - unexpected - immediate error and discard (continue to next )
	I - ignore - discard (continue to next )
	L - leave_queued - (continue to next )
*/
static
char action_model[DFG_Last_State][DFGlast_msg] = {
/* join		deploy_ack	shutdown_contrib	conn_shutdown	flush_reconfig */
  {'H', 	'U',	 	'U', 			'U', 		'U'},/* state Joining */
  {'Q',		'H',		'H',			'Q',		'Q'},/* state Starting */
  {'H',		'U',		'H',			'H',		'H'},/* state Running */
  {'Q',		'H',		'H',			'Q',		'Q'},/* state Reconfiguring */
  {'U',		'U',		'U',			'I',		'U'}/* state Shutting Down */
};

typedef void (*master_msg_handler_func) (EVmaster master, EVmaster_msg_ptr msg);
static master_msg_handler_func master_msg_handler[DFGlast_msg] = {handle_node_join, handle_deploy_ack_wrapper, handle_shutdown_contrib, handle_conn_shutdown, handle_flush_reconfig};
static void dfg_master_msg_handler(CManager cm, CMConnection conn, void *vmsg, 
				   void *client_data, attr_list attrs);

static void
handle_queued_messages(CManager cm, void* vmaster)
{
    /* SHOULD */
    /*  1 - consolidate node failure messages (likely to get several for each node) and handle these first */
    /*  2 -  handle node join messages last */
    /* FOR THE MOMENT */
    /* just do everything in order */
    /* beware the the list might change while we're running a handler */
    EVmaster master = (EVmaster) vmaster;
    EVmaster_msg_ptr next;
    EVmaster_msg_ptr *last_ptr;

    if (master->queued_messages == NULL) return;
    assert(CManager_locked(cm));
    next = master->queued_messages;
    last_ptr = &master->queued_messages;
    while(next != NULL) {
	CMtrace_out(cm, EVdfgVerbose, "EVDFG handle_queued_messages -  master DFG state is %s\n", str_state[master->state]);
	switch (action_model[master->state][next->msg_type]) {
	case 'H':
	    CMtrace_out(cm, EVdfgVerbose, "Master Message is type %s, calling handler\n", master_msg_str[next->msg_type]);
	    *last_ptr = next->next;  /* remove msg from queue */
	    (*master_msg_handler[next->msg_type])(master, next);
	    free_master_msg(next);
	    next = master->queued_messages;   /* start from scratch in case state changed */
	    break;
	case 'U':
	    printf("Master Message is type %s, UNEXPECTED!  Discarding...\n", master_msg_str[next->msg_type]);
	    *last_ptr = next->next;  /* remove msg from queue */
	    free_master_msg(next);
	    next = *last_ptr;
	    break;
	case 'Q':
	    printf("Master Message is type %s, not appropriate now, leaving queued...\n", master_msg_str[next->msg_type]);
	    next = next->next;
	    break;
	default:
	    printf("Unexpected action type '%c', discarding\n", action_model[master->state][next->msg_type]);
	    /* falling through */
	case 'I':
	    *last_ptr = next->next;  /* remove msg from queue */
	    free_master_msg(next);
	    next = *last_ptr;
	    break;
	}
	CMtrace_out(cm, EVdfgVerbose, "EVDFG handle queued end loop -  master DFG state is now %s\n", str_state[master->state]);
    }	    
    CMtrace_out(cm, EVdfgVerbose, "EVDFG handle queued exiting -  master DFG state is now %s\n", str_state[master->state]);
}

static void
handle_queued_messages_lock(CManager cm, void* vmaster)
{
    CManager_lock(cm);
    handle_queued_messages(cm, vmaster);
    CManager_unlock(cm);
}

EVdfg_stone
INT_EVdfg_create_source_stone(EVdfg dfg, char *source_name)
{
    EVdfg_stone tmp;
    size_t len = strlen(source_name) + strlen("source:");
    char *act = malloc(len + 1);
    strcpy(&act[0], "source:");
    strcat(&act[0], source_name);
    tmp = INT_EVdfg_create_stone(dfg, &act[0]);
    free(act);
    return tmp;
}

extern void 
INT_EVdfg_add_sink_action(EVdfg_stone stone, char *sink_name)
{
    size_t len = strlen(sink_name) + strlen("sink:");
    char *act = malloc(len + 1);
    strcpy(&act[0], "sink:");
    strcat(&act[0], sink_name);
    INT_EVdfg_add_action(stone, &act[0]);
    free(act);
}

EVdfg_stone
INT_EVdfg_create_sink_stone(EVdfg dfg, char *sink_name)
{
    EVdfg_stone tmp;
    size_t len = strlen(sink_name) + strlen("sink:");
    char *act = malloc(len + 1);
    strcpy(&act[0], "sink:");
    strcat(&act[0], sink_name);
    tmp = INT_EVdfg_create_stone(dfg, &act[0]);
    free(act);
    return tmp;
}

static int
EVdfg_perform_act_on_state(EVdfg_configuration state, EVdfg_config_action act, int build_queue);
static void fdump_dfg_config_action(FILE* out, EVdfg_config_action a);
/*static void dump_dfg_config_action(EVdfg_config_action a)
{
    fdump_dfg_config_action(stdout, a);
    }*/

static void
EVdfg_add_act_to_queue(EVdfg_configuration state, EVdfg_config_action act)
{
    if (state->pending_action_queue == NULL) {
	state->pending_action_count = 0;
	state->pending_action_queue = malloc(sizeof(state->pending_action_queue[0]));
	state->pending_action_queue[state->pending_action_count++] = act;
	return;
    }
    state->pending_action_queue = realloc(state->pending_action_queue, 
					  sizeof(state->pending_action_queue[0]) * (state->pending_action_count+1));
    if (act.type == ACT_create_bridge) {
	/* insert at beginning */
	memmove(&state->pending_action_queue[1], state->pending_action_queue, 
		sizeof(state->pending_action_queue[0]) * state->pending_action_count);
	state->pending_action_queue[0] = act;
	state->pending_action_count++;
    } else {
	state->pending_action_queue[state->pending_action_count++] = act;
    }
}

void
INT_EVdfg_add_action(EVdfg_stone stone, char *action)
{
    EVdfg_config_action act;
    act.type = ACT_add_action;
    act.stone_id = stone->stone_id;
    act.u.create.action = action ? strdup(action) : NULL;
    EVdfg_perform_act_on_state(stone->dfg->working_state, act, 1 /* add to queue */);
}

EVdfg_stone
INT_EVdfg_create_stone(EVdfg dfg, char *action)
{
    EVdfg_stone stone = malloc(sizeof(struct _EVdfg_stone));
    EVdfg_config_action act;
    stone->dfg = dfg;
    stone->stone_id = 0x80000000 | dfg->stone_count++;
    act.type = ACT_create;
    act.stone_id = stone->stone_id;
    act.u.create.action = action ? strdup(action) : NULL;
    dfg->stones = realloc(dfg->stones, sizeof(dfg->stones[0]) * dfg->stone_count);
    dfg->stones[dfg->stone_count-1] = stone;
    EVdfg_perform_act_on_state(dfg->working_state, act, 1 /* add to queue */);
    return stone;
}

static
EVdfg_stone_state find_stone_state(int stone_id, EVdfg_configuration config)
{
    int i = 0;
    for (i=0; i < config->stone_count ; i++) {
	if (stone_id == config->stones[i]->stone_id) 
	    return config->stones[i];
    }
    return NULL;
}

static void
assign_actions_to_nodes(EVdfg_configuration config, EVmaster master)
{
    int i;
    for (i=0; i < config->pending_action_count; i++) {
	EVdfg_config_action act = config->pending_action_queue[i];
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	int node = stone->node;
	if (master->nodes[node].shutdown_status_contribution == STATUS_FAILED) {
	    if (CMtrace_on(master->cm, EVdfgVerbose)) {
		printf("Skipping action because node failed -> ");
		fdump_dfg_config_action(master->cm->CMTrace_file, act);
	    }
	    config->pending_action_queue[i].type = ACT_no_op;
	} else {
	    if (CMtrace_on(master->cm, EVdfgVerbose)) {
		fprintf(master->cm->CMTrace_file, "Assigning action to node %d -> ", node);
		fdump_dfg_config_action(master->cm->CMTrace_file, act);
	    }
	    config->pending_action_queue[i].node_for_action = node;
	}
    }
}
	
static void
remove_actions_for_node(EVdfg_configuration config, int node, EVmaster master)
{
    int i;
    for (i=0; i < config->pending_action_count; i++) {
	EVdfg_config_action act = config->pending_action_queue[i];
	if (act.node_for_action == node) {
	    if (CMtrace_on(master->cm, EVdfgVerbose)) {
		fdump_dfg_config_action(master->cm->CMTrace_file, act);
	    }
	    config->pending_action_queue[i].type = ACT_no_op;
	}
    }
}
	
static void
build_deploy_msg_for_node_stones(EVdfg_configuration config, int act_num, EVmaster master);

static void
perform_actions_on_nodes(EVdfg_configuration config, EVmaster master)
{
    CManager cm = master->cm;
    int i;
    for (i=0; i < config->pending_action_count; i++) {
	EVdfg_config_action act = config->pending_action_queue[i];
	int local = 0;
	CMConnection conn = NULL;
	if (CMtrace_on(cm, EVdfgVerbose)) {
	    fdump_dfg_config_action(master->cm->CMTrace_file, act);
	}
	if (master->nodes[act.node_for_action].self) {
	    local = 1;
	} else {
	    conn = master->nodes[act.node_for_action].conn;
	}
	switch(act.type) {
	case ACT_no_op: break;
	case ACT_create_bridge:
	case ACT_create:
	    /* build and send deploy msg for this stone and all stones created on it's node */
	    /* kill all actions related to this stone (add_action, link_port, set_attrs, assign_node) */
	    build_deploy_msg_for_node_stones(config, i, master);
	    break;
	case ACT_add_action: {
	    break;
	}
	case ACT_set_auto_period: {
	    break;
	}
	case ACT_link_port: {
	    if (local) {
		INT_EVstone_set_output(cm, act.stone_id, act.u.link.port, act.u.link.dest_id);
	    } else {
		INT_REVstone_set_output(conn, act.stone_id, act.u.link.port, act.u.link.dest_id);
	    }
	    break;
	}
	case ACT_link_dest: {
	    if (local) {
		INT_EVstone_add_split_target(cm, act.stone_id, act.u.link.dest_id);
	    } else {
		INT_REVstone_add_split_target(conn, act.stone_id, act.u.link.dest_id);
	    }
	    break;
	}
	case ACT_unlink_port: {
	    if (local) {
		INT_EVstone_set_output(cm, act.stone_id, act.u.link.port, -1);
	    } else {
		INT_REVstone_set_output(conn, act.stone_id, act.u.link.port, -1);
	    }
	    break;
	}
	case ACT_unlink_dest: {
	    if (local) {
		INT_EVstone_remove_split_target(cm, act.stone_id, act.u.link.dest_id);
	    } else {
		INT_REVstone_remove_split_target(conn, act.stone_id, act.u.link.dest_id);
	    }
	    break;
	}
	case ACT_set_attrs: {
	    if (local) {
		INT_EVset_attr_list(cm, act.stone_id, act.u.attrs.attrs);
	    } else {
		INT_REVset_attr_list(conn, act.stone_id, act.u.attrs.attrs);
	    }
	    break;
	}
	case ACT_assign_node: {
	    break;
	}
	case ACT_destroy: {
	    if (local) {
		INT_EVdestroy_stone(cm, act.stone_id);
	    } else {
		INT_REVdestroy_stone(conn, act.stone_id);
	    }
	    break;
	}
	case ACT_freeze: {
	    if (local) {
		INT_EVfreeze_stone(cm, act.stone_id);
	    } else {
		INT_REVfreeze_stone(conn, act.stone_id);
	    }
	    break;
	}
	case ACT_unfreeze: {
	    if (local) {
		INT_EVunfreeze_stone(cm, act.stone_id);
	    } else {
		INT_REVunfreeze_stone(conn, act.stone_id);
	    }
	    break;
	}
	default:
	    printf("Bad action in perform_action_on_nodes %d\n", act.type);
	    break;
	}
    }
}

int
EVdfg_perform_act_on_state(EVdfg_configuration config, EVdfg_config_action act, int build_queue)
{
    int bridge = 0;
    switch(act.type) {
    case ACT_create_bridge:
	bridge = 1;
	/* falling through */
    case ACT_create: {
	EVdfg_stone_state stone = malloc(sizeof(*stone));
	stone->node = -1;
	stone->bridge_stone = 0;
	stone->attrs = NULL;
	stone->period_secs = -1;
	stone->period_usecs = -1;
	stone->out_count = 0;
	stone->out_links = NULL;
	stone->in_count = 0;
	stone->in_links = NULL;
	stone->action_count = 1;
	stone->extra_actions = NULL;
	stone->condition = EVstone_Undeployed;
	stone->bridge_target = -1;
	stone->fetched_events = NULL;
	if (bridge) {
	    stone->bridge_stone = 1;
	    stone->stone_id = act.stone_id;
	    stone->action = act.u.bridge.action;
	    stone->bridge_target = act.u.bridge.target_id;
	    stone->node = act.node_for_action;
	} else {
	    stone->stone_id = act.stone_id;
	    stone->action = act.u.create.action;
	}
	if (config->stone_count == 0) {
	    config->stones = malloc(sizeof(config->stones[0]));
	} else {
	    config->stones = realloc(config->stones, 
				     sizeof(config->stones[0]) * (config->stone_count + 1));
	}
	config->stones[config->stone_count++] = stone;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_add_action: {
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	if (!stone) {
	    return 0;
	}
	if (stone->action == NULL) {
	    stone->action = act.u.create.action;
	    break;
	}
	if (stone->extra_actions == NULL) {
	    stone->extra_actions = malloc(sizeof(stone->extra_actions[0]));
	} else {
	    stone->extra_actions = realloc(stone->extra_actions,
					   stone->action_count * sizeof(stone->extra_actions[0]));
	}
	stone->extra_actions[stone->action_count - 1] = act.u.create.action;
	stone->action_count++;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_set_auto_period: {
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	if (!stone) {
	    return 0;
	}
	stone->period_secs = act.u.period.secs;
	stone->period_usecs = act.u.period.usecs;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_link_port: {
	EVdfg_stone_state src = find_stone_state(act.stone_id, config);
	EVdfg_stone_state dst = find_stone_state(act.u.link.dest_id, config);
	int i, in_link_found;
	if (!src) {
	    return 0;
	}
	if (src->out_count == 0) {
	    src->out_links = malloc(sizeof(src->out_links[0]) * (act.u.link.port+1));
	    memset(src->out_links, 0, sizeof(src->out_links[0]) * (act.u.link.port+1));
	    src->out_count = act.u.link.port + 1;
	} else if (src->out_count < act.u.link.port + 1) {
	    src->out_links = realloc(src->out_links,
				     sizeof(src->out_links[0]) * (act.u.link.port+1));
	    memset(&src->out_links[src->out_count], -1, sizeof(src->out_links[0]) * (act.u.link.port+1-src->out_count));
	    src->out_count = act.u.link.port + 1;
	}
	in_link_found =  0;
	for (i=0; i < dst->in_count; i++) {
	    if (dst->in_links[i] == act.stone_id) in_link_found = 1;
	}
	if (in_link_found == 0) {
	    if (dst->in_count == 0) {
		dst->in_links = malloc(sizeof(dst->in_links[0]));
		memset(dst->in_links, 0, sizeof(dst->in_links[0]));
		dst->in_count = 1;
	    } else {
		dst->in_links = realloc(dst->in_links,
					sizeof(dst->in_links[0]) * (dst->in_count+1));
		memset(&dst->in_links[dst->in_count], 0, sizeof(dst->in_links[0]) * (dst->in_count+1-dst->in_count));
		dst->in_count++;
	    }
	    dst->in_links[dst->in_count - 1] = act.stone_id;
	}
	if (build_queue) {
	    if (src->out_links[act.u.link.port] != -1) {
		EVdfg_config_action dact;
		dact.type = ACT_unlink_port;
		dact.node_for_action = src->node;
		dact.stone_id = src->stone_id;
		dact.u.link.port = act.u.link.port;
		EVdfg_perform_act_on_state(config, dact, build_queue);
	    }
	}
	src->out_links[act.u.link.port] = act.u.link.dest_id;
	if (build_queue) {
	    if (src->condition == EVstone_Deployed) {
		EVdfg_config_action fact;
		fact.type = ACT_freeze;
		fact.stone_id = src->stone_id;
		fact.node_for_action = src->node;
		EVdfg_add_act_to_queue(config, fact);
		src->condition = EVstone_Frozen;
	    }
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_link_dest: {
	EVdfg_stone_state src = find_stone_state(act.stone_id, config);
	EVdfg_stone_state dst = find_stone_state(act.u.link.dest_id, config);
	int i, in_link_found, out_link_found =  0;
	for (i=0; i < src->out_count; i++) {
	    if (src->out_links[i] == act.u.link.dest_id) out_link_found = 1;
	}
	if (out_link_found) {
	    return 0;
	}
	if (!src) {
	    return 0;
	}
	if (src->out_count == 0) {
	    src->out_links = malloc(sizeof(src->out_links[0]) * (src->out_count+1));
	    memset(src->out_links, 0, sizeof(src->out_links[0]) * (src->out_count+1));
	    src->out_count++;
	} else {
	    src->out_links = realloc(src->out_links,
				     sizeof(src->out_links[0]) * (src->out_count+1));
	    src->out_count++;
	}
	src->out_links[src->out_count-1] = act.u.link.dest_id;
	in_link_found = 0;
	for (i=0; i < dst->in_count; i++) {
	    if (dst->in_links[i] == act.stone_id) in_link_found = 1;
	}
	if (in_link_found == 0) {
	    if (dst->in_count == 0) {
		dst->in_links = malloc(sizeof(dst->in_links[0]));
		memset(dst->in_links, 0, sizeof(dst->in_links[0]));
		dst->in_count = 1;
	    } else {
		dst->in_links = realloc(dst->in_links,
					sizeof(dst->in_links[0]) * (dst->in_count+1));
		memset(&dst->in_links[dst->in_count], 0, sizeof(dst->in_links[0]) * (dst->in_count+1-dst->in_count));
		dst->in_count++;
	    }
	    dst->in_links[dst->in_count - 1] = act.stone_id;
	}
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_unlink_port: {
	EVdfg_stone_state src = find_stone_state(act.stone_id, config);
	EVdfg_stone_state dest;
	if (!src) {
	    return 0;
	}
	if (src->out_count <= act.u.link.port) return 0;
	if (src->out_links[act.u.link.port] == -1) {
	    /* unassigned, error to unlink */
	    return 0;
	}
	dest = find_stone_state(src->out_links[act.u.link.port], config);
	if (!dest) {
	    return 0;
	}
	if (dest->bridge_stone) {
	    /* this stone only exists because of this link, destroy it */
	    EVdfg_config_action dact;
	    dact.type = ACT_destroy;
	    dact.stone_id = dest->stone_id;
	    EVdfg_perform_act_on_state(config, dact, build_queue);
	}
	src->out_links[act.u.link.port] = -1;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_unlink_dest: {
	EVdfg_stone_state src = find_stone_state(act.stone_id, config);
	EVdfg_stone_state dest = find_stone_state(act.u.link.dest_id, config);
	int found = 0;
	int i = 0;
	if (!src) {
	    return 0;
	}
	if (src->out_count <= act.u.link.port) return 0;
	for (i=0; i < src->out_count; i++) {
	    if (src->out_links[i] == dest->stone_id) {
		/* remove this, move remaining down */
		memmove(&src->out_links[i], &src->out_links[i+1], 
			sizeof(src->out_links[0]) * (src->out_count - i - 1));
		found++;
	    } else {
		EVdfg_stone_state possible_bridge = find_stone_state(src->out_links[i], config);
		if (possible_bridge->bridge_stone && (possible_bridge->out_links[0] == dest->stone_id)) {
		    /* there was a bridge stone in between, but we found the link */
		    /* the bridge stone only exists because of this link, destroy it */
		    EVdfg_config_action dact;
		    dact.type = ACT_destroy;
		    dact.stone_id = possible_bridge->stone_id;
		    EVdfg_perform_act_on_state(config, dact, build_queue);
		    /* remove this link, move remaining down */
		    memmove(&src->out_links[i], &src->out_links[i+1], 
			    sizeof(src->out_links[0]) * (src->out_count - i - 1));
		    found++;
		    /* rewrite action to refer to link by position */
		    act.type = ACT_unlink_port;
		    act.u.link.port = i;
		}
	    }
	}
	if (found == 0) {
	    /* not found, error to unlink */
	    return 0;
	} 
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_set_attrs: {
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	if (!stone) {
	    return 0;
	}
	if (stone->attrs != NULL) {
	    free_attr_list(stone->attrs);
	}
	stone->attrs = act.u.attrs.attrs;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_assign_node: {
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	if (!stone) {
	    return 0;
	}
	stone->node = act.u.assign.dest_node;
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    case ACT_destroy: {
	EVdfg_stone_state stone = find_stone_state(act.stone_id, config);
	if (!stone) {
	    return 0;
	}
	if (build_queue) {
	    EVdfg_add_act_to_queue(config, act);
	}
	break;
    }
    default:
	switch(act.type) {
	case ACT_no_op:
	case ACT_freeze:
	case ACT_unfreeze:
	    printf("Bad action in perform_act_on_state %s (%d)\n", ACT_string[act.type], act.type);
	    break;
	default:
	    printf("Bad action in perform_act_on_state %d\n", act.type);
	    break;
	}
	return 0;
    }
    return 1;
}

static void
fdump_dfg_config_action(FILE *out, EVdfg_config_action act)
{
    switch(act.type) {
    case ACT_no_op:
	fprintf(out, "Action is NO_OP\n");
	break;
    case ACT_create_bridge:
	fprintf(out, "Action is CREATE_BRIDGE, stone_id = %x, target %x, action %s\n",
		act.stone_id, act.u.bridge.target_id, act.u.bridge.action);
	break;
    case ACT_create:
	fprintf(out, "Action is CREATE_STONE, stone_id = %x, action %s\n",
		act.stone_id, act.u.create.action);
	break;
    case ACT_add_action:
	fprintf(out, "Action is ADD_ACTION, stone_id = %x, action %s\n",
		act.stone_id, act.u.create.action);
	break;
    case ACT_set_auto_period:
	fprintf(out, "Action is SET_AUTO_PERIOD, stone_id = %x, secs %d, usecs %d\n",
		act.stone_id, act.u.period.secs, act.u.period.usecs);
	break;
    case ACT_link_port:
	fprintf(out, "Action is LINK_PORT, src_stone_id = %x, port %d, dest_stone_id %x\n",
		act.stone_id, act.u.link.port, act.u.link.dest_id);
	break;
    case ACT_link_dest:
	fprintf(out, "Action is LINK_DEST, src_stone_id = %x, dest_stone_id %x\n",
		act.stone_id, act.u.link.dest_id);
	break;
    case ACT_unlink_port:
	fprintf(out, "Action is UNLINK_PORT, src_stone_id = %x, port %d\n",
		act.stone_id, act.u.link.port);
	break;
    case ACT_unlink_dest:
	fprintf(out, "Action is UNLINK_DEST, src_stone_id = %x, dest %x\n",
		act.stone_id, act.u.link.dest_id);
	break;
    case ACT_set_attrs: 
	fprintf(out, "Action is SET_ATTRS, stone_id = %x, attrs ",
		act.stone_id);
	fdump_attr_list(out, act.u.attrs.attrs);
	break;
    case ACT_assign_node:
	fprintf(out, "Action is ASSIGN_NODE, stone_id = %x, node = %d\n",
		act.stone_id, act.u.assign.dest_node);
	break;
    case ACT_destroy:
	fprintf(out, "Action is DESTROY, stone_id = %x\n",
		act.stone_id);
	break;
    case ACT_freeze:
	fprintf(out, "Action is FREEZE, stone_id = %x\n",
		act.stone_id);
	break;
    case ACT_unfreeze:
	fprintf(out, "Action is UNFREEZE, stone_id = %x\n",
		act.stone_id);
	break;
    }
}

extern void 
INT_EVdfg_enable_auto_stone(EVdfg_stone stone, int period_sec, 
			int period_usec)
{
    EVdfg_config_action act;
    act.type = ACT_set_auto_period;
    act.stone_id = stone->stone_id;
    act.u.period.secs = period_sec;
    act.u.period.usecs = period_usec;
    (void)EVdfg_perform_act_on_state(stone->dfg->working_state, act, 1 /* add to queue */);
}


static void fdump_dfg_stone(FILE* out, EVdfg_stone_state s);

extern int
INT_EVdfg_link_port(EVdfg_stone src, int port, EVdfg_stone dest)
{
    EVdfg_config_action act;
    if (port < 0) return 0;
    act.type = ACT_link_port;
    act.stone_id = src->stone_id;
    act.u.link.port = port;
    act.u.link.dest_id = dest->stone_id;
    return EVdfg_perform_act_on_state(src->dfg->working_state, act, 1 /* add to queue */);
}

extern int
INT_EVdfg_link_dest(EVdfg_stone src, EVdfg_stone dest)
{
    EVdfg_config_action act;
    act.type = ACT_link_dest;
    act.stone_id = src->stone_id;
    act.u.link.dest_id = dest->stone_id;
    return EVdfg_perform_act_on_state(src->dfg->working_state, act, 1 /* add to queue */);
}

extern int
INT_EVdfg_unlink_port(EVdfg_stone src, int port)
{
    EVdfg_config_action act;
    if (port < 0) return 0;
    act.type = ACT_unlink_port;
    act.stone_id = src->stone_id;
    act.u.link.port = port;
    return EVdfg_perform_act_on_state(src->dfg->working_state, act, 1 /* add to queue */);
}

extern int
INT_EVdfg_unlink_dest(EVdfg_stone src, EVdfg_stone dest)
{
    EVdfg_config_action act;
    act.type = ACT_unlink_dest;
    act.stone_id = src->stone_id;
    act.u.link.dest_id = dest->stone_id;
    return EVdfg_perform_act_on_state(src->dfg->working_state, act, 1 /* add to queue */);
}

extern int
INT_EVdfg_set_attr_list(EVdfg_stone stone, attr_list attrs)
{
    EVdfg_config_action act;
    act.type = ACT_set_attrs;
    act.stone_id = stone->stone_id;
    act.u.attrs.attrs = attrs;
    add_ref_attr_list(attrs);
    return EVdfg_perform_act_on_state(stone->dfg->working_state, act, 1 /* add to queue */);
}

extern attr_list
INT_EVdfg_get_attr_list(EVdfg_stone stone)
{
    if (stone->dfg->deployed_state) {
	EVdfg_stone_state dstone = find_stone_state(stone->stone_id, stone->dfg->deployed_state);
	if (dstone) {
	    if (dstone->attrs) add_ref_attr_list(dstone->attrs);
	    return dstone->attrs;
	}
    }
    if (stone->dfg->working_state) {
	EVdfg_stone_state wstone = find_stone_state(stone->stone_id, stone->dfg->deployed_state);
	if (wstone) {
	    if (wstone->attrs) add_ref_attr_list(wstone->attrs);
	    return wstone->attrs;
	}
    }	
    return NULL;
}


static void check_all_nodes_registered(EVmaster master);
static void possibly_signal_shutdown(EVmaster master, int value, CMConnection conn);
static int new_shutdown_condition(EVclient client, CMConnection conn);

static void
enable_auto_stones(CManager cm, EVclient client)
{
    int i = 0;
    auto_stone_list *auto_list = client->pending_auto_list;
    client->pending_auto_list = NULL;
    CMtrace_out(cm, EVdfgVerbose, "ENABLING AUTO STONES, list is %p\n", auto_list);
    while (auto_list && auto_list[i].period_secs != -1) {
        /* everyone is ready, enable auto stones */
	CMtrace_out(cm, EVdfgVerbose, "auto stone %d, period %d sec, %d usec\n", auto_list[i].stone, auto_list[i].period_secs, auto_list[i].period_usecs);
	INT_EVenable_auto_stone(cm, auto_list[i].stone, auto_list[i].period_secs, auto_list[i].period_usecs);
	i++;
    }
    if (auto_list) free(auto_list);
}

static void
dfg_ready_handler(CManager cm, CMConnection conn, void *vmsg, 
		  void *client_data, attr_list attrs)
{
    EVclient client = client_data;
    EVready_ptr msg =  vmsg;
    (void) conn;
    (void) attrs;
    client->my_node_id = msg->node_id;
    CManager_lock(cm);
    enable_auto_stones(cm, client);
    if (client->ready_condition != -1) {
	CMtrace_out(cm, EVdfgVerbose, "Client DFG %p Node id %d is ready, signalling %d\n", client, client->my_node_id, client->ready_condition);
	INT_CMCondition_signal(cm, client->ready_condition);
    } else {
	CMtrace_out(cm, EVdfgVerbose, "Client DFG %p Node id %d got ready, reconfig done\n", client, client->my_node_id);
    }	
    CManager_unlock(cm);
}

static void 
handle_conn_shutdown(EVmaster master, EVmaster_msg_ptr msg)
{
    int stone = msg->u.conn_shutdown.stone;
    EVdfg_stone_state reporting_stone = find_stone_state(stone, master->dfg->deployed_state);
    EVdfg dfg = master->dfg;

    /* this stone is automatically frozen by EVPath */
    reporting_stone->condition = EVstone_Frozen;
    master->state = DFG_Reconfiguring;
    
    CMtrace_out(master->cm, EVdfgVerbose, "EVDFG conn_shutdown_handler -  master DFG state is now %s\n", str_state[master->state]);
    if (master->node_fail_handler != NULL) {
	int i;
	int target_stone = -1;
	char *failed_node = NULL;
	char *contact_str = NULL;
	CMtrace_out(master->cm, EVdfgVerbose, "IN CONN_SHUTDOWN_HANDLER\n");
	for (i=0; i< dfg->stone_count; i++) {
	    int j;
	    for (j = 0; j < dfg->deployed_state->stones[i]->out_count; j++) {
		if (dfg->deployed_state->stones[i]->out_links[j] == stone) {
		    int out_stone_id = dfg->deployed_state->stones[i]->out_links[j];
		    EVdfg_stone_state out_stone = find_stone_state(out_stone_id, dfg->deployed_state);
		    CMtrace_out(master->cm, EVdfgVerbose, "Found reporting stone as output %d of stone %d\n",
				j, i);
		    parse_bridge_action_spec(out_stone->action, 
					     &target_stone, &contact_str);
		    CMtrace_out(master->cm, EVdfgVerbose, "Dead stone is %d\n", target_stone);
		}
	    }
	}
	for (i=0; i< dfg->stone_count; i++) {
	    if (dfg->deployed_state->stones[i]->stone_id == target_stone) {
		int node = dfg->deployed_state->stones[i]->node;
		int j;
		CMtrace_out(master->cm, EVdfgVerbose, "Dead node is %d, name %s\n", node,
			    master->nodes[node].canonical_name);
		failed_node = master->nodes[node].canonical_name;
		master->nodes[node].shutdown_status_contribution = STATUS_FAILED;
		for (j=0; j< dfg->stone_count; j++) {
		    if (dfg->deployed_state->stones[j]->node == node) {;
			CMtrace_out(master->cm, EVdfgVerbose, "Dead node is %d, name %s\n", node,
				    master->nodes[node].canonical_name);
			dfg->deployed_state->stones[j]->condition = EVstone_Lost;
		    }
		}
	    }
	}
	CManager_unlock(master->cm);
	master->node_fail_handler(dfg, failed_node, target_stone);
	CManager_lock(master->cm);
	master->reconfig = 1;
	master->sig_reconfig_bool = 1;
	check_all_nodes_registered(master);
    }
}

static void
dfg_shutdown_handler(CManager cm, CMConnection conn, void *vmsg, 
		  void *client_data, attr_list attrs)
{
    EVclient client = client_data;
    EVshutdown_ptr msg =  vmsg;
    (void)cm;
    (void)conn;
    (void)attrs;
    int i = 0;
    CManager_lock(cm);
    /* I'm the client, all is done */
    client->shutdown_value = msg->value;
    client->already_shutdown = 1;
    CMtrace_out(cm, EVdfgVerbose, "Client %d has confirmed shutdown\n", client->my_node_id);
    while (client->shutdown_conditions && (client->shutdown_conditions[i] != -1)){
	CMtrace_out(cm, EVdfgVerbose, "Client %d shutdown signalling %d\n", client->my_node_id, client->shutdown_conditions[i]);
	INT_CMCondition_signal(client->cm, client->shutdown_conditions[i++]);
    }
//GSE  - does client have state?
//    CMtrace_out(cm, EVdfgVerbose, "EVDFG exit shutdown master DFG state is %s\n", str_state[client->state]);
    CManager_unlock(cm);
}

static void
handle_shutdown_contrib(EVmaster master, EVmaster_msg_ptr mmsg)
{

    EVshutdown_contribution_ptr msg =  &mmsg->u.shutdown_contrib;
    possibly_signal_shutdown(master, msg->value, mmsg->conn);
    CMtrace_out(master->cm, EVdfgVerbose, "EVDFG exit shutdown master DFG state is %s\n", str_state[master->state]);
}

static void
dfg_stone_close_handler(CManager cm, CMConnection conn, int stone, 
		  void *client_data)
{
    EVclient client = (EVclient)client_data;
    event_path_data evp = cm->evp;
    int global_stone_id = -1;
    CMFormat conn_shutdown_msg = INT_CMlookup_format(client->cm, EVdfg_conn_shutdown_format_list);
    EVconn_shutdown_msg msg;
    (void)cm;
    (void)conn;
    CManager_lock(cm);
    /* first, freeze the stone so that we don't lose any more data */
    INT_EVfreeze_stone(cm, stone);

    int i;
    for (i=0; i < evp->stone_lookup_table_size; i++ ) {
	if (stone == evp->stone_lookup_table[i].local_id) {
	    global_stone_id = evp->stone_lookup_table[i].global_id;
	}
    }
    if (global_stone_id == -1) {
	CMtrace_out(cm, EVdfgVerbose, "Bad mojo, failed to find global stone id after stone close of stone %d\n", stone);
	CMtrace_out(cm, EVdfgVerbose, "  If the above message occurs during shutdown, this is likely not a concern\n");
	CManager_unlock(cm);
	return;
    }
    msg.stone = global_stone_id;
    if (client->master_connection != NULL) {
	INT_CMwrite(client->master_connection, conn_shutdown_msg, &msg);
    } else {
	queue_master_msg(client->master, (void*)&msg, DFGconn_shutdown, NULL, /*copy*/0);
    }
    CManager_unlock(client->cm);
}

extern int
INT_EVmaster_assign_canonical_name(EVmaster master, char *given_name, char *canonical_name)
{
    int node;
    for (node = 0; node < master->node_count; node++) {
	if (master->nodes[node].name == given_name) {
	    if (master->dfg && (master->dfg->realized == 1)) {
		CMtrace_out(master->cm, EVdfgVerbose, "Reconfigure canonical name assignment, node = %d\n", node);
	    } else {
		CMtrace_out(master->cm, EVdfgVerbose, "Canonical name assignment, node = %d, given name was %s, canonical is %s\n", node, given_name, canonical_name);
	    }
	    master->nodes[node].canonical_name = strdup(canonical_name);
	}
    }
    return 1;
}

static void
handle_flush_reconfig(EVmaster master, EVmaster_msg_ptr mmsg)
{
    EVflush_attrs_reconfig_ptr msg = &mmsg->u.flush_reconfig;
    int i, j;
    EVdfg dfg = master->dfg;
    assert(CManager_locked(master->cm));
    if (((EVflush_attrs_reconfig_ptr)msg)->reconfig) {
	master->state = DFG_Reconfiguring;
    }
    CMtrace_out(master->cm, EVdfgVerbose, "EVDFG flush_attr_reconfig -  master DFG state is now %s\n", str_state[master->state]);
    for (i=0; i < msg->count; i++) {
	/* go through incoming attributes */
	for (j=0; j< dfg->stone_count; j++) {
	    if (dfg->deployed_state->stones[j]->stone_id == msg->attr_stone_list[i].stone) {
		if (dfg->deployed_state->stones[j]->attrs != NULL) {
		    free_attr_list(dfg->deployed_state->stones[j]->attrs);
		}
		dfg->deployed_state->stones[j]->attrs = attr_list_from_string(msg->attr_stone_list[i].attr_str);
		if (dfg->working_state->stones[j]->attrs != NULL) {
		    free_attr_list(dfg->working_state->stones[j]->attrs);
		}
		dfg->working_state->stones[j]->attrs = attr_list_from_string(msg->attr_stone_list[i].attr_str);
		break;
	    }
	}
    }
    if (msg->reconfig) {
	CManager_unlock(master->cm);
	master->node_reconfig_handler(master->dfg);
	CManager_lock(master->cm);
	master->reconfig = 1;
	master->sig_reconfig_bool = 1;
	check_all_nodes_registered(master);
    }
}

static void
handle_node_join(EVmaster master, EVmaster_msg_ptr msg)
{
    char *node_name = msg->u.node_join.node_name;
    char *contact_string = msg->u.node_join.contact_string;
    CMConnection conn = msg->conn;
    int node;
    int new_node = -1;

    assert(CManager_locked(master->cm));

    if (master->state == DFG_Running) {
	master->state = DFG_Reconfiguring;
	CMtrace_out(master->cm, EVdfgVerbose, "EVDFG node_join -  master DFG state is now %s\n", str_state[master->state]);
    }

    if (master->node_join_handler == NULL) {
	/* static node list */
	for (node = 0; node < master->node_count; node++) {
	    if (strcmp(master->nodes[node].name, node_name) == 0) {
		if (conn == NULL) {
		    /* we are the master joining as a client node */
		    master->nodes[node].self = 1;
		    master->client->my_node_id = node;
		} else {
		    INT_CMConnection_add_reference(conn);  /* cause we'll be keeping this */
		    master->nodes[node].conn = conn;
		    master->nodes[node].str_contact_list = strdup(contact_string);
		    master->nodes[node].contact_list = attr_list_from_string(master->nodes[node].str_contact_list);
		    master->nodes[node].shutdown_status_contribution = STATUS_UNDETERMINED;
		}
		new_node = node;
		break;
	    }
	}
	if (new_node == -1) {
	    printf("Registering node \"%s\" not found in node list\n", 
		   node_name);
	    return;
	}
    } else {
	int n;
	
	if (master->dfg && master->dfg->realized == 1 && master->reconfig == 0) {
	    master->reconfig = 1;
	    master->sig_reconfig_bool = 1;
	    master->old_node_count = master->node_count;
	    CMtrace_out(master->cm, EVdfgVerbose, "Reconfigure, contact_string = %s\n", contact_string);
	    CMtrace_out(master->cm, EVdfgVerbose, "node_count = %d, stone_count = %d\n", master->node_count, master->dfg->stone_count);
	}
	n = master->node_count++;
	master->nodes = realloc(master->nodes, (sizeof(master->nodes[0])*master->node_count));
	memset(&master->nodes[n], 0, sizeof(master->nodes[0]));
	master->nodes[n].name = strdup(node_name);
	master->nodes[n].canonical_name = NULL;
	master->nodes[n].shutdown_status_contribution = STATUS_UNDETERMINED;
	new_node = n;
	if (conn == NULL) {
	    master->nodes[n].self = 1;
	    master->client->my_node_id = n;
	} else {
	    INT_CMConnection_add_reference(conn);  /* cause we'll be keeping this */
	    master->nodes[n].self = 0;
	    master->nodes[n].conn = conn;
	    master->nodes[n].str_contact_list = strdup(contact_string);
	    master->nodes[n].contact_list = attr_list_from_string(master->nodes[n].str_contact_list);
	}
    }
    CMtrace_out(master->cm, EVdfgVerbose, "Client \"%s\" has joined DFG, contact %s\n", node_name, master->nodes[new_node].str_contact_list);
    check_all_nodes_registered(master);
}


static void
dfg_deploy_handler(CManager cm, CMConnection conn, void *vmsg, 
		  void *client_data, attr_list attrs)
{
    EVclient client = (EVclient) client_data;
    event_path_data evp = cm->evp;
    (void) conn;
    (void) attrs;
    static int first_time_deploy = 1;
    EVdfg_deploy_ptr msg =  vmsg;
    int i, base = evp->stone_lookup_table_size;
    int auto_stones = 0;
    auto_stone_list *auto_list = malloc(sizeof(auto_stone_list));

    CMtrace_out(cm, EVdfgVerbose, "Client %d getting Deploy message\n", client->my_node_id);

    CManager_lock(cm);
    /* add stones to local lookup table */
    if (evp->stone_lookup_table_size == 0) {
	evp->stone_lookup_table = 
	    malloc(sizeof(evp->stone_lookup_table[0]) * msg->stone_count);
    } else {
	evp->stone_lookup_table = 
	    realloc(evp->stone_lookup_table,
		    sizeof(evp->stone_lookup_table[0]) * (msg->stone_count+base));
    }
    for (i=0; i < msg->stone_count; i++) {
	evp->stone_lookup_table[base + i].global_id = msg->stone_list[i].global_stone_id;
	evp->stone_lookup_table[base + i].local_id = INT_EValloc_stone(cm);
    }
    evp->stone_lookup_table_size = base + i;
    for (i=0; i < msg->stone_count; i++) {
	int local_stone = evp->stone_lookup_table[base + i].local_id;
	int local_list[1024]; /* List of output actions for this stone... better be enough */
	int j;
	if (msg->stone_list[i].attrs != NULL) {
	    attr_list tmp_attrs = attr_list_from_string(msg->stone_list[i].attrs);
	    INT_EVset_attr_list(cm, local_stone, tmp_attrs);
	    free_attr_list(tmp_attrs);
	}
	for (j=0; j < msg->stone_list[i].out_count; j++) {
	    if (msg->stone_list[i].out_links[j] != -1) {
		local_list[j] = lookup_local_stone(evp, msg->stone_list[i].out_links[j]);
		if (local_list[j] == -1) {
		    printf("Didn't found global stone %d\n", msg->stone_list[i].out_links[j]);
		}
	    } else {
		local_list[j] = -1;
	    }
	}
	local_list[msg->stone_list[i].out_count] = -1;
	INT_EVassoc_general_action(cm, local_stone, msg->stone_list[i].action, 
				   &local_list[0]);
	for (j=0; j < msg->stone_list[i].extra_actions; j++) {
	    INT_EVassoc_general_action(cm, local_stone, msg->stone_list[i].xactions[j], 
				       &local_list[0]);
	}	    
	if (msg->stone_list[i].period_secs != -1) {
	    auto_list= realloc(auto_list, sizeof(auto_list[0]) * (auto_stones+2));
	    auto_list[auto_stones].stone = local_stone;
	    auto_list[auto_stones].period_secs = msg->stone_list[i].period_secs;
	    auto_list[auto_stones].period_usecs = msg->stone_list[i].period_usecs;
	    auto_stones++;
	}
	if (action_type(msg->stone_list[i].action) == Action_Terminal) {
	    client->active_sink_count++;
	}
    }    
    auto_list[auto_stones].period_secs = -1;
    if (conn != NULL) {
	CMFormat deploy_ack_msg = INT_CMlookup_format(client->cm, EVdfg_deploy_ack_format_list);
	EVdeploy_ack_msg response_msg;
	response_msg.node_id = msg->canonical_name;
	INT_CMwrite(client->master_connection, deploy_ack_msg, &response_msg);
	CMtrace_out(cm, EVdfgVerbose, "Client %d wrote deploy ack\n", client->my_node_id);
    } else {
      	CMtrace_out(cm, EVdfgVerbose, "Client %d no master conn\n", client->my_node_id);
    }
    if (first_time_deploy) {
	first_time_deploy = 0;
    }
    if (auto_stones == 0) {
	free(auto_list);
	auto_list = NULL;
    }
    client->pending_auto_list = auto_list;
    
    CManager_unlock(cm);
}

static void
free_master(CManager cm, void *vmaster)
{
    EVmaster master = (EVmaster)vmaster;
    int i;
    for (i=0; i < master->node_count; i++) {
	if (master->nodes[i].name) free(master->nodes[i].name);
	if (master->nodes[i].canonical_name) 
	    free(master->nodes[i].canonical_name);
	if (master->nodes[i].contact_list) 
	    free_attr_list(master->nodes[i].contact_list);
	if (master->nodes[i].str_contact_list) 
	    free(master->nodes[i].str_contact_list);
    }
    free(master->nodes);
    if (master->my_contact_str) free(master->my_contact_str);
    free(master);
}

static void
free_client(CManager cm, void *vclient)
{
    EVclient client = (EVclient)vclient;
    if (client->master_connection) 
	INT_CMConnection_close(client->master_connection);
    if (client->master_contact_str) free(client->master_contact_str);
    if (client->shutdown_conditions) free(client->shutdown_conditions);
    free(client);
}

static void
free_dfg(CManager cm, void *vdfg)
{
    EVdfg dfg = vdfg;
    int i;
    for (i=0; i < dfg->stone_count; i++) {
	if (dfg->stones[i]) free(dfg->stones[i]);
    }
    if (dfg->stones) free(dfg->stones);
    if (dfg->deployed_state) free_dfg_state(dfg->deployed_state);
    if (dfg->working_state) free_dfg_state(dfg->working_state);
    free(dfg);
}

static EVflush_attrs_reconfig_ptr
build_attrs_msg(EVclient client)
{
    CManager cm = client->cm;
    event_path_data evp = cm->evp;
    int i = 0, cur_stone;
    EVflush_attrs_reconfig_ptr msg = malloc(sizeof(*msg));
    memset(msg, 0, sizeof(*msg));
    msg->attr_stone_list = malloc(sizeof(EVattr_stone_struct));
    for (cur_stone = evp->stone_base_num; cur_stone < evp->stone_count + evp->stone_base_num; ++cur_stone) {
	stone_type stone = stone_struct(evp, cur_stone);
	if (stone->stone_attrs != NULL) {
	    msg->attr_stone_list[i].stone = lookup_global_stone(evp, stone->local_id);
	    msg->attr_stone_list[i].attr_str = attr_list_to_string(stone->stone_attrs);
	    i++;
	    msg->attr_stone_list = realloc(msg->attr_stone_list, sizeof(EVattr_stone_struct)*(i+1));
	}
    }
    msg->count = i;
    return (msg);
}

static void
free_attrs_msg(EVflush_attrs_reconfig_ptr msg)
{
    int i = 0;
    for (i = 0; i < msg->count; i++) {
	free(msg->attr_stone_list[i].attr_str);
    }
    free(msg->attr_stone_list);
    free(msg);
}

static void
flush_and_trigger(EVclient client, int reconfig)
{
    EVflush_attrs_reconfig_ptr msg = build_attrs_msg(client);
    CMFormat flush_msg = INT_CMlookup_format(client->cm, EVdfg_flush_attrs_reconfig_format_list);
    msg->reconfig = reconfig;
    if (client->master_connection != NULL) {
	/* we are a client, send the reconfig to the master */
	INT_CMwrite(client->master_connection, flush_msg, msg);
	free_attrs_msg(msg);
    } else {
	queue_master_msg(client->master, &flush_msg, DFGflush_reconfig, NULL, /*copy*/ 0);
    }
}


static void
cod_EVdfg_trigger_reconfig(cod_exec_context ec)
{
    CManager cm = get_cm_from_ev_state((void*)cod_get_client_data(ec, 0x34567890));
    event_path_data evp = cm->evp;
    EVclient client = evp->app_stone_close_data;  /* cheating a bit.  We know we store the DFG pointer here */
    flush_and_trigger(client, 1);
}

static void
cod_EVdfg_flush_attrs(cod_exec_context ec)
{
    CManager cm = get_cm_from_ev_state((void*)cod_get_client_data(ec, 0x34567890));
    event_path_data evp = cm->evp;
    EVclient client = evp->app_stone_close_data;  /* cheating a bit.  We know we store the DFG pointer here */
    flush_and_trigger(client, 0);
}

extern EVmaster
INT_EVmaster_create(CManager cm)
{
    EVmaster master = malloc(sizeof(struct _EVmaster));
    attr_list contact_list;

    memset(master, 0, sizeof(struct _EVmaster));
    master->cm = cm;
    master->reconfig = 0;
    master->sig_reconfig_bool = 0;
    master->old_node_count = 1;
    master->no_deployment = 0;
    master->state = DFG_Joining;

    CMtrace_out(cm, EVdfgVerbose, "EVDFG initialization -  master DFG state set to %s\n", str_state[master->state]);
    contact_list = INT_CMget_contact_list(cm);
    master->my_contact_str = attr_list_to_string(contact_list);
    free_attr_list(contact_list);

    /*
     * EVdfg master-sent messages
     */
    INT_CMregister_format(cm, EVdfg_ready_format_list);
    INT_CMregister_format(cm, EVdfg_deploy_format_list);
    INT_CMregister_format(cm, EVclient_shutdown_format_list);

    /*
     * EVdfg master-handled messages
     */
    INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_node_join_format_list),
			   dfg_master_msg_handler, (void*)(((uintptr_t)master)|DFGnode_join));
    INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_deploy_ack_format_list),
			   dfg_master_msg_handler, (void*)(((uintptr_t)master)|DFGdeploy_ack));
    INT_CMregister_handler(INT_CMregister_format(cm, EVclient_shutdown_contribution_format_list),
			   dfg_master_msg_handler, (void*)(((uintptr_t)master)|DFGshutdown_contrib));
    INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_conn_shutdown_format_list),
			   dfg_master_msg_handler, (void*)(((uintptr_t)master)|DFGconn_shutdown));
    INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_flush_attrs_reconfig_format_list),
			   dfg_master_msg_handler, (void*)(((uintptr_t)master)|DFGflush_reconfig));

    INT_CMadd_poll(cm, handle_queued_messages_lock, master);
    INT_CMadd_shutdown_task(cm, free_master, master, FREE_TASK);
    return master;
}

static EVdfg_configuration
new_dfg_configuration(EVdfg master)
{
    EVdfg_configuration ret = malloc(sizeof(*ret));
    memset(ret, 0, sizeof(*ret));
    ret->stone_count = 0;
    ret->stones = NULL;
    return ret;
}

extern EVdfg
INT_EVdfg_create(EVmaster master)
{
    EVdfg dfg = malloc(sizeof(struct _EVdfg));

    memset(dfg, 0, sizeof(struct _EVdfg));
    dfg->master = master;
    dfg->deployed_stone_count = 0;
    dfg->deploy_ack_condition = -1;
    master->dfg = dfg;
    if (master->client) {
	master->client->dfg = dfg;
	dfg->client = master->client;
    }
    master->reconfig = 0;
    master->sig_reconfig_bool = 0;
    master->old_node_count = 1;
    master->state = DFG_Joining;
    CMtrace_out(master->cm, EVdfgVerbose, "EVDFG initialization -  master DFG state set to %s\n", str_state[master->state]);

    dfg->working_state = new_dfg_configuration(dfg);
    dfg->stones = malloc(sizeof(dfg->stones[0]));
    INT_CMadd_shutdown_task(master->cm, free_dfg, dfg, FREE_TASK);
    return dfg;
}


extern char *INT_EVmaster_get_contact_list(EVmaster master)
{
    attr_list contact_list = NULL;
    CManager cm = master->cm;
    char *tmp = NULL;

    /* use enet transport if available */
#if defined(ENET_FOUND) || defined(ZPL_ENET_AVAILABLE)
    atom_t CM_ENET_CONN_TIMEOUT = attr_atom_from_string("CM_ENET_CONN_TIMEOUT");
    atom_t CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    (void) CM_TRANSPORT;
    attr_list listen_list = create_attr_list();
#if defined(ENET_FOUND)
    add_string_attr(listen_list, CM_TRANSPORT, strdup("enet"));
#elif defined(ZPL_ENET_AVAILABLE)
    add_string_attr(listen_list, CM_TRANSPORT, strdup("zplenet"));
#endif
    /* and kick up the connection timeout value.  We can wait 60 secs */
    contact_list = INT_CMget_specific_contact_list(cm, listen_list);
    if(contact_list) {
        add_int_attr(contact_list, CM_ENET_CONN_TIMEOUT, 60000);
    }
    free_attr_list(listen_list);
#endif

    if (contact_list == NULL) {
	contact_list = INT_CMget_contact_list(cm);
	if (contact_list == NULL) {
	    CMlisten(cm);
	    contact_list = INT_CMget_contact_list(cm);
	}
    }
    if(contact_list) {
        tmp = attr_list_to_string(contact_list);
    }
    free_attr_list(contact_list);
    return tmp;
}

static int
max_output_for_action(char *action, int cur_max)
{
    if (cur_max == -1) return -1;
    switch (action_type(action)) {
    case Action_NoAction:
    case Action_Terminal:
    case Action_Bridge:
	return cur_max;
    case Action_Filter:
    case Action_Immediate:
	if (strncmp(action, "Router Action", 13) == 0) {
	    /* router spec, could be lots */
	    return -1;
	}
	if (cur_max < 1) {
	    return 1;
	} else {
	    return cur_max;
	}
    case Action_Multi:
    case Action_Split:
    case Action_Source:
	return -1; /* no max */
    default:
	printf("Didn't expect case in max_output_for_action\n");
	exit(1);
    }
}

static void
check_connectivity(EVdfg_configuration state, EVmaster master)
{
    int i;
    for (i=0; i< state->stone_count; i++) {
	int j;
	int max_output = 0;
	CMtrace_out(master->cm, EVdfgVerbose, "Stone %d - assigned to node %s, action %s\n", i, 
		    master->nodes[state->stones[i]->node].canonical_name, (state->stones[i]->action ? state->stones[i]->action : "NULL"));
	if (state->stones[i]->node == -1) {
	    printf("Warning, stone %d has not been assigned to any node.  This stone will not be deployed.\n", i);
	    printf("    This stones particulars are:\n");
	    fdump_dfg_stone(stdout, state->stones[i]);
	}
	if (state->stones[i]->bridge_stone) {
	    /* we created the bridge stones, assume they're OK */
	    continue;
	}
	if (state->stones[i]->action_count == 0) {
	    printf("Warning, stone %d (assigned to node %s) has no actions registered", i, master->nodes[state->stones[i]->node].canonical_name);
	    continue;
	}
	max_output = max_output_for_action(state->stones[i]->action, max_output);
	for (j=0; j< state->stones[i]->action_count - 1; j++) {
	    max_output = max_output_for_action(state->stones[i]->extra_actions[j], max_output);
	}
	if ((state->stones[i]->out_count == 0) && (max_output > 0)) {
	    printf("Warning, stone %d (assigned to node %s) has no outputs connected to other stones\n", i, master->nodes[state->stones[i]->node].canonical_name);
	    printf("    This stones particulars are:\n");
	    fdump_dfg_stone(stdout, state->stones[i]);
	    continue;
	}
	if ((max_output == 1) && (state->stones[i]->out_count > 1)) {
	    printf("Warning, stone %d (assigned to node %s) has more than one output port linked, but can only support one output\n", i, master->nodes[state->stones[i]->node].canonical_name);
	    printf("    This stones particulars are:\n");
	    fdump_dfg_stone(stdout, state->stones[i]);
	    continue;
	}
	if ((max_output == 1) && (state->stones[i]->out_links[0] == -1)) {
	    printf("Warning, stone %d (assigned to node %s) produces at least one output, but output port 0 is unlinked\n", i, master->nodes[state->stones[i]->node].canonical_name);
	    printf("    This stones particulars are:\n");
	    fdump_dfg_stone(stdout, state->stones[i]);
	}
    }
}

extern int
INT_EVdfg_realize(EVdfg dfg)
{
    check_connectivity(dfg->working_state, dfg->master);
//    check_types(dfg);

    if (dfg->realized == 1) {
	dfg->master->reconfig = 0;
    }
    dfg->realized = 1;
    return 1;
}

extern void
INT_EVmaster_register_node_list(EVmaster master, char **nodes)
{
    int count = 0, i = 0;
    while(nodes[count] != NULL) count++;
    master->node_count = count;
    master->nodes = malloc(sizeof(master->nodes[0]) * count);
    memset(master->nodes, 0, sizeof(master->nodes[0]) * count);
    for (i = 0; i < master->node_count; i++) {
	master->nodes[i].name = strdup(nodes[i]);
	master->nodes[i].canonical_name = strdup(nodes[i]);
	master->nodes[i].shutdown_status_contribution = STATUS_UNDETERMINED;
    }
}

extern int
INT_EVdfg_assign_node(EVdfg_stone stone, char *node_name)
{
    EVdfg dfg = stone->dfg;
    EVmaster master = dfg->master;
    int i, node = -1;
    for (i = 0; i < master->node_count; i++) {
	EVint_node_list n = &master->nodes[i];
	if (n->canonical_name && (strcmp(n->canonical_name, node_name) == 0)) {
	    node = i;
	} else 	if (n->name && (strcmp(n->name, node_name) == 0)) {
	    node = i;
	}

    }
    if (node == -1) {
	printf("Node \"%s\" not found in node list\n", node_name);
	return 0;
    }
	
    if (dfg->realized == 1) {
	CMtrace_out(master->cm, EVdfgVerbose, "assign node, node# = %d\n", node);
    }
    EVdfg_config_action act;
    act.type = ACT_assign_node;
    act.stone_id = stone->stone_id;
    act.u.assign.dest_node = node;
    (void) EVdfg_perform_act_on_state(stone->dfg->working_state, act, 1 /* add to queue */);
    return 1;
}

extern int 
INT_EVclient_ready_wait(EVclient client)
{
    CMtrace_out(client->cm, EVdfgVerbose, "DFG %p wait for ready\n", client);
    INT_CMCondition_wait(client->cm, client->ready_condition);
    client->ready_condition = -1;
    CMtrace_out(client->cm, EVdfgVerbose, "DFG %p ready wait released\n", client);
    return 1;

}

extern int
INT_EVclient_shutdown(EVclient client, int result)
{
    CMFormat shutdown_msg = INT_CMlookup_format(client->cm, EVclient_shutdown_contribution_format_list);
    EVshutdown_contribution_msg msg;
    if (client->already_shutdown) printf("Node %d, already shut down BAD!\n", client->my_node_id);
    msg.value = result;
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d calling client_shutdown\n", client->my_node_id);
    if (client->master_connection == NULL) {
	queue_master_msg(client->master, (void*)&msg, DFGshutdown_contrib, NULL, /*copy*/0);
    } else {
	/* we are a client, tell the master to shutdown */
	INT_CMwrite(client->master_connection, shutdown_msg, &msg);
    }
    if (!client->already_shutdown) {
	CManager_unlock(client->cm);
	CMtrace_out(client->cm, EVdfgVerbose, "Client %d shutdown condition wait\n", client->my_node_id);
	CMCondition_wait(client->cm, new_shutdown_condition(client, client->master_connection));
	CMtrace_out(client->cm, EVdfgVerbose, "Client %d shutdown condition wait DONE!\n", client->my_node_id);
	CManager_lock(client->cm);
    }
    return client->shutdown_value;
}

extern int
INT_EVclient_force_shutdown(EVclient client, int result)
{
    result |= STATUS_FORCE;
    if (client->already_shutdown) printf("Node %d, already contributed to shutdown.  Don't call shutdown twice!\n", client->my_node_id);
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d calling client_FORCE_shutdown\n", client->my_node_id);
    if (client->master_connection != NULL) {
	/* we are a client, tell the master to shutdown */
	CMFormat shutdown_msg = INT_CMlookup_format(client->cm, EVclient_shutdown_contribution_format_list);
	EVshutdown_contribution_msg msg;
	msg.value = result;
	INT_CMwrite(client->master_connection, shutdown_msg, &msg);
	/* and wait until we hear back */
    } else {
	possibly_signal_shutdown(client->master, result, NULL);
    }
    if (!client->already_shutdown) {
	CManager_unlock(client->cm);
	CMtrace_out(client->cm, EVdfgVerbose, "Client %d shutdown condition wait\n", client->my_node_id);
	CMCondition_wait(client->cm, new_shutdown_condition(client, client->master_connection));
	CMtrace_out(client->cm, EVdfgVerbose, "Client %d shutdown condition wait DONE!\n", client->my_node_id);
	CManager_lock(client->cm);
    }
    return client->shutdown_value;
}

extern int
INT_EVclient_active_sink_count(EVclient client)
{
    return client->active_sink_count;
}

extern void
INT_EVclient_ready_for_shutdown(EVclient client)
{
    if (client->already_shutdown) return;
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d ready for shutdown \n", client->my_node_id);
    if (client->master_connection != NULL) {
	/* we are a client, tell the master to shutdown */
	CMFormat shutdown_msg = INT_CMlookup_format(client->cm, EVclient_shutdown_contribution_format_list);
	EVshutdown_contribution_msg msg;
	msg.value = STATUS_NO_CONTRIBUTION;   /* no status contribution */
	INT_CMwrite(client->master_connection, shutdown_msg, &msg);
    } else {
	possibly_signal_shutdown(client->master, STATUS_NO_CONTRIBUTION, NULL);
    }
}

extern int 
INT_EVclient_wait_for_shutdown(EVclient client)
{
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d wait for shutdown \n", client->my_node_id);
    if (client->already_shutdown) return client->shutdown_value;
    INT_CMCondition_wait(client->cm, new_shutdown_condition(client, client->master_connection));
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d wait for shutdown DONE! \n", client->my_node_id);
    return client->shutdown_value;
}

extern int 
INT_EVclient_test_for_shutdown(EVclient client)
{
    CMtrace_out(client->cm, EVdfgVerbose, "Client %d testing for shutdown return %d\n", client->my_node_id,
	client->already_shutdown);
    return client->already_shutdown;
}

extern int INT_EVclient_source_active(EVsource src)
{
    return (src->local_stone_id != -1);
}

extern EVclient_sources
INT_EVclient_register_source(char *name, EVsource src)
{
    CManager cm = src->cm;
    event_path_data evp = cm->evp;
    if (evp->source_count == 0) {
	evp->sources = malloc(sizeof(evp->sources[0]));
    } else {
	evp->sources = realloc(evp->sources,
			       sizeof(evp->sources[0]) * (evp->source_count + 1));
    }
    evp->sources[evp->source_count].name = strdup(name);
    evp->sources[evp->source_count].src = src;
    evp->source_count++;
    return evp->sources;
}

extern EVclient_sinks
INT_EVclient_register_sink_handler(CManager cm, char *name, FMStructDescList list, EVSimpleHandlerFunc handler, void* client_data)
{
    event_path_data evp = cm->evp;
    if (evp->sink_handler_count == 0) {
	evp->sink_handlers = malloc(sizeof(evp->sink_handlers[0]));
    } else {
	evp->sink_handlers = realloc(evp->sink_handlers,
				     sizeof(evp->sink_handlers[0]) * (evp->sink_handler_count + 1));
    }
    evp->sink_handlers[evp->sink_handler_count].name = strdup(name);
    evp->sink_handlers[evp->sink_handler_count].format_list = list;
    evp->sink_handlers[evp->sink_handler_count].handler = handler;
    evp->sink_handlers[evp->sink_handler_count].client_data = client_data;
    evp->sink_handler_count++;
    return evp->sink_handlers;
}

extern EVclient_sinks
INT_EVclient_register_raw_sink_handler(CManager cm, char *name, EVRawHandlerFunc handler, void *client_data)
{
    event_path_data evp = cm->evp;
    if (evp->sink_handler_count == 0) {
	evp->sink_handlers = malloc(sizeof(evp->sink_handlers[0]));
    } else {
	evp->sink_handlers = realloc(evp->sink_handlers,
				     sizeof(evp->sink_handlers[0]) * (evp->sink_handler_count + 1));
    }
    evp->sink_handlers[evp->sink_handler_count].name = strdup(name);
    evp->sink_handlers[evp->sink_handler_count].format_list = NULL;
    evp->sink_handlers[evp->sink_handler_count].handler = (EVSimpleHandlerFunc)handler;
    evp->sink_handlers[evp->sink_handler_count].client_data = client_data;
    evp->sink_handler_count++;
    return evp->sink_handlers;
}

static int
new_shutdown_condition(EVclient client, CMConnection conn)
{
    int cur_count = 0;
    if (client->shutdown_conditions == NULL) {
	client->shutdown_conditions = malloc(2*sizeof(client->shutdown_conditions[0]));
    } else {
	while (client->shutdown_conditions[cur_count++] != -1) ; 
	cur_count--;
	client->shutdown_conditions = realloc(client->shutdown_conditions, 
					   (cur_count+2)*sizeof(client->shutdown_conditions[0]));
    }
    client->shutdown_conditions[cur_count] = INT_CMCondition_get(client->cm, conn);
    client->shutdown_conditions[cur_count+1] = -1;
    return client->shutdown_conditions[cur_count];
}

static char dfg_extern_string[] = "\
	void EVdfg_trigger_reconfiguration(cod_exec_context ec);\n\
	void EVdfg_flush_attrs(cod_exec_context ec);\n\
";

static cod_extern_entry dfg_extern_map[] = {
    {"EVdfg_trigger_reconfiguration", (void *) 0},
    {"EVdfg_flush_attrs", (void *) 0},
    {(void*)0, (void*)0}
};

extern EVclient
dfg_assoc_client(CManager cm, char* node_name, char *master_contact, EVmaster master,
		 EVclient_sources source_capabilities, EVclient_sinks sink_capabilities)
{
    event_path_data evp = cm->evp;
    attr_list master_attrs = NULL;
    CMConnection conn;
    CMFormat register_msg = NULL;
    EVnode_join_msg msg;
    attr_list contact_list = INT_CMget_contact_list(cm);
    char *my_contact_str;
    EVclient client;
    int i;

    register_msg = INT_CMlookup_format(cm, EVdfg_ready_format_list);
    if ((master && master->client) ||
	(!master && register_msg)) {
	fprintf(stderr, "Rejecting attempt to associate a DFG client with another DFG or with the same DFG multiple tiems.\n");
	fprintf(stderr, "Only one call to EVclient_assoc() or EVclient_assoc_local() per CManager allowed.\n");
	return NULL;
    }
    dfg_extern_map[0].extern_value = (void*)(intptr_t)cod_EVdfg_trigger_reconfig;
    dfg_extern_map[1].extern_value = (void*)(intptr_t)cod_EVdfg_flush_attrs;

    INT_EVadd_standard_routines(cm, dfg_extern_string, dfg_extern_map);

    client = malloc(sizeof(*client));
    memset(client, 0, sizeof(*client));
    client->cm = cm;
    client->pending_auto_list = NULL;
    if (master_contact) {
	master_attrs = attr_list_from_string(master_contact);
	client->master_contact_str = strdup(master_contact);
    } else {
	client->master = master;
	client->dfg = master->dfg;
	if (master->dfg)
	    master->dfg->client = client;
	master->client = client;
    }

    client->ready_condition = INT_CMCondition_get(cm, NULL);

    if (contact_list == NULL) {
	INT_CMlisten(cm);
	contact_list = INT_CMget_contact_list(cm);
    }

    my_contact_str = attr_list_to_string(contact_list);
    free_attr_list(contact_list);
    
    msg.node_name = strdup(node_name);
    msg.contact_string = my_contact_str;
    msg.source_count = evp->source_count;
    msg.sources = malloc(msg.source_count * sizeof(msg.sources[0]));
    for (i=0; i < evp->source_count; i++) {
	msg.sources[i].name = strdup(evp->sources[i].name);
	msg.sources[i].FMtype = NULL;
    }
    msg.sink_count = evp->sink_handler_count;
    msg.sinks = malloc(msg.sink_count * sizeof(msg.sinks[0]));
    for (i=0; i < evp->sink_handler_count; i++) {
	msg.sinks[i].name = strdup(evp->sink_handlers[i].name);
	msg.sinks[i].FMtype = NULL;
    }
    INT_EVregister_close_handler(cm, dfg_stone_close_handler, (void*)client);
    
    if (master) {
	/* local master */
	queue_master_msg(master, (void*)&msg, DFGnode_join, NULL, /*copy*/0);
    } else {
	/*
	 * EVdfg client-sent messages
	 */
	register_msg = INT_CMregister_format(cm, EVdfg_node_join_format_list);
	INT_CMregister_format(cm, EVdfg_deploy_ack_format_list);
	INT_CMregister_format(cm, EVclient_shutdown_contribution_format_list);
	INT_CMregister_format(cm, EVdfg_conn_shutdown_format_list);
	INT_CMregister_format(cm, EVdfg_flush_attrs_reconfig_format_list);

	/*
	 * EVdfg client-handled messages
	 */
	INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_ready_format_list),
			       dfg_ready_handler, client);
	INT_CMregister_handler(INT_CMregister_format(cm, EVdfg_deploy_format_list),
			       dfg_deploy_handler, client);
	INT_CMregister_handler(INT_CMregister_format(cm, EVclient_shutdown_format_list),
			       dfg_shutdown_handler, client);

	conn = INT_CMget_conn(cm, master_attrs);
	if (conn == NULL) {
	    fprintf(stderr, "failed to contact Master at %s\n", attr_list_to_string(master_attrs));
	    fprintf(stderr, "Join DFG failed\n");
	    return NULL;
	}
	INT_CMwrite(conn, register_msg, &msg);
	client->master_connection = conn;
	for (i=0; i < evp->source_count; i++) {
	    free(msg.sources[i].name);
	}
	free(msg.sources);
	for (i=0; i < evp->sink_handler_count; i++) {
	    free(msg.sinks[i].name);
	}
	free(msg.sinks);
	free(msg.contact_string);
	free(msg.node_name);
    }
    CMtrace_out(cm, EVdfgVerbose, "DFG %p node name %s\n", client, node_name);
    if (master_attrs) free_attr_list(master_attrs);
    INT_CMadd_shutdown_task(cm, free_client, client, FREE_TASK);
    return client;
}

extern EVclient
INT_EVclient_assoc_local(CManager cm, char* node_name, EVmaster master,
			     EVclient_sources source_capabilities, EVclient_sinks sink_capabilities)
{
    return dfg_assoc_client(cm, node_name, NULL, master, source_capabilities, sink_capabilities);
}

extern EVclient
INT_EVclient_assoc(CManager cm, char* node_name, char *master_contact_str,
    EVclient_sources source_capabilities, EVclient_sinks sink_capabilities)
{
    return dfg_assoc_client(cm, node_name, master_contact_str, NULL, source_capabilities, sink_capabilities);
}

static int
create_bridge_stone(EVdfg dfg, EVdfg_stone_state target, int node)
{
    EVdfg_config_action act, link_act;

    char *contact = dfg->master->nodes[target->node].str_contact_list;
    char *action;
    if (dfg->master->nodes[target->node].self) {
	contact = dfg->master->my_contact_str;
    }
    action = INT_create_bridge_action_spec(target->stone_id, contact);
    act.type = ACT_create_bridge;
    act.u.bridge.action = action;
    act.u.bridge.target_id = target->stone_id;
    act.stone_id = 0x80000000 | dfg->stone_count++;
    dfg->stones = realloc(dfg->stones, sizeof(dfg->stones[0]) * dfg->stone_count);
    dfg->stones[dfg->stone_count-1] = NULL;  /* not a user-visible stone */
    act.node_for_action = node;
    EVdfg_perform_act_on_state(dfg->working_state, act, 1 /* add to queue */);
    link_act.type = ACT_link_dest;
    link_act.stone_id = act.stone_id;
    link_act.node_for_action = node;
    link_act.u.link.dest_id = target->stone_id;
    EVdfg_perform_act_on_state(dfg->working_state, link_act, 0 /* don't add to queue */);
    return act.stone_id;
}

static void
add_bridge_stones(EVdfg dfg, EVdfg_configuration config)
{
    int i;
    for (i=0; i< config->stone_count; i++) {
	int j;
	for (j = 0; j < config->stones[i]->out_count; j++) {
	    EVdfg_stone_state cur = config->stones[i];
	    EVdfg_stone_state target = NULL;
	    int k;
	    for (k = 0; k < config->stone_count; k++) {
		if (config->stones[k]->stone_id == cur->out_links[j]) {
		    target = config->stones[k];
		}
	    }
	    if (target && (!cur->bridge_stone) && (cur->node != target->node)) {
		EVdfg_stone_state bridge = NULL;
		cur->out_links[j] = create_bridge_stone(dfg, target, cur->node);
		/* put the bridge stone where the source stone is */
		for (k = 0; k < config->stone_count; k++) {
		    if (config->stones[k]->stone_id == cur->out_links[j]) {
			bridge = config->stones[k];
		    }
		}
		for (k = 0; k < config->pending_action_count; k++) {
		    EVdfg_config_action act = config->pending_action_queue[k];
		    if (((act.type == ACT_link_port) || (act.type == ACT_link_dest)) &&
			(act.u.link.dest_id == target->stone_id)) {
			config->pending_action_queue[k].u.link.dest_id = bridge->stone_id;
			break;
		    }
		}
		bridge->node = cur->node;
		CMtrace_out(dfg->master->cm, EVdfgVerbose, "Created bridge stone %x, target node is %d, assigned to node %d\n", cur->out_links[j], target->node, cur->node);
	    }
	}
    }
}

static void
add_unfreeze_actions(EVdfg dfg, EVdfg_configuration config)
{
    int i;
    for (i=0; i< config->stone_count; i++) {
	EVdfg_stone_state cur = config->stones[i];
	if (cur->condition == EVstone_Frozen) {
	    EVdfg_config_action fact;
	    fact.type = ACT_unfreeze;
	    fact.stone_id = cur->stone_id;
	    fact.node_for_action = cur->node;
	    EVdfg_add_act_to_queue(config, fact);
	    cur->condition = EVstone_Deployed;
	}
    }
}

static void
add_stone_to_deploy_msg(EVdfg_configuration config, EVdfg_deploy_msg *msg, EVdfg_stone_state dstone) 	    
{
    deploy_msg_stone mstone;
    int k;
    msg->stone_list = realloc(msg->stone_list, (msg->stone_count +1) * sizeof(msg->stone_list[0]));
    memset(&msg->stone_list[msg->stone_count], 0, sizeof(msg->stone_list[0]));
    mstone = &msg->stone_list[msg->stone_count];
    mstone->global_stone_id = dstone->stone_id;
    mstone->attrs = NULL;
    if (dstone->attrs != NULL) {
	mstone->attrs = attr_list_to_string(dstone->attrs);
    }
    mstone->period_secs = dstone->period_secs;
    mstone->period_usecs = dstone->period_usecs;
    if (dstone->bridge_stone) {
	/* bridge stone virtual links should not be deployed, just for bookkeeping */
	mstone->out_count = 0;
	mstone->out_links = NULL;
    } else {
	mstone->out_count = dstone->out_count;
	mstone->out_links = malloc(sizeof(mstone->out_links[0])*mstone->out_count);
	for (k=0; k< dstone->out_count; k++) {
	    if (dstone->out_links[k] != -1) {
		mstone->out_links[k] = dstone->out_links[k];
	    } else {
		mstone->out_links[k] = -1;
	    }
	}
    }
    mstone->action = dstone->action;
    if (dstone->action_count > 1) {
	mstone->extra_actions = dstone->action_count - 1;
	mstone->xactions = malloc(sizeof(mstone->xactions[0])*mstone->extra_actions);
	for (k=0; k < mstone->extra_actions; k++) {
	    mstone->xactions[k] = dstone->extra_actions[k];
	}
    } else {
	mstone->extra_actions = 0;
	mstone->xactions = NULL;
    }
    msg->stone_count++;
}

static void
build_deploy_msg_for_node_stones(EVdfg_configuration config, int act_num, EVmaster master)
{
    EVdfg_config_action orig_act = config->pending_action_queue[act_num];
    int i;
    EVdfg_deploy_msg msg;
    int node = orig_act.node_for_action;
    CMFormat deploy_msg = INT_CMlookup_format(master->cm, EVdfg_deploy_format_list);

    CMtrace_out(master->cm, EVdfgVerbose, "Master in deploy_msg_for_node for client %s, node %d\n",
		master->nodes[node].canonical_name, node);
    memset(&msg, 0, sizeof(msg));
    msg.canonical_name = master->nodes[node].canonical_name;
    msg.stone_count = 0;
    msg.stone_list = malloc(sizeof(msg.stone_list[0]));
    for (i=act_num; i< config->pending_action_count; i++) {
	EVdfg_config_action act = config->pending_action_queue[i];
	switch (act.type) {
	case ACT_no_op: break;
	case ACT_create_bridge:
	case ACT_create:
	    /* build and send deploy msg for this stone and all stones created on it's node */
	    /* kill all actions related to this stone (add_action, link_port, set_attrs, assign_node) */
	    if (act.node_for_action == node) {
		int j;
		for (j=0; j < config->stone_count; j++) {
		    if (config->stones[j]->stone_id == act.stone_id) {
			add_stone_to_deploy_msg(config, &msg, config->stones[j]);
		    }
		}
		for (j=i; j < config->pending_action_count; j++) {
		    if (config->pending_action_queue[j].stone_id == act.stone_id) {
			config->pending_action_queue[j].type = ACT_no_op;
		    }
		}
	    }
	    break;
	case ACT_add_action:
	case ACT_set_auto_period:
	case ACT_link_port:
	case ACT_link_dest:
	case ACT_unlink_port:
	case ACT_unlink_dest:
	case ACT_set_attrs:
	case ACT_assign_node:
	case ACT_destroy:
	case ACT_freeze:
	case ACT_unfreeze:
	    break;
	default:
	    printf("Bad action type in build_deploy_msg_for_nodes (action type %d)\n", act.type);
	    break;
	}
    }
    if (!master->nodes[node].needs_ready) {
	master->nodes[node].needs_ready = 1;
	master->dfg->deploy_ack_count--;
    }
    if (master->nodes[node].conn) {
	INT_CMwrite(master->nodes[node].conn, deploy_msg, &msg);
    } else {
	EVmaster_msg mmsg;
	CManager_unlock(master->cm);
	mmsg.u.deploy_ack.node_id = "master";
	dfg_deploy_handler(master->cm, NULL, &msg, master->client, NULL);
	handle_deploy_ack(master, &mmsg);
	CManager_lock(master->cm);
    }
    for(i=0 ; i < msg.stone_count; i++) {
	free(msg.stone_list[i].out_links);
	if (msg.stone_list[i].attrs) free(msg.stone_list[i].attrs);
	if (msg.stone_list[i].xactions) free(msg.stone_list[i].xactions);
    }
    free(msg.stone_list);
}

static void
deploy_to_node(EVdfg dfg, int node, EVdfg_configuration config)
{
    int i;
    int stone_count = 0;
    EVdfg_deploy_msg msg;
    CMFormat deploy_msg = INT_CMlookup_format(dfg->master->cm, EVdfg_deploy_format_list);

    for (i=dfg->deployed_stone_count; i< dfg->stone_count; i++) {
	if (config->stones[i]->node == node) {
	    stone_count++;
	}
    }
    CMtrace_out(dfg->master->cm, EVdfgVerbose, "Master in deploy_to_node for client %s, node %d, stones to deploy %d\n",
		dfg->master->nodes[node].canonical_name, node, stone_count);
    if (stone_count == 0) {
        dfg->deploy_ack_count++;
	dfg->master->nodes[node].needs_ready = 0;
      	return;
    }
    memset(&msg, 0, sizeof(msg));
    msg.canonical_name = dfg->master->nodes[node].canonical_name;
    msg.stone_count = 0;
    msg.stone_list = malloc(sizeof(msg.stone_list[0]));
    for (i=dfg->deployed_stone_count; i< dfg->stone_count; i++) {
	if (config->stones[i]->node == node) {
	    EVdfg_stone_state dstone = config->stones[i];
	    add_stone_to_deploy_msg(config, &msg, dstone);
	    dstone->condition = EVstone_Deployed;
	}
    }
    dfg->master->nodes[node].needs_ready = 1;
    if (dfg->master->nodes[node].conn) {
	INT_CMwrite(dfg->master->nodes[node].conn, deploy_msg, &msg);
    } else {
	EVmaster_msg mmsg;
	CManager_unlock(dfg->master->cm);
	dfg_deploy_handler(dfg->master->cm, NULL, &msg, dfg->master->client, NULL);
	mmsg.u.deploy_ack.node_id = "master";
	handle_deploy_ack(dfg->master, &mmsg);
	CManager_lock(dfg->master->cm);
    }
    for(i=0 ; i < msg.stone_count; i++) {
	free(msg.stone_list[i].out_links);
	if (msg.stone_list[i].attrs) free(msg.stone_list[i].attrs);
	if (msg.stone_list[i].xactions) free(msg.stone_list[i].xactions);
    }
    free(msg.stone_list);
}

static void
dump_dfg_stone(EVdfg_stone_state s)
{
    fdump_dfg_stone(stdout, s);
}

static void
fdump_dfg_stone(FILE* out, EVdfg_stone_state s)
{
    int i;

    (void)dump_dfg_stone;   /* stop warning aboud dump_dfg_stone, handy to keep around for debugging */

    fprintf(out, "stone %p, node %d, stone_id %x  (current condition %s)\n", s, s->node, s->stone_id,
	    stone_condition_str[s->condition]);
    if (s->bridge_stone) fprintf(out, "      bridge_stone\n");
    fprintf(out, " out_count %d : ", s->out_count);
    for (i=0; i < s->out_count; i++) {
	fprintf(out, "%x, ", s->out_links[i]);
    }
    fprintf(out, "\n action_count %d, action = \"%s\"\n", s->action_count, (s->action ? s->action : "NULL"));
    fprintf(out, "\nbridge_target %x\n", s->bridge_target);
}

static void
free_master_msg(EVmaster_msg *msg)
{
    switch(msg->msg_type) {
    case DFGnode_join: {
	EVnode_join_ptr in = &msg->u.node_join;
	int i;
	free(in->node_name);
	free(in->contact_string);
	for (i=0; i < in->sink_count; i++) {
	    leaf_element *l = &in->sinks[i];
	    if(l->name) free(l->name);
	    if(l->FMtype) free(l->FMtype);
	}
	free(in->sinks);
	for (i=0; i < in->source_count; i++) {
	    leaf_element *l = &in->sources[i];
	    if (l->name) free(l->name);
	    if (l->FMtype) free(l->FMtype);
	}
	free(in->sources);
	break;
    }
    case DFGflush_reconfig: {
	EVflush_attrs_reconfig_ptr in = &msg->u.flush_reconfig;
	int i;
	for (i=0 ; i < in->count; i++) {
	    free(in->attr_stone_list[i].attr_str);
	}
	free(in->attr_stone_list);
	break;
    }
    case DFGdeploy_ack:
    case DFGshutdown_contrib:
    case DFGconn_shutdown:
    default:
	break;
    }
    free(msg);
}

static void
queue_master_msg(EVmaster master, void*vmsg, EVmaster_msg_type msg_type, CMConnection conn, int copy)
{
    EVmaster_msg_ptr msg = malloc(sizeof(EVmaster_msg));
    msg->msg_type = msg_type;
    msg->conn = conn;
    switch(msg_type) {
    case DFGnode_join: {
	EVnode_join_ptr in = (EVnode_join_ptr)vmsg;
	if (!copy) {
	    msg->u.node_join = *in;
	} else {
	    int i;
	    msg->u.node_join.node_name = strdup(in->node_name);
	    msg->u.node_join.contact_string = strdup(in->contact_string);
	    msg->u.node_join.source_count = in->source_count;
	    msg->u.node_join.sink_count = in->sink_count;
	    msg->u.node_join.sinks = (leaf_element*)malloc(sizeof(leaf_element) * in->sink_count);
	    for (i=0; i < in->sink_count; i++) {
		leaf_element *l = &in->sinks[i];
		msg->u.node_join.sinks[i].name = l->name ? strdup(l->name) : NULL;
		msg->u.node_join.sinks[i].FMtype = l->FMtype ? strdup(l->FMtype) : NULL;
	    }
	    msg->u.node_join.sources = (leaf_element*)malloc(sizeof(leaf_element) * in->source_count);
	    for (i=0; i < in->source_count; i++) {
		leaf_element *l = &in->sources[i];
		msg->u.node_join.sources[i].name = l->name ? strdup(l->name) : NULL;
		msg->u.node_join.sources[i].FMtype = l->FMtype ? strdup(l->FMtype) : NULL;
	    }
	}
	break;
    }
    case DFGdeploy_ack: {
	EVdeploy_ack_ptr in = (EVdeploy_ack_ptr)vmsg;
	msg->u.deploy_ack = *in;
	break;
    }
    case DFGshutdown_contrib: {
	EVshutdown_contribution_ptr in = (EVshutdown_contribution_ptr)vmsg;
	msg->u.shutdown_contrib = *in;
	break;
    }
    case  DFGconn_shutdown: {
	EVconn_shutdown_ptr in = (EVconn_shutdown_ptr)vmsg;
	msg->u.conn_shutdown = *in;
	break;
    }
    case DFGflush_reconfig: {
	EVflush_attrs_reconfig_ptr in = (EVflush_attrs_reconfig_ptr)vmsg;
	msg->u.flush_reconfig = *in;
	if (copy) {
	    int i;
	    msg->u.flush_reconfig.attr_stone_list = malloc(sizeof(EVattr_stone_struct) * in->count);
	    for (i=0 ; i < in->count; i++) {
		msg->u.flush_reconfig.attr_stone_list[i].stone = in->attr_stone_list[i].stone;
		msg->u.flush_reconfig.attr_stone_list[i].attr_str = strdup(in->attr_stone_list[i].attr_str);
	    }
	}
	break;
    }
    default:
	printf("MEssage type bad, value is %d  %d\n", msg_type, msg->msg_type);
	assert(FALSE);
    }
    msg->next = NULL;
    if (master->queued_messages == NULL) {
	master->queued_messages = msg;
    } else {
	EVmaster_msg_ptr last = master->queued_messages;
	while (last->next != NULL) last = last->next;
	last->next = msg;
    }
    if (master->cm->control_list->server_thread != 0) {
	CMwake_server_thread(master->cm);
    } else {
	handle_queued_messages(master->cm, master);
    }
}

static void
dfg_master_msg_handler(CManager cm, CMConnection conn, void *vmsg, 
		       void *client_data, attr_list attrs)
{
    EVmaster master = (EVmaster)((uintptr_t)client_data & (~0x7));
    EVmaster_msg_type msg_type = (EVmaster_msg_type) ((uintptr_t)client_data & 0x7);
    queue_master_msg(master, vmsg, msg_type, conn, /*copy*/1);
    /* we'll handle this in the poll handler */
}

static void
handle_deploy_ack_wrapper(EVmaster master, EVmaster_msg_ptr mmsg)
{
    CManager cm = master->cm;
    CManager_unlock(cm);
    handle_deploy_ack(master, mmsg);
    CManager_lock(cm);
}

static void
handle_deploy_ack(EVmaster master, EVmaster_msg_ptr mmsg)
{
    EVdeploy_ack_ptr msg =  &mmsg->u.deploy_ack;
    CManager cm = master->cm;
    EVdfg dfg = master->dfg;
    master->dfg->deploy_ack_count++;
    CMtrace_out(cm, EVdfgVerbose, "Client %s reports deployed, count %d\n", msg->node_id, master->dfg->deploy_ack_count);
    if ((master->dfg->deploy_ack_count == dfg->master->node_count) && (dfg->deploy_ack_condition != -1)) {
	CMtrace_out(cm, EVdfgVerbose, "That was the last one, Signalling %d\n", dfg->deploy_ack_condition);
	CMtrace_out(cm, EVdfgVerbose, "EVDFG exit deploy ack handler -  master DFG state is %s\n", str_state[master->state]);
	CMCondition_signal(cm, master->dfg->deploy_ack_condition);
	master->dfg->deploy_ack_condition = -1;
	assert(master->state == DFG_Starting);
	master->state = DFG_Running;
	CMtrace_out(cm, EVdfgVerbose, "EVDFG  -  master DFG state set to %s\n", str_state[master->state]);
    } else {
      if (master->state == DFG_Reconfiguring) {
	  if (master->dfg->deploy_ack_count == dfg->master->node_count) {
	      master->state = DFG_Running;
	      CMtrace_out(cm, EVdfgVerbose, "EVDFG after reconfiguration -  master DFG state set to %s\n", str_state[master->state]);
	  } else {
	      CMtrace_out(cm, EVdfgVerbose, "EVDFG reconfiguration in progress.  Deploy ack count %d, -  master DFG state set remains %s\n", master->dfg->deploy_ack_count, 
			  str_state[master->state]);
	  }
      }
    }
    CMtrace_out(cm, EVdfgVerbose, "EVDFG exit deploy ack handler -  master DFG state is %s\n", str_state[master->state]);
}

extern void INT_EVdfg_reconfig_transfer_events(EVdfg dfg, int src_stone_index, int src_port, int dest_stone_index, int dest_port) 
{
	
    if (dfg->transfer_events_count == 0) {
	dfg->transfer_events_list = malloc(sizeof(int *));
    } else {
	dfg->transfer_events_list = realloc(dfg->transfer_events_list, (dfg->transfer_events_count + 1) * sizeof(int *));
    }
	
    dfg->transfer_events_list[dfg->transfer_events_count] = malloc(4 * sizeof(int));
	
    dfg->transfer_events_list[dfg->transfer_events_count][0] = src_stone_index;
    dfg->transfer_events_list[dfg->transfer_events_count][1] = src_port;
    dfg->transfer_events_list[dfg->transfer_events_count][2] = dest_stone_index;
    dfg->transfer_events_list[dfg->transfer_events_count][3] = dest_port;
	
    ++dfg->transfer_events_count;
}

static void
perform_deployment(EVdfg dfg)
{
    int i;
    EVmaster master = dfg->master;

    if (dfg->master->sig_reconfig_bool == 0) {
	assert(master->state == DFG_Joining);
	master->state = DFG_Starting;
	CMtrace_out(master->cm, EVdfgVerbose, "EVDFG check all nodes registered -  master DFG state is %s\n", str_state[master->state]);
	add_bridge_stones(dfg, dfg->working_state);
	dfg->deploy_ack_count = 0;
	if (dfg->deploy_ack_condition == -1) {
	    dfg->deploy_ack_condition = INT_CMCondition_get(master->cm, NULL);
	}
	for (i=0; i < dfg->master->node_count; i++) {
	    deploy_to_node(dfg, i, dfg->working_state);
	    dfg->master->nodes[i].needs_ready = 1;   /* everyone needs a ready the first time through */
	}
    } else {
        CMtrace_out(master->cm, EVdfgVerbose, "EVDFG perform_deployment -  master DFG state set to %s\n", str_state[master->state]);
	assert(master->state == DFG_Reconfiguring);
	assign_actions_to_nodes(dfg->working_state, dfg->master);
	add_bridge_stones(dfg, dfg->working_state);
	add_unfreeze_actions(dfg, dfg->working_state);
	if (CMtrace_on(master->cm, EVdfgVerbose)) {
	    fprintf(master->cm->CMTrace_file, "EVDFG pending actions to be done on nodes\n");
	    for (i=0; i < dfg->working_state->pending_action_count; i++) {
		EVdfg_config_action act = dfg->working_state->pending_action_queue[i];
		fdump_dfg_config_action(master->cm->CMTrace_file, act);
	    }
	}
	for (i=dfg->master->old_node_count; i < dfg->master->node_count; i++) {
	    deploy_to_node(dfg, i, dfg->working_state);
	    dfg->master->nodes[i].needs_ready = 1;   /* everyone needs a ready the first time through */
	    remove_actions_for_node(dfg->working_state, i, master);
	}
	perform_actions_on_nodes(dfg->working_state, dfg->master);
	free(dfg->working_state->pending_action_queue);
	dfg->working_state->pending_action_queue = NULL;
    }
}

static void
wait_for_deploy_acks(EVdfg dfg)
{
    if (dfg->deploy_ack_count != dfg->master->node_count) {
	if (dfg->deploy_ack_condition != -1)  {
	    CManager_unlock(dfg->master->cm);
	    CMCondition_wait(dfg->master->cm, dfg->deploy_ack_condition);
	    CManager_lock(dfg->master->cm);
	}
    }
}

static void
signal_ready(EVdfg dfg)
{
    int i;
    CMFormat ready_msg = INT_CMlookup_format(dfg->master->cm, EVdfg_ready_format_list);
    EVready_msg msg;
    CMtrace_out(dfg->master->cm, EVdfgVerbose, "Master signaling DFG %p ready for operation\n",
		dfg);
    for (i=0; i < dfg->master->node_count; i++) {
	if (!dfg->master->nodes[i].needs_ready) {
	    CMtrace_out(dfg->master->cm, EVdfgVerbose, "Master - ready not required for node %d \"%s\"\n", i, 
			dfg->master->nodes[i].name);
	    continue;
	}
	if (dfg->master->nodes[i].conn != NULL) {
	    msg.node_id = i;
	    INT_CMwrite(dfg->master->nodes[i].conn, ready_msg, &msg);
	    CMtrace_out(dfg->master->cm, EVdfgVerbose, "Master - ready sent to node %d \"%s\"\n", i, 
			dfg->master->nodes[i].name);
	} else {
	    if (!dfg->master->nodes[i].self) {
		printf("Failure, no connection, not self, node %d\n", i);
		exit(1);
	    }
	    CManager_unlock(dfg->master->cm);
	    msg.node_id = i;
	    CMtrace_out(dfg->master->cm, EVdfgVerbose, "Master DFG %p is ready, local signalling %d\n", dfg, dfg->client->ready_condition);
	    dfg_ready_handler(dfg->master->cm, NULL, &msg, dfg->client, NULL);
	    CManager_lock(dfg->master->cm);
	}
	dfg->master->nodes[i].needs_ready = 0;
    }
}

static void
possibly_signal_shutdown(EVmaster master, int value, CMConnection conn)
{
    int i;
    CMFormat shutdown_msg = INT_CMlookup_format(master->cm, EVclient_shutdown_format_list);
    EVshutdown_msg msg;
    int status = STATUS_SUCCESS;
    int shutdown = 1;
    int force_shutdown = 0;
    int signal_from_client = -1;
    assert(CManager_locked(master->cm));
    CMtrace_out(master->cm, EVdfgVerbose, "possibly signal shutdown\n");
    for (i=0; i < master->node_count; i++) {
	if ((conn == NULL) && master->nodes[i].self) {
	    /* we're the master and node i */
	    signal_from_client = i;
	} else if (conn == master->nodes[i].conn) {
	    signal_from_client = i;
	}
    }
	
    if ((value >= 0) && ((value & STATUS_FORCE) == STATUS_FORCE)) {
	force_shutdown = 1;
	value ^= STATUS_FORCE;   /* yes, that's xor-assign */
    }
    if (force_shutdown) {
	CMtrace_out(master->cm, EVdfgVerbose, "Client %d signals %d, forces shutdown\n", signal_from_client, value);
    } else {
	CMtrace_out(master->cm, EVdfgVerbose, "Client %d signals %d, See if we're all ready to signal shutdown\n", signal_from_client, value);
    }
    master->nodes[signal_from_client].shutdown_status_contribution = value;
    int contributed_status = 0;
    for (i=0; i < master->node_count; i++) {
	CMtrace_out(master->cm, EVdfgVerbose, "NODE %d status is :", i);
	switch (master->nodes[i].shutdown_status_contribution) {
	case STATUS_UNDETERMINED:
	    CMtrace_out(master->cm, EVdfgVerbose, "NOT READY FOR SHUTDOWN\n");
	    shutdown = 0;
	    break;
	case STATUS_NO_CONTRIBUTION:
	    CMtrace_out(master->cm, EVdfgVerbose, "READY for shutdown, no status\n");
	    break;
	case STATUS_SUCCESS:
	    CMtrace_out(master->cm, EVdfgVerbose, "READY for shutdown, SUCCESS\n");
	    contributed_status = 1;
	    break;
	case STATUS_FAILED:
	    CMtrace_out(master->cm, EVdfgVerbose, "ALREADY FAILED\n");
	    contributed_status = 1;
	    break;
	default:
	    CMtrace_out(master->cm, EVdfgVerbose, "READY for shutdown, FAILURE %d\n",
			master->nodes[i].shutdown_status_contribution);
	    status |= master->nodes[i].shutdown_status_contribution;
	    break;
	}
    }
    if (force_shutdown) {
	shutdown = 1;
	status = value;
	CMtrace_out(master->cm, EVdfgVerbose, "DFG undergoing forced shutdown\n");
    }
    if (!shutdown) {
	CMtrace_out(master->cm, EVdfgVerbose, "DFG not ready for shutdown\n");
	return;
    }
    if (!contributed_status) {
	CMtrace_out(master->cm, EVdfgVerbose, "DFG nobody has contributed status - not ready for shutdown\n");
	return;
    }
    CMtrace_out(master->cm, EVdfgVerbose, "DFG shutdown with value %d\n", status);
    master->state = DFG_Shutting_Down;
    CMtrace_out(master->cm, EVdfgVerbose, "EVDFG possibly signal shutdown -  master DFG state is %s\n", str_state[master->state]);
    msg.value = status;
    for (i=0; i < master->node_count; i++) {
	if (master->nodes[i].conn != NULL) {
	    INT_CMwrite(master->nodes[i].conn, shutdown_msg, &msg);
	    CMtrace_out(master->cm, EVdfgVerbose, "DFG shutdown message sent to client \"%s\"(%d)\n", master->nodes[i].canonical_name, i);
	} else {
	    if (!master->nodes[i].self) {
		printf("Failure, no connection, not self, node %d\n", i);
		exit(1);
	    }
	}
    }
    /* sleep briefly to let the shutdown messages go out */
    INT_CMsleep(master->cm, 1);
    for (i=0; i < master->node_count; i++) {
	if (master->nodes[i].conn != NULL) {
	    if (master->nodes[i].conn->closed) {
		INT_CMConnection_dereference(master->nodes[i].conn);
	    } else {
		INT_CMConnection_close(master->nodes[i].conn);
	    }
	    master->nodes[i].conn = NULL;
	}
    }
    if (master->client) {
	master->client->shutdown_value = status;
	i = 0;
	master->client->already_shutdown = 1;
	while(master->client->shutdown_conditions && (master->client->shutdown_conditions[i] != -1)) {
	    CMtrace_out(master->cm, EVdfgVerbose, "Client %d shutdown signalling %d\n", master->client->my_node_id, master->client->shutdown_conditions[i]);
	    INT_CMCondition_signal(master->cm, master->client->shutdown_conditions[i++]);
	}
	CMtrace_out(master->cm, EVdfgVerbose, "Master DFG shutdown\n");
    }
}

extern void INT_EVmaster_node_join_handler(EVmaster master, EVmasterJoinHandlerFunc func)
{
    master->node_join_handler = func;
}

extern void INT_EVmaster_node_fail_handler(EVmaster master, EVmasterFailHandlerFunc func)
{
    master->node_fail_handler = func;
}

extern void INT_EVmaster_node_reconfig_handler(EVmaster master, EVmasterReconfigHandlerFunc func)
{
    master->node_reconfig_handler = func;
}

static void
free_dfg_state(EVdfg_configuration state)
{
    int i;
    for (i=0; i < state->stone_count; i++) {
    	if (state->stones[i]->out_links) free(state->stones[i]->out_links);
    	if (state->stones[i]->in_links) free(state->stones[i]->in_links);
    	if (state->stones[i]->action) free(state->stones[i]->action);
    	if (state->stones[i]->extra_actions) {
    	    int j;
    	    for (j=0; j < state->stones[i]->action_count-1; j++) {
    		free(state->stones[i]->extra_actions[j]);
    	    }
    	    free(state->stones[i]->extra_actions);
    	}
    	if (state->stones[i]->attrs) free_attr_list(state->stones[i]->attrs);
    	free(state->stones[i]);
    }
    if (state->pending_action_queue) free(state->pending_action_queue);
    free(state->stones);
    free(state);
}

static EVdfg_configuration
copy_dfg_state(EVdfg_configuration state)
{
    int i;
    EVdfg_configuration ret = malloc(sizeof(*ret));
    memset(ret, 0, sizeof(*ret));
    ret->stone_count = state->stone_count;
    ret->stones = malloc(sizeof(ret->stones[0])*state->stone_count);
    for (i=0; i < state->stone_count; i++) {
	ret->stones[i] = malloc(sizeof(*ret->stones[i]));
	*(ret->stones[i]) = *(state->stones[i]);
    	if (state->stones[i]->out_links) {
	    ret->stones[i]->out_links = malloc(sizeof(ret->stones[i]->out_links[0]) * state->stones[i]->out_count);
	    memcpy(ret->stones[i]->out_links, state->stones[i]->out_links,
		   sizeof(ret->stones[i]->out_links[0]) * ret->stones[i]->out_count);
	}
    	if (state->stones[i]->in_links) {
	    ret->stones[i]->in_links = malloc(sizeof(ret->stones[i]->in_links[0]) * state->stones[i]->in_count);
	    memcpy(ret->stones[i]->in_links, state->stones[i]->in_links,
		   sizeof(ret->stones[i]->in_links[0]) * ret->stones[i]->in_count);
	}
    	if (state->stones[i]->action) ret->stones[i]->action = strdup(state->stones[i]->action);
    	if (state->stones[i]->extra_actions) {
    	    int j;
	    ret->stones[i]->extra_actions = malloc(sizeof(char*) * state->stones[i]->action_count);
    	    for (j=0; j < state->stones[i]->action_count-1; j++) {
		if (state->stones[i]->extra_actions[j]) 
		    ret->stones[i]->extra_actions[j] = strdup(state->stones[i]->extra_actions[j]);
    	    }
    	}
    	if (state->stones[i]->attrs) add_ref_attr_list(state->stones[i]->attrs);
    }
    return ret;
}

static void
fdump_dfg_gml(FILE* out, EVdfg_configuration state)
{
    int i;
    char *prefix = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<!-- This file was written by the JAVA GraphML Library.-->\n\
<graphml\n\
 xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n\
 xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n\
 xmlns:y=\"http://www.yworks.com/xml/graphml\"\n\
 xmlns:yed=\"http://www.yworks.com/xml/yed/3\"\n\
 xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd\">\n\
  <key id=\"d0\" for=\"node\" attr.name=\"color\" attr.type=\"string\">\n\
    <default>yellow</default>\n\
  </key>\n\
  <key for=\"node\" id=\"d1\" yfiles.type=\"nodegraphics\"/>\n\
  <graph id=\"G\" edgedefault=\"directed\">\n";

    fprintf(out, "%s", prefix);
    for (i=0; i < state->stone_count; i++) {
	int j;
	fprintf(out, "<node id=\"n%d\" name=\"stone%d\">\n", i, i);
	for (j=0; j< state->stones[i]->out_count; j++) {
	    fprintf(out, "<port name=\"P%d\"/>\n", j);
	}
	fprintf(out, "      <data key=\"d1\">\n        <y:ShapeNode>\n            <y:NodeLabel>S%d</y:NodeLabel>                    <!-- label text -->\n        </y:ShapeNode>\n      </data>\n", i);
	fprintf(out, "</node>\n"); 
	for (j=0; j< state->stones[i]->out_count; j++) {
	    fprintf(out, "<edge id=\"n%de%d\" source=\"n%d\" sourceport=\"P%d\" target=\"n%d\">\n", i, ~0x80000000 & state->stones[i]->out_links[j], i,  j, ~0x80000000 & state->stones[i]->out_links[j]);
	    fprintf(out, "</edge>\n");
	}
    }
    fprintf(out, "</graph>\n</graphml>\n");
}

static void
dump_dfg_gml(EVdfg_configuration state)
{
    fdump_dfg_gml(stdout, state);
}

void
INT_EVdfg_dump_graph(EVdfg_state_type which, EVdfg dfg)
{
    switch(which) {
    case EVdfgWorking:
	dump_dfg_gml(dfg->working_state);
	break;
    case EVdfgDeployed:
	dump_dfg_gml(dfg->deployed_state);
	break;
    }
}
static void
check_all_nodes_registered(EVmaster master)
{
    int i;
    EVdfg dfg = master->dfg;
    if (master->node_join_handler != NULL) {
	EVint_node_list node = &master->nodes[master->node_count-1];
	CManager_unlock(master->cm);
	(master->node_join_handler)(master, node->name, NULL, NULL);
	CManager_lock(master->cm);
	dfg = master->dfg;
	if ((dfg == NULL) || (dfg->realized == 0) || 
	    (dfg->realized == 1 && master->reconfig == 1)) return;
    } else {
	/* must be static node list */
	for(i=0; i<master->node_count; i++) {
	    if (!master->nodes[i].self && (master->nodes[i].conn == NULL)) {
		return;
	    }
	}
    }
	
    if (master->no_deployment == 0) {
	perform_deployment(dfg);
	wait_for_deploy_acks(dfg);
    }
    master->no_deployment = 0;
    if (dfg->deployed_state) {
	free_dfg_state(dfg->deployed_state);
    }
    dfg->deployed_state = dfg->working_state;
    dfg->working_state = copy_dfg_state(dfg->deployed_state);
    
    signal_ready(dfg);
    dfg->deployed_stone_count = dfg->stone_count;
    master->old_node_count = master->node_count;
}
