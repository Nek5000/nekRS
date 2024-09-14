
#include "config.h"
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <time.h>

#include "evpath.h"
#include "cm_internal.h"
#include "response.h"
#include "dlloader.h"

extern 
FMFormat EVregister_format_set(CManager cm, FMStructDescList list);
static void reference_event(event_item *event);
//static void dump_action(stone_type stone, response_cache_element *resp, 
//			int a, const char *indent);
static void dump_stone(stone_type stone);
static void fdump_stone(FILE* out, stone_type stone);
static int is_bridge_stone(CManager cm, EVstone stone_num);
static event_item *
dequeue_event(CManager cm, stone_type stone);

static const char *action_str[] = { "Action_NoAction","Action_Bridge", "Action_Thread_Bridge", "Action_Terminal", "Action_Filter", "Action_Immediate", "Action_Multi", "Action_Decode", "Action_Encode_to_Buffer", "Action_Split", "Action_Store", "Action_Congestion", "Action_Source"};


static void
fprint_stone_identifier(FILE *out, event_path_data evp, int stone_num)
{
    int local_stone_num = -1;
    int global_stone_num = -1;
    if ((stone_num & 0x80000000) == 0x80000000) {
	local_stone_num = lookup_local_stone(evp, stone_num);
	global_stone_num = stone_num;
    } else {
	int i;
	local_stone_num = stone_num;
	for (i=0; i < evp->stone_lookup_table_size; i++) {
	    if (evp->stone_lookup_table[i].local_id == stone_num) {
		global_stone_num = evp->stone_lookup_table[i].global_id;
		break;
	    }
	}
    }
    fprintf(out, "local stone number %x", local_stone_num);
    if (global_stone_num != -1) {
	fprintf(out, " (global %x)", global_stone_num);
    }
}

static void
print_stone_identifier(event_path_data evp, int stone_num)
{
    fprint_stone_identifier(stdout, evp, stone_num);
}

extern int
lookup_local_stone(event_path_data evp, int stone_num)
{
    int i;
    int ret = -1;
    for (i=0; i < evp->stone_lookup_table_size; i++) {
	if (evp->stone_lookup_table[i].global_id == stone_num) {
	    ret = evp->stone_lookup_table[i].local_id;
	    break;
	}
    }
    if (ret == -1) {
	printf("EVPATH: Invalid GLOBAL stone ID %x\n", stone_num);
	return -1;
    }
    return ret;
}

extern void
remove_stone_from_lookup(event_path_data evp, int stone_num)
{
    int i, stone = -1;
    for (i=0; i < evp->stone_lookup_table_size; i++) {
	if (evp->stone_lookup_table[i].global_id == stone_num) {
	    stone = i;
	    break;
	}
    }
    if (stone == -1) return;
    for (i=stone; i < evp->stone_lookup_table_size-1; i++) {
	evp->stone_lookup_table[i] = evp->stone_lookup_table[i+1];
    }
}

extern void
INT_CMadd_stone_to_global_lookup(CManager cm, int stone_num, int global_stone_num)
{
    event_path_data evp = cm->evp;
    if ((global_stone_num & 0x80000000) != 0x80000000) {
        fprintf(stderr, "Global stone num must have 32nd bit set.  Value provided was %x\n", global_stone_num);
        fprintf(stderr, "Ignoring call to CMadd_stone_to_global_lookup for stone %d\n", stone_num);
        return;
    }
    add_stone_to_lookup(evp, stone_num, global_stone_num);
}

extern void
add_stone_to_lookup(event_path_data evp, int stone_num, int global_stone_num)
{
    int base = evp->stone_lookup_table_size;
    if (evp->stone_lookup_table_size == 0) {
	evp->stone_lookup_table = 
	    malloc(sizeof(evp->stone_lookup_table[0]));
    } else {
	evp->stone_lookup_table = 
	    realloc(evp->stone_lookup_table,
		    sizeof(evp->stone_lookup_table[0]) * (1+base));
    }
    evp->stone_lookup_table[base].global_id = global_stone_num;
    evp->stone_lookup_table[base].local_id = stone_num;
    evp->stone_lookup_table_size++;
}

extern int
lookup_global_stone(event_path_data evp, int stone_num)
{
    int i;
    int ret = -1;
    for (i=0; i < evp->stone_lookup_table_size; i++) {
	if (evp->stone_lookup_table[i].local_id == stone_num) {
	    ret = evp->stone_lookup_table[i].global_id;
	    break;
	}
    }
    if (ret == -1) {
	printf("EVPATH: stone ID %x has no global counterpart\n", stone_num);
	return -1;
    }
    return ret;
}

stone_type
stone_struct(event_path_data evp, int stone_num)
{
    int was_global = 0;
    if ((stone_num & 0x80000000) == 0x80000000) {
	stone_num = lookup_local_stone(evp, stone_num);
	was_global++;
    }
	    
    if (evp->stone_count <= stone_num - evp->stone_base_num) {
	printf("EVPATH: Invalid stone ID %x\n", stone_num);
        return NULL;
    }
    if (was_global && ((evp->stone_map[stone_num - evp->stone_base_num] == NULL) ||
		       (evp->stone_map[stone_num - evp->stone_base_num]->local_id == -1))) {
	printf("EVPATH: Invalid stone ID %d (local ID -1)\n", stone_num);
        return NULL;
    }
    return evp->stone_map[stone_num - evp->stone_base_num];
}

void
INT_EVadd_dll_search_dir(char *path_string)
{
    CMdladdsearchdir(path_string);
}

EVstone
INT_EValloc_stone(CManager cm)
{
    event_path_data evp = cm->evp;
    int stone_num = evp->stone_count;
    stone_type stone;

    evp->stone_map = realloc(evp->stone_map, 
			     (evp->stone_count + 1) * sizeof(evp->stone_map[0]));
    evp->stone_map[stone_num] = malloc(sizeof(*stone));
    stone = evp->stone_map[stone_num];
    stone_num += evp->stone_base_num;
    memset(stone, 0, sizeof(*stone));
    stone->local_id = stone_num;
    stone->default_action = -1;
    stone->response_cache_count = 0;
    stone->response_cache = NULL;
    stone->is_frozen = 0;
    stone->is_processing = 0;
    stone->is_outputting = 0;
    stone->is_draining = 0;
    stone->queue = malloc(sizeof(queue_struct));
    stone->queue->queue_tail = stone->queue->queue_head = NULL;
    stone->new_enqueue_flag = 0;
    stone->write_callback = -1;
    stone->proto_actions = NULL;
    stone->stone_attrs = CMcreate_attr_list(cm);
    stone->output_count = 0;
    stone->output_stone_ids = malloc(sizeof(int));
    stone->output_stone_ids[0] = -1;

    stone->queue_size = 0;
    stone->is_stalled = 0;
    stone->stall_from = Stall_None;
    stone->last_remote_source = NULL;
    stone->squelch_depth = 0;
    stone->unstall_callbacks = NULL;
    evp->stone_count++;
    return stone_num;
}

static void
empty_queue(event_path_data evp, queue_ptr queue)
{
    while(queue->queue_head != NULL && queue->queue_tail != NULL) {
	queue_item *tmp = queue->queue_head;
	return_event(evp, queue->queue_head->item);
        if(queue->queue_head == queue->queue_tail) {
	    queue->queue_head = NULL;
	    queue->queue_tail = NULL;
	} else {
	    queue->queue_head = queue->queue_head->next;
	}
	free(tmp);
    }
}

/* {{{ storage_queue_* */
static storage_queue_ptr 
storage_queue_init(CManager cm, storage_queue_ptr queue,
        storage_queue_ops_ptr ops, attr_list attrs) {
    memset(queue, 0, sizeof *queue);
    queue->ops = ops;
    if (queue->ops->init)
        (queue->ops->init)(cm, queue, attrs);
    return queue;
}

static void
storage_queue_cleanup(CManager cm, storage_queue_ptr queue) {
    if (queue->ops->cleanup)
        (queue->ops->cleanup)(cm, queue);
}

static void
storage_queue_empty(CManager cm, storage_queue_ptr queue) {
    (queue->ops->empty)(cm, queue);
}

static void
storage_queue_enqueue(CManager cm, storage_queue_ptr queue, event_item *item) {
    (queue->ops->enqueue)(cm, queue, item);
}

static event_item *
storage_queue_dequeue(CManager cm, storage_queue_ptr queue) {
    return (queue->ops->dequeue)(cm, queue);
}
/* }}} */

static void
free_response_cache(stone_type stone)
{
    int i;
    for(i = 0; i < stone->response_cache_count; i++) {
	response_cache_element *resp = &stone->response_cache[i];
	switch(resp->action_type) {
	case Action_Decode:
	    if (resp->o.decode.context) {
		free_FFSContext(resp->o.decode.context);
		resp->o.decode.context = NULL;
	    }
	    break;
	case Action_Immediate:
	    if (resp->o.imm.free_func) {
		(resp->o.imm.free_func)(resp->o.imm.client_data);
	    }
	    break;
	case Action_Multi:
	    if (resp->o.multi.free_func) {
		(resp->o.multi.free_func)(resp->o.multi.client_data);
	    }
	    break;
	default:
	    break;
	}
    }
    if (stone->response_cache) free(stone->response_cache);
}

void
INT_EVfree_stone(CManager cm, EVstone stone_num)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int i;

    stone = stone_struct(evp, stone_num);
    CMtrace_out(cm, CMFreeVerbose, "Freeing stone %d\n", stone_num);
    if (stone == NULL) return;
    if (stone->local_id == -1) return;
    if (stone->periodic_handle != NULL) {
	INT_CMremove_task(stone->periodic_handle);
	stone->periodic_handle = NULL;
    }
    for(i = 0; i < stone->proto_action_count; i++) {
	proto_action *act = &stone->proto_actions[i];
	if (act->attrs != NULL) {
	    INT_CMfree_attr_list(cm, act->attrs);
	}
	if (act->matching_reference_formats != NULL) 
	    free(act->matching_reference_formats);
	switch(act->action_type) {
	case Action_NoAction:
	case Action_Source:
	case Action_Encode_to_Buffer:
	case Action_Thread_Bridge:
	    break;
	case Action_Bridge:
	    if (act->o.bri.conn) {
		CMtrace_out(cm, CMFreeVerbose, "Closing and dereferencing conn %p in free stone\n", act->o.bri.conn);
		INT_CMConnection_dereference(act->o.bri.conn);
	    }
            if (act->o.bri.remote_contact) {
                free_attr_list(act->o.bri.remote_contact);
                act->o.bri.remote_contact = NULL;
            }
	    if (act->o.bri.remote_path) {
		free(act->o.bri.remote_path);
                act->o.bri.remote_path = NULL;
            }
	    break;
	case Action_Terminal:
	    break;
	case Action_Filter:
	    break;
	case Action_Decode:
	    if (act->o.decode.context) {
		free_FFSContext(act->o.decode.context);
		act->o.decode.context = NULL;
	    }
	    break;
	case Action_Split:
	    break;
	case Action_Immediate:
        case Action_Multi:
        case Action_Congestion:
	    if (act->o.imm.mutable_response_data != NULL) {
		response_data_free(cm, act->o.imm.mutable_response_data);
	    }
	    break;
        case Action_Store:
            storage_queue_cleanup(cm, &act->o.store.queue);
            break; 
	}
    }
    while (stone->queue->queue_head != NULL) {
      event_item *event = dequeue_event(cm, stone);
      return_event(evp, event);
    }
    if (stone->proto_actions != NULL) free(stone->proto_actions);
    if (stone->response_cache != NULL) free_response_cache(stone);
    
    free(stone->queue);
      
    /* XXX unsquelch senders */
    stone->queue = NULL;
    stone->proto_action_count = 0;
    stone->proto_actions = NULL;
    if (stone->stone_attrs != NULL) {
	INT_CMfree_attr_list(cm, stone->stone_attrs);
	stone->stone_attrs = NULL;
    }
    free(stone->output_stone_ids);
    remove_stone_from_lookup(evp, stone_num);
    evp->stone_map[stone->local_id - evp->stone_base_num] = NULL;
    stone->local_id = -1;
    free(stone);
}

EVaction
add_proto_action(CManager cm, stone_type stone, proto_action **act)
{
    int proto_action_num = stone->proto_action_count;
    (void)cm;
    stone->proto_actions = realloc(stone->proto_actions, 
				   (proto_action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[proto_action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    ++stone->proto_action_count;
    *act = &stone->proto_actions[proto_action_num];
    return proto_action_num;
}

static void
clear_response_cache(stone_type stone)
{
    stone->response_cache_count = 0;
    /* GSE  free response entitites */
    if (stone->response_cache) free_response_cache(stone);
    stone->response_cache = NULL;
}

void
fix_response_cache(stone_type stone)
{
    int i;
    for (i=stone->response_cache_count-1; i >= 0; i--) {
	int j;
	FMFormat ref = stone->response_cache[i].reference_format;
	for (j=0; j < i; j++) {
	    if (((stone->response_cache[j].reference_format == ref) ||
		 (stone->response_cache[j].reference_format == NULL)) && 
		(stone->response_cache[j].action_type == Action_NoAction)) {
		/* need to kill this */
		memmove(&stone->response_cache[j], &stone->response_cache[j+1], 
			sizeof(stone->response_cache[0]) * (stone->response_cache_count - j - 1));
		stone->response_cache_count--;
	    }			
	}
    }
}

EVstone
INT_EVcreate_terminal_action(CManager cm, FMStructDescList format_list, 
			     EVSimpleHandlerFunc handler, void *client_data)
{
    EVstone stone = INT_EValloc_stone(cm);
    (void) INT_EVassoc_terminal_action(cm, stone, format_list, 
				       handler, client_data);
    return stone;
}    

EVstone
INT_EVcreate_stone_action(CManager cm, char *action_spec)
{
    EVstone stone = INT_EValloc_stone(cm);
    (void) INT_EVassoc_general_action(cm, stone, action_spec, NULL);
    return stone;
}    

EVaction
INT_EVassoc_terminal_action(CManager cm, EVstone stone_num, 
			    FMStructDescList format_list, EVSimpleHandlerFunc handler,
			    void *client_data)
{
    event_path_data evp = cm->evp;
    int action_num;
    stone_type stone;
    int proto_action_num;

    stone = stone_struct(evp, stone_num);
    proto_action_num = stone->proto_action_count;
    stone->proto_actions = realloc(stone->proto_actions, 
				   (proto_action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[proto_action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    stone->proto_actions[proto_action_num].input_format_requirements =
	format_list;
    stone->proto_actions[proto_action_num].action_type = Action_Terminal;
    stone->proto_actions[proto_action_num].o.term.handler = handler;
    stone->proto_actions[proto_action_num].o.term.client_data = client_data;
    stone->proto_actions[proto_action_num].matching_reference_formats = NULL;
    action_num = stone->response_cache_count;
    stone->response_cache = realloc(stone->response_cache, (action_num + 1) * 
				   sizeof(stone->response_cache[0]));
    memset(&stone->response_cache[action_num], 0, sizeof(stone->response_cache[0]));
    if (format_list != NULL) {
	stone->proto_actions[proto_action_num].data_state = Requires_Decoded;
	stone->proto_actions[proto_action_num].matching_reference_formats = 
	    malloc(2*sizeof(FMFormat));
	stone->proto_actions[proto_action_num].matching_reference_formats[0] = 
	    EVregister_format_set(cm, format_list);
	stone->proto_actions[proto_action_num].matching_reference_formats[1] = NULL;
    } else {
	stone->proto_actions[proto_action_num].data_state = Requires_Contig_Encoded;
	stone->default_action = action_num;
    }
    stone->response_cache[action_num].action_type = Action_Terminal;
    stone->response_cache[action_num].requires_decoded =
	stone->proto_actions[proto_action_num].data_state;
    if (stone->proto_actions[proto_action_num].matching_reference_formats) {
	stone->response_cache[action_num].reference_format = 
	    stone->proto_actions[proto_action_num].matching_reference_formats[0];
    } else {
	stone->response_cache[action_num].reference_format = NULL;
    }
    stone->response_cache[action_num].proto_action_id = proto_action_num;
    stone->proto_action_count++;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding Terminal action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, "\nStone dump->\n");
	fdump_stone(cm->CMTrace_file, stone);
    }
    return action_num;
}
    
EVaction
INT_EVassoc_raw_terminal_action(CManager cm, EVstone stone_num, 
				EVRawHandlerFunc handler,
				void *client_data)
{
    return INT_EVassoc_terminal_action(cm, stone_num, NULL,
				       (EVSimpleHandlerFunc)handler,
				       client_data);
}
    

EVstone
INT_EVcreate_immediate_action(CManager cm, char *action_spec, 
			      EVstone *target_list)
{
    int i = 0;
    EVstone stone = INT_EValloc_stone(cm);
    EVaction action = EVassoc_immediate_action(cm, stone, action_spec, NULL);
    while (target_list && (target_list[i] != 0)) {
	INT_EVaction_set_output(cm, stone, action, i, target_list[i]);
	i++;
    }
    return stone;
}

EVaction
INT_EVassoc_general_action(CManager cm, EVstone stone_num, char*action_spec,
    EVstone *output_list)
{
    event_path_data evp = cm->evp;
    EVaction ret = -1;
    switch (action_type(action_spec))
    {
    case Action_Immediate: {
	int i = 0;
	ret = INT_EVassoc_immediate_action(cm, stone_num, action_spec, NULL);
	while (output_list && (output_list[i] != -1)) {
	    INT_EVaction_set_output(cm, stone_num, ret, i, output_list[i]);
	    i++;
	}
	break;
    }
    case Action_Multi: {
	int i = 0;
	ret = INT_EVassoc_multi_action(cm, stone_num, action_spec, NULL);
	while (output_list && (output_list[i] != -1)) {
	    INT_EVaction_set_output(cm, stone_num, ret, i, output_list[i]);
	    i++;
	}
	break;
    }

    case Action_Bridge:
    {
	EVstone target;
	char *contact;
	attr_list attrs;
	parse_bridge_action_spec(action_spec, &target, &contact);
	attrs = attr_list_from_string(contact);
	ret = INT_EVassoc_bridge_action(cm, stone_num, attrs, target);
	free_attr_list(attrs);
	break;
    }
    case Action_Terminal: {
	char *name = action_spec+5;
	int i;
	for (i=0; i < evp->sink_handler_count; i++) {
	    if (strcmp(name, evp->sink_handlers[i].name) == 0) {
		ret = INT_EVassoc_terminal_action(cm, stone_num, 
						  evp->sink_handlers[i].format_list,
						  evp->sink_handlers[i].handler,
						  evp->sink_handlers[i].client_data);
		break;
	    }
	}
	if (i == evp->sink_handler_count) {
	    printf("Failed to find handler func \"%s\"\n", name);
	}
	break;
    }
    case Action_Split:
	ret = INT_EVassoc_split_action(cm, stone_num, output_list);
	break;
    case Action_Source: {
	char *name = action_spec+7;
	int i;
	for (i=0; i < evp->source_count; i++) {
	    if (strcmp(name, evp->sources[i].name) == 0) {
		evp->sources[i].src->local_stone_id = stone_num;
		ret = INT_EVassoc_split_action(cm, stone_num, output_list);
		break;
	    }
	}
	if (i == evp->source_count) {
	    printf("Failed to find source \"%s\"\n", name);
	}
	break;
    }
    default:
	printf("Missed case\n");
    } 
    return ret;
}

EVaction
INT_EVassoc_immediate_action(CManager cm, EVstone stone_num, 
			 char *action_spec, void *client_data)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    EVaction action_num;
    proto_action *act;
    stone_type stone;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = add_proto_action(cm, stone, &act);
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding Immediate action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, "\naction value is \"%s\"\n", action_spec);
    }
    stone->proto_actions = realloc(stone->proto_actions, (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].data_state = Requires_Decoded;
    stone->proto_actions[action_num].action_type = Action_Immediate;
    stone->proto_actions[action_num].o.imm.mutable_response_data = 
 	install_response_handler(cm, stone_num, action_spec, client_data, 
				 &stone->proto_actions[action_num].matching_reference_formats);
    if (stone->proto_actions[action_num].matching_reference_formats &&
	(stone->proto_actions[action_num].matching_reference_formats[0] == NULL)) {
	stone->default_action = action_num;
	stone->proto_actions[action_num].data_state = Accepts_All;
    }	
    clear_response_cache(stone);
    return action_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

EVaction
INT_EVassoc_multi_action(CManager cm, EVstone stone_num, 
			  char *action_spec, void *client_data)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    int action_num;
    stone_type stone;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;
    action_num = stone->proto_action_count;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding Multi action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, "\nmulti action is \"%s\"\n", action_spec);
    }

    stone->proto_actions = realloc(stone->proto_actions, (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].data_state = Requires_Decoded;
    stone->proto_actions[action_num].action_type = Action_Multi;
    stone->proto_actions[action_num].o.imm.mutable_response_data = 
 	install_response_handler(cm, stone_num, action_spec, client_data, 
				 &stone->proto_actions[action_num].matching_reference_formats);
   stone->proto_action_count++;
    clear_response_cache(stone);
    return action_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

EVaction
INT_EVassoc_congestion_action(CManager cm, EVstone stone_num, 
			      char *action_spec, void *client_data)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    int action_num;
    stone_type stone;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = stone->proto_action_count;
    CMtrace_out(cm, EVerbose, "Adding Congestion action %d to stone %x\n",
		action_num, stone_num);
    stone->proto_actions = realloc(stone->proto_actions, (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].data_state = Requires_Decoded;
    stone->proto_actions[action_num].action_type = Action_Congestion;
    stone->proto_actions[action_num].o.imm.mutable_response_data = 
 	install_response_handler(cm, stone_num, action_spec, client_data, 
				 &stone->proto_actions[action_num].matching_reference_formats);
    stone->proto_action_count++;
    clear_response_cache(stone);
    return action_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

EVstone
INT_EVassoc_filter_action(CManager cm, EVstone stone_num, 
		      FMStructDescList format_list, EVSimpleHandlerFunc handler,
		      EVstone out_stone_num, void *client_data)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    stone_type stone;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    int proto_action_num = stone->proto_action_count;
    stone->proto_actions = realloc(stone->proto_actions, 
				   (proto_action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[proto_action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    stone->proto_actions[proto_action_num].input_format_requirements =
	format_list;
    stone->proto_actions[proto_action_num].action_type = Action_Filter;
    stone->proto_actions[proto_action_num].o.term.handler = handler;
    stone->proto_actions[proto_action_num].o.term.client_data = client_data;
    stone->proto_actions[proto_action_num].o.term.target_stone_id = out_stone_num;
    stone->proto_actions[proto_action_num].data_state = Requires_Decoded;
    stone->proto_actions[proto_action_num].matching_reference_formats = NULL;
    if (format_list != NULL) {
	stone->proto_actions[proto_action_num].matching_reference_formats = 
	    malloc(2*sizeof(FMFormat));
	stone->proto_actions[proto_action_num].matching_reference_formats[0] = 
	    EVregister_format_set(cm, format_list);
	stone->proto_actions[proto_action_num].matching_reference_formats[1] = NULL;
    }	
    stone->proto_action_count++;
    clear_response_cache(stone);
    CMtrace_out(cm, EVerbose, "Adding filter action %d to stone %x\n",
		proto_action_num, stone_num);
    return proto_action_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

static storage_queue_ops storage_queue_default_ops;

EVaction 
INT_EVassoc_store_action(CManager cm, EVstone stone_num, EVstone out_stone,
                            int store_limit)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    proto_action *act;
    int action_num;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = add_proto_action(cm, stone, &act);

    act->data_state = Accepts_All;
    act->action_type = Action_Store;
    act->matching_reference_formats = malloc(sizeof(FMFormat));
    act->matching_reference_formats[0] = NULL; /* signal that we accept all formats */
    storage_queue_init(cm, &act->o.store.queue, &storage_queue_default_ops, NULL);
    act->o.store.target_stone_id = out_stone;
    act->o.store.max_stored = store_limit;
    act->o.store.num_stored = 0;
    clear_response_cache(stone);
    stone->default_action = action_num;

    return action_num;
}

EVstone
INT_EVcreate_store_action(CManager cm, EVstone out_stone, int store_limit)
{
    EVstone stone = INT_EValloc_stone(cm);
    (void) INT_EVassoc_store_action(cm, stone, out_stone, store_limit);
    return stone;
}

void
INT_EVclear_stored(CManager cm, EVstone stone_num, EVaction action_num)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    storage_queue_ptr queue;

    stone = stone_struct(evp, stone_num);
    if (!stone) return;

    queue = &stone->proto_actions[action_num].o.store.queue;

    storage_queue_empty(cm, queue);
}

typedef struct action_tracking_state {
    int last_active_stone;
    int events_in_play;
} *action_state;

static void
raw_enqueue_event(CManager cm, queue_ptr queue, event_item *event)
{
    queue_item *item;
    event_path_data evp = cm->evp;
    if (evp->queue_items_free_list == NULL) {
	item = malloc(sizeof(*item));
    } else {
	item = evp->queue_items_free_list;
	evp->queue_items_free_list = item->next;
    }
    item->item = event;
    item->handled = 0;
    reference_event(event);
    if (queue->queue_head == NULL) {
	queue->queue_head = item;
	queue->queue_tail = item;
	item->next = NULL;
    } else {
	queue->queue_tail->next = item;
	queue->queue_tail = item;
	item->next = NULL;
    }
}

static void backpressure_check(CManager, EVstone);

static void
enqueue_event(CManager cm, int stone_id, event_item *event)
{
    event_path_data evp = cm->evp;
    stone_type stone = stone_struct(evp, stone_id);
    action_state as = evp->as;
    if (as == NULL) {
	as = evp->as = malloc(sizeof(*as));
	memset(as, 0, sizeof(*as));
    }
    raw_enqueue_event(cm, stone->queue, event);
    stone->new_enqueue_flag = 1;
    stone->queue_size++;
    backpressure_check(cm, stone_id);
    as->last_active_stone = stone_id;
    as->events_in_play++;
}

static event_item *
raw_dequeue_event(CManager cm, queue_ptr q)
{
    event_path_data evp = cm->evp;
    queue_item *item = q->queue_head;
    event_item *event = NULL;
    if (item == NULL) return event;
    event = item->item;
    if (q->queue_head == q->queue_tail) {
	q->queue_head = NULL;
	q->queue_tail = NULL;
    } else {
	q->queue_head = q->queue_head->next;
    }   
    item->next = evp->queue_items_free_list;
    evp->queue_items_free_list = item;
    return event;
}

static event_item *
dequeue_event(CManager cm, stone_type stone)
{
    queue_ptr q = stone->queue;
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    event_item *event = NULL;
    event = raw_dequeue_event(cm, q);
    stone->queue_size--;
    as->events_in_play--;
    return event;
}

static event_item *
dequeue_item(CManager cm, stone_type stone, queue_item *to_dequeue)
{
    queue_ptr q = stone->queue;
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    queue_item *item = q->queue_head;
    event_item *event = NULL;
    
    assert(CManager_locked(cm));
    if (to_dequeue == NULL) return event;
    event = to_dequeue->item;
    if (q->queue_head == to_dequeue) {
	if (q->queue_head == q->queue_tail) {
	    q->queue_head = NULL;
	    q->queue_tail = NULL;
	} else {
	    q->queue_head = q->queue_head->next;
	}
    } else {
	queue_item *cur, *last;
	last = q->queue_head;
	cur = q->queue_head->next;
	while (cur != to_dequeue) {
	    last = cur;
	    cur = cur->next;
	}
        assert(cur == to_dequeue);
	item = cur;
	last->next = cur->next;
	if (cur == q->queue_tail) {
	    q->queue_tail = last;
	}
	cur = q->queue_head;
	while (cur != NULL) {
	    cur = cur->next;
	}
    }
    item->next = evp->queue_items_free_list;
    evp->queue_items_free_list = item;
    stone->queue_size--;
    as->events_in_play--;
    return event;
}

void
EVdiscard_queue_item(CManager cm, int s, queue_item *item) {
    stone_type stone = stone_struct(cm->evp, s);
    event_item *event;
    event = dequeue_item(cm, stone, item);
    if (event) return_event(cm->evp, event);
}

static void encode_event(CManager, event_item*);

/* {{{ storage_queue_default_* */
extern void 
ensure_ev_owned(CManager cm, event_item *event)
{
    if (event->contents == Event_App_Owned && !event->free_func) {
        encode_event(cm, event);
        event->decoded_event = NULL;
        event->event_encoded = 1;
        event->contents = Event_CM_Owned;
        assert(event->encoded_event);
    }
}

static void
storage_queue_default_empty(CManager cm, storage_queue_ptr queue) {
    (void)cm;
    empty_queue(cm->evp, &queue->u.queue); 
}

static void
storage_queue_default_enqueue(CManager cm, storage_queue_ptr queue, event_item *item) {
    ensure_ev_owned(cm, item);
    raw_enqueue_event(cm, &queue->u.queue, item);
}

static event_item *
storage_queue_default_dequeue(CManager cm, storage_queue_ptr queue) {
    event_item *ev = raw_dequeue_event(cm, &queue->u.queue);
    return ev;
}

static storage_queue_ops storage_queue_default_ops = {
    /* init    */ NULL,
    /* cleanup */ storage_queue_default_empty,
    /* enqueue */ storage_queue_default_enqueue,
    /* dequeue */ storage_queue_default_dequeue,
    /* empty   */ storage_queue_default_empty
};


extern void
INT_EVassoc_conversion_action(CManager cm, int stone_id, int stage, 
			      FMFormat target_format, FMFormat incoming_format)
{
    response_cache_element *act;
    stone_type stone;
    int a;
    int id_len;
    FFSTypeHandle format;
    char *server_id;

    stone = stone_struct(cm->evp, stone_id);
    if (!stone) return;

    a = stone->response_cache_count;
    server_id = get_server_ID_FMformat(incoming_format, &id_len);

    if (CMtrace_on(cm, EVerbose)) {
	char *target_tmp = global_name_of_FMFormat(target_format);
	char *incoming_tmp = global_name_of_FMFormat(incoming_format);
	fprintf(cm->CMTrace_file, "Adding Conversion action %d to ", a);
	fprint_stone_identifier(cm->CMTrace_file, cm->evp, stone_id);
	fprintf(cm->CMTrace_file, "\n   Incoming format is %s, target %s\n", incoming_tmp, 
	       target_tmp);
    }
    stone->response_cache = realloc(stone->response_cache,
			     sizeof(stone->response_cache[0]) * (a + 1));
    act = & stone->response_cache[a];
    memset(act, 0, sizeof(*act));
    act->requires_decoded = 0;
    act->action_type = Action_Decode;
    act->reference_format = incoming_format;
    act->stage = Immediate;

    act->o.decode.context = create_FFSContext_FM(cm->evp->fmc);
    format = FFSTypeHandle_from_encode(act->o.decode.context, 
				      server_id);
    act->o.decode.decode_format = format;
    act->o.decode.target_reference_format = target_format;
    establish_conversion(act->o.decode.context, format,
			 format_list_of_FMFormat(act->o.decode.target_reference_format));
    stone->response_cache_count++;
}

int
INT_EVaction_set_output(CManager cm, EVstone stone_num, EVaction act_num, 
		    int output_index, EVstone output_stone)
{
    return INT_EVstone_set_output(cm, stone_num, output_index, output_stone);
}

int
INT_EVstone_set_output(CManager cm, EVstone stone_num, 
		       int output_index, EVstone output_stone)
{
    stone_type stone;
    int output_count = 0;

    stone = stone_struct(cm->evp, stone_num);
    if (!stone) return -1;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Setting output %d on ", output_index);
	fprint_stone_identifier(cm->CMTrace_file, cm->evp, stone_num);
	fprintf(cm->CMTrace_file, " to forward to ");
	fprint_stone_identifier(cm->CMTrace_file, cm->evp, output_stone);
	fprintf(cm->CMTrace_file, "\n");
    }
    output_count = stone->output_count;
    if (output_index >= output_count) {
	stone->output_stone_ids = 
	    realloc(stone->output_stone_ids,
                        sizeof(int) * (output_index + 2));
	for ( ; output_count < output_index; output_count++) {
	    stone->output_stone_ids[output_count] = -1;
	}
	stone->output_count = output_index + 1;
    }
    stone->output_stone_ids[output_index] = output_stone;
    return 1;
}

static int
compatible_stages(int real_stage, int cache_stage) {
    return real_stage == cache_stage || (real_stage == Immediate_and_Multi && cache_stage == Immediate);
}

static action_class cached_stage_for_action(proto_action*);

static int
check_response_cache(CManager cm, stone_type stone, action_class stage, event_item *event)
{
    int i;
    for (i=0; i < stone->response_cache_count; i++) {
//	CMtrace_out(cm, EVerbose, "Response cache %d reference_format is %p (%s), Type %s, stage is %d, requires_decoded is %d\n", i, 
//		    stone->response_cache[i].reference_format, 
//		    global_name_of_FMFormat(stone->response_cache[i].reference_format), action_str[stone->response_cache[i].action_type],
//		    stone->response_cache[i].stage, stone->response_cache[i].requires_decoded);
	if (stone->response_cache[i].reference_format == event->reference_format) {
	    /* 
	     * if the event is encoded and the action requires decoded data,
	     * this action won't do.  Scan further for decode action or 
	     * generate one with response_determination().
	     */
	    if ((stone->response_cache[i].action_type == Action_NoAction) &&
		(stage != stone->response_cache[i].stage)) {
		/* don't return NoAction unless we're querying this exact stage */
		continue;
	    }
	    if (!compatible_stages(stage, stone->response_cache[i].stage)) {
		continue;
	    }
	    if (event->event_encoded && stone->response_cache[i].requires_decoded) {
		continue;
	    }
	    if (!event->event_encoded &&
	    	(stone->response_cache[i].action_type == Action_Decode) &&
	    	(stone->response_cache[i].o.decode.target_reference_format == event->reference_format)) {
	    	continue;
	    }
	    return i;
	} else if (stone->response_cache[i].reference_format == NULL &&
                   !stone->response_cache[i].requires_decoded) {
            return i;
        }
    }
    return -1;
}

static int
determine_action(CManager cm, stone_type stone, action_class stage, event_item *event, int recursed_already)
{
    int return_response;
    if (event->reference_format == NULL) {
	CMtrace_out(cm, EVerbose, "Call to determine_action, event reference_format is NULL\n");
    } else {
	CMtrace_out(cm, EVerbose, "Call to determine_action, event reference_format is %p (%s), stage is %d, encoded is %d\n",
		    event->reference_format, global_name_of_FMFormat(event->reference_format),
		    stage, event->event_encoded);
    }
    return_response = check_response_cache(cm, stone, stage, event);

    if (return_response != -1) return return_response;

    if (response_determination(cm, stone, stage, event) == 1) {
	return_response = check_response_cache(cm, stone, stage, event);
	return return_response;
    }

    /* 
     * there was no action for this event, install a dummy so we 
     * don't search again.
     */
    if (stone->response_cache_count == 0) {
	if (stone->response_cache != NULL) free_response_cache(stone);
	stone->response_cache = malloc(sizeof(stone->response_cache[0]));
    } else {
	stone->response_cache = 
	    realloc(stone->response_cache,
		    (stone->response_cache_count + 1) * sizeof(stone->response_cache[0]));
    }
    stone->response_cache[stone->response_cache_count].reference_format =
	event->reference_format;
    stone->response_cache[stone->response_cache_count].action_type =
	Action_NoAction;
    return_response = stone->response_cache_count++;

    if (stone->default_action != -1 
            && compatible_stages(stage, cached_stage_for_action(&stone->proto_actions[stone->default_action]))
        ) {
/*	    printf(" Returning ");
	    dump_action(stone, NULL, stone->default_action, "   ");*/
	response_cache_element *resp = 
	    &stone->response_cache[return_response];
	proto_action *proto = &stone->proto_actions[stone->default_action];
	resp->proto_action_id = stone->default_action;
	resp->action_type = proto->action_type;
	resp->requires_decoded = proto->data_state;
        resp->stage = stage;
	return return_response;
    }
    stone->response_cache[return_response].action_type = Action_NoAction;
    stone->response_cache[return_response].stage = stage;
    stone->response_cache[return_response].requires_decoded = 0;

    return return_response;
}

/*   GSE
 *   Need to seriously augment buffer handling.  In particular, we have to
 *   formalize the handling of event buffers.  Make sure we have all the
 *   situations covered, try to keep both encoded and decoded versions of
 *   events where possible.  Try to augment testing because many cases are
 *   not covered in homogeneous regression testing (our normal mode).

 *   Possible event situations:
 *  Event_CM_Owned:   This is an event that came in from the network, so CM 
 *     owns the buffer that it lives in.   The event maybe encoded or decoded, 
 *     but whatever the situation, we need to do INT_CMreturn_buffer() to 
 *     free it.   ioBuffer should be NULL;  The 'cm' value tells us where to 
 *     do cmreturn_buffer() .
 *  Event_Freeable:   This is an event that came from a higher-level and 
 *     is provided with a free() routine that we should call when we are 
 *     finished with it.  We are free to return control the the application 
 *     while still retaining the event as long as we later call the free 
 *     routine. 
 *  Event_App_Owned:  This is an event that came from a higher-level and 
 *     is *NOT* provided with a free() routine.  We should not return 
 *     until we are done processing it..
 */

event_item *
get_free_event(event_path_data evp)
{
    event_item *event = malloc(sizeof(*event));
    (void)evp;
    memset(event, 0, sizeof(event_item));
    event->ref_count = 1;
    event->event_len = -1;
    event->ioBuffer = NULL;
    return event;
}

extern void
return_event(event_path_data evp, event_item *event)
{
    (void)evp;
    event->ref_count--;
    if (event->ref_count == 0) {
	/* return event memory */
	switch (event->contents) {
	case Event_CM_Owned:
	    if (event->decoded_event) {
		CMtrace_out(event->cm, CMBufferVerbose, "RETURN decoded event %p\n", event->decoded_event);
		INT_CMreturn_buffer(event->cm, event->decoded_event);
	    } else {
		CMtrace_out(event->cm, CMBufferVerbose, "RETURN encoded event %p\n", event->decoded_event);
		INT_CMreturn_buffer(event->cm, event->encoded_event);
	    }
	    break;
	case Event_Freeable:
	    (event->free_func)(event->decoded_event, event->free_arg);
	    break;
	case Event_App_Owned:
	    if (event->free_func) {
		(event->free_func)(event->free_arg, NULL);
	    }
	    break;
	}
	if (event->attrs != NULL) INT_CMfree_attr_list(event->cm, event->attrs);
	if (event->ioBuffer != NULL)
	    free_FFSBuffer(event->ioBuffer);
	free(event);
    }
}

static event_item *
decode_action(CManager cm, event_item *event, response_cache_element *act)
{
    event_path_data evp = cm->evp;
    if (!event->event_encoded) {
	if (event->reference_format == act->o.decode.target_reference_format) {
	    return event;
	}
	assert(0);
    }
	
    switch(event->contents) {
    case Event_CM_Owned:
	if (FFSdecode_in_place_possible(act->o.decode.decode_format)) {
	    void *decode_buffer;
	    if (!FFSdecode_in_place(act->o.decode.context,
				    event->encoded_event, 
				    (void**) (intptr_t) &decode_buffer)) {
		printf("Decode failed\n");
		return 0;
	    }
	    event->decoded_event = decode_buffer;
	    event->encoded_event = NULL;
	    event->event_encoded = 0;
	    event->reference_format = act->o.decode.target_reference_format;
	    return event;
	} else {
	    size_t decoded_length = FFS_est_decode_length(act->o.decode.context, 
						       event->encoded_event,
						       event->event_len);
	    CMbuffer cm_decode_buf = cm_get_data_buf(cm, decoded_length);
	    void *decode_buffer = cm_decode_buf->buffer;
	    CMtrace_out(event->cm, CMBufferVerbose, "Last cm_get_data_buf was for EVPath decode buffer, return was %p\n", cm_decode_buf);
	    if (event->event_len == -1) printf("BAD LENGTH\n");
	    FFSdecode_to_buffer(act->o.decode.context, event->encoded_event, 
				decode_buffer);
	    event->decoded_event = decode_buffer;
	    event->event_encoded = 0;
	    CMtrace_out(event->cm, CMBufferVerbose, "EVPath now returning original, data is %p\n", event->encoded_event);
	    INT_CMreturn_buffer(cm, event->encoded_event);
	    event->encoded_event = NULL;
	    event->reference_format = act->o.decode.target_reference_format;
	    return event;
	}
    case Event_Freeable:
    case Event_App_Owned:
    {
	/* can't do anything with the old event, make a new one */
	size_t decoded_length = FFS_est_decode_length(act->o.decode.context, 
						   event->encoded_event,
						   event->event_len);
	event_item *new_event = get_free_event(evp);
	CMbuffer cm_decode_buf = cm_get_data_buf(cm, decoded_length);
	void *decode_buffer = cm_decode_buf->buffer;
	CMtrace_out(event->cm, CMBufferVerbose, "Last cm_get_data_buf was for EVPath decode buffer2, return was %p\n", cm_decode_buf);
	if (event->event_len == -1) printf("BAD LENGTH\n");
	FFSdecode_to_buffer(act->o.decode.context, 
			    event->encoded_event, decode_buffer);
	new_event->decoded_event = decode_buffer;
	new_event->event_encoded = 0;
	new_event->encoded_event = NULL;
	new_event->event_len = 0;
	new_event->encoded_eventv = NULL;
	new_event->cm = cm;
	new_event->reference_format = act->o.decode.target_reference_format;
	new_event->contents = Event_CM_Owned;
	if (event->attrs) {
	    new_event->attrs = attr_copy_list(event->attrs);
	} else {
	    new_event->attrs = NULL;
	}
	/* new event will take the place of old event */
	return_event(evp, event);
	return new_event;
    }
    }
    return NULL;   /* shouldn't ever happen */
}

static void
encode_event(CManager cm, event_item *event)
{
    (void)cm;
    if (event->event_encoded) {
	return;
    }
	
    if (event->ioBuffer != NULL) {
	return;  /* already encoded */
    }
    event->ioBuffer = create_FFSBuffer();
    event->encoded_event = 
	FFSencode(event->ioBuffer, event->reference_format,
		  event->decoded_event,
		  (size_t*)&event->event_len);
    event->event_encoded = 1;
}

extern void
cod_encode_event(CManager cm, event_item *event)
{
    encode_event(cm, event);
}

static action_class
cached_stage_for_action(proto_action *act) {
    switch (act->action_type) {
    case Action_Congestion:
        return Congestion;
    case Action_Multi:
        return Immediate_and_Multi;
    case Action_Bridge:
        return Bridge;
    case Action_Terminal:
    case Action_Filter:
    case Action_Split:
    case Action_Immediate:
    case Action_Store:
    case Action_Thread_Bridge:
    case Action_NoAction:
        return Immediate;
    default:
	printf("Action_type in cached_stage_for_action is %d\n", act->action_type);
        assert(0);    
    }
    /*NOTREACHED*/
    return Immediate;
}

extern event_item * 
cod_decode_event(CManager cm, int stone_num, int act_num, event_item *event) {
    stone_type stone;
    action_class stage;
    int resp_id;

    assert(!event->decoded_event);

    stone = stone_struct(cm->evp, stone_num);
    stage = cached_stage_for_action(&stone->proto_actions[act_num]);

    resp_id = determine_action(cm, stone, stage, event, 0);
    if (stone->response_cache[resp_id].action_type != Action_Decode) {
	resp_id = determine_action(cm, stone, Immediate, event, 0);
    }
    if (stone->response_cache[resp_id].action_type != Action_Decode) {
	char *tmp;
	printf("Warning!  bad multiq action found for incoming an event on stone %x, stage %d\n",
	       stone->local_id, stage);
	printf("A decode response should be installed into the response cache for event type \"%s\" (%p)\n", tmp = global_name_of_FMFormat(event->reference_format), event->reference_format);
        free(tmp);
	dump_stone(stone);
    }
    return decode_action(cm, event, &stone->response_cache[resp_id]);
}

static void
fdump_proto_action(FILE *out, stone_type stone, int a, const char *indent)
{
    proto_action *proto = &stone->proto_actions[a];
    (void)indent;
    fprintf(out, " Proto-Action %d - %s\n", a, action_str[proto->action_type]);
}

static void
fdump_action(FILE* out, stone_type stone, response_cache_element *resp, int a, const char *indent)
{
    proto_action *act;
    (void)indent;
    if ((resp != NULL) && (resp->action_type == Action_NoAction)) {
	fprintf(out, "NO ACTION REGISTERED\n");
	return;
    }
    act = &stone->proto_actions[a];
    fprintf(out, " Action %d - %s  ", a, action_str[act->action_type]);
    if (act->data_state == Accepts_All) {
	fprintf(out, "accepts any encode state\n");
    } else if (act->data_state == Requires_Decoded) {
	fprintf(out, "requires decoded\n");
    } else if (act->data_state == Requires_Contig_Encoded) {
	fprintf(out, "requires contiguous encoded\n");
    } else if (act->data_state == Requires_Vector_Encoded) {
	fprintf(out, "requires vector encoded\n");
    }
    fprintf(out, "  expects formats ");
    if (act->matching_reference_formats) {
	int i = 0;
	while (act->matching_reference_formats[i] != NULL) {
	    char *tmp;
	    fprintf(out, "\"%s\" (%p), ", tmp = global_name_of_FMFormat(act->matching_reference_formats[i]), act->matching_reference_formats[i]);
	    i++;
	    free(tmp);
	}
    } else {
	fprintf(out, " NULL");
    }
    fprintf(out, "\n");
    switch(act->action_type) {
    case Action_Bridge:
	fprintf(out, "  Target: %s: connection %p, remote_stone_id %d\n",
	       (act->o.bri.remote_path ? act->o.bri.remote_path : "NULL" ),
	       (void*)act->o.bri.conn, act->o.bri.remote_stone_id);
	if (act->o.bri.conn != NULL) fdump_attr_list(out, act->o.bri.conn->attrs);
	if (act->o.bri.conn_failed) fprintf(out, "Connection has FAILED!\n");
	break;
    case Action_Thread_Bridge:
	fprintf(out, "  Target: CManager %p, stone_id %d\n",
	       act->o.thr_bri.target_cm, act->o.thr_bri.target_stone_id);
	if (act->o.thr_bri.target_cm_shutdown) fprintf(out, "TARGET CM HAS SHUTDOWN!\n");
	break;
    case Action_Terminal:
	break;
    case Action_Filter:
/*	printf("  Filter proto action number %d\n",
	act->proto_action_id);*/
	break;
    case Action_Decode:
	fprintf(out, "   Decoding action\n");
	break;
    case Action_Split:
	fprintf(out, "    Split action\n");
	break;
    case Action_Immediate: 
	fprintf(out, "   Immediate action\n");
	dump_mrd(act->o.imm.mutable_response_data);
	break;
    case Action_Store:
        fprintf(out, "   Store action: %d/%d items\n", act->o.store.num_stored, 
                    act->o.store.max_stored);
    case Action_NoAction:
	fprintf(out, "   NoAction\n");
	break;
    case Action_Multi:
	fprintf(out, "   Multi action\n");
	dump_mrd(act->o.imm.mutable_response_data);
	break;
    default:
	assert(FALSE);
    }
}

static void
dump_stone(stone_type stone)
{
    fdump_stone(stdout, stone);
}

static void
fdump_stone(FILE* out, stone_type stone)
{
    int i;
    fprintf(out, "Dump stone ID %d, local addr %p, default action %d\n",
	    stone->local_id, stone, stone->default_action);
    fprintf(out, "       Target Stones:");
    {
	int i;
	for(i = 0; i < stone->output_count; i++) {
	    if (i != stone->output_count - 1) {
		fprintf(out, " %d,", stone->output_stone_ids[i]);
	    } else {
		fprintf(out, " %d\n", stone->output_stone_ids[i]);
	    }
	}
    }
    fprintf(out, "  proto_action_count %d:\n", stone->proto_action_count);
    for (i=0; i< stone->proto_action_count; i++) {
	fdump_proto_action(out, stone, i, "    ");
    }
    fprintf(out, "  proto_action_count %d:\n", stone->proto_action_count);
    for (i=0; i< stone->proto_action_count; i++) {
	fdump_action(out, stone, NULL, i, "    ");
    }
    fprintf(out, "  response_cache_count %d:\n", stone->response_cache_count);
    for (i=0; i< stone->response_cache_count; i++) {
	response_cache_element *resp = &stone->response_cache[i];
	fprintf(out, "Response cache item %d, reference format %p (%s)\n", i, resp->reference_format,
		resp->reference_format ? global_name_of_FMFormat(resp->reference_format) : "<none>");
	fprintf(out, "stage %d, action_type %s, proto_action_id %d, requires_decoded %d\n", resp->stage,
	       action_str[resp->action_type], resp->proto_action_id, resp->requires_decoded);
    }
}

/*static void
dump_action(stone_type stone, response_cache_element *resp, int a, const char *indent)
{
    fdump_action(stdout, stone, resp, a, indent);
}
*/
void
EVdump_stone(CManager cm,  EVstone stone_num)
{
    stone_type stone = stone_struct(cm->evp, stone_num);
    dump_stone(stone);
}

int
internal_path_submit(CManager cm, int local_path_id, event_item *event)
{
    event_item *event_to_submit = event;

    assert(CManager_locked(cm));
    enqueue_event(cm, local_path_id, event_to_submit);
    return 1;
}

static void
update_event_length_sum(CManager cm, proto_action *act, event_item *event)
{
    int eventlength;
    int totallength; 
    static atom_t CM_EVENT_SIZE = -1;
    static atom_t EV_EVENT_LSUM = -1;

    if (CM_EVENT_SIZE == -1) {
	CM_EVENT_SIZE = attr_atom_from_string("CM_EVENT_SIZE");
	EV_EVENT_LSUM = attr_atom_from_string("EV_EVENT_LSUM");
    }
    /*update act->event_length_sum:*/
    if (get_int_attr(event->attrs, CM_EVENT_SIZE, & eventlength)) {
	if (eventlength >= 0 )
	    act->event_length_sum += eventlength; 
	else 
	    act->event_length_sum = -1; /*received event with undecided size, invalidate length_sum*/
    } else {
	/*attr CM_EVENT_SIZE doesn't exist. 
	  two possibilities: 1. something broken. 2. event is 
	  from application (via EVSubmit, so doesn't have a valid 
	  CM_EVENT_SIZE attr (This is not implemented unless it is 
	  really required by someone.)
	*/
	return;
    } 
    /* also update the EV_EVENT_LSUM attrs */
    if(act->attrs == NULL){
	act->attrs = CMcreate_attr_list(cm);
    }
    totallength = (int) act->event_length_sum;/*1024*/ 
    set_int_attr(act->attrs, EV_EVENT_LSUM, totallength);
}

static int
is_immediate_action(response_cache_element *act)
{
    switch(act->action_type) {
    case Action_Terminal:
    case Action_Filter:
    case Action_Split:
    case Action_Immediate:
    case Action_Store:
    case Action_Thread_Bridge:
	return 1;
    default:
	return 0;
    }
}

static int
is_multi_action(response_cache_element *act)
{
    switch(act->action_type) {
    case Action_Multi:
	return 1;
    default:
	return 0;
    }
}

static int
is_congestion_action(response_cache_element *act) {
    return act->action_type == Action_Congestion;
}

static int do_bridge_action(CManager cm, int s);

static void backpressure_set(CManager, EVstone, int stalledp);
static int process_stone_pending_output(CManager, EVstone);

static void
push_activation_record_on_stack(CManager cm, ev_handler_activation_ptr rec)
{
    event_path_data evp = cm->evp;
    rec->thread_id = thr_thread_self();
    if (evp->activation_stack != NULL) {
	evp->activation_stack->prev = rec;
    }
    rec->prev = NULL;
    rec->next = evp->activation_stack;
    evp->activation_stack = rec;
}

static void
pop_activation_record_from_stack(CManager cm, ev_handler_activation_ptr rec)
{
    event_path_data evp = cm->evp;
    ev_handler_activation_ptr tmp;
    thr_thread_t self = thr_thread_self();
    if (!evp->activation_stack) {
	printf("Activation stack inconsistency!  No records!\n"); 
	return;
    }
    if (evp->activation_stack->thread_id == self) {
	evp->activation_stack = evp->activation_stack->next;
	if (evp->activation_stack) {
	    evp->activation_stack->prev = NULL;
	}
	return;
    }
    tmp = evp->activation_stack->next;
    while(tmp) {
	if (tmp->thread_id == self) {
	    tmp->prev->next = tmp->next;
	    if (tmp->next) {
		tmp->next->prev = tmp->prev;
	    }
	    return;
	}
	tmp = tmp->next;
    }
    printf("Activation stack inconsistency!  Record with thread ID now found!\n"); 
}

static ev_handler_activation_ptr
find_activation_record_from_stack(CManager cm)
{
    event_path_data evp = cm->evp;
    ev_handler_activation_ptr tmp;
    thr_thread_t self = thr_thread_self();
    tmp = evp->activation_stack;
    while(tmp) {
	if (tmp->thread_id == self) {
	    return tmp;
	}
	tmp = tmp->next;
    }
    return NULL;
}

extern EVstone
INT_EVexecuting_stone(CManager cm)
{
    ev_handler_activation_ptr rec = find_activation_record_from_stack(cm);
    if (!rec) return -1;
    return rec->stone_id;
}

static event_item *
reassign_memory_event(CManager cm, event_item *event, int do_decode);

static int
process_events_stone(CManager cm, int s, action_class c)
{

/*  Process all possible events on this stone that match the class */

    event_path_data evp = cm->evp;
    stone_type stone = NULL;
    int more_pending = 0;
    int did_something = 0;
    queue_item *item;
    action_state as = evp->as;


    if (s == -1) return 0;
    if (as->last_active_stone == s) as->last_active_stone = -1;
    stone = stone_struct(cm->evp, s);
    if (!stone) return 0;
    if (stone->queue_size == 0) return 0;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Considering events on ");
	fprint_stone_identifier(cm->CMTrace_file, evp, s);
	fprintf(cm->CMTrace_file, "\n");
    }
    if (stone->local_id == -1) return 0;
    if (stone->is_draining == 2) return 0;
    if (stone->is_frozen == 1) return 0;
    if (c == Immediate_and_Multi && stone->pending_output) {
        more_pending += process_stone_pending_output(cm, s);
    }
    if (is_bridge_stone(cm, s) && (c != Bridge) && (c != Congestion)) return 0;
    
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Process events ");
	fprint_stone_identifier(cm->CMTrace_file, evp, s);
	fprintf(cm->CMTrace_file, "\n");
    }
    item = stone->queue->queue_head;
    if (is_bridge_stone(cm, s) && (c == Bridge)) {
	stone->is_processing = 1;
	do_bridge_action(cm, s);
	stone->is_processing = 0;
	return 0;
    }
    while (item != NULL && stone->is_draining == 0 && stone->is_frozen == 0) {
	queue_item *next = item->next;
	response_cache_element *resp;
	response_cache_element *act = NULL;
	if (item->handled) {
	    item = next;
	    continue;
	} else {
	    /* determine what kind of action to take here */
	    int resp_id;
	    event_item *event = item->item;
	    resp_id = determine_action(cm, stone, c, item->item, 0);
            assert(resp_id < stone->response_cache_count);
	    if (stone->response_cache[resp_id].action_type == Action_NoAction
                && c == Immediate_and_Multi) {
                /* ignore event */
                char *tmp = NULL;
                if (event->reference_format)
                    tmp = global_name_of_FMFormat(event->reference_format);
                printf("No action found for event %p submitted to \n", event);
		print_stone_identifier(evp, s);
		printf("\n");
                dump_stone(stone_struct(evp, s));
                if (tmp != NULL) {
                    static int first = 1;
                    printf("    Unhandled incoming event format was \"%s\"\n", tmp);
                    if (first) {
                        first = 0;
                        printf("\n\t** use \"format_info <format_name>\" to get full format information ** \n\n");
                    }
                } else {
                    printf("    Unhandled incoming event format was NULL\n");
                }
                if (tmp) free(tmp);
                event = dequeue_item(cm, stone, item);
                return_event(evp, event);
            }
	    resp = &stone->response_cache[resp_id];
	    if (resp->action_type == Action_Decode) {
		event_item *event_to_submit;
		if (!event->event_encoded) {
		    event_item *old_data_event;
		    CMtrace_out(cm, EVerbose, "Encoding event prior to decode for conversion, action id %d\n", resp_id);
		    old_data_event = reassign_memory_event(cm, event, 0);  /* reassign memory */
		    return_event(evp, old_data_event);
		    event->ref_count++;
		}
		CMtrace_out(cm, EVerbose, "Decoding event, action id %d\n", resp_id);
		event_to_submit = decode_action(cm, event, resp);
		if (event_to_submit == NULL) return more_pending;
		item->item = event_to_submit;
		resp_id = determine_action(cm, stone, c, event_to_submit, 0);
		resp = &stone->response_cache[resp_id];
	    }
	    if (CMtrace_on(cm, EVerbose)) {
		fprintf(cm->CMTrace_file, "next action event %p on ", event);
		fprint_stone_identifier(cm->CMTrace_file, evp, s);
		fprintf(cm->CMTrace_file, " action type is %s, reference_format is %p (%s), stage is %d, requires_decoded is %d\n",
			action_str[resp->action_type], resp->reference_format, 
			resp->reference_format ? global_name_of_FMFormat(resp->reference_format) : "<none>",
			resp->stage, resp->requires_decoded);
		fdump_action(cm->CMTrace_file, stone, resp, resp->proto_action_id, "    ");
	    }
            act = &stone->response_cache[resp_id];
	}
        backpressure_check(cm, s);
        if (!compatible_stages(c, act->stage)) {
	    /* do nothing */
        } else if (is_immediate_action(act)) {
	    event_item *event = dequeue_item(cm, stone, item);
	    ev_handler_activation_rec act_rec;
	    switch(act->action_type) {
	    case Action_Terminal:
	    case Action_Filter: {
		/* the data should already be in the right format */
		int proto = act->proto_action_id;
		proto_action *p = &stone->proto_actions[proto];
		int out;
		struct terminal_proto_vals *term = &(p->o.term);
		EVSimpleHandlerFunc handler = term->handler;
		void *client_data = term->client_data;
		CMtrace_out(cm, EVerbose, "Executing terminal/filter event\n");
		update_event_length_sum(cm, p, event);
		{
		    queue_item *new = malloc(sizeof(queue_item));
		    new->item = event;
		    new->next = cm->evp->current_event_list;
		    cm->evp->current_event_list = new;
		}
		did_something=1;
		stone->is_processing = 1;
		if ((p->data_state == Requires_Contig_Encoded) && 
		    (event->event_encoded == 0)) {
		    encode_event(cm, event);
		}
		act_rec.stone_id = stone->local_id;
		push_activation_record_on_stack(cm, &act_rec);
		CManager_unlock(cm);
		if (event->event_encoded == 0) {
		    out = (handler)(cm, event->decoded_event, client_data,
				    event->attrs);
		} else {
		    EVRawHandlerFunc han = (EVRawHandlerFunc)term->handler;
		    out = (han)(cm, event->encoded_event, event->event_len,
				    client_data, event->attrs);
		}
		CManager_lock(cm);
		stone = stone_struct(cm->evp, s);
		pop_activation_record_from_stack(cm, &act_rec);
		stone->is_processing = 0;
		{
		    if (cm->evp->current_event_list->item == event) {
			queue_item *tmp = cm->evp->current_event_list;
			cm->evp->current_event_list = tmp->next;
			free(tmp);
		    } else {
			queue_item *tmp = cm->evp->current_event_list;
			queue_item *last;
			while(tmp->item != event) {
			    last = tmp;
			    tmp = tmp->next;
			    if (tmp == NULL) break;
			}
			if (tmp == NULL) {
			    printf("Failed to dequeue item from executing events list!\n");
			} else {
			    last->next = tmp->next;
			    free(tmp);
			}
		    }
		}

		if (act->action_type == Action_Filter) {
		    if (out) {
			CMtrace_out(cm, EVerbose, "Filter passed event to stone %x, submitting\n", term->target_stone_id);
			internal_path_submit(cm, 
					     term->target_stone_id,
					     event);
			more_pending++;
		    } else {
			CMtrace_out(cm, EVerbose, "Filter discarded event\n");
		    }			    
		} else {
		    CMtrace_out(cm, EVerbose, "Finish terminal event\n");
		}
		return_event(evp, event);
		break;
	    }
	    case Action_Split: {
		int t = 0;
		proto_action *p = &stone->proto_actions[act->proto_action_id];
		update_event_length_sum(cm, p, event);
		for (t=0; t < stone->output_count; t++) {
		    if (stone->output_stone_ids[t] != -1) {
			internal_path_submit(cm, 
					     stone->output_stone_ids[t],
					     event);
			more_pending++;
		    }
		}
		return_event(evp, event);
		break;
	    }
	    case Action_Immediate: {
		EVImmediateHandlerFunc func;
		void *client_data;
		int *out_stones;
		int out_count;
		int in_play = as->events_in_play;
		/* data is already in the right format */
		func = act->o.imm.handler;
		client_data = act->o.imm.client_data;
		out_stones = stone->output_stone_ids;
		out_count = stone->output_count;
		did_something = 1;
		stone->is_processing = 1;
		func(cm, event, client_data, event->attrs, out_count, out_stones);
		stone = stone_struct(cm->evp, s);
		stone->is_processing = 0;
		return_event(evp, event);
		if (as->events_in_play > in_play)
		    more_pending++;   /* maybe??? */
		break;
	    }
            case Action_Store: {
		proto_action *p = &stone->proto_actions[act->proto_action_id];
                storage_queue_enqueue(cm, &p->o.store.queue, event);
                CMtrace_out(cm, EVerbose, "Enqueued item to store\n");
                if (++p->o.store.num_stored > p->o.store.max_stored
                        && p->o.store.max_stored != -1) {
                    CMtrace_out(cm, EVerbose, "Dequeuing item because of limit\n");
                    event_item *last_event;
                    last_event = storage_queue_dequeue(cm, &p->o.store.queue);
                    assert(last_event->ref_count > 0);
                    --p->o.store.num_stored;
                    internal_path_submit(cm, p->o.store.target_stone_id, last_event);
                    return_event(evp, last_event);
                    more_pending++;
                }
                break;
            }
	    case Action_Thread_Bridge: {
		proto_action *p = &stone->proto_actions[act->proto_action_id];
		if (p->o.thr_bri.target_cm_shutdown) continue;
                thread_bridge_transfer(cm, event, p->o.thr_bri.target_cm,
				       p->o.thr_bri.target_stone_id);
		break;
	    }
	    case Action_Bridge:
	    default:
		assert(FALSE);
	    }
	} else if (is_multi_action(act) || is_congestion_action(act)) {
	    /* event_item *event = dequeue_item(cm, stone->queue, item); XXX */
	    if (stone->new_enqueue_flag) {
		if (!item->handled) {
		    did_something = 1;
		    stone->is_processing = 1;
		    item->handled = 1;
		    if ((act->o.multi.handler)(cm, stone->queue, item,
					       act->o.multi.client_data, stone->output_count,
					       stone->output_stone_ids))
			more_pending++;
		    stone = stone_struct(cm->evp, s);
		    stone->is_processing = 0;
		}
	    }
            break;    
	} 
	item = next;
    }
    if (did_something) {
	stone->new_enqueue_flag = 0;
    }
    return more_pending;
}
	

static
int
process_local_actions(CManager cm)
{
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    int s, more_pending = 0;
/*    CMtrace_out(cm, EVerbose, "Process local actions\n");*/
    if (as == NULL) {
	as = evp->as = malloc(sizeof(*as));
	memset(as, 0, sizeof(*as));
	as->last_active_stone = -1;
    }
 restart:
    if (as->last_active_stone != -1) {
	more_pending = 1;
	while (more_pending) {
	    CMtrace_out(cm, EVerbose, "Process local actions on stone %x\n",
			as->last_active_stone);
	    
	    CMtrace_out(cm, EVerbose, "0 - in-play %d\n", as->events_in_play);
	    more_pending = process_events_stone(cm, as->last_active_stone, Immediate);
	}
    }
    if (as->events_in_play > 0) {
	/* check all stones */
	for (s = evp->stone_base_num; s < evp->stone_count + evp->stone_base_num; s++) {
	    stone_type stone = stone_struct(evp, s);
	    if (!stone) continue;
	    if (stone->local_id == -1) continue;
	    if (stone->is_draining == 2) continue;
	    if (stone->is_frozen == 1) continue;
	    CMtrace_out(cm, EVerbose, "1 - in-play %d\n", as->events_in_play);
	    more_pending += process_events_stone(cm, s, Immediate_and_Multi);
	    if (more_pending && (as->last_active_stone != -1)) goto restart;
	}
    }
    if (as->last_active_stone != -1) {
	CMtrace_out(cm, EVerbose, "Process output actions on stone %x\n",
		    as->last_active_stone);
	CMtrace_out(cm, EVerbose, "2 - in-play %d\n", as->events_in_play);
	more_pending += process_events_stone(cm, as->last_active_stone, Bridge);
    }
    if (as->events_in_play > 0) {
	/* check all stones */
	for (s = evp->stone_base_num; s < evp->stone_count + evp->stone_base_num; s++) {
	    stone_type stone = stone_struct(evp, s);
	    if (!stone) continue;
	    if (stone->local_id == -1) continue;
	    if (stone->is_frozen == 1) continue;
	    CMtrace_out(cm, EVerbose, "3 - in-play %d\n", as->events_in_play);
	    more_pending += process_events_stone(cm, s, Bridge);
	}
    }

    return more_pending;
}

void
INT_EVforget_connection(CManager cm, CMConnection conn)
{
    event_path_data evp = cm->evp;
    int s;
    for (s = evp->stone_base_num; s < evp->stone_count + evp->stone_base_num; ++s)  {
	stone_type stone = stone_struct(evp, s);
	if (!stone) continue;
        if (stone->last_remote_source == conn) {
            stone->last_remote_source = NULL;
            stone->squelch_depth = 0;
        }
    }
}

/*
 *  this handler gets called on a CMConnection write failure, indicating that an event couldn't be sent to a particular remote stone.
 *  Its job is to determine what stone and to call the upper-level stone close handler (probably EV_DFG)
 */
static void
stone_close_handler(CManager cm, CMConnection conn, void *client_data)
{
    event_path_data evp = cm->evp;
    int s = (int)(intptr_t)client_data;  /* stone ID */
    int a = 0;
    stone_type stone;
    EVStoneCloseHandlerFunc handler = NULL;
    CManager_lock(cm);
    stone = stone_struct(evp, s);
    if (!stone) {
	CMtrace_out(cm, EVerbose, "Got a close for connection %p on already free'd stone %x, shutting down\n",
		    conn, s);
	CManager_unlock(cm);
	return;
    }
    CMtrace_out(cm, EVerbose, "Got a close for connection %p on stone %x, shutting down\n",
		conn, s);
    for (a=0 ; a < stone->proto_action_count; a++) {
	proto_action *act = &stone->proto_actions[a];
	if ((act->action_type == Action_Bridge) && 
	    (act->o.bri.conn == conn)) {
	    act->o.bri.conn_failed = 1;
	    act->o.bri.conn = NULL;
	    CMtrace_out(cm, CMFreeVerbose, "Closing and dereferencing conn %p\n", conn);
	    INT_CMConnection_close(conn);   /* dereference the connection */
/*	    while (act->queue->queue_head != NULL) {
		int action_id;
		event_item *event = dequeue_event(cm, act->queue, &action_id);
		return_event(evp, event);
		}*/
	    if (evp->app_stone_close_handler) {
	        handler = evp->app_stone_close_handler;
	    }
	}
    }
    CManager_unlock(cm);
    if (handler) handler(cm, conn, s, evp->app_stone_close_data);
}

#ifdef NOT_DEF
static void
write_callback_handler(CManager cm, CMConnection conn, void *client_data)
{
    event_path_data evp = cm->evp;
    int s = (int)(long) client_data;
    stone_type stone = stone_struct(evp, s);
    CMtrace_out(cm, EVerbose, "In Write callback, write_pending is %d, stone is %d\n", conn->write_pending, s);
    if (conn->write_pending) {
	/* nothing for EVPath level to do yet */
	return;
    }
    assert(CManager_locked(conn->cm));
    assert(stone->write_callback != -1);
    INT_CMunregister_write_callback(conn, stone->write_callback);
    stone->write_callback = -1;
    /* try calling the congestion handler before we write again... */
    process_events_stone(cm, s, Congestion);
    do_bridge_action(cm, s);
}
#endif

/* this is called to make a store-send happen; it just
 * runs the standard loop later.
 */
static void
deferred_process_actions(CManager cm, void *client_data)
{
    (void)client_data;
    CManager_lock(cm);
    if (cm->evp) cm->evp->delay_task_pending = 0;
    while (cm->evp && process_local_actions(cm));
    CManager_unlock(cm);
}

static
int
do_bridge_action(CManager cm, int s)
{
    event_path_data evp = cm->evp;
    proto_action *act = NULL;
    stone_type stone;
    int a;
    stone = stone_struct(evp, s);
    CMtrace_out(cm, EVerbose, "Process output action on stone %x, frozen %d draining %d outputting %d, in_get_conn %d\n", s, stone->is_frozen, stone->is_draining, stone->is_outputting, evp->in_get_conn);

    if (stone->is_frozen || (stone->is_draining == 2)) return 0;
    if (stone->is_outputting) return 0;
    /* if we're sitting in get_conn, don't proceed to make sure we don't try get_conn again */
    if (evp->in_get_conn) {
        if (!evp->delay_task_pending) {
            evp->delay_task_pending = 1;
            (void) INT_CMadd_delayed_task(cm, 0, 100, deferred_process_actions, NULL);
        }
        return 0;
    }
    stone->is_outputting = 1;
    for (a=0 ; a < stone->proto_action_count && stone->is_frozen == 0 && (stone->is_draining != 2); a++) {
	if (stone->proto_actions[a].action_type == Action_Bridge) {
	    act = &stone->proto_actions[a];
	}
    }
    if (act->o.bri.conn_failed) return 0;
    if (act->o.bri.conn == NULL) {
        attr_list contact_list = act->o.bri.remote_contact;
        CMConnection conn;
        evp->in_get_conn++;
        conn = INT_CMget_conn(cm, contact_list);
        evp->in_get_conn--;
	if (act->o.bri.conn != NULL) {
	    /*
	     *  INT_CMget_conn() isn't synchronous.  
	     *  Someone else might do the same CMget_conn() too.  That isn't 
	     *  a problem, *unless* that someone else is this routine and 
	     *  they have done this and stored the result since we checked 
	     *  for o.bri.conn was NULL.  If they have, and it's now not NULL,
	     *  don't overwrite that value, just use it, but kill the extra 
	     *  reference count.  
	     */
	    INT_CMConnection_dereference(conn);
	    conn = act->o.bri.conn;
	} else {
	    act->o.bri.conn = conn;
	}
        if (conn == NULL) {
            if (CMtrace_on(cm, EVWarning)) {
                fprintf(cm->CMTrace_file, "EVassoc_bridge_action - failed to contact host at contact point \n\t");
                if (contact_list != NULL) {
                    fdump_attr_list(cm->CMTrace_file, contact_list);
                } else {
                    fprintf(cm->CMTrace_file, "NULL\n");
                }
                fprintf(cm->CMTrace_file, "Bridge action association failed for stone %x, outputting to remote stone %x\n",
                        s, act->o.bri.remote_stone_id);
                act->o.bri.conn_failed = 1;
            }
            return -1;
        }
        INT_CMconn_register_close_handler(conn, stone_close_handler, 
                                          (void*)(intptr_t)s);
    }
    while (stone->queue->queue_head != NULL) {
	int ret = 1;
	if (act->o.bri.conn && 
	    INT_CMConnection_write_would_block(act->o.bri.conn)) {
/*            queue_item *q = stone->queue->queue_head;*/
	    {
		/* this is temporary, a disabling of congestion handlers */
		INT_CMConnection_wait_for_pending_write(act->o.bri.conn);
		if (stone->queue->queue_head == NULL) {
		    /* all the events disappeared while we were waiting */
		    return 0;
		}
	    }
		

/*	    CMtrace_out(cm, EVerbose, "Would call congestion_handler, new flag %d\n", stone->new_enqueue_flag);
	    if (stone->new_enqueue_flag == 1) {
		stone->new_enqueue_flag = 0;
                while (q != NULL) {q = q->next; i++;}
		CMtrace_out(cm, EVerbose, "Would call congestion_handler, %d items queued\n", i);
		if (stone->write_callback == -1) {
		    stone->write_callback = 
                        INT_CMregister_write_callback(act->o.bri.conn, 
						  write_callback_handler, 
						  (void*)(long)s);
		}
                process_events_stone(cm, s, Congestion);
		return 0;
	    }*/
	}
	event_item *event = dequeue_event(cm, stone);
	size_t event_length = 0;
	if (act->o.bri.conn == NULL) {
	    CMtrace_out(cm, EVerbose, "Bridge stone %x has closed connection\n", s);
	} else {
	    CMtrace_out(cm, EVerbose, "Writing event to remote stone %x\n",
			act->o.bri.remote_stone_id);
	    if (event->format) {
		ret = internal_write_event(act->o.bri.conn, event->format,
					   &act->o.bri.remote_stone_id, 4, 
					   event, event->attrs, &event_length);
		if (ret == 0) {
		  printf("DETECTED FAILED EVENT WRITE\n");
		}
	    } else {
		struct _CMFormat tmp_format;
		if (event->reference_format == NULL) {
		    CMtrace_out(cm, EVWarning, "Tried to output event with NULL reference format.  Event discarded.\n");
		} else {
		    tmp_format.fmformat = event->reference_format;
		    tmp_format.format_name = name_of_FMformat(event->reference_format);
		    tmp_format.registration_pending = 0;
		    ret = internal_write_event(act->o.bri.conn, &tmp_format,
					       &act->o.bri.remote_stone_id, 4, 
					       event, event->attrs, &event_length);
		}
	    }
	}
	return_event(evp, event);
        backpressure_check(cm, s);
	if (ret == 0) {
	    if (CMtrace_on(cm, EVWarning)) {
		fprintf(cm->CMTrace_file, "Warning!  Write failed for output action %d on stone %x, event likely not transmitted\n", a, s);
		fprintf(cm->CMTrace_file, "   -  Bridge stone ");
		fprint_stone_identifier(cm->CMTrace_file, evp, s);
		fprintf(cm->CMTrace_file, " disabled\n");
	    }
	    if (act->o.bri.conn != NULL)  {
		CMtrace_out(cm, CMFreeVerbose, "Closing and dereferencing conn after write failure %p\n", act->o.bri.conn);
                INT_CMConnection_failed(act->o.bri.conn);
		INT_CMConnection_close(act->o.bri.conn);
	    }
	    act->o.bri.conn_failed = 1;
	    act->o.bri.conn = NULL;
	} else {
	    static atom_t EV_EVENT_COUNT = -1;
	    static atom_t EV_EVENT_LSUM = -1;
	    ssize_t length_sum = 0;
	    int event_count = 0;
	    if (EV_EVENT_COUNT == -1) {
		EV_EVENT_COUNT = attr_atom_from_string("EV_EVENT_COUNT");
		EV_EVENT_LSUM = attr_atom_from_string("EV_EVENT_LSUM");
	    }
	    if (!stone->stone_attrs) {
		stone->stone_attrs = create_attr_list();
	    } else {
		get_int_attr(stone->stone_attrs, EV_EVENT_COUNT, &event_count);
		get_long_attr(stone->stone_attrs, EV_EVENT_LSUM, &length_sum);
	    }
	    event_count++;
	    length_sum += event_length;
	    set_int_attr(stone->stone_attrs, EV_EVENT_COUNT, event_count);
	    set_long_attr(stone->stone_attrs, EV_EVENT_LSUM, length_sum);
	}
    }
    stone->is_outputting = 0;
    return 0;
}

void
INT_EVset_store_limit(CManager cm, EVstone stone_num, EVaction action_num, int new_limit)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    proto_action *p;

    stone = stone_struct(evp, stone_num);
    if (!stone) return;

    p = &stone->proto_actions[action_num];
    p->o.store.max_stored = new_limit;
    if (p->o.store.max_stored != -1) {
        while (p->o.store.num_stored > p->o.store.max_stored) {
            event_item *item;
            item = storage_queue_dequeue(cm, &p->o.store.queue);
            if (!item)
                break;
            p->o.store.num_stored--;
            internal_path_submit(cm, p->o.store.target_stone_id, item);
            while (process_local_actions(cm));
            return_event(evp, item);
        }
    }
}

extern EVaction
INT_EVassoc_mutated_imm_action(CManager cm, EVstone stone_id, EVaction act_num,
			       EVImmediateHandlerFunc func, void *client_data, 
			       FMFormat reference_format, int_free_func free_func)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int resp_num;

    stone = stone_struct(evp, stone_id);
    if (!stone) return -1;

    resp_num = stone->response_cache_count;
    stone->response_cache = realloc(stone->response_cache, sizeof(stone->response_cache[0]) * (resp_num + 1));
    response_cache_element *resp = &stone->response_cache[stone->response_cache_count];
    resp->action_type = Action_Immediate;
    resp->requires_decoded = 1;
    resp->proto_action_id = act_num;
    resp->o.imm.handler = func;
    resp->o.imm.client_data = client_data;
    resp->o.imm.free_func = free_func;
    resp->reference_format = reference_format;
    resp->stage = cached_stage_for_action(&stone->proto_actions[act_num]);
    stone->response_cache_count++;
    return resp_num;
}

extern EVaction
INT_EVassoc_anon_multi_action(CManager cm, EVstone stone_id, EVaction act_num,
			       EVMultiHandlerFunc func, void *client_data, FMFormat anon_target)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    stone_type stone;
    int resp_num;

    stone = stone_struct(evp, stone_id);
    resp_num = stone->response_cache_count;
    stone->response_cache = realloc(stone->response_cache, sizeof(stone->response_cache[0]) * (resp_num + 1));
    response_cache_element *resp;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Installing anon action response for multi action %d on ", act_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_id);
	fprintf(cm->CMTrace_file, "\n");
    }
    resp = &stone->response_cache[stone->response_cache_count];
    resp->action_type = stone->proto_actions[act_num].action_type;
    resp->requires_decoded = 0;
    resp->proto_action_id = act_num;
    resp->o.multi.handler = func;
    resp->o.multi.client_data = client_data;
    resp->o.multi.free_func = NULL;
    resp->stage = cached_stage_for_action(&stone->proto_actions[act_num]);
    resp->reference_format = anon_target;
    if (CMtrace_on(cm, EVerbose)) {
	char *tmp;
	if (resp->reference_format) {
	    tmp = global_name_of_FMFormat(resp->reference_format);
	} else {
	    tmp = strdup("<none>");
	}
	fprintf(cm->CMTrace_file, "\tResponse %d for format \"%s\"(%p)", stone->response_cache_count, tmp, resp->reference_format);
	free(tmp);
    }
    stone->response_cache_count += 1;
    fix_response_cache(stone);
    return resp_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

extern EVaction
INT_EVassoc_mutated_multi_action(CManager cm, EVstone stone_id, EVaction act_num,
				  EVMultiHandlerFunc func, void *client_data, 
				 FMFormat *reference_formats, int_free_func free_func)
{
#ifdef HAVE_COD_H
    event_path_data evp = cm->evp;
    stone_type stone;
    int resp_num;
    int queue_count = 0, i;

    stone = stone_struct(evp, stone_id);
    resp_num = stone->response_cache_count;
    while (reference_formats[queue_count] != NULL) queue_count++;
    stone->response_cache = realloc(stone->response_cache, sizeof(stone->response_cache[0]) * (resp_num + queue_count));
    response_cache_element *resp;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Installing %d mutated action responses for multi action %d on ",
		queue_count, act_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_id);
	fprintf(cm->CMTrace_file, "\n");
    }
    for (i=0; i < queue_count; i++) {
	resp = &stone->response_cache[stone->response_cache_count + i];
	resp->action_type = stone->proto_actions[act_num].action_type;
	resp->requires_decoded = 1;
	resp->proto_action_id = act_num;
	resp->o.multi.handler = func;
	resp->o.multi.client_data = client_data;
	resp->o.multi.free_func = free_func;
        resp->stage = cached_stage_for_action(&stone->proto_actions[act_num]);
	resp->reference_format = reference_formats[i];
	if (CMtrace_on(cm, EVerbose)) {
	    char *tmp;
	    if (resp->reference_format) {
		tmp = global_name_of_FMFormat(resp->reference_format);
	    } else {
		tmp = strdup("<none>");
	    }
	    fprintf(cm->CMTrace_file, "\tResponse %d for format \"%s\"(%p)\n", stone->response_cache_count+i, tmp, resp->reference_format);
	    free(tmp);
	}
    }
    stone->response_cache_count += queue_count;
    fix_response_cache(stone);
    return resp_num;
#else
    fprintf(stderr, "No code generation in FFS, action unavailable\n");
    return -1;
#endif    
}

extern EVstone
EVcreate_output_action(CManager cm, attr_list contact_list,
			   EVstone remote_stone)
{
    static int first = 1;
    if (first) {
	first = 0;
	printf("EVassoc_output_action is deprecated.  Please use EVassoc_bridge_action()\n");
    }
    return EVcreate_bridge_action(cm, contact_list, remote_stone);
}

extern EVstone
INT_EVcreate_bridge_action(CManager cm, attr_list contact_list,
			   EVstone remote_stone)
{
    EVstone stone = INT_EValloc_stone(cm);
    INT_EVassoc_bridge_action(cm, stone, contact_list, remote_stone);
    return stone;
}

extern EVstone
INT_EVcreate_thread_bridge_action(CManager cm, CManager target_cm,
				  EVstone target_stone)
{
    EVstone stone = INT_EValloc_stone(cm);
    INT_EVassoc_thread_bridge_action(cm, stone, target_cm, target_stone);
    return stone;
}

static int
is_bridge_stone(CManager cm, EVstone stone_num)
{
    event_path_data evp = cm->evp;
    stone_type stone = stone_struct(evp, stone_num);
    if (stone->default_action == -1) return 0;
    return (stone->proto_actions[stone->default_action].action_type == Action_Bridge);
}

extern EVaction
EVassoc_output_action(CManager cm, EVstone stone_num, attr_list contact_list,
		      EVstone remote_stone)
{
    static int first = 1;
    if (first) {
	first = 0;
	printf("EVassoc_output_action is deprecated.  Please use EVassoc_bridge_action()\n");
    }
    return EVassoc_bridge_action(cm, stone_num, contact_list, remote_stone);
}

extern EVaction
INT_EVassoc_bridge_action(CManager cm, EVstone stone_num, attr_list contact_list,
		      EVstone remote_stone)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int action_num;
    CMConnection conn = NULL;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = stone->proto_action_count;
    add_ref_attr_list(contact_list);
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding bridge action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, " remote stone target is %x\n", remote_stone);
    }
    if (getenv("NoLazyBridge")) {
        conn = INT_CMget_conn(cm, contact_list);
        if (conn == NULL) {
            if (CMtrace_on(cm, EVWarning)) {
                fprintf(cm->CMTrace_file, "EVassoc_bridge_action - failed to contact host at contact point \n\t");
                if (contact_list != NULL) {
                    fdump_attr_list(cm->CMTrace_file, contact_list);
                } else {
                    fprintf(cm->CMTrace_file, "NULL\n");
                }
                fprintf(cm->CMTrace_file, "Bridge action association failed for stone %x, outputting to remote stone %x\n",
                        stone_num, remote_stone);
            }
            return -1;
        }
        INT_CMconn_register_close_handler(conn, stone_close_handler, 
                                          (void*)(intptr_t)stone_num);
    }
    stone->proto_actions = realloc(stone->proto_actions, 
				   (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].action_type = Action_Bridge;
    stone->proto_actions[action_num].o.bri.remote_stone_id = remote_stone;
    stone->proto_actions[action_num].o.bri.remote_contact = contact_list;
    stone->proto_actions[action_num].o.bri.conn = conn;
    stone->default_action = action_num;
    stone->proto_action_count++;
    clear_response_cache(stone);
    return action_num;
}

extern EVaction
INT_EVassoc_thread_bridge_action(CManager cm, EVstone stone_num, 
				 CManager target_cm, EVstone target_stone)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int action_num;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = stone->proto_action_count;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding thread bridge action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, "\n");
    }
    stone->proto_actions = realloc(stone->proto_actions, 
				   (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].action_type = Action_Thread_Bridge;
    stone->proto_actions[action_num].o.thr_bri.target_cm = target_cm;
    stone->proto_actions[action_num].o.thr_bri.target_stone_id = target_stone;
    stone->proto_actions[action_num].o.thr_bri.target_cm_shutdown = 0;
    stone->default_action = action_num;
    stone->proto_action_count++;
    clear_response_cache(stone);
    return action_num;
}

extern EVstone
INT_EVcreate_split_action(CManager cm, EVstone *target_stone_list)
{
    EVstone stone = INT_EValloc_stone(cm);
    (void) INT_EVassoc_split_action(cm, stone, target_stone_list);
    return stone;
}

extern EVaction
INT_EVassoc_split_action(CManager cm, EVstone stone_num, 
		     EVstone *target_stone_list)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int action_num;
    int target_count = 0, i;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    action_num = stone->proto_action_count;
    stone->proto_actions = realloc(stone->proto_actions, 
				   (action_num + 1) * 
				   sizeof(stone->proto_actions[0]));
    memset(&stone->proto_actions[action_num], 0, 
	   sizeof(stone->proto_actions[0]));
    stone->proto_actions[action_num].action_type = Action_Split;
    while (target_stone_list && (target_stone_list[target_count] != -1)) {
	target_count++;
    }
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Adding Split action %d to ", action_num);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, ", %d target stones -> ", target_count);
	for (i=0; i < target_count; i++) {
	    fprintf(cm->CMTrace_file, "%x, ", target_stone_list[i]);
	}
	fprintf(cm->CMTrace_file, "\n");
    }
    for (i=0; i < target_count; i++) {
	INT_EVstone_add_split_target(cm, stone_num, target_stone_list[i]);
    }
    stone->output_count = target_count;
    stone->default_action = action_num;
    stone->proto_action_count++;
    clear_response_cache(stone);
    return action_num;
}

extern int
INT_EVaction_add_split_target(CManager cm, EVstone stone_num, 
			  EVaction action_num, EVstone new_stone_target)
{
    return INT_EVstone_add_split_target(cm, stone_num, new_stone_target);
}

extern int
INT_EVstone_add_split_target(CManager cm, EVstone stone_num, EVstone new_stone_target)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    EVstone *target_stone_list;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    if ((new_stone_target & 0x80000000) == 0x80000000) {
	new_stone_target = lookup_local_stone(evp, new_stone_target);
    }

    target_stone_list = stone->output_stone_ids;
    target_stone_list = realloc(target_stone_list, 
				(stone->output_count + 1) * sizeof(EVstone));
    target_stone_list[stone->output_count++] = new_stone_target;
    stone->output_stone_ids = target_stone_list;

    return 1;
}

extern void
INT_EVaction_remove_split_target(CManager cm, EVstone stone_num, 
			  EVaction action_num, EVstone stone_target)
{
    INT_EVstone_remove_split_target(cm, stone_num, stone_target);
}

extern void
INT_EVstone_remove_split_target(CManager cm, EVstone stone_num, 
				EVstone stone_target)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    EVstone *target_stone_list;
    int target_count = 0;

    stone = stone_struct(evp, stone_num);
    if (!stone) return;

    if ((stone_target & 0x80000000) == 0x80000000) {
	stone_target = lookup_local_stone(evp, stone_target);
    }

    target_stone_list = stone->output_stone_ids;
    if (!target_stone_list) return;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Removing split target %x from stone ", stone_target);
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_num);
	fprintf(cm->CMTrace_file, "\n");
    }
    while ((target_stone_list[target_count] != stone_target) && (target_count < stone->output_count)) {
	target_count++;
	CMtrace_out(cm, EVerbose, "    Found target to remove at position %d\n", target_count);
    }
    for ( ; target_count < stone->output_count-1; target_count++) {
	/* move them down, overwriting target */
	target_stone_list[target_count] = target_stone_list[target_count+1];
    }
    stone->output_count--;
}

void
INT_EVsend_stored(CManager cm, EVstone stone_num, EVaction action_num)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    event_item *item;
    proto_action *p;

    stone = stone_struct(evp, stone_num);
    if (!stone) return;

    p = &stone->proto_actions[action_num];
    while ((item = storage_queue_dequeue(cm, &p->o.store.queue)) != NULL) {
        internal_path_submit(cm, p->o.store.target_stone_id, item);
        p->o.store.num_stored--;
        return_event(evp, item);
        while (process_local_actions(cm));
    }
}

void
INT_EVstore_start_send(CManager cm, EVstone stone_num, EVaction action_num)
{
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    stone_type stone;
    proto_action *act;

    stone = stone_struct(evp, stone_num);
    if (!stone) return;

    act = &stone->proto_actions[action_num];

    if (act->o.store.num_stored == 0) return;
    if (act->o.store.is_sending == 1) return;

    act->o.store.is_sending = 1;
    act->o.store.is_paused = 0;
    as->events_in_play++;
    stone->pending_output++;
    
    /* make sure the local action loop is called soon */
    (void) INT_CMadd_delayed_task(cm, 0, 0, deferred_process_actions, NULL);
}

static int 
process_stone_pending_output(CManager cm, EVstone stone_num) {
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    stone_type stone = stone_struct(evp, stone_num);
    EVaction action_num;
    int found = 0;
    int more_pending = 0;

    for (action_num = 0; action_num < stone->proto_action_count
            && found < stone->pending_output; ++action_num) {
        proto_action *act = &stone->proto_actions[action_num];
        if (act->action_type == Action_Store &&
            act->o.store.is_sending && !act->o.store.is_paused) {
            event_item *item;
            ++found;
            item = storage_queue_dequeue(cm, &act->o.store.queue);
            assert(item->ref_count > 0);
            assert(!stone_struct(evp, act->o.store.target_stone_id)->is_stalled);
            internal_path_submit(cm, act->o.store.target_stone_id, item);
            /* printf("submitted one <%d / %d>\n", as->events_in_play, as->last_active_stone); */
            act->o.store.num_stored--;
            /* return_event(evp, item); */
            if (!act->o.store.num_stored) {
                /* printf("stopping sending\n"); */
                act->o.store.is_sending = act->o.store.is_paused = 0;
                as->events_in_play--;
                stone->pending_output--;
                found--;
            } else {
                ++more_pending;
            }
        }
    }
    return more_pending;
}

int
INT_EVstore_is_sending(CManager cm, EVstone stone_num, EVaction action_num)
{
    event_path_data evp = cm->evp;
    stone_type stone;;
    proto_action *p;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    p = &stone->proto_actions[action_num];
    return p->o.store.is_sending;
}

int
INT_EVstore_count(CManager cm, EVstone stone_num, EVaction action_num)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    proto_action *p;

    stone = stone_struct(evp, stone_num);
    if (!stone) return -1;

    p = &stone->proto_actions[action_num];
    return p->o.store.num_stored;
}

struct source_info {
    EVstone to_stone;
    void *user_data;
    enum { SOURCE_ACTION, SOURCE_REMOTE } type;
    EVstone stone;
    union {
        struct {
            EVaction action;
            int would_recurse;
        } action;
        struct {
            CMConnection conn;
        } remote;
    }u;
};

typedef void (*ForeachSourceCB)(CManager, struct source_info *);

static void
foreach_source_inner(CManager cm, EVstone to_stone, char *seen,
        ForeachSourceCB cb, struct source_info *info) {
    event_path_data evp = cm->evp;
    EVstone cur_stone;
    for (cur_stone = evp->stone_base_num; cur_stone < evp->stone_count + evp->stone_base_num; ++cur_stone) {
        EVaction cur_action;
        stone_type stone = stone_struct(evp, cur_stone);
        if (seen[cur_stone - evp->stone_base_num]) continue;
	if (!stone) continue;
        if (stone->local_id == -1) continue;
        if (cur_stone == to_stone) {
            if (stone->last_remote_source != NULL) {
                info->type = SOURCE_REMOTE;
                info->stone = cur_stone;
                info->u.remote.conn = stone->last_remote_source;
                cb(cm, info);
            }
        } else {
            int was_stalled = stone->is_stalled;
	    int i;
	    int matches = 0;
	    int matches_recursive = 0;
	    for (i = 0; i < stone->output_count; ++i) {
		if (to_stone == stone->output_stone_ids[i]) {
		    ++matches;
		    ++matches_recursive;
		}
	    }
            for (cur_action = 0; cur_action < stone->proto_action_count; ++cur_action) {
                proto_action *act;
                act = &stone->proto_actions[cur_action];
                switch (act->action_type) {
                case Action_Store:
                    if (act->o.store.target_stone_id == to_stone) {
                        ++matches;
                    }
                    break;
                case Action_Split:
                    break;
                case Action_Filter:
                    if (to_stone == act->o.term.target_stone_id) {
                        ++matches;
                        ++matches_recursive;
                    }
                    break;
                case Action_Immediate:
                    break;
                case Action_Terminal:
                    /* nothing to do */
                    break;
                default: ;
                    /* printf("source searching: unhandled type %d\n", act->action_type); */
                    /* TODO unhandled cases
                            - Source handles?
                     */
                }
                if (matches) {
                    info->type = SOURCE_ACTION;
                    info->stone = cur_stone;
                    info->u.action.action = cur_action;
                    info->u.action.would_recurse = matches_recursive;
                    cb(cm, info);
                }
                /* If a stone is stalled, conceptually no traffic is going to it.
                 * Note that the only way we should be trying to recurse through
                 * a stalled stone is if it were manually stalled.
                 */
                if (!was_stalled && matches_recursive) {
                    /* avoid infinite recursion */
                    seen[cur_stone - evp->stone_base_num] = 1;
                    foreach_source_inner(cm, cur_stone, seen, cb, info);
                    seen[cur_stone - evp->stone_base_num] = 0; /* XXX */
                }
            }
        }
    }
}

static void
foreach_source(CManager cm, EVstone to_stone, ForeachSourceCB cb, void *user_data) {
    char* seen = calloc(1, cm->evp->stone_count); /* XXX try to keep static */
    struct source_info info;
    info.user_data = user_data;
    info.to_stone = to_stone;
    foreach_source_inner(cm, to_stone, seen, cb, &info);
    free(seen);
}

static void
backpressure_transition(CManager cm, EVstone s, stall_source src, int new_value);

enum { CONTROL_SQUELCH, CONTROL_UNSQUELCH };

static void
register_backpressure_callback(CManager cm, EVstone s, EVSubmitCallbackFunc cb, void *user_data) {
    stall_callback *new_cb = INT_CMmalloc(sizeof(stall_callback));
    stone_type stone = stone_struct(cm->evp, s);
    assert(CManager_locked(cm));
    new_cb->cb = cb;
    new_cb->user_data = user_data;
    new_cb->next = stone->unstall_callbacks;
    stone->unstall_callbacks = new_cb;
}

static void
do_backpressure_unstall_callbacks(CManager cm, EVstone stone_id) {
    stone_type stone = stone_struct(cm->evp, stone_id);
    stall_callback *cur = stone->unstall_callbacks;
    assert(CManager_locked(cm));
    if (!cur) return;
    stone->unstall_callbacks = NULL;
    CManager_unlock(cm);
    while (cur) {
        stall_callback *next = cur->next;
        (cur->cb)(cm, stone_id, cur->user_data);
        INT_CMfree(cur);
        cur = next;
    }
    CManager_lock(cm);
}


static void
backpressure_set_one(CManager cm, struct source_info *info)
{
    event_path_data evp = cm->evp;
    action_state as = evp->as;
    assert(as->events_in_play >= 0);
    int s = info->to_stone;
    stone_type to_stone = stone_struct(evp,s);
    stone_type stone = stone_struct(evp, info->stone);
    switch (info->type) {
    case SOURCE_ACTION: {
            proto_action *act = &stone->proto_actions[info->u.action.action];
            
            if (info->u.action.would_recurse) {
                /* If we might be stalling the stone that has sources, we mark it as stalled before
                 * calling * backpressure_transition() because we are already recursing to its sources
                 * and performing actions as if the stone is stalled. If the stone unstalls
                 * as a result of the transition, we do want to do the recursive search
                 * since we do not recurse through stalled stones (though we do set would_recurse
                 * for them). 
                 */
                if (to_stone->is_stalled) {
                    printf("recurse stall %d\n", info->stone);
                    stone->is_stalled = 1;
                } else {
                    printf("recurse unstall %d\n", info->stone);
                    do_backpressure_unstall_callbacks(cm, info->stone);
                }
                /* TODO for, e.g., split stones check that we should actually unstall it 
                 * (maybe just count our upstream stall depth?) 
                 */
                backpressure_transition(cm, info->stone, Stall_Upstream, to_stone->is_stalled);
            }
            switch (act->action_type) {
            case Action_Store:
                {
                    struct storage_proto_vals *store = &act->o.store;
                    if (store->is_paused != to_stone->is_stalled) {
                        store->is_paused = to_stone->is_stalled;
                        if (store->is_sending) {
                            if (store->is_paused) {
                                --as->events_in_play;
                                --stone->pending_output;
                            } else {
                                ++as->events_in_play;
                                ++stone->pending_output;
                                (void) INT_CMadd_delayed_task(cm, 0, 0, deferred_process_actions, NULL);
                            }
                        }
                    }
                }
                break;
            default:;
                /* printf("unhandled action %d\n", act->action_type); */
                /* TODO more cases? */
            }
        }
        break;
    case SOURCE_REMOTE: {
            if (to_stone->is_stalled) {
                if (stone->squelch_depth++ == 0) {
                    INT_CMwrite_evcontrol(info->u.remote.conn, CONTROL_SQUELCH, info->stone);
                }
            } else if (0 == --stone->squelch_depth) {
                INT_CMwrite_evcontrol(info->u.remote.conn, CONTROL_UNSQUELCH, info->stone);
            }
        }
        break;
    default:
        /* XXX */
        ;
    }
}

static void
backpressure_set(CManager cm, EVstone to_stone, int stalledp) {
    event_path_data evp = cm->evp;
    stone_type stone = stone_struct(evp, to_stone);
    assert(cm->evp->use_backpressure);
    if (stone->is_stalled == stalledp) {
        return;
    }
    /*
    if (stalledp) {
        printf("backpressure stalled %d\n", (int) to_stone);
    } else {
        printf("backpressure unstalled %d\n", (int) to_stone);
    }
    */
    stone->is_stalled = stalledp;
    if (!stalledp) {
        do_backpressure_unstall_callbacks(cm, to_stone);
    }
    foreach_source(cm, to_stone, backpressure_set_one, NULL);
}

#if defined (__INTEL_COMPILER)
//  Allow extern declarations with no prior decl
#  pragma warning (disable: 188)
#endif
static void
backpressure_transition(CManager cm, EVstone s, stall_source src, int new_value) {
    stone_type stone = stone_struct(cm->evp, s);
    assert(cm->evp->use_backpressure);
    if (new_value) {
        stone->stall_from |= src;
    } else {
        stone->stall_from &= ~src;
    }
    backpressure_set(cm, s, stone->stall_from ? 1 : 0);
}

static void
backpressure_check(CManager cm, EVstone s) {
    if (cm->evp->use_backpressure) {
        stone_type stone = stone_struct(cm->evp, s);
        int old_stalled = stone->is_stalled;
        int low_threshold = 50, high_threshold = 200;
        if (stone->stone_attrs) {
	    static atom_t EV_BACKPRESSURE_HIGH = -1;
	    static atom_t EV_BACKPRESSURE_LOW = -1;
	    if (EV_BACKPRESSURE_HIGH == -1) {
		EV_BACKPRESSURE_HIGH = attr_atom_from_string("EV_BACKPRESSURE_HIGH");
		EV_BACKPRESSURE_LOW = attr_atom_from_string("EV_BACKPRESSURE_LOW");
	    }
            get_int_attr(stone->stone_attrs, EV_BACKPRESSURE_HIGH, &high_threshold);
            get_int_attr(stone->stone_attrs, EV_BACKPRESSURE_LOW, &low_threshold);
        }
        backpressure_transition(cm, s, Stall_Overload,
            (stone->queue_size > (old_stalled ? low_threshold : high_threshold)));
    }
}

int
INT_EVsubmit_or_wait(EVsource source, void *data, attr_list attrs,
                    EVSubmitCallbackFunc cb, void *user_data) {
    CManager cm = source->cm;
    stone_type stone;

    stone = stone_struct(cm->evp, source->local_stone_id);
    if (!stone) return -1;

    if (stone->is_stalled) {
        register_backpressure_callback(source->cm, source->local_stone_id, cb, user_data);
        return 0;
    } else {
        INT_EVsubmit(source, data, attrs);
        return 1;
    }
}

int
INT_EVsubmit_encoded_or_wait(CManager cm, EVstone s, void *data, int data_len, 
                    attr_list attrs, EVSubmitCallbackFunc cb, void *user_data) {
    stone_type stone;
    stone = stone_struct(cm->evp, s);
    if (!stone) return -1;

    if (stone->is_stalled) {
        register_backpressure_callback(cm, s, cb, user_data);
        return 0;
    } else {
        INT_EVsubmit_encoded(cm, s, data, data_len, attrs);
        return 1;
    }
}

void
INT_EVstall_stone(CManager cm, EVstone s) {
    backpressure_transition(cm, s, Stall_Requested, 1);
}

void
INT_EVunstall_stone(CManager cm, EVstone s) {
    backpressure_transition(cm, s, Stall_Requested, 0);
}

void
INT_EVhandle_control_message(CManager cm, CMConnection conn, unsigned char type, int arg) {
    event_path_data evp = cm->evp;
    switch (type) {
    case CONTROL_SQUELCH:
    case CONTROL_UNSQUELCH: {
            int s;
            stone_type stone;
            for (s = evp->stone_base_num; s < evp->stone_count + evp->stone_base_num; ++s) {
                stone = stone_struct(evp, s);
                if (is_bridge_stone(cm, s) &&  stone->proto_actions[stone->default_action].o.bri.conn == conn
                        && stone->proto_actions[stone->default_action].o.bri.remote_stone_id == arg) {
                    backpressure_transition(cm, s, Stall_Squelch, type == CONTROL_SQUELCH);
                }
            }
        }
        break;
    default:
        assert(FALSE);
    }
}

extern void
do_local_actions(CManager cm)
{
    while (process_local_actions(cm));
}

extern FMFormat
EVregister_format_set(CManager cm, FMStructDescList list)
{
    FMFormat format = NULL;

    if (list[0].format_name != NULL) {
	format = register_data_format(cm->evp->fmc, list);
    }
    return format;
}

static void
EVauto_submit_func(CManager cm, void* vstone)
{
    int stone_num = (int)(intptr_t) vstone;
    event_item *event;
    CManager_lock(cm);
    event = get_free_event(cm->evp);
    event->event_encoded = 0;
    event->decoded_event = NULL;
    event->reference_format = NULL;
    event->format = NULL;
    event->free_func = NULL;
    event->attrs = NULL;
    event->cm = cm;
    internal_path_submit(cm, stone_num, event);
    while (process_local_actions(cm));
    return_event(cm->evp, event);
    CManager_unlock(cm);
}

extern EVstone
INT_EVcreate_auto_stone(CManager cm, int period_sec, int period_usec, 
			char *action_spec, EVstone out_stone)
{
    EVstone stone = INT_EValloc_stone(cm);
    EVaction action = INT_EVassoc_immediate_action(cm, stone, action_spec, NULL);
    INT_EVaction_set_output(cm, stone, action, 0, out_stone);
    INT_EVenable_auto_stone(cm, stone, period_sec, period_usec);
    return stone;
}

extern void
INT_EVenable_auto_stone(CManager cm, EVstone stone_num, int period_sec, 
		    int period_usec)
{
    CMTaskHandle handle;
    stone_type stone;
    int acceptable_action = 0;
    int action_num;

    stone = stone_struct(cm->evp, stone_num);
    if (!stone) return;
    for (action_num = 0; action_num < stone->proto_action_count; action_num ++) {
	switch (stone->proto_actions[action_num].action_type) {
	case Action_Filter: case Action_Immediate: case Action_Multi:
	    acceptable_action++;
	    break;
	default:
	    break;
	}
    }
    if (acceptable_action == 0) {
	printf("Warning!  Enabling auto events on ");
	fprint_stone_identifier(cm->CMTrace_file, cm->evp, stone_num);
	printf(", but no acceptable actions found!\n");
    }
    handle = INT_CMadd_periodic_task(cm, period_sec, period_usec,
				     EVauto_submit_func, 
				     (void*)(intptr_t)stone_num);
    stone->periodic_handle = handle;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Enabling auto events on ");
	fprint_stone_identifier(cm->CMTrace_file, cm->evp, stone_num);
	fprintf(cm->CMTrace_file, "\n");
    }
}



extern EVsource
INT_EVcreate_submit_handle(CManager cm, EVstone stone, FMStructDescList data_format)
{
    EVsource source = malloc(sizeof(*source));
    memset(source, 0, sizeof(*source));
    source->local_stone_id = stone;
    source->cm = cm;
    source->preencoded = 0;
    if (data_format != NULL) {
	source->format = INT_CMregister_format(cm, data_format);
	source->reference_format = EVregister_format_set(cm, data_format);
    };
    return source;
}

extern EVsource
INT_EVcreate_submit_handle_free(CManager cm, EVstone stone, 
			    FMStructDescList data_format, 
			    EVFreeFunction free_func, void *free_data)
{
    EVsource source = malloc(sizeof(*source));
    memset(source, 0, sizeof(*source));
    source->local_stone_id = stone;
    source->cm = cm;
    source->format = INT_CMregister_format(cm, data_format);
    source->reference_format = EVregister_format_set(cm, data_format);
    source->free_func = free_func;
    source->free_data = free_data;
    source->preencoded = 0;
    return source;
}

extern void
INT_EVfree_source(EVsource source)
{
    free(source);
}

static void
reference_event(event_item *event)
{
    event->ref_count++;
}

extern void
internal_cm_network_submit(CManager cm, CMbuffer cm_data_buf, 
			   attr_list attrs, CMConnection conn, 
			   void *buffer, size_t length, int stone_id)
{
    event_path_data evp = cm->evp;
    event_item *event = get_free_event(evp);
    stone_type stone;
    (void)cm_data_buf;
    FMFormat reference_format = FMformat_from_ID(evp->fmc, buffer);
    if (reference_format == NULL) {
	printf("FFS failure format not found, incoming data incomprehensible, ignored\n");
	fprintf(cm->CMTrace_file, "Buffer format is ");
	fprint_server_ID(cm->CMTrace_file, buffer);
	fprintf(cm->CMTrace_file, "\n");

	printf("  This could be a FFS format server issue, a CMSelfFormats issue, a transport corruption issue, or something else...\n");
	return;
    }
    event->contents = Event_CM_Owned;
    event->event_encoded = 1;
    event->event_len = length;
    event->encoded_event = buffer;
    event->reference_format = reference_format;
    event->attrs = CMadd_ref_attr_list(cm, attrs);
    event->cm = cm;
    event->format = NULL;
    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Event coming in from network to ");
	fprint_stone_identifier(cm->CMTrace_file, evp, stone_id);
	fprintf(cm->CMTrace_file, "\n");
    }
    if (CMtrace_on(conn->cm, EVerbose)) {
	static int dump_char_limit = 256;
	static int warned = 0;
	static int size_set = 0;
	int r;
	if (size_set == 0) {
	    char *size_str = getenv("CMDumpSize");
	    size_set++;
	    if (size_str != NULL) {
		dump_char_limit = atoi(size_str);
	    }
	}
	fprintf(cm->CMTrace_file, "CM - record type %s, contents are:\n  ", global_name_of_FMFormat(event->reference_format));
	r = FMfdump_encoded_data(cm->CMTrace_file, event->reference_format,
				event->encoded_event, dump_char_limit);
	if (r && !warned) {
	    fprintf(cm->CMTrace_file, "\n\n  ****  Warning **** CM record dump truncated\n");
	    fprintf(cm->CMTrace_file, "  To change size limits, set CMDumpSize environment variable.\n\n\n");
	    warned++;
	}
    }
    INT_CMtake_buffer(cm, buffer);
    event->cm = cm;
    stone = stone_struct(evp, stone_id);
    if (stone->squelch_depth == 0) {
        stone->last_remote_source = conn;
    }
    internal_path_submit(cm, stone_id, event);
    return_event(evp, event);
    while (process_local_actions(cm));
}

static void free_ioBuffer(void *event_data, void *client_data)
{
    free_FFSBuffer((FFSBuffer) client_data);
}    

static event_item *
reassign_memory_event(CManager cm, event_item *event, int do_decode)
{
    /* 
     *  The old event item is enqueued on stones elsewhere in EVPath.  We must make the data memory 
     *  it references available application reuse.  So, we are going to:
     *   1 - encode and then decode the event (to get new, clean CM-owned data)
     *   2 - modify the old event structure to reference only the new data
     *   3 - create and return tmp_event, which references only the previous data and 
     *       has a reference count of 1 (so that return_event() does what it needs to do).
     */
    event_item *tmp_event = get_free_event(cm->evp);
    FFSContext tmp_context;
    FFSTypeHandle format;
    void *decode_buffer;

    CMtrace_out(cm, EVerbose, "Doing deep copy to free up event before returning from EVsubmit()\n");
    *tmp_event = *event;
    tmp_event->ref_count = 1;   /* we're going to make sure the enqueued events don't reference anything that this does */
    tmp_event->attrs = CMadd_ref_attr_list(cm, event->attrs);
    event->free_func = NULL;   /* if these were present, no longer applicable */
    event->free_arg = NULL;
    event->cm = cm;
    encode_event(cm, event);  /* Copy all data to an FFS encode buffer in event->ioBuffer */
    event->decoded_event = NULL;      /* we're not touching this anymore */
    event->contents = Event_Freeable; /* it's all our data */
    event->free_arg = event->ioBuffer;
    event->ioBuffer = NULL;
    event->free_func = free_ioBuffer;
    /*
     * The new event type is a bit unique in EVPath.  event->ioBuffer is meant to hold only encoded data, but we 
     * want to decode that stuff, so we don't leave our data there.   Rather, the FFSbuffer created by this will be freed 
     * via the Event_Freeable mechanism.
     */
    if (do_decode) {
	tmp_context = create_FFSContext_FM(cm->evp->fmc);
	format = FFSTypeHandle_from_encode(tmp_context, event->encoded_event);
	establish_conversion(tmp_context, format, format_list_of_FMFormat(event->reference_format));

	if (!FFSdecode_in_place(tmp_context, event->encoded_event, &decode_buffer)) {
	    printf("Decode failed\n");
	    return 0;
	}
	event->decoded_event = decode_buffer;
	event->encoded_event = NULL;
	event->event_encoded = 0;
	free_FFSContext(tmp_context);
    }
    event->ref_count--;   /* we've essentially split the event.  tmp_event will be dereferenced */
    return tmp_event;
}

extern void
INT_EVsubmit_general(EVsource source, void *data, EVFreeFunction free_func, 
		 attr_list attrs)
{
    event_item *event = get_free_event(source->cm->evp);
    event->contents = Event_App_Owned;
    event->event_encoded = 0;
    event->decoded_event = data;
    event->reference_format = source->reference_format;
    event->format = source->format;
    event->free_func = free_func;
    event->free_arg = data;
    event->cm = source->cm;
    event->attrs = CMadd_ref_attr_list(source->cm, attrs);
    internal_path_submit(source->cm, source->local_stone_id, event);
    while (process_local_actions(source->cm));
    return_event(source->cm->evp, event);
}
    
void
INT_EVsubmit(EVsource source, void *data, attr_list attrs)
{
    event_path_data evp = source->cm->evp;
    event_item *event;
    if (source->local_stone_id == -1) return;  /* not connected */
    event = get_free_event(evp);
    if (source->free_func != NULL) {
	event->contents = Event_Freeable;
    } else {
	event->contents = Event_App_Owned;
    }
    event->cm = source->cm;
    if (source->preencoded) {
	event->event_encoded = 1;
	event->encoded_event = data;
	event->reference_format = FMFormat_of_original(FFSTypeHandle_from_encode(evp->ffsc, 
							    data));
    } else {
	event->event_encoded = 0;
	event->decoded_event = data;
	event->reference_format = source->reference_format;
	event->format = source->format;
    }
    event->free_func = source->free_func;
    event->free_arg = source->free_data;
    event->attrs = CMadd_ref_attr_list(source->cm, attrs);
    internal_path_submit(source->cm, source->local_stone_id, event);
    while (process_local_actions(source->cm));
    if (event->ref_count != 1 && (event->contents == Event_App_Owned)) {
	event = reassign_memory_event(source->cm, event, 1);  /* reassign memory */
    }
    return_event(source->cm->evp, event);
}

void
INT_EVsubmit_encoded(CManager cm, EVstone stone, void *data, size_t data_len, attr_list attrs)
{
    event_path_data evp = cm->evp;
    event_item *event = get_free_event(evp);
    if (stone_struct(evp, stone) == NULL) return;

    event->contents = Event_App_Owned;
    event->event_encoded = 1;
    event->encoded_event = data;
    event->cm = cm;
    event->event_len = data_len;
    event->reference_format = FMFormat_of_original(FFSTypeHandle_from_encode(evp->ffsc, 
							data));

    event->attrs = CMadd_ref_attr_list(cm, attrs);
    internal_path_submit(cm, stone, event);
    while (process_local_actions(cm));
    return_event(cm->evp, event);
}

static void
free_evp(CManager cm, void *not_used)
{
    event_path_data evp = cm->evp;
    int s;
    (void)not_used;
    CMtrace_out(cm, CMFreeVerbose, "Freeing evpath information, evp %p\n", evp);
    for (s = 0 ; s < evp->stone_count; s++) {
	INT_EVfree_stone(cm, s + evp->stone_base_num);
    }
    cm->evp = NULL;
    if (evp == NULL) return;
    free(evp->stone_map);
    free(evp->as);

    free_FFSContext(evp->ffsc);
/*    free_FMcontext(evp->fmc);   unnecessary?*/
    while (evp->queue_items_free_list != NULL) {
	queue_item *tmp = evp->queue_items_free_list->next;
	free(evp->queue_items_free_list);
	evp->queue_items_free_list = tmp;
    }
    if (evp->sources) {
	int i;
	for (i=0; i<evp->source_count; i++) {
	    if (evp->sources[i].name) free(evp->sources[i].name);
	}
	free(evp->sources);
    }
    if (evp->sink_handlers) {
	int i;
	for (i=0; i<evp->sink_handler_count; i++) {
	    if (evp->sink_handlers[i].name) free(evp->sink_handlers[i].name);
	}
	free(evp->sink_handlers);
    }
    if (evp->stone_lookup_table) free(evp->stone_lookup_table);
    if (evp->externs) free(evp->externs);
    thr_mutex_free(evp->lock);
    free(evp);
}

void
EVPinit(CManager cm)
{
    static int first_evpinit = 1;
    cm->evp = INT_CMmalloc(sizeof( struct _event_path_data));
    memset(cm->evp, 0, sizeof( struct _event_path_data));
    cm->evp->ffsc = cm->FFScontext;
    cm->evp->fmc = FMContext_from_FFS(cm->evp->ffsc);
    cm->evp->stone_base_num = 0;
    if (!first_evpinit) {
	/* 
	 * after the first evpinit, start stones at random base number,
	 * just so that we're more likely to catch mitmatched stone/CM
	 * combos in threaded situations.
	*/
	srand((int)time(NULL));
	while (cm->evp->stone_base_num == 0) {
	    cm->evp->stone_base_num = rand() & 0xffff;
	}
    }
    CMtrace_out(cm, EVerbose, "INITATED EVPATH, base stone num is %x\n", 
		cm->evp->stone_base_num);
    first_evpinit = 0;
    cm->evp->queue_items_free_list = NULL;
    thr_mutex_init(cm->evp->lock);
    internal_add_shutdown_task(cm, free_evp, NULL, FREE_TASK);
    {
        char *backpressure_env;
        backpressure_env = getenv("EVBackpressure");
        if (backpressure_env && atoi(backpressure_env) != 0) {
            cm->evp->use_backpressure = 1;
        } else {
            cm->evp->use_backpressure = 0;
        }
    }
    INT_CMadd_poll(cm, deferred_process_actions, NULL);
    REVPinit(cm);
}

extern void
INT_EVregister_close_handler(CManager cm, EVStoneCloseHandlerFunc handler, void *client_data)
{
    cm->evp->app_stone_close_handler = handler;
    cm->evp->app_stone_close_data = client_data;
}
    
extern int
INT_EVtake_event_buffer(CManager cm, void *event)
{
    queue_item *item;
    event_item *cur = NULL;
    queue_item *running_events = cm->evp->current_event_list;
    event_path_data evp = cm->evp;

    while(running_events != NULL) {
	cur = running_events->item;
	if (!(((cur->decoded_event <= event) &&
	       ((char *) event <= ((char *) cur->decoded_event + cur->event_len))) ||
	      ((cur->encoded_event <= event) &&
	       ((char *) event <= ((char *) cur->encoded_event + cur->event_len))))) {
	    cur = NULL;
	} else {
	    break;
	}
	running_events = running_events->next;
    }
    if (cur == NULL) {
	fprintf(stderr,
		"Event address (%p) in INT_EVtake_event_buffer does not match currently executing event on this CM.\n",
		event);
	return 0;
    }

/*    if (cur->block_rec == NULL) {
	static int take_event_warning = 0;
	if (take_event_warning == 0) {
	    fprintf(stderr,
		    "Warning:  INT_EVtake_event_buffer called on an event submitted with \n    INT_EVsubmit_event(), INT_EVsubmit_typed_event() or INT_EVsubmit_eventV() .\n    This violates ECho event data memory handling requirements.  See \n    http://www.cc.gatech.edu/systems/projects/ECho/event_memory.html\n");
	    take_event_warning++;
	}
	return 0;
    }
*/
    if (evp->queue_items_free_list == NULL) {
	item = malloc(sizeof(*item));
    } else {
	item = evp->queue_items_free_list;
	evp->queue_items_free_list = item->next;
    }
    item->item = cur;
    reference_event(cur);
    item->next = evp->taken_events_list;
    evp->taken_events_list = item;
    return 1;
}

void
INT_EVreturn_event_buffer(CManager cm, void *event)
{
    event_path_data evp = cm->evp;
    queue_item *tmp, *last = NULL;
    /* search through list for event and then dereference it */
    (void)cm;
    (void)event;

    tmp = evp->taken_events_list;
    while (tmp != NULL) {
	if (((tmp->item->decoded_event <= event) &&
	    ((char *) event <= ((char *) tmp->item->decoded_event + tmp->item->event_len))) ||
	    ((tmp->item->encoded_event <= event) &&
	     ((char *) event <= ((char *) tmp->item->encoded_event + tmp->item->event_len)))) {
	    if (last == NULL) {
		evp->taken_events_list = tmp->next;
	    } else {
		last->next = tmp->next;
	    }
	    return_event(cm->evp, tmp->item);
	    tmp->next = evp->queue_items_free_list;
	    evp->queue_items_free_list = tmp;
	    return;
	}
	last = tmp;
	tmp = tmp->next;
    }
    fprintf(stderr, "Event %p not found in taken events list\n",
	    event);
}

extern FMFormat
INT_EVget_src_ref_format(EVsource source)
{
    return source->reference_format;
}

extern int
INT_EVfreeze_stone(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    stone = stone_struct(evp, stone_id);
    if (!stone) return -1;
    stone->is_frozen = 1;
    return 1;	
}

extern int
INT_EVunfreeze_stone(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    CMTaskHandle handle;
    stone = stone_struct(evp, stone_id);
    if (!stone) return -1;

    stone->is_frozen = 0;
    /* ensure that we run the process_actions loop soon so the stone's
       pending events (or pending output) won't be ignored */
    handle = INT_CMadd_delayed_task(cm, 0, 0, deferred_process_actions, NULL);
    free(handle); /* we don't need this */
    return 1;	
}

extern int
INT_EVdrain_stone(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    int count = 0;

    stone = stone_struct(evp, stone_id);
    if (!stone) return -1;

    stone->is_draining = 1;
    while(stone->is_processing || stone->is_outputting ||
	  (stone->queue->queue_head != NULL)) {
	if (count++ > 20) {
	    /* if not drained in 10 seconds, fail */
	    return 0;
	}
	INT_CMusleep(cm, 500000);
    }
    stone->is_draining = 2;
    return 1;
}

extern EVevent_list
INT_EVextract_stone_events(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    EVevent_list list = malloc(sizeof(list[0]));

    list[0].length = (size_t)-1;
    stone = stone_struct(evp, stone_id);
    if (!stone) return NULL;
    list = extract_events_from_queue(cm, stone->queue, list);
    return list;
}

extern int
INT_EVtransfer_events(CManager cm, EVstone src_stone_id, EVstone dest_stone_id)
{
    event_path_data evp = cm->evp;
    stone_type src_stone, dest_stone;
    int count = 0;
    queue_item * item;
    event_item * event;

    src_stone = stone_struct(evp, src_stone_id);
    if (!src_stone) return -1;
    dest_stone = stone_struct(evp, dest_stone_id);
    if (!dest_stone) return -1;
    item = src_stone->queue->queue_head;
    while (item != NULL) {
	queue_item *next = item->next;
	event = dequeue_item(cm, src_stone, item);
	internal_path_submit(cm, dest_stone_id, event);
	return_event(evp, event);
	count++;
	item = next;
    }
    return count;
}

extern attr_list
INT_EVextract_attr_list(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    stone = stone_struct(evp, stone_id);
    if (!stone) return NULL;
    return(stone->stone_attrs);
}

extern void
INT_EVset_attr_list(CManager cm, EVstone stone_id, attr_list stone_attrs)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    stone = stone_struct(evp, stone_id);
    if (!stone) return;
    if (stone->stone_attrs) free_attr_list(stone->stone_attrs);
    stone->stone_attrs = stone_attrs;
    add_ref_attr_list(stone_attrs);
}

EVevent_list
extract_events_from_queue(CManager cm, queue_ptr que, EVevent_list list)
{
    EVevent_list current_entry = NULL;
    queue_item *first = NULL, *last = NULL;
    int num_of_elements = 0;
    
    first = que->queue_head;
    last = que->queue_tail;
                
    while (list[num_of_elements].length != (size_t)-1) num_of_elements++;
    while(first != NULL && last != NULL) {
	list = (EVevent_list) realloc (list, (num_of_elements + 2) * sizeof(list[0]));
	current_entry = &list[num_of_elements];
	if((first->item->event_encoded) || (first->item->ioBuffer != NULL)) {
	    current_entry->length = first->item->event_len;
	    current_entry->buffer = first->item->encoded_event;
	} else {
	    encode_event(cm, first->item);
	    current_entry->length = first->item->event_len;
	    current_entry->buffer = first->item->encoded_event;
	}
	num_of_elements++;
	first = first->next;
    }
    list[num_of_elements].length = (size_t)-1;
    return list;
}

extern int
INT_EVdestroy_stone(CManager cm, EVstone stone_id)
{
    event_path_data evp = cm->evp;
    stone_type stone;
    stone = stone_struct(evp, stone_id);
    if (!stone) return -1;
    INT_EVdrain_stone(cm, stone_id);
    empty_queue(evp, stone->queue);
    INT_EVfree_stone(cm, stone_id);  
    return 1;      
} 
    


