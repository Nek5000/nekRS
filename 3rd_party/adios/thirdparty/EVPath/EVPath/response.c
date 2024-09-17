#include "config.h"
#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "evpath.h"
#ifdef HAVE_COD_H
#include "cod.h"
#else
#define cod_code void*
#define cod_exec_context void*
#define cod_assoc_client_data(x, y, z) 0
#define cod_get_client_data(x, y) NULL
#define cod_create_exec_context(x) 0
#define cod_exec_context_free(x) 0
#define cod_code_free(x) 0
#define cod_assoc_externs(x,y) 0
#define cod_parse_for_context(x,y) 0
#define cod_set_closure(x,y,z) 0
#define cod_add_int_constant_to_parse_context(name, i, context) 0
#endif
#include "cm_internal.h"
#include "dlloader.h"

typedef enum {Response_Filter, Response_Transform, Response_Router, Response_Multityped} response_types;

struct terminal_spec {
    FMStructDescList format_list;
    void *handler;
    void *client_data;
};

struct filter_spec {
    FMStructDescList format_list;
    char *function;
    void *client_data;
    FMFormat reference_format;
};

struct transform_spec {
    FMStructDescList in_format_list;
    FMStructDescList out_format_list;
    char *function;
    void *client_data;
    FMFormat reference_input_format;
    FMFormat reference_output_format;
    EVsource source_handle;
    int output_base_struct_size;
};

struct multityped_spec {
    FMStructDescList *struct_list;
    char *function;
    void *client_data;
    int accept_anonymous;
    FMFormat *reference_input_format_list;
};

typedef struct response_spec {
    response_types response_type;
    union {
	struct terminal_spec term;
	struct filter_spec filter;
	struct transform_spec transform;
	struct multityped_spec multityped;
    }u;
} *handler_list;

struct filter_instance {
    int (*func_ptr)(void *, attr_list);
    cod_code code;
    cod_exec_context ec;
    void *client_data;
};

struct transform_instance {
    int (*func_ptr)(void *, void*, attr_list, attr_list);
    cod_code code;
    cod_exec_context ec;
    int out_size;
    void *client_data;
    FMFormat out_format;
};

struct queued_instance {
    int ref_count;
    cod_code code;
    cod_exec_context ec;
    void *client_data;
    FMFormat *formats;
};

typedef struct response_instance {
    response_types response_type;
    int stone;
    int proto_action_id;
    union {
	struct filter_instance filter;
	struct transform_instance transform;
	struct queued_instance queued;
    }u;
} *response_instance;


static char *
add_FMfieldlist_to_string(char *str, FMStructDescRec *f)
{
    int index, field_count = 0;
    FMFieldList list = f->field_list;
    int len = (int) strlen(str);
    char *tmp_str;
    len += (int)strlen(f->format_name) + 5 + 35 + 20;
    str = realloc(str, len);
    while(list && (list[field_count].field_name != NULL)) field_count++;
    tmp_str = str + strlen(str);
    sprintf(tmp_str, "FMFormat \"%s\" StructSize %d FieldCount %d\n",
	    f->format_name, f->struct_size, field_count);
    for (index = 0; index < field_count; index++) {
	len += (int)strlen(list[index].field_name) + (int) strlen(list[index].field_type) + 50;
	str = realloc(str, len);
	tmp_str = str + strlen(str);
	sprintf(tmp_str, "    FMField \"%s\" \"%s\" %d %d\n",
		list[index].field_name, list[index].field_type,
		list[index].field_size, list[index].field_offset);
    }
    return str;
}

/*static char *
add_IOformat_to_string(char *str, IOFormat ioformat)
{
    return add_IOfieldlist_to_string(str, name_of_IOformat(ioformat),
				     field_list_of_IOformat(ioformat));
}*/

static char *
get_str(char *str, const char **name_p)
{
    int name_len = 0;
    char *name = malloc(1);
    while (*str != '"') {
	name = realloc(name, (name_len + 2));
	name[name_len++] = *(str++);
    }
    name[name_len] = 0;
    str++;
    *name_p = name;
    return str;
}

static char *
parse_FMformat_from_string(char *str, FMStructDescRec *f)
{
    char *name;
    FMFieldList list;
    int struct_size;
    f->format_name = NULL;
    f->field_list = NULL;
    f->struct_size = 0;
    f->opt_info = NULL;
    if (strncmp(str, "FMFormat \"", 10) == 0) {
	int field_count;
	int index = 0;
	str += 10;
	str = get_str(str, (const char **)&name);
	str += strlen(" StructSize ");
	if (sscanf(str, "%d", &struct_size) == 1) {
	    while(isdigit((int)*str)) str++;
	}
	str += strlen(" FieldCount ");
	if (sscanf(str, "%d", &field_count) == 1) {
	    while(isdigit((int)*str)) str++;
	}
	str++;
	list = malloc(sizeof(*list) * (field_count + 1));
	for (index = 0; index < field_count; index++) {
	    str += strlen("    FMField \"");
	    str = get_str(str, &(list[index].field_name));
	    str += 2;
	    str = get_str(str, &(list[index].field_type));
	    str++;
	    if (sscanf(str, "%d", &list[index].field_size) == 1) {
		while(isdigit((int)*str)) str++;
	    }
	    str++;
	    if (sscanf(str, "%d", &list[index].field_offset) == 1) {
		while(isdigit((int)*str)) str++;
	    }
	    str = strchr(str, '\n') + 1;
	}
	list[field_count].field_name = NULL;
	list[field_count].field_type = NULL;
	list[field_count].field_size = 0;
	list[field_count].field_offset = 0;
	if (field_count == 0) {
	    free(list);
	    list = NULL;
	}
	f->format_name = name;
	f->field_list = list;
	f->struct_size = struct_size;
    }
    return str;
}

void *
install_response_handler(CManager cm, int stone_id, char *response_spec,
			 void *local_data, FMFormat **ref_ptr)
{
    char *str = response_spec;
    (void)stone_id;
    if (strncmp("Terminal Action", str, strlen("Terminal Action")) == 0) {
	int format_count, i;
	FMStructDescList list;
	str += strlen("Terminal Action") + 1;
	sscanf(str, "  Format Count %d\n", &format_count);
	str = strchr(str, '\n') + 1;
	list = malloc(sizeof(list[0]) * (format_count + 1));
	for (i=0; i < format_count; i++) {
	    str = parse_FMformat_from_string(str, &list[i]);
	}
	list[format_count].format_name = NULL;
/*	INT_EVassoc_terminal_action(cm, stone_id, list, local_data, NULL);*/
    }
    if (strncmp("Filter Action", str, strlen("Filter Action")) == 0) {
	struct response_spec *response = malloc(sizeof(struct response_spec));
	int format_count, i;
	char *function;
	FMStructDescList list;
	str += strlen("Filter Action") + 1;
	sscanf(str, "  Format Count %d\n", &format_count);
	str = strchr(str, '\n') + 1;
	list = malloc(sizeof(list[0]) * (format_count + 1));
	for (i=0; i < format_count; i++) {
	    str = parse_FMformat_from_string(str, &list[i]);
	}
	list[format_count].format_name = NULL;
	function = malloc(strlen(str) + 1);
	strcpy(function, str);
	response->response_type = Response_Filter;
	response->u.filter.format_list = list;
	response->u.filter.function = function;
	response->u.filter.client_data = local_data;
	response->u.filter.reference_format =
	    EVregister_format_set(cm, list);
	if (ref_ptr) {
	    FMFormat *formats = malloc(2*sizeof(FMFormat));
	    formats[1] = NULL;
	    formats[0] = response->u.filter.reference_format;

	    *ref_ptr = formats;
	}
	return (void*)response;
    }
    if (strncmp("Router Action", str, strlen("Router Action")) == 0) {
	struct response_spec *response = malloc(sizeof(struct response_spec));
	int format_count, i;
	char *function;
	FMStructDescList list;
	str += strlen("Router Action") + 1;
	sscanf(str, "  Format Count %d\n", &format_count);
	str = strchr(str, '\n') + 1;
	list = malloc(sizeof(list[0]) * (format_count + 1));
	for (i=0; i < format_count; i++) {
	    str = parse_FMformat_from_string(str, &list[i]);
	}
	list[format_count].format_name = NULL;
	function = malloc(strlen(str) + 1);
	strcpy(function, str);
	response->response_type = Response_Router;
	response->u.filter.format_list = list;
	response->u.filter.function = function;
	response->u.filter.client_data = local_data;
	response->u.filter.reference_format =
	    EVregister_format_set(cm, list);
	if (ref_ptr) {
	    FMFormat *formats = malloc(2*sizeof(FMFormat));
	    formats[1] = NULL;
	    formats[0] = response->u.filter.reference_format;

	    *ref_ptr = formats;
	}
	return (void*)response;
    }
    if (strncmp("Transform Action", str, strlen("Transform Action")) == 0) {
	struct response_spec *response = malloc(sizeof(struct response_spec));
	int format_count, i;
	char *function;
	FMStructDescList in_list, out_list;
	str += strlen("Transform Action") + 1;
	sscanf(str, "  Input Format Count %d\n", &format_count);
	str = strchr(str, '\n') + 1;
	in_list = malloc(sizeof(in_list[0]) * (format_count + 1));
	for (i=0; i < format_count; i++) {
	    str = parse_FMformat_from_string(str, &in_list[i]);
	}
	in_list[format_count].format_name = NULL;
	in_list[format_count].field_list = NULL;
	if (sscanf(str, "  Output Format Count %d\n", &format_count) != 1) {
	    printf("output format parse failed\n");
	    return 0;
	}
	str = strchr(str, '\n') + 1;
	out_list = malloc(sizeof(out_list[0]) * (format_count + 1));
	for (i=0; i < format_count; i++) {
	    str = parse_FMformat_from_string(str, &out_list[i]);
	}
	out_list[format_count].format_name = NULL;
	out_list[format_count].field_list = NULL;
	function = malloc(strlen(str) + 1);
	strcpy(function, str);
	response->response_type = Response_Transform;
	response->u.transform.in_format_list = in_list;
	response->u.transform.out_format_list = out_list;
	response->u.transform.function = function;
	response->u.transform.client_data = local_data;
	response->u.transform.reference_input_format = NULL;
	if (in_list[0].format_name != NULL)
	    response->u.transform.reference_input_format =
		EVregister_format_set(cm, in_list);
	if (ref_ptr) {
	    FMFormat *formats = malloc(2*sizeof(FMFormat));
	    formats[1] = NULL;
	    formats[0] = response->u.transform.reference_input_format;
	    *ref_ptr = formats;
	}
	if (out_list[0].format_name != NULL)
	    response->u.transform.reference_output_format =
		EVregister_format_set(cm, out_list);
	response->u.transform.output_base_struct_size = out_list[0].struct_size;
	return (void*)response;
    }
    if (strncmp("Multityped Action", str, strlen("Multityped Action")) == 0) {
	struct response_spec *response = malloc(sizeof(struct response_spec));
	int list_count, j;
	char *function;
	FMStructDescList *struct_list;
	int accept_anonymous = 0;

	str += strlen("Multityped Action") + 1;
	sscanf(str, "  List Count %d\n", &list_count);
	str = strchr(str, '\n') + 1;
	struct_list = malloc(sizeof(struct_list[0]) * (list_count + 1));
	for (j = 0; j < list_count; j++) {
	    int format_count2, k;
	    FMStructDescList in_list;
	    sscanf(str, "Next format   Subformat Count %d\n", &format_count2);
	    str = strchr(str, '\n') + 1;

	    in_list = malloc(sizeof(in_list[0]) * (format_count2 + 1));
	    for (k=0; k < format_count2; k++) {
		str = parse_FMformat_from_string(str, &in_list[k]);
	    }
	    in_list[format_count2].format_name = NULL;
	    in_list[format_count2].field_list = NULL;
	    struct_list[j] = in_list;
	    if (struct_list[j]->field_list == NULL) {  /* anonymous */
		free((void*)struct_list[j]->format_name);
		free(in_list);
		struct_list[j] = NULL;
		list_count--;
		j--;
		accept_anonymous++;
	    }
	}
	struct_list[list_count] = NULL;
	function = malloc(strlen(str) + 1);
	strcpy(function, str);
	response->response_type = Response_Multityped;
	response->u.multityped.struct_list = struct_list;
	response->u.multityped.function = function;
	response->u.multityped.client_data = local_data;
	response->u.multityped.accept_anonymous = accept_anonymous;
	response->u.multityped.reference_input_format_list =
	    malloc((list_count +1) * sizeof(FMFormat));
	for (j = 0; j < list_count; j++) {
	    if ((struct_list[j])[0].format_name != NULL) {
		response->u.multityped.reference_input_format_list[j] =
		    EVregister_format_set(cm, struct_list[j]);
	    }
	}
	if (ref_ptr) {
	    FMFormat *formats = malloc((list_count + 1)*sizeof(FMFormat));
	    int k = 0;
	    for (k=0; k < list_count; k++) {
		formats[k] = response->u.multityped.reference_input_format_list[k];
	    }
	    formats[list_count] = NULL;
	    *ref_ptr = formats;
	}
	return (void*)response;
    }
    printf("Unparsed action : %s\n", str);
    return NULL;
}


char *
create_terminal_action_spec(FMStructDescList format_list)
{
    int format_count = 0;
    int i;
    char *str;
    while(format_list[format_count].format_name != NULL) format_count++;
    str = malloc(50);
    sprintf(str, "Terminal Action   Format Count %d\n", format_count);

    for (i = 0 ; i < format_count; i++) {
	str = add_FMfieldlist_to_string(str, &format_list[i]);
    }
    return str;
}

char *
INT_create_bridge_action_spec(int stone_id, char *contact)
{
    int size = (int) strlen(contact);
    char *output;
    size += (int) strlen("Bridge Action") + 20;
    output = malloc(size);
    sprintf(output, "Bridge Action %d %s", stone_id, contact);
    return output;
}

void
parse_bridge_action_spec(char *action_spec, int *target, char **contact)
{
    action_spec += strlen("Bridge Action ");
    sscanf(action_spec, "%d", target);
    while(*action_spec != ' ') action_spec++;
    action_spec++;
    *contact = action_spec;
}

action_value
action_type(char *action_spec)
{
    if (action_spec == NULL)
	return Action_Split;
    if (strncmp(action_spec, "Bridge Action", 13) == 0)
	return Action_Bridge;
    if (strncmp(action_spec, "Filter Action", 13) == 0)
	return Action_Immediate;
    if (strncmp(action_spec, "Router Action", 13) == 0)
	return Action_Immediate;
    if (strncmp(action_spec, "Transform Action", 16) == 0)
	return Action_Immediate;
    if (strncmp(action_spec, "Multityped Action", 17) == 0)
	return Action_Multi;
    if (strncmp(action_spec, "sink:", 5) == 0)
	return Action_Terminal;
    if (strncmp(action_spec, "source:", 7) == 0)
	return Action_Source;
    if (strncmp(action_spec, "Split Action", 7) == 0)
	return Action_Split;
    return Action_NoAction;
}

char *
INT_create_filter_action_spec(FMStructDescList format_list, char *function)
{
    int format_count = 0;
    int i;
    char *str;
    while(format_list && (format_list[format_count].format_name != NULL)) format_count++;
    str = malloc(50);
    sprintf(str, "Filter Action   Format Count %d\n", format_count);

    for (i = 0 ; i < format_count; i++) {
	str = add_FMfieldlist_to_string(str, &format_list[i]);
    }
    str = realloc(str, strlen(str) + strlen(function) + 1);
    strcpy(&str[strlen(str)], function);
    return str;
}

char *
INT_create_router_action_spec(FMStructDescList format_list, char *function)
{
    int format_count = 0;
    int i;
    char *str;
    while(format_list && (format_list[format_count].format_name != NULL)) format_count++;
    str = malloc(50);
    sprintf(str, "Router Action   Format Count %d\n", format_count);

    for (i = 0 ; i < format_count; i++) {
	str = add_FMfieldlist_to_string(str, &format_list[i]);
    }
    str = realloc(str, strlen(str) + strlen(function) + 1);
    strcpy(&str[strlen(str)], function);
    return str;
}

char *
INT_create_transform_action_spec(FMStructDescList format_list, FMStructDescList out_format_list, char *function)
{
    int format_count = 0;
    int i;
    char *str;
    while(format_list && format_list[format_count].format_name != NULL)
	format_count++;
    str = malloc(50);
    sprintf(str, "Transform Action   Input Format Count %d\n", format_count);

    for (i = 0 ; i < format_count; i++) {
	str = add_FMfieldlist_to_string(str, &format_list[i]);
    }

    format_count = 0;
    while(out_format_list[format_count].format_name != NULL) format_count++;
    str = realloc(str, strlen(str) + 30);
    sprintf(str + strlen(str), "  Output Format Count %d\n", format_count);

    for (i = 0 ; i < format_count; i++) {
	str = add_FMfieldlist_to_string(str, &out_format_list[i]);
    }
    str = realloc(str, strlen(str) + strlen(function) + 1);
    strcpy(&str[strlen(str)], function);
    return str;
}

extern char *
INT_create_multityped_action_spec(FMStructDescList *input_format_lists, char *function)
{
    int list_count = 0;
    int l, i;
    char *str;
    while(input_format_lists && input_format_lists[list_count] != NULL)
	list_count++;

    str = malloc(50);
    sprintf(str, "Multityped Action   List Count %d\n", list_count);

    for (l = 0; l < list_count; l++) {
	int format_count = 0;
	FMStructDescList format_list = input_format_lists[l];
	while(format_list && format_list[format_count].format_name != NULL)
	    format_count++;
	str = realloc(str, strlen(str) + 50);
	sprintf(str + strlen(str), "Next format   Subformat Count %d\n",
		format_count);
	for (i = 0 ; i < format_count; i++) {
	    str = add_FMfieldlist_to_string(str, &format_list[i]);
	}
    }

    str = realloc(str, strlen(str) + strlen(function) + 1);
    strcpy(&str[strlen(str)], function);
    return str;
}

struct ev_state_data {
    CManager cm;
    struct _event_item *cur_event;
    int stone;
    int proto_action_id;
    int out_count;
    int *out_stones;
    queue_item *item;
    struct _queue *queue;
    response_instance instance;
    int did_output;
};

extern CManager
get_cm_from_ev_state(void *vevstate)
{
    struct ev_state_data *evstate = (struct ev_state_data *)vevstate;
    return evstate->cm;
}

static int
filter_wrapper(CManager cm, struct _event_item *event, void *client_data,
	       attr_list attrs, int out_count, int *out_stones)
{
    response_instance instance = (response_instance)client_data;
    int ret;
    cod_exec_context ec = instance->u.filter.ec;
    struct ev_state_data ev_state;

    ev_state.cm = cm;
    ev_state.cur_event = event;
    ev_state.out_count = out_count;
    ev_state.out_stones = out_stones;
    if (ec != NULL) {
#ifdef HAVE_COD_H
	cod_assoc_client_data(ec, 0x34567890, (intptr_t)&ev_state);

	ret = ((int(*)(cod_exec_context, void *, attr_list))instance->u.filter.code->func)(ec, event->decoded_event, attrs);
#endif
    } else {
	/* DLL-based handler */
	ret = ((int(*)(void *, attr_list))instance->u.filter.func_ptr)(event->decoded_event, attrs);
    }
    if (ret) {
	CMtrace_out(cm, EVerbose, "Filter function returned %d, submitting further to stone %d\n", ret, out_stones[0]);
	internal_path_submit(cm, out_stones[0], event);
    } else {
	CMtrace_out(cm, EVerbose, "Filter function returned %d, NOT submitting\n", ret);
    }
    return ret;
}
static int
router_wrapper(CManager cm, struct _event_item *event, void *client_data,
	       attr_list attrs, int out_count, int *out_stones)
{
    response_instance instance = (response_instance)client_data;
    int ret;
    if (instance->u.filter.func_ptr) {
	ret = ((int(*)(void *, attr_list))instance->u.filter.func_ptr)(event->decoded_event, attrs);
    } else {
#ifdef HAVE_COD_H
	int (*func)(cod_exec_context, void *, attr_list) =
	    (int(*)(cod_exec_context, void *, attr_list))instance->u.filter.code->func;
	cod_exec_context ec = instance->u.filter.ec;
	struct ev_state_data ev_state;

	ev_state.cm = cm;
	ev_state.cur_event = event;
	ev_state.out_count = out_count;
	ev_state.out_stones = out_stones;
	cod_assoc_client_data(ec, 0x34567890, (intptr_t)&ev_state);
	ret = (func)(ec, event->decoded_event, attrs);
#endif
    }
    if (ret >= 0) {
	if (ret >= out_count) {
	    CMtrace_out(cm, EVerbose, "Router function returned %d, larger than the number of associated outputs\n", ret);
	} else if (out_stones[ret] == -1) {
	    CMtrace_out(cm, EVerbose, "Router function returned %d, which has not been set with EVaction_set_output()\n", ret);
	} else {
	    CMtrace_out(cm, EVerbose, "Router function returned %d, submitting further to stone %d\n", ret, out_stones[ret]);
	    internal_path_submit(cm, out_stones[ret], event);
	}
    } else {
	CMtrace_out(cm, EVerbose, "Router function returned %d, NOT submitting\n", ret);
    }
    return ret;
}

static void
transform_free_wrapper(void *data, void *free_data)
{
    FMFormat out_format = (FMFormat)free_data;
    FMfree_var_rec_elements(out_format, data);
    free(data);
}

static int
transform_wrapper(CManager cm, struct _event_item *event, void *client_data,
		  attr_list attrs, int out_count, int *out_stones)
{
    response_instance instance = (response_instance)client_data;
    int ret;
    void *out_event = malloc(instance->u.transform.out_size);
    int(*func)(cod_exec_context, void *, void*, attr_list, attr_list) = NULL;
    cod_exec_context ec = instance->u.transform.ec;
    struct ev_state_data ev_state;
    attr_list output_attrs = create_attr_list();

    ev_state.cm = cm;
    ev_state.cur_event = event;
    ev_state.stone = instance->stone;
    ev_state.proto_action_id = instance->proto_action_id;
    ev_state.out_count = out_count;
    ev_state.out_stones = out_stones;

    if (CMtrace_on(cm, EVerbose)) {
	fprintf(cm->CMTrace_file, "Input Transform Event is :\n");
	if (event->reference_format) {
	    FMfdump_data(cm->CMTrace_file, event->reference_format, event->decoded_event, 10240);
	} else {
	    fprintf(cm->CMTrace_file, "       ****  UNFORMATTED  ****\n");
	}
    }
    memset(out_event, 0, instance->u.transform.out_size);
    if (ec != NULL) {
#ifdef HAVE_COD_H
	func = (int(*)(cod_exec_context, void *, void*, attr_list, attr_list))instance->u.transform.code->func;
	cod_assoc_client_data(ec, 0x34567890, (intptr_t)&ev_state);
	ret = func(ec, event->decoded_event, out_event, attrs, output_attrs);
#endif
    } else {
	/* DLL-based handler */
	ret = ((int(*)(void *, void *, attr_list, attr_list))instance->u.transform.func_ptr)(event->decoded_event, out_event, attrs, output_attrs);
    }

    if (ret && (out_stones[0] == -1)) {
	printf("Transform output stone ID not set, event discarded\n");
	ret = 0;
    }
    if (ret) {
	struct _EVSource s;
	if (CMtrace_on(cm, EVerbose)) {
	    FMFormat f = instance->u.transform.out_format;
	    fprintf(cm->CMTrace_file, " Transform function returned %d, submitting further\n", ret);
	    FMfdump_data(cm->CMTrace_file, f, out_event, 10240);
	}
	s.local_stone_id = out_stones[0];
	s.cm = cm;
	s.format = NULL;
	s.reference_format = instance->u.transform.out_format;
	s.free_func = transform_free_wrapper;
	s.free_data = instance->u.transform.out_format;
	s.preencoded = 0;
	INT_EVsubmit(&s, out_event, output_attrs);
    } else {
	CMtrace_out(cm, EVerbose, "Transform function returned %d, NOT submitting\n", ret);
	transform_free_wrapper(out_event, instance->u.transform.out_format);
    }
    free_attr_list(output_attrs);
    return ret;
}

/* {{{ cod_find_index */
static queue_item *queue_find_index(queue_item *item, int i, FMFormat format) {
    for (;;) {
        if (!item) {
            return NULL;
        }
        if (!format || (item->item->reference_format == format)) {
            if (i == 0)
                return item;
            --i;
        }
        item = item->next;
    }
}

static queue_item *queue_find_anonymous(queue_item *item, int i, FMFormat *formats) {
    for (;;) {
	int known = 0;
	int j = 0;
        if (!item) {
            return NULL;
        }
	while(formats[j]) {
	    if (item->item->reference_format == formats[j]) known++;
	    j++;
	}
        if (known == 0) {
            if (i == 0)
                return item;
            --i;
        }
        item = item->next;
    }
}

static queue_item *cod_find_index_rel(struct ev_state_data *ev_state, int queue, int index)
{
    if (queue != -2) {
	return queue_find_index(
	    ev_state->queue->queue_head, index,
	    queue < 0 ?  NULL : ev_state->instance->u.queued.formats[queue]);
    } else {
	return queue_find_anonymous(
	    ev_state->queue->queue_head, index, ev_state->instance->u.queued.formats);
    }
}

static queue_item *cod_find_index_abs(struct ev_state_data *ev_state, int queue, int index) {
    queue_item *ret;
    ret = queue_find_index(ev_state->queue->queue_head, index, NULL);
    if (!ret)
        return NULL;
    if (queue < 0 || ret->item->reference_format ==
            ev_state->instance->u.queued.formats[queue])
        return ret;
    else
        return NULL;
}

static queue_item *cod_find_index(int absp, struct ev_state_data *ev_state, int queue, int index) {
    if (absp)
        return cod_find_index_abs(ev_state, queue, index);
    else
        return cod_find_index_rel(ev_state, queue, index);
}

/* }}} */

static void cod_ev_discard(cod_exec_context ec, int absp, int queue, int index) 
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    CManager cm = ev_state->cm;
    queue_item *item;

    item = cod_find_index(absp, ev_state, queue, index);

    assert(item);

    EVdiscard_queue_item(cm, ev_state->stone, item);
}

static void cod_ev_discard_rel(cod_exec_context ec, int queue, int index) {
    cod_ev_discard(ec, 0, queue, index);
}

#ifdef NOT_DEF
static void cod_ev_discard_abs(cod_exec_context ec, int queue, int index) {
    cod_ev_discard(ec, 1, queue, index);
}
#endif

static EVstone
port_to_stone(struct ev_state_data *evstate, int port)
{
    if (port >= evstate->out_count) {
	fprintf(stderr, "Stone has %d outbound ports, port %d invalid\n",
		evstate->out_count, port);
	return -1;
    }
    if (evstate->out_stones[port] == -1) {
	fprintf(stderr, "Stone port %d target has not been set\n",
		port);
    }
    return evstate->out_stones[port];
}

static void cod_ev_discard_and_submit(cod_exec_context ec,
        int absp, int port, int queue, int index) 
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    CManager cm = ev_state->cm;
    queue_item *item;
    EVstone target_stone = port_to_stone(ev_state, port);

    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return;
    }

    item = cod_find_index(absp, ev_state, queue, index);

    if (item == NULL) {
        printf("Item %x not found on queue %d, stone %d\n", index, queue, ev_state->stone);
	return;
    }

    item->handled = 0;

    internal_path_submit(cm, target_stone, item->item);

    ev_state->did_output++;
    EVdiscard_queue_item(cm, ev_state->stone, item);
}

static void cod_ev_submit(cod_exec_context ec,
        int absp, int port, int queue, int index) 
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    CManager cm = ev_state->cm;
    queue_item *item;
    EVstone target_stone = port_to_stone(ev_state, port);

    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return;
    }

    item = cod_find_index(absp, ev_state, queue, index);

    if (item == NULL) {
        printf("Item %x not found on queue %d, stone %d\n", index, queue, ev_state->stone);
	return;
    }

    item->handled = 0;

    internal_path_submit(cm, target_stone, item->item);

    ev_state->did_output++;
}

static void cod_ev_submit_rel(cod_exec_context ec,int port, int queue, int index) 
{
    cod_ev_submit(ec, 0, port, queue, index);
}

#ifdef NOT_DEF
static void cod_ev_submit_abs(cod_exec_context ec,int port, int queue, int index) 
{
    cod_ev_submit(ec, 1, port, queue, index);
}
#endif

static int cod_ev_get_port(cod_exec_context ec, int queue)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    int port = (ev_state->out_stones[queue]);
    
    return  port;
}

static int cod_ev_target_size(cod_exec_context ec, int stone_num)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);

    CManager cm = ev_state->cm;
    stone_type stone = stone_struct(cm->evp, stone_num);
    if (!stone) return -1;
    return stone->queue_size;
}



static void cod_ev_discard_and_submit_rel(cod_exec_context ec, int port, int queue,
        int index) {
    struct ev_state_data *ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    EVstone target_stone = port_to_stone(ev_state, port);
    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return;
    }

    cod_ev_discard_and_submit(ec, 0, port, queue, index);
}

#ifdef NOT_DEF
static void cod_ev_discard_and_submit_abs(cod_exec_context ec, int port, int queue,
        int index) {
    struct ev_state_data *ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    EVstone target_stone = port_to_stone(ev_state, port);
    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return;
    }

    cod_ev_discard_and_submit(ec, 1, port, queue, index);
}
#endif

static void *cod_ev_get_data(cod_exec_context ec, int absp, int queue, int index)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    queue_item *item;
    item = cod_find_index(absp, ev_state, queue, index);

    if (!item) {
	return NULL;
    }
    assert(item->item);

    if (!item->item->decoded_event) {
        item->item = cod_decode_event(ev_state->cm, ev_state->stone,
				      ev_state->proto_action_id, item->item);
    }
    assert(item->item->decoded_event);

    return item->item->decoded_event;
}

static void *cod_ev_get_data_rel(cod_exec_context ec, int queue, int index) {
    return cod_ev_get_data(ec, 0, queue, index);
}

static void *cod_ev_get_data_abs(cod_exec_context ec, int queue, int index) {
    return cod_ev_get_data(ec, 1, queue, index);
}

static int cod_ev_conforms(cod_exec_context ec, int queue, int index) {
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    return cod_find_index_abs(ev_state, queue, index) != NULL;
}

static int cod_ev_present(cod_exec_context ec, int queue, int index) {
    struct ev_state_data *ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    return cod_find_index_rel(ev_state, queue, index) != NULL;
}

static int cod_ev_count(cod_exec_context ec, int queue) {
    struct ev_state_data *ev_state;
    FMFormat type;
    queue_item *item;
    int count = 0;

    /*    queue == -1 RETURNS total event count */
    /*    queue == -2 returns anonymous event count (I.E. count of events not in the queued format list) */

    int format_count = 0;

    ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    while(ev_state->instance->u.queued.formats[format_count]) format_count++;

    if (format_count <= queue) {
	printf("Error, queue parameter(%d) to EVCount is larger than queue count (%d)\n",
	       queue, format_count);
	return -1;
    }
    if (queue == -2) {
	item = ev_state->queue->queue_head;
	while (item) {
	    int i;
	    for (i =0; i < format_count; i++) {
		/* on match break out of loop */
		if (item->item->reference_format == ev_state->instance->u.queued.formats[i]) break;
	    }
	    /* if we got to format_count without matching anything, increment count */
	    if (i == format_count) ++count;
	    item = item->next;
	}
	return count;
    }
    type = queue < 0 ? NULL :
        ev_state->instance->u.queued.formats[queue];
    item = ev_state->queue->queue_head;
    while (item) {
        if (!type || item->item->reference_format == type)
            ++count;
        item = item->next;
    }

    return count;
}

static attr_list cod_ev_get_attrs(cod_exec_context ec, int queue, int index) {
    struct ev_state_data *ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    queue_item *item = cod_find_index_rel(ev_state, queue, index);
    attr_list *pattr;

    if (NULL == item) {
	printf("No item at index %d on queue %d\n", index, queue);

	return NULL;
    }
    pattr = &item->item->attrs;
    if (!*pattr) {
        *pattr = CMcreate_attr_list(ev_state->cm);
    }
    return *pattr;
}

static attr_list cod_ev_get_stone_attrs(cod_exec_context ec, char *stone_name) {
    struct ev_state_data *ev_state = (void*) cod_get_client_data(ec, 0x34567890);
    CManager cm = ev_state->cm;
    event_path_data evp = cm->evp;
    attr_list ret_list = NULL;
    int cur_stone;
    static atom_t STONE_NAME_ATOM = -1;
    if (STONE_NAME_ATOM == -1) {
	STONE_NAME_ATOM = attr_atom_from_string("EVP_STONE_NAME");
    }
    for (cur_stone = evp->stone_base_num; cur_stone < evp->stone_count + evp->stone_base_num; ++cur_stone) {
	stone_type stone = stone_struct(evp, cur_stone);
	if (stone && (stone->stone_attrs != NULL)) {
	    char *this_stone_name = NULL;
	    if (get_string_attr(stone->stone_attrs, STONE_NAME_ATOM, &this_stone_name)) {
		if (stone_name && (strcmp(this_stone_name, stone_name) == 0)) {
		    if (ret_list) printf("Warning, duplicate stone name \"%s\" found during attr query\n", stone_name);
		    ret_list = stone->stone_attrs;
		}
	    }
	}
    }
    return ret_list;
}

static int
queued_wrapper(CManager cm, struct _queue *queue, queue_item *item,
                void *client_data, int out_count, int *out_stones)
{
#ifdef HAVE_COD_H
    response_instance instance = (response_instance)client_data;
    int(*func)(cod_exec_context) =  /* XXX wrong type */
	(int(*)(cod_exec_context))instance->u.queued.code->func;
    cod_exec_context ec = instance->u.queued.ec;
    struct ev_state_data ev_state;

    ev_state.cm = cm;
    ev_state.cur_event = NULL;
    ev_state.stone = instance->stone;
    ev_state.proto_action_id = instance->proto_action_id;
    ev_state.out_count = out_count;
    ev_state.out_stones = out_stones;
    ev_state.queue = queue;
    ev_state.item = item;
    ev_state.instance = instance;
    ev_state.did_output = 0;
    cod_assoc_client_data(ec, 0x34567890, (intptr_t)&ev_state);

    func(ec);

    return ev_state.did_output;
#else
    return 0;
#endif
}

static response_instance
generate_filter_code(CManager cm, struct response_spec *mrd, stone_type stone,
		     FMFormat format);
static response_instance
generate_multityped_code(CManager cm, struct response_spec *mrd, stone_type stone,
			  FMFormat *formats);

extern
void
free_struct_list(FMStructDescList list)
{
    int format_count = 0;
    int format;

    while(list[format_count].format_name != NULL) format_count++;

    for (format = 0; format < format_count; format++) {
	free((void*)list[format].format_name);
	free_FMfield_list(list[format].field_list);
    }
    free(list);
}

static FMFormat
localize_format(CManager cm, FMFormat format)
{
    FMFormat ret;
    FMStructDescList local_formats = get_localized_formats(format);
    ret = EVregister_format_set(cm, local_formats);
    free_struct_list(local_formats);
    return ret;
}

void
dump_mrd(void *mrdv)
{
    struct response_spec *mrd = (struct response_spec *) mrdv;
    switch (mrd->response_type) {
    case Response_Filter:
	printf("Response Filter, code is %s\n",
	       mrd->u.filter.function);
	break;
    case Response_Router:
	printf("Response Router, code is %s\n",
	       mrd->u.filter.function);
	break;
    case Response_Transform:
	printf("Response Transform, code is %s\n",
	       mrd->u.transform.function);
	break;
    case Response_Multityped:
	printf("Multityped Action, code is %s\n",
	       mrd->u.transform.function);
	break;
    }
}

static int
proto_action_in_stage(proto_action *act, action_class stage) {
    switch (stage) {
    case Immediate_and_Multi:
        if (act->action_type == Action_Multi) return 1;
        /* fallthrough */
    case Immediate:
        switch (act->action_type) {
        case Action_Terminal:
        case Action_Filter:
        case Action_Split:
        case Action_Immediate:
        case Action_Store:
            return 1;
        default:
            return 0;
        }
    case Bridge:
        return act->action_type == Action_Bridge;
    case Congestion:
        return act->action_type == Action_Congestion;
    default:
        assert(0);
    }
    return 0;
}

extern int
FMformat_compat_cmp2(FMFormat format, FMFormat *formatList,
		     int listSize, FMcompat_formats * older_format);

static void
free_multi_response(void *client_data)
{
    response_instance resp = (response_instance) client_data;
    resp->u.queued.ref_count--;
    if (resp->u.queued.ref_count != 0) return;
    if (resp->u.queued.code) cod_code_free(resp->u.queued.code);
    if (resp->u.queued.ec) cod_exec_context_free(resp->u.queued.ec);
    free(resp);
}

static void
free_imm_response(void *client_data)
{
    response_instance resp = (response_instance) client_data;
    switch (resp->response_type) {
    case Response_Filter:
    case Response_Router:
	if (resp->u.filter.code) cod_code_free(resp->u.filter.code);
	if (resp->u.filter.ec) cod_exec_context_free(resp->u.filter.ec);
	break;
    case Response_Transform:
	if (resp->u.transform.code) cod_code_free(resp->u.transform.code);
	if (resp->u.transform.ec) cod_exec_context_free(resp->u.transform.ec);
	break;
    default:
	break;
    }
    free(resp);
}

int
response_determination(CManager cm, stone_type stone, action_class stage, event_item *event)
{
    int nearest_proto_action = -1;
    int return_value = 0;
    FMFormat conversion_target_format = NULL;
    FMFormat matching_format = NULL;
    int i, format_count = 0;
    FMFormat * formatList;
    int *format_map;
    FMcompat_formats older_format = NULL;

    formatList =
	(FMFormat *) malloc((stone->proto_action_count + 1) * sizeof(FMFormat));
    format_map = (int *) malloc((stone->proto_action_count + 1) * sizeof(int));
    for (i = 0; i < stone->proto_action_count; i++) {
	int j = 0;
        if (!proto_action_in_stage(&stone->proto_actions[i], stage)) {
            continue;
        }
	while (stone->proto_actions[i].matching_reference_formats &&
	       (stone->proto_actions[i].matching_reference_formats[j] != NULL)) {
	    if (strcmp(name_of_FMformat(event->reference_format), name_of_FMformat(stone->proto_actions[i].matching_reference_formats[j])) == 0 ) {
		formatList = (FMFormat *) realloc(formatList, (format_count + 2) * sizeof(FMFormat));
		format_map = realloc(format_map, (format_count + 2) * sizeof(int));
		formatList[format_count] = stone->proto_actions[i].matching_reference_formats[j];
		format_map[format_count] = i;
		format_count++;
	    }
	    j++;
	}
    }
    formatList[format_count] = NULL;
    if (event->reference_format == NULL) {
	/* special case for unformatted input */
	for (i=0 ; i < stone->proto_action_count ; i++) {
            if (!proto_action_in_stage(&stone->proto_actions[i], stage))
		continue;
	    if ((stone->proto_actions[i].matching_reference_formats == NULL) ||
		(stone->proto_actions[i].matching_reference_formats[0] == NULL))
		nearest_proto_action = i;
	}
    } else {
	int map_entry = FMformat_compat_cmp2(event->reference_format,
						    formatList,
						    format_count,
						    &older_format);
	if (map_entry != -1) {
            nearest_proto_action = format_map[map_entry];
            matching_format = formatList[map_entry];
        }
    }
    if (nearest_proto_action == -1) {
        /* special case for accepting anything */
        for (i=0; i < stone->proto_action_count; i++) {
            if (!proto_action_in_stage(&stone->proto_actions[i], stage)) continue;
            if (((stone->proto_actions[i].matching_reference_formats == NULL) ||
		 (stone->proto_actions[i].matching_reference_formats[0] == NULL))
                && stone->proto_actions[i].data_state != Requires_Decoded) {
                nearest_proto_action = i;
            }
            if (stone->proto_actions[i].action_type == Action_Multi) {
		struct response_spec *mrd;

		mrd =
		    stone->proto_actions[i].o.imm.mutable_response_data;
		if (mrd->u.multityped.accept_anonymous) {
		    nearest_proto_action = i;
		}
            }
        }
    }
    free(formatList);
    free(format_map);
    if (nearest_proto_action != -1) {
	int action_generated = 0;
	proto_action *proto = &stone->proto_actions[nearest_proto_action];
	if (proto->action_type == Action_Immediate) {
	    /* must be immediate action */
	    response_instance instance;
	    struct response_spec *mrd;
	    mrd =
		proto->o.imm.mutable_response_data;
	    switch(mrd->response_type) {
	    case Response_Filter:
	    case Response_Router:
		if (event->event_encoded) {
		    conversion_target_format =
			localize_format(cm, event->reference_format);
		} else {
		    conversion_target_format = event->reference_format;
		}
		break;
	    case Response_Transform:
		conversion_target_format = mrd->u.transform.reference_input_format;
		break;
	    case Response_Multityped:
		assert(FALSE);
                break;
	    }

	    instance = generate_filter_code(cm, mrd, stone, conversion_target_format);
	    if (instance == NULL) return 0;
	    instance->stone = stone->local_id;
	    instance->proto_action_id = nearest_proto_action;
	    action_generated++;
	    switch(mrd->response_type) {
	    case Response_Filter:
		INT_EVassoc_mutated_imm_action(cm, stone->local_id, nearest_proto_action,
					       filter_wrapper, instance,
					       conversion_target_format, free_imm_response);
		break;
	    case Response_Router:
		INT_EVassoc_mutated_imm_action(cm, stone->local_id, nearest_proto_action,
					       router_wrapper, instance,
					       conversion_target_format, free_imm_response);
		break;
	    case Response_Transform:
		INT_EVassoc_mutated_imm_action(cm, stone->local_id, nearest_proto_action,
					       transform_wrapper, instance,
					       conversion_target_format, free_imm_response);
		break;
            default:
		assert(FALSE);
		break;
	    }
	    return_value = 1;
	} else 	if (proto->action_type == Action_Multi || proto->action_type == Action_Congestion) {
	    response_instance instance;
	    struct response_spec *mrd;

	    mrd =
		proto->o.imm.mutable_response_data;
	    instance = generate_multityped_code(cm, mrd, stone,
						 proto->matching_reference_formats);
	    if (instance == 0) {
                return 0;
            }
	    instance->stone = stone->local_id;
	    instance->proto_action_id = nearest_proto_action;
	    action_generated++;
	    INT_EVassoc_mutated_multi_action(cm, stone->local_id, nearest_proto_action,
					     queued_wrapper, instance,
					     proto->matching_reference_formats, free_multi_response);
	    if (mrd->u.multityped.accept_anonymous && (matching_format == NULL)) {
		/* we're accepting this as an anonymous target */
		INT_EVassoc_anon_multi_action(cm, stone->local_id, nearest_proto_action, queued_wrapper, instance,
					      event->reference_format);
	    }
            if (event->event_encoded) {
                conversion_target_format = matching_format;
            }
            return_value = 1;
	} else {
	    response_cache_element *resp;

	    conversion_target_format = NULL;
	    if (proto->matching_reference_formats) {
		conversion_target_format = proto->matching_reference_formats[0];
	    }

	    /* we'll install the conversion later, first map the response */
	    if (stone->response_cache_count == 0) {
		if (stone->response_cache != NULL) free(stone->response_cache);
		stone->response_cache = malloc(sizeof(stone->response_cache[0]));
	    } else {
		stone->response_cache =
		    realloc(stone->response_cache,
			    (stone->response_cache_count + 1) * sizeof(stone->response_cache[0]));
	    }
	    resp = &stone->response_cache[stone->response_cache_count++];
	    proto_action *proto2 = &stone->proto_actions[nearest_proto_action];
	    if (conversion_target_format) {
		resp->reference_format = conversion_target_format;
	    } else {
		resp->reference_format = event->reference_format;
	    }
	    resp->proto_action_id = nearest_proto_action;
	    resp->action_type = proto2->action_type;
	    resp->requires_decoded = (proto2->data_state == Requires_Decoded);
            resp->stage = stage;
	}
	if (conversion_target_format != NULL) {
	    if (event->event_encoded) {
		/* create a decode action */
		INT_EVassoc_conversion_action(cm, stone->local_id, stage,
					      conversion_target_format,
					      event->reference_format);
		return_value = 1;
	    } else {
		if (event->reference_format != conversion_target_format) {
		    /* 
		     * create a decode action anyway, the event will be
		     * encoded to a buffer and then decoded into the target
		     * format.  Doing this more efficiently is difficult. 
		     */
		    INT_EVassoc_conversion_action(cm, stone->local_id, stage,
						  conversion_target_format,
						  event->reference_format);
		    return_value = 1;
		} else {
		    return_value = 1;
		}
	    }
	} else {
            return_value = 1;
        }
    }
    fix_response_cache(stone);
    return return_value;
}

void
response_data_free(CManager cm, void *resp_void)
{
    struct response_spec *resp = (struct response_spec*)resp_void;
    switch(resp->response_type) {
    case Response_Filter:
    case Response_Router:
	free_struct_list(resp->u.filter.format_list);
        free(resp->u.filter.function);
	break;
    case Response_Transform:
	free_struct_list(resp->u.transform.in_format_list);
	free_struct_list(resp->u.transform.out_format_list);
        free(resp->u.transform.function);
	break;
    case Response_Multityped:
      {
	  int i = 0;
	  while(resp->u.multityped.struct_list[i] != NULL) {
	      FMStructDescList list = resp->u.multityped.struct_list[i];
	      int j = 0;
	      while (list[j].format_name != NULL) {
		  free((void*)list[j].format_name);
		  free_FMfield_list(list[j].field_list);
		  j++;
	      }
	      free(list);
	      i++;
	  }
      }
      free(resp->u.multityped.struct_list);
      free(resp->u.multityped.reference_input_format_list);
      free(resp->u.multityped.function);
      break;
    default:
	break;
    }
    free(resp);
}

#ifdef NOT_DEF
static void
cod_free_wrapper(void *data, void *free_data)
{
    event_item *event = (event_item *)free_data;
    FMfree_var_rec_elements(event->reference_format,
			    data);
}
#endif

extern void 
INT_EVadd_standard_routines(CManager cm, char *extern_string,
			    cod_extern_entry *externs)
{
    event_path_data evp = cm->evp;
    int count = 0;
    if (evp->externs == NULL) {
	evp->externs = malloc(sizeof(evp->externs[0]) * 2);
    } else {
	while(evp->externs[count].extern_decl != NULL) count++;
	evp->externs = realloc(evp->externs, 
			       sizeof(evp->externs[0]) * (count + 2));
    }
    evp->externs[count].extern_decl = extern_string;
    evp->externs[count].externs = externs;
    evp->externs[count + 1].extern_decl = NULL;
    evp->externs[count + 1].externs = NULL;
}

extern void 
INT_EVadd_standard_structs(CManager cm, FMStructDescList *lists)
{
    event_path_data evp = cm->evp;
    int count = 0, new = 0, i;
    
    while (lists[new] != NULL) new++;

    if (evp->extern_structs == NULL) {
	evp->extern_structs = malloc(sizeof(evp->extern_structs[0]) * (new+1));
    } else {
	while(evp->extern_structs[count] != NULL) count++;
	evp->extern_structs = realloc(evp->extern_structs, 
			       sizeof(evp->extern_structs[0]) * (count + new + 1));
    }
    for (i=0; i<= new; i++) {
	evp->extern_structs[count + i] = lists[i];
    }
}

static void
cod_ffs_write(cod_exec_context ec, FFSFile fname,  int queue, int index)
{
    FMFormat ref_format, file_format;
    struct ev_state_data * ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    queue_item * my_item = cod_find_index(0, ev_state, queue, index);
    FMContext fmc;
    attr_list *temp_attr;
    FMStructDescList format_list;

    if(!my_item) {
	fprintf(stderr, "No corresponding item in the queue\n");
	return;
    }

    ref_format = my_item->item->reference_format;
    fmc = FMContext_of_file(fname);

    format_list = format_list_of_FMFormat(ref_format);
    file_format = FMregister_data_format(fmc, format_list);
    
    temp_attr = &my_item->item->attrs;
    if(!*temp_attr) {
	printf("There is no attr for: %s\n", format_list->format_name);
    }
    
    if(my_item->item->event_encoded) {
	fprintf(stderr, "Event is encoded, have not handled this case.  Can not write to file\n");
	return;
    } else {
	void * temp_data = my_item->item->decoded_event;
	if(!write_FFSfile_attrs(fname, file_format, temp_data, *temp_attr))
	    fprintf(stderr, "Error in writing FFS_file!\n");
    }
    return;
}

static void
cod_ffs_read(cod_exec_context ec, FFSFile fname, void * data, attr_list * temp, int queue)
{
    FMFormat ref_format;
    struct ev_state_data * ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    FFSTypeHandle temp_type;
    FFSContext fmc = FFSContext_of_file(fname);
    FMStructDescList format_list;

    ref_format = ev_state->instance->u.queued.formats[queue];
    format_list = format_list_of_FMFormat(ref_format);
    temp_type = FFSset_fixed_target(fmc, format_list);
    (void)temp_type;
    FFSread_attr(fname, data, temp);

    return;
}


static
int
cod_max_output(cod_exec_context ec)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    return ev_state->out_count;
}

static int
cod_target_stone_on_port(cod_exec_context ec, int port, void *data, void *type_info, attr_list attrs)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    EVstone target_stone = port_to_stone(ev_state, port);

    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return -1;
    }
    return target_stone;
}

struct delayed_event {
    EVstone to_stone;
    event_item *event;
};

extern void do_local_actions(CManager cm);

static void
EVdelayed_submit_func(CManager cm, void* vdelayed)
{
    struct delayed_event *delayed = (struct delayed_event *)vdelayed;
    int stone_num = delayed->to_stone;
    event_item *event = delayed->event;
    free(delayed);
    CManager_lock(cm);
    internal_path_submit(cm, stone_num, event);
    do_local_actions(cm);
    return_event(cm->evp, event);
    CManager_unlock(cm);
}

static void
internal_cod_submit_general(cod_exec_context ec, int port, void *data, void *type_info, attr_list attrs, struct timeval *tp)
{
    struct ev_state_data *ev_state = (void*)cod_get_client_data(ec, 0x34567890);
    CManager cm = ev_state->cm;
    event_path_data evp = ev_state->cm->evp;
    event_item *event;
    EVstone target_stone = port_to_stone(ev_state, port);

    if (target_stone == -1) {
        printf("Port %d on stone %d invalid\n", port, ev_state->stone);
	return;
    }

    assert(CManager_locked(cm));
    ev_state->did_output++;
    if (ev_state->cur_event && data == ev_state->cur_event->decoded_event) {
	CMtrace_out(cm, EVerbose,
		    "Internal COD submit, resubmission of current input event to stone %d\n",
		    target_stone);
	if (tp) {
	    /* delayed event */
	    struct delayed_event *ev = malloc(sizeof(struct delayed_event));
	    ev->to_stone = target_stone;
	    ev->event = ev_state->cur_event;
	    ev_state->cur_event->ref_count++;
	    INT_CMadd_delayed_task(cm, tp->tv_sec, tp->tv_usec, EVdelayed_submit_func, (void*)ev);
	} else {
	    internal_path_submit(ev_state->cm, target_stone, ev_state->cur_event);
	}
    } else {
	FMFormat event_format = NULL;
	CMtrace_out(cm, EVerbose,
		    "Internal COD submit, submission of new data to stone %d\n",
		    target_stone);
	if (event_format == NULL) {
	    event_format = EVregister_format_set(cm, (FMStructDescList) type_info);
	    if (event_format == NULL) {
		printf("Bad format information on submit\n");
		return;
	    }
	}
	event = get_free_event(evp);
	event->event_encoded = 0;
	event->decoded_event = data;
	event->reference_format = event_format;
	event->format = NULL;
/*	event->free_func = cod_free_wrapper;*/
	event->free_func = NULL;
	event->free_arg = event;
	event->attrs = add_ref_attr_list(attrs);
	event->cm = cm;
	cod_encode_event(cm, event);  /* map to memory we trust */
	event->event_encoded = 1;
	event->decoded_event = NULL;  /* lose old data */
	if (tp) {
	    /* delayed event */
	    struct delayed_event {
		EVstone to_stone;
		event_item *event;
	    };
	    CMTaskHandle handle;
	    struct delayed_event *ev = malloc(sizeof(struct delayed_event));
	    ev->to_stone = target_stone;
	    ev->event = event;
	    handle = INT_CMadd_delayed_task(cm, tp->tv_sec, tp->tv_usec, EVdelayed_submit_func, (void*)ev);
	    free(handle);
	} else {
	    internal_path_submit(cm, target_stone, event);
	    return_event(cm->evp, event);
	}
    }
}

static void
internal_cod_submit_attr(cod_exec_context ec, int port, void *data, void *type_info, attr_list attrs)
{
    internal_cod_submit_general(ec, port, data, type_info, attrs, NULL);
}

static void
internal_cod_submit(cod_exec_context ec, int port, void *data, void *type_info)
{
    internal_cod_submit_general(ec, port, data, type_info, NULL, NULL);
}

#ifdef _MSC_VER
static long lrand48()
{
    return rand();
}

static double drand48() {
    return (double)(rand()) / (double)(RAND_MAX);
}
static void
sleep(int t)
{
    Sleep(t * 1000);
}
#endif
static void
add_standard_routines(stone_type stone, cod_parse_context context)
{
    static char extern_string[] = "\
		int printf(string format, ...);\n\
		void *malloc(int size);\n\
		void sleep(int seconds);\n\
		void free(void *pointer);\n\
		long lrand48();\n\
		double drand48();\n\
		int EVmax_output(cod_exec_context ec);\n\
		int EVtarget_stone_on_port(cod_exec_context ec, int port);\n\
		void EVsubmit(cod_exec_context ec, int port, void* d, cod_type_spec dt);\n\
		void EVsubmit_attr(cod_exec_context ec, int port, void* d, cod_type_spec dt, attr_list list);\n\
		void EVsubmit_delayed(cod_exec_context ec, int port, void* d, cod_type_spec dt, attr_list list, timeval *tp);\n\
        	attr_list EVget_stone_attrs(cod_exec_context ec, char *stone_name);\n \
		attr_list stone_attrs;\n";
		//time_t time(time_t *timer);\n";

    static cod_extern_entry externs[] = {
	{"printf", (void *) 0},
	{"malloc", (void*) 0},
	{"free", (void*) 0},
	{"lrand48", (void *) 0},
	{"drand48", (void *) 0},
	{"stone_attrs", (void *) 0},
	{"EVsubmit", (void *) 0},
	{"EVsubmit_attr", (void *) 0},
	{"EVsubmit_delayed", (void *) 0},
	{"sleep", (void*) 0},
	{"EVmax_output", (void*)0},
	{"EVtarget_stone_on_port", (void*)0},
        {"EVget_stone_attrs",  (void *)0},
	{(void *) 0, (void *) 0}
    };

    //{"time", (void*) 0},

    /*
     * some compilers think it isn't a static initialization to put this
     * in the structure above, so do it explicitly.
     */

    externs[0].extern_value = (void *) (intptr_t) printf;
    externs[1].extern_value = (void *) (intptr_t) malloc;
    externs[2].extern_value = (void *) (intptr_t) free;
    externs[3].extern_value = (void *) (intptr_t) lrand48;
    externs[4].extern_value = (void *) (intptr_t) drand48;
    externs[5].extern_value = (void *) (intptr_t) &stone->stone_attrs;
    externs[6].extern_value = (void *) (intptr_t) &internal_cod_submit;
    externs[7].extern_value = (void *) (intptr_t) &internal_cod_submit_attr;
    externs[8].extern_value = (void *) (intptr_t) &internal_cod_submit_general;
    externs[9].extern_value = (void *) (intptr_t) &sleep;
    externs[10].extern_value = (void *) (intptr_t) &cod_max_output;
    externs[11].extern_value = (void *) (intptr_t) &cod_target_stone_on_port;
    externs[12].extern_value = (void *) (intptr_t) &cod_ev_get_stone_attrs;

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);
}

static void
add_typed_queued_routines(cod_parse_context context, int index, const char *fmt_name)
{
    char *extern_string;
    char *data_extern_string;
    static char *extern_string_fmt =
        "void EVdiscard_%s(cod_exec_context ec, cod_closure_context type, int index);\n"
        "int EVcount_%s(cod_exec_context ec, cod_closure_context type);\n"
        "int EVpresent_%s(cod_exec_context ec, cod_closure_context queue, int index);\n"
        "void EVdiscard_and_submit_%s(cod_exec_context ec, int target, cod_closure_context queue, int index);\n"
        "void EVsubmit_%s(cod_exec_context ec, int target, cod_closure_context queue, int index);\n"
        "attr_list EVget_attrs_%s(cod_exec_context ec, cod_closure_context queue, int index);\n"
	"void write_%s(cod_exec_context ec, ffs_file fname, cod_closure_context type, int index);\n"
	"void read_%s(cod_exec_context ec, ffs_file fname, void * data, attr_list * attr_data, cod_closure_context queue);\n";
    static char *data_extern_string_fmt =
        "%s *EVdata_%s(cod_exec_context ec, cod_closure_context type, int index);\n"
        "%s *EVdata_full_%s(cod_exec_context ec, cod_closure_context type, int index);\n";
    static cod_extern_entry externs_fmt[] = {
        {"EVdiscard_%s", (void *) 0},
        {"EVcount_%s", (void *) 0},
        {"EVpresent_%s", (void *) 0},
        {"EVdiscard_and_submit_%s", (void *) 0},
        {"EVget_attrs_%s", (void *) 0},
        {"EVsubmit_%s", (void *) 0},
	{"write_%s", (void *) 0},
	{"read_%s", (void *) 0},
        {NULL, (void *) 0}
    };
    static cod_extern_entry data_externs_fmt[] = {
        {"EVdata_%s", (void *) 0},
        {"EVdata_full_%s", (void *) 0},
        {NULL, (void *) 0}
    };
    cod_extern_entry *cur;
    cod_extern_entry *externs;
    cod_extern_entry *data_externs;

    extern_string = malloc(strlen(fmt_name) * 9 + strlen(extern_string_fmt));
    assert(extern_string);
    data_extern_string = malloc(strlen(fmt_name) * 9 + strlen(data_extern_string_fmt));

    sprintf(extern_string, extern_string_fmt,
	    fmt_name, fmt_name, fmt_name, fmt_name,
	    fmt_name, fmt_name, fmt_name, fmt_name);
    sprintf(data_extern_string, data_extern_string_fmt,
	    fmt_name, fmt_name, fmt_name, fmt_name);
    externs = malloc(sizeof(externs_fmt));
    assert(externs);
    memcpy(externs, externs_fmt, sizeof(externs_fmt));
    externs[0].extern_value = (void*) cod_ev_discard_rel;
    externs[1].extern_value = (void*) cod_ev_count;
    externs[2].extern_value = (void*) cod_ev_present;
    externs[3].extern_value = (void*) cod_ev_discard_and_submit_rel;
    externs[4].extern_value = (void*) cod_ev_get_attrs;
    externs[5].extern_value = (void*) cod_ev_submit_rel;
    externs[6].extern_value = (void*) cod_ffs_write;
    externs[7].extern_value = (void*) cod_ffs_read;

    data_externs = malloc(sizeof(externs_fmt));
    assert(data_externs);
    memcpy(data_externs, data_externs_fmt, sizeof(data_externs_fmt));
    data_externs[0].extern_value = (void*) cod_ev_get_data_rel;
    data_externs[1].extern_value = (void*) cod_ev_get_data_abs;

    for (cur = externs; cur->extern_name; ++cur) {
        char *real_name = malloc(strlen(cur->extern_name) + strlen(fmt_name));
        assert(real_name);
        sprintf(real_name, cur->extern_name, fmt_name);
        cur->extern_name = real_name;
    }

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);
    for (cur = externs; cur->extern_name; ++cur) {
	/* 
	 * the index here is the index of the queue itself, 
	 * while the index in the calls above is the index of the referenced queue item 
	 */
	cod_set_closure(cur->extern_name, (void*)(intptr_t)index, context);
        free(cur->extern_name);
    }
    free(externs);
    free(extern_string);

    if (index >= 0) {
	for (cur = data_externs; cur->extern_name; ++cur) {
	    char *real_name = malloc(strlen(cur->extern_name) + strlen(fmt_name));
	    assert(real_name);
	    sprintf(real_name, cur->extern_name, fmt_name);
	    cur->extern_name = real_name;
	}
	cod_assoc_externs(context, data_externs);
	cod_parse_for_context(data_extern_string, context);
	for (cur = data_externs; cur->extern_name; ++cur) {
	    /* 
	     * the index here is the index of the queue itself, 
	     * while the index in the calls above is the index of the referenced queue item 
	     */
	    cod_set_closure(cur->extern_name, (void*)(intptr_t)index, context);
	    free(cur->extern_name);
	}
    }
    free(data_externs);
    free(data_extern_string);
}

static void
add_queued_routines(cod_parse_context context, FMFormat *formats)
{
    static char extern_string[] = "\
        int EVconforms(cod_exec_context ec, int queue, int index);\n\
        void EVdiscard(cod_exec_context ec, int queue, int index);\n\
        void EVdiscard_full(cod_exec_context ec, cod_closure_context queue, int index);\n\
        void EVdiscard_and_submit(cod_exec_context ec, int target,\
                    int queue, int index);\n\
        void EVdiscard_and_submit_full(cod_exec_context ec, int target,\
                    cod_closure_context queue, int index);\n			       \
        void *EVdata(cod_exec_context ec, int queue, int index);\n\
        void *EVdata_full(cod_exec_context ec, cod_closure_context queue, int index);\n\
        int EVcount(cod_exec_context ec, int queue);\n\
        int EVcount_full(cod_exec_context ec, cod_closure_context type);\n\
        int EVpresent(cod_exec_context ec, int queue, int index);\n\
	int EVget_port(cod_exec_context ec, int queue);\n\
    	int EVtarget_size(cod_exec_context ec, int outstone);\n\
        attr_list EVget_attrs(cod_exec_context ec, int queue, int index);\n	\
        attr_list EVget_attrs_full(cod_exec_context ec, cod_closure_context queue, int index);\n\
";

    static cod_extern_entry externs[] = {
        {"EVconforms", (void *)0},  //0 
        {"EVdiscard", (void *)0},  //1 
        {"EVdiscard_full",  (void *)0},  //2 
        {"EVdiscard_and_submit", (void *)0},  //3 
        {"EVdiscard_and_submit_full", (void *)0},  //4 
        {"EVdata", (void *)0},  //5 
        {"EVdata_full", (void *)0},  //6 
        {"EVcount", (void *)0},  //7 
        {"EVcount_full", (void *)0},  //8 
        {"EVpresent", (void *)0},  //9 
        {"EVget_port", (void *)0},  //10 
        {"EVtarget_size", (void *)0},  //11 
        {"EVget_attrs", (void *)0},  //12 
        {"EVget_attrs_full",  (void *)0},  //13 
        {(void *)0, (void *)0}
    };
    int i;
    FMFormat *cur;

    externs[0].extern_value = (void*)cod_ev_conforms;
    externs[1].extern_value = (void*)cod_ev_discard_rel;
    externs[2].extern_value = (void*)cod_ev_discard_rel;
    externs[3].extern_value = (void*)cod_ev_discard_and_submit_rel;
    externs[4].extern_value = (void*)cod_ev_discard_and_submit_rel;
    externs[5].extern_value = (void*)cod_ev_get_data_rel;
    externs[6].extern_value = (void*)cod_ev_get_data_rel;
    externs[7].extern_value = (void*)cod_ev_count;
    externs[8].extern_value = (void*)cod_ev_count;
    externs[9].extern_value = (void*)cod_ev_present;
    externs[10].extern_value = (void*)cod_ev_get_port;
    externs[11].extern_value = (void*)cod_ev_target_size;
    externs[12].extern_value = (void*)cod_ev_get_attrs;
    externs[13].extern_value = (void*)cod_ev_get_attrs;

    cod_assoc_externs(context, externs);
    cod_parse_for_context(extern_string, context);
    cod_set_closure("EVdiscard_full", (void*)(intptr_t)-1, context);
    cod_set_closure("EVdiscard_and_submit_full", (void*)(intptr_t)-1, context);
    cod_set_closure("EVget_attrs_full", (void*)(intptr_t)-1, context);
    cod_set_closure("EVdata_full", (void*)(intptr_t)-1, context);
    cod_set_closure("EVcount_full", (void*)(intptr_t)-1, context);

    for (cur = formats, i = 0; *cur; ++cur, ++i) {
        add_typed_queued_routines(context, i, name_of_FMformat(*cur));
    }
    add_typed_queued_routines(context, -2, "anonymous");
}

static void
add_queued_constants(cod_parse_context context, FMFormat *formats)
{
    FMFormat *cur_format;
    int i = 0;
    for (cur_format = formats; *cur_format; ++cur_format, ++i) {
        const char *fmt_name = name_of_FMformat(*cur_format);
        char *name = malloc(4 + strlen(fmt_name));
        sprintf(name, "%s_ID", fmt_name);
        cod_add_int_constant_to_parse_context(name, i, context);
	free(name);
    }
}

#ifdef HAVE_COD_H
extern sm_ref
cod_build_type_node(const char *name, FMFieldList field_list);
extern sm_ref
cod_build_param_node(const char *id, sm_ref typ, int param_num);
extern void
cod_add_decl_to_parse_context(const char *name, sm_ref item, cod_parse_context context);
extern void
cod_add_param(const char *id, const char *typ, int param_num,
	      cod_parse_context context);

static void
add_param(cod_parse_context parse_context, char *name, int param_num,
	  FMFormat format)
{
    FMStructDescList list = format_list_of_FMFormat(format);
    int i = 1;
    sm_ref type, param;
    while (list[i].format_name != NULL) {
	FMFieldList fl = list[i].field_list;
	/* step through input formats */
	cod_add_simple_struct_type(list[i].format_name, fl, parse_context);
	i++;
    }
    type = cod_build_type_node(list[0].format_name, list[0].field_list);
    cod_add_decl_to_parse_context(list[0].format_name, type, parse_context);

    param = cod_build_param_node(name, type, param_num);

    cod_add_decl_to_parse_context(name, param, parse_context);
}

static void
add_type(cod_parse_context parse_context, FMFormat format)
{
    FMStructDescList list = format_list_of_FMFormat(format);
    for (; list->format_name; ++list) {
	cod_add_simple_struct_type(list->format_name, list->field_list, parse_context);
    }
}

#if 0
static void
add_param_list(cod_parse_context parse_context, char *name, int param_num,
	  FMStructDescList list)
{
    char *tname = malloc(strlen(name) + strlen("_type") +1);
    sm_ref type, param;
    int i = 0;
    while (list[i].format_name != NULL) {
	sm_ref typ;
	/* step through input formats */
	typ = cod_build_type_node(list[i].format_name,
				  list[i].field_list);
	cod_add_decl_to_parse_context(list[i].format_name, typ,
				      parse_context);
	i++;
    }
    sprintf(tname, "%s_type", name);
    type = cod_build_type_node(tname, list[i-1].field_list);
    cod_add_decl_to_parse_context(tname, type, parse_context);

    param = cod_build_param_node(name, type, 0);

    cod_add_decl_to_parse_context(name, param, parse_context);
}
#endif

static int
dll_prefix_present(char *filter)
{

    if (filter[0] == 'd' && filter[1] == 'l' && filter[2] == 'l' && filter[3] == ':') {
	return 1;
    }
    return 0;
}

static char *
extract_dll_path(char *filter)
{
    char *copy = strdup(filter);
    char *temp;
    char *path;


    temp = strtok(copy, ":");
    if (strcmp(temp, "dll")) {
	free(copy);
	return NULL;
    }
    temp = strtok(NULL, ":");

    if (temp == NULL) {
	free(copy);
	return NULL;
    }

    path = strdup(temp);
    free(copy);

    return path;
}

static char *
extract_symbol_name(char *filter)
{

    char *copy = strdup(filter);
    char *temp;
    char *symbol;

    temp = strtok(copy, ":");
    if (strcmp(temp, "dll")) {
	free(copy);
	return NULL;
    }
    temp = strtok(NULL, ":");
    temp = strtok(NULL, ":");

    if (temp == NULL) {
	free(copy);
	return NULL;
    }

    symbol = strdup(temp);
    free(copy);

    return symbol;
}

static void*
load_dll_symbol(CManager cm, char *path, char *symbol_name)
{
    lt_dlhandle handle;

    handle = CMdlopen(cm->CMTrace_file, path, 0);
    if (!handle) {
    	fprintf(stderr, "failed opening %s\n", path);
	    return NULL;
    }
    return lt_dlsym(handle, symbol_name);
}

extern void
add_metrics_routines(stone_type stone, cod_parse_context context);

static response_instance
generate_filter_code(CManager cm, struct response_spec *mrd, stone_type stone,
		     FMFormat format)
{
    response_instance instance = malloc(sizeof(*instance));

    cod_code code;
    cod_parse_context parse_context = new_cod_parse_context();
    /*    sm_ref conn_info_data_type, conn_info_param;*/

    memset(instance, 0, sizeof(*instance));
    add_standard_routines(stone, parse_context);
    add_metrics_routines(stone, parse_context);
    if (cm->evp->extern_structs) {
	int count = -1;
	while(cm->evp->extern_structs[++count] != NULL) {
	    cod_add_struct_type(cm->evp->extern_structs[count], parse_context);
	}
    }
	
    if (cm->evp->externs) {
	int count = -1;
	while (cm->evp->externs[++count].extern_decl != NULL) {
	    cod_assoc_externs(parse_context, cm->evp->externs[count].externs);
	    cod_parse_for_context(cm->evp->externs[count].extern_decl, 
				  parse_context);
	}
    }

    switch (mrd->response_type) {
    case Response_Filter:
    case Response_Router:
    case Response_Transform:
	cod_add_param("ec", "cod_exec_context", 0, parse_context);
	if (format) {
	    add_param(parse_context, "input", 1, format);
	} else {
	    cod_add_param("input", "int", 1, parse_context);
	}
	if (mrd->response_type == Response_Transform) {
	    add_param(parse_context, "output", 2,
		      mrd->u.transform.reference_output_format);
	    cod_add_param("event_attrs", "attr_list", 3, parse_context);
	    cod_add_param("output_attrs", "attr_list", 4, parse_context);
	} else {
	    cod_add_param("event_attrs", "attr_list", 2, parse_context);
	}
	break;
    case Response_Multityped:
        /* this should call generate_multityped_code() */
        assert(FALSE);
	break;
    }

/*    conn_info_data_type = cod_build_type_node("output_conn_info_type",
					      output_conn_field_list);
    cod_add_decl_to_parse_context("output_conn_info_type",
				  conn_info_data_type, parse_context);
    conn_info_param = cod_build_param_node("output_conn_info",
					   conn_info_data_type, 3);
    cod_add_decl_to_parse_context("output_conn_info", conn_info_param,
				  parse_context);
*/
    switch(mrd->response_type) {
    case Response_Filter:
    case Response_Router:
	if (dll_prefix_present(mrd->u.filter.function)) {
	    /* it is a dll */
	    char *path = NULL;
	    char *symbol_name = NULL;

	    path = extract_dll_path(mrd->u.filter.function);
	    symbol_name = extract_symbol_name(mrd->u.filter.function);
	    if (!path || !symbol_name) {
		fprintf(stderr, "could not parse string \"%s\" for dll path and symbol information\n", mrd->u.filter.function);
		free(instance);
		return NULL;
	    }
	    instance->u.filter.func_ptr = (int(*)(void*,attr_list)) load_dll_symbol(cm, path, symbol_name);
	    if (instance->u.filter.func_ptr == NULL) {
		fprintf(stderr, "Failed to load symbol \"%s\" from file \"%s\"\n",
			symbol_name, path);
		free(instance);
		free(path);
		free(symbol_name);
		return NULL;
	    }
	    free(symbol_name);
	    free(path);
	    instance->u.filter.code = NULL;
	} else {
	    code = cod_code_gen(mrd->u.filter.function, parse_context);
	    instance->response_type = mrd->response_type;
	    instance->u.filter.code = code;
	    if (code)
		instance->u.filter.ec = cod_create_exec_context(code);

	    instance->u.filter.func_ptr = NULL;
	}
	break;
    case Response_Transform:
	if (dll_prefix_present(mrd->u.transform.function)) {
	    /* it is a dll */
	    char *path = NULL;
	    char *symbol_name = NULL;

	    path = extract_dll_path(mrd->u.transform.function);
	    symbol_name = extract_symbol_name(mrd->u.transform.function);
	    if (!path || !symbol_name) {
		fprintf(stderr, "could not parse string \"%s\" for dll path and symbol information\n", mrd->u.transform.function);
		free(instance);
		return NULL;
	    }
	    instance->u.transform.func_ptr =
		(int(*)(void*,void*,attr_list,attr_list)) load_dll_symbol(cm, path, symbol_name);
	    if (instance->u.transform.func_ptr == NULL) {
		fprintf(stderr, "Failed to load symbol \"%s\" from file \"%s\"\n",
			symbol_name, path);
		free(instance);
		free(path);
		free(symbol_name);
		return NULL;
	    }
	    instance->u.transform.code = NULL;
	    free(path);
	    free(symbol_name);
	} else {
	    code = cod_code_gen(mrd->u.transform.function, parse_context);
	    instance->response_type = Response_Transform;
	    instance->u.transform.code = code;
	    if (code)
		instance->u.transform.ec = cod_create_exec_context(code);
	}
	instance->u.transform.out_size =
	    mrd->u.transform.output_base_struct_size;
	instance->u.transform.out_format =
	    mrd->u.transform.reference_output_format;
	break;
    case Response_Multityped:
	break;
    }
    cod_free_parse_context(parse_context);

    return instance;
}

#ifdef NOT_DEF
static int
verify_multityped_code(CManager cm, struct response_spec *mrd, stone_type stone,
			 FMFormat *formats)
{
    FMFormat *cur_format;
    int ret;

/*    cod_code code;*/
    cod_parse_context parse_context = new_cod_parse_context();
    /*    sm_ref conn_info_data_type, conn_info_param;*/

    for (cur_format = formats; *cur_format; ++cur_format) {
        add_type(parse_context, *cur_format);
    }

    add_standard_routines(stone, parse_context);
    add_metrics_routines(stone, parse_context);
    add_queued_routines(parse_context, formats);
    add_queued_constants(parse_context, formats);
    if (cm->evp->extern_structs) {
	int count = -1;
	while(cm->evp->extern_structs[++count] != NULL) {
	    cod_add_struct_type(cm->evp->extern_structs[count], parse_context);
	}
    }
	
    if (cm->evp->externs) {
	int count = -1;
	while (cm->evp->externs[++count].extern_decl != NULL) {
	    cod_assoc_externs(parse_context, cm->evp->externs[count].externs);
	    cod_parse_for_context(cm->evp->externs[count].extern_decl, 
				  parse_context);
	}
    }



    assert(mrd->response_type == Response_Multityped);
    cod_add_param("ec", "cod_exec_context", 0, parse_context);
    ret = cod_code_verify(mrd->u.multityped.function, parse_context);
    return ret;
}
#endif

static response_instance
generate_multityped_code(CManager cm, struct response_spec *mrd, stone_type stone,
			 FMFormat *formats)
{
    response_instance instance = malloc(sizeof(*instance));
    FMFormat *cur_format;
    int format_count = 0;

    cod_code code;
    cod_parse_context parse_context = new_cod_parse_context();
    /*    sm_ref conn_info_data_type, conn_info_param;*/

    memset(instance, 0, sizeof(*instance));

    for (cur_format = formats; *cur_format; ++cur_format) {
	add_type(parse_context, *cur_format);
	format_count++;
    }

    add_standard_routines(stone, parse_context);
    add_metrics_routines(stone, parse_context);
    add_queued_routines(parse_context, formats);
    add_queued_constants(parse_context, formats);
    if (cm->evp->extern_structs) {
	int count = -1;
	while(cm->evp->extern_structs[++count] != NULL) {
	    cod_add_struct_type(cm->evp->extern_structs[count], parse_context);
	}
    }
	
    if (cm->evp->externs) {
	int count = -1;
	while (cm->evp->externs[++count].extern_decl != NULL) {
	    cod_assoc_externs(parse_context, cm->evp->externs[count].externs);
	    cod_parse_for_context(cm->evp->externs[count].extern_decl, 
				  parse_context);
	}
    }



    assert(mrd->response_type == Response_Multityped);
    cod_add_param("ec", "cod_exec_context", 0, parse_context);
/*    if (format) {
	add_param(parse_context, "input", 1, format);
    } else {
	cod_add_param("input", "int", 1, parse_context);
	}*/
    cod_set_return_type("void", parse_context);
    code = cod_code_gen(mrd->u.multityped.function, parse_context);
    instance->response_type = mrd->response_type;
    instance->u.queued.ref_count = format_count;
    instance->u.queued.formats = formats;
    instance->u.queued.code = code;
    if (code)
	instance->u.queued.ec = cod_create_exec_context(code);

    cod_free_parse_context(parse_context);

    if (!instance->u.queued.ec) {
        free(instance);
        return NULL;
    }

    return instance;
}

#else
static response_instance
generate_filter_code(CManager cm, struct response_spec *mrd, stone_type stone,
		     FMFormat format){return NULL;}
static response_instance
generate_multityped_code(CManager cm, struct response_spec *mrd, stone_type stone,
			 FMFormat *formats){return NULL;}
#endif
