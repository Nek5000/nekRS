#include "config.h"
#if !NO_DYNAMIC_LINKING
#include "dlloader.h"
#endif
#include <stdio.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>
#include <ctype.h>
#include <math.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <limits.h>
#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#define __ANSI_CPP__
#else
#include <netinet/in.h>
#include <arpa/inet.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <inttypes.h>
#include <ffs.h>
#include <atl.h>
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"
#include "cm_transport.h"

#if NO_DYNAMIC_LINKING
struct select_data;
extern void libcmselect_LTX_add_select(CMtrans_services svc, struct select_data *sdp, int fd, 
				       select_list_func func, void *arg1, void *arg2);
extern void libcmselect_LTX_write_select(CMtrans_services svc, struct select_data *sdp, int fd, 
			     select_list_func func, void *arg1, void *arg2);
extern periodic_task_handle libcmselect_LTX_add_periodic(CMtrans_services svc, struct select_data *sdp, int interval_sec, int interval_usec, select_list_func func, void *arg1, void *arg2);
extern periodic_task_handle libcmselect_LTX_add_delayed_task(CMtrans_services svc, 
struct select_data *sdp, int delay_sec, int delay_usec, select_list_func func,void *arg1, 
void *arg2);
extern void libcmselect_LTX_remove_periodic(CMtrans_services svc, struct select_data *sdp, periodic_task_handle handle);
extern void libcmselect_LTX_remove_select(CMtrans_services svc, struct select_data *sdp, int fd);
extern void libcmselect_LTX_wake_function(CMtrans_services svc, struct select_data *sdp);
extern void libcmselect_LTX_blocking_function(CMtrans_services svc, void *client_data);
extern void libcmselect_LTX_polling_function(CMtrans_services svc,void *client_data);
extern void libcmselect_LTX_select_initialize(CMtrans_services svc,CManager cm,void *client_data);
extern void libcmselect_LTX_select_shutdown(CMtrans_services svc,CManager cm,void *client_data);
extern void libcmselect_LTX_select_free(CMtrans_services svc,CManager cm,void *client_data);
extern void libcmselect_LTX_select_stop(CMtrans_services svc,void *client_data);
#endif


static void CMinitialize (CManager cm);

static atom_t CM_TRANSPORT = -1;
static atom_t CM_NETWORK_POSTFIX = -1;
static atom_t CM_CONN_BLOCKING = -1;
atom_t CM_REBWM_RLEN = -1;
atom_t CM_REBWM_REPT = -1;
atom_t CM_BW_MEASURE_INTERVAL = -1;
atom_t CM_BW_MEASURE_TASK = -1;
atom_t CM_BW_MEASURED_VALUE = -1;
atom_t CM_BW_MEASURED_COF = -1;
atom_t CM_BW_MEASURE_SIZE = -1;
atom_t CM_BW_MEASURE_SIZEINC = -1;
static atom_t CM_EVENT_SIZE = -1;
static atom_t CM_INCOMING_CONNECTION = -1;
static atom_t CM_TRANSPORT_RELIABLE = -1;
static atom_t CM_IP_INTERFACE = -1;

void wait_for_pending_write(CMConnection conn);
static void cm_wake_any_pending_write(CMConnection conn);
static void transport_wake_any_pending_write(CMConnection conn);
static void cm_set_pending_write(CMConnection conn);
static int drop_CM_lock(CManager cm, const char *file, int line);
static int acquire_CM_lock(CManager cm, const char *file, int line);
static int return_CM_lock_status(CManager cm, const char *file, int line);
static void add_buffer_to_pending_queue(CManager cm, CMConnection conn, CMbuffer buf, long length);
static void cond_wait_CM_lock(CManager cm, void *cond, char *file, int line);

struct CMtrans_services_s CMstatic_trans_svcs = {INT_CMmalloc, INT_CMrealloc, INT_CMfree, 
						 INT_CM_fd_add_select, 
						 CM_fd_write_select, 
						 CM_fd_remove_select, 
						 CMtransport_trace,
						 CMtransport_verbose,
						 CMConnection_create,
						 INT_CMadd_shutdown_task,
						 INT_CMadd_periodic_task,
						 INT_CMremove_periodic,
						 INT_CMadd_poll,
						 cm_get_data_buf,
						 cm_return_data_buf,
						 internal_connection_close,
						 cm_create_transport_buffer,
						 cm_create_transport_and_link_buffer,
						 INT_CMget_transport_data,
						 cm_set_pending_write,
						 transport_wake_any_pending_write,
						 drop_CM_lock,
						 acquire_CM_lock,
						 return_CM_lock_status,
						 cond_wait_CM_lock,
						 add_buffer_to_pending_queue,
						 INT_CMConnection_dereference,
						 INT_CMConnection_add_reference,
						 INT_CMConnection_failed,
						 CMwake_server_thread,
						 INT_CMCondition_signal
};
static void INT_CMControlList_close(CMControlList cl, CManager cm);
static int CMcontrol_list_poll(CMControlList cl);
int CMdo_non_CM_handler(CMConnection conn, int header,
			      char *buffer, size_t length);
void CMdo_performance_response(CMConnection conn, size_t length,
					    int func, int byte_swap,
					    char *buffer);

void CMhttp_handler(CMConnection conn, char* buffer, size_t length);
static void CM_init_select(CMControlList cl, CManager cm);

static void cond_wait_CM_lock(CManager cm, void *vcond, char *file, int line)
{
    thr_condition_t *cond = vcond;
    CMtrace_out(cm, CMLowLevelVerbose, "CManager Condition wait at \"%s\" line %d\n",
		file, line);
    cm->locked--;
    thr_condition_wait(*cond, cm->exchange_lock);
    CMtrace_out(cm, CMLowLevelVerbose, "CManager Condition wake at \"%s\" line %d\n",
		file, line);
    cm->locked++;
}

static int drop_CM_lock(CManager cm, const char *file, int line)
{
    int ret = cm->locked;
    IntCManager_unlock(cm, file, line);
    return ret;
}

static int acquire_CM_lock(CManager cm, const char *file, int line)
{
    IntCManager_lock(cm, file, line);
    return cm->locked;
}

static int return_CM_lock_status(CManager cm, const char *file, int line)
{
    (void) file;
    (void) line;
    return cm->locked;
}

static void
CMpoll_forever(CManager cm)
{
    CMControlList cl = cm->control_list;
    int should_exit = 0;
    CManager_lock(cm);
    if (!cm->control_list->select_initialized) {
	CM_init_select(cm->control_list, cm);
    }
    if (cl->has_thread > 0 && cl->server_thread == thr_thread_self()) {
	/* 
	 * if we're actually the server thread here, do a thread exit when
	 * we're done
	 */
	should_exit++;
    }
    while(!cl->closed) {
	CMtrace_out(cm, CMLowLevelVerbose, "CM Poll Forever - thread %zx doing wait\n", (size_t)thr_thread_self());
	if (CMcontrol_list_wait(cl) == -1) {
	    CMtrace_out(cm, CMLowLevelVerbose, "CM Poll Forever - doing close and exit\n");
	    /* 
	     * error.  others will free the CM too, add to the ref count 
	     * here so we can close.
	     */
	    cm->reference_count++;
	    CManager_unlock(cm);
	    CManager_close(cm);
	    exit(1);
	}
    }
    CMtrace_out(cm, CMLowLevelVerbose, "CM Poll Forever - doing close\n");
    CManager_unlock(cm);
    CManager_close(cm);
    if (should_exit != 0) thr_thread_exit(NULL);
}

static void CManager_free(CManager cm);

static void
server_thread_func(CManager cm)
{
    CMpoll_forever(cm);
    CManager_free(cm);
}

extern void
INT_CMrun_network(CManager cm)
{
    if (!cm->control_list->select_initialized) {
	CM_init_select(cm->control_list, cm);
    }
    if ((cm->control_list->server_thread != 0) &&
	(cm->control_list->server_thread != thr_thread_self())) {
	/* What?  We're polling, but we're not the server thread? */
	fprintf(stderr, "Warning:  CMrun_network() called when another thread may already be handling the network\n");
	fprintf(stderr, "          This situation may result in unexpected I/O blocking.\n");
	fprintf(stderr, "          Server thread set to %zx.\n", (size_t) thr_thread_self());
    }
    cm->control_list->server_thread = thr_thread_self();
    cm->control_list->has_thread = 1;
    CManager_unlock(cm);
    CMpoll_forever(cm);
}

static int
CM_test_thread_func()
{
    return 1;
}

static thr_thread_t 
thr_fork(void*(*func)(void*), void *arg)
{
    thr_thread_t new_thread = 0;
    int err = thr_thread_create(&new_thread, NULL, (void*(*)(void*))func, arg);
    if (err != 0) {
	return (thr_thread_t) (intptr_t)NULL;
    } else {
	return (thr_thread_t) new_thread;
    }
}

int
INT_CMfork_comm_thread(CManager cm)
{
    /* if we're on a kernel-level-threads package, for the thread and 
       return 1, else return 0; */
    if (!cm->control_list->select_initialized) {
	CM_init_select(cm->control_list, cm);
    }
    if (cm->control_list->has_thread == 0) {
	if (cm->control_list->network_blocking_function.func) {
	    thr_thread_t server_thread = 
		thr_fork((void*(*)(void*))server_thread_func, 
			 (void*)cm);
	    CMtrace_out(cm, CMLowLevelVerbose,
			"CM - Forked comm thread %p\n", (void*)(intptr_t)server_thread);
	    if (server_thread ==  (thr_thread_t)(intptr_t) NULL) {
		return 0;
	    }
	    cm->control_list->server_thread = server_thread;
	    cm->control_list->has_thread = 1;
	    cm->reference_count++;
	    CMtrace_out(cm, CMFreeVerbose, "Forked - CManager %p ref count now %d\n", 
			cm, cm->reference_count);
	    cm->control_list->cl_reference_count++;
	    cm->control_list->free_reference_count++;
	} else {
	    /*
	     *  Can't start a server thread yet, but lets see 
	     *  if we can fork anything successfully.
	     */
	    thr_thread_t test_thread = 
		thr_fork((void*(*)(void*))CM_test_thread_func, 
			 (void*)cm);
	    if (test_thread ==  (thr_thread_t)(intptr_t) NULL) {
		/* No.  Say we can't. */
		CMtrace_out(cm, CMLowLevelVerbose,
			    "CM - Test fork failed, no comm thread\n");
		return 0;
	    }
	    /* OK, we'll fork it later. */
	    CMtrace_out(cm, CMLowLevelVerbose,
			"CM - Will fork comm thread later\n");
	    cm->control_list->has_thread = -1; /* should fork one */
	}
    }
    return 1;
}

extern
void
CMControlList_set_blocking_func(CMControlList cl, CManager cm, 
				CMPollFunc bfunc, CMPollFunc pfunc,
				void *client_data)
{
}

extern void
INT_CMpoll_network(CManager cm)
{
    CMControlList cl = cm->control_list;
    CMtrace_out(cm, CMLowLevelVerbose, "CM Poll Network\n");
    cl->network_polling_function.func((void*)&CMstatic_trans_svcs,
				      cl->network_polling_function.client_data);
    CMcontrol_list_poll(cl);
}

static void
add_contact_list(CManager cm, attr_list attrs)
{
    int list_size = 0;
    if (cm->contact_lists == NULL) {
	cm->contact_lists = INT_CMmalloc(sizeof(attr_list) *2);
	list_size = 0;
    } else {
	while(cm->contact_lists[list_size] != NULL) list_size++;
	cm->contact_lists = INT_CMrealloc(cm->contact_lists, 
				      sizeof(attr_list) * (list_size + 2));
    }
    cm->contact_lists[list_size] = attrs;
    cm->contact_lists[list_size+1] = NULL;
}

void
INT_CM_insert_contact_info(CManager cm, attr_list attrs)
{
    attr_merge_lists(cm->contact_lists[0], attrs);
}

attr_list
INT_CMget_contact_list(CManager cm)
{
    if (cm->contact_lists == NULL) return NULL;
    CMadd_ref_attr_list(cm, cm->contact_lists[0]);
    return (cm->contact_lists[0]);
}

extern attr_list
INT_CMderef_and_copy_list(CManager cm, attr_list attrs)
{
  // done inside the CM lock, so a safe way to convert a shared list to an owned list
  attr_list ret = attr_copy_list(attrs);
  free_attr_list(attrs);
  return ret;
}

extern attr_list
INT_CMget_specific_contact_list(CManager cm, attr_list attrs)
{
    char *chosen_transport = NULL, *chosen_net = NULL, *chosen_interface = NULL;
    char *freeable_transport = NULL;
    int i = 0;

    if (attrs != NULL) {
	get_string_attr(attrs, CM_TRANSPORT, &chosen_transport);
    }
    if (chosen_transport && (strchr(chosen_transport, ':') != NULL)) {
	freeable_transport = strdup(chosen_transport);
	*(strchr(freeable_transport, ':')) = 0;
	chosen_transport = freeable_transport;
    }
    if (attrs != NULL) {
	get_string_attr(attrs, CM_NETWORK_POSTFIX, &chosen_net);
    }
    if (attrs != NULL) {
	get_string_attr(attrs, CM_IP_INTERFACE, &chosen_interface);
    }
    if ((chosen_transport == NULL) && (chosen_net == NULL) && (chosen_interface == NULL)) {
	CMadd_ref_attr_list(cm, cm->contact_lists[0]);
	return cm->contact_lists[0];
    }
    /* specific transport or interface chosen */
    i = 0;
    while (cm->contact_lists && (cm->contact_lists[i] != NULL)) {
	char *this_transport = NULL, *this_postfix = NULL, *this_interface = NULL;

	get_string_attr(cm->contact_lists[i], CM_TRANSPORT, &this_transport);
	get_string_attr(cm->contact_lists[i], CM_NETWORK_POSTFIX, &this_postfix);
	get_string_attr(cm->contact_lists[i], CM_IP_INTERFACE, &this_interface);
	if (this_transport == NULL) {
	    this_transport = "sockets";
	}
	if (chosen_transport == NULL) {
	    chosen_transport = "sockets";
	}
	if (strcmp(this_transport, chosen_transport) == 0) {
	    if ((chosen_net != NULL) || (this_postfix != NULL)) {
		/* one is not null */
		if (chosen_net && this_postfix) {
		    if (strcmp(chosen_net, this_postfix) != 0) {
			i++;
			continue;
		    }
		} else {
		    i++;
		    continue;
		}
	    }
	    if ((chosen_interface != NULL) || (this_interface != NULL)) {
		/* one is not null */
		if (chosen_interface && this_interface) {
		    if (strcmp(chosen_interface, this_interface) != 0) {
			i++;
			continue;
		    }
		} else {
		    i++;
		    continue;
		}
	    }
	    CMadd_ref_attr_list(cm, cm->contact_lists[i]);
	    if (freeable_transport) free(freeable_transport);
	    return cm->contact_lists[i];
	}
	i++;
    }
    /* chosen transport not listened? */
    CMinternal_listen(cm, attrs, /* try others*/ 0);
    /* try again */
    i = 0;
    while (cm->contact_lists && (cm->contact_lists[i] != NULL)) {
	char *this_transport = NULL, *this_postfix = NULL, *this_interface = NULL;

	get_string_attr(cm->contact_lists[i], CM_TRANSPORT, &this_transport);
	get_string_attr(cm->contact_lists[i], CM_NETWORK_POSTFIX, 
			&this_postfix);
	get_string_attr(cm->contact_lists[i], CM_IP_INTERFACE, &this_interface);
	if (this_transport == NULL) {
	    this_transport = "sockets";
	}
	if (chosen_transport == NULL) {
	    chosen_transport = "sockets";
	}
	if (strcmp(this_transport, chosen_transport) == 0) {
	    if ((chosen_net != NULL) || (this_postfix != NULL)) {
		/* one is not null */
		if (chosen_net && this_postfix) {
		    if (strcmp(chosen_net, this_postfix) != 0) {
			i++;
			continue;
		    }
		} else {
		    i++;
		    continue;
		}
	    }
	    if ((chosen_interface != NULL) || (this_interface != NULL)) {
		/* one is not null */
		if (chosen_interface && this_interface) {
		    if (strcmp(chosen_interface, this_interface) != 0) {
			i++;
			continue;
		    }
		} else {
		    i++;
		    continue;
		}
	    }
	    CMadd_ref_attr_list(cm, cm->contact_lists[i]);
	    if (freeable_transport) free(freeable_transport);
	    return cm->contact_lists[i];
	}
	i++;
    }
    /* maybe it failed to load */
    if (freeable_transport) free(freeable_transport);
    return NULL;
}

int
INT_CMlisten(CManager cm)
{
  return INT_CMlisten_specific (cm, NULL);
}

static attr_list
split_transport_attributes(attr_list list)
{
    char *chosen_transport = NULL;
    if (list) {
	get_string_attr(list, CM_TRANSPORT, &chosen_transport);
    }
    if (chosen_transport && (strchr(chosen_transport, ':') != NULL)) {
	attr_list new_list = attr_copy_list(list);
	atom_t atom;
	char *old_transport, *params;
	char *next_param;
	get_string_attr(new_list, CM_TRANSPORT, &old_transport);
	params = strchr(old_transport, ':');
	*(params++) = 0;
	set_string_attr(new_list, CM_TRANSPORT, strdup(old_transport));
	while (params != NULL) {
	    char *equal, *end;
	    next_param = strchr(params, ',');
	    if (next_param) *(next_param++) = 0;  /* kill comma */
	    if ((equal = strchr(params, '=')) != NULL) {
		/* there's an equal sign */
		*(equal++) = 0;
		/* we'll deal with this later */
	    }
	    while (isspace(*params)) params++;  /* skip white */
	    end = params + strlen(params) - 1;
	    while(end > params && isspace(*end)) end--;
	    // Write new null terminator
	    *(end+1) = 0;
	    atom = attr_atom_from_string(params);
	    if (equal == NULL) {
		set_int_attr(new_list, atom, 1);
	    } else {
		char *tail;
		long value;
		while (isspace(*equal)) equal++;  /* skip white */
		end = equal + strlen(equal) - 1;
		while(end > equal && isspace(*end)) end--;
		// Write new null terminator
		*(end+1) = 0;
		value = strtol(equal, &tail, 10);
		if (strlen(tail) == 0) {
		    /* valid integer! */
		    if ((value < INT_MAX) && (value > INT_MIN)) {
			set_int_attr(new_list, atom, (int)value);
		    } else {
			set_long_attr(new_list, atom, value);
		    }
		} else {
		    /* string... */
		    set_string_attr(new_list, atom, strdup(equal));
		}
	    }
	    params = next_param;
	}
	free(old_transport);  /* not free'd by replace */
	free_attr_list(list);
	list = new_list;
    }
    return list;
}
extern int
CMinternal_listen(CManager cm, attr_list listen_info, int try_others)
{
    int success = 0;
    transport_entry *trans_list;
    char *chosen_transport = NULL;
    char *iface = NULL;

    if (listen_info) {
        listen_info = split_transport_attributes(attr_copy_list(listen_info));
	get_string_attr(listen_info, CM_TRANSPORT, &chosen_transport);
	get_string_attr(listen_info, CM_IP_INTERFACE, &iface);
    }
    if (chosen_transport != NULL) {
        CMtrace_out(cm, CMConnectionVerbose,
		    "CM - Listening only on transport \"%s\"\n",
		    chosen_transport);
	if (load_transport(cm, chosen_transport, 1) == 0) {
	    CMtrace_out(cm, CMConnectionVerbose,
			"Failed to load transport \"%s\".  Revert to default.\n",
					chosen_transport);
	    CMtrace_out(cm, CMTransportVerbose,
			"Failed to load transport \"%s\".  Revert to default.\n",
			chosen_transport);
	    if (!try_others) {
		if (listen_info) free_attr_list(listen_info);
		return success;
	    }
	    chosen_transport = NULL;
	}
    }
    trans_list = cm->transports;
    while ((trans_list != NULL) && (*trans_list != NULL)) {
	attr_list attrs;
	if ((chosen_transport == NULL) || 
	    (strcmp((*trans_list)->trans_name, chosen_transport) == 0)) {
	    attrs = (*trans_list)->listen(cm, &CMstatic_trans_svcs,
					  *trans_list,
					  listen_info);
	    if (iface) {
		add_string_attr(attrs, CM_IP_INTERFACE, strdup(iface));
	    }
	    add_contact_list(cm, attrs);
	    if (CMtrace_on(cm, CMConnectionVerbose)) {
		fprintf(cm->CMTrace_file, "Adding contact list -> ");
		fdump_attr_list(cm->CMTrace_file, attrs);
	    }
	    if (attrs != NULL) {
		success++;
	    }
	}
	trans_list++;
    }
    if (listen_info) free_attr_list(listen_info);
    return success;
}

int
INT_CMlisten_specific(CManager cm, attr_list listen_info)
{
    int success = 0;
    if (!cm->initialized) CMinitialize(cm);
    success = CMinternal_listen(cm, listen_info, /* try others*/ 1);
    return (success != 0);
}

#ifdef CM_DEFAULT_TRANSPORT
static char *CMglobal_default_transport = CM_DEFAULT_TRANSPORT;
#else 
static char *CMglobal_default_transport = NULL;
#endif

static char *CMglobal_alternate_transports[] = {NULL};

static void 
CMinitialize(CManager cm)
{
    char **transport_names = CMglobal_alternate_transports;
    char *def = getenv("CMDefaultTransport");
    if (def != NULL) CMglobal_default_transport = def;
    if (CMglobal_default_transport) {
	if (load_transport(cm, CMglobal_default_transport, 0) == 0) {
	    fprintf(stderr, "Failed to initialize default transport.  Exiting.\n");
	    exit(1);
	}
    }
    while ((transport_names != NULL) && (transport_names[0] != NULL)) {
	load_transport(cm, transport_names[0], 1);
	transport_names++;
    }
    cm->initialized++;
}

static int
CMcontrol_list_poll(CMControlList cl)
{
    func_entry *poll_list = cl->polling_function_list;
    while ((poll_list != NULL) && (poll_list->func != NULL)){
	int consistency_number = cl->cl_consistency_number;

	CManager_unlock(poll_list->cm);
	poll_list->func(poll_list->cm, poll_list->client_data);
	CManager_lock(poll_list->cm);
	/* do function */
	if (consistency_number != cl->cl_consistency_number) {
	    return 1;
	}
	poll_list++;
    }
    return 1;
}

static 
void
CMControlList_add_poll(CMControlList cl, CManager cm, CMPollFunc func,
		       void *client_data)
{
    func_entry *poll_list;
    int count = 0;
    poll_list = cl->polling_function_list;
    while ((poll_list != NULL) && (poll_list[count].func != NULL)) {
	count++;
    }
    /*
     *  We're going to navigate the poll list without locks.  This is
     *  somewhat dangerous, but safe enough if we only have a couple of
     *  functions so we never realloc.  If we ever hit the realloc() below,
     *  there's a chance of some bad data references.  At this point, using
     *  that many poll functions seems unlikely, but later, maybe not...
     */
    if (poll_list != NULL) {
	if (cl->pflist_size < count - 2) {
	    cl->pflist_size *= 2;
	    poll_list = INT_CMrealloc(poll_list, sizeof(func_entry) * (cl->pflist_size));
	}
    } else {
	poll_list = INT_CMmalloc(sizeof(func_entry)*10);
	cl->pflist_size = 10;
    }
    poll_list[count].cm = cm;
    poll_list[count].func = func;
    poll_list[count].client_data = client_data;
    poll_list[count+1].func = NULL;
    cl->polling_function_list = poll_list;
}
    
extern
void
INT_CMadd_poll(CManager cm, CMPollFunc func, void *client_data)
{
    CMControlList_add_poll(cm->control_list, cm, func, client_data);
}

extern
int
CMcontrol_list_wait(CMControlList cl)
{
    /* associated CM should be locked */
    if ((cl->server_thread != 0) &&
	(cl->server_thread != thr_thread_self())) {
	/* What?  We're polling, but we're not the server thread? */
	fprintf(stderr, "Warning:  Multiple threads calling CMnetwork_wait\n");
	fprintf(stderr, "          This situation may result in unexpected I/O blocking.\n");
	fprintf(stderr, "          Server thread set to %zx.\n", (size_t) thr_thread_self());
    }
    cl->server_thread = thr_thread_self();
    if (cl->network_blocking_function.func != NULL) {
	cl->network_blocking_function.func((void*)&CMstatic_trans_svcs,
					   cl->network_blocking_function.client_data);
    }
    CMcontrol_list_poll(cl);
    return 1;
}

static CMControlList CMControlList_create();

static thr_mutex_t atl_mutex;
static int atl_mutex_initialized = 0;
static void process_pending_queue(CManager cm, void *junk);

extern
CManager
INT_CManager_create()
{
    return INT_CManager_create_control(NULL);
}

static void
atl_mutex_lock(void* lock)
{
    thr_mutex_lock(*((thr_mutex_t*)lock));
}
static void
atl_mutex_unlock(void* lock)
{
    thr_mutex_unlock(*((thr_mutex_t*)lock));
}
extern
CManager
INT_CManager_create_control(char *control_module)
{
    CManager cm = (CManager) INT_CMmalloc(sizeof(CManager_s));
    int atom_init = 0;

    if (!atl_mutex_initialized) {
	atl_mutex_initialized++;
	thr_mutex_init(atl_mutex);
	atl_install_mutex_funcs((atl_lock_func)atl_mutex_lock, (atl_lock_func)atl_mutex_unlock, 
				&atl_mutex);
    }
    if (cm == NULL)
	return NULL;
    memset(cm, 0, sizeof(CManager_s));

    if (atom_init == 0) {
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
	CM_CONN_BLOCKING = attr_atom_from_string("CM_CONN_BLOCKING");
	CM_REBWM_RLEN = attr_atom_from_string("CM_REG_BW_RUN_LEN");
	CM_REBWM_REPT = attr_atom_from_string("CM_REG_BW_REPEAT_CNT");
	CM_BW_MEASURE_INTERVAL = attr_atom_from_string("CM_BW_MEASURE_INTERVAL");
	CM_BW_MEASURE_TASK = attr_atom_from_string("CM_BW_MEASURE_TASK");
	CM_BW_MEASURED_VALUE = attr_atom_from_string("CM_BW_MEASURED_VALUE");
	CM_BW_MEASURED_COF = attr_atom_from_string("CM_BW_MEASURED_COF");
	CM_BW_MEASURE_SIZE = attr_atom_from_string("CM_BW_MEASURE_SIZE");
	CM_BW_MEASURE_SIZEINC = attr_atom_from_string("CM_BW_MEASURE_SIZEINC");
	CM_EVENT_SIZE = attr_atom_from_string("CM_EVENT_SIZE");
	CM_INCOMING_CONNECTION = attr_atom_from_string("CM_INCOMING_CONNECTION");
	CM_TRANSPORT_RELIABLE = attr_atom_from_string("CM_TRANSPORT_RELIABLE");
	CM_IP_INTERFACE = attr_atom_from_string("IP_INTERFACE");
    }

    /* initialize data structs */
    cm->transports = NULL;
    cm->initialized = 0;
    cm->reference_count = 1;

    char *tmp;
    if ((tmp = getenv("CMControlModule"))) {
	control_module = tmp;
    }

    if (control_module != NULL) {
	char *tmp = strdup(control_module);
	char *c;
	for (c = tmp; *c; ++c) *c = tolower(*c);
#ifdef HAVE_SYS_EPOLL_H
	if (strcmp(tmp, "epoll") == 0) {
	    cm->control_module_choice = "epoll";
	} else
#endif
	if (strcmp(tmp, "select") == 0) {
	    cm->control_module_choice = "select";
	} else {
	    fprintf(stderr, "Warning:  Specified CM/EVPath control module \"%s\" unknown or not built.\n", control_module);
	    /* force to default */
	    control_module = NULL;
	}
	free(tmp);
    }	
    if (control_module == NULL) {
#ifdef HAVE_SYS_EPOLL_H
	cm->control_module_choice = "epoll";
#else
	cm->control_module_choice = "select";
#endif
    }
    cm->control_list = CMControlList_create();
    thr_mutex_init(cm->exchange_lock);

    cm->locked = 0;
    cm->closed = 0;
    cm->abort_read_ahead = 0;
    cm->CMTrace_file = NULL;
    CMtrace_init(cm, EVerbose);
    CMinit_local_formats(cm);
    thr_mutex_init(cm->context_lock);

    cm->in_format_count = 0;
    cm->in_formats = INT_CMmalloc(1);

    cm->reg_format_count = 0;
    cm->reg_formats = INT_CMmalloc(1);

    cm->pending_request_max = 1;
    cm->pbio_requests = INT_CMmalloc(sizeof(struct _pending_format_requests));
    cm->pbio_requests[0].server_id = NULL;
    cm->pbio_requests[0].id_length = 0;
    cm->pbio_requests[0].condition = 0;
    cm->pbio_requests[0].top_request = 0;

    cm->connection_count = 0;
    cm->connections = INT_CMmalloc(1);
    cm->reg_user_format_count = 0;
    cm->reg_user_formats = INT_CMmalloc(1);
    cm->cm_buffer_list = NULL;
    cm->pending_data_queue = NULL;

    cm->contact_lists = NULL;
    cm->shutdown_functions = NULL;
    cm->perf_upcall = NULL;
    INT_CMadd_poll(cm, process_pending_queue, NULL);
#ifdef EV_INTERNAL_H
    CManager_lock(cm);
    EVPinit(cm);
    CManager_unlock(cm);
#endif
    return cm;
}

static void CMControlList_free(CManager cm, CMControlList cl);
static void remove_conn_from_CM(CManager cm, CMConnection conn);

static void
CManager_free(CManager cm)
{
    int i;
    CMbuffer list = NULL;

    INT_CMfree(cm->transports);
    cm->transports = NULL;
/*    free_FFSContext(cm->FFScontext);*/
    cm->FFScontext = NULL;
    INT_CMfree(cm->in_formats);

    for (i=0 ; i < cm->reg_format_count; i++) {
	INT_CMfree(cm->reg_formats[i]->format_name);
	INT_CMfree(cm->reg_formats[i]);
    }
    INT_CMfree(cm->reg_formats);

    /*
     *  Applications are expected to free the user contexts that 
     *  they request.  If they do this, there will be no user formats to
     *  free at this point.  (Doing so might result in double freeing.)
     */
    INT_CMfree(cm->reg_user_formats);

    INT_CMfree(cm->pbio_requests);

    INT_CMfree(cm->connections);

    thr_mutex_free(cm->exchange_lock);

    thr_mutex_free(cm->context_lock);

    if (cm->contact_lists != NULL) {
	i = 0;
	while(cm->contact_lists[i] != NULL) {
	    INT_CMfree_attr_list(cm, cm->contact_lists[i]);
	    i++;
	}
	INT_CMfree(cm->contact_lists);
    }
    list = cm->cm_buffer_list;
    i=0;
    while (list != NULL) {
	CMbuffer next = list->next;
	CMtrace_out(cm, CMBufferVerbose, "Final buffer disposition buf %d, %p, size %zd, ref_count %d\n", i++, list, list->size, list->ref_count);
	if (list->return_callback) {
	    (list->return_callback)(list->return_callback_data);
	} else {
	    INT_CMfree(list->buffer);
	}
	INT_CMfree(list);
	list = next;
    }
    cm->cm_buffer_list = NULL;
    if (cm->shutdown_functions) INT_CMfree(cm->shutdown_functions);
    INT_CMfree(cm->avail);
     INT_CMfree(cm);
 }

 extern void
 INT_CMinstall_perf_upcall(CManager cm, CMperf_upcall upcall)
 {
     cm->perf_upcall = upcall;
 }

 extern void
 INT_CManager_close(CManager cm)
 {
     CMControlList cl = cm->control_list;

     CMtrace_out(cm, CMFreeVerbose, "CManager %p closing, ref count %d\n", cm,
		 cm->reference_count);

     CMtrace_out(cm, CMFreeVerbose, "CMControlList close CL=%p current reference count will be %d, sdp = %p\n", 
		 cl, cl->cl_reference_count - 1, cl->select_data);
     INT_CMControlList_close(cl, cm);
     CMtrace_out(cm, CMFreeVerbose, "CMControlList CL=%p is closed\n", cl);

     while (cm->connection_count != 0) {
	 /* connections are moved down as they are closed... */
	 CMtrace_out(cm, CMFreeVerbose, "CManager in close, closing connection %p , ref count %d\n", cm->connections[0],
		     cm->connections[0]->conn_ref_count);
	 internal_connection_close(cm->connections[0]);
	 INT_CMConnection_failed(cm->connections[0]);
     }

     if (cm->shutdown_functions != NULL) {
	 func_entry *shutdown_functions = cm->shutdown_functions;
	 int i = 0;

	 while (shutdown_functions[i].func != NULL) {
	     if (shutdown_functions[i].task_type == SHUTDOWN_TASK) {
		 CMtrace_out(cm, CMFreeVerbose, "CManager calling shutdown function SHUTDOWN %d, %p\n", i, shutdown_functions[i].func);
		 shutdown_functions[i].func(cm, shutdown_functions[i].client_data);
		 shutdown_functions[i].task_type = NO_TASK;
	     }
	     i++;
	 }
     }
     cm->reference_count--;
     CMtrace_out(cm, CMFreeVerbose, "CManager %p ref count now %d\n", 
		 cm, cm->reference_count);
     if (cm->reference_count == 0) {
	 if (cm->shutdown_functions != NULL) {
	     int i = 0;
	     func_entry *shutdown_functions = cm->shutdown_functions;
	     cm->shutdown_functions = NULL;

	     while (shutdown_functions[i].func != NULL) {
		 i++;
	     }
	     i--;
	     for ( ; i >= 0; i--) {
		 if (shutdown_functions[i].task_type == FREE_TASK) {
		     CMtrace_out(cm, CMFreeVerbose, "CManager calling shutdown function FREE %d, %p\n", i, shutdown_functions[i].func);
		     shutdown_functions[i].func(cm, shutdown_functions[i].client_data);
		     shutdown_functions[i].func = NULL;
		 }
	     }
	     INT_CMfree(shutdown_functions);
	 }
	 CMtrace_out(cm, CMFreeVerbose, "Freeing CManager %p\n", cm);
	 cl->free_reference_count = 1;
	 CMControlList_free(cm, cl);
	 CManager_unlock(cm);
	 CManager_free(cm);
     } else {
	 CManager_unlock(cm);
     }
 }

 extern void
 internal_add_shutdown_task(CManager cm, CMPollFunc func, void *client_data, int task_type)
 {
     int func_count = 0;
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     if (cm->shutdown_functions == NULL) {
	 cm->shutdown_functions = 
	     INT_CMmalloc(sizeof(cm->shutdown_functions[0]) * 2);
     } else {
	 while (cm->shutdown_functions[func_count].func != NULL) {
	     func_count++;
	 }
	 cm->shutdown_functions = 
	     INT_CMrealloc(cm->shutdown_functions,
		       sizeof(cm->shutdown_functions[0]) * (func_count +2));
     }
     cm->shutdown_functions[func_count].func = func;
     cm->shutdown_functions[func_count].task_type = task_type;
     cm->shutdown_functions[func_count].client_data = client_data;
     func_count++;
     cm->shutdown_functions[func_count].func = NULL;
 }

 extern void
 INT_CMadd_shutdown_task(CManager cm, CMPollFunc func, void *client_data, int task_type)
 {
     internal_add_shutdown_task(cm, func, client_data, task_type);
 }

 static void
 add_conn_to_CM(CManager cm, CMConnection conn)
 {
     cm->connections = 
	 INT_CMrealloc(cm->connections, 
		   (cm->connection_count + 1) * sizeof(cm->connections[0]));
     cm->connections[cm->connection_count] = conn;
     INT_CMConnection_add_reference(conn);
     cm->connection_count++;
 }

 static void
 remove_conn_from_CM(CManager cm, CMConnection conn)
 {
     int i;
     int found = 0;
     for (i=0; i < cm->connection_count; i++) {
	 if (cm->connections[i] == conn) {
	     found++;
	     INT_CMConnection_dereference(conn);
	 } else if (found) {
	     /* copy down */
	     cm->connections[i-1] = cm->connections[i];
	 }
     }
     if (found == 0) {
	 fprintf(stderr, "Internal error, remove_conn_from_CM.  Not found\n");
     } else {
	 cm->connection_count--;
	 cm->abort_read_ahead = 1;
     }
 }

 static CMControlList
 CMControlList_create()
 {
     CMControlList new_list = (CMControlList) INT_CMmalloc(sizeof(CMControlList_s));
     new_list->select_initialized = 0;
     new_list->select_data = NULL;
     new_list->add_select = NULL;
     new_list->remove_select = NULL;
     new_list->server_thread =  (thr_thread_t)(intptr_t) NULL;
     new_list->network_blocking_function.func = NULL;
     new_list->network_polling_function.func = NULL;
     new_list->polling_function_list = NULL;
     new_list->cl_consistency_number = 1;
     new_list->cl_reference_count = 1;
     new_list->free_reference_count = 1;
     thr_mutex_init(new_list->list_mutex);
     new_list->condition_list = NULL;
     new_list->next_condition_num = 1;
     new_list->closed = 0;
     new_list->locked = 0;
     new_list->has_thread = 0;
     new_list->cond_polling = 0;
     return new_list;
 }

 CMConnection
 CMConnection_create(transport_entry trans, void *transport_data, 
		     attr_list conn_attrs)
 {
     static int first = 1;
     static int non_block_default = 0;
     static int read_thread_default = 0;
     int blocking_on_conn;
     CMConnection conn = INT_CMmalloc(sizeof(struct _CMConnection));
     if (first) {
	 char *value = getenv("CMNonBlockWrite");
	 first = 0;
	 if (value != NULL) {
	     sscanf(value, "%d", &non_block_default);
	     CMtrace_out(trans->cm, CMConnectionVerbose, "CM default blocking %d\n",
			 non_block_default);
	 }
	 value = getenv("CMReadThread");
	 if (value != NULL) {
	     sscanf(value, "%d", &read_thread_default);
	     CMtrace_out(trans->cm, CMConnectionVerbose, "CM default read thread %d\n",
			 read_thread_default);
	 }
     }
     conn->cm = trans->cm;
     conn->trans = trans;
     conn->transport_data = transport_data;
     conn->conn_ref_count = 1;
     conn->closed = 0;
     conn->failed = 0;
     conn->preloaded_formats = NULL;
     conn->remote_format_server_ID = 0;
     conn->remote_CManager_ID = 0;
     conn->handshake_condition = -1;
     conn->io_out_buffer = create_FFSBuffer();
     conn->close_list = NULL;
     conn->write_callback_len = 0;
     conn->write_callbacks = NULL;
     if (conn_attrs) {
	 CMadd_ref_attr_list(conn->cm, conn_attrs);
     }
     conn->attrs = conn_attrs;
     conn->attr_encode_buffer = create_AttrBuffer();

     conn->message_buffer = NULL;
     conn->buffer_full_point = 0;
     conn->buffer_data_end = 0;

     conn->characteristics = NULL;
     conn->write_pending = 0;
     conn->do_non_blocking_write = non_block_default;
     conn->XML_output = 0;
     conn->use_read_thread = read_thread_default; 

     if (get_int_attr(conn_attrs, CM_CONN_BLOCKING, &blocking_on_conn)) {
	 conn->do_non_blocking_write = !blocking_on_conn;
     }
     add_conn_to_CM(trans->cm, conn);
     CMtrace_out(trans->cm, CMFreeVerbose, "CMConnection_create %p \n",
		 conn);
     return conn;
 }

 extern attr_list
 INT_CMConnection_get_attrs(CMConnection conn)
 {
     return conn->attrs;
 }

 typedef struct {
     int size; /*stream length to probe bw. */
     int size_inc;/*size=size+size_inc if previous size is not adequate*/
     int successful_run; /*if we have at least 5 successful runs, 
			   size is set properly*/
     int failed_run;/*After parameters are set properly, if we have 3
		      contineous runs failed, we need to set the size again 
		      for the changed network condition*/
     CMConnection conn;
     attr_list attrs;
 } *bw_measure_data;

 static void
 do_bw_measure(CManager cm, void *client_data)
 {
     double bw;
     (void)cm;
     bw_measure_data data = (bw_measure_data) client_data;
     CManager_lock(cm);
     bw=INT_CMregressive_probe_bandwidth(data->conn, data->size, data->attrs);
     CManager_unlock(cm);

     /*Initialization phase*/
     if(bw<0 && data->successful_run<5){
	 data->size+=data->size_inc;
	 data->successful_run=0;
     }
     if(bw>=0 && data->successful_run<5) data->successful_run++; /*if measured correctly in several back-to-back measurements, increase data->successful_run up to 5. 5 indicates parameters are well tuned now*/

     /*After initialization, if network condition changed, and previously tuned parameters come to not adequate:*/
     if(bw<0 && data->successful_run>=5 && data->failed_run<5)
	 data->failed_run++;
     if(bw>=0)  data->failed_run=0; /*guard for contineous failures. */
     if(data->failed_run>=5){ /*need to tune parameters again now. */
	 data->successful_run=0;
	 data->failed_run=0;
     }

     CMtrace_out(data->conn->cm, CMLowLevelVerbose,"successful run: %d, failed run: %d, size: %d, bw: %f\n", data->successful_run, data->failed_run, data->size, bw);
     return;
 }    

 extern int
 INT_CMConnection_set_character(CMConnection conn, attr_list attrs)
 {
     ssize_t interval_value;
     if (attrs == NULL) return 0;
     if (get_long_attr(attrs, CM_BW_MEASURE_INTERVAL, &interval_value)) {
	 bw_measure_data data;
	 int previous_interval;
	 CMTaskHandle task = NULL;
	 if ((interval_value <= 1) || (interval_value > 60*60*8)) {
	     printf("Bad CM_BW_MEASURE_INTERVAL, %zd seconds\n",
		    interval_value);
	     return 0;
	 }

	 CMtrace_out(conn->cm, CMLowLevelVerbose,"CM_BW_MEASURE_INTERVAL set, interval is %zd\n", interval_value);
	 if (conn->characteristics && 
	     (get_int_attr(conn->characteristics, CM_BW_MEASURE_INTERVAL,
			   &previous_interval) != 0)) {
	     CMTaskHandle prior_task = NULL;
	     if (interval_value >= previous_interval) {
		 CMtrace_out(conn->cm, CMLowLevelVerbose,"CM_BW_MEASURE_INTERVAL prior interval is %d, no action.\n", previous_interval);
		 return 1;
	     }
	     CMtrace_out(conn->cm, CMLowLevelVerbose,"CM_BW_MEASURE_INTERVAL prior interval is %d, killing prior task.\n", previous_interval);
	     get_long_attr(conn->characteristics, CM_BW_MEASURE_TASK,
			   (ssize_t*)(intptr_t)&prior_task);
	     if (prior_task) {
		 INT_CMremove_task(prior_task);
		 set_long_attr(conn->characteristics, CM_BW_MEASURE_TASK, (intptr_t)0);
	     }
	 }
	 data = malloc(sizeof(*data));
	 data->size=data->size_inc=-1;

	 /*Get attr about size, size_inc from attributes. */
	 get_int_attr(attrs, CM_BW_MEASURE_SIZE, &(data->size));
	 if(data->size<1024)
	     data->size=1024;
	 get_int_attr(attrs, CM_BW_MEASURE_SIZEINC, &(data->size_inc));
	 if(data->size_inc<1024)
	     data->size_inc=1024;

	 data->successful_run=0;
	 data->failed_run=0;

	 /*app set attr about N and repeat time. store in data->attrs automically and pass to regressive_probe_bandwidth. 
	  */
	 data->conn = conn;
	 data->attrs = CMattr_copy_list(conn->cm, attrs);
	 /* do one task almost immediately */
	 task = INT_CMadd_delayed_task(conn->cm, 0, 1000, do_bw_measure, 
				   (void*)data);
	 free(task);
	 /* schedule tasks periodically */
	 task = INT_CMadd_periodic_task(conn->cm, (int)interval_value, 0, 
				    do_bw_measure, (void*)data);
	 if (conn->characteristics == NULL) {
	     conn->characteristics = CMcreate_attr_list(conn->cm);
	 }
	 set_int_attr(conn->characteristics, CM_BW_MEASURE_INTERVAL,
		      (int)interval_value);
	 set_long_attr(conn->characteristics, CM_BW_MEASURE_TASK, (intptr_t)task);

	 return 1;
     }
     return 0;
 }

 extern CMConnection 
 INT_CMget_indexed_conn(CManager cm, int i)
 {
     if (i>=0 && i<cm->connection_count) {
	 if (cm->connections[i] != NULL) {
	     return cm->connections[i];
	 } else {
	     CMtrace_out(cm, CMConnectionVerbose,
			 "cm->connection[%d] is NULL. INT_CMget_indexed_conn\n", i);
	     return NULL;
	 }
     } else {
	 CMtrace_out(cm, CMConnectionVerbose,
		     "Invalid index. i=%d. INT_CMget_indexed_conn\n", i);
	 return NULL;
     }
 }

 void
 INT_CMConnection_dereference(CMConnection conn)
 {
     conn->conn_ref_count--;
     if (conn->conn_ref_count > 0) {
	 CMtrace_out(conn->cm, CMFreeVerbose, "CM - Dereference connection %p, ref count now %d\n",
		     (void*)conn, conn->conn_ref_count);
	 return;
     }
     if (conn->conn_ref_count < 0) {
	 CMtrace_out(conn->cm, CMFreeVerbose, "CM - connection reference count less than 0, conn %p\n", conn);
	 return;   /*  BAD! */
     }
     CMtrace_out(conn->cm, CMFreeVerbose, "CM - Shut down connection %p\n",
		 (void*)conn);
     if (conn->write_pending) {
	 wait_for_pending_write(conn);
     }
     conn->closed = 1;
     if (conn->failed == 0) {
	 CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection failed with no dereference %p\n", conn);
	 INT_CMConnection_failed(conn);
     }
     CMtrace_out(conn->cm, CMFreeVerbose, "CM - Dereference connection %p FREEING\n", (void*)conn);
     if (conn->write_callbacks) INT_CMfree(conn->write_callbacks);
     INT_CMfree(conn->preloaded_formats);
     INT_CMfree_attr_list(conn->cm, conn->attrs);
     free_FFSBuffer(conn->io_out_buffer);
     free_AttrBuffer(conn->attr_encode_buffer);
 #ifdef EV_INTERNAL_H
     INT_EVforget_connection(conn->cm, conn);
 #endif
     INT_CMfree(conn);
 }

void
INT_CMConnection_failed(CMConnection conn)
{
    CMTaskHandle prior_task = NULL;
    if (conn->failed) return;
    conn->failed = 1;
    transport_wake_any_pending_write(conn);

    assert(CManager_locked(conn->cm));
    CMtrace_out(conn->cm, CMFreeVerbose, "CMConnection failed conn=%p\n", 
		conn);
    CMconn_fail_conditions(conn);
    conn->trans->shutdown_conn(&CMstatic_trans_svcs, conn->transport_data);
    get_long_attr(conn->characteristics, CM_BW_MEASURE_TASK, 
		  (ssize_t*)(intptr_t)&prior_task);
    if (prior_task) {
	INT_CMremove_task(prior_task);
	set_long_attr(conn->characteristics, CM_BW_MEASURE_TASK, (intptr_t)0);
    }
    if (conn->close_list) {
	CMCloseHandlerList list = conn->close_list;
	conn->close_list = NULL;
	while (list != NULL) {
	    CMCloseHandlerList next = list->next;
	    if (! conn->closed ) {
		CMtrace_out(conn->cm, CMConnectionVerbose, 
			    "CM - Calling close handler %p for connection %p\n",
			    (void*) list->close_handler, (void*)conn);
		CManager_unlock(conn->cm);
		list->close_handler(conn->cm, conn, list->close_client_data);
		CManager_lock(conn->cm);
	    }
	    INT_CMfree(list);
	    list = next;
	}
    }
    conn->closed = 1;
    remove_conn_from_CM(conn->cm, conn);
}

 void
 internal_connection_close(CMConnection conn)
 {
     CMtrace_out(conn->cm, CMFreeVerbose, "internal_connection_close conn=%p ref count is %d\n", 
		 conn, conn->conn_ref_count);
     conn->closed = 1;
 }

 void
 INT_CMConnection_close(CMConnection conn)
 {
     internal_connection_close(conn);
     CMtrace_out(conn->cm, CMFreeVerbose, "User CMConnection close conn=%p ref count will be %d\n", 
		 conn, conn->conn_ref_count - 1);
     INT_CMConnection_dereference(conn);
 }

 void
 INT_CMconn_register_close_handler(CMConnection conn, CMCloseHandlerFunc func, 
				   void *client_data)
 {
     CMCloseHandlerList *lastp = &conn->close_list;
     CMCloseHandlerList entry = INT_CMmalloc(sizeof(*entry));
     while (*lastp != NULL) lastp = &((*lastp)->next);
     entry->close_handler = func;
     entry->close_client_data = client_data;
     entry->next = NULL;
     *lastp = entry;
 }

 static void
 INT_CMControlList_close(CMControlList cl, CManager cm)
 {
     void *status;
     cl->cl_reference_count--;
     cl->closed = 1;

     (cl->stop_select)((void*)&CMstatic_trans_svcs, &cl->select_data);
     if ((cl->has_thread > 0) && (cl->server_thread != thr_thread_self())){
	     (cl->wake_select)((void*)&CMstatic_trans_svcs,
			       &cl->select_data);
     }	
     if ((cl->has_thread > 0) && (cl->server_thread != thr_thread_self())){
	 (cl->stop_select)((void*)&CMstatic_trans_svcs,
			   &cl->select_data);

	 (cl->wake_select)((void*)&CMstatic_trans_svcs,
			   &cl->select_data);
	 CManager_unlock(cm);
	 thr_thread_join(cl->server_thread, &status);
	 CManager_lock(cm);
	 cl->has_thread = 0;
     }
 }


 void
 CMwake_server_thread(CManager cm)
 {
     CMControlList cl = cm->control_list;
     (cl->wake_select)((void*)&CMstatic_trans_svcs, &cl->select_data);
 }

 extern void
 internal_condition_free(CMControlList cl);

 static void
 CMControlList_free(CManager cm, CMControlList cl)
 {
     cl->free_reference_count--;
     if (CMtrace_val[CMFreeVerbose]) {
	 fprintf(cm->CMTrace_file, "CMControlList_free, %p, ref count now %d\n", cl,
		cl->free_reference_count);
     }
     if(cl->free_reference_count == 0) {
	 if (CMtrace_val[CMFreeVerbose]) {
	     fprintf(cm->CMTrace_file, "CMControlList_free freeing %p\n", cl);
	 }
	 if (cl->polling_function_list != NULL) {
	     INT_CMfree(cl->polling_function_list);
	 }
	 thr_mutex_free(cl->list_mutex);
	 internal_condition_free(cl);
	 INT_CMfree(cl);
     }
 }

extern void
INT_CMget_qual_hostname(CManager cm, char *buf, int len)
{
    get_IP_config(buf, len, NULL, NULL, NULL, NULL, NULL,
		  CMstatic_trans_svcs.trace_out, cm);
}

extern void
INT_CMget_port_range(CManager cm, int *high_bound, int *low_bound)
{
    get_IP_config(NULL, 0, NULL, low_bound, high_bound, NULL, NULL,
		  CMstatic_trans_svcs.trace_out, cm);
}

extern int
INT_CMget_self_ip_addr(CManager cm)
{
    int IP;
    get_IP_config(NULL, 0, &IP, NULL, NULL, NULL, NULL,
		  CMstatic_trans_svcs.trace_out, cm);
    return IP;
}

extern char *
INT_CMget_ip_config_diagnostics(CManager cm)
{
    return IP_get_diagnostics(cm, CMstatic_trans_svcs.trace_out);
}

 #define CURRENT_HANDSHAKE_VERSION 1

 static 
 int
 transport_is_reliable(CMConnection conn)
 {
     attr_list list;
     int ret = 0;
     if (conn->trans->get_transport_characteristics == NULL) {
	 return 0; /* don't know */
     }
     list = conn->trans->get_transport_characteristics(conn->trans, &CMstatic_trans_svcs, 
						       conn->trans->trans_data);
     (void) get_int_attr(list, CM_TRANSPORT_RELIABLE, &ret);

     free_attr_list(list);

     return ret;
 }

 static
 void
 send_and_maybe_wait_for_handshake(CManager cm, CMConnection conn)
 {
     struct FFSEncodeVec tmp_vec[1];
     int reliable = transport_is_reliable(conn);
     int msg[5], actual;
     if (!cm->FFSserver_identifier) cm->FFSserver_identifier = -1;
     msg[0] = 0x434d4800;  /* CMH\0 */
     msg[1] = (CURRENT_HANDSHAKE_VERSION << 24) + sizeof(msg);
     msg[2] = cm->FFSserver_identifier;
     msg[3] = 5;  /* not implemented yet */
     msg[4] = 0;  /* not implemented yet */
     if (conn->remote_format_server_ID != 0) {
	 /* set high bit if we already have his ID */
	 msg[3] |= 0x80000000;
     }
     tmp_vec[0].iov_base = &msg;
     tmp_vec[0].iov_len = sizeof(msg);
     CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - sending handshake\n");
     if ((conn->remote_format_server_ID == 0) && reliable) {
	 /* we will await his respone */
	 conn->handshake_condition = INT_CMCondition_get(cm, conn);
     }
     actual = conn->trans->writev_func(&CMstatic_trans_svcs, 
				       conn->transport_data, 
				       &tmp_vec[0], 1, NULL);

     CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - after handshake, pending is %d\n", conn->write_pending);
     if (conn->write_pending) {
	 wait_for_pending_write(conn);
     }
     if (actual != 1) {
	 printf("handshake write failed\n");
     }
     if ((conn->remote_format_server_ID == 0) && reliable) {
	 CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - waiting for handshake response\n");
	 INT_CMCondition_wait(cm, conn->handshake_condition);
     }
 }

static
void
timeout_conn(CManager cm, void *client_data)
{
    INT_CMCondition_fail(cm, (int)(intptr_t) client_data);
}

 static
 CMConnection
 try_conn_init(CManager cm, transport_entry trans, attr_list attrs)
 {
     CMConnection conn = NULL;
     if (trans->initiate_conn) {
	 conn = trans->initiate_conn(cm, &CMstatic_trans_svcs,
				     trans, attrs);
     } else if (trans->initiate_conn_nonblocking) {
	 int result;
	 long wait_condition = INT_CMCondition_get(cm, NULL);
	 CMTaskHandle task = INT_CMadd_delayed_task(cm, 5, 0, timeout_conn,
						    (void*)(intptr_t)wait_condition);
	 if (CMtrace_on(cm, CMConnectionVerbose)) {
	     char *attr_str = attr_list_to_string(attrs);
	     CMtrace_out(cm, CMConnectionVerbose, 
			 "CM - Try to establish connection %p - %s, wait condition %ld\n", (void*)conn,
			 attr_str, wait_condition);
	     INT_CMfree(attr_str);
	 }
	 void *client_data = trans->initiate_conn_nonblocking(cm, &CMstatic_trans_svcs,
						 trans, attrs, wait_condition);
	 // Upon wake, condition will have either been signaled (return 1) or failed (return 0)
	 result = INT_CMCondition_wait(cm, wait_condition);
	 CMtrace_out(cm, CMConnectionVerbose, 
		     "CM - CMConnection wait returned, result %d\n", result);
	 if (result == 1) {
	     INT_CMremove_task(task);
	 }
	 conn = trans->finalize_conn_nonblocking(cm,  &CMstatic_trans_svcs,
						 trans, client_data, result);
     } else {
	 assert(0);
     }
     if (conn != NULL) {
	 if (CMtrace_on(conn->cm, CMConnectionVerbose)) {
	     char *attr_str = attr_list_to_string(attrs);
	     CMtrace_out(conn->cm, CMConnectionVerbose, 
			 "CM - Establish connection %p - %s\n", (void*)conn,
			 attr_str);
	     INT_CMfree(attr_str);
	 }
	 if (conn->use_read_thread) {
	     INT_CMstart_read_thread(conn);
	 }
	 send_and_maybe_wait_for_handshake(cm, conn);
     }
     return conn;
 }

 CMConnection
 CMinternal_initiate_conn(CManager cm, attr_list attrs)
 {
     transport_entry *trans_list;
     char *chosen_transport = NULL;

     assert(CManager_locked(cm));

     if (attrs) {
	 attrs = split_transport_attributes(attr_copy_list(attrs));
	 get_string_attr(attrs, CM_TRANSPORT, &chosen_transport);
     }
     if (chosen_transport != NULL) {
	 if (load_transport(cm, chosen_transport, 1) == 0) {
	     CMtrace_out(cm, CMConnectionVerbose,
			 "Failed to load transport \"%s\".  Revert to default.\n",
			 chosen_transport);
	     chosen_transport = NULL;
	 }
     }
     trans_list = cm->transports;
     if (chosen_transport == NULL) {
	 CMtrace_out(cm, CMConnectionVerbose,
		     "INT_CMinitiate_conn no transport attr found\n");

	 while ((trans_list != NULL) && (*trans_list != NULL)) {
	     CMConnection conn;
	     conn = try_conn_init(cm, *trans_list, attrs);
	     if (conn != NULL) {
		 if (attrs) free_attr_list(attrs);
		 return conn;
	     }
	     trans_list++;
	 }
     } else {
	 CMtrace_out(cm, CMConnectionVerbose,
		     "INT_CMinitiate_conn looking for transport \"%s\"\n", 
		     chosen_transport);
	 while ((trans_list != NULL) && (*trans_list != NULL)) {
	     if (strcmp((*trans_list)->trans_name, chosen_transport) == 0) {
		 CMConnection conn = try_conn_init(cm, *trans_list, attrs);
		 if (attrs) free_attr_list(attrs);
		 return conn;
	     }
	     trans_list++;
	 }
	 CMtrace_out(cm, CMConnectionVerbose,
		     "INT_CMinitiate_conn transport \"%s\" not found - no connection\n", 
		     chosen_transport);
	 if (attrs) free_attr_list(attrs);
	 return NULL;
     }

     if (attrs) free_attr_list(attrs);
     return NULL;
 }

 static void
 fdump_CMConnection(FILE *out, CMConnection conn)
 {
     if (conn == NULL) {
	 fprintf(out, "CMConnection NULL\n");
	 return;
     }
     fprintf(out, "CMConnection %p, reference count %d, closed %d\n\tattrs : ", 
	    conn, conn->conn_ref_count, conn->closed);
     fdump_attr_list(out, conn->attrs);
     fprintf(out, "\tbuffer_full_point %zd, current buffer_end %zd\n", 
	     conn->buffer_full_point, conn->buffer_data_end);
     fprintf(out, "\twrite_pending %d\n", conn->write_pending);
 }

 /*static void
 dump_CMConnection(CMConnection conn)
 {
     fdump_CMConnection(stdout, conn);
     }*/

 CMConnection
 INT_CMinitiate_conn(CManager cm, attr_list attrs)
 {
     CMConnection conn;
     if (!cm->initialized) CMinitialize(cm);
     if (CMtrace_on(cm, CMConnectionVerbose)) {
	 fprintf(cm->CMTrace_file,"Doing CMinitiate_conn\n");
     }
     conn = CMinternal_initiate_conn(cm, attrs);
     if (CMtrace_on(cm, CMConnectionVerbose)) {
	 if (conn != NULL) {
	     fdump_CMConnection(cm->CMTrace_file, conn);
	 } else {
	     fprintf(cm->CMTrace_file, "NULL\n");
	 }
     }
     return conn;
 }

 void
 INT_CMConnection_add_reference(CMConnection conn)
 {
     conn->conn_ref_count++;
     CMtrace_out(conn->cm, CMFreeVerbose, "Add reference to connection %p, value is now %d\n", conn, conn->conn_ref_count);
 }

 CMConnection
 CMinternal_get_conn(CManager cm, attr_list attrs)
 {
     int i;
     CMConnection conn = NULL;
     assert(CManager_locked(cm));
     if (CMtrace_on(cm, CMConnectionVerbose)) {
	 fprintf(cm->CMTrace_file, "In CMinternal_get_conn, attrs ");
	 if (attrs) fdump_attr_list(cm->CMTrace_file, attrs); else fprintf(cm->CMTrace_file, "\n");
     }
     for (i=0; i<cm->connection_count; i++) {
	 CMConnection tmp = cm->connections[i];
	 if (tmp->closed || tmp->failed) continue;
	 if (tmp->trans->connection_eq(cm, &CMstatic_trans_svcs,
					tmp->trans, attrs,
					tmp->transport_data)) {

	     CMtrace_out(tmp->cm, CMFreeVerbose, "internal_get_conn found conn=%p ref count will be %d\n", 
			 tmp, tmp->conn_ref_count +1);
	     CMtrace_out(tmp->cm, CMConnectionVerbose, "internal_get_conn found conn=%p ref count will be %d\n", 
			 tmp, tmp->conn_ref_count +1);
	     tmp->conn_ref_count++;
	     conn = tmp;
	     break;
	 }
     }
     if (conn == NULL) {
	 if (CMtrace_on(cm, CMConnectionVerbose)) {
	     fprintf(cm->CMTrace_file, "In CMinternal_get_conn, no existing connection found, initiating\n");
	 }
	 conn = CMinternal_initiate_conn(cm, attrs);
	 if (conn) {
	     CMtrace_out(conn->cm, CMFreeVerbose, "internal_get_conn initiated connection %p ref count now %d\n", 
			 conn, conn->conn_ref_count);
	 }
     }
     if (conn != NULL) {
	 CMtrace_out(conn->cm, CMFreeVerbose, "internal_get_conn returning conn=%p ref count %d\n", 
		     conn, conn->conn_ref_count);
     }
     if (CMtrace_on(cm, CMConnectionVerbose)) {
	 fprintf(cm->CMTrace_file, "CMinternal_get_conn returning ");
	 if (conn != NULL) {
	     fdump_CMConnection(cm->CMTrace_file, conn);
	 } else {
	     fprintf(cm->CMTrace_file, "NULL\n");
	 }
     }
     return conn;
 }

 CMConnection
 INT_CMget_conn(CManager cm, attr_list attrs)
 {
     CMConnection conn;
     if (!cm->initialized) CMinitialize(cm);
     conn = CMinternal_get_conn(cm, attrs);
     return conn;
 }

 int
 INT_CMcontact_self_check(CManager cm, attr_list attrs)
 {
     transport_entry *trans_list;
     if (!cm->initialized) CMinitialize(cm);
     trans_list = cm->transports;
     while ((trans_list != NULL) && (*trans_list != NULL)) {
	 int result = 0;
	 result = (*trans_list)->self_check(cm, &CMstatic_trans_svcs, 
					    *trans_list, attrs);
	 if (result) return result;
	 trans_list++;
     }
     return 0;
 }

 extern CMbuffer
 cm_create_transport_buffer(CManager cm, void *buffer, ssize_t length)
 {
     CMbuffer tmp;
     (void)cm;
     tmp = INT_CMmalloc(sizeof(*tmp));
     memset(tmp, 0, sizeof(*tmp));
     tmp->buffer = buffer;
     tmp->size = length;
     tmp->ref_count = 1;
     CMtrace_out(cm, CMBufferVerbose, "Creating buffer %p, ref_count is %d\n", tmp, tmp->ref_count);
 //   This should just return the buffer... not update the link list.  That's handled in the calling routine.
 //    tmp->next = cm->cm_buffer_list;
 //    cm->cm_buffer_list = tmp;
     return tmp;
 }

 extern CMbuffer
 cm_create_transport_and_link_buffer(CManager cm, void *buffer, ssize_t length)
 {
     CMbuffer tmp;
     tmp = INT_CMmalloc(sizeof(*tmp));
     memset(tmp, 0, sizeof(*tmp));
     tmp->buffer = buffer;
     tmp->size = length;
     tmp->ref_count = 1;
     CMtrace_out(cm, CMBufferVerbose, "Create and link buffer %p, ref_count is %d\n", tmp, tmp->ref_count);
     tmp->next = cm->cm_buffer_list;
     cm->cm_buffer_list = tmp;
     return tmp;
 }

 /* alloc temporary buffer for CM use */
 extern CMbuffer
 cm_get_data_buf(CManager cm, ssize_t length)  
 {
     int buffer_count = 0;
     CMbuffer tmp = cm->cm_buffer_list;

     CMtrace_out(cm, CMBufferVerbose, "cm_get_data_buf called with len %zu\n",
		 length);
     while (tmp != NULL) {
	 CMtrace_out(cm, CMBufferVerbose, "  buffer %d %p, size is %zd, data %p, ref_count %d\n",
		     buffer_count, tmp, tmp->size, tmp->buffer, tmp->ref_count);
	 buffer_count++;
	 tmp = tmp->next;
     }

     tmp = cm->cm_buffer_list;
     buffer_count = 0;
     /* first baseline consistency check */
     while (tmp != NULL) {
	 if (tmp->ref_count < 0) {
	     CMtrace_out(cm, CMBufferVerbose, "cm_get_data_buf buffer %p, ref_count is %d, should not be negative\n", tmp, tmp->ref_count);
	 }
	 buffer_count++;
	 tmp = tmp->next;
     }
     /* first look for a buffer big enough, but not too big */
     tmp = cm->cm_buffer_list;
     while (tmp != NULL) {
	 if (tmp->ref_count <= 0) {
	     if ((tmp->size >= (size_t)length) && ((tmp->size/10) < (size_t)length)) {
		 CMtrace_out(cm, CMBufferVerbose, "cm_get_data_buf called len %zu, return existing %p, next %p, count %d\n",
			     length, tmp, tmp->next, buffer_count);
		 tmp->ref_count = 1;
		 return tmp;
	     }
	 }
	 tmp = tmp->next;
     }		
     /* ok, we'll settle for way too big, but realloc it down */
     tmp = cm->cm_buffer_list;
     while (tmp != NULL) {
	 if (tmp->ref_count <= 0) {
	     if ((tmp->size >= (size_t)length)) {
		 char *t = INT_CMrealloc(tmp->buffer, length);
		 if (t == NULL) {
		     return NULL;
		 }
		 tmp->buffer = t;
		 tmp->size = length;
		 tmp->ref_count = 1;
		 CMtrace_out(cm, CMBufferVerbose, "      cm_get_data_buf resizing down!  return is %p\n", tmp);
		 return tmp;
	     }
	 }
	 tmp = tmp->next;
     }		
     tmp = cm->cm_buffer_list;
     /* well, look for a small one to realloc up */
     while (tmp != NULL) {
	 if (tmp->ref_count <= 0) {
	     if (tmp->size <= (size_t)length) {
		 char *t = INT_CMrealloc(tmp->buffer, length);
		 if (t == NULL) {
		     return NULL;
		 }
		 tmp->buffer = t;
		 tmp->size = length;
		 tmp->ref_count = 1;
		 CMtrace_out(cm, CMBufferVerbose, "      cm_get_data_buf resizingup!  return is %p\n", tmp);
		 return tmp;
	     }
	 }
	 tmp = tmp->next;
     }
     tmp = cm_create_transport_buffer(cm, INT_CMmalloc(length), length);
     tmp->ref_count = 1;
     tmp->next = cm->cm_buffer_list;
     cm->cm_buffer_list = tmp;
     CMtrace_out(cm, CMBufferVerbose, "cm_get_data_buf create new len %zu, return %p, count %d\n",
		 length, tmp, buffer_count);
     return tmp;
 }

 /* realloc temporary buffer for CM use */
 extern CMbuffer
 cm_extend_data_buf(CManager cm, CMbuffer tmp, size_t length)  
 {
     (void)cm;
     if (tmp->size < length) {
	 char *t = INT_CMrealloc(tmp->buffer, length);
	 if (t == NULL) {
	     return NULL;
	 }
	 tmp->buffer = t;
	 tmp->size = length;
     }
     return tmp;
 }

 /* CM says that it is done with temporary buffer */
 extern void
 cm_return_data_buf(CManager cm, CMbuffer cmb)
 {
     cmb->ref_count--;
     CMtrace_out(cm, CMBufferVerbose, "cm_return_data_buf buffer %p, callback %p, ref_count is now %d\n", cmb, cmb->return_callback, cmb->ref_count);
     if ((cmb->ref_count == 0) && (cmb->return_callback != NULL)) {
	 CMbuffer last = NULL, tmp = cm->cm_buffer_list;
	 /* UNLINK */
	 CMtrace_out(cm, CMBufferVerbose, "cm_return_data_buf --- Unlinking %p cmb\n", cmb);
	 while (tmp != NULL) {
	     if (tmp != cmb) {
		 last = tmp;
		 tmp = tmp->next;
		 continue;
	     }
	     /* remove the buffer from CM's list */
	     if (last == NULL) {
		 cm->cm_buffer_list = tmp->next;
	     } else {
		 last->next = tmp->next;
	     }
	     tmp = tmp->next;
	     (cmb->return_callback)(cmb->return_callback_data);
	     free(cmb);
	     break;
	 }
     }
 }

 static CMbuffer
 cm_buffer_lookup(CManager cm, void *buffer)
 {
     CMbuffer tmp = cm->cm_buffer_list;
     while (tmp != NULL) {
	 if ((tmp->buffer <= buffer) && 
	     ((char*)buffer < ((char*)tmp->buffer + tmp->size))){
	     return tmp;
	 }
	 tmp = tmp->next;
     }
     return NULL;
 }

 static CMbuffer
 cm_buffer_dump_list(CManager cm)
 {
     CMbuffer tmp = cm->cm_buffer_list;
     printf("Known CM buffers are:\n");
     while (tmp != NULL) {
	 printf("Buffer begin %p, size %zd, end %p\n",
		tmp->buffer, tmp->size, (char*)tmp->buffer + tmp->size);
	 tmp = tmp->next;
     }
     return NULL;
 }

 int (*cm_write_hook)(size_t) = (int (*)(size_t)) NULL;
 int (*cm_preread_hook)(size_t,char*) = (int (*)(size_t, char*)) NULL;
 void (*cm_postread_hook)(size_t,char*) = (void (*)(size_t, char*)) NULL;
 void (*cm_last_postread_hook)() = (void (*)()) NULL;
 static size_t CMact_on_data(CMConnection conn, CMbuffer cm_buffer, char *buffer, size_t length);

 static void process_pending_queue(CManager cm, void *junk)
 {
     /* shortcircuit if no data, no lock */
     if (!cm->pending_data_queue) return;

     CManager_lock(cm);
     while (cm->pending_data_queue) {
	 pending_queue entry = cm->pending_data_queue;
	 size_t result;
	 cm->pending_data_queue = entry->next;
	 result = CMact_on_data(entry->conn, entry->buffer, entry->buffer->buffer, entry->length);
	 if (result != 0) {
	     printf("in process pending, CMact_on_data returned %zd\n", result);
	 }
	 cm_return_data_buf(cm, entry->buffer);
	 free(entry);
     }
     CManager_unlock(cm);
 }

 static void add_buffer_to_pending_queue(CManager cm, CMConnection conn, CMbuffer buf, long length)
 {
     assert(CManager_locked(cm));
     pending_queue entry = malloc(sizeof(struct pending_queue_entry));
     entry->next = NULL;
     entry->conn = conn;
     entry->buffer = buf;
     entry->length = length;
     if (!cm->pending_data_queue) {
	 cm->pending_data_queue = entry;
     } else {
	 pending_queue last = cm->pending_data_queue;
	 pending_queue tmp = last->next;
	 while (tmp) {
	     last = tmp;
	     tmp = tmp->next;
	 }
	 last->next = entry;
     }
     CMwake_server_thread(cm);
 }

 static
 CMbuffer
 fill_cmbuffer(CManager cm, char *buf, size_t length)
 {
     CMbuffer ret = cm_get_data_buf(cm, length);
     memcpy(ret->buffer, buf, length);
     return ret;
 }

 extern void CMDataAvailable(transport_entry trans, CMConnection conn)
 {
     CManager cm = conn->cm;
     int do_read = 1;
     int read_msg_count = 0;
     size_t read_byte_count = 0;
     size_t result;
     static int first = 1;
     static int read_ahead_msg_limit = 50;
     static size_t read_ahead_byte_limit = 1024*1024*1024;
     static int use_blocking_reads = 1;
     int first_four = 0;
     char *tmp_message_buffer = NULL;
     ssize_t data_length;
     CMbuffer message_buffer;
     size_t buffer_full_point, buffer_data_end;

     /* called from the transport, grab the locks */
     if (first) {
	 char *tmp;
	 first = 0;
	 tmp = getenv("CMReadAheadMsgLimit");
	 if (tmp != NULL) {
	     if (sscanf(tmp, "%d", &read_ahead_msg_limit) != 1) {
		 printf("Read ahead msg limit \"%s\" not parsed\n", tmp);
	     }
	 }
	 tmp = getenv("CMReadAheadByteLimit");
	 if (tmp != NULL) {
	     if (sscanf(tmp, "%zd", &read_ahead_byte_limit) != 1) {
		 printf("Read ahead byte limit \"%s\" not parsed\n", tmp);
	     }
	 }
	 tmp = getenv("CMBlockingReads");
	 if (tmp != NULL) {
	     use_blocking_reads = atoi(tmp);
	 }
     }

     if (conn->buffer_full_point == 0) {
	 /* nothing saved from last invocation */
	 message_buffer = NULL;
	 buffer_full_point = 0;
	 buffer_data_end = 0;
     } else {
	 /* message_buffer, buffer_full_point and buffer_data_end cached from last call to CMDataAvailable */
	 message_buffer = conn->message_buffer;
	 conn->message_buffer = NULL;
	 buffer_full_point = conn->buffer_full_point;
	 conn->buffer_full_point = 0;
	 buffer_data_end = conn->buffer_data_end;
	 conn->buffer_data_end = 0;
     }
  start_read:
     if (conn->failed) {
	 return;
     }
     if (buffer_full_point == 0) {
	 buffer_full_point = 4; /* read first 4 bytes first thing */
	 buffer_data_end = 0; /* no bytes read yet */
	 first_four = 1;
     }
     if (trans->read_to_buffer_func) {
	 CMtrace_out(cm, CMLowLevelVerbose, "CMdata continuing read, already have %zd bytes, trying to read total %zd\n", buffer_data_end, buffer_full_point);
     } else {
	 CMtrace_out(cm, CMLowLevelVerbose, "CMdata block read beginning\n");
     }	
     if (buffer_full_point < HEADER_BUFFER_SIZE) {
	 tmp_message_buffer = &conn->header_buffer[0];
     } else {
	 if (message_buffer == NULL) {
	     /* we had data in the header buffer, but need more space */
	     message_buffer = cm_get_data_buf(cm, buffer_full_point);
	     memcpy(message_buffer->buffer, &conn->header_buffer[0],
		    buffer_data_end);
	     tmp_message_buffer = message_buffer->buffer;
	 } else {
	     /* make sure buffer is big enough */
	     cm_extend_data_buf(cm, message_buffer, buffer_full_point);
	     tmp_message_buffer = message_buffer->buffer;
	 }
     }
     if (cm_preread_hook) {
	 do_read = cm_preread_hook(buffer_full_point - buffer_data_end, tmp_message_buffer);
     }
     if (do_read) {
	 if (trans->read_to_buffer_func) {
	 /* 
	  * read_to_buffer functionality present means we can read the transport directly to any memory area
	  * This gives us the most flexibility in managing our data.  
	  *   At this point, tmp_message_buffer is the char* buffer to which we want to read
	  *   If non-NULL, message_buffer is the CMBuffer that holds this memory
	  *   buffer_data_end is how much data is already in the buffer
	  *   buffer_full_point is how much data we think we need for the message to be complete
	  */
	     size_t read_len = buffer_full_point - buffer_data_end;
	     char *read_target_buf = tmp_message_buffer + buffer_data_end;
	     /* 
	      * non blocking is True only if :
	      *    - we're reading the first four bytes (first_four is true) and use_read_thread is false
	      *    - or use_blocking_reads is false and use_read_thread is false
	      */
	     int non_blocking = first_four || !use_blocking_reads;
	     ssize_t actual;
	     if (conn->use_read_thread) {
		 non_blocking = 0;
		 CManager_unlock(cm);
	     }
	     if (read_len == 0) {
		 /* should never happen that we are here but don't need more data */
		 printf("Seriously bad shit\n");
	     }
	     actual = trans->read_to_buffer_func(&CMstatic_trans_svcs, 
						 conn->transport_data, 
						 read_target_buf, read_len, non_blocking);
	     if (conn->use_read_thread) {
		 CManager_lock(cm);
	     }
	     if (actual == -1) {
		 CMtrace_out(cm, CMLowLevelVerbose, 
			     "CMdata read failed, actual %zu, failing connection %p\n", actual, conn);
		 CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection failed read_len with dereference %p\n", conn);
		 INT_CMConnection_failed(conn);
		 return;
	     }
	     buffer_data_end += actual;
	     if (actual < (ssize_t)read_len) {
		 /* partial read, we know we don't have enough data now, roll on */
		 CMtrace_out(cm, CMLowLevelVerbose, 
			     "CMdata read partial, got %zu\n", actual);
		 /* save state and return */
		 conn->message_buffer = message_buffer;
		 conn->buffer_full_point = buffer_full_point;
		 conn->buffer_data_end = buffer_data_end;
		 return;
	     }
	     data_length = buffer_data_end;
	 } else {
	     ssize_t offset;
	     message_buffer = trans->read_block_func(&CMstatic_trans_svcs, 
						     conn->transport_data,
						     &data_length, &offset);
	     if (message_buffer == NULL) {
		 CMtrace_out(cm, CMLowLevelVerbose, 
			     "CMdata NULL return from read_block_func");
		 return;
	     }
	     message_buffer->ref_count++;
	     CMtrace_out(cm, CMBufferVerbose, "Received buffer %p from transport read_block_func, increment ref count, now is %d\n", message_buffer, message_buffer->ref_count);
	     tmp_message_buffer = (char *)message_buffer->buffer + offset;
	     buffer_data_end = data_length;
	     cm->abort_read_ahead = 1;

	     if (data_length == -1) {
		 CMtrace_out(cm, CMLowLevelVerbose, 
			     "CMdata read failed, actual %zu, failing connection %p\n", data_length, conn);
		 CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection failed with dereference, data length %p\n", conn);
		 INT_CMConnection_failed(conn);
		 return;
	     }
	     if (data_length == 0) {
		 return;
	     }
	     if (tmp_message_buffer == NULL) {
		 CMtrace_out(cm, CMLowLevelVerbose, "CMdata read_block failed, failing connection %p\n", conn);
		 CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection failed read_block withdereference %p\n", conn);
		 INT_CMConnection_failed(conn);
		 return;
	     }
	     CMtrace_out(cm, CMLowLevelVerbose, "CMdata read_block returned %zu bytes of data\n", data_length);
	 }
	 if (cm_postread_hook) {
	     cm_postread_hook(data_length, tmp_message_buffer);
	 }
     }
     result = CMact_on_data(conn, message_buffer, tmp_message_buffer, data_length);
     if (result != 0) {
	 /* whoops, reevaluate.  We need more data */
	 buffer_full_point += result;
	 goto start_read;
     }
     /* OK, we consumed all data */
     buffer_full_point = 0;
     buffer_data_end = 0;
     if (message_buffer) {
	 cm_return_data_buf(cm, message_buffer);
	 message_buffer = NULL;
     }
     read_msg_count++;
     read_byte_count += data_length;

     /* try read-ahead */
     if (cm->abort_read_ahead == 1) {
	 cm->abort_read_ahead = 0;
	 CMtrace_out(cm, CMDataVerbose, 
		     "CM - readahead not tried, aborted for condition signal\n");
	 return;
     }	
     if ((read_msg_count > read_ahead_msg_limit) || 
	 (read_byte_count > read_ahead_byte_limit)) {
	 CMtrace_out(cm, CMDataVerbose, 
		     "CM - readahead not tried, fairness, read %d msgs, %zd bytes\n",
		     read_msg_count, read_byte_count);
	 return;
     } else {
	 goto start_read;
     }
 }

 static void
 CMdo_handshake(CMConnection conn, int handshake_version, int byte_swap, char *base)
 {
     int do_send = 1;
     int remote_format_server_ID;
     int remote_CManager_ID;
     if (byte_swap) {
	 ((char*)&remote_format_server_ID)[0] = base[3];
	 ((char*)&remote_format_server_ID)[1] = base[2];
	 ((char*)&remote_format_server_ID)[2] = base[1];
	 ((char*)&remote_format_server_ID)[3] = base[0];
	 ((char*)&remote_CManager_ID)[0] = base[7];
	 ((char*)&remote_CManager_ID)[1] = base[6];
	 ((char*)&remote_CManager_ID)[2] = base[5];
	 ((char*)&remote_CManager_ID)[3] = base[4];
     } else {
	 remote_format_server_ID = ((int *) base)[0];
	 remote_CManager_ID = ((int *) base)[1];
     }

     CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - Received CONN handshake message\n");
     if ((remote_CManager_ID & 0x80000000) == 0x80000000) {
	 /* the other fellow already has our ID */
	 do_send = 0;
	 remote_CManager_ID ^= 0x80000000;  /* kill high bit */
     }
     if (conn->remote_format_server_ID != 0) {
	 if (conn->remote_format_server_ID != remote_format_server_ID) {
	     printf("Gaak.  Got a second handshake on connection 0x%p, with a different format server ID %x vs. %x\n",
		    conn, conn->remote_format_server_ID, remote_format_server_ID);
	 } else {
	     printf("Less Gaak.  Got a second handshake on connection 0x%p, remote id %x\n",
		    conn, conn->remote_format_server_ID);
	 }
     } else {
	 conn->remote_format_server_ID = remote_format_server_ID;
	 conn->remote_CManager_ID = remote_CManager_ID;
	 CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - CONN handshake condition %d\n", conn->handshake_condition);
	 if (conn->handshake_condition != -1) {
	     INT_CMCondition_signal(conn->cm, conn->handshake_condition);
	     conn->handshake_condition = -1;
	 }
     }
     if (do_send) {
	 CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - Sending CONN handshake message\n");
	 send_and_maybe_wait_for_handshake(conn->cm, conn);
     } else {
	 CMtrace_out(conn->cm, CMLowLevelVerbose, "CM - *NOT* Sending CONN handshake message\n");
     }
 }

 static size_t
 CMact_on_data(CMConnection conn, CMbuffer cm_buffer, char *buffer, size_t length)
 {
     char *base = buffer;
     char *check_sum_base = buffer;
     int byte_swap = 0;
     int get_attrs = 0;
     int skip = 0;
     int performance_msg = 0, event_msg = 0, evcontrol_msg = 0, handshake = 0;
     int performance_func = 0, handshake_version = 0;
     CMbuffer cm_decode_buf = NULL;
     attr_list attrs = NULL;
     int64_t data_length, decoded_length;
     int attr_length = 0, i;
     int header_len;
     int stone_id;
     char *decode_buffer = NULL, *data_buffer;
     FFSTypeHandle local_format, original_format;
     CManager cm = conn->cm;
     CMincoming_format_list cm_format = NULL;
     int message_key;
     unsigned char checksum;
     int short_length = 1;

     if (length < 4) {
	 return 4 - length;
     }
     message_key = 0x00ffff00 & *((int*)buffer);
     switch (message_key) {  /* assume 4-byte int */
     case 0x00444d00: /* \0DMC reversed byte order */
	 byte_swap = 1;
     case 0x004d4400:  /* CMD\0 */
	 break;
     case 0x00444d01: /* \1DMC reversed byte order long msg*/
	 byte_swap = 1;
     case 0x004d4401:  /* CMD\1 long msg*/
	 short_length = 0;
	 break;
     case 0x00414d00: /* \0AMC reversed byte order */
	 byte_swap = 1;
     case 0x004d4100:  /* CMA\0 */
	 get_attrs = 1;
	 break;
     case 0x00005645: /* \0CVE reversed byte order */
	 byte_swap = 1;
     case 0x45560000: /* EVC\0 */
	 evcontrol_msg = 1;
	 break;
     case 0x00504d00: /* \0PMC reversed byte order */
	 byte_swap = 1;
     case 0x004d5000:  /* CMP\0 */
	 performance_msg = 1;
	 short_length = 0;
	 break;
     case 0x004c4d00: /* \0LMC reversed byte order */
	 byte_swap = 1;
     case 0x004d4c00:  /* CML\0 */
	 event_msg = 1;
	 get_attrs = 1;
	 break;
     case 0x00484d00: /* \0HMC reversed byte order - handshake */
	 byte_swap = 1;
     case 0x004d4800:  /* CMH\0 - handshake */
	 handshake = 1;
	 get_attrs = 0;
	 break;
     default:
	 /*  non CM message */
	 /*  lookup registered message prefixes and try to find handler */
	 /*  otherwise give up */
       {
	   int ret;
	   CMbuffer local = NULL;
 #ifdef EVER_HAVE_HANDLERS_OTHER_THAN_PBIO
	   if (cm_buffer == NULL) {
	       local = fill_cmbuffer(cm, buffer, length);
	       buffer = local->buffer;
	   }
 #endif
	   CManager_unlock(cm);
	   switch (*(int*)buffer) {
	   case 0x5042494f:
	   case 0x4f494250:  /* incoming FFS format protocol message */
	     {
	       extern int CM_pbio_query(CMConnection conn, CMTransport trans,
					char *buffer, size_t length);
	       
	       size_t ret = CM_pbio_query(conn, conn->trans, buffer, length);
	       CManager_lock(cm);
	       return ret;
	     }
	   }
	   ret = CMdo_non_CM_handler(conn, *(int*)buffer, buffer, length);
	   CManager_lock(cm);
	   if (local) cm_return_data_buf(cm, local);
	   if (ret == -1) {
	       printf("Unknown message on connection %p, failed %d, closed %d, %x\n", conn, conn->failed, conn->closed, *(int*)buffer);
	       CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection unknown message failed with dereference %p\n", conn);
	       INT_CMConnection_failed(conn);
	   }
	   return 0;
       }
     }

     if (get_attrs == 1) {
	 if (!event_msg) {
	     header_len = 16;/* magic plus two 4-byte sizes (attrs + data) */
	     skip = 4;
	 } else {
	     header_len = 16;
	 }
	 if (!short_length) {
	     header_len += 4; /* extra data length bytes */
	 }
     } else {
	 if (short_length) {
	     header_len = 8; /* magic plus 4-byte size */
	 } else {
	     header_len = 12; /* magic plus 8-byte size */
	 }
     }

     if (length < header_len) {
	 return header_len - length;
     }
     base = buffer + 4 + skip; /* skip used data */
     if (short_length) {
	 if (byte_swap) {
	     int tmp;
 #ifdef WORDS_BIGENDIAN	    
	     checksum = (unsigned char) check_sum_base[0];
 #else
	     checksum = (unsigned char) check_sum_base[3];
 #endif
	     ((char*)&tmp)[0] = base[3];
	     ((char*)&tmp)[1] = base[2];
	     ((char*)&tmp)[2] = base[1];
	     ((char*)&tmp)[3] = base[0];
	     data_length = tmp;
	     if (header_len != 8) {
		 ((char*)&attr_length)[0] = base[7];
		 ((char*)&attr_length)[1] = base[6];
		 ((char*)&attr_length)[2] = base[5];
		 ((char*)&attr_length)[3] = base[4];
	     }
	 } else {
 #ifdef WORDS_BIGENDIAN	    
	     checksum = (unsigned char) check_sum_base[3];
 #else
	     checksum = (unsigned char) check_sum_base[0];
 #endif
	     data_length = ((int *) base)[0];
	     if (header_len != 8) {
		 attr_length = ((int *) base)[1];
	     }
	 }
     } else {
	 if (byte_swap) {
 #ifdef WORDS_BIGENDIAN	    
	     checksum = (unsigned char) check_sum_base[0];
 #else
	     checksum = (unsigned char) check_sum_base[3];
 #endif
	     int tmp;
	     ((char*)&tmp)[0] = base[3];
	     ((char*)&tmp)[1] = base[2];
	     ((char*)&tmp)[2] = base[1];
	     ((char*)&tmp)[3] = base[0];
	     data_length = ((int64_t)tmp) << 32;
	     ((char*)&tmp)[0] = base[7];
	     ((char*)&tmp)[1] = base[6];
	     ((char*)&tmp)[2] = base[5];
	     ((char*)&tmp)[3] = base[4];
	     if (header_len != 12) {
		 ((char*)&attr_length)[0] = base[11];
		 ((char*)&attr_length)[1] = base[10];
		 ((char*)&attr_length)[2] = base[9];
		 ((char*)&attr_length)[3] = base[8];
	     }
	 } else {
 #ifdef WORDS_BIGENDIAN	    
	     checksum = (unsigned char) check_sum_base[3];
 #else
	     checksum = (unsigned char) check_sum_base[0];
 #endif
	     data_length = ((int64_t)(((int *) base)[0])) << 32;
	     data_length += ((int *) base)[1];
	     if (header_len != 12) {
		 attr_length = ((int *) base)[1];
	     }
	 }
     }
     if (performance_msg || evcontrol_msg) {
	 performance_func = 0xff & (data_length >> 56);
	 data_length &= 0xffffffffffffff;
	 data_length -= 12;  /* subtract off header size */
     }
     if (handshake) {
	 handshake_version = 0xff & (data_length >> 24);
	 data_length &= 0xffffff;
	 data_length -= 8;  /* subtract off header size */
     }
     if (event_msg) {
	 char *stone_base = base;
	 if (header_len == 16) {
	     if (byte_swap) {
		 ((char*)&stone_id)[0] = stone_base[11];
		 ((char*)&stone_id)[1] = stone_base[10];
		 ((char*)&stone_id)[2] = stone_base[9];
		 ((char*)&stone_id)[3] = stone_base[8];
	     } else {
		 stone_id = ((int *) stone_base)[2];
	     }
	 } else if (header_len == 20) {
	     if (byte_swap) {
		 ((char*)&stone_id)[0] = stone_base[15];
		 ((char*)&stone_id)[1] = stone_base[14];
		 ((char*)&stone_id)[2] = stone_base[13];
		 ((char*)&stone_id)[3] = stone_base[12];
	     } else {
		 stone_id = ((int *) stone_base)[3];
	     }
	 }
     }

     if ((ssize_t)length < header_len + data_length + attr_length) {
	 return header_len + data_length + attr_length - 
	     length;
     }
     /* At this point, the message is accepted.  Determine processing */
     base = buffer + header_len;
     if (handshake) {
	 CMdo_handshake(conn, handshake_version, byte_swap, base);
	 return 0;
     }
     if (checksum != 0) {
	 unsigned char calculated_checksum = 0;
	 for (i=4; i < length; i++) {
	     calculated_checksum += ((unsigned char *)buffer)[i];
	 }
	 if (calculated_checksum != checksum) {
	     printf("Discarding incoming message because of corruption.  Checksum mismatch got %x, expected %x\n",
		    calculated_checksum, checksum);
	     printf("Message was : ");
	     for (i=0 ; i < length; i++) {
		 printf(" %02x",  ((unsigned char *)buffer)[i]);
	     }
	     printf("\n");
	     return 0;
	 }
     }
     if (performance_msg) {
	 CMdo_performance_response(conn, data_length, performance_func, byte_swap,
				   base);
	 return 0;
     } else if (evcontrol_msg) {
	 int arg;
	 if (byte_swap) {
	     ((char*)&arg)[0] = base[3];
	     ((char*)&arg)[1] = base[2];
	     ((char*)&arg)[2] = base[1];
	     ((char*)&arg)[3] = base[0];
	 } else {
	     arg = ((int *) base)[0];
	 }
 #ifdef EV_INTERNAL_H
	 INT_EVhandle_control_message(conn->cm, conn, (unsigned char) performance_func, arg);
 #endif
	 return 0;
     }
     data_buffer = base + attr_length;
     if (attr_length != 0) {
	 attrs = CMdecode_attr_from_xmit(conn->cm, base);
	 if (CMtrace_on(conn->cm, CMDataVerbose)) {
	     fprintf(cm->CMTrace_file, "CM - Incoming read attributes -> ");
	     fdump_attr_list(cm->CMTrace_file, attrs);
	 }
     }
     if (event_msg) {
	 CMbuffer local = NULL;
	 CMtrace_out(cm, CMDataVerbose, "CM - Receiving event message data len %"  PRId64 ", attr len %d, stone_id %x\n",
		     data_length, attr_length, stone_id);
	 if (attrs == NULL){
	     attrs = CMcreate_attr_list(cm);
	 }
	 set_int_attr(attrs, CM_EVENT_SIZE, (int)data_length);
	 set_long_attr(attrs, CM_INCOMING_CONNECTION, (intptr_t)conn);

	 if (cm_buffer == NULL) {
	     local = fill_cmbuffer(cm, buffer, length);
	     data_buffer = (data_buffer - buffer) + (char*)local->buffer;
	     buffer = local->buffer;
	     cm_buffer = local;
	 }
 #ifdef EV_INTERNAL_H	
	 internal_cm_network_submit(cm, cm_buffer, attrs, conn, data_buffer,
				    data_length, stone_id);
 #endif
	 if (local) cm_return_data_buf(cm, local);
	 free_attr_list(attrs);
	 return 0;
     }
     {
	 FMFormat format = FMformat_from_ID(FMContext_from_FFS(cm->FFScontext), data_buffer);
	 char *incoming_name;
	 if (format == NULL) {
	     printf("BAD INCOMING DATA\n");
	     return 0;
	 }
	 incoming_name = name_of_FMformat(format);

	 for (i = 0; i < cm->reg_format_count; i++) {
	     if ((cm->reg_formats[i]->registration_pending) && 
		 (strcmp(incoming_name, cm->reg_formats[i]->format_name) == 0)) {
		 CMcomplete_format_registration(cm->reg_formats[i], 0);
	     }
	 }
     }
     local_format = FFS_target_from_encode(conn->cm->FFScontext, data_buffer);
     original_format = FFSTypeHandle_from_encode(conn->cm->FFScontext, data_buffer);
     if (local_format == NULL) {
	 if (conn->cm->unregistered_format_handler) {
	     conn->cm->unregistered_format_handler(conn, name_of_FMformat(FMFormat_of_original(original_format)));
	 } else {
	     fprintf(stderr, "No conversion found for incoming CM message\n");
	 }
	 return 0;
     }
     CMtrace_out(cm, CMDataVerbose, "CM - Receiving record of type %s, FFSformat %p\n", 
		 name_of_FMformat(FMFormat_of_original(original_format)), original_format);
     for (i=0; i< cm->in_format_count; i++) {
	 if (cm->in_formats[i].format == local_format) {
	     cm_format = &cm->in_formats[i];
	     CMtrace_out(cm, CMDataVerbose, "CM - Found incoming cm_format %p, matching FFSformat %p\n", 
			 cm_format, local_format);
	 }
     }
     if (cm_format == NULL) {
	 cm_format = CMidentify_rollbackCMformat(cm, data_buffer);
	 CMtrace_out(cm, CMDataVerbose, "CM - Created cm_format %p, matching FFSformat %p\n", 
		     cm_format, original_format);
	 if(cm_format) {
	     CMtrace_out(cm, CMDataVerbose, "CM - Calling CMcreate_conversion type %s, format %p\n", 
			 name_of_FMformat(FMFormat_of_original(original_format)), original_format);
	     CMcreate_conversion(cm, cm_format);
	     CMtrace_out(cm, CMDataVerbose, "CM - after CMcreate_conversion format %p, has_conversion is %d\n", 
			 original_format, FFShas_conversion(original_format));
	 }
     }

     if ((cm_format == NULL) || (cm_format->handler == NULL)) {
	 fprintf(stderr, "CM - No handler for incoming data of this version of format \"%s\"\n",
		 name_of_FMformat(FMFormat_of_original(original_format)));
	 return 0;
     } else if (!FFShas_conversion(original_format)) {
	 CMcreate_conversion(cm, cm_format);
     }
     assert(FFShas_conversion(original_format));

     if (FFSdecode_in_place_possible(original_format)) {
	 if (!FFSdecode_in_place(cm->FFScontext, data_buffer, 
					(void**) (intptr_t) &decode_buffer)) {
	     printf("Decode failed\n");
	     return 0;
	 }
     } else {
	 decoded_length = FFS_est_decode_length(cm->FFScontext, data_buffer, data_length);
	 cm_decode_buf = cm_get_data_buf(cm, decoded_length);
	 decode_buffer = cm_decode_buf->buffer;
	 FFSdecode_to_buffer(cm->FFScontext, data_buffer, decode_buffer);
     }
     if(cm_format->older_format) {
 #ifdef EVOL
	 if(!process_old_format_data(cm, cm_format, &decode_buffer, &cm_decode_buf)){
	     return 0;
	 }
 #endif
     }
     if (CMtrace_on(conn->cm, CMDataVerbose)) {
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
	 fprintf(cm->CMTrace_file, "CM - record type %s, contents are:\n  ", name_of_FMformat(FMFormat_of_original(cm_format->format)));
	 r = FMfdump_data(cm->CMTrace_file, FMFormat_of_original(cm_format->format), decode_buffer, dump_char_limit);
	 if (!r && !warned) {
	     printf("\n\n  ****  Warning **** CM record dump truncated\n");
	     printf("  To change size limits, set CMDumpSize environment variable.\n");
	     warned++;
	 }
	 fprintf(cm->CMTrace_file, "\n=======\n");
     }
     if (attrs == NULL) {
	 attrs = CMcreate_attr_list(cm);
	 CMattr_merge_lists(cm, attrs, conn->attrs);
     } else {
	 CMattr_merge_lists(cm, attrs, conn->attrs);
     }

     /* 
      *  Handler may recurse, so clear these structures first
      */
     CMtrace_out(cm, CMFreeVerbose, "CM - add reference connection %p - handler\n", conn);
     INT_CMConnection_add_reference(conn);
     {
         CMHandlerFunc handler = cm_format->handler;
	 void *client_data = cm_format->client_data;
	 CMbuffer local = NULL;
	 if ((cm_buffer == NULL) && (decode_buffer == NULL)) {
	     local = fill_cmbuffer(cm, buffer, length);
	     buffer = local->buffer;
	 }
	 CManager_unlock(cm);
	 handler(cm, conn, decode_buffer, client_data, attrs);
	 CManager_lock(cm);
	 if (local) cm_return_data_buf(cm, local);
	 CMtrace_out(cm, CMFreeVerbose, "CM - delete reference connection %p - handler\n", conn);
	 INT_CMConnection_dereference(conn);
     }
     if (cm_decode_buf) {
	 cm_return_data_buf(cm, cm_decode_buf);
	 cm_decode_buf = NULL;
     }
     if (attrs) {
	 INT_CMfree_attr_list(cm, attrs);
	 attrs = NULL;
     }
     CMtrace_out(cm, CMDataVerbose, "CM - Finish processing - record of type %s\n", 
		 name_of_FMformat(FMFormat_of_original(original_format)));
     return 0;
 }

 void *
 INT_CMtake_buffer(CManager cm, void *data)
 {
     CMbuffer buf = cm_buffer_lookup(cm, data);
     if (buf == NULL) {
	 fprintf(stderr, "Error: INT_CMtake_buffer called with record %p not associated with cm\n", data);
	 cm_buffer_dump_list(cm);
	 return NULL;
     } else {
	 buf->ref_count++;
	 CMtrace_out(cm, CMBufferVerbose, "CMtake_buffer, data %p found buffer %p, ref_count incremented, now %d\n", data, buf, buf->ref_count);
	 return data;
     }
 }

 extern void
 INT_CMreturn_buffer(CManager cm, void *data)
 {
     CMbuffer buf = cm_buffer_lookup(cm, data);
     if (buf == NULL) {
	 fprintf(stderr, "Error: INT_CMreturn_buffer called with record %p not associated with cm\n", data);
	 cm_buffer_dump_list(cm);
	 return;
     }
     CMtrace_out(cm, CMBufferVerbose, "CMreturn_buffer, data %p found buffer %p, ref_count now %d, calling cm_return_data_buf\n", data, buf, buf->ref_count);
     cm_return_data_buf(cm, buf);
 }

void
INT_CMregister_handler(CMFormat format, CMHandlerFunc handler,
		       void *client_data)
{
    CManager cm = format->cm;
    int i;
    format->handler = handler;
    format->client_data = client_data;

    for (i=0; i< cm->in_format_count; i++) {
	if (strcmp(name_of_FMformat(FMFormat_of_original(cm->in_formats[i].format)), format->format_name) == 0) {
	    if (format->registration_pending) {
	        CMcomplete_format_registration(format, 1);
	    }
	    if (cm->in_formats[i].format == format->ffsformat) {
	        if (!cm->in_formats[i].handler) {
		    cm->in_formats[i].handler = handler;
		    cm->in_formats[i].client_data = client_data;
		} else if ((cm->in_formats[i].handler != handler) ||
			   (cm->in_formats[i].client_data != client_data)) {
		    fprintf(stderr, "Warning, CMregister_handler() called multiple times for the same format with different handler or client_data\n");
		    fprintf(stderr, "Repeated calls will be ignored\n");
		}
	    }
	}
    }
}

extern void
INT_CMregister_invalid_message_handler(CManager cm, CMUnregCMHandler handler)
{
    cm->unregistered_format_handler = handler;
}

 extern int
 INT_CMwrite(CMConnection conn, CMFormat format, void *data)
 {
     return INT_CMwrite_attr(conn, format, data, NULL);
 }

 extern void CMWriteQueuedData(transport_entry trans, CMConnection conn)
 {
     attr_list attrs = NULL;  /* GSE fix */
     CMtrace_out(conn->cm, CMLowLevelVerbose, "CMWriteQueuedData, conn %p, header %zd, attr %zd\n", 
		 conn, conn->queued_data.rem_header_len, 
		 conn->queued_data.rem_attr_len);
     if (conn->queued_data.rem_header_len != 0) {
	 struct FFSEncodeVec tmp_vec[1];
	 ssize_t actual;
	 tmp_vec[0].iov_base = conn->queued_data.rem_header;
	 tmp_vec[0].iov_len = conn->queued_data.rem_header_len;
	 actual = trans->NBwritev_func(&CMstatic_trans_svcs,
					    conn->transport_data,
					    &tmp_vec[0], 1,
					    attrs);
	 if (actual == -1) {
	     goto failed;
	 }
	 if (actual < (ssize_t)conn->queued_data.rem_header_len) {
	     conn->queued_data.rem_header_len -= actual;
	     memmove(&conn->queued_data.rem_header[0],
		     &conn->queued_data.rem_header[actual],
		     conn->queued_data.rem_header_len);
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "CMWriteQueuedData, conn %p, remaining header %zd\n", 
			 conn, conn->queued_data.rem_header_len);
	     return;
	 }
     }
     if (conn->queued_data.rem_attr_len != 0) {
	 struct FFSEncodeVec tmp_vec[1];
	 ssize_t actual;
	 tmp_vec[0].iov_base = conn->queued_data.rem_attr_base;
	 tmp_vec[0].iov_len = conn->queued_data.rem_attr_len;
	 actual = trans->NBwritev_func(&CMstatic_trans_svcs,
					    conn->transport_data,
					    &tmp_vec[0], 1,
					    attrs);
	 if (actual == -1) {
	     goto failed;
	 }
	 if (actual < (ssize_t)conn->queued_data.rem_attr_len) {
	     conn->queued_data.rem_attr_len -= actual;
	     conn->queued_data.rem_attr_base += actual;
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "CMWriteQueuedData, conn %p, remaining attr %zd\n", 
			 conn, conn->queued_data.rem_attr_len);
	     return;
	 }
     }
     if (conn->queued_data.vector_data) {
	 size_t vec_count = 0;
	 size_t length = 0;
	 FFSEncodeVector vec = conn->queued_data.vector_data;
	 ssize_t actual = 0;

	 while(vec[vec_count].iov_base != NULL) {
	     length += vec[vec_count].iov_len;
	     vec_count++;
	 }
	 actual = trans->NBwritev_func(&CMstatic_trans_svcs,
					    conn->transport_data,
					    vec, vec_count,
					    attrs);
	 if (actual == -1) {
	     goto failed;
	 }
	 if (actual < (ssize_t)length) {
	     int i = 0;
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "Continued partial pending write, %zu bytes sent\n", actual);
	     while (actual > (ssize_t)vec[i].iov_len) {
		 actual -= vec[i].iov_len;
		 i++;
		 vec_count--;
	     }
	     vec[i].iov_len -= actual;
	     vec[i].iov_base = (char*)vec[i].iov_base + actual;
	     conn->queued_data.vector_data = &vec[i];
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "CMWriteQueuedData, conn %p, %zu remaining data vectors\n", 
			 conn, vec_count);
	     return;
	 }
     }
     if (conn->queued_data.buffer_to_free) {
	 cm_return_data_buf(conn->cm, conn->queued_data.buffer_to_free);
     }
     conn->write_pending = 0;
     conn->trans->set_write_notify(conn->trans, &CMstatic_trans_svcs, 
				   conn->transport_data, 0);

     if(!CManager_locked(conn->cm)) {
	 printf("Not LOCKED in write queued data!\n");
     }
     cm_wake_any_pending_write(conn);
     return;
  failed:
     CMtrace_out(conn->cm, CMFreeVerbose, "Calling write failed connection failed with dereference %p\n", conn);
     INT_CMConnection_failed(conn);
     if (conn->queued_data.buffer_to_free) {
	 cm_return_data_buf(conn->cm, conn->queued_data.buffer_to_free);
	 conn->queued_data.buffer_to_free = NULL;
     }
     conn->write_pending = 0;
     conn->trans->set_write_notify(conn->trans, &CMstatic_trans_svcs, 
				   conn->transport_data, 0);
     cm_wake_any_pending_write(conn);
     return;
 }

 static void
 transport_wake_any_pending_write(CMConnection conn)
 {

     conn->write_pending = 0;
     CMtrace_out(conn->cm, CMTransportVerbose, "UNSet Pending write for conn %p\n", conn);
     cm_wake_any_pending_write(conn);
 }

 static void
 cm_wake_any_pending_write(CMConnection conn)
 {
     if (conn->write_callbacks) {
	 int i = 0;
	 CMConnHandlerListEntry callbacks[16];
	 size_t callback_len = conn->write_callback_len;
	 assert(conn->write_callback_len <= 16);
	 memcpy(callbacks, conn->write_callbacks, sizeof(callbacks[0]) * conn->write_callback_len);
	 for (i = 0; i < callback_len; ++i) {
	     if (callbacks[i].func) {
		 (callbacks[i].func)(conn->cm, conn, callbacks[i].client_data);
	     }
	 }
	 CMtrace_out(conn->cm, CMTransportVerbose, "Completed pending write, did %d notifications\n", i);
     } else {
	 CMtrace_out(conn->cm, CMTransportVerbose, "Completed pending write, No notifications\n");
     }
     CMwake_server_thread(conn->cm);
 }

 static void
 cm_set_pending_write(CMConnection conn)
 {
     assert(CManager_locked(conn->cm));
     CMtrace_out(conn->cm, CMTransportVerbose, "Set Pending write for conn %p\n", conn);
     conn->write_pending = 1;
 }

 static void
 queue_remaining_write(CMConnection conn, FFSEncodeVector tmp_vec, 
		       FFSEncodeVector pbio_vec, int vec_count,
		       attr_list attrs, size_t actual_bytes_written,
		       int attrs_present)
 {
     int i = 0, j = 0;
     size_t total_bytes = 0;
     size_t remaining_bytes = 0;

     for (i=0; i < vec_count; i++) {
	 total_bytes += tmp_vec[i].iov_len;
     }
     remaining_bytes = total_bytes - actual_bytes_written;
     i = 0;
     while (actual_bytes_written > tmp_vec[i].iov_len) {
	 actual_bytes_written -= tmp_vec[i].iov_len;
	 i++;
     }
     tmp_vec[i].iov_len -= actual_bytes_written;
     tmp_vec[i].iov_base = (char*)tmp_vec[i].iov_base + actual_bytes_written;
     if (pbio_vec) {
	 if (attrs_present) {
	     pbio_vec[i-2] = tmp_vec[i];
	 } else {
	     pbio_vec[i-1] = tmp_vec[i];
	 }
     }
     actual_bytes_written = 0;
 /*
  *    Data is either:
  *    vec 0 = header
  *    vec 1 = ffs
  *    vec 2 = ffs

  *    vec 0 = header
  *    vec 1 = attrs
  *    vec 2 = ffs
  *    vec 2 = ffs

  *    vec 0 = header
  *    vec 1 = notffs
  *    vec 2 = notffs

  *    vec 0 = header
  *    vec 1 = attrs
  *    vec 2 = notffs
  *    vec 2 = notffs
  */

     conn->queued_data.buffer_to_free = NULL;

     if (i == 0) {
	 /* didn't even write the 8 or 12 or 16 or 24 byte header */
	 assert(sizeof(conn->queued_data.rem_header) >= tmp_vec[0].iov_len);
	 memcpy(&conn->queued_data.rem_header, tmp_vec[0].iov_base, 
		tmp_vec[0].iov_len);
	 conn->queued_data.rem_header_len = tmp_vec[0].iov_len;
     } else {
	 conn->queued_data.rem_header_len = 0;
     }
     if (attrs_present && (i <= 1)) {
	 /* got stuck in encoded attributes */
	 conn->queued_data.rem_attr_base = tmp_vec[1].iov_base;
	 conn->queued_data.rem_attr_len = tmp_vec[1].iov_len;
	 i++;
     } else {
	 conn->queued_data.rem_attr_len = 0;
     }


     /* fixup pbio_vec */
     j = i - 1;  /* how far into the pbio vector are we? */
     if (attrs_present) j--;
     if (j == -1) {
	 /* no PBIO data */
	 conn->queued_data.vector_data = NULL;
	 return;
     }
     if (pbio_vec == NULL) {
	 /* DATA NOT FFS, don't do optimized copying */
	 /* Errr, something smarter here */
	 size_t data_length = remaining_bytes - conn->queued_data.rem_attr_len - conn->queued_data.rem_header_len;
	 size_t length = data_length + (sizeof(tmp_vec[0])*2);
	 char *ptr;
	 CMbuffer buf = cm_get_data_buf(conn->cm, length);
	 FFSEncodeVector vec = buf->buffer;
	 vec[0].iov_len = data_length;
	 vec[0].iov_base = ((char*)buf->buffer) + (sizeof(tmp_vec[0])*2);
	 vec[1].iov_len = 0;
	 vec[1].iov_base = NULL;
	 conn->queued_data.buffer_to_free = buf;
	 conn->queued_data.vector_data = vec;
	 if (i == 0) i=1;
	 if (attrs_present && i <= 1) i=2;
	 ptr = vec[0].iov_base;
	 while (i < vec_count) {
	     memcpy(ptr, tmp_vec[i].iov_base, tmp_vec[i].iov_len);
	     ptr += tmp_vec[i].iov_len;
	     data_length -= tmp_vec[i].iov_len;
	     i++;
	 }
	 return;
     } else {
	 conn->queued_data.buffer_to_free = NULL;
     }
     if (j >= 0) {
	 CMtrace_out(conn->cm, CMLowLevelVerbose, "Removing from pbio_vec at offset %d\n", (int) j);
	 pbio_vec[j].iov_len -= actual_bytes_written;
	 pbio_vec[j].iov_base = (char*)pbio_vec[j].iov_base + 
	     actual_bytes_written;
     } else {
	 j = 0;  /* nothing written */
     }

     /* 
      * copy application data (which had been left in place) into temporary
      * PBIO buffer as well.
      */
     conn->queued_data.vector_data = 
	 copy_all_to_FFSBuffer(conn->io_out_buffer, &pbio_vec[j]);
     tmp_vec = conn->queued_data.vector_data;
     i = 0; 
 }

 static void
 remove_pending_write_callback_by_id(CMConnection conn, SOCKET ids) {
     int id = (int)(intptr_t)ids;
     assert(id < conn->write_callback_len && id >= 0);
     conn->write_callbacks[id].func = NULL;
 }

 static void
 remove_pending_write_callback(CMConnection conn, CMWriteCallbackFunc handler,
			       void *client_data)
 {
     int i = 0;
     while (conn->write_callbacks[i].func != handler
	      && conn->write_callbacks[i].client_data != client_data) i++;
     conn->write_callbacks[i].func = NULL;
 }

 static int
 add_pending_write_callback(CMConnection conn, CMWriteCallbackFunc handler, 
			    void* client_data)
 {
     int count = 0;
     while (conn->write_callbacks && count < conn->write_callback_len &&
	    (conn->write_callbacks[count].func != NULL)) count++;
     if (count + 1 > conn->write_callback_len) {
	 if (conn->write_callbacks == NULL) {
	     conn->write_callbacks = malloc(sizeof(conn->write_callbacks[0]));
	     conn->write_callback_len = 1;
	 } else {
	     conn->write_callbacks = 
		 realloc(conn->write_callbacks,
			 sizeof(conn->write_callbacks[0])*(count+1));
	     conn->write_callback_len = count+1;
	 }
     }
     conn->write_callbacks[count].func = handler;
     conn->write_callbacks[count].client_data = client_data;
     return count;
 }


 static void
 wake_pending_write(CManager cm, CMConnection conn, void *param)
 {
     int cond = (int)(intptr_t)param;
     remove_pending_write_callback(conn, wake_pending_write, param);
     INT_CMCondition_signal(cm, cond);
 }

 void
 wait_for_pending_write(CMConnection conn)
 {
     CMControlList cl = conn->cm->control_list;
     assert(CManager_locked(conn->cm));
     CMtrace_out(conn->cm, CMLowLevelVerbose, "Wait for pending write for conn %p\n", conn);

     if ((!cl->has_thread) || (thr_thread_self() == cl->server_thread)) {
	 /* single thread working, just poll network */
	 while(conn->write_pending && !conn->closed) {
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "Control list wait for conn %p\n", conn);
	     CMcontrol_list_wait(cl);
	 }
     } else {
	 /* other thread is handling the network wait for it to wake us up */
	 while (conn->write_pending && !conn->closed) {
	     int cond = INT_CMCondition_get(conn->cm, conn);
	     add_pending_write_callback(conn, wake_pending_write, 
					(void*) (intptr_t)cond);
	     CMtrace_out(conn->cm, CMLowLevelVerbose, "Condition wait for conn %p\n", conn);
	     if (INT_CMCondition_wait(conn->cm, cond) == 0) {
		 /* condition wait failed, connection is dead */
		 conn->write_pending = 0;
	     }
	 }
     }	    
     CMtrace_out(conn->cm, CMLowLevelVerbose, "Done waiting for pending write for conn %p\n", conn);
 }

 void
 INT_CMConnection_wait_for_pending_write(CMConnection conn)
 {
     wait_for_pending_write(conn);
 }

 int
 INT_CMwrite_raw(CMConnection conn, FFSEncodeVector full_vec, FFSEncodeVector data_vec,
		 long vec_count, size_t byte_count, attr_list attrs, int data_vec_stack)
 {
     return INT_CMwrite_raw_notify(conn, full_vec, data_vec, vec_count, byte_count, attrs, data_vec_stack,
				   NULL, NULL);
 }

 /* Returns 1 if successful, -1 if deferred, 0 on error */
 int
 INT_CMwrite_raw_notify(CMConnection conn, FFSEncodeVector full_vec, FFSEncodeVector data_vec,
			long vec_count, size_t byte_count, attr_list attrs, int data_vec_stack,
			CMcompletion_notify_func notify_func, void *notify_client_data)
 {
     size_t actual = 0;
     unsigned char checksum = 0;
     int i, j, start;
     size_t count = 0;
     size_t length = 0;
     if (conn->closed || conn->failed) return 0;

     if (conn->write_pending) {
	 wait_for_pending_write(conn);
     }
     for (i=0; i < vec_count; i++) {
	 length += full_vec[i].iov_len;
     }
     start = 4;
     if (length < 10240) {
	 /* do checksum for small messages */
	 for (i=0; i < vec_count; i++) {
	     count += full_vec[i].iov_len - start;
	     for (j=start; j< full_vec[i].iov_len; j++) {
		 checksum += ((unsigned char*)full_vec[i].iov_base)[j];
	     }
	     start = 0;
	 }
     }
     ((int*)full_vec[0].iov_base)[0] = 
	 (((int*)full_vec[0].iov_base)[0] & 0xffffff00) | (unsigned char) checksum;
     if ((conn->do_non_blocking_write == 1) && (conn->trans->NBwritev_func)) {
	 ssize_t actual_bytes;
	 actual_bytes = 
	     conn->trans->NBwritev_func(&CMstatic_trans_svcs, 
					     conn->transport_data, 
					     full_vec, vec_count, attrs);
	 if (actual_bytes < 0) {
	     CMtrace_out(conn->cm, CMFreeVerbose, "Calling write failed connection failed with dereference %p\n", conn);
	     INT_CMConnection_failed(conn);
	     if (conn->queued_data.buffer_to_free) {
		 cm_return_data_buf(conn->cm, conn->queued_data.buffer_to_free);
		 conn->queued_data.buffer_to_free = NULL;
	     }
	     conn->write_pending = 0;
	     conn->trans->set_write_notify(conn->trans, &CMstatic_trans_svcs, 
					   conn->transport_data, 0);
	     cm_wake_any_pending_write(conn);
	 }
	 if (actual_bytes < (ssize_t)length) {
	     /* copy remaining and send it later */
	     if (actual_bytes < 0 ) actual_bytes = 0;
	     if (data_vec_stack) {
		 data_vec = copy_vector_to_FFSBuffer(conn->io_out_buffer, data_vec);
	     }
	     queue_remaining_write(conn, full_vec, data_vec, vec_count, 
				   attrs, actual_bytes, attrs != NULL);
	     conn->trans->set_write_notify(conn->trans, &CMstatic_trans_svcs, conn->transport_data, 1);
	     conn->write_pending = 1;
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Partial write, queued %zd bytes\n",
			 length - actual_bytes);
	     return 1;
	 }
	 actual = vec_count;  /* set actual for success */
     } else if (conn->trans->writev_complete_notify_func && notify_func) {
	 actual = conn->trans->writev_complete_notify_func(&CMstatic_trans_svcs, 
							   conn->transport_data, 
							   full_vec, vec_count, 
							   attrs, notify_func, notify_client_data);
     } else {
	 actual = conn->trans->writev_func(&CMstatic_trans_svcs, 
					   conn->transport_data, 
					   full_vec, vec_count, attrs);
	 if (actual <= 0) {
	     CMtrace_out(conn->cm, CMFreeVerbose, "Calling write failed connection failed with dereference %p\n", conn);
	     INT_CMConnection_failed(conn);
	 }
	 if (notify_func) {
	     (notify_func)(notify_client_data);
	 }
     }
     return actual == vec_count ? 1 : 0;
 }

 extern int
 INT_CMwrite_evcontrol(CMConnection conn, unsigned char type, int argument) {
     int evcontrol_header[2] = {0x45564300, 0};
     struct FFSEncodeVec static_vec[3];
     int success;
     FFSEncodeVector vec = &static_vec[0];
     assert(sizeof(int) == 4);
     vec[0].iov_base = evcontrol_header;
     vec[0].iov_len = sizeof(evcontrol_header);
     vec[1].iov_base = &argument; /* XXX int size */
     vec[1].iov_len = sizeof(int);
     vec[2].iov_base = NULL;
     vec[2].iov_len = 0;
     evcontrol_header[1] = type << 24 | (sizeof(evcontrol_header) + sizeof(int));
     success = INT_CMwrite_raw(conn, vec, vec + 1, 2, evcontrol_header[1] & 0xffffff, NULL, 1) != 0;
     return success;
 }

 extern int
 INT_CMwrite_attr(CMConnection conn, CMFormat format, void *data, 
		  attr_list attrs)
 {
     void *header_ptr = NULL;
     int header_len = 0;
     int no_attr_header[2] = {0x434d4400, 0};  /* CMD\0 in first entry */
//  not yet impl     int no_attr_long_header[4] = {0x434d4401, 0x434d4401, 0, 0};  /* CMD\1 in first entry, pad to 16 */
     int attr_header[4] = {0x434d4100, 0x434d4100, 0, 0};  /* CMA\0 in first entry, pad to 16 */
     int attr_long_header[4] = {0x434d4101, 0, 0, 0};  /* CMA\1 in first entry */
     FFSEncodeVector vec;
     size_t length = 0, vec_count = 0, actual;
     int do_write = 1;
     void *encoded_attrs = NULL;
     int attrs_present = 0;
     int long_message = 0;
     CManager cm = conn->cm;

     /* ensure conn is open */
     if (conn->closed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to closed connection\n");
	 return 0;
     }
     if (conn->failed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to failed connection\n");
	 return 0;
     }
     if (conn->write_pending) {
	 wait_for_pending_write(conn);
     }
     if (conn->closed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to closed connection\n");
	 return 0;
     }
     if (format->registration_pending) {
	 CMcomplete_format_registration(format, 1);
     }
     if (format->fmformat == NULL) {
	 printf("Format registration has failed for format \"%s\" - write aborted\n",
		format->format_name);
	 return 0;
     }
     if (conn->closed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to closed connection\n");
	 return 0;
     }
     CMformat_preload(conn, format);

     if (conn->closed != 0) {
	 return 0;
     }
     if (CMtrace_on(conn->cm, CMDataVerbose)) {
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
	 fprintf(cm->CMTrace_file, "CM - Writing record of type %s\n",
		name_of_FMformat(format->fmformat));
	 if (attrs != NULL) {
	     fprintf(cm->CMTrace_file, "CM - write attributes are:");
	     fdump_attr_list(cm->CMTrace_file, attrs);
	 }
	 fprintf(cm->CMTrace_file, "CM - record type %s, contents are:\n  ", name_of_FMformat(format->fmformat));
	 r = FMfdump_data(cm->CMTrace_file, format->fmformat, data, dump_char_limit);
	 if (!r && !warned) {
	     fprintf(cm->CMTrace_file, "\n\n  ****  Warning **** CM record dump truncated\n");
	     fprintf(cm->CMTrace_file, "  To change size limits, set CMDumpSize environment variable.\n");
	     warned++;
	 }
	 fprintf(cm->CMTrace_file, "\n=======\n");
     }

     /* encode data with CM context */
     vec = FFSencode_vector(conn->io_out_buffer, format->fmformat, data);
     while(vec[vec_count].iov_base != NULL) {
	 length += vec[vec_count].iov_len;
	 vec_count++;
     }
     if ((length & 0x7fffffff) == 0) {
	 long_message = 1;
     }
     if (attrs != NULL) {
	 attrs_present++;
     }
     if (!long_message) {
	 if (attrs_present) {
	     attr_header[2] = (int)length;
	     header_ptr = &attr_header;
	     header_len = sizeof(attr_header);
	 } else {
	     no_attr_header[1] = (int)length;
	     header_ptr = &no_attr_header;
	     header_len = sizeof(no_attr_header);
	 }
     } else {
	 if (attrs_present) {
	     memcpy((void*) &attr_long_header[1], &length, sizeof(length));
	     header_ptr = &attr_long_header;
	     header_len = sizeof(attr_long_header);
	 } else {
	     memcpy((void*) &attr_long_header[2], &length, sizeof(length));
	     header_ptr = no_attr_header;
	     header_len = sizeof(no_attr_header);
	 }
     }	 
     if (attrs_present) {
	 encoded_attrs = encode_attr_for_xmit(attrs, conn->attr_encode_buffer,
					      &attr_header[3]);
	 attr_header[3] = (attr_header[3] +7) & -8;  /* round up to even 8 */
	 attr_long_header[3] = (attr_header[3] +7) & -8;  /* round up to even 8 */
     }
     CMtrace_out(conn->cm, CMDataVerbose, "CM - Total write size is %zu bytes data + %d bytes attrs\n", length, attr_header[3]);
     if (cm_write_hook != NULL) {
	 do_write = cm_write_hook(length);
     }
     if (do_write) {
	 struct FFSEncodeVec static_vec[100];
	 FFSEncodeVector tmp_vec = &static_vec[0];
	 size_t byte_count = length;/* sum lengths */
	 if (vec_count >= sizeof(static_vec)/ sizeof(static_vec[0])) {
	     tmp_vec = INT_CMmalloc((vec_count+1) * sizeof(*tmp_vec));
	 }
	 tmp_vec[0].iov_base = header_ptr;
	 tmp_vec[0].iov_len = header_len;
	 if (attrs == NULL) {
	     memcpy(&tmp_vec[1], vec, sizeof(*tmp_vec) * vec_count);
	     vec_count++;
	     byte_count += sizeof(no_attr_header);
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writing %zu vectors, total %zu bytes in writev\n", 
			 vec_count, byte_count);
	 } else {
	     tmp_vec[1].iov_base = encoded_attrs;
	     tmp_vec[1].iov_len = attr_header[3];
	     memcpy(&tmp_vec[2], vec, sizeof(*tmp_vec) * vec_count);
	     byte_count += sizeof(attr_header) + attr_header[3];
	     vec_count += 2;
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writing %zu vectors, total %zu bytes (including attrs) in writev\n", 
			 vec_count, byte_count);
	 }

	 actual = INT_CMwrite_raw(conn, tmp_vec, vec, (int)vec_count, byte_count, attrs, 0);
	 if (tmp_vec != &static_vec[0]) {
	     INT_CMfree(tmp_vec);
	 }
	 if (actual == 0) {
	     /* fail */
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writev failed\n");
	     return 0;
	 }
     }
     CMtrace_out(conn->cm, CMLowLevelVerbose, "Writev success\n");
     return 1;
 }

 #ifdef EV_INTERNAL_H
 extern int
 internal_write_event(CMConnection conn, CMFormat format, void *remote_path_id,
		      int path_len, event_item *event, attr_list attrs, size_t *event_len_p)
 {
     FFSEncodeVector vec;
     struct FFSEncodeVec preencoded_vec[2];
     size_t data_length = 0, actual;
     long vec_count = 0;
     int attr_len = 0;
     int do_write = 1;
     void *encoded_attrs = NULL;
     int attrs_present = 0;
     CManager cm = conn->cm;

     /* ensure conn is open */
     if (conn->closed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to closed connection\n");
	 return 0;
     }
     if (conn->failed != 0) {
	 CMtrace_out(conn->cm, CMDataVerbose, "Not writing data to failed connection\n");
	 return 0;
     }
     if (conn->write_pending) {
	 wait_for_pending_write(conn);
     }
     if (format->registration_pending) {
	 CMcomplete_format_registration(format, 1);
     }
     if (format->fmformat == NULL) {
	 printf("Format registration has failed for format \"%s\" - write aborted\n",
		format->format_name);
	 return 0;
     }
     CMformat_preload(conn, format);

     if (CMtrace_on(conn->cm, CMDataVerbose)) {
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
	 fprintf(cm->CMTrace_file, "CM - Writing EVENT record %p of type %s\n", event,
		name_of_FMformat(format->fmformat));
	 if (attrs != NULL) {
	     fprintf(cm->CMTrace_file, "CM - write attributes are:");
	     fdump_attr_list(cm->CMTrace_file, attrs);
	 } else {
	     fprintf(cm->CMTrace_file, "CM - write attrs NULL\n");
	 }
	 fprintf(cm->CMTrace_file, "CM - record type %s, contents ", name_of_FMformat(format->fmformat));
	 if (event->decoded_event) {
	     fprintf(cm->CMTrace_file, "DECODED are:\n  ");
	     r = FMfdump_data(cm->CMTrace_file, format->fmformat, event->decoded_event,
			      dump_char_limit);
	 } else {
	     fprintf(cm->CMTrace_file, "ENCODED are:\n  ");
	     r = FMfdump_encoded_data(cm->CMTrace_file, format->fmformat,
				      event->encoded_event, dump_char_limit);
	 }	    
	 if (!r && !warned) {
	     fprintf(cm->CMTrace_file, "\n\n  ****  Warning **** CM record dump truncated\n");
	     fprintf(cm->CMTrace_file, "  To change size limits, set CMDumpSize environment variable.\n");
	     warned++;
	 }
	 fprintf(cm->CMTrace_file, "\n=======\n");
     }

     if (!event->encoded_event) {
	 /* encode data with CM context */
	 vec = FFSencode_vector(conn->io_out_buffer,
				format->fmformat, event->decoded_event);
	 while(vec[vec_count].iov_base != NULL) {
	     data_length += vec[vec_count].iov_len;
	     vec_count++;
	 }
     } else {
	 vec = &preencoded_vec[0];
	 preencoded_vec[0].iov_base = event->encoded_event;
	 preencoded_vec[0].iov_len = event->event_len;
	 preencoded_vec[1].iov_base = NULL;
	 preencoded_vec[1].iov_len = 0;
	 vec_count = 1;
	 data_length = event->event_len;
     }
     if (attrs != NULL) {
	 attrs_present++;
	 encoded_attrs = encode_attr_for_xmit(attrs, conn->attr_encode_buffer,
					      &attr_len);
	 attr_len = (attr_len +7) & -8;  /* round up to even 8 */
     }
     CMtrace_out(conn->cm, CMDataVerbose, "CM - Total write size is %zd bytes data + %d bytes attrs\n", data_length, attr_len);
     if (cm_write_hook != NULL) {
	 do_write = cm_write_hook(data_length);
     }
     if (do_write) {
	 struct FFSEncodeVec static_vec[100];
	 FFSEncodeVector tmp_vec = &static_vec[0];
	 size_t byte_count = data_length;/* sum lengths */
	 int header[4] = {0x434d4C00, 0, 0, 0};  /* CML\0 in first entry */
	 if (vec_count >= sizeof(static_vec)/ sizeof(static_vec[0])) {
	     tmp_vec = INT_CMmalloc((vec_count+3) * sizeof(*tmp_vec));
	 }
	 header[1] = (int)data_length;
	 if (path_len != 4) {
	     header[0] = 0x434d4700;
	     header[3] = (path_len + 7) & -8;
	 } else {
	     header[0] = 0x434d4C00;  /* 4 byte chan ID */
	     header[3] = *((int*)remote_path_id);
	 }
	 if (attrs == NULL) {
	     FFSEncodeVector assign_vec = tmp_vec;
	     header[2] = 0;
	     if (path_len != 4) {
		 tmp_vec[1].iov_base = remote_path_id;
		 tmp_vec[1].iov_len = (path_len + 7) & -8;
		 byte_count += tmp_vec[1].iov_len;
		 assign_vec++;
	     }
	     tmp_vec[0].iov_base = &header;
	     tmp_vec[0].iov_len = sizeof(header);
	     memcpy(&assign_vec[1], vec, sizeof(*tmp_vec) * vec_count);
	     vec_count++;
	     byte_count += sizeof(header);
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writing %lu vectors, total %zu bytes in writev\n", 
			 vec_count, byte_count);
	 } else {
	     tmp_vec[0].iov_base = &header;
	     tmp_vec[0].iov_len = sizeof(header);
	     tmp_vec[1].iov_base = encoded_attrs;
	     header[2] = attr_len;
	     tmp_vec[1].iov_len = attr_len;
	     memcpy(&tmp_vec[2], vec, sizeof(*tmp_vec) * vec_count);
	     byte_count += sizeof(header) + header[2];
	     vec_count += 2;
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writing %lu vectors, total %zu bytes (including attrs) in writev\n", 
			 vec_count, byte_count);
	 }
	 actual = INT_CMwrite_raw(conn, tmp_vec, vec, vec_count, byte_count, attrs,
				     vec == &preencoded_vec[0]);
	 if (tmp_vec != &static_vec[0]) {
	     INT_CMfree(tmp_vec);
	 }
	 if (actual <= 0) {
	     /* fail */
	     CMtrace_out(conn->cm, CMFreeVerbose, "Calling connection (write failed) failed with dereference %p\n", conn);
	     INT_CMConnection_failed(conn);
	     CMtrace_out(conn->cm, CMLowLevelVerbose, 
			 "Writev failed\n");
	     return 0;
	 }
     }
     if (event_len_p) *event_len_p = data_length;
     CMtrace_out(conn->cm, CMLowLevelVerbose, "Writev success\n");
     return 1;
 }
 #endif

 static void
 init_non_blocking_conn(CMConnection conn)
 {
     /* default */
     conn->do_non_blocking_write = 0;

     if (conn->trans->NBwritev_func == NULL) return;
     if (conn->trans->set_write_notify == NULL) return;

     /* only if we make it this far should we try non blocking writes */
     conn->do_non_blocking_write = 1;
 }

 extern int
 INT_CMConnection_write_would_block(CMConnection conn)
 {
     if (conn->do_non_blocking_write == -1) {
	 init_non_blocking_conn(conn);
     }
     return conn->write_pending;
 }

 extern int 
 INT_CMregister_write_callback(CMConnection conn, CMWriteCallbackFunc handler,
			       void *client_data)
 {
     if (conn->do_non_blocking_write == -1) {
	 init_non_blocking_conn(conn);
     }
     return add_pending_write_callback(conn, handler, client_data);
 }

 extern void
 INT_CMunregister_write_callback(CMConnection conn, SOCKET id)
 {
     remove_pending_write_callback_by_id(conn, id);
 }

 extern void
 INT_CM_fd_add_select(CManager cm, SOCKET fd, select_list_func handler_func,
		      void *param1, void *param2)
 {
     if (!handler_func) {
	 CMtrace_out(cm, EVWarning, "INT_CM_fd_add_select called with bogus notification function; ignored\n");
	 return;
     }
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     cm->control_list->add_select(&CMstatic_trans_svcs,
				  &cm->control_list->select_data, fd,
				  handler_func, param1, param2);
 }

 extern void
 CM_fd_write_select(CManager cm, SOCKET fd, select_list_func handler_func,
		    void *param1, void *param2)
 {
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     cm->control_list->write_select(&CMstatic_trans_svcs,
				    &cm->control_list->select_data, fd,
				    handler_func, param1, param2);
 }

 extern void
 CM_fd_remove_select(CManager cm, SOCKET fd)
 {
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     cm->control_list->remove_select(&CMstatic_trans_svcs,
				     &cm->control_list->select_data, fd);
 }

 extern CMTaskHandle
 INT_CMadd_periodic(CManager cm, long period, CMPollFunc func,
		    void *client_data)
 {
     CMTaskHandle handle = INT_CMmalloc(sizeof(*handle));
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     handle->cm = cm;
     handle->task = 
	 cm->control_list->add_periodic(&CMstatic_trans_svcs,
					&cm->control_list->select_data,
					0, period, (select_list_func)func, 
					(void*)cm, client_data);
     if (handle->task == NULL) {
	 free(handle);
	 return NULL;
     }
     return handle;
 }

 extern CMTaskHandle
 INT_CMadd_periodic_task(CManager cm, int period_sec, int period_usec,
			 CMPollFunc func, void *client_data)
 {
     CMTaskHandle handle = INT_CMmalloc(sizeof(*handle));
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     handle->cm = cm;
     handle->task = 
	 cm->control_list->add_periodic(&CMstatic_trans_svcs,
					&cm->control_list->select_data,
					period_sec, period_usec, 
					(select_list_func)func, 
					(void*)cm, client_data);
     if (handle->task == NULL) {
	 free(handle);
	 return NULL;
     }
     return handle;
 }

 extern void
 INT_CMremove_periodic(CMTaskHandle handle)
 {
     CManager cm = handle->cm;
     cm->control_list->remove_periodic(&CMstatic_trans_svcs,
				       &cm->control_list->select_data, 
				       handle->task);
     free(handle);
 }

 extern void
 INT_CMremove_task(CMTaskHandle handle)
 {
     CManager cm = handle->cm;
     cm->control_list->remove_periodic(&CMstatic_trans_svcs,
				       &cm->control_list->select_data, 
				       handle->task);
     free(handle);
 }

 extern CMTaskHandle
 INT_CMadd_delayed_task(CManager cm, int delay_sec, int delay_usec, 
			CMPollFunc func, void *client_data)
 {
     CMTaskHandle handle = INT_CMmalloc(sizeof(*handle));
     if (!cm->control_list->select_initialized) {
	 CM_init_select(cm->control_list, cm);
     }
     handle->cm = cm;
     handle->task = 
	 cm->control_list->add_delayed_task(&CMstatic_trans_svcs,
					    &cm->control_list->select_data,
					    delay_sec, delay_usec,
					    (select_list_func)func, 
					    (void*)cm, client_data);
     if (handle->task == NULL) {
	 free(handle);
	 return NULL;
     }
     return handle;
 }

 static void
 select_shutdown(CManager cm, void *shutdown_funcv)
 {
     SelectInitFunc shutdown_function = (SelectInitFunc)shutdown_funcv;
     CMtrace_out(cm, CMFreeVerbose, "calling select shutdown function sdp%p\n", cm->control_list->select_data);
     shutdown_function(&CMstatic_trans_svcs, cm, &cm->control_list->select_data);
 }

 static void
 select_free(CManager cm, void *task_datav)
 {
     void **task_data = (void**)task_datav;
     SelectInitFunc select_free_function = (SelectInitFunc)task_data[0];
     CMtrace_out(cm, CMFreeVerbose, "calling select FREE function, %p\n", task_data[1]);
     select_free_function(&CMstatic_trans_svcs, cm, &task_data[1]);
#if !NO_DYNAMIC_LINKING
     CMdlclose(task_data[2]);
#endif
     free(task_data);
 }

#ifdef HAVE_SYS_EPOLL_H
extern void
libcmepoll_init_sel_item(struct _select_item *sel_item);
#endif
extern void
libcmselect_init_sel_item(struct _select_item *sel_item);

static void
CM_init_select(CMControlList cl, CManager cm)
{
    SelectInitFunc init_function;
    SelectInitFunc shutdown_function;
    SelectInitFunc select_free_function;
    void *dlhandle = NULL;
    struct _select_item sel_item;
    char *select_module = cm->control_module_choice;
     
    CMtrace_out(cm, CMControlVerbose, "Loading CMselect module %s\n", select_module);
#if !NO_DYNAMIC_LINKING
    char *libname;
    lt_dlhandle handle;	
    lt_dladdsearchdir(EVPATH_MODULE_BUILD_DIR);
    lt_dladdsearchdir(EVPATH_MODULE_INSTALL_DIR);
    libname = malloc(strlen("lib" CM_LIBRARY_PREFIX "cm") + strlen(select_module) + strlen(MODULE_EXT) + 1);
#ifndef HAVE_WINDOWS_H
    strcpy(libname, "lib" CM_LIBRARY_PREFIX "cm");
#else
    strcpy(libname, CM_LIBRARY_PREFIX "cm");
#endif
    strcat(libname, select_module);
    strcat(libname, MODULE_EXT);
    handle = CMdlopen(cm->CMTrace_file, libname, 0);
    dlhandle = handle;
    free(libname);
    if (!handle) {
	fprintf(stderr, "Failed to load requested libcm%s dll.\n", select_module);
	fprintf(stderr, "Search path includes '.', '%s', '%s' and any default search paths supported by ld.so\n", EVPATH_MODULE_BUILD_DIR, 
		EVPATH_MODULE_INSTALL_DIR);
	fprintf(stderr, "Consider setting LD_LIBRARY_PATH or otherwise modifying module search paths.\n");
	exit(1);
    }
    sel_item.add_select = (CMAddSelectFunc)lt_dlsym(handle, "add_select");  
    sel_item.remove_select = (CMRemoveSelectFunc)lt_dlsym(handle, "remove_select");  
    sel_item.write_select = (CMAddSelectFunc)lt_dlsym(handle, "write_select");  
    sel_item.add_periodic = (CMAddPeriodicFunc)lt_dlsym(handle, "add_periodic");  
    sel_item.add_delayed_task = 
	(CMAddPeriodicFunc)lt_dlsym(handle, "add_delayed_task");  
    sel_item.remove_periodic = (CMRemovePeriodicFunc)lt_dlsym(handle, "remove_periodic");  
    sel_item.wake_function = (CMWakeSelectFunc)lt_dlsym(handle, "wake_function");
    sel_item.blocking_function = (CMPollFunc)lt_dlsym(handle, "blocking_function");
    sel_item.polling_function = (CMPollFunc)lt_dlsym(handle, "polling_function");;
    sel_item.initialize = (SelectInitFunc)lt_dlsym(handle, "select_initialize");
    sel_item.shutdown = (SelectInitFunc)lt_dlsym(handle, "select_shutdown");
    sel_item.free = (SelectInitFunc)lt_dlsym(handle, "select_free");
    sel_item.stop = (CMWakeSelectFunc)lt_dlsym(handle, "select_stop");
#else

#ifdef HAVE_SYS_EPOLL_H
    if (strcmp(select_module, "epoll") == 0) {
	libcmepoll_init_sel_item(&sel_item);
    }
#endif
    if (strcmp(select_module, "select") == 0) {
	libcmselect_init_sel_item(&sel_item);
    }

#endif
     cl->add_select = sel_item.add_select;
     cl->remove_select = sel_item.remove_select;
     cl->write_select = sel_item.write_select;
     cl->add_periodic = sel_item.add_periodic;
     cl->add_delayed_task = sel_item.add_delayed_task;
     cl->remove_periodic = sel_item.remove_periodic;
     cl->wake_select = sel_item.wake_function;
     cl->network_blocking_function.func = sel_item.blocking_function;
     cl->network_polling_function.func = sel_item.polling_function;
     init_function = sel_item.initialize;
     shutdown_function = sel_item.shutdown;
     select_free_function = sel_item.free;
     cl->stop_select = sel_item.stop;

    cl->network_blocking_function.client_data = (void*)&(cl->select_data);
    cl->network_blocking_function.cm = NULL;
    cl->network_polling_function.client_data = (void*)&(cl->select_data);
    cl->network_polling_function.cm = NULL;
     if ((cl->add_select == NULL) || (cl->remove_select == NULL) || 
	 (cl->network_blocking_function.func == NULL) || (cl->add_periodic == NULL) ||
	 (cl->remove_periodic == NULL)) {
	 printf("Select failed to load properly\n");
	 exit(1);
     }
     init_function(&CMstatic_trans_svcs, cm, &cm->control_list->select_data);
    if (cl->has_thread == -1) {
	thr_thread_t server_thread = 
	    thr_fork((void*(*)(void*))server_thread_func, 
		     (void*)cm);
	if (server_thread ==  (thr_thread_t) (intptr_t)NULL) {
	    return;
	}
	CMtrace_out(cm, CMLowLevelVerbose,
		    "CM - Forked comm thread %p\n", (void*)(intptr_t)server_thread);
	cm->control_list->server_thread = server_thread;
	cm->control_list->cl_reference_count++;
	cm->control_list->free_reference_count++;
	cl->has_thread = 1;
	cm->reference_count++;
	CMtrace_out(cm, CMFreeVerbose, "Forked - CManager %p ref count now %d\n", 
		    cm, cm->reference_count);
    }
     cl->select_initialized = 1;
     CMtrace_out(cm, CMFreeVerbose, "CManager adding select shutdown function, %p\n",shutdown_function);
     internal_add_shutdown_task(cm, select_shutdown, (void*)shutdown_function, SHUTDOWN_TASK);
     {
	 void ** data = malloc(3 * sizeof(void*));
	 data[0] = select_free_function;
	 data[1] = cm->control_list->select_data;
	 data[2] = dlhandle;
	 internal_add_shutdown_task(cm, select_free, (void*)data, FREE_TASK);
     }
 }

 static void
 wake_function(CManager cm, void *cond)
 {
     CManager_lock(cm);
     INT_CMCondition_signal(cm, (int)(intptr_t)cond);
     CManager_unlock(cm);
 }

 extern void
 INT_CMsleep(CManager cm, int sec)
 {
     int cond = INT_CMCondition_get(cm, NULL);
     CMTaskHandle handle = 
	 INT_CMadd_delayed_task(cm, sec, 0, wake_function, (void*)(intptr_t)cond);
     INT_CMfree(handle);
     INT_CMCondition_wait(cm, cond);
 }

 extern void
 INT_CMusleep(CManager cm, int usec)
 {
     int cond = INT_CMCondition_get(cm, NULL);
     CMTaskHandle handle = 
	 INT_CMadd_delayed_task(cm, 0, usec, wake_function, (void*)(intptr_t)cond);
     INT_CMfree(handle);
     INT_CMCondition_wait(cm, cond);
 }

typedef struct foreign_handler_struct {
    int header;
    CMNonCMHandler handler;
} *handler_list;

static handler_list foreign_handler_list;
static int foreign_handler_count = 0;

/* static void */
/* clear_foreign_handlers() */
/* { */
/*     if (foreign_handler_count == 0) return; */
/*     free(foreign_handler_list); */
/* } */

 extern void
 INT_CMregister_non_CM_message_handler(int header, CMNonCMHandler handler)
 {
     if (foreign_handler_count > 0) {
	 foreign_handler_list = INT_CMrealloc(foreign_handler_list, 
					  sizeof(foreign_handler_list[0]) * 
					  (foreign_handler_count + 1));
     } else {
	 foreign_handler_list = INT_CMmalloc(sizeof(foreign_handler_list[0]));
/*	 atexit(clear_foreign_handlers);*/
     }
     foreign_handler_list[foreign_handler_count].header = header;
     foreign_handler_list[foreign_handler_count].handler = handler;
     foreign_handler_count++;
 }

 int
 CMdo_non_CM_handler(CMConnection conn, int header, char *buffer, size_t length)
 {
     int i = 0;
     while (i < foreign_handler_count) {
	 if (foreign_handler_list[i].header == header) {
	     return foreign_handler_list[i].handler(conn, conn->trans, buffer, 
					     length);
	 }
	 i++;
     }
     return -1;
 }

 extern CMtrans_services
 INT_CMget_static_trans_services()
 {
   return &CMstatic_trans_svcs;
 }

 extern void*
 INT_CMget_transport_data (CMConnection conn)
 {
   return conn->transport_data;
 }

static 
int offset_compare(const void* lhsv, const void* rhsv)
 {
     CMavail_period_ptr lhs = (CMavail_period_ptr) lhsv;
     CMavail_period_ptr rhs = (CMavail_period_ptr) rhsv;
     if (lhs->offset.tv_sec < rhs->offset.tv_sec)
	 return -1;
     if (lhs->offset.tv_sec > rhs->offset.tv_sec)
	 return 1;
     return lhs->offset.tv_usec - rhs->offset.tv_usec;
 }

#ifdef _MSC_VER
static inline void timeradd(struct timeval *a, struct timeval *b,
							struct timeval *res)
{
	res->tv_sec = a->tv_sec + b->tv_sec;
	res->tv_usec = a->tv_usec + b->tv_usec;
	if (res->tv_usec >= 1000000)
	{
		res->tv_usec -= 1000000;
		res->tv_sec++;
	}
}
#endif

extern int
INT_CMinstall_pull_schedule(CManager cm, struct timeval *base_time, 
			    struct timeval *period, CMavail_period_ptr avail)
{
    int i = 0, count = 0;
    struct timeval zero = {0,0}, last_end = {0,0};
    CMavail_period_ptr sorted;
    while (timercmp(&avail[count].offset, &zero, !=) ||
	   timercmp(&avail[count].duration, &zero, !=)) {
	if (avail[count].offset.tv_sec < 0) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail sec offset is negative.  Rejected\n");
	    return 0;
	}
	if (avail[count].offset.tv_usec < 0) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail usec offset is negative.  Rejected\n");
	    return 0;
	}
	if (avail[count].duration.tv_sec < 0) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail sec duration is negative.  Rejected\n");
	    return 0;
	}
	if (avail[count].duration.tv_usec < 0) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail usec duration is negative.  Rejected\n");
	    return 0;
	}
	count++;
    }
    sorted = malloc(sizeof(avail[0]) * (count+1));
    memcpy(sorted, avail, sizeof(avail[0]) * count);
    qsort(sorted, count, sizeof(avail[0]), offset_compare);
    for (i = 0; i < count; i++) {
	struct timeval end;
	timeradd(&avail[i].offset, &avail[i].duration, &end);
	if (timercmp(&end, period, >)) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail region %d rejected, extends beyond period\n", i);
	    free(sorted);
	    return -1;
	}
	if (timercmp(&avail[i].offset, &last_end, <)) {
	    fprintf(stderr, "CMinstall_pull_schedule(), avail regions overlap. Rejected\n");
	    free(sorted);
	    return -1;
	}
	last_end = end;
    }

    cm->base_time = *base_time;
    cm->period = *period;
    cm->avail = sorted;
    transport_entry *trans_list;
    trans_list = cm->transports;
    CMtrace_out(cm, CMTransportVerbose, "CM installed pull schedule with period %ld secs, %zd usecs\n", period->tv_sec, (size_t) period->tv_usec);
    while ((trans_list != NULL) && (*trans_list != NULL)) {
	if ((*trans_list)->install_pull_schedule_func) {
	    (*trans_list)->install_pull_schedule_func(&CMstatic_trans_svcs,
						      *trans_list,
						      base_time, period,
						      cm->avail);
	    CMtrace_out(cm, CMTransportVerbose, "CM installed pull schedule to transport %s\n", (*trans_list)->trans_name);
	}
	trans_list++;
    }
    return 0;
}
