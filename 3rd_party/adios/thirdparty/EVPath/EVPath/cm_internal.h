#ifndef __I_O__
#include <ffs.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include "sys/time.h"
#endif
#ifndef _CM_SCHEDULE_H
#include "cm_schedule.h"
#endif
#ifndef _CM_CONFIG_H
#include "config.h"
#endif

#ifndef _MSC_VER
#include <pthread.h>
#define thr_mutex_t pthread_mutex_t
#define thr_thread_t pthread_t
#define thr_condition_t pthread_cond_t
#define thr_thread_self() pthread_self()
#define thr_thread_exit(status) pthread_exit(status);
#define thr_thread_detach(thread) pthread_detach(thread);
#define thr_thread_yield() sched_yield()
#define thr_thread_join(t, s) pthread_join(t, s)
#define thr_mutex_init(m) pthread_mutex_init(&m, NULL);
#define thr_mutex_lock(m) pthread_mutex_lock(&m);
#define thr_mutex_unlock(m) pthread_mutex_unlock(&m);
#define thr_mutex_free(m) pthread_mutex_destroy(&m);
#define thr_condition_init(c) pthread_cond_init(&c, NULL);
#define thr_condition_wait(c, m) pthread_cond_wait(&c, &m);
#define thr_condition_signal(c) pthread_cond_signal(&c);
#define thr_condition_broadcast(c) pthread_cond_broadcast(&c);
#define thr_condition_free(c) pthread_cond_destroy(&c);
#define thr_thread_create(w,x,y,z) pthread_create(w,x,y,z);
#else
//#include <mutex>
#include <Windows.h>
#define thr_mutex_t HANDLE
#define thr_thread_t DWORD
#define thr_condition_t HANDLE
#define thr_thread_create(w,x,y,z) 0
#define thr_thread_self() GetCurrentThreadId()
#define thr_thread_exit(status) 
#define thr_thread_detach(thread) 
#define thr_thread_yield() 
#define thr_thread_join(t, s) (void)s
#define thr_mutex_init(m) 
#define thr_mutex_lock(m)
#define thr_mutex_unlock(m)
#define thr_mutex_free(m)
#define thr_condition_init(c)
#define thr_condition_wait(c, m)
#define thr_condition_signal(c)
#define thr_condition_broadcast(c) 
#define thr_condition_free(c) 
#endif

#include <ev_internal.h>


#if defined (__INTEL_COMPILER)
/*  Allow extern declarations with no prior decl */
#  pragma warning (disable: 1418)
/*  Allow extern declarations in primary source files. */
#  pragma warning (disable: 1419)
/*  Don't tell me about unspecified order of evaluation */
#  pragma warning (disable: 981)
/*  Don't tell me about floating point equality tests */
#  pragma warning (disable: 1572)
/* Assert warning */
#  pragma warning (disable: 181)
#endif
#ifndef HAVE_COD_H
struct _ecl_code_struct;
typedef void *cod_parse_context;
#ifndef EV_INTERNAL_H
typedef struct extern_entry {
    /*! the textual name of the external entry */
    char *extern_name;
    /*! the address of the external entry */
    void *extern_value;
} cod_extern_entry;
typedef cod_extern_entry *cod_extern_list;
#endif
#endif


typedef struct _DelaySizeMtx {
/************
AveRTDelay[i] is the Average RTT to tranferring data of 
size MsgSize[i].  i=0,1,...,MsgNum. 
************/
    int MsgNum;			/* num of msgs */
    double *AveRTDelay;
    int *MsgSize;
}DelaySizeMtx;


typedef struct _CMincoming_format {
    FFSTypeHandle format;
    CMHandlerFunc handler;
    void *client_data;
    FMcompat_formats older_format;
    FMFormat local_prior_format;
    FMContext local_iocontext;
    CMFormat f2_format;
    int f1_struct_size;
    struct _cod_code_struct *code;
} *CMincoming_format_list;

struct _CMControlList;
typedef struct _CMControlList *CMControlList;

struct _pending_format_requests {
    char *server_id;
    int id_length;
    int condition;
    int top_request;
};

typedef struct func_entry {
    CMPollFunc func;
    CManager cm;
    void *client_data;
    int task_type;
} func_entry;

#include "cm_transport.h"

typedef struct pending_queue_entry {
    CMConnection conn;
    CMbuffer buffer;
    long length;
    struct pending_queue_entry *next;
} *pending_queue;

typedef struct _CManager {
    transport_entry *transports;
    int initialized;
    int reference_count;
    char *control_module_choice;  /* this is static, doesn't need to be free'd */

    CMControlList control_list;	/* the control list for this DE */

    int in_format_count;
    CMincoming_format_list in_formats;
    
    int reg_format_count;
    CMFormat *reg_formats;

    int reg_user_format_count;
    CMFormat *reg_user_formats;
  
    int pending_request_max;
    struct _pending_format_requests *pbio_requests;

    int connection_count;
    CMConnection *connections;

    thr_mutex_t exchange_lock;
    int locked;
    int closed;
    int abort_read_ahead;

    FFSContext FFScontext;	/* FFS context for data encoding */
    int FFSserver_identifier;	/* identifier for what FFS server we're talking to */
    thr_mutex_t context_lock;

    CMbuffer cm_buffer_list;
    pending_queue pending_data_queue;

    attr_list *contact_lists;

    func_entry *shutdown_functions;
    CMperf_upcall perf_upcall;
    CMUnregCMHandler unregistered_format_handler;

    struct _event_path_data *evp;
    FILE * CMTrace_file;

    /* pull schedule entries */
    struct timeval base_time;
    struct timeval period;
    CMavail_period_ptr avail;
} CManager_s;

typedef struct _CMCondition *CMCondition;

typedef void (*INT_CMfree_func)(void *block);

typedef enum _CMControlStyle {
    CMSingleThreaded, CMDedicatedServerThread, CMOccasionalPolling
} CMControlStyle;

typedef struct free_block_rec {
    int ref_count;
    CManager cm;
    void *block;
    INT_CMfree_func free_func;
    CManager locking_cm;
} *free_block_rec_p;

typedef void (*CMNetworkFunc)(void *svcs, void *client_data);


struct _CMTaskHandle {
    CManager cm;
    periodic_task_handle task;
};

typedef struct _CMControlList {
    func_entry network_blocking_function;
    func_entry network_polling_function;
    func_entry *polling_function_list;
    int pflist_size;
    int cl_consistency_number;

    int select_initialized;
    void *select_data;
    CMAddSelectFunc add_select;
    CMRemoveSelectFunc remove_select;
    CMAddSelectFunc write_select;
    CMAddPeriodicFunc add_periodic;
    CMAddPeriodicFunc add_delayed_task;
    CMRemovePeriodicFunc remove_periodic;
    CMWakeSelectFunc stop_select;
    CMWakeSelectFunc wake_select;
    /* 
     * CLs can be used by multiple DEs, close it when ref count reaches
     * zero 
     */
    int cl_reference_count;
    int free_reference_count;

    CMCondition condition_list;
    int next_condition_num;

    thr_mutex_t list_mutex;
    int locked;

    int closed;
    int has_thread;
    int cond_polling;
    thr_thread_t server_thread;
} CMControlList_s;

struct queued_data_rec {
    char rem_header[32];
    size_t rem_header_len;
    char *rem_attr_base;
    size_t rem_attr_len;
    FFSEncodeVector vector_data;
    CMbuffer buffer_to_free;
};

typedef struct _CMCloseHandlerList {
    CMCloseHandlerFunc close_handler;
    void *close_client_data;
    struct _CMCloseHandlerList *next;
} *CMCloseHandlerList;

typedef struct _CMConnHandlerList {
    CMCloseHandlerFunc func;
    void *client_data;
} *CMConnHandlerList, CMConnHandlerListEntry;

#define HEADER_BUFFER_SIZE 20
#ifndef CHR_TIME_H
/* to avoid including this everywhere it isn't needed, struct chr_time just needs to be big-ish.  Actually implementation can vary */
  typedef struct chr_time {
    double d1;
    double d2;
    double d3;
  } chr_time;
#endif

struct _CMConnection {
    CManager cm;
    /* remote contact info */

    transport_entry trans;
    void *transport_data;
    int conn_ref_count;
    FFSBuffer io_out_buffer;
    int closed;
    int failed;

    FMFormat *preloaded_formats;
    int remote_format_server_ID;   /* ID for the FFS format server in use by the peer */
    int remote_CManager_ID;   /* random, unique ID */
    int handshake_condition;

    CMCloseHandlerList close_list;

    size_t write_callback_len;
    CMConnHandlerList write_callbacks;
    AttrBuffer attr_encode_buffer;

    char header_buffer[HEADER_BUFFER_SIZE];		/* holds data until we know final size */
    CMbuffer message_buffer;    /* final destination of buffer */
    size_t buffer_full_point;	/* data required for buffer to be full */
    size_t buffer_data_end;	/* last point with valid data in buffer */

    attr_list characteristics;
    chr_time bandwidth_start_time;
    chr_time regressive_bandwidth_start_time; /*ztcai*/
    attr_list attrs;
    struct queued_data_rec queued_data;
    int write_pending;
    int do_non_blocking_write;
    int XML_output;
    int use_read_thread;
};

struct _CMFormat {
    CManager cm;
    char *format_name;
    FMFormat fmformat;
    FFSTypeHandle ffsformat;
    void *format_list_addr;
    CMHandlerFunc handler;
    void *client_data;
    FMStructDescList format_list;
    int registration_pending;
};

#define CManager_lock(cm) IntCManager_lock(cm, __FILE__, __LINE__)
#define CManager_unlock(cm) IntCManager_unlock(cm, __FILE__, __LINE__)
extern void IntCManager_lock(CManager cm, const char *file, int line);
extern void IntCManager_unlock(CManager cm, const char *file, int line);
extern int CManager_locked(CManager cm);
extern void CMControlList_lock(CMControlList cl);
extern void CMControlList_unlock(CMControlList cl);
extern int CMControlList_locked(CMControlList cl);

#define CMConn_write_lock(cm) IntCMConn_write_lock(cm, __FILE__, __LINE__)
#define CMConn_write_unlock(cm) IntCMConn_write_unlock(cm, __FILE__, __LINE__)
extern void IntCMConn_write_lock(CMConnection cl, char *file, 
				       int line);
extern void IntCMConn_write_unlock(CMConnection cl, char *file,
					 int line);
extern int CMConn_write_locked(CMConnection cl);

extern void 
CMtrace_out(CManager cm, CMTraceType trace_type, char *format, ...);

extern int
CMtrace_on(CManager cm, CMTraceType trace_type);

extern void 
CMDataAvailable(transport_entry trans, CMConnection conn);

extern void 
CMWriteQueuedData(transport_entry trans, CMConnection conn);

extern CMincoming_format_list
CMidentify_CMformat(CManager cm, FFSTypeHandle format);

extern void CMtransport_trace(CManager cm, const char *format, ...);
extern void CMtransport_verbose(CManager cm, CMTraceType trace, const char *format, ...);

extern void
CM_fd_add_select(CManager cm, SOCKET fd, select_list_func handler_func,
		       void *param1, void *param2);

extern void
CM_fd_write_select(CManager cm, SOCKET fd, select_list_func handler_func,
			 void *param1, void *param2);

extern void CM_fd_remove_select(CManager cm, SOCKET fd);

extern CMConnection
CMConnection_create(transport_entry trans, void *transport_data,
			  attr_list conn_attrs);

extern void free_CMFormat(CMFormat format);

extern void CMcomplete_format_registration(CMFormat format, int lock);
extern int CMcontrol_list_wait(CMControlList cl);
extern int load_transport(CManager cm, const char *trans_name, int quiet);
extern transport_entry add_transport_to_cm(CManager cm, transport_entry trans);

extern int CMinternal_listen(CManager cm, attr_list listen_info, int try_others);
extern CMConnection CMinternal_get_conn(CManager cm, attr_list attrs);
extern void CMconn_fail_conditions(CMConnection conn);
extern int CMpbio_send_format_preload(FMFormat ioformat, CMConnection conn);
extern void CMformat_preload(CMConnection conn, CMFormat format);
extern void CMinit_local_formats(CManager cm);

extern CMbuffer cm_get_data_buf(CManager cm, ssize_t length);
extern void cm_return_data_buf(CManager cm, CMbuffer cmb);
extern CMbuffer cm_create_transport_buffer(CManager cmb, void* buffer, ssize_t length);
extern CMbuffer cm_create_transport_and_link_buffer(CManager cmb, void* buffer, ssize_t length);

extern CMincoming_format_list CMidentify_rollbackCMformat 
(CManager cm, char *data_buffer);
extern void
CMcreate_conversion(CManager cm, CMincoming_format_list cm_format);
extern int
process_old_format_data(CManager cm, CMincoming_format_list cm_format,
	   	char **decode_buff, CMbuffer *cm_decode_buffer);
extern void
internal_add_shutdown_task(CManager cm, CMPollFunc func, void *client_data, int task_type);
extern void
internal_cm_network_submit(CManager cm, CMbuffer cm_data_buf, 
			   attr_list attrs, CMConnection conn, 
			   void *buffer, size_t length, int stone_id);
#define CMcreate_attr_list(cm) CMint_create_attr_list(cm, __FILE__, __LINE__)
#define INT_CMfree_attr_list(cm, l) CMint_free_attr_list(cm, l, __FILE__, __LINE__)
#define CMadd_ref_attr_list(cm, l) CMint_add_ref_attr_list(cm, l, __FILE__, __LINE__)
#define CMattr_copy_list(cm, l) CMint_attr_copy_list(cm, l, __FILE__, __LINE__)
#define CMattr_merge_lists(cm, l1, l2) CMint_attr_merge_lists(cm, l1, l2, __FILE__, __LINE__)
#define CMdecode_attr_from_xmit(cm, l) CMint_decode_attr_from_xmit(cm, l, __FILE__, __LINE__)

extern attr_list CMint_create_attr_list(CManager cm, char *file, int line);
extern void CMint_free_attr_list(CManager cm, attr_list l, char *file, int line);
extern attr_list CMint_add_ref_attr_list(CManager cm, attr_list l, char *file, int line);
extern attr_list CMint_attr_copy_list(CManager cm, attr_list l, char *file, int line);
extern void CMint_attr_merge_lists(CManager cm, attr_list l1, attr_list l2, 
					char *file, int line);
extern attr_list CMint_decode_attr_from_xmit(CManager cm, void * buf, char *file, int line);
extern void* INT_CMrealloc(void *ptr, size_t size);
extern void* INT_CMmalloc(size_t size);
#define malloc(x) INT_CMmalloc(x)
#define realloc(ptr, size) INT_CMrealloc(ptr, size)
extern void INT_CMfree(void *ptr);
extern void INT_CMadd_shutdown_task(CManager cm, CMPollFunc func, void *client_data, int task_type);
extern void INT_CManager_close(CManager cm);
extern CManager INT_CManager_create ();
extern CManager INT_CManager_create_control ( char *control_module);
extern int INT_CMlisten_specific(CManager cm, attr_list listen_info);
extern void INT_CMConnection_close(CMConnection conn);
extern void internal_connection_close(CMConnection conn);
extern void INT_CMremove_task(CMTaskHandle handle);
extern CMTaskHandle INT_CMadd_periodic(CManager cm, long period, 
					     CMPollFunc func, void *client_data);
extern CMTaskHandle
INT_CMadd_periodic_task(CManager cm, int period_sec, int period_usec, 
			  CMPollFunc func, void *client_data);
extern double
INT_CMregressive_probe_bandwidth(CMConnection conn, long size, attr_list attrs);
extern CMTaskHandle
INT_CMadd_delayed_task(CManager cm, int secs, int usecs, CMPollFunc func,
			     void *client_data);
extern int
INT_CMwrite_attr(CMConnection conn, CMFormat format, void *data, 
		       attr_list attrs);
extern int
INT_CMwrite_evcontrol(CMConnection conn, unsigned char type, int arg);
int INT_CMCondition_get(CManager cm, CMConnection dep);
void INT_CMCondition_signal(CManager cm, int condition);
void INT_CMCondition_set_client_data(CManager cm, int condition,
				       void *client_data);
void *INT_CMCondition_get_client_data(CManager cm, int condition);
int INT_CMCondition_wait(CManager cm, int condition);
extern void INT_CMCondition_fail(CManager cm, int condition);
extern attr_list INT_CMget_contact_list(CManager cm);
extern attr_list INT_CMderef_and_copy_list(CManager cm, attr_list attrs);
extern void INT_CMregister_non_CM_message_handler(int header, CMNonCMHandler handler);
extern void *INT_CMtake_buffer(CManager cm, void *data);
extern void INT_CMreturn_buffer(CManager cm, void *data);
extern CMConnection INT_CMget_conn(CManager cm, attr_list contact_list);
extern CMFormat INT_CMregister_format(CManager cm, FMStructDescList format_list);
extern CMFormat INT_CMregister_simple_format(CManager cm, char *format_name, FMFieldList field_list, int struct_size);
extern void
INT_EVforget_connection(CManager, CMConnection);
extern void
INT_EVhandle_control_message(CManager, CMConnection, unsigned char type, int arg);

extern void
INT_CMregister_handler(CMFormat format, CMHandlerFunc handler, 
			void *client_data);
extern long INT_CMprobe_latency(CMConnection conn, int msg_size,
				  attr_list attrs);
extern int
INT_CMwrite(CMConnection conn, CMFormat format, void *data);
extern CMConnection
INT_CMget_indexed_conn(CManager cm, int i);
extern int
INT_CMcontact_self_check(CManager cm, attr_list attrs);
extern int INT_CMtry_return_buffer(CManager cm, void *data);
extern FMFormat INT_CMget_IOformat_by_name(CManager cm, FMContext context,
					     char *name);
extern 
void INT_CMpoll_network(CManager cm);
extern 
void INT_CMrun_network(CManager cm);
extern void*
INT_CMget_transport_data(CMConnection conn);

extern int INT_CMCondition_has_failed(CManager cm, int condition);
extern int
INT_EVtake_event_buffer(CManager cm, void *event);
extern void
INT_EVPsubmit(CManager cm, int local_path_id, void *data, FMFormat format);
extern int INT_CMlisten(CManager cm);
extern char *
INT_create_filter_action_spec(FMStructDescList format_list, char *function);
extern char *
INT_create_bridge_action_spec(int stone_id, char *contact_string);
extern char *
INT_create_router_action_spec(FMStructDescList format_list, char *function);
extern int INT_CMfork_comm_thread(CManager cm);
extern int
INT_CMregister_write_callback(CMConnection conn, 
				CMWriteCallbackFunc handler,
				void *client_data);
extern void
INT_CMunregister_write_callback(CMConnection conn, SOCKET id);
extern void
INT_CMadd_poll(CManager cm, CMPollFunc func, void *client_data);
extern void
INT_EVPsubmit_encoded(CManager cm, int local_path_id, void *data, int len);
extern CMFormat INT_CMlookup_format(CManager cm, FMStructDescList format_list);
extern char *
INT_create_transform_action_spec(FMStructDescList format_list, FMStructDescList out_format_list, char *function);
extern char *
INT_create_multityped_action_spec(FMStructDescList *input_format_lists, char *function);

extern int INT_CMCondition_has_signaled(CManager cm, int condition);

extern attr_list
INT_CMget_specific_contact_list(CManager cm, attr_list attrs);

extern CMtrans_services
INT_CMget_static_trans_services();

extern void INT_CMsleep(CManager cm, int secs);
extern int INT_CMget_self_ip_addr(CManager cm);
extern void INT_CMget_port_range(CManager cm, int *high, int *low);
extern char *INT_CMget_ip_config_diagnostics(CManager cm);
extern void INT_CMget_qual_hostname(CManager cm, char *buf, int len);
extern attr_list INT_CMConnection_get_attrs(CMConnection conn);
extern void * INT_CMcreate_compat_info(CMFormat format, char *xform_code,
			int *len_p);
extern FMContext INT_CMget_user_type_context(CManager cm);
extern FFSTypeHandle INT_CMget_format_app_IOcontext(CManager cm, FFSContext context,
					     void *buffer, void *app_context);
extern FFSTypeHandle INT_CMget_format_IOcontext(CManager cm, FFSContext context,
					     void *buffer);
extern CMConnection
INT_CMinitiate_conn(CManager cm, attr_list contact_list);
extern void
INT_CMconn_register_close_handler(CMConnection conn, 
				    CMCloseHandlerFunc func, 
				    void *client_data);
extern void
INT_EVreturn_event_buffer(CManager cm, void *event);
extern void
INT_CMConnection_add_reference(CMConnection conn);
extern int
INT_CMConnection_set_character(CMConnection conn, attr_list attrs);
extern void
INT_CMremove_periodic(CMTaskHandle handle);
extern void INT_CMfree_user_type_context(CManager cm, FMContext context);
extern double
INT_CMprobe_bandwidth(CMConnection conn, long size, attr_list attrs);
extern int INT_CMConnection_write_would_block(CMConnection conn);
extern void INT_CMusleep(CManager cm, int usecs);
extern void INT_CM_insert_contact_info(CManager cm, attr_list attrs);
extern void INT_CM_fd_add_select(CManager cm, SOCKET fd, select_func handler_func, void *param1, void *param2);
extern void INT_CMstart_read_thread(CMConnection conn);
extern void INT_EVadd_standard_routines(CManager cm, char *extern_string, cod_extern_entry *externs);
extern void INT_EVadd_standard_structs(CManager cm, FMStructDescList *lists);
extern void INT_EVregister_close_handler(CManager cm, EVStoneCloseHandlerFunc handler, void *client_data );
extern void CMwake_server_thread(CManager cm);
extern int CMtrace_val[];
extern int CMtrace_timing;
extern int CMtrace_PID;
extern int CMtrace_init(CManager cm, CMTraceType t);
extern void INT_CMTrace_file_id(int ID);
#define CMtrace_on(cm, trace_type)  ((cm->CMTrace_file == NULL) ? CMtrace_init(cm, trace_type) : CMtrace_val[trace_type])
#ifdef _MSC_VER
#include <process.h>
#include <time.h>
#define CLOCK_MONOTONIC 1
static int clock_gettime(int cl, struct timespec* spec)
{
    __int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
    wintime -= 116444736000000000i64;  //1jan1601 to 1jan1970
    spec->tv_sec = (long)(wintime / 10000000i64);           //seconds
    spec->tv_nsec = wintime % 10000000i64 * 100;      //nano-seconds
    return 0;
}
#define HAVE_CLOCK_GETTIME
#endif
#ifdef HAVE_CLOCK_GETTIME
#define TRACE_TIME_DECL struct timespec ts
#define TRACE_TIME_GET clock_gettime(CLOCK_MONOTONIC, &ts)
#define TRACE_TIME_PRINTDETAILS "%lld.%.9ld - ", (long long)ts.tv_sec, ts.tv_nsec
#else
#define TRACE_TIME_DECL	struct timeval tv
#define TRACE_TIME_GET gettimeofday(&tv, NULL)
#define TRACE_TIME_PRINTDETAILS "%lld.%.6ld - ", (long long)tv.tv_sec, (long)tv.tv_usec
#endif

#define CMtrace_out(cm, trace_type, ...) {TRACE_TIME_DECL ; (CMtrace_on(cm,trace_type) ? (CMtrace_PID ? fprintf(cm->CMTrace_file, "P%lxT%lx - ", (long) getpid(), (long)thr_thread_self()) : 0) , CMtrace_timing? TRACE_TIME_GET,fprintf(cm->CMTrace_file, TRACE_TIME_PRINTDETAILS):0, fprintf(cm->CMTrace_file, __VA_ARGS__) : 0);fflush(cm->CMTrace_file);}
extern void CMdo_performance_response(CMConnection conn, size_t length, int func,
				      int byte_swap, char *buffer);
extern int
INT_CMwrite_raw(CMConnection conn, FFSEncodeVector full_vec, FFSEncodeVector data_vec,
                long vec_count, size_t byte_count, attr_list attrs, int data_vec_stack);
int
INT_CMwrite_raw_notify(CMConnection conn, FFSEncodeVector full_vec, FFSEncodeVector data_vec,
		       long vec_count, size_t byte_count, attr_list attrs, int data_vec_stack,
		       CMcompletion_notify_func notify_func, void *notify_client_data);
extern void
INT_CMConnection_dereference(CMConnection conn);
extern void INT_CMConnection_failed (CMConnection conn);

extern FMContext INT_CMget_FMcontext(CManager cm);
extern void INT_CMinstall_perf_upcall(CManager cm, CMperf_upcall upcall);
extern attr_list INT_CMtest_transport(CMConnection conn, attr_list how);
extern void INT_CMConnection_wait_for_pending_write(CMConnection conn);
extern EVstone INT_EVexecuting_stone(CManager cm);
extern void wait_for_pending_write(CMConnection conn);
extern int
INT_CMinstall_pull_schedule(CManager cm, struct timeval *base_time, 
			    struct timeval *period, CMavail_period_ptr avail);
extern void
INT_CMregister_invalid_message_handler(CManager cm, CMUnregCMHandler handler);

