#ifndef __CM_TRANSPORT_H__
#define __CM_TRANSPORT_H__

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif
#ifdef HAVE_SYS_TIME_H
#include "sys/time.h"
#endif
#ifndef _CM_SCHEDULE_H
#include "cm_schedule.h"
#endif

#include <stddef.h>
#include "sys/types.h"
#ifdef _MSC_VER
#include <BaseTsd.h>
    typedef SSIZE_T ssize_t;
#else
#include <unistd.h>
#ifndef SOCKET
#define SOCKET int
#endif
#endif

typedef struct _transport_item *transport_entry;
typedef struct _transport_item *CMTransport;

typedef struct _CMbuffer {
    void *buffer;
    size_t size;
    int ref_count;
    struct _CMbuffer *next;
    void (*return_callback)(void *);
    void *return_callback_data;
} *CMbuffer;

typedef enum _CMTraceType {
    CMAlwaysTrace, CMControlVerbose, CMConnectionVerbose, CMLowLevelVerbose, CMDataVerbose, CMTransportVerbose, CMFormatVerbose, CMFreeVerbose, CMAttrVerbose, CMBufferVerbose, EVerbose, EVWarning, CMSelectVerbose, EVdfgVerbose, 
    CMLastTraceType /* add before this one */
} CMTraceType;

typedef void *(*CMTransport_malloc_func)(size_t);
typedef void *(*CMTransport_realloc_func)(void*, size_t);
typedef void (*CMTransport_free_func)(void*);
typedef void (*CMTransport_wake_comm_thread_func)(CManager cm);
typedef void (*CMTransport_condition_signal_func)(CManager cm, int condition);

typedef void (*select_list_func)(void *, void*);

typedef void (*CMAddSelectFunc)(void *svcs, void *select_data, SOCKET fd,
				select_list_func func,
				void *param1, void *param2);

typedef void (*CMTransport_fd_add_select)(CManager cm, SOCKET fd, select_list_func handler_func,
					  void *param1, void *param2);
typedef void (*CMTransport_fd_remove_select)(CManager cm, SOCKET fd);
typedef void (*CMTransport_trace)(CManager cm, const char *format, ...);
typedef void (*CMTransport_verbose)(CManager cm, CMTraceType trace, const char *format, ...);
typedef CMConnection (*CMTransport_conn_create)(transport_entry trans,
						void *transport_data,
						attr_list conn_attrs);
typedef void (*CMTransport_add_shut_task)(CManager cm, CMPollFunc func, 
					  void *client_data, int task_type);
typedef CMTaskHandle (*CMTransport_add_period_task)(CManager cm, 
						    int period_sec,
						    int period_usec,
						    CMPollFunc func,
						    void *client_data);
typedef void (*CMTransport_remove_periodic)(CMTaskHandle cmt);
typedef void (*CMTransport_add_poll)(CManager cm, CMPollFunc func,
				     void *client_data);
typedef CMbuffer (*CMTransport_get_data_buffer)(CManager cm, ssize_t length);
typedef void (*CMTransport_return_data_buffer)(CManager cm, CMbuffer cmb);
typedef void (*CMTransport_connection_close)(CMConnection conn);
typedef void *(*CMTransport_get_transport_data)(CMConnection conn);
typedef void (*CMTransport_action_pending_write)(CMConnection conn);
typedef CMbuffer (*CMTransport_create_data_buffer)(CManager cm, void *buffer, ssize_t length);
typedef int (*CMTransport_modify_global_lock)(CManager cm, const char *file, int line);
typedef void (*CMTransport_add_buffer_to_pending_queue)(CManager cm, CMConnection conn, CMbuffer buf, long length);
typedef void (*CMTransport_cond_wait_CM_lock)(CManager cm, void *cond, char *file, int line);
typedef void (*CMRemoveSelectFunc)(void *svcs, void *select_data, SOCKET fd);
typedef struct _periodic_task *periodic_task_handle;

typedef periodic_task_handle (*CMAddPeriodicFunc) 
   (void *svcs, void *select_data, int period_sec, int period_usec,
	  select_list_func func, void *param1, void *param2);

typedef void (*CMRemovePeriodicFunc)(void *svcs, void *select_data, 
					   periodic_task_handle handle);

typedef void (*CMWakeSelectFunc)(void *svcs, void *select_data);

typedef struct CMtrans_services_s {
    CMTransport_malloc_func malloc_func;
    CMTransport_realloc_func realloc_func;
    CMTransport_free_func free_func;
    CMTransport_fd_add_select fd_add_select;
    CMTransport_fd_add_select fd_write_select;
    CMTransport_fd_remove_select fd_remove_select;
    CMTransport_trace trace_out;
    CMTransport_verbose verbose;
    CMTransport_conn_create connection_create;
    CMTransport_add_shut_task add_shutdown_task;
    CMTransport_add_period_task add_periodic_task;
    CMTransport_remove_periodic remove_periodic;
    CMTransport_add_poll add_poll;
    CMTransport_get_data_buffer get_data_buffer;
    CMTransport_return_data_buffer return_data_buffer;
    CMTransport_connection_close connection_close;
    CMTransport_create_data_buffer create_data_buffer;
    CMTransport_create_data_buffer create_data_and_link_buffer;
    CMTransport_get_transport_data get_transport_data;
    CMTransport_action_pending_write set_pending_write;
    CMTransport_action_pending_write wake_any_pending_write;
    CMTransport_modify_global_lock drop_CM_lock;
    CMTransport_modify_global_lock acquire_CM_lock;
    CMTransport_modify_global_lock return_CM_lock_status;
    CMTransport_cond_wait_CM_lock cond_wait_CM_lock;
    CMTransport_add_buffer_to_pending_queue add_buffer_to_pending_queue;
    CMTransport_connection_close connection_deref;
    CMTransport_connection_close connection_addref;
    CMTransport_connection_close connection_fail;
    CMTransport_wake_comm_thread_func wake_comm_thread;
    CMTransport_condition_signal_func condition_signal;
} *CMtrans_services;
#define DROP_CM_LOCK(svc, cm) (svc)->drop_CM_lock((cm), __FILE__, __LINE__)
#define ACQUIRE_CM_LOCK(svc, cm) (svc)->acquire_CM_lock((cm), __FILE__, __LINE__)
#define CM_LOCKED(svc, cm) (svc)->return_CM_lock_status((cm), __FILE__, __LINE__)

typedef void *(*CMTransport_func)(CManager cm, CMtrans_services svc, transport_entry trans);
typedef attr_list (*CMTransport_listen_func)(CManager cm,
					     CMtrans_services svc,
					     transport_entry trans,
					     attr_list listen_info);
typedef void *(*CMTransport_read_block_func)(CMtrans_services svc,
					     void *conn_data,
					     ssize_t *actual, ssize_t *offset);
typedef int (*CMTransport_read_to_buffer_func)(CMtrans_services svc,
					       void *conn_data,
					       void *buffer,
					       ssize_t len, int block_flag);
typedef int (*CMTransport_writev_func)(CMtrans_services svc,
				       void *transport_data,
				       void *buffer, ssize_t len,
				       attr_list attrs);

typedef void (*CMcompletion_notify_func)(void *client_data);
typedef int (*CMTransport_writev_complete_notify_func)(CMtrans_services svc,
						       void *transport_data,
						       void *buffer, ssize_t len,
						       attr_list attrs, CMcompletion_notify_func func,
						       void *client_data);
typedef void (*CMTransport_shutdown_conn_func)(CMtrans_services svc,
					       void *conn_data);

typedef CMConnection (*CMTransport_conn_func)(CManager cm, 
					      CMtrans_services svc,
					      transport_entry trans, 
					      attr_list attrs);

typedef CMConnection (*CMTransport_NBconn_func)(CManager cm, 
					        CMtrans_services svc,
					        transport_entry trans, 
                                                attr_list attrs,
                                                int condition);

typedef CMConnection (*CMTransport_NBconn_final_func)(CManager cm, 
                                                      CMtrans_services svc,
                                                      transport_entry trans, 
                                                      void *client_data,
                                                      int result);

typedef int (*CMTransport_self_check_func)(CManager cm,
					   CMtrans_services svc,
					   transport_entry trans,
					   attr_list attrs);

typedef int (*CMTransport_connection_eq_func)(CManager cm,
					      CMtrans_services svc,
					      transport_entry trans,
					      attr_list attrs,
					      void *conn_data);

typedef int (*CMTransport_set_write_notify_func) (transport_entry trans, 
						  CMtrans_services svc, 
						  void *conn_data, int enable);

typedef attr_list (*CMTransport_get_transport_characteristics) (transport_entry trans, 
								CMtrans_services svc, 
								void *conn_data);
typedef int (*CMTransport_install_pull_schedule)(CMtrans_services svc,
						 void *transport_data,
						 struct timeval *base_time, 
						 struct timeval *period, 
						 CMavail_period_ptr avail);

typedef void (*DataAvailableCallback)(transport_entry trans, CMConnection conn);
typedef void (*WritePossibleCallback)(transport_entry trans, CMConnection conn);
typedef void (*SelectInitFunc)(CMtrans_services svc, CManager cm, void *client_data);

struct _transport_item {
    char *trans_name;
    CManager cm;
    void *dlhandle;
    DataAvailableCallback data_available;
    WritePossibleCallback write_possible;
    CMTransport_func  transport_init;
    CMTransport_listen_func  listen;
    CMTransport_conn_func  initiate_conn;
    CMTransport_NBconn_func  initiate_conn_nonblocking;
    CMTransport_NBconn_final_func  finalize_conn_nonblocking;
    CMTransport_self_check_func  self_check;
    CMTransport_connection_eq_func  connection_eq;
    CMTransport_shutdown_conn_func  shutdown_conn;
    CMTransport_read_to_buffer_func read_to_buffer_func;
    CMTransport_read_block_func read_block_func;
    CMTransport_writev_func writev_func;
    CMTransport_writev_func NBwritev_func; /* non blocking */
    CMTransport_writev_complete_notify_func writev_complete_notify_func;
    CMTransport_set_write_notify_func set_write_notify;
    void *trans_data;
    CMTransport_get_transport_characteristics get_transport_characteristics;
    CMTransport_install_pull_schedule install_pull_schedule_func;
};

struct _select_item {
    CMAddSelectFunc add_select;
    CMRemoveSelectFunc remove_select;
    CMAddSelectFunc write_select;
    CMAddPeriodicFunc add_periodic;
    CMAddPeriodicFunc add_delayed_task;
    CMRemovePeriodicFunc remove_periodic;
    CMWakeSelectFunc wake_function;
    CMPollFunc blocking_function;
    CMPollFunc polling_function;
    SelectInitFunc initialize;
    SelectInitFunc shutdown;
    SelectInitFunc free;
    CMWakeSelectFunc stop;
};

extern void
get_IP_config(char *hostname_buf, int len, int* IP_p, int *port_range_low_p, int *port_range_high_p, 
	      int *use_hostname_p, attr_list attrs, CMTransport_trace trace_func, void *trace_data);

extern char *
IP_get_diagnostics(CManager cm, CMTransport_trace trace_out);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif


#endif
/*
 *  Documentation on Transport interfaces.

 *  There are two sets of interfaces, the set exported by the transport and 
 *  the set of upcalls provided by CM for the transport to use.  (Upcalls are 
 *  necessarily function pointers because DLLs can't resolve symbols
 *  in the main program on many platforms).
 *
 *  Calls that can or should be exported by the transport:
 *   (in the transport source, all calls must have names of the form
 *   "libcm<transport>_LTX_<routine>" where <transport> is the name of
 *   the transport and <routine> is the name of the subroutine.) 
 *
 *  - void *initialize(CManager cm, CMtrans_services svc, transport_entry trans);
 *      The initialize routine will be called once only when the
 *      transport is loaded.  It is passed in the cm value in use, a
 *      pointer to the set of upcalls (CMtrans_services) and a link to
 *      its own entry in the transport list.  The return value here is
 *      a void* that should be a pointer to a structure of
 *      transport-private data.  CM will not examine this pointer, but
 *      will pass it back to later routines as "transport_data".
 *      There are currently no provisions for unloading transports, so
 *      there is no provision to free() the transport-private data (aside
 *      from that done by the transport with 'add_shutdown_task').
 *  - attr_list non_blocking_listen(CManager cm, CMtrans_services svc,
 *                             transport_entry trans, attr_list listen_info);
 *      This routine will be called in response to a CMlisten() or
 *      CMlisten_specific() invocation.  It should cause the transport
 *      to listen for future incoming connections (without blocking the
 *      current thread).  The listen_info list, if non-NULL, contains
 *      attributes which should specify what 'address' or 'port' (or
 *      other transport-specific binding) should be listened on.  If
 *      unspecified, the listen binding can be arbitrary.  The return
 *      attr_list should contain sufficient information to initiate a
 *      connection to this process when passed to the same transport
 *      running in a different process/host.  The routine should also
 *      perform whatever tasks are necessary to service the listen
 *      'port', to respond appropriate to connection requests and to
 *      establish the necessary CMConnections when connections are
 *      successful. 
 *  - CMConnection initiate_conn(CManager cm, CMtrans_services svc,
 *                               transport_entry trans, attr_list attrs);
 *      This routine should initiate a connection to the host/process
 *      specified by the attrs parameters.  The return value is a
 *      CMConnection whose private data will be specific to this
 *      particular connection (which will be provided to routines
 *      below as 'conn_data').  The routine should also perform
 *      whatever tasks are necessary for servicing this connection
 *      (e.g. adding the appropriate FD to the select() list,
 *      establishing a periodic task that will check for data, etc.)
 *      Generally, when data is available on a connection, a call to
 *      trans->data_available() should be performed.
 *  - CMConnection initiate_conn_nonblocking(CManager cm, CMtrans_services svc,
 *                               transport_entry trans, attr_list attrs, int condition);
 *      This routine should initiate a connection to the host/process
 *      specified by the attrs parameters.  The return value is a
 *      CMConnection whose private data will be specific to this
 *      particular connection (which will be provided to routines
 *      below as 'conn_data').  The routine should also perform
 *      whatever tasks are necessary for servicing this connection
 *      (e.g. adding the appropriate FD to the select() list,
 *      establishing a periodic task that will check for data, etc.)
 *      Generally, when data is available on a connection, a call to
 *      trans->data_available() should be performed.
 *  - int self_check(CManager cm, CMtrans_services svc,
 *                   transport_entry trans, attr_list attrs);
 *      Because only the individual CMtransports can fully interpret
 *      the attribute lists that comprise connection information,
 *      layers above sometimes don't know if a particular attribute
 *      list is actually a reference to itself.  This routine asks
 *      the question "Am I the host/process referenced by the contact
 *      list 'attr'?".  Return value is 1 for true and 0 for false.
 *  - int connection_eq(CManager cm, CMtrans_services svc,
 *                      transport_entry trans, attr_list attrs,
 *                      void *conn_data);
 *      This routine is similar to self_check, but is used to avoid
 *      creating duplicate connections between communicating
 *      processes.  It asks the question "If I were to use 'attrs' to
 *      initiate a connection, would I be connecting to a destination
 *      already represented by this connection?"  Generally, CM tracks
 *      what attribute lists had been used to create outbound
 *      connections, so this routine is most necessary to identify
 *      incoming connections that have been established passively.
 *      (I.E. as the result of a listen/accept.)  (CM assumes that
 *      CMConnections are bidirectional.)  Note that many transports
 *      will require an exchange of contact information at the time
 *      that a connection is initiated in order to support this call.
 *  - void shutdown_conn(CMtrans_services svc, void *conn_data);
 *      This routine should close the connection associated with the
 *      transport-private data conn_data.  This includes removing it
 *      from service and deallocating all resources associated with it
 *      (including the conn_data structure).
 *  - int read_to_buffer(CMtrans_services svc, void *conn_data, void *buffer, 
 *                       int len, int block_flag);
 *      There are two basic "read" calls that might be provided by a
 *      transport, read_to_buffer() and read_block().
 *      read_to_buffer() is designed for use by streaming transports
 *      such as TCP which do not have transport-imposed message
 *      boundaries.  In this case, CM will provide overall buffer
 *      management (allocating, extending and managing the message
 *      buffer), and the transport's responsibilities are simply to
 *      drop a specified number of bytes at the address provided.
 *      Reading a complete message often requires a number of calls to
 *      read_to_buffer().  The block_flag parameter specifies whether
 *      or not the call should block waiting on the specified amount
 *      of data, and the return value indicates the number of bytes
 *      actually read and placed in the buffer.  A return value of -1
 *      indicates a fatal error and initiates a shutdown of the
 *      connection.
 *  -  CMBuffer read_block(CMtrans_services svc, void *conn_data, int *actual, int *offset);
 *      This "read" call is designed for use by transports which have
 *      a stronger sense of message boundaries and may need to manage
 *      their own buffers.  Generally this is called by trans->
 *      data_available(), returns a pointer to a transport-managed
 *      memory region and sets the integer pointed to by actual to the
 *      number of valid bytes in the message.  As with
 *      read_to_buffer(), a length of -1 or a return value of NULL,
 *      indicates a fatal error and connection shutdown is initated.
 *      A length of 0 indicates that a complete message has not yet
 *      been received.  The offset value should be the offset of the start
 * 	of data (in case the transport requires a header). The return value 
 *      should be a CMBuffer value containing the data (or partial data) 
 *	available on the connection.
 *  -  int writev(CMtrans_services svc, void *conn_data, void *iov, 
 *                int iovcnt, attr_list attrs);
 *      writev() is the basic vector write function (similar to Posix
 *      writev()).  It takes a vector of buffers and a count of
 *      vectors and is expected to transfer all across the
 *      connection.  This attr_list is designed to specify characteristics
 *      of the transport of this message (such as priority, reliability,
 *      etc.)  The parameter may be NULL and may be ignored by transport
 *      that do not support such characterstics.  This call is blocking and
 *      should write all bytes and return the number of complete vectors
 *      written.  Writing less than the requested number of vectors
 *      indicates a fatal error and will likely initiate the shutdown
 *      of the connection.
 * - int NBwritev(CMtrans_services svc, void *conn_data, void *iov, 
 *                int iovcnt, attr_list attrs);
 *      NBwritev() is the non-blocking version of writev().
 *      CM's non-blocking write support is experimental, is not
 *      currently well tested and is not enabled by default.  This
 *      routine differs from writev() in that its return value is
 *      the number of bytes (not vectors) written, and that a return
 *      of less than the requested count is not a fatal error.
 *      Instead, the remaining bytes will be copied and queued for a
 *      later write.
 *  -  int writev_complete_notify(CMtrans_services svc, void *conn_data, void *iov, 
 *                int iovcnt, attr_list attrs, CMcompletion_notify_func func, void*client_data);
 *      writev_complete_notify() differs from writev() in that it does *not*
 *      enforce the semantic that the write is essentially complete
 *      (I.E. data can be immeditately overwritten) when the function
 *      returns.  Instead, the higher level is assuring the transport that
 *      the data will remain valid until such time as the transport
 *      indicates that the write is complete by calling the
 *      CMcompletion_notify_func that was specified in the call).  Providing
 *      this call is OPTIONAL for transports, but CM and EVPath will attempt
 *      to use it where possible.
 * - int set_write_notify(transport_entry trans, CMtrans_services svc,
 *                        void *conn_data, int enable);
 *      This routine is used if a non-blocking write fails to write
 *      all of it's data and some is queued.  This routine should
 *      enable notifications on the connection so that the routine
 *      trans->write_possible() is called when more data can be
 *      successfully written on the connection.  Both
 *      set_write_notify() and NBwritev() must be exported for
 *      non-blocking writes to be possible with a given transport.
 * - attr_list get_transport_characteristics(transport_entry trans, 
 *                        CMtrans_services svc, void *conn_data);
 *      This routine is used to provide CM with characteristics of a 
 *      transport that are not required for connections, such as reliability, 
 *      multicast nature, etc.
 *
 * 	In addition to the calls above that are exported by the transport
 * DLL, there are several other entries in the transport_item data
 * structure.  The "trans_name" and "cm" entries are obvious.  "trans_data"
 * is a void pointer to the transport data structure (as returned by
 * initialize).  The other two entries are "data_available" and
 * "write_possible".  These function pointers which are really upcalls
 * provided by CM to let the transport notify CM of certain conditions on
 * specific connections:  I.E. that data is a available on a specific
 * connection, or that a write is now possible (after a non-blocking write
 * has failed to write all data).  These function profiles are:
 *  void (*data_available)(transport_entry trans, CMConnection conn);
 *  void (*write_possible)(transport_entry trans, CMConnection conn);
 * 
 *
 *
 *    The CMtrans_services data structure, passed in to all the transport
 *    functions above as the 'svc' parameter, contains functions that
 *    transports can use to manipulate CM-level data structures, schedule
 *    tasks, etc.  These are:
 *  - void* (*malloc_function)(int size);
 *  - void* (*realloc_function)(void *ptr, int size);
 *  - void (*free_function)(void *ptr);
 *    The malloc/realloc/free function set duplicate the standard library
 *    functions.   (No strong reason at the moment to use these instead of
 *    calling directly.)
 *  - void (*fd_add_select)(CManager cm, int fd, select_list_func func, 
 *                     void *param1, void *param2);
 *    This upcall adds a function to CM's select() list.  The function
 *    'func' will be called with 'param1' and 'param2' when select() detects
 *    that 'fd' has data available.  The most common use of this function is
 *    to pass trans->data_available as 'func', the 'trans' entry as param1
 *    and the CM-level connection structure 'conn' as param2.
 *  - void (*fd_remove_select)(Cmanager cm, int fd)
 *    This upcall clears the select list entry for 'fd'.
 *  - void (*fd_write_select)(Cmanager cm, int fd)
 *    This upcall adds a function to CM's select() list for the 'write
 *    possible' condition.  The function 'func' will be called with 'param1'
 *    and 'param2' when select() detects that 'fd' is in a 'write possible'
 *    condition.  The most common use of this function is to pass
 *    trans->write_possible as 'func', the 'trans' entry as param1 
 *    and the CM-level connection structure 'conn' as param2.
 *  - (*trace_out)(CManager cm, char *format, ...)
 *    This is a printf()-like call that will result in debugging output IFF
 *    the CMTransportVerbose environment variable is set.
 *  - CMConnection (*connection_create)(transport_entry trans, 
 *                                      void *conn_data, attr_list conn_list)
 *    This upcall creates a CM-level CMConnection value to be associated
 *    with a new transport-level structure.  The 'conn_data' provided will
 *    be passed as 'conn_data' to other calls.  The conn_list set of
 *    attributes will be returned by CMConnection_get_attrs().  This call is
 *    typically used by the transport 'initate_conn' function, creating a
 *    CMConnection value that will be returned directly from the function.
 *    It is also employed when 'accepting' a new connection.  In this later
 *    case, the CMConnection value is not returned directly to any caller as
 *    accepting a connection is an asynchronous operation in CM.
 *  - void (*add_shutdown_task)(CManager cm, CMPollFunc func, void *data, int task_type)
 *    Add a function that will be called when the CManager is shutdown. The
 *    function is called with 'cm' and 'data' as its parameters.  Mostly
 *    useful for freeing transport data.  Task type should be FREE_TASK or 
 *    SHUTDOWN task.  SHUTDOWN tasks will be called at first close.  FREE tasks
 *    will be called after ref counts go to zero.
 *  - void (*add_periodic_task)(CManager cm, int period_sec, 
 *                              int period_usec, CMPollFunc func, void *data)
 *    Add a function that will be called at regular intervals, specified by
 *    period_sec and period_usec. The function is called with 'cm' and
 *    'data' as its parameters.  Note that the periodicity is not
 *    guaranteed.  This is serviced by the network handler thread, or by
 *    non-threaded applications calling some function that causes CM to
 *    service the network.  If those threads are otherwise engaged, the
 *    function won't be invoked until the network is to be serviced again.
 *  - void (*connection_close)(CMConnection conn)
 *    Call CMConnection_close() on the particular connection.  Generally,
 *    transports should not call this, but should let the higher levels
 *    close a connection in response to a reported error.  However, if a
 *    transport comes to know about a failure outside of a read/write
 *    situation, this can be called to close the connection.
 *  - CMbuffer (*get_data_buffer)(CManager cm, int length)
 *      There are several buffer management calls for use by CMtransports.
 *    Generally stream-oriented transports like TCP sockets don't need these
 *    calls at all.  If a transport can read a few bytes at a type from an
 *    incoming data stream and efficiently deposit those bytes into whatever
 *    memory destination memory is needed, the transport can respond to
 *    read_to_buffer() requests and let CM manage the memory for incoming
 *    messages.  If not, then the transport should export a 'read_buffer'
 *    interface and return a CMBuffer value.  CMBuffer values generally must
 *    be managed by CM so that reference counts, deallocation,
 *    CMtake_buffer, etc. all work, so the buffer management calls allow a
 *    transport to work with CM to create appropriate CMbuffer values.  The
 *    simplest of these calls is 'get_data_buffer'.  This allows a transport
 *    to create a CMBuffer value of a specific size, using CM-allocated
 *    memory, and give up control of that buffer once it has been filled
 *    with data.  Here the transport controls only the size of the
 *    allocation (as opposed to allowing CM to realloc() the buffer manually
 *    as data comes in).  The get_data_buffer() call is currently used by
 *    the cmudp transport.
 *  - CMbuffer (*create_data_buffer)(CManager cm, void *buffer, int length)
 *  - CMbuffer (*create_data_and_link_buffer)(CManager cm, void *buffer, int length)
 */
