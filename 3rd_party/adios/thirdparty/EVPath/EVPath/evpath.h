
#ifndef __EVPATH__H__
#define __EVPATH__H__
/*! \file */

#include "ffs.h"
#include "atl.h"
#ifdef	__cplusplus
extern "C" {
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif




struct _CManager;
struct _CMConnection;
struct _CMFormat;

/*!
 * CManager is the root of control flow and message handling in a CM program.
 *
 * CManager is an opaque handle.  
 */
typedef struct _CManager *CManager;

/*!
 * CMConnection is a handle to a communications link.
 *
 * CManager is an opaque handle.  
 */
typedef struct _CMConnection *CMConnection;

/*!
 * CMFormat is a handle to a registered native format.
 *
 * CMFormat is an opaque handle.  It is the return value from
 * CMregister_format() and is used both to identify data for writing (in
 * CMwrite() and CMwrite_attr() and to register handlers for incoming data
 * (with CMregister_handler()).
 */
typedef struct _CMFormat *CMFormat;

/*!
 * CMTaskHandle is a handle to a delayed or periodic task.
 */
typedef struct _CMTaskHandle *CMTaskHandle;

/*!
 * buf_entry is a structure used to return the lengths of encoded events
 * and a pointer to their locations.
 */
typedef struct buf_entry {
    size_t length;	/*!< length of the encoded buffer */
    void *buffer;	/*!< base address of the encoded buffer */
} *EVevent_list;

/*!
 * The prototype for a CM data handling function.
 *
 * CM allows application-routines matching this prototype to be registered
 * as data handlers.
 * \param cm The CManager with which this handler was registered.
 * \param conn The CMConnection upon which the message arrived.
 * \param message A pointer to the incoming data, cast to void*.  The real
 * data is formatted to match the fields of with which the format was
 * registered. 
 * \param client_data This value is the same client_data value that was
 * supplied in the CMregister_handler() call.  It is not interpreted by CM,
 * but instead can be used to maintain some application context.
 * \param attrs The attributes (set of name/value pairs) that this message
 * was delivered with.  These are determined by the transport and may
 * include those specified in CMwrite_attr() when the data was written.
 */
typedef void (*CMHandlerFunc) (CManager cm, 
			       CMConnection conn,
			       void *message, void *client_data,
			       attr_list attrs);

/*!
 * The prototype for a CM polling handler (and others).
 *
 * Functions matching of this prototype can be registered with CMadd_poll(),
 * CMadd_periodic_task(), CMadd_delayed_task() and CMadd_shutdown_task().
 * \param cm The CManager with which this handler was registered.
 * \param client_data This value is the same client_data value that was
 * supplied in the CMadd_poll() call.  It is not interpreted by CM,
 * but instead can be used to maintain some application context.
 */
typedef void (*CMPollFunc) (CManager cm, void *client_data);

/*!
 * The prototype for a CM connection close handler.
 *
 * Functions matching of this prototype can be registered with
 * CMregister_close_handler(). 
 * \param cm The CManager with which this handler was registered.
 * \param conn The CMConnection which is being closed.
 * \param client_data This value is the same client_data value that was
 * supplied in the CMregister_close_handler() call.  It is not interpreted 
 * by CM, but instead can be used to maintain some application context.
 */
typedef void (*CMCloseHandlerFunc) (CManager cm, CMConnection conn,
				    void *client_data);

/*!
 * The prototype for a CM write possible callback.
 *
 * Functions matching of this prototype can be registered with
 * CMregister_write_callback(). 
 * \param cm The CManager with which this callback function was registered.
 * \param conn The CMConnection upon which a non-blocking write is now (potentially) possible. 
 * \param client_data This value is the same client_data value that was
 * supplied in the CMregister_write_callback() call.  It is not interpreted 
 * by CM, but instead can be used to maintain some application context.
 */
typedef void (*CMWriteCallbackFunc) (CManager cm, CMConnection conn,
				     void *client_data);



/*!
 * create a CManager.
 *
 * CManager is the root of control flow and message handling in a CM program.
 */
/*NOLOCK*/
extern CManager CManager_create();

/*!
 * create a CManager.
 *
 * CManager is the root of control flow and message handling in a CM program.
 * This version of a call selects a specific control_module for managing waiting 
 * network and control events, currently either "select" or "epoll".
 */
/*NOLOCK*/
extern CManager CManager_create_control(char *control_module);

/*!
 * close a CManager
 *
 * the close operation shuts down all connections and causes the
 * termination of any network handling thread.
 * \param cm The CManager to be shut down.
 */
extern void CManager_close (CManager cm);

/*!
 * specify a numerical identifier to be used as part of a trace output filename
 *
 * \param ID    the numerical identifier;
 *
 *
 */
/*NOLOCK*/
extern void CMTrace_file_id (int ID);

/*!
 * fork a thread to handle the network input operations of a CM.
 *
 * \param cm The CManager whose input should be handled by a new thread.
 * \return 
 * - 0 if a communications manager thread cannot be forked
 * - 1 success
 * \warning Only one thread should be handling input for a CM.  If this call
 * is used then no other thread should call CMpoll_network() or
 * CMrun_network(). 
 * \warning If this call is to be used (or indeed if any threading is to be
 * used), one of the gen_thread init routines should be called <b>before</b>
 * the call to CManager_create().  Otherwise bad things will happen.
 */
extern int CMfork_comm_thread (CManager cm);

/*!
 * Tell CM to listen for incoming network connections.
 *
 * \param cm The CManager which should listen.
 * \return the number of transports which successfully initiated connection
 * listen operations (by reporting contact attributes).
 * \note CMlisten() is identical to calling CMlisten_specific() with a NULL
 * value for the listen_info parameter.
 * \note The listening transports will add their contact information to the
 * list returned by CMget_contact_list().
 */
extern int CMlisten (CManager cm);

/*!
 * Tell CM to listen for incoming network connections with 
 * specific characteristics.
 *
 * \param cm The CManager which should listen.
 * \param listen_info An attribute list to be passed to the
 * transport-provided listen operations of all currently-loaded transports.
 * \return the number of transports which successfully initiated connection
 * listen operations (by reporting contact attributes).
 * \note The listen_info value is interpreted by each individual transport.
 * An incomplete (and probably dated) list of transports that use this include: 
 * - the <b>sockets</b> transport which uses the CM_IP_PORT attribute to control
 *   which port it listens on.  If this attribute is not present it listens
 *   on any available port. 
 * - the <b>rudp</b> transport which uses the CM_UDP_PORT attribute to control
 *   which port it listens on.  If this attribute is not present it listens
 *   on any available port. 
 * - the <b>enet</b> transport which uses the CM_ENET_PORT attribute to control
 *   which port it listens on.  If this attribute is not present it listens
 *   on any available port. 
 * - the <b>udp</b> transport - a raw unreliable UDP transport.
 * - the <b>multicast</b> transport  - a raw unreliable Multicast transport.
 * - the <b>nnti</b> transport  - a multi-network RDMA transport.
 * - the <b>ib</b> transport - kind-of-functional InfiniBand transport.
 */
extern int CMlisten_specific (CManager cm, attr_list listen_info);

/*!
 * get the contact information for this CM.
 *
 * This call returns the set of attributes that define the contact
 * information for the network transports that have performed listen
 * operations.  
 * \param cm the CManager for which to return contact information.
 * \return the contact list.
 */
extern attr_list
CMget_contact_list (CManager cm);

/*!
 * insert contact information into this CM.
 *
 * This call adds to the set of attributes that define the contact
 * information for the network transports that have performed listen
 * operations.  
 * \param cm the CManager for which to add contact information.
 * \param attrs the information to add.
 */
extern void
CM_insert_contact_info (CManager cm, attr_list attrs);

/*!
 * get a specific subset of the contact information for this CM.
 *
 * This call returns the set of attributes that define the contact
 * information for a particular network transport.  If no listen operation
 * has been performed on that transport, one will be done with a NULL attr
 * list. 
 * \param cm the CManager for which to return contact information.
 * \param attrs the attribute list specifying the transport.
 * \return the contact list.
 */
extern attr_list
CMget_specific_contact_list (CManager cm, attr_list attrs);

/*!
 * get a thread-owned version of a shared attribute list in a thread-safe way
 * \param cm  the CManager which owns the attribute list.
 * \param attrs the shared attribute list (from CMget_contact list, etc.)
 * \return a single-owner contact list
 */
extern attr_list
CMderef_and_copy_list(CManager cm, attr_list attrs);

/*!
 * check to see if this is contact information for <b>this</b> CM.
 *
 * Since attribute lists are generally opaque, it is not necessarily obvious
 * when you have a set of contract attributes that is actually your own
 * contact list.  This call is designed to answer that question.
 *
 * \param cm The CManager whose contact information should be compared.
 * \param attrs The contact list to compare.
 * \return 1 if for some loaded transport the attrs list matches the contact
 * information in the cm. 0 otherwise.
 */
extern int
CMcontact_self_check (CManager cm, attr_list attrs);

/*!
 * acquire a (possibly existing) connection to another CM process.
 *
 * \param cm The CManager in which to make the connection.
 * \param contact_list The attribute list specifying contact information for
 * the remote CM.
 * \return A CMConnection value, or NULL in the event of failure.
 *
 * CMget_conn() attempts to determine if the contact attribute match any
 * existing connection (Using the transport connection_eq() method).  If a
 * connection matches, that connection's reference count is incremented and
 * its value is returned.  If no connection matches, a CMinitiate_conn() is
 * performed using the contact list and its result value is returned.
 */
extern CMConnection
CMget_conn (CManager cm, attr_list contact_list);

/*!
 * initiate connection to another CM process.
 *
 * \param cm The CManager in which to make the connection.
 * \param contact_list The attribute list specifying contact information for
 * the remote CM.
 * \return A CMConnection value, or NULL in the event of failure.
 *
 * CMinitiate_conn() will attempt to initiate a connection using each of the
 * currently loaded transports.  It will return the first that succeeds, or
 * NULL if none succeed.
 * \note If the contact list contains a CM_TRANSPORT attribute with a string
 * value, CMinitiate_conn() will attempt to load that transport, then if that
 * succeeds will attempt to initiate a connection using only that transport.
 */
extern CMConnection
CMinitiate_conn (CManager cm, attr_list contact_list);

/*!
 * kill and potentially deallocate a connection.
 *
 * \param conn the CMConnection to kill
 *
 * CMConnection_close decrements the reference count of a connection.  If
 * the resulting reference count is zero, then the connection is shut down.
 * All resources associated with the connection are free'd, the close
 * handler is called and the CMConnection structure itself is free'd.
 * \warning CMConnection values should not be used after
 * CMConnection_close().  CMConnection_close() should only be used on 
 * CMConnection values created with CMget_conn() or CMinitiate_conn(), not 
 * with connections that are passively created (accepted through CMlisten()).
*/
extern void
CMConnection_close (CMConnection conn);

/*!
 * manually increment the reference count of a connection.
 *
 * \param conn the CMConnection whose count should be incremented.
 * \note  Used if some mechanism other than CMget_conn() is used to "share"
 * connection in multiple contexts so that it is closed only when all users
 * have closed it.
*/
extern void
CMConnection_add_reference (CMConnection conn);

/*!
 * manually decrement the reference count of a connection.
 *
 * \param conn the CMConnection whose count should be decremented.
 * \note  Used if some mechanism other than CMget_conn() is used to "share"
 * connection in multiple contexts so that it is closed only when all users
 * have closed it.
*/
extern void
CMConnection_dereference (CMConnection conn);

/*!
 * register a function to be called when a connection is closed.
 *
 * \param conn the connection with which the function is associated.
 * \param func the function to be called when the connection closes.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * \note There can be only one close handler per connection.  Multiple
 * registrations overwrite each other.
 */
extern void
CMconn_register_close_handler (CMConnection conn, 
			       CMCloseHandlerFunc func, 
			       void *client_data);
/*!
 * return the list of attributes associated with a connection.
 *
 * \param conn the connection for which to return the attributes.
 * \return an attr_list value containing connection attributes.
 */
extern attr_list 
CMConnection_get_attrs (CMConnection conn);

/*!
 * modify the characteristics of a connection.
 *
 * \param conn the connection for to modify characteristics.
 * \param attrs the characteristics to apply (specific to CM and transport).
 * \return a true/false failure value.
 */
extern int
CMConnection_set_character (CMConnection conn, attr_list attrs);

/*!
 * return connection 'i' associated with a CM value.
 *
 * \param cm the cmanager from which to return a connection.
 * \param i the index value into the CManager's list of connections.
 * \return a CMConnection value associated with connection 'i'
 */
extern CMConnection
CMget_indexed_conn (CManager cm, int i);

/*!
 * register a format with CM.
 *
 * \param cm  The CManager in which to register the format.
 * \param format_list The FM format list which describes the structure.
 * It should contain the transitive closure of all data types necessary to
 * specify the message representation.  The list is terminated with a
 * <tt>{NULL, NULL, 0, NULL}</tt> value.
 *
 * Registering a format is a precursor to sending a message or registering a
 * handler for incoming messages.
 */
extern CMFormat
CMregister_format (CManager cm, FMStructDescList format_list);

/*!
 * register a simple (no internal structures) format with CM.
 *
 * \param cm  The CManager in which to register the format.
 * \param format_name The name of the data/message structure being registered
 * \param field_list The FM field list which describes the structure, listing all 
 * structure fields (their names, data types, offsets and sizes).
 * As with all FMFieldLists, the list is terminated with a <tt>{NULL, NULL, 0, 0}</tt> value.
 * \param struct_size The (padded if necessary) size of the structure, typically as returned by sizeof().
 *
 * Registering a format is a precursor to sending a message or registering a
 * handler for incoming messages.  This call is the equivalent to calling CMregister_format(), 
 * specifying a single entry in the format_list parameter with the format_name, field_list and 
 * struct_size and the opt_info specified as NULL.
 */
extern CMFormat
CMregister_simple_format (CManager cm, char *format_name, FMFieldList field_list, int struct_size);

/*!
 * lookup the CMFormat associated with a particular FMStructDescList
 *
 * \param cm The CManager in which the format was registered.
 * \param format_list The format list which was used in the registration.
 *
 * CMLookup_format() addresses a specific problem particular to libraries.
 * CMwrite() requires a CMFormat value that results from a
 * CMregister_format() call.  Efficiency would dictate that the
 * CMregister_format() be performed once and the CMFormat value used
 * repeatedly for multiple writes.  However, libraries which want to avoid
 * the use of static variables, or which wish to support multiple CM values
 * per process have no convenient way to store the CMFormat values for
 * reuse.   CMlookup_format() exploits the fact that field_lists are
 * often constants with fixed addresses (I.E. their memory is not reused for
 * other field lists later in the application).  This call quickly looks up
 * the CMFormat value in a CM by searching for a matching field_list
 * address. 
 * \warning You should *not* use this convenience routine if you cannot
 * guarantee that all field lists used to register formats have a unique
 * address. 
 */
extern CMFormat CMlookup_format (CManager cm, FMStructDescList format_list);

/*!
 * return the FMContext value used by a CM.
 *
 * \param cm The CManager from which to return the FMcontext;
 *
 * CMget_FMcontext() returns the FMcontext value in use by the CManager
 */
/*NOLOCK*/
extern FMContext CMget_FMcontext(CManager cm);

/*!
 * send a message on a connection.
 *
 * \param conn The CMConnection upon which to send the message.
 * \param format The CMFormat value returned by CMregister_format().
 * \param data The unencoded message, cast to a <tt>void*</tt> value.
 * \return
 * - 1 if the write was successful.
 * - 0 if the write did not complete successfully.
 * \note CMwrite() is equivalent to CMwrite_attr() with a NULL value 
 * passed for the attrs parameter.
 */
extern int
CMwrite (CMConnection conn, CMFormat format, void *data);

/*!
 * send a message on a connection with particular attributes.
 *
 * \param conn The CMConnection upon which to send the message.
 * \param format The CMFormat value returned by CMregister_format().
 * \param data The unencoded message, cast to a <tt>void*</tt> value.
 * \param attrs The set of name/value attributes with which to write the data.
 * \return
 * - 1 if the write was successful.
 * - 0 if the write did not complete successfully.
 * \note The values in the attrs parameter serve two purposes.  First, 
 * they may be interpreted by CM or the CM transport layers on either 
 * the writing or reading sides to customize event delivery.  Second, 
 * they are made available (perhaps with additional transport-specific 
 * attributes) to the read-side handler in the attrs argument to the 
 * CMHandlerFunc that handles the message.
 * \note CMwrite_attr() with a NULL value for the attrs parameter is 
 * equivalent to CMwrite().
 */
extern int
CMwrite_attr (CMConnection conn, CMFormat format, void *data, 
	      attr_list attrs);

/*!
 * register a function to be called when message matching a particular 
 * format arrives. 
 *
 * \param format The CMFormat value returned by CMregister_format()
 * \param handler The function to be called to handle the message.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 */
extern void
CMregister_handler (CMFormat format, CMHandlerFunc handler, 
		    void *client_data);

/*!
 * register a function to be called when a write is again possible on a particular CMconnection.
 *
 * \param conn The CMConnection upon which to send the message.
 * \param handler The function to be called to handle the message.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 *
 * 
 */
extern void
CMregister_write_callback (CMConnection conn, 
			   CMWriteCallbackFunc handler,
			   void *client_data);

/*!
 * test whether a write to a particular connection would block
 *
 * \param conn The CMConnection to test
 * \return boolean TRUE(1) if the write would certainly block and 
 *   FALSE(0) if it may not.  
 * At the network level, we likely only know that <B>something</B> can be
 * sent, not how much.  So even if this returns false, a write may still
 * block at the network level.  If this happens, CM will copy the remaining
 * bytes and allow the CMwrite() to return, finishing the send
 * asynchronously.  However, if a CMwrite() is initiated when
 * write_would_block() is already TRUE, the write <b>will block</b> until
 * the blocking condition is gone (I.E. CMConnection_write_would_block() is
 * again FALSE).)
 */
extern int
CMConnection_write_would_block (CMConnection conn);

/*!
 * assume control over a incoming buffer of data.
 *
 * This call is designed to be used inside a CMHandlerFunc.  Normally data
 * buffers are recycled and CM only guarantees that the data delivered to a
 * CMHandlerFunc will be valid for the duration of the call.  In that
 * circumstance, a handler that wanted to preserve the data for longer than
 * its own duration (to pass it to a thread or enter it into some other data
 * structure for example) would have to copy the data.  To avoid that
 * inefficiency, CMtake_buffer() allows the handler to take control of the
 * buffer holding its incoming data.  The buffer will then not be recycled
 * until it is returned to CM with CMreturn_buffer().
 * \param cm The CManager in which the handler was called.
 * \param data The base address of the data (I.E. the message parameter to
 * the CMHandlerFunc).
 * \return NULL on error, otherwise returns the data parameter. 
*/
extern void *CMtake_buffer (CManager cm, void *data);

/*!
 * return a buffer of incoming data.
 *
 * This call recycles a data buffer that the application has taken control
 * of through CMtake_buffer().
 * \param cm The CManager in which the handler was called.
 * \param data The base address of the data (I.E. same value that was passed
 * to CMtake_buffer().
*/
extern void CMreturn_buffer (CManager cm, void *data);

#include "cm_transport.h"
/*!
 * The prototype for a non-CM message handler.
 *
 * Functions matching of this prototype can be registered with
 * CMregister_non_CM_message_handler().
 * \param conn The CMConnection on which the message is available.
 * \param header The first 4 bytes of the message, encoded as an integer.
 * \return 0 if the message was completely handled.  
 *   Otherwise return the number of additional bytes necessary.
 */
typedef int (*CMNonCMHandler) (CMConnection conn, CMTransport transport,
			       char *buffer, size_t length);

/*!
 * register a handler for raw (non-CM) messages.
 *
 * CM, like may other message layers, embeds a unique value in the first 4
 * bytes of the incoming message to identify it as a CM message.  CM
 * actually has several sets of identifying 4-byte values that it recognizes
 * as CM-internal messages.  This interface can be used to add to that set
 * to include non-CM messages (such as IIOP, HTTP, etc.).  
 * \param header The 4 bytes that identify (encoded as an integer) that
 * identify the messages to be handled.
 * \param handler The handler to be called when messages with this header
 * arrive. 
 * \note Registration is not CManager-specific, but apply to all CManagers
 * in the address space (ugly).
 * \warning Don't do this at home kids!  This API is not complete enough to
 * actually implement something like IIOP externally, but it's the thought
 * that counts.
 */
/*NOLOCK*/
extern void
CMregister_non_CM_message_handler (int header, CMNonCMHandler handler);

/*!
 * register a handler for CM messages that don't match registered handlers.
 *
 */
typedef void (*CMUnregCMHandler) (CMConnection conn, char *format_name);

/*NOLOCK*/
extern void
CMregister_invalid_message_handler (CManager cm, CMUnregCMHandler handler);

/*!
 * return the pointer to the static transport services structure.
 *
 * All transports share a set of basic services provided by CM, whose function
 * pointers are available through this structure.
 *
 * \return returns the pointer to the services structure.
 */
/*NOLOCK*/
extern CMtrans_services
CMget_static_trans_services ();

  /*!
   * return the pointer to a CMConnection's transport data.
   *
   * Think of this structure as the cross-product of a transport and CMConnection.
   * Transport functions use this structure to store per-connection data.
   *
   * \return returns the pointer to the transport data structure.
   */
extern void*
CMget_transport_data (CMConnection conn);

/*!
 * add a task (function) to be executed occasionally.
 *
 * \param cm The CManager to which the task is added.
 * \param func The function to be called occasionally.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * CM poll functions are called after each round of message delivery.  I.E. 
 * once per call to CMpoll_network() if that function is used.
 */
extern void
CMadd_poll (CManager cm, CMPollFunc func, void *client_data);

/*!
 * add a task (function) to be executed with a specified periodicity.
 *
 * \param cm The CManager to which the task is added.
 * \param period_sec The number of whole seconds of the period.
 * \param period_usec The number of additional microseconds of the period.
 * \param func The function to be called periodically.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * \return a CMTaskHandle which can be used to remove the task.
 * \note CM does not guarantee a particular periodicity, it merely applies
 * its best efforts.  I.E. It will not block wait in select() past the
 * timeout period for the next task.  However handlers may run long and I/O
 * may intervene to delay the task execution.  The task will be executed
 * when the first opportunity arises after it is scheduled.  After execution
 * is complete, the next execution will be scheduled based upon the actual
 * execution time of the current invocation (not when it was scheduled to be
 * executed). 
 */
extern CMTaskHandle
CMadd_periodic_task (CManager cm, int period_sec, int period_usec, 
		     CMPollFunc func, void *client_data);

/*!
 * add a task (function) to be executed at a later time.
 *
 * \param cm The CManager to which the task is added.
 * \param secs The number of whole seconds to delay the task.
 * \param usecs The number of additional microseconds to delay the task.
 * \param func The function to be called after the delay.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * \return a CMTaskHandle which can be used to remove the task (only before
 * it executes).
 * \note CM does not guarantee a particular delay, it merely applies
 * its best efforts.  I.E. It will not block wait in select() past the
 * timeout period for the next task.  However handlers may run long and I/O
 * may intervene to delay the task execution.  The task will be executed
 * when the first opportunity arises after it is scheduled.  
 */
extern CMTaskHandle
CMadd_delayed_task (CManager cm, int secs, int usecs, CMPollFunc func,
		    void *client_data);

/*!
 * remove a registered periodic or delayed task.
 *
 * \param handle The handle to the task to remove.
 */
extern void
CMremove_task (CMTaskHandle handle);

/*!
 * add a task (function) to be called when the CM is shut down.
 *
 * \param cm The CManager to which a shutdown task is added.
 * \param func The function to be called upon shutdown.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * \param task_type Should be either SHUTDOWN_TASK or FREE_TASK
 *
 * Multiple shutdown tasks can be added to the same CM and they are called
 * in the order registered.  There is currently no API for removing them.
 */
extern void
CMadd_shutdown_task (CManager cm, CMPollFunc func, void *client_data, int task_type);

/*!
 * task type for CMadd_shutdown_task.  NO_TASK is used internally 
 */
#define NO_TASK 0
/*!
 * task type for CMadd_shutdown_task.  SHUTDOWN_TASKs will be called at first close.
 */
#define SHUTDOWN_TASK 1
/*!
 * task type for CMadd_shutdown_task.  FREE_TASKs will be called when the CM reference count reaches zero.
 */
#define FREE_TASK 2

/*!
 * add a task to be executed with a particular periodicity.
 *
 * \param cm The CManager which should execute the task.
 * \param period The period of the task in microseconds.
 * \param func The function to be called.
 * \param client_data An uninterpreted value that is passed to the function
 * when it is called.
 * \return a CMTaskHandle which can be used to remove the task.
 * \deprecated Use CMadd_periodic_task().
 */
extern CMTaskHandle
CMadd_periodic (CManager cm, long period, CMPollFunc func,
		void *client_data);

/*!
 * remove a registered periodic task.
 *
 * \param handle The handle to the task to remove.
 * \deprecated Use CMremove_task()
 */
extern void
CMremove_periodic (CMTaskHandle handle);

/*!
 * sleep for a given number of seconds.
 *
 * Unlike system sleep() calls, CMsleep() will continue to handle network
 * messages during the sleep time.  In particular, if CMsleep is called by
 * the network handler thread or in a single threaded program, then it will
 * enter a network handling loop until the time has elapsed.  If called by
 * other than the network handler thread in a multithread application, then
 * it will suspend on a thread condition wait until the time has elapsed.
 * \param cm The CManager upon which to sleep.
 * \param secs The number of seconds for which to sleep.
 */
extern void
CMsleep (CManager cm, int secs);

/*!
 * sleep for a given number of microseconds.
 *
 * Unlike system sleep() calls, CMusleep() will continue to handle network
 * messages during the sleep time.  In particular, if CMusleep is called by
 * the network handler thread or in a single threaded program, then it will
 * enter a network handling loop until the time has elapsed.  If called by
 * other than the network handler thread in a multithread application, then
 * it will suspend on a thread condition wait until the time has elapsed.
 * \param cm The CManager upon which to sleep.
 * \param usecs The number of microseconds for which to sleep.
 */
extern void
CMusleep (CManager cm, int usecs);

/*!
 * handle one round of network events
 *
 * \param cm The CManager for which to handle events.
 * CMpoll_network()} is one of the basic <b>network event</b> handling calls
 * in CM.  A CM network event is a basic communications occurrence, such as
 * a connection request or message arrival. The routine CMpoll_network()
 * essentially polls the network and handles some pending messages before
 * returning.  
 * \note Not all pending messages will be handled, but generally one message
 * will be handled for each connection upon which input is pending.
 */
extern 
void CMpoll_network (CManager cm);

/*!
 * handle network events until shutdown.
 *
 * \param cm The CManager for which to handle events.
 * CMrun_network()} is one of the basic <b>network event</b> handling calls
 * in CM.  A CM network event is a basic communications occurrence, such as
 * a connection request or message arrival. The routine CMrun_network()
 * essentially handles network events until the CManager is shutdown.
 */
extern 
void CMrun_network (CManager cm);

/*!
 * The prototype for a CM select handling function.
 *
 * CM allows application-routines matching this prototype to be registered
 * and called when a particular file descriptor has data available to read.
 * \param param1 is the value specified as param1 in the CM_fd_add_select() call.
 * \param param2 is the value specified as param2 in the CM_fd_add_select() call.
 */
typedef void (*select_func) (void *param1, void*param2);

/*!
 * Register an application routine to be called with a particular 
 * file descriptor has data available for read.
 *
 * \param cm The CManager with which this handler is to be registered.
 * \param fd The file descriptor to be monitored.
 * \param handler_func The function to be called when data is available.
 * \param param1 The value to be passed as param1 to the handler_func.
 * \param param2 The value to be passed as param2 to the handler_func.
 */
extern void
CM_fd_add_select (CManager cm, SOCKET fd, select_func handler_func,
		  void *param1, void *param2);

/*!
 * allocate a new CM condition value.
 *
 * \param cm the CManager value in which to allocate the condition.
 * \param dep the CMConnection value upon which the condition depends.
 * \return an integer value representing a CM condition.
 * \note CM condition values are used to cause a thread or program to wait
 * for a particular situation, usually for a message response to arrive.
 * In this case the condition value is acquired before sending the request
 * message, integer condition value is sent as part of the request and
 * returned in the response.  The response handler then does a
 * CMCondition_signal() as part of its operation.
 * \note The dep CMConnection value is used in error handling.  In
 * particular, if that connection dies or is closed, the condition will be
 * marked as <b>failed</b> and the corresponding CMCondition_wait() will
 * return.  Thus if the situation in which the condition is used relies upon
 * the continued operation of a connection (such as waiting for a response),
 * then that connection should be specified as the dep parameter in this
 * call.  If there is no such reliance, dep can be NULL.
 */
extern int CMCondition_get (CManager cm, CMConnection dep);

/*!
 * wait for a CM condition value.
 *
 * \param cm the CManager value in which the condition was allocated.
 * \param condition the condition upon which to wait.
 * \return 
 * - 1 if the condition was signalled normally.
 * - 0 if the CMConnection specified as dep in the CMCondition_get()
 *	        call was closed.
 * \note CM condition values are used to cause a thread or program to wait
 * for a particular situation, usually for a message response to arrive.
 * \note CMCondition_wait() is useful because it does the "right thing" in
 * both single-threaded and multi-threaded applications.  In single-threaded
 * applications it enters a network-handling loop until the condition has
 * been signaled.  In applications with a network handler thread, it checks
 * to see if it is being executed by that handler thread.  If it is *not*,
 * then it does a thread condition wait to suspect the thread.  If it is
 * being executed by the network handler thread, then it also enters a
 * network-handling loop until the condition has been signaled.
 * \warning The condition value is considered 'free'd upon return from
 * CMCondition_wait() and should not be used in any subsequent call
 * (including calls to CMCondition_get_client_data(), etc.).
 */
extern int CMCondition_wait (CManager cm, int condition);

/*!
 * signal a CM condition value.
 *
 * \param cm the CManager value in which the condition was allocated.
 * \param condition the condition to be signaled.
 * \note CM condition values are used to cause a thread or program to wait
 * for a particular situation, usually for a message response to arrive.
 * \note CMCondition_signal() notifies CM that the situation needed to
 * satisfy a particular condition variable has occurred and any waiting
 * thread should awaken.
 */
extern void CMCondition_signal (CManager cm, int condition);

/*!
 * set the client_data associated with a condition value.
 *
 * \param cm the CManager value in which the condition is allocated.
 * \param condition the condition with which the client_data should be
 * associated. 
 * \param client_data the value to be associated with the condition.
 * \note The client_data value is not interpreted by CM, but instead
 * provides a mechanism through which information can be conveyed between
 * the requesting thread and response handler.  In a typical usage, the
 * requesting site sets the client_data to the address of storage for a
 * return value.  The response handler then uses
 * CMCondition_get_client_data() to access that address and store the return
 * value in the appropriate location.
 * \warning Calls to CMCondition_set_client_data() should occur between the
 * call to CMCondition_alloc() and CMCondition_wait().  The condition value
 * is considered 'free'd upon return from CMCondition_wait() and should not
 * be used in any subsequent call.  To avoid possible race conditions, calls
 * to CMCondition_set_client_data() should also occur before the CMwrite of
 * the request to ensure that the response doesn't arrive before the client
 * data is set.
 */
extern void CMCondition_set_client_data (CManager cm, int condition,
				       void *client_data);
/*!
 * get the client_data associated with a condition value.
 *
 * \param cm the CManager value in which the condition is allocated.
 * \param condition the condition to query for client_data.
 * \return the client_data value associated with the condition.
 * \note The client_data value is not interpreted by CM, but instead
 * provides a mechanism through which information can be conveyed between
 * the requesting thread and response handler.  In a typical usage, the
 * requesting site sets the client_data to the address of storage for a
 * return value.  The response handler then uses
 * CMCondition_get_client_data() to access that address and store the return
 * value in the appropriate location.
 * \warning Calls to CMCondition_get_client_data() should generally occur
 * in the response handler (as opposed to after CMCondition_wait()).  The
 * condition value is considered 'free'd upon return from CMCondition_wait()
 * and should not be used in any subsequent call.
 */
extern void *CMCondition_get_client_data (CManager cm, int condition);

/*!
 * test whether or not a particular condition has been signaled.
 *
 * \param cm the CManager value in which the condition is allocated.
 * \param condition the condition to test.
 * \return boolean value representing whether or not the condition has been
 * signaled. 
 * This call essentially provides a mechanism of examining the state of a
 * condition without blocking on CMCondition_wait().
 * \warning This call should not be used on a condition after
 * a CMCondition_wait() has been performed.
 */
extern int CMCondition_has_signaled (CManager cm, int condition);
/*!
 * test whether or not a particular condition has failed.
 *
 * \param cm the CManager value in which the condition is allocated.
 * \param condition the condition to test.
 * \return boolean value representing whether or not the condition has 
 * failed (I.E. its dependent connection has been closed.)
 * This call essentially provides a mechanism of examining the state of a
 * condition without blocking on CMCondition_wait().
 * \warning This call should not be used on a condition after
 * a CMCondition_wait() has been performed.
 */
extern int CMCondition_has_failed (CManager cm, int condition);

/** @defgroup malloc CM memory allocation functions
 *
 * This group of functions is used to manage CM-returned memory.
 * They are provided to handle the eventuality when CM uses its own memory
 * manager.  That hasn't happened yet, so these are identical to realloc,
 * malloc and free.
 */

/*!
 * reallocate a chunk of memory
 *
 * \param ptr the memory to reallocate
 * \param size the new size
 * \return a pointer to the new block
 */
/*NOLOCK*/
extern void* CMrealloc (void *ptr, long size);
/*!
 * allocate a chunk of memory
 *
 * \param size the requested size
 * \return a pointer to the new block
 */
/*NOLOCK*/
extern void* CMmalloc (long size);
/*!
 * free a chunk of memory
 *
 * \param ptr the memory to free
 */
/*NOLOCK*/
extern void CMfree (void *ptr);

/** @defgroup perf Performance-query functions
 * These functions intrusively test the characteristics of a connection,
 * measuring available bandwidth and current round-trip latency.
 * @{
 */
/*!
 * Probe the approximate round-trip latency on a particular connection by
 * sending a burst of data.
 *
 * This is an intrusive probe.
 * \param conn The CMConnection to be tested.
 * \param msg_size The size of message to be sent in the test.  (Latency
 * varies dramatically with the message size.)
 * \param attrs Currently this parameter is ignored, but it *should* allow
 * control over the number of messages sent.
 * \return The return value is in units of microseconds.
 * \note CM measures latency by sending a message and waiting for a
 * response.  This round-trip is called a "ping".  In the current
 * implementation, CM performs 2 ping operations to "warm up" the 
 * connection.  It then performs 5 additional ping operations, measuring the
 * time required for each.  The return value is the average of these final
 * operations. 
*/
extern long CMprobe_latency (CMConnection conn, long msg_size,
				  attr_list attrs);

/*!
 * Probe the available bandwidth on a particular connection by sending a
 * burst of data.
 *
 * This is an intrusive probe.
 * \param conn The CMConnection to be tested.
 * \param size The size of message to be sent in the test.  (Bandwidth
 * varies dramatically with the message size.)
 * \param attrs Currently this parameter is ignored, but it *should* allow
 * control over the number of messages sent.
 * \return The return value is in units of MBytes per second.  
 * \note In the current implementation, CM sends \f$N\f$ messages to probe
 * available bandwidth, where \f$N\f$ is calculated as \f$100000/size\f$.
 * That is, CMprobe_bandwidth sends about 100Kbytes of data.
*/
extern double
CMprobe_bandwidth (CMConnection conn, long size, attr_list attrs);

/*!
 * Probe the available bandwidth on a particular connection by sending several streams
 * and do a linear regression.
 *
 * This is an intrusive probe.
 * \param conn The CMConnection to be tested.
 * \param size The size of message to be sent in the test.  (Bandwidth
 * varies dramatically with the message size.)
 * \param attrs Currently this parameter is ignored, but it *should* allow
 * control over the number of messages sent.
 * \return The return value is in units of KBytes per second.  
 * \note In the current implementation, CM sends \f$N\f$ messages to probe
 * available bandwidth, where \f$N\f$ is calculated as \f$100000/size\f$.
 * That is, CMprobe_bandwidth sends about 100Kbytes of data.
*/
extern double
CMregressive_probe_bandwidth (CMConnection conn, long size, attr_list attrs);

/*@}*/
/*!
 * Try to return the IP address of the current host as an integer.
 */
/*NOLOCK*/
extern int
CMget_self_ip_addr(CManager cm);

/*@}*/
/*!
 * Try to return the qualified name of the current host.
 */
/*NOLOCK*/
extern void
CMget_qual_hostname(CManager cm, char *buf, int len);

/*@}*/
/*!
 * Try to return the port range in use.
 */
/*NOLOCK*/
extern void
CMget_port_range(CManager cm, int *high_bound, int *low_bound);

/*!
 * Return a caller-owned string that give information about network interfaces, etc.
 */
/*NOLOCK*/
extern char *
CMget_ip_config_diagnostics(CManager cm);

/** @defgroup evpath EVPath functions and types
 * @{
 */
struct _EVStone;
struct _EVSource;
/*!
 * EVStone a stone is an elementary building block of paths
 *
 * EVStone is an integer-typed opaque handle.  Its only external use is 
 * to act as an external stone identifier for remote operations (such as 
 * specifying the remote target stone in a bridge action)
 */
typedef int EVstone;
/*!
 * EVaction actions, associated with stones, are the mechanisms through 
 * which data flow operations are defined.
 *
 * EVaction is an opaque integer-typed handle.  An EVaction handle is 
 * interpreted in the context of the stone it is associated with and is 
 * not unique across stones.
 */
typedef int EVaction;
/*!
 * EVsource an EVsource is a source handle used to submit events to EVpath.
 * An EVsource specifies both the (local) target stone and the format 
 * (fully-specified structured data type) of the data that will be submitted 
 * using this handle.  
 *
 * EVsource is an opaque handle.
 */
typedef struct _EVSource *EVsource;

/*!
 * The prototype for a EV submit callback function. 
 *
 * Used by EVsubmit_or_wait() 
 * \param cm The CManager with which this callback function was registered.
 * \param target The target stone that can now submit without stalling.
 * \param client_data This value is the same client_data value that was
 * supplied in the call.
 */
typedef void (*EVSubmitCallbackFunc) (CManager cm, EVstone target, 
					  void *client_data);

/*!
 * The prototype for an EVPath terminal handler function.
 *
 * EVPath allows application-routines matching this prototype to be 
 * registered as sinks on stones.
 * \param cm The CManager with which this handler was registered.
 * \param message A pointer to the incoming data, cast to void*.  The real
 * data is formatted to match the fields of with which the format was
 * registered. 
 * \param client_data This value is the same client_data value that was
 * supplied in the EVassoc_terminal_action() call.  It is not interpreted by CM,
 * but instead can be used to maintain some application context.
 * \param attrs The attributes (set of name/value pairs) that this message
 * was delivered with.  These are determined by the transport and may
 * include those specified in CMwrite_attr() when the data was written.
 */
typedef int (*EVSimpleHandlerFunc) (CManager cm, void *message, void *client_data,
				    attr_list attrs);

/*!
 * The prototype for an EVPath raw terminal handler function.
 *
 * EVPath allows application-routines matching this prototype to be 
 * registered as raw sinks (receiving encoded FFS data) on stones.
 * \param cm The CManager with which this handler was registered.
 * \param message A pointer to the incoming data, cast to void*.  The real
 * data is formatted to match the fields of with which the format was
 * registered. 
 * \param msg_len The length in bytes of the message block
 * \param client_data This value is the same client_data value that was
 * supplied in the EVassoc_terminal_action() call.  It is not interpreted by CM,
 * but instead can be used to maintain some application context.
 * \param attrs The attributes (set of name/value pairs) that this message
 * was delivered with.  These are determined by the transport and may
 * include those specified in CMwrite_attr() when the data was written.
 */
typedef int (*EVRawHandlerFunc) (CManager cm, void *message, size_t msg_len, void *client_data,
				 attr_list attrs);

/*!
 * The prototype for a EVPath bridge stone close handler.
 *
 * Functions matching of this prototype can be registered with
 * EVregister_close_handler(). 
 * \param cm The CManager with which this handler was registered.
 * \param conn The CMConnection which is being closed.
 * \param the stone ID of the bridge stone which has closed
 * \param client_data This value is the same client_data value that was
 * supplied in the EVregister_close_handler() call.  It is not interpreted 
 * by CM, but instead can be used to maintain some application context.
 */
typedef void (*EVStoneCloseHandlerFunc) (CManager cm, CMConnection conn, int stone, void *client_data);

struct _event_item;

/*!
 * Allocate a stone.
 *
 * Stones are the basic abstraction of EVPath, the entity to which events
 * are submitted and with which actions are associated.  The value returned
 * from EValloc_stone() is actually a simple integer which may be transmitted
 * to remote locations (for example for use in remote bridge actions).
 * \param cm The CManager which will manage the control for this stone.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls to associate actions with the stone.
 */
/*REMOTE*/
extern EVstone
EValloc_stone(CManager cm);

/*!
 * Free a stone.
 *
 * This call also free's all actions and data associated with a stone, 
 * including enqueued events if any.
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to free.
 */
/*REMOTE*/
extern void
EVfree_stone(CManager cm, EVstone stone);

/*!
 * Associate a terminal action (sink) with a stone.
 *
 * The specified handler will be called when data matching the 
 * format_list arrives at the stone.  The event data supplied may not 
 * remain valid after the handler call returns.  EVtake_event_buffer() may 
 * be used to ensure longer-term validity of the event data.  The 
 * parameters to the handler are those of EVSimpleHandlerFunc.
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param format_list The list of formats which describe the event data 
 * structure that the function accepts.
 * \param handler The handler function that will be called with data arrives.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
/*REMOTE*/
extern EVaction
EVassoc_terminal_action(CManager cm, EVstone stone, FMStructDescList format_list, 
			EVSimpleHandlerFunc handler, void* client_data);

/*!
 * Associate a raw terminal action (sink) with a stone.
 *
 * The specified handler will be called when any data.  Data is delivered in
 * FFS-encoded form using the EVRawHandlerFunc interface. The event data
 * supplied may not remain valid after the handler call returns.
 * EVtake_event_buffer() may be used to ensure longer-term validity of the
 * event data.  The parameters to the handler are those of
 * EVRawHandlerFunc.
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param handler The handler function that will be called with data arrives.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
extern EVaction
EVassoc_raw_terminal_action(CManager cm, EVstone stone, 
			    EVRawHandlerFunc handler, void* client_data);

/*!
 * Associate a terminal action (sink) with a new stone.
 *
 * The specified handler will be called when data matching the 
 * format_list arrives at the stone.  The event data supplied may not 
 * remain valid after the handler call returns.  EVtake_event_buffer() may 
 * be used to ensure longer-term validity of the event data.  The 
 * parameters to the handler are those of EVSimpleHandlerFunc.  This 
 * function differs from the previous function only in that it creates
 * a stone rather than using an existing stone.
 * \param cm The CManager from which this stone was allocated.
 * \param format_list The list of formats which describe the event data 
 * structure that the function accepts.
 * \param handler The handler function that will be called with data arrives.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 */
/*REMOTE*/
extern EVstone
EVcreate_terminal_action(CManager cm, FMStructDescList format_list, 
			EVSimpleHandlerFunc handler, void* client_data);

/*!
 * Associate a multiple-event action with a stone.
 *
 * EVassoc_multi_action() can be used to install handlers which
 * potentially consume multiple events as input and will not necessarily
 * consume any event immediately.  In evpath, these handlers are run
 * whenever a new event arrives in their queue, but that handler must
 * explicitly dequeue an event to consume it.  There is no declaration of
 * conditions under which a handler is to run.  It is expected to inspect
 * its queue upon activation and simply return without action if conditions
 * are not suitable for it to act (I.E. it needs more events, a different
 * set of events, etc.).  While the stone has only one real queue, events in
 * that queue are segregated into the types specified in the action_spec and
 * presented as being in multiple queues.
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param action_spec An action specification of the sort created by
 * create_multityped_action_spec().
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
/*REMOTE*/
extern EVaction
EVassoc_multi_action(CManager cm, EVstone stone, char *action_spec, 
		     void *client_data);

/*!
 * Associate an immediate non-terminal action with a stone.
 *
 * EVassoc_immediate_action() can be used to install handlers which
 * take only a single event as input and can therefore run and "consume"
 * their data immediately.  In particular, they are distinct from actions
 * which may leave their input data enqueued for some time (typically
 * handlers which might require more than one event to act).  The current
 * EVPath implementation supports only immediate actions with one input and
 * one output, but multiple output actions will be implemented soon.  
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param action_spec An action specification of the sort created by
 * create_filter_action_spec() or create_transform_action_spec().
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
/*REMOTE*/
extern EVaction
EVassoc_immediate_action(CManager cm, EVstone stone, char *action_spec, 
		      void *client_data);

/*!
 * Associate an immediate non-terminal action with a new stone.
 *
 * EVassoc_immediate_action() can be used to install handlers which
 * take only a single event as input and can therefore run and "consume"
 * their data immediately.  In particular, they are distinct from actions
 * which may leave their input data enqueued for some time (typically
 * handlers which might require more than one event to act).  The current
 * EVPath implementation supports only immediate actions with one input and
 * one output, but multiple output actions will be implemented soon.  This 
 * function differs from the previous function only in that it creates
 * a stone rather than using an existing stone.
 * \param cm The CManager from which this stone was allocated.
 * \param action_spec An action specification of the sort created by
 * create_filter_action_spec() or create_transform_action_spec().
 * \param target_list A -1 terminated list of stones to which outgoing
 * data is to be sent.  This initial list can be NULL (or merely have
 * an initial 0) to specify no targets at action initialization time.  
 * Values are filled in later with EVaction_set_output().
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 */
/*REMOTE*/
extern EVstone
EVcreate_immediate_action(CManager cm, char *action_spec, EVstone *target_list);

/*!
 * Direct the output of a stone action to another local target stone
 *
 * Immediate and queued actions have one or more outputs from which data
 * will emerge.  EVaction_set_output() is used to assign each of these
 * outputs to a local stone.  (It is NOT used with output stones.)
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which the action is registered.
 * \param action The action whose output is to be assigned.
 * \param output_index The zero-based index of the output to assign.
 * \param target_stone The stone to which the specified output should be
 * directed. 
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 *
 * \deprecated Early EVPath associated output ports with actions
 * individually, with the result that the topology of the stone graph was
 * potentially multiplex for each action, rather than being a property of
 * the stone.  This is now changed, and all actions share the same set of
 * output ports.  When this routine is called, the action parameter is
 * ignored, but new code should use EVstone_set_output().
 */
/*REMOTE*/
extern int
EVaction_set_output(CManager cm, EVstone stone, EVaction action, 
		    int output_index, EVstone target_stone);

/*!
 * Direct a particular output port of a stone to another local target stone
 *
 * Immediate and queued actions have one or more outputs from which data
 * will emerge.  EVaction_set_output() is used to assign each of these
 * outputs to a local stone.  (It is NOT used with output stones.)
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which the action is registered.
 * \param output_index The zero-based index of the output to assign.
 * \param target_stone The stone to which the specified output should be
 * directed. 
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
/*REMOTE*/
extern int
EVstone_set_output(CManager cm, EVstone stone, int output_index, EVstone target_stone);

/*!
 * Associate an immediate non-ECL filter action with a stone.
 *
 * EVassoc_filter_action() is similar to EVassoc_immediate_action() called
 * with an action spec generated by create_filter_action_spec(), except that
 * a function pointer is provided directly instead of having the function
 * generated by ECL.
 * 
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param incoming_format_list The list of formats which describe the event data 
 * structure that the function accepts.
 * \param handler The handler function that will be called with data arrives.
 * \param out_stone The local stone to which output should be directed.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 *
 * \deprecated  This function needs to go away and instead the functionality
 * should be integrated into a new create_*_action_spec() call that would
 * then be passed to EVassoc_immediate_action().
 */
/*REMOTE*/
extern EVaction
EVassoc_filter_action(CManager cm, EVstone stone, 
		      FMStructDescList incoming_format_list, 
		      EVSimpleHandlerFunc handler, EVstone out_stone,
		      void* client_data);

/*!
 * Associate a bridge action with a stone.
 *
 * Bridge actions perform network data transmission between address spaces.
 * EVassoc_bridge_action will acquire a CM-level connection to the remote
 * process specified by the \b contact_list parameter.  Data delivered to
 * the local stone specified by \b stone will be encoded, sent over the 
 * network link and delivered to \b remote_stone in the target address space.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The local stone to which to register the action.
 * \param contact_list A CM-level contact list (such as from
 * CMget_contact_list()) specifying the remote address space to connect to. 
 * \param remote_stone The stone ID in the remote address space to which
 * data is to be delivered.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 *
 * Bridge actions are associated with the default action of a stone and are
 * non-specific as far as input data, encoding and transmitting any event
 * presented to the action.  Bridge actions may not be modified after
 * association. 
 */
/*REMOTE*/
extern EVaction
EVassoc_bridge_action(CManager cm, EVstone stone, attr_list contact_list, 
		      EVstone remote_stone);

/*!
 * Associate a bridge action with a new stone.
 *
 * Bridge actions perform network data transmission between address spaces.
 * EVcreate_bridge_action will acquire a CM-level connection to the remote
 * process specified by the \b contact_list parameter.  Data delivered to
 * the local stone specified by \b stone will be encoded, sent over the 
 * network link and delivered to \b remote_stone in the target address space.
 * This function differs from the previous function only in that it creates
 * a stone rather than using an existing stone.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param contact_list A CM-level contact list (such as from
 * CMget_contact_list()) specifying the remote address space to connect to. 
 * \param remote_stone The stone ID in the remote address space to which
 * data is to be delivered.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 *
 * Bridge actions are associated with the default action of a stone and are
 * non-specific as far as input data, encoding and transmitting any event
 * presented to the action.  Bridge actions may not be modified after
 * association. 
 */
/*REMOTE*/
extern EVstone
EVcreate_bridge_action(CManager cm, attr_list contact_list, 
		       EVstone remote_stone);

/*!
 * Associate a thread bridge action with a stone.
 *
 * Thread bridge actions transfer events between CM control domains in the
 * same address space.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The local stone to which to register the action.
 * \param target_cm Another CManager in the same address space to which data should be transferred.
 * \param target_stone The stone ID associated with the target CM to which
 * data is to be delivered.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 *
 * Bridge actions are associated with the default action of a stone and are
 * non-specific as far as input data, generally transferring events without
 * encoding or copying.  Thread bridge actions may not be modified after
 * association. 
 */
extern EVaction
EVassoc_thread_bridge_action(CManager cm, EVstone stone, CManager target_cm,
			     EVstone target_stone);

/*!
 * Associate a thread bridge action with a stone.
 *
 * Thread bridge actions transfer events between CM control domains in the
 * same address space.
 *
 * \param cm The CManager from which a stone should be allocated.
 * \param target_cm Another CManager in the same address space to which data should be transferred.
 * \param target_stone The stone ID associated with the target CM to which
 * data is to be delivered.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 *
 * Bridge actions are associated with the default action of a stone and are
 * non-specific as far as input data, generally transferring events without
 * encoding or copying.  Thread bridge actions may not be modified after
 * association. 
 */
extern EVstone
EVcreate_thread_bridge_action(CManager cm, CManager target_cm,
			      EVstone target_stone);

/*!
 * Associate a split action with a stone.
 *
 * Split actions replicate an incoming event to multiple output target
 * stones.  All output paths receive every incoming event. (Reference counts
 * are updated, the event is not actually copied.)  Split actions may be
 * modified after association by using EVaction_add/remote_split_target() to
 * modify the target list.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The local stone to which to register the action.
 * \param target_list A '-1' terminated list of stones to which incoming
 * data is to be replicated.  This initial list can be NULL (or merely have
 * an initial '-1') to specify no targets at action initialization time.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
/*REMOTE*/
extern EVaction
EVassoc_split_action(CManager cm, EVstone stone, EVstone *target_list);

/*!
 * Associate a split action with a new stone.
 *
 * Split actions replicate an incoming event to multiple output target
 * stones.  All output paths receive every incoming event. (Reference counts
 * are updated, the event is not actually copied.)  Split actions may be
 * modified after association by using EVaction_add/remote_split_target() to
 * modify the target list.  This function differs from the previous function 
 * only in that it creates a stone rather than using an existing stone.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param target_list A '-1' terminated list of stones to which incoming
 * data is to be replicated.  This initial list can be NULL (or merely have
 * an initial -1) to specify no targets at action initialization time.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 */
/*REMOTE*/
extern EVstone
EVcreate_split_action(CManager cm, EVstone *target_list);

/*!
 * Add a target to a split action.
 *
 * This call adds a new target stone to the list of stones to which a split
 * action will replicate data.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The split stone.
 * \param action The split action ID (as returned by EVassoc_split_action()).
 * \param target_stone The target stone to add to the list.
 * \return Returns 1 on success, 0 on failure (fails if there is not a split
 * action on the specified stone).
 */
/*REMOTE*/
extern int
EVaction_add_split_target(CManager cm, EVstone stone, EVaction action,
			  EVstone target_stone);

/*!
 * Remove a target from a split action.
 *
 * This call removes a target stone from the list of stones to which a split
 * action will replicate data.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The split stone.
 * \param action The split action ID (as returned by EVassoc_split_action()).
 * \param target_stone The target stone to remove from the list.
 */
/*REMOTE*/
extern void
EVaction_remove_split_target(CManager cm, EVstone stone, EVaction action,
			  EVstone target_stone);

/*!
 * Add a particular target to any split action associated with a stone.
 *
 * This call adds a target stone to the list of stones to which a split
 * action will replicate data.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The split stone.
 * \param target_stone The target stone to add to the list.
 */
/*REMOTE*/
extern void
EVstone_add_split_target(CManager cm, EVstone stone, EVstone target_stone);

/*!
 * Remove a particular target from any split action associated with a stone.
 *
 * This call removes a target stone from the list of stones to which a split
 * action will replicate data.
 *
 * \param cm The CManager from which this stone was allocated.
 * \param stone The split stone.
 * \param target_stone The target stone to remove from the list.
 */
/*REMOTE*/
extern void
EVstone_remove_split_target(CManager cm, EVstone stone, EVstone target_stone);

/*!
 * Associate a congestion-event action with an output stone.
 *
 * EVassoc_congestion_action() can be used to install a handler which
 * should run when the events queued on an output stone cannot be
 * transmitted immediately (usually because of network congestion).
 * Like the multi_action, there is only real queue for output stones, but
 * for congestion actions events in that queue are segregated into the types
 * specified in the action_spec and presented as being in multiple queues.
 * \param cm The CManager from which this stone was allocated.
 * \param stone The stone to which to register the action.
 * \param action_spec An action specification of the sort created by
 * create_multityped_action_spec().
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return An action identifier, an integer EVaction value, which can be used
 * in subsequent calls to modify or remove the action.
 */
extern EVaction
EVassoc_congestion_action(CManager cm, EVstone stone, char *action_spec,
			  void* client_data);

/*!
 * Create a new storage stone.
 *
 * Store stone buffer up to a limit of events temporarily, and send
 * the event to out_stone when done.
 *
 * \param cm The CManager from which this stone is allocated.
 * \param out_stone The target stone for output
 * \param store_limit The maximum number of events to buffer
 *          (when this limit is reached events will be passed through)
 * \return The newly allocated stone.
 */
/*REMOTE*/
extern EVstone
EVcreate_store_action(CManager cm, EVstone out_stone, int store_limit);

/*!
 * Create a new stone with a general action.
 *
 * This is the non-specific action stone creation call
 *
 * \param cm The CManager from which this stone is allocated.
 * \param action_spec The spec returned from some create_*_action_spec() call
 * \return The newly allocated stone.
 */
/*REMOTE*/
extern EVstone
EVcreate_stone_action(CManager cm, char *action_spec);

/*!
 * Create a new storage action.
 *
 * Storage actions implement storage stones internally and will
 * accept any type of event.
 *
 * \param cm The CManager from which the stone is allocated
 * \param stone_num The stone to add the action to.
 * \param out_stone The stone buffered data will be sent to
 *          when done.
 * \param store_limit The maximum number of events to store. Set to -1
 *                      for no limit. Stone will ordinarily act like a
 *                      a buffer with this maximum size when the limit is
 *                      reached.
 * \return The number of the newly created action.
 */
/*REMOTE*/
extern EVaction 
EVassoc_store_action(CManager cm, EVstone stone_num, EVstone out_stone,
                        int store_limit);

/*!
 * Clear the contents stored in the specified storage action.
 * 
 * \param cm The CManager from which the stone is allocated
 * \param stone_num The stone the action is attached to
 * \param action_num The action created
 */
/*REMOTE*/ /* XXX??? */
extern void
EVclear_stored(CManager cm, EVstone stone_num, EVaction action_num);

/*!
 * Send the contents stored in the specified storage action.
 * The storage will be empty when this function returns.
 *
 * \param cm The CManager from which the stone is allocated
 * \param stone_num The stone the action is attached to
 * \param action_num The action created
 */
extern void
EVsend_stored(CManager cm, EVstone stone_num, EVaction action_num);

/*!
 * Count the number of items stored in a storage action.
 * 
 * \param cm The CManager from which the stone is allocated
 * \param stone_num The stone the action is attached to
 * \param action_num The action created (return value of EVassoc_store_action())
 * \return Number of events stored
 */
extern int 
EVstore_count(CManager cm, EVstone stone_num, EVaction action_num); 

/*!
 * Return whether we are sending from this storage stone.
 * 
 * \param cm The CManager from which this stone is allocated.
 * \param stone_num The stone the action is attached to
 * \param action_num The action created
 * \return True iff currently sending
 */
extern int
EVstore_is_sending(CManager cm, EVstone stone_num, EVaction action_num);

/*!
 * Start sending from a storage stone. Will not stop until done sending.
 *
 * \param cm The CManager from which this stone is allocated.
 * \param stone_num The stone the action is attached to
 * \param action_num The action created
 */
extern void
EVstore_start_send(CManager cm, EVstone stone_num, EVaction action_num); 

/*!
 * Set the maximum number of items stored in a storage stone,
 * when using it like a buffer. Excess items will be flushed to the
 * next stone in line.
 *
 * \param cm The CManager from which the stone is allocated
 * \param stone_num Which stone the action is attached to
 * \param action_num The storage action
 * \param store_limit The maximum number of events to keep buffered 
 */
extern void
EVset_store_limit(CManager cm, EVstone stone_num, EVaction action_num,
    int store_limit);

/*!
 * Create a submission handle (EVsource).
 *
 * EVpath is optimized for repetitive event streams.  Rather than specifying
 * the characteristics of data and the stone to which it is to be submitted
 * on every event submission, we use associate those characteristics with
 * EVsource handles.  These handles serve as a cache for internal information.
 *
 * \param cm The CManager associated with the stone.
 * \param stone The stone to which data is to be submitted.
 * \param data_format The FMStructDescList describing the representation of the
 * data. 
 * \return An EVsource handle for use in later EVsubmit() calls.
 */
extern EVsource
EVcreate_submit_handle(CManager cm, EVstone stone, FMStructDescList data_format);

/*!
 * Free a source.
 *
 * This call free's the resources associated with an EVsource handle..
 * \param source  The source to free.
 */
extern void
EVfree_source(EVsource source);

/*!
 * The prototype for a function which will free the memory associated with
 * an event.
 *
 * Normally, the EVpath event submission functions do not return until
 * it is safe for the application to destroy the submitted data (I.E. until
 * EVpath is finished with it).  However, if a "free" function is associated
 * with the event through the EVsource, EVpath will return sooner if there
 * is another thread of control available to prosecute the actions on the
 * event.  EVpath will then call the application-supplied free function to
 * free the event when the event data is no longer required.
 * Application-supplied event free functions must satisfy this profile. 
 * \param event_data  The address of the event data, expressed as a void*.
 * \param client_data The parameter is used to supply the free function with
 * the same client_data value that was specified in the
 * EVcreate_submit_handle_free() call.
 */
typedef void (*EVFreeFunction) (void *event_data, void *client_data);

/*!
 * Create a submission handle (EVsource), specifying a free function for the
 * event. 
 *
 * EVpath is optimized for repetitive event streams.  Rather than specifying
 * the characteristics of data and the stone to which it is to be submitted
 * on every event submission, we use associate those characteristics with
 * EVsource handles.  These handles serve as a cache for internal information.
 * This version of the call allows an EVFreeFunction to be associated with
 * the handle.  EVpath will take ownership of the submitted data, calling
 * the free function when processing is finished.  
 *
 * \param cm The CManager associated with the stone.
 * \param stone The stone to which data is to be submitted.
 * \param data_format The FMStructDescList describing the representation of the
 * data. 
 * \param free_func  The EVFreeFunction to call when EVPath has finished
 * processing the submitted data.
 * \param client_data The parameter is supplied to the free function and can
 * be used to supply it with additional information.
 * \return An EVsource handle for use in later EVsubmit() calls.
 */
extern EVsource
EVcreate_submit_handle_free(CManager cm, EVstone stone, FMStructDescList data_format,
			    EVFreeFunction free_func, void *client_data);

/*!
 * Submit an event for processing by EVPath.
 *
 * EVsubmit submits an event for processing by EVPath.  The format of the
 * submitted data must match the description given by the \b data_format
 * parameter when the EVsource handle was created.  The \b attrs parameter
 * specifies the attributes (name/value pairs) that the event is submitted
 * with.  These attributes will be delivered to the final terminal, as well
 * as being available at intermediate processing points.  Some attributes
 * may affect the processing or transmission of data, depending upon the
 * specific transport or processing agents.
 * \param source The EVsource handle through which data is to be submitted.
 * \param data The data to be submitted, represented as a void*.
 * \param attrs The attribute list to be submitted with the data.
 */
extern void
EVsubmit(EVsource source, void *data, attr_list attrs);

/*!
 * Submit an event for processing by EVPath.
 *
 * EVsubmit submits an event for processing by EVPath.  The format of the
 * submitted data must match the description given by the \b data_format
 * parameter when the EVsource handle was created.  The \b attrs parameter
 * specifies the attributes (name/value pairs) that the event is submitted
 * with.  These attributes will be delivered to the final terminal, as well
 * as being available at intermediate processing points.  Some attributes
 * may affect the processing or transmission of data, depending upon the
 * specific transport or processing agents.
 * \param source The EVsource handle through which data is to be submitted.
 * \param data The data to be submitted, represented as a void*.
 * \param free_func  The EVFreeFunction to call when EVPath has finished
 * processing the submitted data.
 * \param attrs The attribute list to be submitted with the data.
 *
 * \deprecated  This function is used to underly ECho, which allows the free
 * function to be specified with the submit.  New applications should
 * specify the free function in the submit handle.
 */
extern void
EVsubmit_general(EVsource source, void *data, EVFreeFunction free_func,
		 attr_list attrs);

/*!
 * Submit a pre-encoded event for processing by EVPath.
 *
 * EVsubmit submits a pre-encoded event for processing by EVPath.  The event 
 * must be a contiguous FFS-encoded block of data.  The \b attrs parameter
 * specifies the attributes (name/value pairs) that the event is submitted
 * with.  These attributes will be delivered to the final terminal, as well
 * as being available at intermediate processing points.  Some attributes
 * may affect the processing or transmission of data, depending upon the
 * specific transport or processing agents.
 * \param cm The CManager associated with the stone.
 * \param stone The stone to which data is to be submitted.
 * \param data The pre-encoded data to be submitted, represented as a void*.
 * \param data_len The length of the pre-encoded data block.
 * \param attrs The attribute list to be submitted with the data.
 *
 */
extern void
EVsubmit_encoded(CManager cm, EVstone stone, void *data, size_t data_len,
		 attr_list attrs);

/*!
 * Submit a pre-encoded event for processing by EVPath.
 *
 * EVtransfer_events will dequeue events in the incoming data queue 
 * associated with the \b src_stone parameter and re-enqueue them on
 * the dest_stone.
 * \param cm The CManager associated with the stone.
 * \param src_stone The stone from which events are to be removed
 * \param dest_stone The stone to which events are to transfered
 * \return -1 on error, otherwise count of events transferred on success;
 */
/*REMOTE*/
extern int
EVtransfer_events(CManager cm, EVstone src_stone, EVstone dest_stone);

/*!
 * Assume control over a incoming buffer of data.
 *
 * This call is designed to be used inside a EVSimpleHandlerFunc.  Normally
 * data buffers are recycled and EVPath only guarantees that the data
 * data delivered to an EVSimpleHandlerFunc will be valid for the duration
 * data of the call.  In that circumstance, a handler that wanted to
 * data preserve the data for longer than its own duration (to pass it to a
 * data thread or enter it into some other data structure for example) would
 * data have to copy the data.  To avoid that inefficiency, 
 * EVtake_event_buffer() allows the handler to take control of the 
 * buffer holding its incoming data.  The buffer will then not be recycled
 * until it is returned to CM with EVreturn_event_buffer().
 * \param cm The CManager in which the handler was called.
 * \param event The base address of the data (I.E. the message parameter to
 * the EVSimpleHandlerFunc).
 * \return 0 on error, 1 on success;
*/
extern int
EVtake_event_buffer (CManager cm, void *event);

/*!
 * Return a buffer of incoming data.
 *
 * This call recycles a data buffer that the application has taken control
 * of through EVtake_event_buffer().
 * \param cm The CManager in which the handler was called.
 * \param event The base address of the data (I.E. same value that was passed
 * to EVtake_event_buffer().
*/
extern void
EVreturn_event_buffer (CManager cm, void *event);

/*!
 * return the FFSDataHandle associated with an EVsource handle.
 *
 * Some middleware may find it useful to access the FMFormat that is
 * produced when the FMStructDescList associated with a source is registered
 * with FFS.  This call merely gives access to that information to save a
 * reregistration step.
 * \param source The EVsource value for which to retrieve the associated
 * FMFormat.
 */
extern FMFormat
EVget_src_ref_format(EVsource source);

/*!
 * Enable periodic auto-submits of NULL events on a stone.
 *
 * \param cm The CManager in which the stone is registered.
 * \param stone_num The stone which should receive auto-submits.
 * \param period_sec The period at which submits should occur, seconds portion.
 * \param period_usec The period at which submits should occur, microseconds
 * portion.
 */
/*REMOTE*/
extern void
EVenable_auto_stone(CManager cm, EVstone stone_num, int period_sec, 
		    int period_usec);

/*!
 * Enable periodic auto-submits of NULL events on a stone. This 
 * function differs from the previous function only in that it creates
 * a stone rather than using an existing stone.
 *
 * \param cm The CManager in which the stone is registered.
 * \param period_sec The period at which submits should occur, seconds portion.
 * \param period_usec The period at which submits should occur, microseconds
 * portion.
 * \param action_spec An action specification of the sort created by
 * create_filter_action_spec() or create_transform_action_spec().
 * \param out_stone The local stone to which output should be directed.
 * \return The stone identifier, an integer EVstone value, which can be used
 * in subsequent calls.
 */
/*REMOTE*/
extern EVstone
EVcreate_auto_stone(CManager cm, int period_sec, int period_usec, 
		    char *action_spec, EVstone out_stone);


/*!
 * Cause a stone to become "stalled" explicitly. In this state, the stone
 * will continue processing events as usual, but will propagate backpressure
 * as if it were overloaded.
 *
 * A stone marked as stalled with EVstall_stone will remain stalled
 * even as other sources of stalling (remote squelching, too many
 * queued but unprocessed events) change.
 *
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone to be marked as stalled
 */
/*REMOTE*/
extern void
EVstall_stone(CManager cm, EVstone stone_id);

/*! 
 * Undo EVstall_stone(), allowing the stone to become unstalled when other
 * reasons for stalling are not present.
 *
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone to unmark as stalled, which should have been marked
 *                 as stalled with EVstall_stone()
 */
/*REMOTE*/
extern void
EVunstall_stone(CManager cm, EVstone stone_id); 

/*!
 * If the stone pointed to by the source handle is not stalled, submit normally
 * and return true; otherwise, return false and call the supplied callback when it
 * is no longer stalled.
 * \param source The EVsource handle through which data is to be submitted.
 * \param data The data to be submitted, represented as a void*.
 * \param attrs The attribute list to be submitted with the data. 
 * \param cb The function to call if the submit cannot be performed now.
 * \param user_data Passed as a parameter to the callback.
 */
extern int
EVsubmit_or_wait(EVsource source, void *data, attr_list attrs, EVSubmitCallbackFunc cb, void *user_data);

/*!
 * As EVsubmit_or_wait, but as if calling EVsubmit_encoded. 
 * \param cm The CManager with which the target stone is registered
 * \param stone The target stone id
 * \param data The data to be submitted, represented as a void*.
 * \param data_len The length of the pre-encoded data.
 * \param attrs The attribute list to be submitted with the data. 
 * \param cb The function to call if the submit cannot be performed now.
 * \param user_data Passed as a parameter to the callback.
 */
extern int
EVsubmit_encoded_or_wait(CManager cm, EVstone stone, void *data, int data_len, attr_list attrs,
                            EVSubmitCallbackFunc cb, void *user_data);

/*!
 * Cause a stone to suspend operation
 *
 * This function causes a stone to enter a "suspended" state in which
 * incoming data will simply be queued, rather than submitted to any actions
 * which might be registered.  In the case of an output stone, will allow
 * the stone to finish the output action it is currently executing and then
 * prevent the output stone from sending any more data to the target stone.
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone which is to be frozen
 * \return Returns 1 on success, 0 on failure
 */ 
/*REMOTE*/
extern int
EVfreeze_stone(CManager cm, EVstone stone_id);

/*!
 * Cause a stone to resume operation
 *
 * This function causes a frozen stone (via EVfreeze_stone()) to resume
 * operation.  Pending data will be submitted to actions during the next
 * action processing phase. 
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone to unfreeze
 * \return Returns 1 on success, 0 on failure
 */ 
/*REMOTE*/
extern int
EVunfreeze_stone(CManager cm, EVstone stone_id);

/*!
 * Drain a stone
 *
 * This function is a blocking call that suspends the caller until all
 * events queued on a stone are processed (if processing is possible, it
 * might not be for events that require the presence of other events).
 * The function is typically used after upstream stones have been frozen
 * with EVfreeze_stone() during a reconfiguration action.  EVdrain_stone()
 * then makes sure a stone is as empty as possible prior to event extraction
 * and destruction.
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone which is to be drained
 * \return Returns 1 on success, 0 on failure
 */
/*REMOTE*/
extern int
EVdrain_stone(CManager cm, EVstone stone_id);

/*!
 * Return the EVstone id inside an executing terminal handler
 *
 * If called from within an executing terminal handler or a function called
 * by a terminal handler, this function will return the EVstone ID of the
 * stone with home that handler was registered.  Generally, the client_data
 * parameter in EVcreate/assoc_terminal_action() is the best way to pass
 * information to a handler, but in EVdfg a handler is registered once by
 * name and may be assigned to multiple stones.  This function provides a
 * way to differentiate between those stones during execution and to access
 * items such as the attribute list associated with those stones.
 *
 * \param cm The CManager in which the function is currently executing
 * \return The stone ID of the stone in whose handler we're executing (or -1 on error)
 */
extern EVstone
EVexecuting_stone(CManager cm);

/*!
 * Return the queued events associated with a stone and its actions.
 * 
 * This function will be called by EVdrain_stone. It will form an array of
 * structures where each structure will contain the size of the encoded 
 * event and a pointer to the encoded event. The array will contain an entry 
 * for each event, associated with the stone or its actions.   
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone whose associated events are to be extracted
 * \return  Returns an array of structures (EVevent_list) containing the
 * lengths of events and pointers to the encoded versions of events
 */
extern EVevent_list
EVextract_stone_events(CManager cm, EVstone stone_id);

/*!
 * Return the attribute list associated with a stone.
 *
 * This function is used to extract the set of attributes associated with a
 * stone.  It is normally used during a reconfiguration operation to
 * recreate a stone elsewhere. 
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone whose attributes are to be extracted
 * \return attr_list Returns the attribute list associated with the stone
 */
/*REMOTE*/
extern attr_list
EVextract_attr_list(CManager cm, EVstone stone_id);

/*!
 * Set the attribute list associated with a stone.
 *
 * This function is used to set the attributes associated with a
 * stone.  It is normally used during a reconfiguration operation to
 * recreate a stone elsewhere. 
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone whose attributes are to be extracted
 * \param list The attribute list to be associated with the stone
 */
/*REMOTE*/
extern void
EVset_attr_list(CManager cm, EVstone stone_id, attr_list list);

/*!
 * Free a stone after it has been drained.
 *
 * This function will wait till a stone is drained. Then it will free all the
 * data and events associated with the stone.
 * \param cm The CManager in which the stone is registered
 * \param stone_id The stone which is to be destroyed
 * \return Returns 1 on success, 0 on failure
 */ 
/*REMOTE*/
extern int
EVdestroy_stone(CManager cm, EVstone stone_id);

/*!
 * create an action specification for a filter function.
 *
 * 
 * \param format_list A description of the incoming event data that the
 * filter expects. 
 * \param function The filter function itself.  A zero return value means
 * that the data should be discarded. 
 */
/*NOLOCK*/
extern char *
create_filter_action_spec(FMStructDescList format_list, char *function);

/*!
 * create an action specification for a router function.
 *
 * 
 * \param format_list A description of the incoming event data that the
 * router function expects. 
 * \param function The router function itself.  A negative return value means
 * that the data should be discarded.  A positive value less than the number
 * of output values that have been set with EVaction_set_output() indicates
 * which of the output paths the input data should be submitted to.  Return
 * values larger than the number of output paths have undefined behavior.
 */
/*NOLOCK*/
extern char *
create_router_action_spec(FMStructDescList format_list, char *function);

/*!
 * create an action specification that transforms event data.
 *
 * \param format_list A description of the incoming event data that the
 * transformation expects. 
 * \param out_format_list A description of the outgoing event data that the
 * transformation will produce. 
 * \param function The processing that will perform the transformation.  A
 * zero return value means that the output data should be ignored/discarded.
 */
/*NOLOCK*/
extern char *
create_transform_action_spec(FMStructDescList format_list, FMStructDescList out_format_list, char *function);

/*!
 * create an action specification that operates on multiple queues of events
 *
 * \param input_format_lists A null-terminated list of null-terminated lists
 *  of descriptions of  the incoming event data types that the transformation
 *  expects. 
 * \param function The processing that will perform the transformation.  Not outputs are explicit, but must be performed by some EVsubmit* action.
 */
/*NOLOCK*/
extern char *
create_multityped_action_spec(FMStructDescList *input_format_lists, char *function);

#ifdef __COD__H__
/*!
 * Add a set of routines that will be visible in COD.
 *
 * \param cm The CManager in which the routines should be visible
 * \param extern_string A string that declares the routines, C-style 
 * \param externs A NULL-terminated structure of type cod_extern_entry.
 * This structure consists of name/address pairs that give transfer
 * addresses for each routine declared by the extern string.
 * You must include cod.h for this routine to be visible.
 */
extern void
EVadd_standard_routines(CManager cm, char *extern_string, 
			cod_extern_entry *externs);
#endif

/*!
 * Add a directory to search for DLL-based functions
 *
 * \param path_string The name of a directory to include for searches
 */
/*NOLOCK*/
extern void
EVadd_dll_search_dir(char *path_string);

/*!
 * Add a set of structure types that will be visible in COD.
 *
 * \param cm The CManager in which the routines should be visible
 * \param lists A NULL-terminated list of FMStructDescLists.
 */
extern void
EVadd_standard_structs(CManager cm, FMStructDescList *lists);

/*!
 * Register a handler to be called when a bridge stone is closed
 *
 * \param cm The CManager managing the bridge stones
 * \param handler The routine to be called
 * \param client_data This parameter will be supplied unmodified to the handler routine upon close.
 */
extern void
EVregister_close_handler(CManager cm, EVStoneCloseHandlerFunc handler, void *client_data);

/*!
 * Print a description of stone status to standard output.
 *
 * A simple dump function that can be used for debugging.
 * \param cm The CManager to which the stone is registered.
 * \param stone_num  The stone to dump.
 */
void
EVdump_stone(CManager cm,  EVstone stone_num);

/*!
 * The prototype of a specific immediate handler function.
 *
 * This function prototype is used by the EVPath internal "response"
 * interface.  At some point, the response interface will likely become
 * external so that EVPath's response to unknown data can be customized.
 * However, at the moment this is an internal interface.
 */
typedef int (*EVImmediateHandlerFunc) (CManager cm, struct _event_item *event, 
				       void *client_data, attr_list attrs, 
				       int out_count, int *out_stones);
/*!
 * Associate a conversion action.
 *
 * This function is used by the EVPath internal "response" interface.  At
 * some point, the response interface will likely become external so that
 * EVPath's response to unknown data can be customized.  However, at the
 * moment this is an internal interface.
 */
extern void
EVassoc_conversion_action(CManager cm, int stone_id, int stage, FMFormat target_format,
			  FMFormat incoming_format);

/*!
 * Perform a transport-level bandwidth test
 *
 * This function requests that EVPath perform a transport-level performance
 * check (currently of bandwidth) by sending a series of messages of a fixed
 * size, broken up internally into a fixed number of vectors (a la
 * writev()).  Th nature of the test is controlled by the 'how' attribute
 * list parameter.  In particular, the CM_TRANS_TEST_SIZE attribute controls
 * the message size in bytes, the CM_TRANS_TEST_VECS attribute specifies how
 * many vectors the message is to be broken up into (each vector will be
 * roughly size/vectors in bytes), and CM_TRANS_TEST_REPEAT gives the total
 * number of messages to send.  The total number of bytes sent will be size
 * * repeat_count.  The return attribute list value will contain at least
 * the following double-valued attributes, CM_TRANS_TEST_DURATION -
 * containing the duration of the test in seconds, and CM_TRANS_MEGABITS_SEC
 * - containing the estimated bandwidth in megabits per sec.  I.E. size *
 * repeat_count * 8 / 1000 * 1000 * secs.  Note the "bits", not bytes, and
 * that we use 1000 rather than 1024 to more closely match networking
 * conventions.  
 * 
 * \param conn The CMConnection on which to perform the bandwidth test 
 * \param how Parameters to control the bandwidth test
 */
extern attr_list
CMtest_transport(CMConnection conn, attr_list how);

/*!
 * Callback function to support transport tests
 * 
 * CMperf_upcall is the type of a callback function that is used as part of
 * CM transport testing.  This function is called on the receiving side of
 * the transport for every message that arrives.  The data of the message
 * and it's length are provided in parameters 'buffer' and 'length'
 * respectively.  The 'type' parameter varies with the message that arrives.
 * Current types are '0' for test initiation (I.E. the first short header of
 * a bandwidth tests), '1' for "body" messages (I.E. the N-1 repeats of a
 * bandwidth test), and '2' for the "final" message (I.E. the last message
 * of a bandwidth test).  The return value is only read for the final
 * message, and it is returned to caller of CMtest_transport(), so should
 * contain the final values of testing.
 *
 */
typedef attr_list (*CMperf_upcall)(CManager cm, void *buffer, size_t length, int type, attr_list list);

/*!
 * install the callback function to support transport testing
 *
 * CMinstall_perf_upcall() sets up an upcall to be used when
 * CMtest_transport() initiates a transport-level performance test.
 * \param cm The CManager in which to register the upcall.
 * \param upcall  The function to be called when a performance message arrives.
 */
extern void
CMinstall_perf_upcall(CManager cm, CMperf_upcall upcall);

/*!
 * Add a global stone ID to EVpath's internal lookup table
 *
 * CMadd_stone_to_global_lookup() inserts a "translation" of a global stone
 * ID (which must have the high bit set in it's 32-bit value) to a local
 * stone ID.  The mechanism for cleanly allocating global stone IDs is not
 * provided here, so this is an incomplete interface.  Exported for limited
 * use.
 * \param cm The CManager in which to register the translation.
 * \param stone_num the actual stone ID
 * \param global_stone_num the 32-bit ID with high bit set that will be used as a "well known" ID
 */
extern void
CMadd_stone_to_global_lookup(CManager cm, int stone_num, int global_stone_num);
/* @}*/

#ifdef	__cplusplus
}
#endif

#endif
