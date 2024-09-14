#ifndef __EV_DFG__H__
/*! \file */

#include "evpath.h"

#ifdef	__cplusplus
extern "C" {
#endif
/** @defgroup ev_dfg EVdfg functions and types
 * @{
 */
/*
**  Basic approach:
**  Create a DFG of "virtual stones" to be later deployed.
**  - Associate actions (created with our usual create action routines from StructDecls and text actions)
**  - Associate attributes
**  Actual deployment:
**  - first version:
**    * Centralized
**    * Server or distinguished participant is informed of all participants, he gives them DFG 
**      segments and gets reports back of stone IDs, necessary stone IDs are then distributed to 
**      neighbors and the cohort is released.
**    * require all participants to be present before continuing
**    * manage "ready messages" so that DFG is realized without race conditions
**    * Do "validation" of DFG pre-realization.  Check for disconnected vertices, look at data 
**      types to see if they will be handled, etc.
**    * deployment master may be distinguished participant or separate server.
**  - optimization
**    * eliminate actual stone IDs.  As long as the number of stones per CM is relatively small,
**      we can use DFG-assigned stone IDs in the event-message and do a table lookup to determine
**      the right local stone when the message arrives.
**  - next version
**    * add nodes to (fixed) DFG as they come up.  Provide early-arrivers with the contact lists of 
**      late-arriving neighbors as they join.
*/

/*!
 * EVdfg is a handle to a DFG.
 *
 * EVdfg is an opaque handle
 */
typedef struct _EVdfg *EVdfg;

/*!
 * EVmaster is a handle to a DFG master, the distinguished process who creates the DFGs and manages clients.
 *
 * EVmaster is an opaque handle
 */
typedef struct _EVmaster *EVmaster;

/*!
 * EVclient is a handle to a client process which can host DFG components.
 *
 * EVclient is an opaque handle
 */
typedef struct _EVclient *EVclient;

/*!
 * EVdfg_stone is a handle to virtual EVpath stone.
 *
 * EVdfg_stone is an opaque handle.  
 */
typedef struct _EVdfg_stone *EVdfg_stone;

/*!
 * EVclient_sources is a handle to the set of client event source capabilities.
 *
 * This opaque handle is returned by EVclient_register_source() and is then
 * specified to EVclient_assoc() or EVclient_assoc_local().  It represents
 * the capabilities of this client to host event source stones for any DFG
 * that might be created.
 */
typedef struct _EVclient_sources *EVclient_sources;

/*!
 * EVclient_sources is a handle to the set of client event handler capabilities.
 *
 * This opaque handle is returned by EVclient_register_sink_handler() or
 * EVclient_register_raw_sink_handler() and is then specified to
 * EVclient_assoc() or EVclient_assoc_local().  It represents the handlers
 * available on this client for use by event sink stones for any DFG that
 * might be created.
 */
typedef struct _EVclient_sinks *EVclient_sinks;

/*!
 * Special StoneID used in EVcreate_submit_handle()
 *
 * DFG_SOURCE is a special StoneID used in calls to EVcreate_submit_handle() to create 
 * EVsource handles that can be later assigned to DFG source nodes.
 */
#define DFG_SOURCE -1

/*!
 * Create a DFG master
 *
 * This call is used in master process to create a set of communicating
 * EVdfg processes.  The master is the unique point of control and
 * adminstration for client processes and DFGs.
 *
 * \param cm The CManager with which to associate the DFG
 * \return An EVmaster handle, to be used in later calls.
 */
extern EVmaster EVmaster_create(CManager cm);

/*!
 * Get the contact list from an EVmaster handle
 *
 * This call is used to extract contact information from an
 * EVmaster handle.  This call is made on the Master side of EVdfg
 * and the contact information is then provided to the remote Clients
 * for use in EVclient_assoc() calls.  The result of this call is a
 * null-terminated string to be owned by the caller.  (I.E. you should
 * free the string memory when you're done with it.)
 *
 * \param master The EVmaster handle for which to create contact information.
 * \return A null-terminated string representing contact information for this EVdfg
 */
extern char *EVmaster_get_contact_list(EVmaster master);

/*!
 *  Associate this process as a client to an EVmaster
 *
 *  This call is used to join the set of communicating processes as a
 * client.  The master_contact string should be the same one that came from
 * EVmaster_get_contact_list() on the master.  This call cannot be used by
 * the master process to participate in the DFG itself.  In that
 * circumstance, EVclient_assoc_local() should be used.
 *
 * The source_capabilities and sink_capabilities parameters identify the set
 * of sinks and sources that this client is capable of hosting.  Those
 * capabilities may or may not be utilized in the DFGs created by the
 * master.  All calls to register sinks and sources should be done *before*
 * the EVclient_assoc() call.
 *
 * \param cm The CManager with which to associate the DFG client
 * \param node_name The name with which the client can be identified.  This
 *  should be unique among the joining nodes in static joining mode
 *  (I.E. using EVmaster_register_node_list().  In dynamic mode (I.E. where
 *  EVmaster_node_join_handler() is used), then this name is presented to the
 *  registered join handler, but it need not be unique.  EVdfg copies this
 *  string, so it can be free'd after use.
 * \param master_contact The string contact information for the master
 *  process.  This is not stored by EVdfg.
 * \param source_capabilities The last value returned by
 *  EVclient_register_source() or NULL
 * \param sink_capabilities The last value returned by
 *  EVclient_register_sink_handle() or EVclient_register_raw_sink_handle, or NULL
 */
extern EVclient EVclient_assoc(CManager cm, char *node_name, char *master_contact, 
			       EVclient_sources source_capabilities, EVclient_sinks sink_capabilities);

/*!
 *  Associate a client to an EVmaster in the same process
 *
 *  This call is for the process which hosts the EVdfg master to also
 *  participate as an EVclient.
 *
 * The source_capabilities and sink_capabilities parameters identify the set
 * of sinks and sources that this client is capable of hosting.  Those
 * capabilities may or may not be utilized in the DFGs created by the
 * master.  All calls to register sinks and sources should be done *before*
 * the EVclient_assoc() call.
 *
 * \param cm The CManager with which to associate the DFG client
 * \param node_name The name with which the client can be identified.  This
 *  should be unique among the joining nodes in static joining mode
 *  (I.E. using EVmaster_register_node_list().  In dynamic mode (I.E. where
 *  EVmaster_node_join_handler() is used), then this name is presented to the
 *  registered join handler, but it need not be unique.  EVdfg copies this
 *  string, so it can be free'd after use.
 * \param master The handle for the EVmaster.
 * \param source_capabilities The last value returned by
 *  EVclient_register_source() or NULL
 * \param sink_capabilities The last value returned by
 *  EVclient_register_sink_handle() or EVclient_register_raw_sink_handle, or NULL
 * \return a handle to EVclient upon success
 */
extern EVclient EVclient_assoc_local(CManager cm, char *node_name, EVmaster master,
				     EVclient_sources source_capabilities, EVclient_sinks sink_capabilities);


/*!
 *  Supply a static list of client names to the EVmaster.
 *
 * This call is used in static joining mode, that is when the set of nodes
 * that will join the DFG client set is known upfront and each has a predefined unique 
 * name. 
 * \param master The EVmaster handle for which the set of clients is to be
 * specified. 
 * \param list A NULL-terminated list of NULL-terminated strings.  These
 * names must be unique, and each must be used in an EVclient_assoc() or
 * EVclient_assoc_local() call before the ensemble will be complete.  EVdfg
 * copies this list, so it should be free'd by the calling application if
 * dynamic.
 *
 */
extern void EVmaster_register_node_list(EVmaster master, char** list);

/*!
 * The prototype for an EVmaster client-join handling function.
 *
 * In "dynamic join" mode, as opposed to static-client-list mode (See
 * EVmaster_register_node_list()), the EVdfg master calls a registered join
 * handler each time a new client node joins the DFG.  This handler should
 * 1) potentially assign a canonical name to the client node (using
 * EVmaster_assign_canonical_name()), and 2) when the expected nodes have
 * all joined the handler should create the virtual DFG and then call
 * EVdfg_realize() in order to instantiate it.
 * 
 * This call happens in the Master process only.
 * 
 * \param dfg The EVdfg handle with which this handler is associated
 * \param identifier This null-terminated string is the value specified in
 * the corresponding EVclient_assoc() call by a client.
 * \param available_sources This parameter is currently a placeholder for
 * information about what sources the client is capable of hosting.
 * \param available_sinks This parameter is currently a placeholder for
 * information about what sinks the client is capable of hosting.
 */
typedef void (*EVmasterJoinHandlerFunc) (EVmaster master, char *identifier, void* available_sources, void *available_sinks);

/*!
 * The prototype for an EVdfg client-fail handling function.
 *
 * If the an EVmasterFailHandlerFunc has been registered to the DFG, this
 * function will be called in the event that some node has failed.
 * Generally EVdfg becomes aware that a node has failed when some other
 * client tries to do a write to a bridge stone connected to that node and
 * that write fails.  Three things to be cognizant of:
 *  - This call happens in the Master process only.
 *  - If it's the Master that has failed, you're not going to get a call.
 *  - If multiple nodes notice the same failure, you might get multiple
 * calls to this function resulting from a single failure.
 * 
 * \param dfg The EVdfg handle with which this handler is associated
 * \param failed_client This null-terminated string is the value that was
 * specified in EVclient_assoc() in the failed client.
 * \param failed_stone This is local ID of the failed stone (perhaps not
 * useful to anyone - should change or eliminate this).
 */
typedef void (*EVmasterFailHandlerFunc) (EVdfg dfg, char *failed_client, int failed_stone);

/*!
 * The prototype for an EVdfg voluntary-reconfiguration handling function.
 *
 * If the an EVmasterReconfigHandlerFunc has been registered to the DFG, this
 * function will be called in the event that some stone has executed a call
 * to EVdfg_trigger_reconfiguration() inside a CoD-based event handler.
 * This allows a client application to monitor local conditions and to
 * trigger action by the master in the event a reconfiguration is called
 * for.  The EVdfg_trigger_reconfiguration() call and the
 * EVdfg_flush_attrs() call both cause local stone attributes to be flushed
 * back to the Master's virtual stones so that they can be examined using
 * EVdfg_get_attr_list().
 * 
 * \param dfg The EVdfg handle with which this handler is associated
 */
typedef void (*EVmasterReconfigHandlerFunc) (EVdfg dfg);

/*!
 * Register a node join handler function to an EVmaster.
 *
 * \param master The EVmaster handle with which to associate this handler
 * \param func The handler function to associate
 */
extern void EVmaster_node_join_handler (EVmaster master, EVmasterJoinHandlerFunc func);

/*!
 * Register a node fail handler function to an EVmaster
 *
 * \param master The EVmaster handle with which to associate this handler
 * \param func The handler function to associate
 */
extern void EVmaster_node_fail_handler (EVmaster master, EVmasterFailHandlerFunc func);

/*!
 * Register a voluntary reconfiguration handler function to an EVmaster
 *
 * \param master The EVmaster handle with which to associate this handler
 * \param func The handler function to associate
 */
extern void EVmaster_node_reconfig_handler (EVmaster master, EVmasterReconfigHandlerFunc func);

/*!
 * Cause the instantiation of a virtual DFG.
 *
 *  This call is performed by the master to signal the end of the creation
 *  or reorganization of a DFG.  In static-client-list mode, it is generally
 *  called by master after the virtual DFG has been created.  In
 *  dynamic-join mode, creating the virtual DFG and calling EVdfg_realize()
 *  is how the join handler signals EVdfg that all expected nodes have
 *  joined.  In a node fail handler, a node reconfig handler, or when the
 *  join handler is called after the first realization of the DFG, a further
 *  call of EVdfg_realize() represents the finalization of a reconfiguration
 *  of the DFG.
 *
 * \param dfg The handle of the EVdfg to be realized.
 */
extern void EVdfg_realize(EVdfg dfg);

/*!
 * Wait for a DFG to be ready to run
 *
 *  This call is performed by participating clients (including the master
 *  process operating as a client), in order to wait for the deployment of a
 *  DFG.  This call will only return after a DFG has been completely
 *  instantiated and is running.  In initial deployment, all participating
 *  clients will exit this call at roughly the same time.  In the case of
 *  client nodes that join after a DFG is already running, this call will
 *  return after the master has been notified of the join and the client has
 *  been incorporated into the DFG (I.E. stones potentially deployed.)
 *
 * \param client The EVclient handle upon which to wait.
 * \return 1 on success, 0 on failure
 */
extern int EVclient_ready_wait(EVclient client);

/*!
 * Associate a name with a source handle
 *
 *  This call is performed by client nodes in order to register their
 *  ability to host a source with the given name.  If the name appears in
 *  a virtual source stone that is mapped to the client, then that source
 *  handle is `active`. 
 *
 * \param name The name to associate with the source handle.  EVdfg does not
 *  take ownership of this string and the application should ensure that the
 *  string remains valid for the duration of EVdfg operation.
 * \param src The EVsource to be associated with the name.  Source/name
 *  association is actually an EVPath-level operation, so there is no EVdfg
 *  parameter in this call.
 * \return The last value returned by EVclient_register_source() should be passed 
 *  to EVclient_assoc() or EVclient_assoc_local().  Prior return values can be ignored.
 */
extern EVclient_sources
EVclient_register_source(char *name, EVsource src);

/*!
 * Associate a name with a sink handle
 *
 *  This call is performed by client nodes in order to register their
 *  ability to host a sink stone with a handler with the given name.  
 *
 * \param cm The client's CManager
 * \param name The name to associate with the sink handle.  EVdfg does not
 *  take ownership of this string and the application should ensure that the
 *  string remains valid for the duration of EVdfg operation.
 * \param list The FMStructDescList describing the data that this handler
 *  function expects.  EVdfg does not take ownership of this data structure
 *  and the application should ensure that the structure remains valid for
 *  the duration of EVdfg operation.
 * \param handler The EVSimpleHandlerFunc to be associated with the
 *  name/data type.  Sink-handle/name association is actually an
 *  EVPath-level operation, so there is a CM parameter in this call, but no
 *  EVdfg param.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return The last value returned by EVclient_register_sink_handler() or 
 *  EVclient_register_raw_sink_handler(), should be passed to EVclient_assoc() or 
 *  EVclient_assoc_local().  Prior return values can be ignored.
 */
extern EVclient_sinks
EVclient_register_sink_handler(CManager cm, char *name, FMStructDescList list, EVSimpleHandlerFunc handler, void* client_data);

/*!
 * Associate a name with a raw sink handle
 *
 *  This call is performed by client nodes in order to register their
 *  ability to host a sink stone with a raw handler with the given name.  
 *
 * \param cm The client's CManager
 * \param name The name to associate with the sink handle.  EVdfg does not
 *  take ownership of this string and the application should ensure that the
 *  string remains valid for the duration of EVdfg operation.
 * \param handler The EVRawHandlerFunc to be associated with the
 *  name/data type.  Sink-handle/name association is actually an
 *  EVPath-level operation, so there is a CM parameter in this call, but no
 *  EVdfg param.
 * \param client_data An uninterpreted value that is passed to the handler
 * function when it is called.
 * \return The last value returned by EVclient_register_sink_handler() or 
 *  EVclient_register_raw_sink_handler(), should be passed to EVclient_assoc() or 
 *  EVclient_assoc_local().  Prior return values can be ignored.
 */
extern EVclient_sinks
EVclient_register_raw_sink_handler(CManager cm, char *name, EVRawHandlerFunc handler, void *client_data);

/*!
 *  return a boolean describing whether the source has been assigned in EVdfg
 *
 *  This call is performed by client nodes in order to determine if a
 *  particular local EVsource registered with EVdfg has actually been made
 *  active by having a virtual source stone associated with it.
 *
 * \param src The source to test
 * \return true if the source has been assigned, false otherwise
 */
extern int EVclient_source_active(EVsource src);

/*!
 *  return a count of EVdfg sink stones assigned to the current client
 *
 *  This call is performed by client nodes in order to determine how many
 *  sink stones have been assigned to them.
 *
 * \param client The EVclient under consideration.
 * \return the number of active sink stones assigned to the current client
 */
extern unsigned int EVclient_active_sink_count(EVclient client);

/*!
 *  Assign a unique, canonical name to a client of a particular given_name.
 *
 *  This call is performed by the master, typically in the
 *  EVmasterJoinHandlerFunc, in order to assign a unique name to clients who may
 *  not have one previously.  The canonical name is the name to be used in
 *  EVdfg_assign_node().  The only constraint on the name is that not have been 
 *  assigned to any already-participating node.
 *
 * \param master The master the client is joining.
 * \param given_name The original name of the client.
 * \param canonical_name The canonical name to be assigned to the client.
 * \return true on success, false if the name was not unique
 */
extern int EVmaster_assign_canonical_name(EVmaster master, char *given_name, char *canonical_name);

/*!
 * Create a DFG
 *
 * This call is used in the master process to create a DFG that will be
 * deployed on the client processes.
 * \param master The EVmaster with which to associate the DFG
 * \return An EVdfg handle, to be used in later calls.
 */
extern EVdfg EVdfg_create(EVmaster master);


/*!
 *  Create an EVdfg stone with a specific action associated with it.
 *
 *  This call is performed by the master during the Stone Creation and
 *  Assignment process in order to create an EVdfg_stone.
 *
 * \param dfg The DFG under consideration.
 * \param action_spec An action specifier string such as is created by the
 *  EVpath create_*_action_spec() routines.  This parameter may be NULL if
 *  no action is to be initially assigned to the stone.  If non-NULL, EVdfg
 *  takes ownership of the action_spec string and will free it upon DFG
 *  termination. 
 * \return Function returns an EVdfg_stone handle that can be used in
 *  subsequent calls.
 */
extern EVdfg_stone EVdfg_create_stone(EVdfg dfg, char *action_spec);

/*!
 *  Add an action to an existing EVdfg stone.
 *
 *  This call is performed by the master during the Stone Creation and
 *  Assignment process in order to add an action to an EVdfg_stone. 
 *
 * \param stone The EVdfg_stone to which to add the action.
 * \param action_spec An action specifier string such as is created by the
 *  EVpath create_*_action_spec() routines.  EVdfg takes ownership of the
 *  action_spec string and will free it upon DFG termination. 
 */
extern void EVdfg_add_action (EVdfg_stone stone, char *action_spec);

/*!
 *  Create an EVdfg stone that will act as an Event source.
 *
 *  This call is performed by the master during the Stone Creation and
 *  Assignment process in order to create an EVdfg source stone.  No other
 *  actions should be assigned to an EVdfg source stone.
 *
 * \param dfg The DFG under consideration.
 * \param source_name This value must match some value which has been
 *  specified in EVclient_register_source() on the node to which this stone is
 *  eventually mapped.  (EVdfg can't detect a mismatch until EVdfg_realize()
 *  is called.)  The source_name string is *not* owned by EVdfg and should
 *  be free'd by the application if dynamic.
 * \return Function returns an EVdfg_stone handle that can be used in
 *  subsequent calls.
 */
extern EVdfg_stone EVdfg_create_source_stone(EVdfg dfg, char *source_name);

/*!
 *  Create an EVdfg stone that will act as an Event sink (terminal stone).
 *
 *  This call is performed by the master during the Stone Creation and
 *  Assignment process in order to create an EVdfg sink stone.
 *
 * \param dfg The DFG under consideration.
 * \param handler_name This value must match some value which has been
 *  specified in EVclient_register_sink_handler() on the node to which this
 *  stone is eventually mapped.  (EVdfg can't detect a mismatch until
 *  EVdfg_realize() is called.)  The handler_name string is *not* owned by
 *  EVdfg and should be free'd by the application if dynamic.
 * \return Function returns an EVdfg_stone handle that can be used in
 *  subsequent calls.
 */
extern EVdfg_stone EVdfg_create_sink_stone(EVdfg dfg, char *handler_name);

/*!
 *  Add a sink action to an existing EVdfg stone.
 *
 *  This call is performed by the master during the Stone Creation and
 *  Assignment process in order to add a sink (terminal) action to an EVdfg
 *  stone. 
 *
 * \param stone The EVdfg_stone to which to add the sink action.
 * \param handler_name This value must match some value which has been
 *  specified in EVclient_register_sink_handler() on the node to which this
 *  stone is eventually mapped.  (EVdfg can't detect a mismatch until
 *  EVdfg_realize() is called.)  The handler_name string is *not* owned by
 *  EVdfg and should be free'd by the application if dynamic.
 */
extern void EVdfg_add_sink_action(EVdfg_stone stone, char *handler_name);

/*!
 * Link a particular output port of one stone (the source) to a destination
 * (target) stone.
 * 
 * This function is roughly the analog of the EVstone_set_output function,
 * but at the EVdfg level.  All non-terminal stones have one or more output
 * ports from which data will emerge.  EVdfg_link_port() is used to assign
 * each of these outputs to another EVdfg_stone stone.  For EVPath actions
 * which have a single output, such as 'filter' and 'transform', those
 * outputs appear on 'port_index' 0.  Other actions, such as router and
 * multityped actions allow submission to any port, so care must be taken
 * that the code in those actions uses appropriate port numbers as assigned
 * here.  The complement of EVdfg_link_port() is EVdfg_unlink_port().  EVdfg
 * source stones, and stones with no assigned actions, act as EVPath 'split'
 * stones and their outputs are best managed with EVdfg_link_dest() and
 * EVdfg_unlink_dest().
 *
 * \param source The EVdfg_stone whose ports are to be assigned.
 * \param output_index The zero-based index of the output which should be assigned.
 * \param destination The EVdfg_stone which is to receive those events.
 */
extern void EVdfg_link_port(EVdfg_stone source, int output_index, 
			    EVdfg_stone destination);

/*!
 * Link an anonymous output port of one stone (the source) to a destination
 * (target) stone.
 * 
 * This function is roughly the analog of the EVstone_add_split_target, but
 * at the EVdfg level.  EVdfg stones with no additional actions associated
 * with them are implicitly EVPath split stones and events presented to
 * their inputs are replicated to all outputs.  Those outputs can be
 * maintained by their explicit numbers (with EVdfg_link_port()), or, since
 * the specific output port to which a target stone is assigned is
 * unimportant, the EVdfg_link_dest() call will simply assign the target
 * stone to the next available port.  The complement to this call is
 * EVdfg_unlink_dest().
 *
 * \param source The EVdfg_stone whose ports are to be assigned.
 * \param destination The EVdfg_stone which is to receive those events.
 */
extern void EVdfg_link_dest(EVdfg_stone source, EVdfg_stone destination);

/*!
 * Unlink a particular output port of one stone (the source).
 * 
 * The output port should have been previously set with EVdfg_link_port.
 *
 * \param source The EVdfg_stone whose ports are to be assigned.
 * \param output_index The zero-based index of the output which should be assigned.
 * \return returns true on success, false if the particular port was not previously set
 */
extern int EVdfg_unlink_port(EVdfg_stone source, int output_index);

/*!
 * Unlink the output port of one stone (the source) so that it no longer
 * directs events to a particular destination (target) stone.
 * 
 * \param source The EVdfg_stone whose ports are to be assigned.
 * \param destination The EVdfg_stone which should no longer receive those events.
 * \return returns true on success, false if the destination stone was not
 * found on the output links of the source
 */
extern int EVdfg_unlink_dest(EVdfg_stone source, EVdfg_stone destination);

/*!
 * Assign a particular EVdfg_stone to a particular client node.
 * 
 * This function assigns a particular EVdfg_stone to be instantiated upon a
 *  particular client node.  The node parameter must match either 1) a
 *  string specified in EVmaster_register_node_list() (in static joining mode),
 *  or 2) a canonical name that has been assigned to a node (in dynamic
 *  joining mode).
 *
 * \param stone The EVdfg_stone to be assigned to a particular node.
 * \param node The node to which it is to be assigned.  EVdfg does not take
 *  ownership of this string and it should be free'd by the application if dynamic.
 * \return returns true on success, false if the node name given was not found in the set of nodes.
 */
extern int EVdfg_assign_node(EVdfg_stone stone, char *node);

/*!
 * Set virtual stone properties to enable periodic auto-submits of NULL
 * events on an EVdfg_stone.  This becomes active upon DFG realization and
 * activation.
 * 
 * This function is the analog of the EVenable_auto_stone() function, but at
 * the EVdfg level.
 *
 * \param stone The EVdfg_stone which should receive auto-submits.
 * \param period_sec The period at which submits should occur, seconds portion.
 * \param period_usec The period at which submits should occur, microseconds
 * portion.
 *  
 * Autosubmits are initiated on each node just as it is about to return from
 * EVclient_ready_wait(), I.E. with DFG activation.  
 */
extern void EVdfg_enable_auto_stone(EVdfg_stone stone, int period_sec, 
				    int period_usec);

/*!
 * Set the attribute list that will be visible as "stone_attrs" inside EVdfg
 *  dynamic functions.
 * 
 * \param stone The EVdfg_stone affected.
 * \param attrs The attribute list to be assigned.  EVdfg does an
 *  add_ref_attr_list() on this list.  Because lists are only free'd when
 *  the reference count goes to zero, the application should generally call
 *  free_attr_list() on the list as well.
 */
extern void EVdfg_set_attr_list(EVdfg_stone stone, attr_list attrs);

/*!
 * Query the attribute list that is/was visible as "stone_attrs" inside EVdfg
 *  dynamic functions.
 * 
 * \param stone The EVdfg_stone involved.
 * \return The "stone_attrs" attribute list.  EVdfg does an
 *  add_ref_attr_list() on this list before returning it.  Because lists are
 *  only free'd when the reference count goes to zero, the application
 *  should generally call free_attr_list() when it is finished with it.  
 *  Unless dynamic reconfiguration is in play, the stone_attrs value
 *  reported here does not reflect updates from the instantiated stones on
 *  the client nodes (even if local to the master).  However, for voluntary
 *  reconfiguration, instance stone attributes are flushed to the master and
 *  can be interrogated with this call.
 */
extern attr_list EVdfg_get_attr_list(EVdfg_stone stone);

/*!
 *  a value used in EVclient shutdown calls
 *
 * DFG_STATUS_SUCCESS is a value used in EVclient_shutdown* calls
 */
#define DFG_STATUS_SUCCESS 0

/*!
 *  a value used in EVclient shutdown calls
 *
 * DFG_STATUS_FAILURE is a value used in EVclient_shutdown* calls
 */
#define DFG_STATUS_FAILURE 1

/*!
 *  Vote that this node is ready for shutdown, provide it's contribution
 *  to the return value and wait for the actual shutdown to occur.
 *
 *  One of EVclient_shutdown() or EVclient_ready_for_shutdown() should be called
 *  by every participating node for normal shutdown.  The return value from
 *  all calls (on all nodes) to EVclient_shutdown() and
 *  EVclient_wait_for_shutdown() will all be the same.
 *
 * \param client The EVclient handle to which a return value is contributed
 * \param result this node's contribution to the DFG-wide shutdown value.
 *  In order to facilitate this value being passed directly to exit(), the
 *  convention for DFG exit values follows the Unix process exit semantics
 *  where 0 indicates success and non-zero indicates failure.  In
 *  particular, the value provided should be either DFG_STATUS_SUCCESS (0)
 *  or DFG_STATUS_FAILURE (1)
 * \return the DFG-wide shutdown value
 */
extern int EVclient_shutdown(EVclient client, int result);

/*!
 *  Override voting and force shutdown with specific return value.
 *
 *  One of EVclient_shutdown() or EVclient_ready_for_shutdown() should be called
 *  by every participating node for normal shutdown.  However, this call
 *  forces an abnormal shutdown. The return value from all calls (on all
 *  nodes) to EVclient_shutdown() and EVclient_wait_for_shutdown() will all the
 *  value of this result parameter.
 *
 * \param client The EVclient handle for which shutdown should be forced
 * \param result this node's contribution to the DFG-wide shutdown value.
 *  In order to facilitate this value being passed directly to exit(), the
 *  convention for DFG exit values follows the Unix process exit semantics
 *  where 0 indicates success and non-zero indicates failure.  In
 *  particular, the value provided should be either DFG_STATUS_SUCCESS (0)
 *  or DFG_STATUS_FAILURE (1)
 * \return the DFG-wide shutdown value (actually equal to result here)
 *
 */
extern int EVclient_force_shutdown(EVclient client, int result);

/*!
 *  Vote that this node is ready for shutdown and is providing no specific contribution
 *  to the return value.
 *
 *  One of EVclient_shutdown() or EVclient_ready_for_shutdown() should be called
 *  by every participating node for normal shutdown.  The return value from
 *  EVclient_wait_for_shutdown() will be the same on each node and will depend
 *  upon every node's contribution to the shutdown status.
 *
 * \param client The EVclient handle for which shutdown is indicated.
 *
 */
extern void EVclient_ready_for_shutdown(EVclient client);

/*!
 *  Wait for EVdfg to determine that the coordinated shutdown time has arrived.
 *
 *  One of EVclient_shutdown() or EVclient_ready_for_shutdown() should be called
 *  by every participating node for normal shutdown.  The return value from
 *  EVclient_wait_for_shutdown() will be the same on each node and will depend
 *  upon every node's contribution to the shutdown status.
 *
 * \param client The EVclient handle for which shutdown is indicated.
 * \return the DFG-wide shutdown value.  In order to facilitate this value
 * being passed directly to exit(), the convention for DFG exit values
 * follows the Unix process exit semantics where 0 indicates success and
 * non-zero indicates failure.  In particular, the value returned here should be
 * either DFG_STATUS_SUCCESS (0) or DFG_STATUS_FAILURE (1)
 */
extern int EVclient_wait_for_shutdown(EVclient client);

/*!
 *  test to see if the coordinated shutdown time has arrived.
 *
 *  This call essentially tests to see if EVclient_wait_for_shutdown() will
 * return immediately.  If false, EVclient_wait_for_shutdown() would block
 * as the coordinated shutdown time has not yet arrived.  Clients using this
 * call must still use EVclient_wait_for_shutdown() to exit cleanly.
 *
 * \param client The EVclient handle for which shutdown may be imminent
 * \return boolean True if the DFG-wide shutdown time has arrived
 */
extern int EVclient_test_for_shutdown(EVclient client);

typedef enum _EVdfg_state_type {
    EVdfgWorking,
    EVdfgDeployed
} EVdfg_state_type;

/*!
 *  Dump a graphml version of the stone graph
 *
 * As a means of discerning whether the programmer has created the stone
 * graph that is intended, this call provides a mechanism to print (to
 * stdout) a GraphML representation of the stones and links as graph
 * vertices and directed edges.  At present, the output is best viewed in
 * the yEd graph editor.  Copy the output from the screen, paste to a file,
 * open in yEd and choose a layout scheme (try auto-layout).
 *
 * Currently we set a node label (specific to yEd) that is just S%d based on
 * the stone id.  Links also have the appropriate output port property set,
 * though this does not affect their graphical display.
 *
 * \param which Set to EVdfgDeployed to display a realized graph, or
 * EVdfgWorking to display an unrealized or in progress reconfiguration
 * \param dfg The dfg to display
 */
extern void
EVdfg_dump_graph(EVdfg_state_type which, EVdfg dfg);

#ifdef	__cplusplus
}
#endif

#define __EV_DFG__H__
#endif
