#include "config.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <atl.h>
#include <evpath.h>
#include "chr_time.h"
#if !NO_DYNAMIC_LINKING
#include "dlloader.h"
#endif
#undef NDEBUG
#include "assert.h"
#include <cm_internal.h>
#include <cm_transport.h>

extern struct CMtrans_services_s CMstatic_trans_svcs;
/* const lt_dlsymlist lt_preloaded_symbols[1] = { { 0, 0 } }; */

static transport_entry *global_transports = NULL;

/* static void */
/* free_global_transports() */
/* { */
/*     int i = 0; */
/*     transport_entry *tmp; */
/*     if (global_transports == NULL) { */
/* 	return; */
/*     } else { */
/* 	tmp = global_transports; */
/* 	global_transports = NULL; */
/*     } */
/*     while (tmp[i]) { */
/* #if !NO_DYNAMIC_LINKING */
/* 	CMdlclose(tmp[i]->dlhandle); */
/* #endif */
/* 	free(tmp[i]->trans_name); */
/* 	free(tmp[i]); */
/* 	i++; */
/*     } */
/*     free(tmp); */
/* #if !NO_DYNAMIC_LINKING */
/*     CMdlclearsearchlist(); */
/* #endif */
/* } */

int
find_transport_in_cm(CManager cm, const char *trans_name)
{
    int i = 0;
    if (cm->transports == NULL)
	return 0;
    while (cm->transports[i] != NULL) {
	if (strcmp(cm->transports[i]->trans_name, trans_name) == 0)
	    return 1;
	i++;
    }
    return 0;
}

transport_entry
add_transport_to_cm(CManager cm, transport_entry transport)
{
    int num_trans;
    if (cm->transports == NULL) {
	cm->transports = INT_CMmalloc(sizeof(transport_entry) * 2);
	num_trans = 0;
    } else {
	num_trans = 0;
	while (cm->transports[num_trans] != NULL)
	    num_trans++;
	cm->transports = INT_CMrealloc(cm->transports,
				       sizeof(transport_entry) *
				       (num_trans + 2));
    }
    cm->transports[num_trans] =
	INT_CMmalloc(sizeof(struct _transport_item));
    *(cm->transports[num_trans]) = *transport;
    cm->transports[num_trans + 1] = NULL;
    transport = cm->transports[num_trans];
    transport->cm = cm;
    return transport;
}

int
load_transport(CManager cm, const char *trans_name, int quiet)
{
    transport_entry *trans_list = global_transports;
    transport_entry transport = NULL;
    int i = 0;
#if !NO_DYNAMIC_LINKING
    transport_entry transport_to_free = NULL;
    char *libname;
    lt_dlhandle handle;
#endif


    if (find_transport_in_cm(cm, trans_name))
	return 1;

    while ((trans_list != NULL) && (*trans_list != NULL)) {
	if (strcmp((*trans_list)->trans_name, trans_name) == 0) {
	    transport_entry trans = add_transport_to_cm(cm, *trans_list);
	    if (trans->transport_init) {
		trans->trans_data =
		    trans->transport_init(cm, &CMstatic_trans_svcs, trans);
	    }
	    return 1;
	}
	trans_list++;
	i++;
    }
#if !NO_DYNAMIC_LINKING
    libname =
	INT_CMmalloc(strlen(trans_name) +
		     strlen("lib" CM_LIBRARY_PREFIX "cm") +
		     strlen(MODULE_EXT)
		     + 1);

#ifndef HAVE_WINDOWS_H
    strcpy(libname, "lib" CM_LIBRARY_PREFIX "cm");
#else
    strcpy(libname, CM_LIBRARY_PREFIX "cm");
#endif
    strcat(libname, trans_name);
    strcat(libname, MODULE_EXT);

    lt_dladdsearchdir(EVPATH_MODULE_BUILD_DIR);
    lt_dladdsearchdir(EVPATH_MODULE_INSTALL_DIR);
    handle = CMdlopen(cm->CMTrace_file, libname, 0);
    if (!handle) {
	if (!quiet)
	    fprintf(stderr, "Failed to load required '%s' dll.\n",
		    trans_name);
	if (!quiet)
	    fprintf(stderr,
		    "Search path includes '.', '%s', '%s' and any default search paths supported by ld.so\n",
		    EVPATH_MODULE_BUILD_DIR, EVPATH_MODULE_INSTALL_DIR);
    }
    if (!handle) {
	return 0;
    }
    INT_CMfree(libname);
    transport = INT_CMmalloc(sizeof(struct _transport_item));
    transport_to_free = transport;
    transport->trans_name = strdup(trans_name);
    transport->cm = cm;
    transport->dlhandle = handle;
    transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
    transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
    transport->transport_init = (CMTransport_func)
	lt_dlsym(handle, "initialize");
    transport->listen = (CMTransport_listen_func)
	lt_dlsym(handle, "non_blocking_listen");
    transport->initiate_conn = (CMTransport_conn_func)
	lt_dlsym(handle, "initiate_conn");
    transport->initiate_conn_nonblocking = (CMTransport_NBconn_func)
	lt_dlsym(handle, "initiate_conn_nonblocking");
    transport->finalize_conn_nonblocking = (CMTransport_NBconn_final_func)
	lt_dlsym(handle, "finalize_conn_nonblocking");
    transport->self_check = (CMTransport_self_check_func) lt_dlsym(handle, "self_check");
    transport->connection_eq =
	(CMTransport_connection_eq_func) lt_dlsym(handle, "connection_eq");
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)
	lt_dlsym(handle, "shutdown_conn");
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)
	lt_dlsym(handle, "read_to_buffer_func");
    transport->read_block_func = (CMTransport_read_block_func)
	lt_dlsym(handle, "read_block_func");
    transport->writev_func = (CMTransport_writev_func)
	lt_dlsym(handle, "writev_func");
    transport->writev_complete_notify_func =
	(CMTransport_writev_complete_notify_func)
	lt_dlsym(handle, "writev_complete_notify_func");
    transport->install_pull_schedule_func =
	(CMTransport_install_pull_schedule)
	lt_dlsym(handle, "install_pull_schedule");
    transport->NBwritev_func = (CMTransport_writev_func)
	lt_dlsym(handle, "NBwritev_func");
    transport->set_write_notify = (CMTransport_set_write_notify_func)
	lt_dlsym(handle, "set_write_notify");
    transport->get_transport_characteristics =
	(CMTransport_get_transport_characteristics)
	lt_dlsym(handle, "get_transport_characteristics");
    if (transport->transport_init) {
	transport->trans_data =
	    transport->transport_init(cm, &CMstatic_trans_svcs, transport);
    }
    transport = add_transport_to_cm(cm, transport);
    free(transport_to_free);	/* returned transport is a copy */
#else
    if (strcmp(trans_name, "sockets") == 0) {
	extern transport_entry cmsockets_add_static_transport(CManager cm,
							      CMtrans_services
							      svc);
	transport =
	    cmsockets_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
    if (strcmp(trans_name, "udp") == 0) {
	extern transport_entry cmudp_add_static_transport(CManager cm,
							  CMtrans_services
							  svc);
	transport = cmudp_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#ifdef NNTI_FOUND
    if (strcmp(trans_name, "nnti") == 0) {
	extern transport_entry cmnnti_add_static_transport(CManager cm,
							   CMtrans_services
							   svc);
	transport = cmnnti_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
#ifdef LIBFABRIC_FOUND
    if (strcmp(trans_name, "fabric") == 0) {
	extern transport_entry cmfabric_add_static_transport(CManager cm,
							     CMtrans_services
							     svc);
	transport =
	    cmfabric_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
#ifdef IB_FOUND
    if (strcmp(trans_name, "ib") == 0) {
	extern transport_entry cmib_add_static_transport(CManager cm,
							 CMtrans_services
							 svc);
	transport = cmib_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
#ifdef ENET_FOUND
    if (strcmp(trans_name, "enet") == 0) {
	extern transport_entry cmenet_add_static_transport(CManager cm,
							   CMtrans_services
							   svc);
	transport = cmenet_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
#ifdef ZPL_ENET_AVAILABLE
    if (strcmp(trans_name, "zplenet") == 0) {
	extern transport_entry cmzplenet_add_static_transport(CManager cm,
							   CMtrans_services
							   svc);
	transport = cmzplenet_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
#ifdef UDT4_FOUND
    if (strcmp(trans_name, "udt4") == 0) {
	extern transport_entry cmudt4_add_static_transport(CManager cm,
							   CMtrans_services
							   svc);
	transport = cmudt4_add_static_transport(cm, &CMstatic_trans_svcs);
	transport->data_available = CMDataAvailable;	/* callback
							 * pointer */
	transport->write_possible = CMWriteQueuedData;	/* callback
							 * pointer */
	(void) add_transport_to_cm(cm, transport);
    }
#endif
    if (!transport)
	return 0;
#endif
    CMtrace_out(cm, CMTransportVerbose, "Loaded transport %s.\n",
		trans_name);
    if (global_transports != NULL) {
	global_transports = INT_CMrealloc(global_transports,
					  sizeof(global_transports) * (i +
								       2));
    } else {
	global_transports =
	    INT_CMmalloc(sizeof(global_transports) * (i + 2));
/*	atexit(free_global_transports);*/
    }
    global_transports[i] = transport;
    global_transports[i + 1] = NULL;
    if (transport->install_pull_schedule_func && cm->avail) {
	transport->install_pull_schedule_func(&CMstatic_trans_svcs,
					      transport, &cm->base_time,
					      &cm->period, cm->avail);
	CMtrace_out(cm, CMTransportVerbose,
		    "CM installed pull schedule to transport %s\n",
		    transport->trans_name);
    }

    return 1;
}
