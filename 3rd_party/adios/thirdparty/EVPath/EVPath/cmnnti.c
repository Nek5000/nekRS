/***** Includes *****/
#include "config.h"
#include <sys/types.h>
#ifdef ENET_FOUND
#include <enet/enet.h>
#endif

#include <stdio.h>
#include <errno.h>
#include <pthread.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#undef NDEBUG
#include <assert.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#include <arpa/inet.h>

#include <atl.h>
#include "evpath.h"
#include "cm_transport.h"
#include <Trios_nnti.h>

#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 177)
#  pragma warning (disable: 2259)
#  pragma warning (disable: 981)
#endif

typedef enum {MAP_EACH_MESSAGE_PULL, MAP_EACH_MESSAGE_PUSH,
	      PERSISTENT_MAPPED_BUFFER_PULL, PERSISTENT_MAPPED_BUFFER_PUSH,
	      LAST_PROTOCOL} messaging_protocols;

extern attr_list
libcmnnti_LTX_non_blocking_listen(CManager cm, CMtrans_services svc, 
				  transport_entry trans, attr_list listen_info);
struct nnti_connection_data;

static atom_t CM_PEER_IP = -1;
static atom_t CM_PEER_LISTEN_PORT = -1;
static atom_t CM_NNTI_PORT = -1;
static atom_t CM_NNTI_PARAMS = -1;
static atom_t CM_NNTI_ADDR = -1;
static atom_t CM_NNTI_ENET_CONTROL = -1;
static atom_t CM_NNTI_IMMEDIATE_PULL_WAIT = -1;
static atom_t CM_IP_HOSTNAME = -1;
static atom_t CM_TRANSPORT = -1;
static atom_t CM_NNTI_TRANSPORT = -1;
static atom_t CM_ENET_PORT = -1;
static atom_t CM_ENET_ADDR = -1;
static atom_t CM_TRANSPORT_RELIABLE = -1;

char *NNTI_result_string[] = {
    "NNTI_OK",
    "NNTI_EIO",
    "NNTI_EMSGSIZE",
    "NNTI_ECANCELED",
    "NNTI_ETIMEDOUT",
    "NNTI_EINVAL",
    "NNTI_ENOMEM",
    "NNTI_ENOENT",
    "NNTI_ENOTSUP",
    "NNTI_EEXIST",
    "NNTI_EBADRPC",
    "NNTI_ENOTINIT",
    "NNTI_EPERM",
    "NNTI_EAGAIN"};
#define NNTI_ERROR_STRING(error) ((error == 0)? NNTI_result_string[error] : NNTI_result_string[error-1000])

/* requet queue for rdma pull scheduling */
struct pull_request;
struct pull_request_queue;
struct pull_sched_context;
struct client_message;

typedef int (* nnti_pull_sched_func) (struct pull_request *, struct pull_sched_context *, void *callback_data);

typedef void (* nnti_pull_completion_func) (struct pull_request *, struct pull_sched_context *, void *callback_data);

typedef struct nnti_transport_data {
    CManager cm;
    CMtrans_services svc;
    int socket_fd;
    int self_ip;
    int self_port;
    attr_list listen_attrs;
    char *self_hostname;
    char* incoming;
    NNTI_buffer_t  mr_recvs;
    NNTI_work_request_t  wr_recvs;
    NNTI_transport_t trans_hdl;
    pthread_t listen_thread;
    pthread_t listen_thread2;
    pthread_cond_t cond;
    int use_enet;
    int immediate_pull_wait;
    int shutdown_listen_thread;
    struct nnti_connection_data *connections;
    int cache_maps;

    /* for scheduling */
    struct pull_request_queue *pull_req_queue;
    struct pull_request_queue *ongoing_req_queue;
    nnti_pull_sched_func nnti_pull_scheduler;
    void *nnti_pull_sched_data;
    nnti_pull_completion_func nnti_pull_completion;
    void *nnti_pull_completion_data;
    int request_size;
    struct nnti_connection_data **request_list;

#ifdef ENET_FOUND
    /* enet support */
    ENetHost *enet_server;
#endif
    int enet_listen_port;
    attr_list characteristics;

} *nnti_transport_data_ptr;

typedef struct nnti_connection_data {
    char *peer_hostname;
    int nnti_port;
    char* nnti_params;
    CMbuffer read_buffer;
    int read_buf_len;
    nnti_transport_data_ptr ntd;
    CMConnection conn;
    attr_list attrs;
    struct nnti_connection_data *next;

    NNTI_peer_t peer_hdl;
    long size;
    uint64_t cksum;
    uint64_t raddr;
    int acks_offset;
    NNTI_buffer_t mr_send; // registered memory region to send a message to client
    char *send_buffer;
    NNTI_buffer_t res_addr;
    NNTI_buffer_t buf_addr;
    int piggyback_size_max;
    void *outgoing_mapped_region;
    long outgoing_mapped_size;
    NNTI_buffer_t outgoing_mapped_mr;
    CMbuffer mapped_write_buffer;
    NNTI_buffer_t incoming_mapped_region;
    long incoming_mapped_size;
    NNTI_buffer_t incoming_mapped_mr;
    void *pending_request_msg_info;
    long pending_request_size;

    /* enet support */
    char *remote_host;
    int remote_IP;
    int remote_contact_port;
    int use_enet;
#ifdef ENET_FOUND
    ENetPeer *peer;
    ENetPacket *packet;
#endif
    NNTI_buffer_t mr_pull;
    NNTI_work_request_t wr_pull;
} *nnti_conn_data_ptr;


#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif

#define ENET_PIGGYBACK_SIZE 102500

typedef enum {enet, nnti} control_transport;

typedef struct _send_handle {
    control_transport t;
    int size;
#ifdef ENET_FOUND
    ENetPacket *packet;
#endif
    nnti_conn_data_ptr ncd;
} send_handle;

static int buffer_pack(NNTI_transport_t *trans_hdl, void *input, void **output, uint64_t *output_size)
{
    NNTI_dt_sizeof(trans_hdl, input, output_size);
    *output=(char*)malloc(*output_size);
    NNTI_dt_pack(trans_hdl, input, *output, *output_size);

    return(0);
}

//static int buffer_free(NNTI_transport_t *trans_hdl, void *input)
//{
//    NNTI_dt_free(trans_hdl, input);
//
//    return(0);
//}

static int buffer_unpack(NNTI_transport_t *trans_hdl, void *input, int32_t input_size, void *output)
{
    NNTI_dt_unpack(trans_hdl, output, input, input_size);

    return(0);
}

send_handle
get_control_message_buffer(nnti_conn_data_ptr ncd, struct client_message **mp,
			   int size)
{
    send_handle ret;
    memset(&ret, 0, sizeof(ret));
    if (ncd->use_enet) {
	ret.t = enet;
#ifdef ENET_FOUND
	/* Create a reliable packet of the right size */
	ret.packet = enet_packet_create (NULL, size,
					 ENET_PACKET_FLAG_RELIABLE);
	*mp = (struct client_message *) ret.packet->data;
	memset(ret.packet->data, 0, size);
#endif
    } else {
	assert(size < NNTI_REQUEST_BUFFER_SIZE);
	ret.t = nnti;
	*mp = (struct client_message*)ncd -> send_buffer;
    }	
    ret.ncd = ncd;
    ret.size = size;
    return ret;
}

int
send_control_message(send_handle h)
{
    CManager cm = h.ncd->ntd->cm;
    CMtrans_services svc = h.ncd->ntd->svc;
    if (h.t == enet) {
#ifdef ENET_FOUND
        svc->trace_out(cm, "CMNNTI/ENET control write of %d bytes on peer %p",
		       h.size, h.ncd->peer);
	/* Send the packet to the peer over channel id 0. */
	if (enet_peer_send (h.ncd->peer, 0, h.packet) == -1) {
	    svc->trace_out(cm, "CMNNTI/ENET control write failed.");
	    return 0;
	}
	enet_host_flush (h.ncd->ntd->enet_server);
#ifdef NOT_DEF
	static time_t last_flush_call = 0;
	if (last_flush_call == 0) {
	    enet_host_flush (h.ncd->ntd->enet_server);
	    last_flush_call = time(NULL);
	} else {
	    time_t now = time(NULL);
	    if (now > last_flush_call) {
		last_flush_call = now;
		enet_host_flush (h.ncd->ntd->enet_server);
	    }
	}
#endif
#endif
    } else {
        svc->trace_out(cm, "CMNNTI control write of %d bytes",
		       h.size);
	NNTI_status_t               status;
	int timeout = 1000;
	NNTI_work_request_t send_wr;
	int err = NNTI_send(&h.ncd->peer_hdl, &h.ncd->mr_send, NULL, &send_wr);
	
	if (err != NNTI_OK) {
	    fprintf (stderr, "Error: NNTI_send() returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
	    return 1;
	}
	svc->trace_out(cm, "    NNTI_send() returned. Call wait... ");
	
	/* Wait for message to be sent */
	timeout = 1000;
      again:
	err = NNTI_wait(&send_wr, timeout, &status);
	if (err == NNTI_ETIMEDOUT) {
	  timeout *=2;
	  if (h.ncd->ntd->shutdown_listen_thread) return 0;
	  goto again;
	}
	if (err != NNTI_OK) {
	  fprintf (stderr, "Error: NNTI_wait() for sending returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
	  return 1;
	}
	svc->trace_out(cm, "    NNTI_wait() of send request returned... ");
	
    }
    return 1;
}

#ifdef ENET_FOUND
extern void
enet_non_blocking_listen(CManager cm, CMtrans_services svc,
			 transport_entry trans, attr_list listen_info);
#endif

static nnti_conn_data_ptr
create_nnti_conn_data(svc)
CMtrans_services svc;
{
    nnti_conn_data_ptr nnti_conn_data =
	svc->malloc_func(sizeof(struct nnti_connection_data));
    memset(nnti_conn_data, 0, sizeof(struct nnti_connection_data));
    nnti_conn_data->read_buffer = NULL;
    nnti_conn_data->nnti_port = -1;
    nnti_conn_data->peer_hostname = NULL;
    nnti_conn_data->next = NULL;
    return nnti_conn_data;
}

static void
add_connection(nnti_transport_data_ptr ntd, nnti_conn_data_ptr ncd)
{
    nnti_conn_data_ptr tmp = ntd->connections;
    ntd->connections = ncd;
    ncd->next = tmp;
}

static void
unlink_connection(nnti_transport_data_ptr ntd, nnti_conn_data_ptr ncd)
{
    if (ntd->connections == ncd) {
	ntd->connections = ncd->next;
	ncd->next = NULL;
    } else {
	nnti_conn_data_ptr tmp = ntd->connections;
	while (tmp != NULL) {
	    if (tmp->next == ncd) {
		tmp->next = ncd->next;
		ncd->next = NULL;
		return;
	    }
	}
	printf("Serious internal error, NNTI unlink_connection, connection not found\n");
    }
}


enum {CMNNTI_CONNECT=1, CMNNTI_PIGGYBACK=2, CMNNTI_PULL_REQUEST=3, CMNNTI_PULL_COMPLETE=4, CMNNTI_LAST_MSG_TYPE=7};
char *msg_type_name[] = 
{"NO MESSAGE", "CMNNTI_CONNECT", "CMNNTI_PIGGYBACK", "CMNNTI_PULL_REQUEST", "CMNNTI_PULL_COMPLETE"};

struct connect_message {
    short message_type;
    int nnti_port;
    int enet_port;
    uint32_t enet_ip;
    uint32_t name_len;
    char name[1];
};

static int
initiate_nnti_conn(cm, svc, trans, attrs, nnti_conn_data, conn_attr_list)
CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
nnti_conn_data_ptr nnti_conn_data;
attr_list conn_attr_list;
{
    int int_port_num;
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    char *host_name, *nnti_transport, *params;
    char server_url[256];
    struct connect_message *cmsg;
    NNTI_work_request_t send_wr;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & host_name)) {
	svc->trace_out(cm, "NNTI transport found no NNTI_HOST attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "NNTI transport connect to host %s", host_name);
    }
    if (host_name == NULL)
	return -1;

    if (!query_attr(attrs, CM_NNTI_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (long) &int_port_num)) {
	svc->trace_out(cm, "CMNNTI transport found no NNTI_PORT attribute");
	return -1;
    } else {
	svc->trace_out(cm, "CMNNTI transport connect to port %d", int_port_num);
    }

    if (!query_attr(attrs, CM_NNTI_PARAMS, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (long) &params)) {
	svc->trace_out(cm, "CMNNTI transport found no NNTI_PARAMS attribute");
	params = strdup("");
    } else {
	svc->trace_out(cm, "CMNNTI transport connect with params %s", params);
    }

    if (!get_string_attr(attrs, CM_NNTI_TRANSPORT, &nnti_transport)) {
	svc->trace_out(cm, "NNTI transport found no NNTI_TRANSPORT attribute");

	return -1;
    } else {
        svc->trace_out(cm, "NNTI transport connect using transport %s", nnti_transport);
    }

    sprintf(server_url, "%s://%s:%d/%s", nnti_transport, host_name, int_port_num, params);

    if (ntd->self_port == -1) {
        libcmnnti_LTX_non_blocking_listen(cm, svc, trans, NULL);
    }

#ifdef ENET_FOUND
    if (ntd->enet_server == NULL) {
	enet_non_blocking_listen(cm, svc, trans, NULL);
    }
#endif

    int timeout = 500;
    ntd->svc->trace_out(trans->cm, "Connecting to URL \"%s\"", server_url);
    DROP_CM_LOCK(svc, trans->cm);
    int err = NNTI_connect(&ntd->trans_hdl, server_url, timeout, 
			   &nnti_conn_data->peer_hdl);
    ACQUIRE_CM_LOCK(svc, trans->cm);
    if (err != NNTI_OK) {
        fprintf (stderr, "Error: NNTI_connect() returned non-zero: %d, %s\n", err, NNTI_ERROR_STRING(err));
        return 1;
    }

    /* register memory regions */
    err = NNTI_alloc(&ntd->trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1,
                 NNTI_SEND_SRC, &nnti_conn_data->mr_send);

    
    if (err != NNTI_OK) {
        fprintf (stderr, "Error: NNTI_alloc(SEND_SRC) for client message returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
        return 1;
    }

    cmsg = (void*)NNTI_BUFFER_C_POINTER(&nnti_conn_data->mr_send);
    cmsg->message_type = CMNNTI_CONNECT;
    cmsg->nnti_port = ntd->self_port;
    cmsg->enet_port = ntd->enet_listen_port;
    cmsg->enet_ip = ntd->self_ip;
    cmsg->name_len = strlen(ntd->self_hostname);
    strcpy(&cmsg->name[0], ntd->self_hostname);

    svc->trace_out(cm, "CMNNTI initiate sending connect message to remote host");
    err = NNTI_send(&nnti_conn_data->peer_hdl, &nnti_conn_data->mr_send, 
		    NULL, &send_wr);
    if (err != NNTI_OK) {
        fprintf (stderr, "Error: NNTI_send() returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
        return 1;
    }

    svc->trace_out(trans->cm, " NNTI_send() returned. Call wait... ");
    nnti_conn_data->nnti_port = int_port_num;
    nnti_conn_data->peer_hostname = strdup(host_name);
    nnti_conn_data->ntd = ntd;

    /* Wait for message to be sent */
    NNTI_status_t               status;
    timeout = 500;
 again:
    err = NNTI_wait(&send_wr, timeout, &status);
    if ((err == NNTI_ETIMEDOUT) || (err == NNTI_EAGAIN)) {
	if (nnti_conn_data->ntd->shutdown_listen_thread) return 0;
	timeout *=2;
	goto again;
    }
    if (err != NNTI_OK) {
        fprintf (stderr, "Error: NNTI_wait() for sending returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
        return 1;
    }
    svc->trace_out(trans->cm, " NNTI_wait() of send CONNECT request returned... ");

    svc->trace_out(cm, "--> NNTI Connection established");

    nnti_conn_data->nnti_port = int_port_num;
    nnti_conn_data->nnti_params = params;
    nnti_conn_data->ntd = ntd;
    nnti_conn_data->send_buffer = (char*)NNTI_BUFFER_C_POINTER(&nnti_conn_data->mr_send);;
    return 1;
}

#ifdef ENET_FOUND
static void nnti_enet_service_network(CManager cm, void *void_trans);

extern void
enet_non_blocking_listen(CManager cm, CMtrans_services svc,
			 transport_entry trans, attr_list listen_info)
{
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    ENetAddress address;
    ENetHost * server;

    svc->trace_out(cm, "CMnnti begin enet listen\n");

    address.host = ENET_HOST_ANY;

    if (ntd->enet_server != NULL) {
	/* we're already listening */
	return;
    }
    long seedval = time(NULL) + getpid();
    /* port num is free.  Constrain to range 26000 : 26100 */
    int low_bound = 26000;
    int high_bound = 26100;
    int size = high_bound - low_bound;
    int tries = 10;
    srand48(seedval);
    while (tries > 0) {
	int target = low_bound + size * drand48();
	address.port = target;
	svc->trace_out(cm, "CMnnti trying to bind enet port %d", target);
	
	server = enet_host_create (& address /* the address to bind the server host to */, 
				   0     /* allow dynamic clients (This is supported by the GaTech mod ENET only) */,
				   1      /* allow up to 2 channels to be used, 0 and 1 */,
				   0      /* assume any amount of incoming bandwidth */,
				   0      /* assume any amount of outgoing bandwidth */);
	tries--;
	if (server != NULL) tries = 0;
	if (tries == 5) {
	    /* try reseeding in case we're in sync with another process */
	    srand48(time(NULL) + getpid());
	}
    }
    if (server == NULL) {
	fprintf(stderr, "Failed after 5 attempts to bind to a random port.  Lots of undead servers on this host?\n");
	return;
    }
    ntd->enet_server = server;
    ntd->enet_listen_port = address.port;

    svc->fd_add_select(cm, enet_host_get_sock_fd (server), 
		       (select_list_func) nnti_enet_service_network, (void*)cm, (void*)trans);
    return;
}

static int
initiate_enet_link(CManager cm, CMtrans_services svc, transport_entry trans,
		   nnti_conn_data_ptr nnti_conn_data, char *host_name, int int_port_num)
{
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    ENetAddress address;
    ENetEvent event;
    ENetPeer *peer;
    struct in_addr sin_addr;
    enet_address_set_host (& address, host_name);
    sin_addr.s_addr = address.host;

    svc->trace_out(cm, "Attempting ENET RUDP connection, host=\"%s\", IP = %s, port %d",
		   host_name == 0 ? "(unknown)" : host_name, 
		   inet_ntoa(sin_addr),
		   int_port_num);

    enet_address_set_host (& address, host_name);
    address.port = (unsigned short) int_port_num;

    if (ntd->enet_server == NULL) {
	enet_non_blocking_listen(cm, svc, trans, NULL);
    }

    /* Initiate the connection, allocating the two channels 0 and 1. */
    peer = enet_host_connect (ntd->enet_server, & address, 1, 0);    
    
    if (peer == NULL) {
       fprintf (stderr, 
                "No available peers for initiating an ENet connection.\n");
       exit (EXIT_FAILURE);
    }
    
    /* Wait up to 5 seconds for the connection attempt to succeed. */
    if (enet_host_service (ntd->enet_server, & event, 5000) > 0 &&
        event.type == ENET_EVENT_TYPE_CONNECT) {
	svc->trace_out(cm, "Connection to %s:%d succeeded.\n", inet_ntoa(sin_addr), address.port);

    } else {
        /* Either the 5 seconds are up or a disconnect event was */
        /* received. Reset the peer in the event the 5 seconds   */
        /* had run out without any significant event.            */
        enet_peer_reset (peer);

        printf ("Connection to %s:%d failed.", inet_ntoa(sin_addr), address.port);
	return 0;
    }

    svc->trace_out(cm, "--> Enet Connection established");
    nnti_conn_data->remote_host = host_name == NULL ? NULL : strdup(host_name);
    nnti_conn_data->remote_IP = address.host;
    nnti_conn_data->remote_contact_port = int_port_num;
    nnti_conn_data->ntd = ntd;
    nnti_conn_data->peer = peer;
    peer->data = nnti_conn_data;
    return 1;
}
#endif

#ifdef ENET_FOUND
static int
initiate_enet_conn(CManager cm, CMtrans_services svc, transport_entry trans,
	      attr_list attrs, nnti_conn_data_ptr nnti_conn_data,
	      attr_list conn_attr_list)
{
    int int_port_num;
    char *host_name;
    static int host_ip = 0;
    (void)conn_attr_list;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & host_name)) {
	svc->trace_out(cm, "Cmnnti/Enet transport found no CM_IP_HOSTNAME attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "Cmnnti/Enet transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_ENET_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & host_ip)) {
	svc->trace_out(cm, "Cmnnti/Enet transport found no CM_ENET_ADDR attribute");
	/* wasn't there */
	host_ip = 0;
    } else {
        svc->trace_out(cm, "Cmnnti/Enet transport connect to host_IP %lx", host_ip);
    }
    if ((host_name == NULL) && (host_ip == 0))
	return -1;

    if (!query_attr(attrs, CM_ENET_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & int_port_num)) {
	svc->trace_out(cm, "Cmnnti/Enet transport found no CM_ENET_PORT attribute");
	return -1;
    } else {
        svc->trace_out(cm, "Cmnnti/Enet transport connect to port %d", int_port_num);
    }

    return initiate_enet_link(cm, svc, trans, nnti_conn_data, host_name, int_port_num);
}
#endif

/* 
 * Initiate a connection to a nnti group.
 */
extern CMConnection
libcmnnti_LTX_initiate_conn(cm, svc, trans, attrs)
CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
{
    nnti_conn_data_ptr nnti_conn_data = create_nnti_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;

    if (initiate_nnti_conn(cm, svc, trans, attrs, nnti_conn_data, conn_attr_list) != 1) {
	return NULL;
    }

    /* this size might be overridden */
    nnti_conn_data->piggyback_size_max = NNTI_REQUEST_BUFFER_SIZE;

#ifdef ENET_FOUND
    int enet_conn_status;
    sleep(1);

    enet_conn_status = initiate_enet_conn(cm, svc, trans, attrs, nnti_conn_data, conn_attr_list);
    switch (enet_conn_status) {
    case -1:
      svc->trace_out(cm, "NO ENET Contact info provided.");
      nnti_conn_data->use_enet = 0;
      break;
    case 0:
      svc->trace_out(cm, "ENET Connection failed.");
      return NULL;
    default:
      nnti_conn_data->use_enet = 1;
      nnti_conn_data->piggyback_size_max = ENET_PIGGYBACK_SIZE;
    }
#endif

    add_attr(conn_attr_list, CM_IP_HOSTNAME, Attr_String,
	     (attr_value) strdup(nnti_conn_data->peer_hostname));
    add_attr(conn_attr_list, CM_NNTI_PORT, Attr_Int4,
	     (attr_value) (long)nnti_conn_data->nnti_port);
    add_attr(conn_attr_list, CM_NNTI_PARAMS, Attr_String,
	     (attr_value) (long)nnti_conn_data->nnti_params);

    conn = svc->connection_create(trans, nnti_conn_data, conn_attr_list);
    add_connection(nnti_conn_data->ntd, nnti_conn_data);
    nnti_conn_data->conn = conn;
    nnti_conn_data->attrs = conn_attr_list;
    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 */
extern int
libcmnnti_LTX_self_check(cm, svc, trans, attrs)
CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
{
    nnti_transport_data_ptr ntd = trans->trans_data;
    int int_port_num;
    char *host_name;
    char my_host_name[256];
    static int IP = 0;

    get_IP_config(my_host_name, sizeof(host_name), &IP, NULL, NULL, NULL,
		  NULL, svc->trace_out, (void *)cm);

    if (IP == 0) {
	IP = ntohl(INADDR_LOOPBACK);
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & host_name)) {
	svc->trace_out(cm, "CMself check NNTI transport found no IP_HOST attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_NNTI_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & int_port_num)) {
	svc->trace_out(cm, "CMself check NNTI transport found no NNTI_PORT attribute");
	return 0;
    }

    if (host_name && (strcmp(host_name, my_host_name) != 0)) {
	svc->trace_out(cm, "CMself check - Hostnames don't match");
	return 0;
    }
    if (int_port_num != ntd->self_port) {
	svc->trace_out(cm, "CMself check - Ports don't match");
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

extern int
libcmnnti_LTX_connection_eq(cm, svc, trans, attrs, ncd)
CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
nnti_conn_data_ptr ncd;
{

    int int_port_num;
    char *host_name = NULL;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(long) & host_name)) {
	svc->trace_out(cm, "NNTI transport found no NNTI_HOST attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "NNTI transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_NNTI_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (long) &int_port_num)) {
	svc->trace_out(cm, "Conn Eq CMNnti transport found no NNTI_PORT attribute");
	return 0;
    }
    svc->trace_out(cm, "CMNnti Conn_eq comparing host/ports %s/%d and %s/%d",
		   ncd->peer_hostname, ncd->nnti_port,
		   host_name, int_port_num);

    if ((ncd->nnti_port == int_port_num) &&
	(strcmp(ncd->peer_hostname, host_name) == 0)) {
	svc->trace_out(cm, "CMNnti Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "CMNnti Conn_eq returning FALSE");
    return 0;
}


struct piggyback {
  unsigned int size;          // size of message payload
  char payload[1];
};

struct pull_request {
    unsigned long size;          /* size of message to pull */
    void *msg_info;	       /* returned unchanged */
    int reuse_prior_mapping;
    char packed_src_buf[1];
};

struct pull_complete_notify {
  void *msg_info;
};

struct client_message {
  short message_type;
  union {
    struct piggyback pig;
    struct pull_request pull;
    struct pull_complete_notify pull_complete;
  };
};

enum pull_request_state {
  GET_SCHED_QUEUED,             // in queue, ready to be scheduled
  GET_SCHED_GET_ISSUED,         // GET command has been issued
  GET_SCHED_COMPLETED,          // GET command has finished successfully
  GET_SCHED_ERROR,              // some error/failure has happened
};

struct pull_sched_context {
  nnti_conn_data_ptr ncd;       // from which connection
  int target_stone_id;          // local target stone
  void *context_data;           // other context informaiton such as application phases
  // TODO: add source process id and data format id
};

/* queue structure for delayed pull requests */
struct pull_request_queue {
  enum pull_request_state state;
  int num_errors;
  struct pull_sched_context context;
  struct pull_request_queue *next;
  struct pull_request request;
};

static void
handle_pull_request_message(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
			    struct client_message *m);
static void
handle_pull_complete_message(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
			     struct client_message *m);

static void
handle_control_request(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans, 
		       struct client_message *m);

static int nnti_conn_count = 0;
#define WAIT_MODE_BOTH 1
#define WAIT_MODE_REQUEST 2
#define WAIT_MODE_RDMA 3
typedef struct listen_struct {
    nnti_transport_data_ptr ntd;
    CMtrans_services svc;
    transport_entry trans;
    int wait_mode;
} *listen_struct_p;

static void
handle_request_buffer_event(listen_struct_p lsp, NNTI_status_t *wait_status)
{
    nnti_transport_data_ptr ntd = lsp->ntd;
    CMtrans_services svc = lsp->svc;
    transport_entry trans = lsp->trans;
    struct connect_message *cm = (struct connect_message *)(wait_status->start+wait_status->offset);
	    
    nnti_conn_data_ptr ncd = ntd->connections;
    while (ncd != NULL) {
	if (strcmp(&wait_status->src.url[0], &ncd->peer_hdl.url[0]) == 0) {
	    ntd->svc->trace_out(trans->cm, "NNTI data available on existing connection, from host %s, port %d, type %s (%d)", 
				ncd->peer_hostname, ncd->nnti_port,
				msg_type_name[cm->message_type], cm->message_type);
	    break;
	}
	ncd = ncd->next;
    }
    if ((cm->message_type <= 0) || (cm->message_type > CMNNTI_LAST_MSG_TYPE)) {
        printf("BAD INCOMING SHORT MESSAGE! (value was %d)\n", cm->message_type);
	exit(1);
    }
    switch (cm->message_type){
    case CMNNTI_CONNECT: 
    {
	int err;
	attr_list conn_attr_list = NULL;
		
	ntd->svc->trace_out(trans->cm, "  client %s:%d  (enet %d, %x) is connecting",
			    cm->name, cm->nnti_port, cm->enet_port, cm->enet_ip);
	      
	assert(ncd == NULL);
	ncd = create_nnti_conn_data(svc);
	ncd->ntd = ntd;
	ncd->peer_hdl = wait_status->src;
	ncd->peer_hostname = strdup(cm->name);
	ncd->nnti_port = cm->nnti_port;
	ncd->remote_contact_port = cm->enet_port;
	ncd->remote_IP = cm->enet_ip;
	ncd->piggyback_size_max = NNTI_REQUEST_BUFFER_SIZE;

	ncd->size = NNTI_REQUEST_BUFFER_SIZE;
	ncd->raddr = 0;
	ncd->cksum = 0;
	ncd->acks_offset=(nnti_conn_count++)*NNTI_REQUEST_BUFFER_SIZE;
	//            ncd->buf_addr = cm->buf_addr;
	//            ncd->res_addr = cm->res_addr;
		
	/* register memory region for sending acknowledgement to this client */
	ntd->svc->trace_out(trans->cm, "   ID %d register small send buffer to client %d (offset=%d)...",
			    cm->nnti_port, nnti_conn_count-1, ncd->acks_offset);
	
	conn_attr_list = create_attr_list();
	ncd->conn = svc->connection_create(trans, ncd, conn_attr_list);
	ncd->attrs = conn_attr_list;
	add_connection(ntd, ncd);
	/* alloc memory regions */
	err = NNTI_alloc(&ntd->trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1,
                 NNTI_SEND_SRC, &ncd->mr_send);

    
	if (err != NNTI_OK) {
	    fprintf (stderr, "Error: NNTI_alloc(SEND_SRC) for client message returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
	    return;
	}

	ncd->send_buffer = (void*)NNTI_BUFFER_C_POINTER(&ncd->mr_send);
    }
    break;
    case CMNNTI_PIGGYBACK:
    {
	struct client_message *m = (struct client_message *)cm;
	CMbuffer read_buffer;
	if (ncd == NULL) {
	    printf("Incoming message failed to match connection!\n");
	}
	read_buffer = ntd->svc->get_data_buffer(trans->cm, (int)m->pig.size);
	
	memcpy(&((char*)read_buffer->buffer)[0], &(m->pig.payload[0]), m->pig.size);
	
	ntd->svc->add_buffer_to_pending_queue(trans->cm, ncd->conn, read_buffer, m->pig.size);
    }
    break;
    default:
    {
	struct client_message *m = (struct client_message *)(wait_status->start+wait_status->offset);
	if (ncd == NULL) {
	    printf("Default handler: incoming message failed to match connection!\n");
	}
	handle_control_request(ncd, svc, trans, m);
    }
    }
}

int
listen_thread_func(void *vlsp)
{
    int timeout = 1000;
    listen_struct_p lsp = vlsp;
    nnti_transport_data_ptr ntd = lsp->ntd;
    CMtrans_services svc = lsp->svc;
    transport_entry trans = lsp->trans;
    int err;
    NNTI_status_t               wait_status;

    ACQUIRE_CM_LOCK(svc, trans->cm);
    int buf_size = 10;  // initial size
    NNTI_work_request_t **wr_list = malloc(sizeof(wr_list[0]) * buf_size);

    while (1) {
	unsigned int which;
	int wr_count = 1;
	if (lsp->wait_mode != WAIT_MODE_RDMA) {
	    wr_list[0] = &ntd->wr_recvs;
	} else {
	    wr_count = 0;
	}
	if (lsp->wait_mode != WAIT_MODE_REQUEST) {
	    if (buf_size < ntd->request_size) {
		wr_list = realloc(wr_list, sizeof(wr_list[0]) * ntd->request_size);
		buf_size = ntd->request_size;
	    }
	    for(which=1; which <= ntd->request_size; which++) {
		if (ntd->request_list[which-1] != NULL) {
		    wr_list[wr_count] = &ntd->request_list[which-1]->wr_pull;
		    wr_count++;
		}
	    }
	}
	//err = NNTI_wait(&ntd->mr_recvs, NNTI_RECV_QUEUE, timeout, &wait_status);
	if (wr_count > 0) {
	    DROP_CM_LOCK(svc, trans->cm);
	    err = NNTI_waitany(wr_list, wr_count, timeout, &which, &wait_status);
	    ACQUIRE_CM_LOCK(svc, trans->cm);
	} else {
	    svc->cond_wait_CM_lock((trans->cm),&(ntd->cond),  __FILE__, __LINE__);
	    /* nothing to do but go back and wait for other stuff */
	    continue;
	}	    
	if (ntd->shutdown_listen_thread) {
	    DROP_CM_LOCK(svc, trans->cm);
	    return 0;
	}
	if ((err == NNTI_ETIMEDOUT) || (err==NNTI_EAGAIN)) {
	    //ntd->svc->trace_out(trans->cm, "NNTI_wait() on receiving result %d timed out...", err);

            // TODO: check if it's time to schedule Puts and check progress of outstanding requests

	    continue;
        } else if (err != NNTI_OK) {
            fprintf (stderr, "Error: NNTI_wait() on receiving result returned non-zero: %d %s  which is %d, status is %d\n", err, NNTI_ERROR_STRING(err), which, wait_status.result);
	    DROP_CM_LOCK(svc, trans->cm);
            return 1;
        } else {
	    ntd->svc->trace_out(trans->cm, "  message arrived: msg wait_status=%d size=%lu offs=%lu addr=%lu, offset was %ld",
		     wait_status.result, wait_status.length, wait_status.offset, wait_status.start+wait_status.offset, wait_status.offset);
        }

        if (wait_status.result == NNTI_OK) {
	    if ((which == 0) && (lsp->wait_mode != WAIT_MODE_RDMA)) {
		handle_request_buffer_event(lsp, &wait_status);
		ntd->svc->trace_out(trans->cm, "Clearing request buffer work request");
		NNTI_destroy_work_request(&ntd->wr_recvs);
		ntd->svc->trace_out(trans->cm, "init new request buffer work request");
		NNTI_create_work_request(&ntd->mr_recvs, &ntd->wr_recvs);
	    } else {
		/* handle a completed pull request */
		send_handle h;
		struct client_message *r;
		nnti_conn_data_ptr ncd = NULL;
		int i;
		for (i=0; i < ntd->request_size; i++) {
		    if ((ntd->request_list[i] != NULL) &&
			(&ntd->request_list[i]->wr_pull == wr_list[which])) {
			ncd = ntd->request_list[i];
			ntd->request_list[i] = NULL;
		    }
		}
		if (ncd == NULL) {
		    printf("Ugh.  Didn't find connection data for completed pull.  Dying.\n");
		}
		h = get_control_message_buffer(ncd, &r, sizeof(*r));
		r->message_type = CMNNTI_PULL_COMPLETE;
		r->pull_complete.msg_info = ncd->pending_request_msg_info;
		svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET done with pull, returning control message, type %d, local_info %p", r->message_type, r->pull_complete.msg_info);
		if (send_control_message(h) == 0) {
		    svc->trace_out(ncd->ntd->cm, "--- control message send failed!");
		}	  
		
		/* kick this upstairs */
		svc->add_buffer_to_pending_queue(trans->cm, ncd->conn, ncd->read_buffer, ncd->pending_request_size);
//		if (!ncd->ntd->cache_maps) {
//		    ncd->read_buffer = NULL;
//		    ntd->svc->trace_out(trans->cm, "unregistering memory");
//		    NNTI_unregister_memory(&ncd->mr_pull);
//		}
		ntd->svc->trace_out(trans->cm, "Clearing recvs work request\n");
		NNTI_clear_work_request(wr_list[which]);
	    }
	}
        // TODO: check if it's time to schedule Puts and check progress of outstanding requests
    }
}

#ifdef ENET_FOUND
static void *
enet_accept_conn(nnti_transport_data_ptr ntd, transport_entry trans, 
		 ENetAddress *address);

static
void
nnti_enet_service_network(CManager cm, void *void_trans)
{
    transport_entry trans = (transport_entry) void_trans;
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    CMtrans_services svc = ntd->svc;
    ENetEvent event;
    
    if (!ntd->enet_server) {
	printf("nnti_enet_service network returning\n");
	return;
    }

    /* Wait up to 1 milliseconds for an event. */
    while (enet_host_service (ntd->enet_server, & event, 1) > 0) {
        switch (event.type) {
	case ENET_EVENT_TYPE_NONE:
	    break;
        case ENET_EVENT_TYPE_CONNECT: {
	    void *nnti_connection_data;
	    svc->trace_out(cm, "A new client connected from %x:%u.\n", 
			   event.peer -> address.host,
			   event.peer -> address.port);

	    nnti_connection_data = enet_accept_conn(ntd, trans, &event.peer->address);

	    ((nnti_conn_data_ptr)nnti_connection_data)->use_enet = 1;
	    ((nnti_conn_data_ptr)nnti_connection_data)->piggyback_size_max = ENET_PIGGYBACK_SIZE;
            /* Store any relevant client information here. */
            event.peer -> data = nnti_connection_data;
	    ((nnti_conn_data_ptr)nnti_connection_data)->peer = event.peer;

            break;
	}
        case ENET_EVENT_TYPE_RECEIVE: {
	    nnti_conn_data_ptr ncd = event.peer->data;
	    struct client_message *m = (struct client_message *) event.packet->data;
	    svc->trace_out(cm, "An ENET packet of length %u was received on channel %u, message type %s(%d)",
			   (unsigned int) event.packet -> dataLength,
			   (unsigned int) event.channelID,
			   msg_type_name[m->message_type], m->message_type);
	    if ((m->message_type <= 0) || (m->message_type > CMNNTI_LAST_MSG_TYPE)) {
		printf("BAD INCOMING ENET SHORT MESSAGE!\n");
		exit(1);
	    }
	    if (m->message_type == CMNNTI_PIGGYBACK){
		int piggyback_size;
		ncd->packet = event.packet;
	      
		CMbuffer read_buffer = ntd->svc->get_data_buffer(trans->cm, (int)m->pig.size);
	      
		piggyback_size = m->pig.size;
		memcpy(&((char*)read_buffer->buffer)[0], &(m->pig.payload[0]), m->pig.size);
	      

		enet_packet_destroy(event.packet);
		/* kick this upstairs */
		svc->trace_out(cm, "We received piggybacked data of size %d %x.",
			       piggyback_size, piggyback_size);
		ncd->read_buffer = read_buffer;
		ncd->read_buf_len = piggyback_size;
		trans->data_available(trans, ncd->conn);
		svc->return_data_buffer(trans->cm, read_buffer);
//		svc->add_buffer_to_pending_queue(trans->cm, ncd->conn, read_buffer, piggyback_size);
	    } else {
		handle_control_request(ncd, svc, trans, m);
		enet_packet_destroy (event.packet);
	    }
            break;
	}           
        case ENET_EVENT_TYPE_DISCONNECT: {
	    nnti_conn_data_ptr nnti_conn_data = event.peer -> data;
	    svc->trace_out(nnti_conn_data->ntd->cm, "Got a disconnect on connection %p\n",
		event.peer -> data);

            nnti_conn_data = event.peer -> data;
	    nnti_conn_data->peer = NULL;
	    /*	    nnti_conn_data->read_buffer_len = -1;*/

        }
	}
    }
}
static
void
nnti_enet_service_network_lock(CManager cm, void *void_trans)
{
    transport_entry trans = (transport_entry) void_trans;
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    CMtrans_services svc = ntd->svc;
    ACQUIRE_CM_LOCK(svc, cm);
    nnti_enet_service_network(cm, void_trans);
    DROP_CM_LOCK(svc, cm);
}
#endif

static void
handle_control_request(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
		       struct client_message *m)
{

  switch (m->message_type) {
  case CMNNTI_PULL_REQUEST:{ 
      handle_pull_request_message(ncd, svc, trans, m);
      break;
  }
  case CMNNTI_PULL_COMPLETE:{
      handle_pull_complete_message(ncd, svc, trans, m);
      break;
  }
  default:
      printf("Bad control request! type %d\n", m->message_type);
  }

}

#ifdef ENET_FOUND
/* 
 * Accept enet connection
 */
static void *
enet_accept_conn(nnti_transport_data_ptr ntd, transport_entry trans, 
		 ENetAddress *address)
{
    CMtrans_services svc = ntd->svc;
    nnti_conn_data_ptr ncd = ntd->connections;
    int verbose = -1;
    attr_list conn_attr_list = NULL;;

 restart:
    ncd = ntd->connections;
    if (verbose >=1) printf("NCD is %p\n", ncd);
    
    while (ncd && (ncd->remote_IP != address->host) && 
	   (ncd->remote_contact_port != address->port)) {
	if (verbose >=1) {
	    printf("NCD remote IP %x, address->host %x\n", ncd->remote_IP, address->host);
	    printf("NCD remote contact port %d, address->port %d\n", ncd->remote_contact_port, address->port);
	}
	ncd = ncd->next;
    }
    if (ncd == NULL) {
      printf("Waiting...\n");
      sleep(1);
      verbose++;
      goto restart;
    }

    conn_attr_list = ncd->attrs;

    add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4, (void*)(long)address->host);
    ncd->remote_IP = address->host;
    ncd->remote_contact_port = address->port;

    if (ncd->remote_host != NULL) {
	svc->trace_out(ntd->cm, "Accepted NNTI/ENET RUDP connection from host \"%s\"",
		       ncd->remote_host);
    } else {
	svc->trace_out(ntd->cm, "Accepted NNTI/ENET RUDP connection from UNKNOWN host");
    }
    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (long)ncd->remote_contact_port);
    svc->trace_out(ntd->cm, "Remote host (IP %x) is listening at port %d\n",
		   ncd->remote_IP,
		   ncd->remote_contact_port);
    return ncd;
}

static void
setup_enet_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_list)
{
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) trans->trans_data;
    static int socket_global_init = 0;

    ENetAddress address;
    ENetHost * server;
    long seedval = time(NULL) + getpid();
    /* port num is free.  Constrain to range 26000 : 26100 */
    int low_bound, high_bound;
    int size, tries;
    
    get_IP_config(NULL, 0, NULL, &low_bound, &high_bound,
		  NULL, listen_list, svc->trace_out, (void *)cm);
    size = high_bound - low_bound;

    tries = 50;

    if (socket_global_init++ == 0) {
	if (enet_initialize () != 0) {
	    fprintf (stderr, "An error occurred while initializing ENet.\n");
	    //return EXIT_FAILURE;
	}
    }
    svc->add_periodic_task(cm, 0, 30*1000, nnti_enet_service_network_lock, (void*)trans);
    svc->trace_out(cm, "CMNNTI begin ENET listen\n");

    srand48(seedval);
    address.host = ENET_HOST_ANY;
    while (tries > 0) {
	int target = low_bound + size * drand48();
	address.port = target;
	svc->trace_out(cm, "Cmnnti/Enet trying to bind port %d", target);
	
	server = enet_host_create (& address /* the address to bind the server host to */, 
				   0     /* allow up to 4095 clients and/or outgoing connections */,
				   1      /* allow up to 2 channels to be used, 0 and 1 */,
				   0      /* assume any amount of incoming bandwidth */,
				   0      /* assume any amount of outgoing bandwidth */);
	tries--;
	if (server != NULL) tries = 0;
	if (tries == 5) {
	    /* try reseeding in case we're in sync with another process */
	    srand48(time(NULL) + getpid());
	}
    }
    if (server == NULL) {
	fprintf(stderr, "Failed after %d attempts to bind to a random port.  Lots of undead servers on this host?\n", tries);
	return;
    }
    ntd->enet_server = server;
    ntd->enet_listen_port = address.port;
    svc->trace_out(cm, "CMNNTI  ENET listen at port %d, server %p", address.port, server);
    svc->fd_add_select(cm, enet_host_get_sock_fd (server), 
		       (select_list_func) nnti_enet_service_network, (void*)cm, (void*)trans);
    add_attr(listen_list, CM_ENET_PORT, Attr_Int4,
	     (attr_value) (long)address.port);
}
#else
static void
setup_enet_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_list) {}
#endif

static void
setup_nnti_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_list)
{
    static NNTI_transport_t trans_hdl;
    nnti_transport_data_ptr ntd = trans->trans_data;
    static int initialized = 0;
    char url[256];
    char *last_colon, *first_colon, *last_slash;
    char *hostname;
    int incoming_size = 100;
    int int_port_num = 0;
    int err;

    /* hope to eliminate this at some point */
    setenv("TRIOS_NNTI_USE_RDMA_TARGET_ACK", "FALSE", 1);

    if (!initialized) {
	NNTI_result_t rc;
	DROP_CM_LOCK(svc, trans->cm);
	rc = NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);
	ACQUIRE_CM_LOCK(svc, trans->cm);
	if (rc != NNTI_OK) {
	    printf("Failed to initialize NNTI transport, exiting\n");
	    exit(1);
	}
	initialized++;
    }
    NNTI_get_url(&trans_hdl, url, sizeof(url));
    if (getenv("NNTI_INCOMING_SIZE")) {
	char *size_str = getenv("NNTI_INCOMING_SIZE");
	int incoming_tmp;
	if (sscanf(size_str, "%d", &incoming_tmp) != 1) {
	    ntd->svc->trace_out(trans->cm, "Failed to parse NNTI_INCOMING_SIZE, \"%s\"", size_str);
	} else {
	    if (incoming_tmp < incoming_size) {
		ntd->svc->trace_out(trans->cm, "NNTI_INCOMING_SIZE %d ignored, smaller than default \"%d\"", incoming_tmp, incoming_size);
	    } else {
		ntd->svc->trace_out(trans->cm, "NNTI_INCOMING_SIZE set to %d", incoming_tmp);
		incoming_size = incoming_tmp;
	    }
	}
    }
		
    ntd->svc->trace_out(trans->cm, "NNTI_init succeeded, listening on url %s", url);
    last_colon = rindex(url, ':');
    *last_colon = 0;
    first_colon = index(url, ':');
    *first_colon = 0;
    hostname = (first_colon + 1);
    while(hostname[0] == '/') hostname++;
    last_slash = index(last_colon+1, '/');
    *(last_slash++) = 0;
    sscanf((last_colon + 1), "%d", &int_port_num);

    add_attr(listen_list, CM_IP_HOSTNAME, Attr_String,
	     (attr_value) strdup(hostname));
    add_attr(listen_list, CM_NNTI_PORT, Attr_Int4,
	     (attr_value) (long) int_port_num);
    if (strlen(last_slash)) {
	add_string_attr(listen_list, CM_NNTI_PARAMS, strdup(last_slash));
    }
    add_attr(listen_list, CM_TRANSPORT, Attr_String,
	     (attr_value) strdup("nnti"));
    add_attr(listen_list, CM_NNTI_TRANSPORT, Attr_String,
	     (attr_value) strdup(url));

    /* register memory regions */
    err = NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, incoming_size,
			       NNTI_RECV_QUEUE, &ntd->mr_recvs);

    NNTI_create_work_request(&ntd->mr_recvs, &ntd->wr_recvs);

    if (err != NNTI_OK) {
	fprintf (stderr, "Error: NNTI_alloc(NNTI_RECV_QUEUE) for client messages returned non-zero: %d %s\n", err, NNTI_ERROR_STRING(err));
	return;
    } else {
	ntd->svc->trace_out(trans->cm, "Successfully registered memory on listen side incoming %p", (void*)NNTI_BUFFER_C_POINTER(&ntd->mr_recvs));
    }
    ntd->incoming = (void*)NNTI_BUFFER_C_POINTER(&ntd->mr_recvs);
    
    ntd->self_port = int_port_num;
    ntd->trans_hdl = trans_hdl;
    listen_struct_p lsp = malloc(sizeof(*lsp));
    lsp->svc = svc;
    lsp->trans = trans;
    lsp->ntd = ntd;
    lsp->wait_mode = WAIT_MODE_REQUEST;
    ntd->shutdown_listen_thread = 0;
    ntd->listen_thread = 0;
    ntd->self_hostname = strdup(hostname);
    pthread_cond_init(&(ntd->cond), NULL);   
    err = pthread_create(&ntd->listen_thread, NULL, (void*(*)(void*))listen_thread_func, lsp);
    lsp = malloc(sizeof(*lsp));
    lsp->svc = svc;
    lsp->trans = trans;
    lsp->ntd = ntd;
    lsp->wait_mode = WAIT_MODE_RDMA;
    err = pthread_create(&ntd->listen_thread2, NULL, (void*(*)(void*))listen_thread_func, lsp);
}

extern attr_list
libcmnnti_LTX_non_blocking_listen(cm, svc, trans, listen_info)
CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list listen_info;
{
    attr_list listen_list;
    nnti_transport_data_ptr ntd = trans->trans_data;
    int use_enet = 0;
    char *enet = getenv("NNTI_ENET");
    int use_nnti = 1;
    int immediate_pull = 1;
    if (ntd->listen_attrs != NULL) {
	return ntd->listen_attrs;
    }
    if (enet) {
      sscanf(enet, "%d", &use_enet);
    }
    if (listen_info) {
	get_int_attr(listen_info, CM_NNTI_ENET_CONTROL, &use_enet);
	get_int_attr(listen_info, CM_NNTI_IMMEDIATE_PULL_WAIT, &immediate_pull);
    }
    ntd->immediate_pull_wait = immediate_pull;

    listen_list = create_attr_list();
    if (use_nnti) {
	setup_nnti_listen(cm, svc, trans, listen_list);
    }
    if (use_enet) {
	setup_enet_listen(cm, svc, trans, listen_list);
    }
    ntd->listen_attrs = listen_list;
    ntd->use_enet = use_enet;

    return listen_list;

}

#ifdef NEED_IOVEC_DEFINE
struct iovec {
    void *iov_base;
    long iov_len;
};

#endif

/* 
 *  This function will not be used unless there is no read_to_buffer function
 *  in the transport.  It is an example, meant to be copied in transports 
 *  that are more efficient if they allocate their own buffer space.
 */
extern void *
libcmnnti_LTX_read_block_func(svc, ncd, actual_len, offset_ptr)
CMtrans_services svc;
nnti_conn_data_ptr ncd;
int *actual_len;
int *offset_ptr;
{
    CMbuffer cb;

    if (ncd->read_buf_len == -1) return NULL;

    *actual_len = ncd->read_buf_len;
    *offset_ptr = 0;
    cb = ncd->read_buffer;
    ncd->read_buf_len = 0;
    ncd->read_buffer = NULL;
    return cb;
}

typedef struct {
    int send_id;
    int size;
    NNTI_buffer_t mr;
    CMbuffer write_buffer;
    CMcompletion_notify_func notify_func;
    void *notify_client_data;
} *nnti_message_info;

static int send_request = 32;

static int
copy_full_buffer_and_send_pull_request(CMtrans_services svc, nnti_conn_data_ptr ncd,
				       struct iovec *iov, int iovcnt, attr_list attrs,
				       CMcompletion_notify_func notify_func,
				       void *notify_client_data)
{
    send_handle h;
    struct client_message *m;
    CMbuffer write_buffer = NULL;
    NNTI_buffer_t mr;
    char *data = NULL;
    char *ptr;
    int i, write_size = 0;
    int err = NNTI_OK;
    long register_size = 0;
    nnti_message_info local_message_info = malloc(sizeof(*local_message_info));
    void *packed = NULL;
    uint64_t packed_size;

    static void **mapped_segments = NULL;
    static uint64_t *mapped_lengths = NULL;
    static int mapped_count = -1;
    int can_reuse_mapping = 0;

    for(i=0; i<iovcnt; i++) write_size+= iov[i].iov_len;
//    for(i=0; i<iovcnt; i++) printf("NNTI WRITE VEC[%d], %p for %ld\n",
//				   i, iov[i].iov_base, iov[i].iov_len);

/* choices at this point:
  - if there is a notify_func, we could leave the data where it is 
   and tell the application when we're done with it.
  - if no notify func, we must copy the data here.  We can see if the 
   last mapped memory is acceptable.  If so, we can reuse without remapping.
  - else we map new memory.
*/
    
    if (notify_func) {
	can_reuse_mapping = 1;
	/* OK, we're not going to copy the data */
	if (mapped_count == iovcnt) {
	    int i;
	    for(i=0; i < mapped_count; i++) {
		if ((iov[i].iov_len != mapped_lengths[i]) ||
		    (iov[i].iov_base != mapped_segments[i])) {
		    can_reuse_mapping = 0;
		    svc->trace_out(ncd->ntd->cm, "CMNNTI already mapped data, doesn't match write, buf %d, %p vs. %p, %d vs. %d",
				   i, iov[i].iov_base, mapped_segments[i], iov[i].iov_len, mapped_lengths[i]);
		    break;
		}
	    }
	} else {
	    svc->trace_out(ncd->ntd->cm, "CMNNTI either no already mapped data, or wrong buffer count");
	    can_reuse_mapping = 0;
	}
	if (can_reuse_mapping) {
	    /* who'd have guessed the data is in exactly the same place!  
	       Amazing!  We just have to send a pull request. */
	    svc->trace_out(ncd->ntd->cm, "CMNNTI reusing mapping of already mapped user-level data");
	    mr = ncd->outgoing_mapped_mr;
	} else {
	    void **segments = malloc(iovcnt * sizeof(segments[0]));
	    uint64_t *lengths = malloc(iovcnt * sizeof(lengths[0]));
	    if (ncd->outgoing_mapped_region != NULL) {
		ncd->outgoing_mapped_region = NULL;
		NNTI_unregister_memory(&ncd->outgoing_mapped_mr);
		svc->trace_out(ncd->ntd->cm, "CMNNTI unregistering previously mapped region at %p, size %d",
			       ncd->outgoing_mapped_region, ncd->outgoing_mapped_size);
	    }
	    for(i=0; i<iovcnt; i++) {
		segments[i] = iov[i].iov_base;
		lengths[i] = iov[i].iov_len;
	    }

	    svc->trace_out(ncd->ntd->cm, "CMNNTI registering application buffers in place, size %d", write_size);
	    err = NNTI_register_segments(&ncd->ntd->trans_hdl, (char**)segments,
					 lengths, iovcnt, NNTI_GET_SRC, &mr);

//	    for(i=0; i<iovcnt; i++) {
//		printf("registered segment[%d] %ld bytes at %p\n", i, lengths[i], segments[i]);
//	    }
	    if (mapped_segments) free(mapped_segments);
	    if (mapped_lengths) free(mapped_lengths);
	    mapped_segments = segments;
	    mapped_lengths = lengths;
	    mapped_count = iovcnt;
	}
    } else {
	/* we have to copy the data.  Have we got a mapped buffer 
	   we can copy into?  If so, reuse. */
	CMbuffer last_write_buffer = ncd->mapped_write_buffer;
	int register_size = 0;
	if (last_write_buffer) register_size = last_write_buffer->size;
	if (register_size < write_size) {
	    if (ncd->outgoing_mapped_region != NULL) {
		ncd->outgoing_mapped_region = NULL;
		NNTI_unregister_memory(&ncd->outgoing_mapped_mr);
		svc->trace_out(ncd->ntd->cm, "CMNNTI unregistering previously mapped region at %p, size %d",
			       ncd->outgoing_mapped_region, ncd->outgoing_mapped_size);
	    }
	    write_buffer = svc->get_data_buffer(ncd->ntd->cm, write_size);
	    register_size = write_buffer->size;
	    svc->trace_out(ncd->ntd->cm, "CMNNTI getting a new copy region %p, size %d",
			   write_buffer->buffer, register_size);
	    data = write_buffer->buffer;
	    err = NNTI_register_memory(&ncd->ntd->trans_hdl, data, register_size, 1, NNTI_GET_SRC, &mr);
	} else {
	    write_buffer = last_write_buffer;
	    svc->trace_out(ncd->ntd->cm, "CMNNTI reusing old copy region %p, size %d",
			   write_buffer->buffer, write_size);
	}
	write_size = 0;
	data = write_buffer->buffer;
	for(i=0; i<iovcnt; i++) {
	    memcpy(&data[write_size], iov[i].iov_base, iov[i].iov_len);
	    write_size += iov[i].iov_len;
	}
    }	    
    if (err != NNTI_OK) {
	printf ("  CMNNTI: NNTI_register_memory() for message returned non-zero: %d %s\n",
		err, NNTI_ERROR_STRING(err));
    }
    ncd->outgoing_mapped_region = data;
    ncd->outgoing_mapped_size = register_size;
    ncd->outgoing_mapped_mr = mr;

    buffer_pack(&ncd->ntd->trans_hdl, &mr, &packed, &packed_size);

    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
	svc->trace_out(ncd->ntd->cm, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");        
	exit(1);
    }

    int offset_of_packed_src_buf = ((int) (((char *) (&(((struct client_message *)NULL)->pull.packed_src_buf))) - ((char *) NULL)));
    h = get_control_message_buffer(ncd, &m, offset_of_packed_src_buf + packed_size + sizeof(packed_size));
    ptr = &m->pull.packed_src_buf[0];
    memcpy(ptr, &packed_size, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(ptr, packed, packed_size);
    ptr += packed_size;

    m->message_type = CMNNTI_PULL_REQUEST;
    m->pull.size = write_size;
    m->pull.reuse_prior_mapping = can_reuse_mapping;
    m->pull.msg_info = local_message_info;
    local_message_info->send_id = send_request++;
    local_message_info->size = write_size;
    local_message_info->mr = mr;
    local_message_info->write_buffer = write_buffer;
    ncd->mapped_write_buffer = write_buffer;
    local_message_info->notify_func = notify_func;
    local_message_info->notify_client_data = notify_client_data;
    svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET sending pull request for %d bytes", write_size);
    svc->set_pending_write(ncd->conn);
    if (send_control_message(h) == 0) return 0;

//    buffer_free(&ncd->ntd->trans_hdl, packed);

    return 1;
}

/* a simple rate limiting scheduler */
struct rate_limit_sched_data {
    int max_ongoing_reqs;
    int ongoing_req_count;
};

/*
 * Rate limiting scheduler. it only allows fetching the request iff the outstanding RDMA Gets is within
 * the limit.
 * Return 1 for ok to issue RDMA Get; 0 for delay the Get operation; -1 for error.
 */
int nnti_rate_limit_scheduler(struct pull_request *request, 
                              struct pull_sched_context *sched_context, 
                              void *callback_data
                             )
{
    struct rate_limit_sched_data *sched_data = (struct rate_limit_sched_data *) callback_data;
    if (sched_data->ongoing_req_count == sched_data->max_ongoing_reqs) {
        return 0;
    } else if (sched_data->ongoing_req_count < sched_data->max_ongoing_reqs) {
        sched_data->ongoing_req_count ++;
        return 1;
    } else {
        return -1;
    }
}

/*
 * This callback is called when a pull request is completed successfully
 */
void nnti_rate_limit_on_completion(struct pull_request *request, 
                                   struct pull_sched_context *sched_context, 
                                   void *callback_data
                                  )
{
    struct rate_limit_sched_data *sched_data = (struct rate_limit_sched_data *) callback_data;
    sched_data->ongoing_req_count --;
}

void
add_ncd_to_request_list(nnti_conn_data_ptr ncd)
{
    int i;
    for (i=0; i< ncd->ntd->request_size; i++) {
	if (ncd->ntd->request_list[i] == NULL) {
	    ncd->ntd->request_list[i] = ncd;
	    pthread_cond_signal(&(ncd->ntd->cond));
	    return;
	}
    }
    /* above found no empty slots */
    ncd->ntd->request_size ++;
    ncd->ntd->request_list = realloc(ncd->ntd->request_list, sizeof(ncd->ntd->request_list[0]) * ncd->ntd->request_size);
    ncd->ntd->request_list[ncd->ntd->request_size-1] = ncd;
    pthread_cond_signal(&(ncd->ntd->cond));
}

static void
free_mapped_buffer(void* ptr)
{
    int err;
    err = NNTI_free(ptr);
    if (err != NNTI_OK) {
	printf ("  CMNNTI: NNTI_free() for client returned non-zero: %d %s\n",
		err, NNTI_ERROR_STRING(err));
    }
}

/*
 * issue RDMA Get for the request. Return 0 for success and -1 for error.
 */
int perform_pull_request_message(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
               struct pull_request *request)
{
    int err;
    int offset = 0;
    int pullsize = request->size;
    char *data;
    NNTI_status_t               status;
    NNTI_buffer_t get_src_mr;
    void *packed = NULL;
    void *ptr = NULL;
    uint64_t packed_size;

    svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET Received pull request, pulling %d bytes", request->size);

    memset(&get_src_mr, 0, sizeof(get_src_mr));
    ptr = &request->packed_src_buf[0];

    memcpy(&packed_size, ptr, sizeof(packed_size));
    ptr += sizeof(packed_size);
    packed = malloc(packed_size);
    memcpy(packed, ptr, packed_size);

    buffer_unpack(&ncd->ntd->trans_hdl, ptr, packed_size, &get_src_mr);

    if (ncd->read_buffer) {
	svc->trace_out(ncd->ntd->cm, "Considering reuse of buffer %p, ref_count %d\n", ncd->read_buffer, ncd->read_buffer->ref_count);
    }
    if ((request->reuse_prior_mapping) && ncd->read_buffer && (ncd->read_buffer->ref_count == 1)) {
	/* no need to reregister!  We'll reuse!*/
        svc->trace_out(ncd->ntd->cm, "CMNNTI reusing already mapped region at %p, size %d",
		       ncd->incoming_mapped_region, ncd->incoming_mapped_size);
        svc->trace_out(ncd->ntd->cm, "CMNNTI alloced region at %p, size %d",
		       ncd->read_buffer->buffer, ncd->read_buffer->size);

	data = ncd->read_buffer->buffer;
    } else {
	if (ncd->read_buffer != NULL) {
	    svc->return_data_buffer(trans->cm, ncd->read_buffer);
	    svc->trace_out(ncd->ntd->cm, "CMNNTI freeing previously mapped region at %p, size %d",
			   ncd->incoming_mapped_region, ncd->incoming_mapped_size);
	    ncd->read_buffer = NULL;
	}

	err = NNTI_alloc(&ncd->ntd->trans_hdl, request->size, 1, NNTI_GET_DST, &ncd->mr_pull);

	data = NNTI_BUFFER_C_POINTER(&ncd->mr_pull);
	ncd->read_buffer = svc->create_data_and_link_buffer(ncd->ntd->cm, data, request->size);
	ncd->read_buffer->return_callback = free_mapped_buffer;
	ncd->read_buffer->return_callback_data = &ncd->mr_pull;
        svc->trace_out(ncd->ntd->cm, "CMNNTI alloced region at %p, size %d",
		       data, ncd->read_buffer->size);

	if (err != NNTI_OK) {
	    printf ("  CMNNTI: NNTI_register_memory() for client returned non-zero: %d %s\n",
		    err, NNTI_ERROR_STRING(err));
            return -1;
	}
	ncd->incoming_mapped_size = ncd->read_buffer->size;
        memcpy(&ncd->incoming_mapped_region, &get_src_mr, sizeof(NNTI_buffer_t));
    }

    err = NNTI_get (&get_src_mr,
		    offset,  // get from this remote buffer+offset
		    pullsize,      // this amount of data
		    &ncd->mr_pull,
		    offset, &ncd->wr_pull); // into this buffer+offset
    
    if (err != NNTI_OK) {
	printf ("  THREAD: Error: NNTI_get() for client returned non-zero: %d %s\n",
		err, NNTI_ERROR_STRING(err));
    //    conns[which%nc].status = 1; // failed status
        return -1;
    }

    if (ncd->ntd->immediate_pull_wait) {
	int timeout = 500;
	err = NNTI_ETIMEDOUT;
	DROP_CM_LOCK(svc, trans->cm);
	while (err == NNTI_ETIMEDOUT ) {
	    err = NNTI_wait(&ncd->wr_pull, timeout, &status);
	    timeout ++;
	}
	ACQUIRE_CM_LOCK(svc, trans->cm);

	if (err != NNTI_OK) {
	    fprintf (stderr, "  THREAD: Error: pull from client failed. NNTI_wait returned: %d %s, wait_status.result = %ds\n",err, NNTI_ERROR_STRING(err), status.result);
	    /* bad shit */
	} else {
	    // completed a pull here
	    send_handle h;
	    struct client_message *r;
	    h = get_control_message_buffer(ncd, &r, sizeof(*r));
	    r->message_type = CMNNTI_PULL_COMPLETE;
	    r->pull_complete.msg_info = request->msg_info;
	    svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET done with pull, returning control message, type %d, local_info %p", r->message_type, r->pull_complete.msg_info);
	    if (send_control_message(h) == 0) {
		svc->trace_out(ncd->ntd->cm, "--- control message send failed!");
	    }	  
	    
	    /* kick this upstairs */
	    if (!getenv("HANDLE_MESSAGES_IN_NNTI_THREAD")) {
		svc->add_buffer_to_pending_queue(trans->cm, ncd->conn, ncd->read_buffer, request->size);
	    } else {
		ncd->read_buf_len = request->size;
		trans->data_available(trans, ncd->conn);
	    }
//	    if (!ncd->ntd->cache_maps) {
//		NNTI_unregister_memory(&ncd->mr_pull);
//		ncd->read_buffer = NULL;
//	    }
	}
	return 0;
    } else {
	// instead of waiting for current Get operation to finish, go on to pull for any other requests
	// add destination buffer to ntd->request_list to wait
	ncd->pending_request_msg_info = request->msg_info;
	ncd->pending_request_size = request->size;
	add_ncd_to_request_list(ncd);
	return 0;
    }
}


static void
handle_pull_request_message(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
		       struct client_message *m)
{
    // put the current message into pull request queue
    int32_t size;
    memcpy(&size, &m->pull.packed_src_buf[0], sizeof(size));
    int pull_rqst_size = sizeof(struct pull_request_queue) + size + 4;
    int copy_size = sizeof(struct pull_request) + size + 4;
    struct pull_request_queue *request = (struct pull_request_queue *)
        svc->malloc_func (pull_rqst_size);
    
    if (!request) {
        svc->trace_out(ncd->ntd->cm, "CMNNTI: Error: cannot allocate memory %s:%d", __FUNCTION__, __LINE__);
        return;
    }
    memcpy (&request->request, &m->pull, copy_size);
    request->state = GET_SCHED_QUEUED;
    request->num_errors = 0;

#define NO_SCHEDULER
#ifdef NO_SCHEDULER
    perform_pull_request_message(ncd, svc, trans, &request->request);
    free(request);
#else
    /* max number of errors for scheduling a request */
    static int max_num_errors_scheduling = 10;


    request->next = ncd->ntd->pull_req_queue;
    ncd->ntd->pull_req_queue = request;
    // go over the request queue: schedule and pull requests
    struct pull_request_queue *req = ncd->ntd->pull_req_queue;
    struct pull_request_queue *prev_req = NULL;
    while(req) {
        // first consult the scheduler if we should pull this request.
        int should_get = (* ncd->ntd->nnti_pull_scheduler)(&req->request, &req->context, ncd->ntd->nnti_pull_sched_data);
        if (should_get == 0) { // do not pull right now
             prev_req = req;
             req = req->next;
             continue;
	} else if (should_get == -1) { // error took place
            request->state = GET_SCHED_ERROR;
            request->num_errors ++;
            if (request->num_errors == max_num_errors_scheduling) {
                svc->trace_out(ncd->ntd->cm, "CMNNTI: Errror in scheduling pull %s:%d", __FUNCTION__, __LINE__);
                // delete the request
                if (prev_req == NULL) {
                    ncd->ntd->pull_req_queue->next = req->next;
                } else {
                    prev_req->next = req->next;
                }
                // TODO: send a notification to sender side
            }
            prev_req = req;
            req = req->next;
            continue;
        }

        // issue RDMA GET for the request
        int rc = perform_pull_request_message(ncd, svc, trans, &req->request);
        if (rc == 0) { // success
            // move the request to ongoing_req_queue
            struct pull_request_queue *temp = req->next;
            if (prev_req == NULL) {
                ncd->ntd->pull_req_queue->next = req->next;
            } else {
                prev_req->next = req->next;
            }
            req->state = GET_SCHED_GET_ISSUED;
            req->next = ncd->ntd->ongoing_req_queue;
            ncd->ntd->ongoing_req_queue = req;
            req = temp;
        } else {
            request->state = GET_SCHED_ERROR;
            request->num_errors ++;
            if (request->num_errors == max_num_errors_scheduling) {
                svc->trace_out(ncd->ntd->cm, "CMNNTI: Errror in scheduling pull %s:%d", __FUNCTION__, __LINE__);
                // delete the request
                if (prev_req == NULL) {
                    ncd->ntd->pull_req_queue->next = req->next;
                } else {
                    prev_req->next = req->next;
                }
                // TODO: send a notification to sender side
            }
            prev_req = req;
            req = req->next;
        }
    }
#endif
}

static void
handle_pull_complete_message(nnti_conn_data_ptr ncd, CMtrans_services svc, transport_entry trans,
		       struct client_message *m)
{
    nnti_message_info local_message_info = m->pull_complete.msg_info;
    svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET received pull complete message, freeing resources, unblocking any pending writes");
    if (local_message_info->notify_func) {
	(local_message_info->notify_func)(local_message_info->notify_client_data);
    }
    if (!ncd->ntd->cache_maps) {
	/* unregister IF we're not caching, OR the memory wasn't ours */
	ncd->outgoing_mapped_region = NULL;
	NNTI_unregister_memory(&local_message_info->mr);
	if (local_message_info->write_buffer) {
	    /* release message buffer if there is one */
	    svc->return_data_buffer(ncd->ntd->cm, local_message_info->write_buffer);
	    ncd->mapped_write_buffer = NULL;
	    local_message_info->write_buffer = NULL;
	}

    }
    svc->wake_any_pending_write(ncd->conn);
    free(local_message_info);
}

extern int
libcmnnti_LTX_writev_complete_notify_func(CMtrans_services svc, 
					  nnti_conn_data_ptr ncd,
					  void *iovs,
					  int iovcnt,
					  attr_list attrs,
					  CMcompletion_notify_func notify_func,
					  void *notify_client_data)
{
    struct iovec * iov = (struct iovec*) iovs;
    int size= 0, i= 0;
    int client_header_size = ((int) (((char *) (&(((struct client_message *)NULL)->pig.payload[0]))) - ((char *) NULL)));
    for(i=0; i<iovcnt; i++) size+= iov[i].iov_len;

    if (size <= ncd->piggyback_size_max) {
        send_handle h;
	struct client_message *m;
	h = get_control_message_buffer(ncd, &m, size + client_header_size);
	m->message_type = CMNNTI_PIGGYBACK;
	m->pig.size = size;
	size = 0;
	for(i=0; i<iovcnt; i++) {
	  memcpy(&m->pig.payload[size], iov[i].iov_base, iov[i].iov_len);
	  size += iov[i].iov_len;
	}
        svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET outbound piggybacking %d bytes of data in control message", size);
	if (send_control_message(h) == 0) return 0;
    } else {
	svc->trace_out(ncd->ntd->cm, "CMNNTI/ENET outbound size is %d > %d, doing pull request", size, ncd->piggyback_size_max);
        if (copy_full_buffer_and_send_pull_request(svc, ncd, iov, iovcnt, attrs, notify_func, notify_client_data) == 0)
	    return 0;
    }
    return iovcnt;
}

extern int
libcmnnti_LTX_writev_func(CMtrans_services svc, nnti_conn_data_ptr ncd, void *iovs, 
			int iovcnt, attr_list attrs)
{
    return libcmnnti_LTX_writev_complete_notify_func(svc, ncd, iovs, iovcnt, 
						     attrs, NULL, NULL);
}

static void
free_nnti_data(CManager cm, void *ntdv)
{
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) ntdv;
    CMtrans_services svc = ntd->svc;
    ntd->shutdown_listen_thread = 1;
    DROP_CM_LOCK(svc, ntd->cm);
    pthread_join(ntd->listen_thread, NULL);
    ACQUIRE_CM_LOCK(svc, ntd->cm);
    NNTI_unregister_memory(&ntd->mr_recvs);
    ntd->shutdown_listen_thread = 0;
    svc->free_func(ntd->nnti_pull_sched_data);
    if (ntd->nnti_pull_sched_data != ntd->nnti_pull_completion_data) {
        svc->free_func(ntd->nnti_pull_completion_data);
    }
    svc->free_func(ntd);
}


extern void *
libcmnnti_LTX_initialize(cm, svc)
CManager cm;
CMtrans_services svc;
{
    static int atom_init = 0;
    nnti_transport_data_ptr nnti_data;
    svc->trace_out(cm, "Initialize CMNnti transport");
    if (getenv("NNTI_LOGGING")) {   
        extern int logger_init(const int debug_level, const char *file);
        char nnti_log_filename[256];
	sprintf(nnti_log_filename, "nnti_log_%x", getpid());
        logger_init(5, nnti_log_filename);
    }
    if (atom_init == 0) {
	CM_NNTI_PORT = attr_atom_from_string("NNTI_PORT");
	CM_NNTI_ADDR = attr_atom_from_string("NNTI_ADDR");
	CM_NNTI_PARAMS = attr_atom_from_string("NNTI_PARAMS");
	CM_NNTI_ENET_CONTROL = attr_atom_from_string("NNTI_ENET_CONTROL");
	CM_PEER_LISTEN_PORT = attr_atom_from_string("PEER_LISTEN_PORT");
	CM_NNTI_IMMEDIATE_PULL_WAIT = attr_atom_from_string("NNTI_IMMEDIATE_PULL_WAIT");
	CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
	CM_ENET_ADDR = attr_atom_from_string("CM_ENET_ADDR");
	CM_ENET_PORT = attr_atom_from_string("CM_ENET_PORT");
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	CM_NNTI_TRANSPORT = attr_atom_from_string("CM_NNTI_TRANSPORT");
	CM_PEER_IP = attr_atom_from_string("PEER_IP");
	CM_TRANSPORT_RELIABLE = attr_atom_from_string("CM_TRANSPORT_RELIABLE");
	atom_init++;
    }
    nnti_data = svc->malloc_func(sizeof(struct nnti_transport_data));
    memset(nnti_data, 0, sizeof(struct nnti_transport_data));
    nnti_data->cm = cm;
    nnti_data->svc = svc;
    nnti_data->socket_fd = -1;
    nnti_data->self_ip = 0;
    nnti_data->self_port = -1;
    nnti_data->connections = NULL;
    nnti_data->listen_attrs = NULL;
    nnti_data->enet_listen_port = -1;
    nnti_data->characteristics = create_attr_list();
    nnti_data->cache_maps = 1;
    add_int_attr(nnti_data->characteristics, CM_TRANSPORT_RELIABLE, 1);
    svc->add_shutdown_task(cm, free_nnti_data, (void *) nnti_data, FREE_TASK);

    nnti_data->pull_req_queue = NULL;
    nnti_data->ongoing_req_queue = NULL;

    /* hardcode to use rate limiting scheduler */
    struct rate_limit_sched_data *sched_data = (struct rate_limit_sched_data *)
        svc->malloc_func(sizeof(struct rate_limit_sched_data));
    sched_data->max_ongoing_reqs = 5; /* get this from outside environement */
    sched_data->ongoing_req_count = 0;
    nnti_data->nnti_pull_scheduler = nnti_rate_limit_scheduler;
    nnti_data->nnti_pull_sched_data = sched_data;
    nnti_data->nnti_pull_completion = nnti_rate_limit_on_completion;
    nnti_data->nnti_pull_completion_data = sched_data;

    return (void *) nnti_data;
}

extern void
libcmnnti_LTX_shutdown_conn(svc, ncd)
CMtrans_services svc;
nnti_conn_data_ptr ncd;
{
    unlink_connection(ncd->ntd, ncd);
    NNTI_unregister_memory(&ncd->mr_send);
    free(ncd->peer_hostname);
    //    free_attr_list(ncd->attrs);
    free(ncd);
}

extern attr_list
libcmnnti_LTX_get_transport_characteristics(transport_entry trans, CMtrans_services svc,
					       void* ntdv)
{
    nnti_transport_data_ptr ntd = (nnti_transport_data_ptr) ntdv;
    return ntd->characteristics;
}


extern transport_entry
cmnnti_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport = svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("nnti");
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)libcmnnti_LTX_initialize;
    transport->listen = (CMTransport_listen_func)libcmnnti_LTX_non_blocking_listen;
    transport->initiate_conn = (CMConnection(*)())libcmnnti_LTX_initiate_conn;
    transport->self_check = (int(*)())libcmnnti_LTX_self_check;
    transport->connection_eq = (int(*)())libcmnnti_LTX_connection_eq;
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)libcmnnti_LTX_shutdown_conn;
    transport->read_block_func = (CMTransport_read_block_func)libcmnnti_LTX_read_block_func;
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)NULL;
    transport->writev_func = (CMTransport_writev_func)libcmnnti_LTX_writev_func;
    transport->writev_complete_notify_func = (CMTransport_writev_complete_notify_func)libcmnnti_LTX_writev_complete_notify_func;
    transport->get_transport_characteristics = (CMTransport_get_transport_characteristics) libcmnnti_LTX_get_transport_characteristics;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}
