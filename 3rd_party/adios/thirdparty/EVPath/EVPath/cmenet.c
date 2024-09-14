/***** Includes *****/
#include "config.h"
#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 1418)
#endif

#undef NDEBUG
#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <windows.h>
#include <process.h>
#include <time.h>
#define getpid()	_getpid()
#else
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#include <sys/socket.h>
#ifdef HAVE_SYS_SOCKIO_H
#include <sys/sockio.h>
#endif
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#endif
#ifdef HAVE_HOSTLIB_H
#include "hostLib.h"
#endif
#ifdef HAVE_STREAMS_UN_H
#include <streams/un.h>
#endif
#include <netinet/in.h>
#include <arpa/inet.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#endif
#include <stdio.h>
#include <fcntl.h>
#ifndef HAVE_WINDOWS_H
#include <net/if.h>
#include <sys/ioctl.h>
#include <errno.h>
#else
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#include <ws2tcpip.h>
#endif
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

#include <atl.h>
#include "evpath.h"
#include "cm_transport.h"
#include "ev_select.h"

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

#ifdef USE_ZPL_ENET
#define ENET_IMPLEMENTATION
#define USE_IPV6
#define MAX_CLIENTS 4095
#include <netinet/in.h>
#include <zpl-enet/include/enet.h>
    /*  extra function to access the UDP socket FD */
    ENET_API enet_uint32 enet_host_get_sock_fd(ENetHost *);
    enet_uint32 enet_host_get_sock_fd(ENetHost *host) {
        return host->socket;
    }

extern void ZPLENETdummy() {  // for warning suppression
     (void) enet_initialize_with_callbacks(0, NULL);
     (void) enet_deinitialize() ;
     (void) enet_linked_version() ;
     (void) enet_socket_listen(0, 0) ;
     (void) enet_socket_accept(0, NULL) ;
     (void) enet_socket_connect(0, NULL) ;
     (void) enet_socket_get_option(0, (ENetSocketOption) 1, NULL) ;
     (void) enet_socket_shutdown(0, (ENetSocketShutdown) 1) ;
     (void) enet_socketset_select(0, NULL, NULL, 0) ;
     (void) enet_address_get_host(NULL, NULL, 0 ) ;
     (void) enet_host_get_peers_count( NULL) ;
     (void) enet_host_get_packets_sent(NULL) ;
     (void) enet_host_get_packets_received(NULL) ;
     (void) enet_host_get_bytes_sent(NULL) ;
     (void) enet_host_get_bytes_received(NULL) ;
     (void) enet_host_get_received_data(NULL, NULL) ;
     (void) enet_host_get_mtu(NULL) ;
     (void) enet_peer_get_id(NULL) ;
     (void) enet_peer_get_ip(NULL, NULL, 0) ;
     (void) enet_peer_get_port(NULL) ;
     (void) enet_peer_get_rtt(NULL) ;
     (void) enet_peer_get_packets_sent(NULL) ;
     (void) enet_peer_get_packets_lost(NULL) ;
     (void) enet_peer_get_bytes_sent(NULL) ;
     (void) enet_peer_get_bytes_received(NULL) ;
     (void) enet_peer_get_state(NULL) ;
     (void) enet_peer_get_data(NULL) ;
     (void) enet_peer_set_data(NULL, NULL) ;
     (void) enet_packet_get_data(NULL) ;
     (void) enet_packet_get_length(NULL) ;
     (void) enet_packet_set_free_callback(NULL, NULL) ;
     (void) enet_packet_create_offset(NULL, 0, 0, 0) ;
     (void) enet_crc32(NULL, 0) ;
     (void) enet_host_check_events(NULL, NULL);
     (void) enet_host_send_raw(NULL, NULL, NULL, 0) ;
     (void) enet_host_send_raw_ex(NULL, NULL, NULL, 0, 0) ;
     (void) enet_host_set_intercept(NULL, NULL) ;
     (void) enet_host_broadcast(NULL, 0, NULL) ;
     (void) enet_host_compress(NULL, NULL);
     (void) enet_host_channel_limit(NULL, 0);
     (void) enet_host_bandwidth_limit(NULL, 0,0);
     (void) enet_peer_ping_interval(NULL, 0) ;
     (void) enet_peer_disconnect_now(NULL, 0) ;
     (void) enet_peer_disconnect_later(NULL, 0) ;
     (void) enet_peer_throttle_configure(NULL, 0, 0, 0) ;
}

#define TPORT "CMZplEnet"
#define TRANSPORT_STRING "zplenet"
#define INTERFACE_NAME(NAME) libcmzplenet_LTX_ ## NAME
#else
#define TPORT "CMEnet"
#define MAX_CLIENTS 0
#define TRANSPORT_STRING "enet"
#define INTERFACE_NAME(NAME) libcmenet_LTX_ ## NAME
#include <enet/enet.h>
#endif
#ifndef _MSC_VER
#include <arpa/inet.h>
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
 
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif



typedef struct _queued_data {
    struct _queued_data *next;
    struct enet_connection_data *econn_d;
    ENetPacket *packet;
} *queued_data;

typedef struct func_list_item {
    select_list_func func;
    void *arg1;
    void *arg2;
} FunctionListElement;

typedef struct enet_client_data {
    CManager cm;
    char *hostname;
    int listen_port;
    CMtrans_services svc;
    ENetHost *server;
    queued_data pending_data;
    int wake_write_fd;
    int wake_read_fd;
    enet_uint32 last_host_service_zero_return;
    CMTaskHandle periodic_handle;
    thr_mutex_t enet_lock;
    int enet_locked;
    struct enet_connection_data *pending_connections;
} *enet_client_data_ptr;

typedef struct enet_connection_data {
    char *remote_host;
#ifndef USE_IPV6
    int remote_IP;     /* in host byte order */
#else
    struct in6_addr remote_IP;
    int remote_IPv4;     /* in host byte order */
#endif
    int remote_contact_port;
    ENetPeer *peer;
    CMbuffer read_buffer;
    int read_buffer_len;
    ENetPacket *packet;
    enet_client_data_ptr ecd;
    CMConnection conn;
    attr_list conn_attr_list;
    int connect_condition;
    struct enet_connection_data *next_pending;
} *enet_conn_data_ptr;

static atom_t CM_PEER_IP = -1;
static atom_t CM_PEER_LISTEN_PORT = -1;
static atom_t CM_NETWORK_POSTFIX = -1;
static atom_t CM_ENET_PORT = -1;
static atom_t CM_ENET_HOSTNAME = -1;
static atom_t CM_ENET_ADDR = -1;
static atom_t CM_ENET_CONN_TIMEOUT = -1;
static atom_t CM_ENET_CONN_REUSE = -1;
static atom_t CM_TRANSPORT = -1;

static enet_uint32 enet_host_service_warn_interval = 0;

#ifdef __cplusplus
extern "C"
#else
extern
#endif
attr_list
INTERFACE_NAME(non_blocking_listen)(CManager cm, CMtrans_services svc,
                                    transport_entry trans, attr_list listen_info);

#define ENETlock(lock) IntENET_lock(lock, __FILE__, __LINE__)
#define ENETunlock(lock) IntENET_unlock(lock, __FILE__, __LINE__)
#define gettid() pthread_self()

static void
IntENET_lock(enet_client_data_ptr ecd, char *file, int line)
{
//    if (file) printf("(PID %lx, TID %lx) Trying ENET Lock at %s, line %d\n", (long) getpid(), (long)gettid(), file, line);
    thr_mutex_lock(ecd->enet_lock);
//    if (file) printf("GOT ENET Lock at %s, line %d\n", file, line);
    ecd->enet_locked++;
}

static void
IntENET_unlock(enet_client_data_ptr ecd, char *file, int line)
{
//    if (file) printf("(PID %lx, TID %lx) ENET Unlock at %s, line %d\n", (long) getpid(), (long)gettid(), file, line);
    ecd->enet_locked--;
    thr_mutex_unlock(ecd->enet_lock);
}

static int
check_host(char *hostname, void *sin_addr)
{
    (void)hostname; (void)sin_addr;
#ifdef HAS_STRUCT_HOSTEN
    struct hostent *host_addr;
    host_addr = gethostbyname(hostname);
    if (host_addr == NULL) {
	struct in_addr addr;
	if (inet_aton(hostname, &addr) == 0) {
	    /* 
	     *  not translatable as a hostname or 
	     * as a dot-style string IP address
	     */
	    return 0;
	}
	assert(sizeof(int) == sizeof(struct in_addr));
	*((int *) sin_addr) = *((int*) &addr);
    } else {
	memcpy(sin_addr, host_addr->h_addr, host_addr->h_length);
    }
    return 1;
#endif
    printf("Check host called, unimplemented\n");
    return 0;
}

static enet_conn_data_ptr 
create_enet_conn_data(CMtrans_services svc)
{
    enet_conn_data_ptr enet_conn_data = (enet_conn_data_ptr) svc->malloc_func(sizeof(struct enet_connection_data));
    enet_conn_data->remote_host = NULL;
    enet_conn_data->remote_contact_port = -1;
    enet_conn_data->read_buffer = NULL;
    enet_conn_data->read_buffer_len = 1;
    return enet_conn_data;
}

static void *
enet_accept_conn(enet_client_data_ptr ecd, transport_entry trans, 
		 ENetAddress *address);

static void free_func(void *packet)
{
    /* Clean up the packet now that we're done using it. */
//    ENETlock(enet_data);
    enet_packet_destroy ((ENetPacket*)packet);
//    ENETunlock(enet_data);
}

static void
handle_packet(CManager cm, CMtrans_services svc, transport_entry trans, enet_conn_data_ptr econn_d, ENetPacket *packet)
{
    CMbuffer cb;
    svc->trace_out(cm, "A packet of length %u was received.\n",
                   (unsigned int) packet->dataLength);
    econn_d->read_buffer_len = (int) packet->dataLength;
    cb = svc->create_data_and_link_buffer(cm, 
                                          packet->data, 
                                          econn_d->read_buffer_len);
    econn_d->read_buffer = cb;
    cb->return_callback = free_func;
    cb->return_callback_data = packet;
    econn_d->packet = packet;

    /* kick this upstairs */
    trans->data_available(trans, econn_d->conn);
    svc->return_data_buffer(trans->cm, cb);
}
	    
static
void
enet_service_network(CManager cm, void *void_trans)
{
    transport_entry trans = (transport_entry) void_trans;
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    CMtrans_services svc = ecd->svc;
    ENetEvent event;
    
    if (!ecd->server) return;
    if (!(CM_LOCKED(svc, ecd->cm))) {
	printf("Enet service network, CManager not locked\n");
    }

    while (ecd->pending_data) {
        svc->trace_out(cm, "ENET Handling pending data\n");
        queued_data entry = ecd->pending_data;
        ecd->pending_data = entry->next;
        handle_packet(cm, svc, trans, (enet_conn_data_ptr) entry->econn_d, entry->packet);
        free(entry);
    }

    while (ecd->server) {
        IntENET_lock(ecd, NULL, 0);
        int ret = enet_host_service (ecd->server, & event, 0);
        if (enet_host_service_warn_interval && 
            (enet_time_get() > (ecd->last_host_service_zero_return + enet_host_service_warn_interval))) {
            fprintf(stderr, "WARNING, time between zero return for enet_host_service = %d msecs\n",
                    enet_time_get() - ecd->last_host_service_zero_return);
        }
        IntENET_unlock(ecd, NULL, 0);
        if (ret <= 0) {
            break;
        }
        switch (event.type) {
	case ENET_EVENT_TYPE_NONE:
	    break;
        case ENET_EVENT_TYPE_CONNECT: {
	    enet_conn_data_ptr enet_connection_data = NULL;
            if (event.peer->data) {
                enet_conn_data_ptr last = NULL;
                enet_connection_data = ecd->pending_connections;
                while (enet_connection_data) {
                    if (enet_connection_data->peer == event.peer) {
                        if (last) {
                            last->next_pending = enet_connection_data->next_pending;
                        } else {
                            ecd->pending_connections = enet_connection_data->next_pending;
                        }                            
                        enet_connection_data->next_pending = NULL;
                        break;
                    }
                    enet_connection_data = enet_connection_data->next_pending;
                }
            }
            if (enet_connection_data) {
                svc->condition_signal(cm, enet_connection_data->connect_condition);
                break;
            }
#ifndef USE_IPV6
            struct in_addr addr;
            addr.s_addr = event.peer->address.host;
	    svc->trace_out(cm, "A new client connected from %s:%u.\n", 
			   inet_ntoa(addr),
			   event.peer->address.port);
#else
            char straddr[INET6_ADDRSTRLEN];
            inet_ntop(AF_INET6, &event.peer->address.host, straddr,
                      sizeof(straddr));
	    svc->trace_out(cm, "A new client connected from %s:%u.\n", 
			   &straddr[0],
			   event.peer->address.port);
            struct in_addr addr;
            if (ntohl(((enet_uint32 *)&event.peer->address.host.s6_addr)[3]) == 0xffff) {;
                addr.s_addr = ((enet_uint32 *)&event.peer->address.host.s6_addr)[3];
                svc->trace_out(cm, "That was IPV4 address %s\n", inet_ntoa(addr));
            }
#endif            

	    enet_connection_data = enet_accept_conn(ecd, trans, &event.peer->address);

            /* Store any relevant client information here. */
            svc->trace_out(cm, "ENET ========   Assigning peer %p has data %p\n", event.peer, enet_connection_data);
            enet_peer_timeout(event.peer, 0, 0, 200);
            event.peer->data = enet_connection_data;
	    ((enet_conn_data_ptr)enet_connection_data)->peer = event.peer;

            break;
	}
        case ENET_EVENT_TYPE_RECEIVE: {
	    enet_conn_data_ptr econn_d = (enet_conn_data_ptr) event.peer->data;
            if (econn_d) {
                handle_packet(cm, svc, trans, econn_d, event.packet);
            } else {
#ifndef USE_IPV6
                struct in_addr addr;
                addr.s_addr = event.peer->address.host;
                svc->trace_out(cm, "ENET  ====== virgin peer, address is %s, port %u.\n", inet_ntoa(addr), event.peer->address.port);
#else
                char straddr[INET6_ADDRSTRLEN];
                inet_ntop(AF_INET6, &event.peer->address.host, straddr,
                          sizeof(straddr));
                svc->trace_out(cm, "ENET  ====== virgin peer, address is %s, port %u.\n", 
                               &straddr[0],
                               event.peer->address.port);
#endif               
                svc->trace_out(cm, "ENET  ====== DISCARDING DATA\n");
            }
            break;
	}           
#ifdef USE_ZPL_ENET
        case ENET_EVENT_TYPE_DISCONNECT_TIMEOUT:
#endif
        case ENET_EVENT_TYPE_DISCONNECT: {
	    enet_conn_data_ptr enet_conn_data = (enet_conn_data_ptr) event.peer->data;
	    svc->trace_out(cm, "Got a disconnect on connection %p\n", event.peer->data);

            enet_conn_data = (enet_conn_data_ptr) event.peer->data;
	    enet_conn_data->read_buffer_len = -1;
            if (enet_conn_data->conn) {
                svc->connection_fail(enet_conn_data->conn);
            }
            break;
        }
        default:
            printf("UNKNOWN EVENT TYPE! %d\n", event.type);
            break;
        }
    }
    ecd->last_host_service_zero_return = enet_time_get();
}

static
void
enet_service_network_lock(CManager cm, void *void_trans)
{
    transport_entry trans = (transport_entry) void_trans;
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    CMtrans_services svc = ecd->svc;
    ACQUIRE_CM_LOCK(svc, cm);
    enet_service_network(cm, void_trans);
    DROP_CM_LOCK(svc, cm);
}

static
void
read_wake_fd_and_service(CManager cm, void *void_trans)
{
    transport_entry trans = (transport_entry) void_trans;
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;

    char buffer;
    int fd = ecd->wake_read_fd;

#ifdef HAVE_WINDOWS_H
    recv(fd, &buffer, 1, 0);
#else
    if (read(fd, &buffer, 1) != 1) {
	perror("wake read failed\n");
    }
#endif

    enet_service_network(cm, void_trans);
}

#ifdef NOTDEF
static
void 
dump_sockaddr(who, sa)
char *who;
struct sockaddr_in *sa;
{
    unsigned char *addr;

    addr = (unsigned char *) &(sa->sin_addr.s_addr);

    printf("%s: family=%d port=%d addr=%d.%d.%d.%d\n",
	   who,
	   ntohs(sa->sin_family),
	   ntohs(sa->sin_port),
	   addr[0], addr[1], addr[2], addr[3]);
}

static
void 
dump_sockinfo(msg, fd)
char *msg;
int fd;
{
    int nl;
    struct sockaddr_in peer, me;

    printf("Dumping sockinfo for fd=%d: %s\n", fd, msg);

    nl = sizeof(me);
    getsockname(fd, (struct sockaddr *) &me, &nl);
    dump_sockaddr("Me", &me);

    nl = sizeof(peer);
    getpeername(fd, (struct sockaddr *) &peer, &nl);
    dump_sockaddr("Peer", &peer);
}

#endif

static int conn_reuse = 1;

/* 
 * Accept enet connection
 */
static void *
enet_accept_conn(enet_client_data_ptr ecd, transport_entry trans, 
		 ENetAddress *address)
{
    CMtrans_services svc = ecd->svc;
    enet_conn_data_ptr enet_conn_data;

    CMConnection conn;
    attr_list conn_attr_list = NULL;;

    enet_conn_data = create_enet_conn_data(svc);
    enet_conn_data->ecd = ecd;
    conn_attr_list = create_attr_list();
    conn = svc->connection_create(trans, enet_conn_data, conn_attr_list);
    enet_conn_data->conn = conn;

#ifndef USE_IPV6
    add_int_attr(conn_attr_list, CM_PEER_IP, ntohl(address->host));
    enet_conn_data->remote_IP = ntohl(address->host); 
#else
    enet_conn_data->remote_IP = address->host;
    enet_conn_data->remote_IPv4 = ntohl(((enet_uint32 *)&address->host.s6_addr)[3]);
    add_int_attr(conn_attr_list, CM_PEER_IP, enet_conn_data->remote_IPv4);  /* remote_IPv4 is in host byte order */
#endif
    if (!conn_reuse) {
        enet_conn_data->remote_contact_port = -1;
    } else {
        enet_conn_data->remote_contact_port = address->port;
    }

    if (enet_conn_data->remote_host != NULL) {
	svc->trace_out(trans->cm, "Accepted ENET RUDP connection from host \"%s\"",
		       enet_conn_data->remote_host);
    } else {
	svc->trace_out(trans->cm, "Accepted ENET RUDP connection from UNKNOWN host");
    }
    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)enet_conn_data->remote_contact_port);
#ifndef USE_IPV6
    struct in_addr addr;
    addr.s_addr = htonl(enet_conn_data->remote_IP);
    svc->trace_out(trans->cm, "Remote host (IP %s) is listening at port %d\n",
        inet_ntoa(addr),
        enet_conn_data->remote_contact_port);
#else
    char straddr[INET6_ADDRSTRLEN];
    inet_ntop(AF_INET6, &enet_conn_data->remote_IP, straddr,
              sizeof(straddr));
    svc->trace_out(trans->cm, "Remote host (IP %s) is listening at port %d\n",
        &straddr[0],
        enet_conn_data->remote_contact_port);
#endif

    free_attr_list(conn_attr_list);

    return enet_conn_data;
}

#ifdef __cplusplus
extern "C"
#else
extern
#endif
void
INTERFACE_NAME(shutdown_conn)(CMtrans_services svc, enet_conn_data_ptr scd)
{
    svc->connection_deref(scd->conn);
    if (scd->remote_host) free(scd->remote_host);
    free(scd);
}


static void *
enet_initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans,
	      attr_list attrs, enet_conn_data_ptr enet_conn_data,
	      attr_list conn_attr_list)
{
    int int_port_num;
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    char *host_name;
    int host_ip = 0;
#ifdef USE_IPV6
    struct in6_addr host_ipv6;
#endif
    struct in_addr sin_addr;
    (void)conn_attr_list;
    int timeout = 200;   /* connection time out default 100 milliseconds */

    if (!(CM_LOCKED(svc, ecd->cm))) {
	printf("Enet service network, CManager not locked in enet_initiate_conn\n");
    }

    if (!query_attr(attrs, CM_ENET_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, TPORT " transport found no CM_ENET_HOSTNAME attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, TPORT " transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_ENET_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_ip)) {
	svc->trace_out(cm, "CMEnet transport found no CM_ENET_ADDR attribute");
	/* wasn't there */
	host_ip = 0;
    } else {
        svc->trace_out(cm, "CMEnet transport connect to host_IP %lx", host_ip);
    }
    /* HOST_IP is in HOST BYTE ORDER */
    if ((host_name == NULL) && (host_ip == 0)) {
	printf("No host no IP\n");
	return 0;
    }

    if (!query_attr(attrs, CM_ENET_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "CMEnet transport found no CM_ENET_PORT attribute");
	return 0;
    } else {
        svc->trace_out(cm, "CMEnet transport connect to port %d", int_port_num);
    }

    if (!query_attr(attrs, CM_ENET_CONN_TIMEOUT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & timeout)) {
	svc->trace_out(cm, "CMEnet transport found no CM_ENET_CONN_TIMEOUT attribute");
    } else {
        svc->trace_out(cm, "CMEnet transport connection timeout set to %d msecs", timeout);
    }
    if (!query_attr(attrs, CM_ENET_CONN_REUSE, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & conn_reuse)) {
	svc->trace_out(cm, "CMEnet transport found no CM_ENET_CONN_REUSE attribute");
    } else {
        svc->trace_out(cm, "CMEnet transport connection reuse set to %d", conn_reuse);
    }

    /* ENET connection, host_name is the machine name */
    ENetAddress address;
    ENetPeer *peer;
    sin_addr.s_addr = htonl(host_ip);

    if (host_name) {
	enet_address_set_host (& address, host_name);
#ifndef USE_IPV6
	sin_addr.s_addr = address.host;
	svc->trace_out(cm, "Attempting ENET RUDP connection, USING host=\"%s\", IP = %s, port %d",
		       host_name == 0 ? "(unknown)" : host_name, 
		       inet_ntoa(sin_addr),
		       int_port_num);
#else
        char straddr[INET6_ADDRSTRLEN];
        inet_ntop(AF_INET6, &address.host, straddr,
                  sizeof(straddr));
	svc->trace_out(cm, "Attempting ENET RUDP connection, USING host=\"%s\", IP = %s, port %d",
		       host_name == 0 ? "(unknown)" : host_name, 
		       &straddr[0],
		       int_port_num);
#endif
    } else {
#ifndef USE_IPV6
	address.host = ntohl(host_ip);
	svc->trace_out(cm, "Attempting ENET RUDP connection, USING IP = %s, port %d",
		       inet_ntoa(sin_addr),
		       int_port_num);
#else
        char straddr[INET6_ADDRSTRLEN];
        *((enet_uint32 *)&host_ipv6.s6_addr) = 0;
        *(((enet_uint32 *)&host_ipv6.s6_addr) + 1) = 0;
        *(((enet_uint32 *)&host_ipv6.s6_addr) + 2) = htonl(0xffff);
        *(((enet_uint32 *)&host_ipv6.s6_addr) + 3) = htonl(host_ip);
        inet_ntop(AF_INET6, &host_ipv6, straddr,
                  sizeof(straddr));

	svc->trace_out(cm, "Attempting ENET RUDP connection, USING host=\"%s\", IP = %s, port %d",
		       host_name == 0 ? "(unknown)" : host_name, 
		       &straddr[0],
		       int_port_num);
        sin_addr.s_addr = htonl(host_ip);
	svc->trace_out(cm, "Attempting ENET RUDP connection, USING IPv4 = %s\n",
		       inet_ntoa(sin_addr));
        memcpy(&address.host, &host_ipv6, sizeof(struct in6_addr));
#endif        
    }
    address.port = (unsigned short) int_port_num;

    if (ecd->server == NULL) {
	attr_list l = INTERFACE_NAME(non_blocking_listen)(cm, svc, trans, NULL);
	if (l) free_attr_list(l);
    }

    /* Initiate the connection, allocating the two channels 0 and 1. */
    ENETlock(ecd);
    peer = enet_host_connect (ecd->server, & address, 1, 0);    
    if (peer == NULL)
    {
       fprintf (stderr, 
                "No available peers for initiating an ENet connection, count at initiation was %d.\n", MAX_CLIENTS);
       exit (EXIT_FAILURE);
    }
    
    enet_peer_timeout(peer, 0, 0, 200);
    ENETunlock(ecd);
    peer->data = enet_conn_data;
    enet_conn_data->remote_host = host_name == NULL ? NULL : strdup(host_name);
#ifndef USE_IPV6
    enet_conn_data->remote_IP = htonl(host_ip);
#else
    memcpy(&enet_conn_data->remote_IP, &host_ipv6, sizeof(enet_conn_data->remote_IP));
    enet_conn_data->remote_IPv4 = htonl(((enet_uint32 *)&host_ipv6.s6_addr)[3]);
#endif
    enet_conn_data->remote_contact_port = int_port_num;
    enet_conn_data->ecd = ecd;
    enet_conn_data->peer = peer;
    peer->data = enet_conn_data;
    svc->trace_out(cm, "ENET ========   On init Assigning peer %p has data %p moving to wait phase\n", peer, enet_conn_data);

    return enet_conn_data;
}

/* 
 * Initiate a ENET RUDP connection with another CM.
 */
#ifdef __cplusplus
extern "C"
#else
extern
#endif
void *
INTERFACE_NAME(initiate_conn_nonblocking)(CManager cm, CMtrans_services svc,
                                     transport_entry trans, attr_list attrs, int connect_condition)
{
    enet_conn_data_ptr enet_conn_data = create_enet_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    enet_conn_data_ptr ret;

    enet_conn_data->conn_attr_list = conn_attr_list;
    enet_conn_data->connect_condition = connect_condition;
    ret = enet_initiate_conn(cm, svc, trans, attrs, enet_conn_data, conn_attr_list);
    if (ret) {
        enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
        ret->next_pending = ecd->pending_connections;
        ecd->pending_connections = ret;
    }
    return (void*)ret;
}

/* 
 * Initiate a ENET RUDP connection with another CM.
 */
#ifdef __cplusplus
extern "C"
#else
extern
#endif
CMConnection
INTERFACE_NAME(finalize_conn_nonblocking)(CManager cm, CMtrans_services svc,
                                          transport_entry trans, void *client_data, int result)
{
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    enet_conn_data_ptr final_conn_data = (enet_conn_data_ptr) client_data;
    enet_conn_data_ptr last = NULL, enet_conn_data = ecd->pending_connections;
    CMConnection conn;
    attr_list conn_attr_list = final_conn_data->conn_attr_list;

    if (!result) {
        while (enet_conn_data) {
            if (enet_conn_data == final_conn_data) {
                if (last) {
                    last->next_pending = enet_conn_data->next_pending;
                } else {
                    ecd->pending_connections = enet_conn_data->next_pending;
                }                            
                enet_conn_data->next_pending = NULL;
                break;
            }
            last = enet_conn_data;
            enet_conn_data = enet_conn_data->next_pending;
        }

        free_attr_list(conn_attr_list);
        free(enet_conn_data);
        return NULL;
    }

    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)final_conn_data->remote_contact_port);
    conn = svc->connection_create(trans, final_conn_data, conn_attr_list);
    final_conn_data->conn = conn;
    free_attr_list(conn_attr_list);
    final_conn_data->conn_attr_list = NULL;
    svc->connection_addref(conn);  /* one ref count went to CM, 
				the other to the user */

    return conn;
}


/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * For enet, this involves checking to see if the host name is the 
 * same as ours and if the CM_ENET_PORT matches the one we are listening on.
 */
#ifdef __cplusplus
extern "C"
#else
extern
#endif
int
INTERFACE_NAME(self_check)(CManager cm, CMtrans_services svc, 
			 transport_entry trans, attr_list attrs)
{

    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    int host_addr;
    int int_port_num;
    char *host_name;
    char my_host_name[256];
    static int IP = 0;   /* always in host byte order */

    get_IP_config(my_host_name, sizeof(host_name), &IP, NULL, NULL, NULL,
		  NULL, svc->trace_out, (void *)cm);

    if (IP == 0) {
	IP = ntohl(INADDR_LOOPBACK);
    }
    if (!query_attr(attrs, CM_ENET_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "CMself check CMEnet transport found no CM_ENET_HOSTNAME attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_ENET_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_addr)) {
	svc->trace_out(cm, "CMself check CMEnet transport found no CM_ENET_ADDR attribute");
	if (host_name == NULL) return 0;
	host_addr = 0;
    }
    if (!query_attr(attrs, CM_ENET_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "CMself check CMEnet transport found no CM_ENET_PORT attribute");
	return 0;
    }
    //get_qual_hostname(my_host_name, sizeof(my_host_name), svc, NULL, NULL);

    if (host_name && (strcmp(host_name, my_host_name) != 0)) {
	svc->trace_out(cm, "CMself check - Hostnames don't match");
	return 0;
    }
    if (host_addr && (IP != host_addr)) {
	svc->trace_out(cm, "CMself check - Host IP addrs don't match, %lx, %lx", IP, host_addr);
	return 0;
    }
    if (int_port_num != ecd->listen_port) {
	svc->trace_out(cm, "CMself check - Ports don't match, %d, %d", int_port_num, ecd->listen_port);
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

#ifdef __cplusplus
extern "C"
#else
extern
#endif
int
INTERFACE_NAME(connection_eq)(CManager cm, CMtrans_services svc,
			    transport_entry trans, attr_list attrs,
			    enet_conn_data_ptr ecd)
{

    int int_port_num;
    int requested_IP = -1;
    char *host_name = NULL;

    (void) trans;
    if (!query_attr(attrs, CM_ENET_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "CMEnet transport found no CM_ENET_HOST attribute");
    }
    if (!query_attr(attrs, CM_ENET_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "Conn Eq CMenet transport found no CM_ENET_PORT attribute");
	return 0;
    }
    if (!query_attr(attrs, CM_ENET_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & requested_IP)) {
	svc->trace_out(cm, "CMENET transport found no CM_ENET_ADDR attribute");
    }
    if (requested_IP == -1) {
	check_host(host_name, (void *) &requested_IP);
	requested_IP = ntohl(requested_IP);
        struct in_addr addr;
        addr.s_addr = htonl(requested_IP);
	svc->trace_out(cm, "IP translation for hostname %s is %s", host_name,
		       inet_ntoa(addr));
    }
    /* requested IP is in host byte order */
    if (ecd->peer->state != ENET_PEER_STATE_CONNECTED) {
        svc->trace_out(cm, "ENET Conn_eq returning FALSE, peer not connected");
        return 0;
    }
    struct in_addr addr1, addr2;
#ifdef USE_IPV6
    addr1.s_addr = htonl(ecd->remote_IPv4);
#else
    addr1.s_addr = htonl(ecd->remote_IP);
#endif
    addr2.s_addr = htonl(requested_IP);
    svc->trace_out(cm, "ENET Conn_eq comparing IP/ports %s/%d and %s/%d",
		   inet_ntoa(addr1), ecd->remote_contact_port,
                   inet_ntoa(addr2), int_port_num);
#ifdef USE_IPV6
    if ((ecd->remote_IPv4 == requested_IP) &&    /* both in host byte order */
#else
    if ((ecd->remote_IP == requested_IP) &&    /* both in host byte order */
#endif
	(ecd->remote_contact_port == int_port_num)) {
	svc->trace_out(cm, "ENET Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "ENET Conn_eq returning FALSE");
    return 0;
}

static attr_list
build_listen_attrs(CManager cm, CMtrans_services svc, enet_client_data_ptr ecd,
		   attr_list listen_info, int int_port_num)
{
    char host_name[256];
    attr_list ret_list;
    int IP;
    int use_hostname = 0;
    
    svc->trace_out(cm, "CMEnet listen succeeded on port %d",
		       int_port_num);
    get_IP_config(host_name, sizeof(host_name), &IP, NULL, NULL,
		  &use_hostname, listen_info, svc->trace_out, (void *)cm);

    ret_list = create_attr_list();

    if (ecd) {
	ecd->hostname = strdup(host_name);
	ecd->listen_port = int_port_num;
    }
    if ((IP != 0) && !use_hostname) {
	add_attr(ret_list, CM_ENET_ADDR, Attr_Int4,
		 (attr_value) (intptr_t)IP);
    }
    if ((getenv("CMEnetsUseHostname") != NULL) || 
	use_hostname) {
	add_attr(ret_list, CM_ENET_HOSTNAME, Attr_String,
		 (attr_value) strdup(host_name));
    } else if (IP == 0) {
        add_int_attr(ret_list, CM_ENET_ADDR, INADDR_LOOPBACK);
    }
    add_attr(ret_list, CM_ENET_PORT, Attr_Int4,
	     (attr_value) (intptr_t)int_port_num);
    
    add_attr(ret_list, CM_TRANSPORT, Attr_String,
	     (attr_value) strdup(TRANSPORT_STRING));
    return ret_list;
}

static void
wake_enet_server_thread(enet_client_data_ptr ecd)
{
    static char buffer = 'W';  /* doesn't matter what we write */
    if (ecd->wake_write_fd != -1) {
#ifdef HAVE_WINDOWS_H
	send(ecd->wake_write_fd, &buffer, 1, 0);
#else
	if (write(ecd->wake_write_fd, &buffer, 1) != 1) {
	    printf("Whoops, wake write failed\n");
	}
#endif
    }
}

/* 
 * Create an IP socket for connection from other CMs
 */
#ifdef __cplusplus
extern "C"
#else
extern
#endif
attr_list
INTERFACE_NAME(non_blocking_listen)(CManager cm, CMtrans_services svc,
				  transport_entry trans, attr_list listen_info)
{
    enet_client_data_ptr ecd = (enet_client_data_ptr) trans->trans_data;
    ENetAddress address;
    ENetHost * server;


    int attr_port_num = 0;
    u_short port_num = 0;

    if (!(CM_LOCKED(svc, cm))) {
	printf("ENET non_blocking listen, CManager not locked\n");
    }
    /* 
     *  Check to see if a bind to a specific port was requested
     */
    if (listen_info != NULL
	&& !query_attr(listen_info, CM_ENET_PORT,
		       NULL, (attr_value *)(intptr_t) & attr_port_num)) {
	port_num = 0;
    } else {
	if (attr_port_num > USHRT_MAX || attr_port_num < 0) {
	    fprintf(stderr, "Requested port number %d is invalid\n", attr_port_num);
	    return NULL;
	}
	port_num = attr_port_num;
    }

    svc->trace_out(cm, "CMEnet begin listen, requested port %d", attr_port_num);

    address.host = ENET_HOST_ANY;

    if (ecd->server != NULL) {
	/* we're already listening */
        if (port_num == 0) {
	    /* not requesting a specific port, return what we have */
	    return build_listen_attrs(cm, svc, NULL, listen_info, ecd->listen_port);
	} else {
	    printf("CMlisten_specific() requesting a specific port follows other Enet operation which initiated listen at another port.  Only one listen allowed, second listen fails.\n");
	    return NULL;
	}
    }
    
    if (port_num != 0) {
	/* Bind the server to the default localhost.     */
	/* A specific host address can be specified by   */
	/* enet_address_set_host (& address, "x.x.x.x"); */

	address.port = port_num;

	svc->trace_out(cm, "CMEnet trying to bind selected port %d", port_num);
        ENETlock(ecd);
	server = enet_host_create (& address /* the address to bind the server host to */, 
                                   MAX_CLIENTS,     /* max 4095 connections */
				   1      /* allow up to 2 channels to be used, 0 and 1 */,
				   0      /* assume any amount of incoming bandwidth */,
				   0      /* assume any amount of outgoing bandwidth */);
        ENETunlock(ecd);
	if (server == NULL) {
	    fprintf (stderr, 
		     "An error occurred while trying to create an ENet server host.\n");
	    return NULL;
	}
    } else {
	int low_bound, high_bound;
	get_IP_config(NULL, 0, NULL, &low_bound, &high_bound,
		      NULL, listen_info, svc->trace_out, (void *)cm);
        if (high_bound == -1) {
            /* unconstrained port */
            address.port = 0;

            svc->trace_out(cm, "CMEnet trying to bind to any available port");
            ENETlock(ecd);
            server = enet_host_create (& address /* the address to bind the server host to */, 
                                       MAX_CLIENTS,     /* max 4095 connections */
                                       1      /* allow up to 2 channels to be used, 0 and 1 */,
                                       0      /* assume any amount of incoming bandwidth */,
                                       0      /* assume any amount of outgoing bandwidth */);
            ENETunlock(ecd);
            if (server == NULL) {
                fprintf (stderr, 
                         "An error occurred while trying to create an ENet server host.\n");
                return NULL;
            }
            address.port = server->address.port;
            svc->trace_out(cm, "CMEnet is listening on port %d\n", address.port);
        } else {
            /* specified port range */
            /* port num is free.  Constrain to range 26000 : 26100 */
            int size;
            int tries;
            srand48(time(NULL) + getpid());

        restart:
            size = high_bound - low_bound;
            tries = 10;
            while (tries > 0) {
                int target = low_bound + (int)(size * drand48());
                address.port = target;
                
                svc->trace_out(cm, "CMEnet trying to bind port %d", target);
                
                ENETlock(ecd);
                server = enet_host_create (& address /* the address to bind the server host to */, 
                                           MAX_CLIENTS     /* 0 means dynamic alloc clients and/or outgoing connnections */,
                                           1      /* allow up to 2 channels to be used, 0 and 1 */,
                                           0      /* assume any amount of incoming bandwidth */,
                                           0      /* assume any amount of outgoing bandwidth */);
                ENETunlock(ecd);
                tries--;
                if (server != NULL) tries = 0;
                if (tries == 5) {
                    /* try reseeding in case we're in sync with another process */
                    srand48(time(NULL) + getpid());
                }
            }
            if (server == NULL) {
                high_bound += 100;
                goto restart;
            }
        }
    }
    ecd->server = server;
    svc->fd_add_select(cm, enet_host_get_sock_fd (server), 
		       (select_list_func) enet_service_network, (void*)cm, (void*)trans);

    ecd->periodic_handle = svc->add_periodic_task(cm, 0, 100, (CMPollFunc) enet_service_network_lock, (void*)trans);

    svc->trace_out(ecd->cm, "CMENET Adding read_wake_fd as action on fd %d",
		   ecd->wake_read_fd);

    svc->fd_add_select(cm, ecd->wake_read_fd, (select_list_func)read_wake_fd_and_service, 
                       (void*)cm, (void*)trans);

    return build_listen_attrs(cm, svc, ecd, listen_info, address.port);
}

#if defined(HAVE_WINDOWS_H) && !defined(NEED_IOVEC_DEFINE)
#define NEED_IOVEC_DEFINE
#endif

#ifdef NEED_IOVEC_DEFINE
struct iovec {
    void *iov_base;
    size_t iov_len;
};

#endif

#ifdef __cplusplus
extern "C"
#else
extern
#endif
void *
INTERFACE_NAME(read_block_func)(CMtrans_services svc,
                                enet_conn_data_ptr conn_data, ssize_t *actual_len,
                                ssize_t *offset_ptr)
{
    CMbuffer cb;

    if (conn_data->read_buffer_len == -1) return NULL;

    *actual_len = conn_data->read_buffer_len;
    *offset_ptr = 0;
    cb = conn_data->read_buffer;
    conn_data->read_buffer_len = 0;
    conn_data->read_buffer = NULL;
    return cb;
}

#ifdef CURRENT_UTC_TIME_NEEDED
static
void current_utc_time(struct timespec *ts)
{
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, ts);
#endif
}
#endif

#ifdef TIME_DIFF_USED
static struct timespec time_diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
	temp.tv_sec = end.tv_sec-start.tv_sec-1;
	temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
	temp.tv_sec = end.tv_sec-start.tv_sec;
	temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}
#endif

#ifdef __cplusplus
extern "C"
#else
extern
#endif
int
INTERFACE_NAME(writev_func)(CMtrans_services svc, enet_conn_data_ptr ecd,
			  struct iovec *iov, size_t iovcnt, attr_list attrs)
{
    size_t i;
    size_t length = 0;

    (void) attrs;
    for (i = 0; i < iovcnt; i++) {
	length += iov[i].iov_len;
    }

    svc->trace_out(ecd->ecd->cm, "CMENET vector write of %d bytes on peer %p",
		   length, ecd->peer);

   /* Create a reliable packet of the right size */
    if (!(CM_LOCKED(svc, ecd->ecd->cm))) {
	printf("ENET writev, CManager not locked\n");
    }
    ENETlock(ecd->ecd);
    ENetPacket * packet = enet_packet_create (NULL, length, 
					      ENET_PACKET_FLAG_RELIABLE);
    ENETunlock(ecd->ecd);

    length = 0;
    /* copy in the data */
    for (i = 0; i < iovcnt; i++) {
	memcpy(packet->data + length, iov[i].iov_base, iov[i].iov_len);
	length += iov[i].iov_len;
    }

    /* Send the packet to the peer over channel id 0. */
    ENETlock(ecd->ecd);
    if (enet_peer_send (ecd->peer, 0, packet) == -1) {
        enet_packet_destroy(packet);
        svc->trace_out(ecd->ecd->cm, "ENET  ======  failed to send a packet to peer %p, state %d\n", ecd->peer, ecd->peer->state);
	return -1;
    }
    ENETunlock(ecd->ecd);

    wake_enet_server_thread(ecd->ecd);

    return (int)iovcnt;
}


static int enet_global_init = 0;

static void
free_enet_data(CManager cm, void *ecdv)
{
    enet_client_data_ptr ecd = (enet_client_data_ptr) ecdv;
    CMtrans_services svc = ecd->svc;
    (void)cm;
    if (ecd->hostname != NULL)
	svc->free_func(ecd->hostname);
    svc->free_func(ecd);
}

static void
shutdown_enet_thread
(CManager cm, void *ecdv)
{
    enet_client_data_ptr ecd = (enet_client_data_ptr) ecdv;
    CMtrans_services svc = ecd->svc;
    (void)cm;
    if (ecd->server != NULL) {
	ENetHost * server = ecd->server;
        ENETlock(ecd);
	enet_host_flush(ecd->server);
        ENETunlock(ecd);
	svc->fd_remove_select(cm, enet_host_get_sock_fd (server));
	svc->remove_periodic(ecd->periodic_handle);
	ecd->server = NULL;
        ENETlock(ecd);
	enet_host_destroy(server);
        ENETunlock(ecd);
    }
}

#ifdef HAVE_WINDOWS_H
static char*
WSAerror_str(err)
int err;
{
    switch(err) {
    case WSAEINTR: return "WSAEINTR";
    case WSAEBADF: return "WSAEBADF";
    case WSAEACCES: return "WSAEACCES";
    case WSAEFAULT: return "WSAEFAULT";
    case WSAEINVAL: return "WSAEINVAL";
    case WSAEMFILE: return "WSAEMFILE";
    case WSAEWOULDBLOCK: return "WSAEWOULDBLOCK";
    case WSAEINPROGRESS: return "WSAEINPROGRESS";
    case WSAEALREADY: return "WSAEALREADY";
    case WSAENOTSOCK: return "WSAENOTSOCK";
    case WSAEDESTADDRREQ: return "WSAEDESTADDRREQ";
    case WSAEMSGSIZE: return "WSAEMSGSIZE";
    case WSAEPROTOTYPE: return "WSAEPROTOTYPE";
    case WSAENOPROTOOPT: return "WSAENOPROTOOPT";
    case WSAEPROTONOSUPPORT: return "WSAEPROTONOSUPPORT";
    case WSAESOCKTNOSUPPORT: return "WSAESOCKTNOSUPPORT";
    case WSAEOPNOTSUPP: return "WSAEOPNOTSUPP";
    case WSAEPFNOSUPPORT: return "WSAEPFNOSUPPORT";
    case WSAEAFNOSUPPORT: return "WSAEAFNOSUPPORT";
    case WSAEADDRINUSE: return "WSAEADDRINUSE";
    case WSAEADDRNOTAVAIL: return "WSAEADDRNOTAVAIL";
    case WSAENETDOWN: return "WSAENETDOWN";
    case WSAENETUNREACH: return "WSAENETUNREACH";
    case WSAENETRESET: return "WSAENETRESET";
    case WSAECONNABORTED: return "WSAECONNABORTED";
    case WSAECONNRESET: return "WSAECONNRESET";
    case WSAENOBUFS: return "WSAENOBUFS";
    case WSAEISCONN: return "WSAEISCONN";
    case WSAENOTCONN: return "WSAENOTCONN";
    case WSAESHUTDOWN: return "WSAESHUTDOWN";
    case WSAETOOMANYREFS: return "WSAETOOMANYREFS";
    case WSAETIMEDOUT: return "WSAETIMEDOUT";
    case WSAECONNREFUSED: return "WSAECONNREFUSED";
    case WSAELOOP: return "WSAELOOP";
    case WSAENAMETOOLONG: return "WSAENAMETOOLONG";
    case WSAEHOSTDOWN: return "WSAEHOSTDOWN";
    case WSAEHOSTUNREACH: return "WSAEHOSTUNREACH";
    case WSAENOTEMPTY: return "WSAENOTEMPTY";
    case WSAEPROCLIM: return "WSAEPROCLIM";
    case WSAEUSERS: return "WSAEUSERS";
    case WSAEDQUOT: return "WSAEDQUOT";
    case WSAESTALE: return "WSAESTALE";
    case WSAEREMOTE: return "WSAEREMOTE";
    case WSAEDISCON: return "WSAEDISCON";
    case WSASYSNOTREADY: return "WSASYSNOTREADY";
    case WSAVERNOTSUPPORTED: return "WSAVERNOTSUPPORTED";
    case WSANOTINITIALISED: return "WSANOTINITIALISED";
    default: return "Unknown Winsock error";
    }
}
/*
 *  Note.  Unfortunately, the _pipe() function on WinNT 
 *  produces FDs that you can't use in select().  This ruins what we want
 *  this pipe for, which is to wake up a thread sleeping in select().
 *  So, we need to introduce a pipe function that returns two socket FDs.
 *  NT Sux.
 */

static int
pipe(filedes)
SOCKET filedes[2];
{
    
    int length;
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    SOCKET sock1, sock2, conn_sock;
    unsigned long block = TRUE;
    int delay_value = 1;
   
    conn_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (conn_sock == SOCKET_ERROR) {
	fprintf(stderr, "Cannot open INET socket\n");
	return -1;
    }
    sock_addr.sin_family = PF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = 0;
    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
	     sizeof sock_addr) == SOCKET_ERROR) {
	fprintf(stderr, "Cannot bind INET socket\n");
	return -1;
    }
    length = sizeof sock_addr;
    if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
	fprintf(stderr, "Cannot get socket name\n");
	return -1;
    }
    /* begin listening for conns */
    if (listen(conn_sock, FD_SETSIZE)) {
	fprintf(stderr, "listen failed\n");
	return -1;
    }

/* send sock */
    if ((sock1 = socket(AF_INET, SOCK_STREAM, 0)) == SOCKET_ERROR) {
	return -1;
    }
    sock_addr.sin_addr.s_addr = 0x0100007f;  /* loopback */
    sock_addr.sin_family = PF_INET;
    if (ioctlsocket(sock1, FIONBIO, &block) != 0) {
	printf("ioctl failed\n");
    }
    if (connect(sock1, (struct sockaddr *) &sock_addr,
		sizeof sock_addr) == SOCKET_ERROR) {
	int err = WSAGetLastError();
	if (err != WSAEWOULDBLOCK) {
	    printf("unexpected error from connect, %s\n", WSAerror_str(err));
	}
    }

    if ((sock2 = accept(conn_sock, (struct sockaddr *) 0, (int *) 0)) == SOCKET_ERROR) {
	    int err = WSAGetLastError();
	    printf("err was %s\n", WSAerror_str(err));
    }
    
    setsockopt(sock2, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));
    {
	fd_set stXcptFDS,stWriteFDS;
	struct timeval stTimeOut;	/* for select() timeout (none) */
	int wRet;

	EVPATH_FD_ZERO((fd_set FAR*)&(stXcptFDS));
	EVPATH_FD_ZERO((fd_set FAR*)&(stWriteFDS));
	FD_SET(sock1, (fd_set FAR*)&(stWriteFDS));
	FD_SET(sock1, (fd_set FAR*)&(stXcptFDS));
	stTimeOut.tv_sec  = 10;
	stTimeOut.tv_usec = 0;
	wRet = select(-1, NULL, 
		      (fd_set FAR*)&(stWriteFDS),
		      (fd_set FAR*)&(stXcptFDS), 
		      NULL);
	if (wRet == SOCKET_ERROR) {
	    int err = WSAGetLastError();
	    printf("err was %s\n", WSAerror_str(err));
	}
    }
    setsockopt(sock1, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));

    filedes[0] = sock1;
    filedes[1] = sock2;
    return 0;
}
#endif

#ifdef __cplusplus
extern "C"
#else
extern
#endif
void *
INTERFACE_NAME(initialize)(CManager cm, CMtrans_services svc,
			 transport_entry trans, attr_list attrs)
{
    static int atom_init = 0;
    int filedes[2];
    char *env = getenv("ENET_HOST_SERVICE_WARN_INTERVAL");

    enet_client_data_ptr enet_data;
    (void)attrs;
    svc->trace_out(cm, "Initialize ENET reliable UDP transport built in %s",
		   EVPATH_MODULE_BUILD_DIR);
    if (enet_global_init == 0) {
	if (enet_initialize () != 0) {
	    fprintf (stderr, "An error occurred while initializing ENet.\n");
	    //return EXIT_FAILURE;
	}
#ifndef USE_ZPL_ENET
        enet_time_set(0);   /* rollover in 50 days, old ENET only */
#endif
    }
    if (atom_init == 0) {
	CM_ENET_HOSTNAME = attr_atom_from_string("CM_ENET_HOST");
	CM_ENET_PORT = attr_atom_from_string("CM_ENET_PORT");
	CM_ENET_ADDR = attr_atom_from_string("CM_ENET_ADDR");
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	CM_PEER_IP = attr_atom_from_string("PEER_IP");
	CM_PEER_LISTEN_PORT = attr_atom_from_string("PEER_LISTEN_PORT");
	CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
	CM_ENET_CONN_TIMEOUT = attr_atom_from_string("CM_ENET_CONN_TIMEOUT");
	CM_ENET_CONN_REUSE = attr_atom_from_string("CM_ENET_CONN_REUSE");
        (void)CM_NETWORK_POSTFIX;
	atom_init++;
    }
    if (env) {
        sscanf(env, "%d", &enet_host_service_warn_interval);
        fprintf(stderr, "DEBUG: Setting enet_host_service_warn_interval to %d\n", enet_host_service_warn_interval);
    }
    enet_data = (enet_client_data_ptr) svc->malloc_func(sizeof(struct enet_client_data));
    memset(enet_data, 0, sizeof(struct enet_client_data));
    thr_mutex_init(enet_data->enet_lock);
    enet_data->enet_locked = 0;
    enet_data->cm = cm;
    enet_data->hostname = NULL;
    enet_data->listen_port = -1;
    enet_data->svc = svc;
    enet_data->server = NULL;
    enet_data->pending_data = NULL;

    if (pipe(filedes) != 0) {
	perror("Pipe for wake not created.  ENET wake mechanism inoperative.");
	return NULL;
    }
    enet_data->wake_read_fd = filedes[0];
    enet_data->wake_write_fd = filedes[1];
    svc->add_shutdown_task(cm, shutdown_enet_thread, (void *) enet_data, SHUTDOWN_TASK);
    svc->add_shutdown_task(cm, free_enet_data, (void *) enet_data, FREE_TASK);
    return (void *) enet_data;
}

#ifdef USE_ZPL_ENET
extern transport_entry cmzplenet_add_static_transport(CManager cm, CMtrans_services svc)
#else
extern transport_entry cmenet_add_static_transport(CManager cm, CMtrans_services svc)
#endif
{
    transport_entry transport;
    transport = (transport_entry) svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
#ifndef USE_ZPL_ENET
    transport->trans_name = strdup("enet");
#else
    transport->trans_name = strdup("zplenet");
#endif
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)INTERFACE_NAME(initialize);
    transport->listen = (CMTransport_listen_func)INTERFACE_NAME(non_blocking_listen);
    transport->initiate_conn = NULL;
    transport->initiate_conn_nonblocking = (CMTransport_NBconn_func)INTERFACE_NAME(initiate_conn_nonblocking);
    transport->finalize_conn_nonblocking = (CMTransport_NBconn_final_func)INTERFACE_NAME(finalize_conn_nonblocking);
    transport->self_check = (CMTransport_self_check_func)INTERFACE_NAME(self_check);
    transport->connection_eq = (CMTransport_connection_eq_func)INTERFACE_NAME(connection_eq);
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)INTERFACE_NAME(shutdown_conn);
    transport->read_block_func = (CMTransport_read_block_func)INTERFACE_NAME(read_block_func);
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)NULL;
    transport->writev_func = (CMTransport_writev_func)INTERFACE_NAME(writev_func);
    transport->get_transport_characteristics = NULL;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}

