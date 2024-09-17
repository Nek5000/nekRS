/***** Includes *****/
#include "config.h"
#include <sys/types.h>

#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <windows.h>
#define getpid()	_getpid()
#else
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

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif

#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 177)
#  pragma warning (disable: 2259)
#  pragma warning (disable: 981)
#endif

struct udp_connection_data;

static atom_t CM_UDP_PORT = -1;
static atom_t CM_UDP_ADDR = -1;
static atom_t CM_IP_HOSTNAME = -1;
static atom_t CM_TRANSPORT = -1;
static atom_t CM_TRANSPORT_RELIABLE = -1;

typedef struct udp_transport_data {
    CManager cm;
    CMtrans_services svc;
    SOCKET socket_fd;
    int self_ip;
    int self_port;
    attr_list characteristics;
    struct udp_connection_data *connections;
} *udp_transport_data_ptr;

#define MSGBUFSIZE 25600

typedef struct udp_connection_data {
    int udp_IP;
    int udp_port;
    struct sockaddr_in dest_addr;
    CMbuffer read_buffer;
    size_t read_buf_len;
    udp_transport_data_ptr utd;
    CMConnection conn;
    attr_list attrs;
    struct udp_connection_data *next;
} *udp_conn_data_ptr;

#ifdef WSAEWOULDBLOCK
#define read(fd, buf, len) recv(fd, buf, len, 0)
#define write(fd, buf, len) send(fd, buf, len, 0)
#endif

static udp_conn_data_ptr
create_udp_conn_data(CMtrans_services svc)
{
    udp_conn_data_ptr udp_conn_data =
	svc->malloc_func(sizeof(struct udp_connection_data));
    udp_conn_data->read_buffer = NULL;
    udp_conn_data->udp_port = -1;
    udp_conn_data->next = NULL;
    return udp_conn_data;
}

static void
add_connection(udp_transport_data_ptr utd, udp_conn_data_ptr ucd)
{
    udp_conn_data_ptr tmp = utd->connections;
    utd->connections = ucd;
    ucd->next = tmp;
}

static void
unlink_connection(udp_transport_data_ptr utd, udp_conn_data_ptr ucd)
{
    if (utd->connections == ucd) {
	utd->connections = ucd->next;
	ucd->next = NULL;
    } else {
	udp_conn_data_ptr tmp = utd->connections;
	while (tmp != NULL) {
	    if (tmp->next == ucd) {
		tmp->next = ucd->next;
		ucd->next = NULL;
		return;
	    }
	}
	printf("Serious internal error, UDP unlink_connection, connection not found\n");
    }
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

extern void
libcmudp_LTX_shutdown_conn(CMtrans_services svc, udp_conn_data_ptr ucd)
{
    unlink_connection(ucd->utd, ucd);
    svc->connection_deref(ucd->conn);
    free_attr_list(ucd->attrs);
    free(ucd);
}

#include "qual_hostname.c"

#ifdef _MSC_VER
static int inet_aton(const char* cp, struct in_addr* addr)
{
    addr->s_addr = inet_addr(cp);
    return (addr->s_addr == INADDR_NONE) ? 0 : 1;
}
#endif

static int
check_host(char *hostname, void *sin_addr)
{
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
}

static int
initiate_udp_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, udp_conn_data_ptr udp_conn_data, attr_list conn_attr_list)
{
    int int_port_num;
    udp_transport_data_ptr utd = (udp_transport_data_ptr) trans->trans_data;
    char *host_name;
    static int udp_ip = 0;
    struct sockaddr_in dest_addr;
    char *network_string;

    memset(&dest_addr, 0, sizeof(dest_addr));
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "UDP transport found no UDP_HOST attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "UDP transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_UDP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &udp_ip)) {
	svc->trace_out(cm, "CMUDP transport found no UDP_ADDR attribute");
	/* wasn't there */
	udp_ip = 0;
    } else {
	svc->trace_out(cm, "CMUDP transport connect to UDP_IP %lx", udp_ip);
    }
    if ((host_name == NULL) && (udp_ip == 0))
	return -1;

    if (!query_attr(attrs, CM_UDP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &int_port_num)) {
	svc->trace_out(cm, "CMUDP transport found no UDP_PORT attribute");
	return -1;
    } else {
	svc->trace_out(cm, "CMUDP transport connect to port %d", int_port_num);
    }

    if (((network_string = getenv("CM_NETWORK")) != NULL) &&
	(host_name != NULL)) {
	size_t name_len = strlen(host_name) + 2 + strlen(network_string);
	char *new_host_name = svc->malloc_func(name_len);
	char *first_dot = strchr(host_name, '.');
	memset(new_host_name, 0, name_len);
	if (first_dot == NULL) {
	    strcpy(new_host_name, host_name);
	    strcat(new_host_name, network_string);
	} else {
	    strncpy(new_host_name, host_name, first_dot - host_name);
	    strcat(new_host_name, network_string);
	    strcat(new_host_name, first_dot);
	}
	if (check_host(new_host_name, (void *) &dest_addr.sin_addr) == 0) {
	    /* host has no NETWORK interface */
	    if (check_host(host_name, (void *) &dest_addr.sin_addr) == 0) {
		svc->trace_out(cm, "--> Host not found \"%s\"",
			       host_name);
	    }
	} else {
	    svc->trace_out(cm, "--> Using non default network interface with hostname %s",
			   new_host_name);
	}
	svc->free_func(new_host_name);
    } else {
	if (host_name != NULL) {
	    if (check_host(host_name, (void *) &dest_addr.sin_addr) == 0) {
		if (udp_ip == 0) {
		    svc->trace_out(cm, "CMSocket connect FAILURE --> Host not found \"%s\", no IP addr supplied in contact list", host_name);
		} else {
		    svc->trace_out(cm, "CMSOCKET --> Host not found \"%s\", Using supplied IP addr %x",
				   host_name == NULL ? "(unknown)" : host_name,
				   udp_ip);
		    dest_addr.sin_addr.s_addr = ntohl(udp_ip);
		}
	    }
	} else {
	    dest_addr.sin_addr.s_addr = ntohl(udp_ip);
	}
    }
    dest_addr.sin_family = AF_INET;
    dest_addr.sin_port = htons(int_port_num);

    svc->trace_out(cm, "--> Connection established");

    udp_conn_data->udp_IP = udp_ip;
    udp_conn_data->udp_port = int_port_num;
    udp_conn_data->dest_addr = dest_addr;
    udp_conn_data->utd = utd;

    return 1;
}

static void
libcmudp_data_available(void *vtrans, void *vinput)
{
    transport_entry trans = vtrans;
    SOCKET input_fd = (SOCKET) (intptr_t) vinput;
    ssize_t nbytes;
    udp_transport_data_ptr utd = (udp_transport_data_ptr) trans->trans_data;
    udp_conn_data_ptr ucd = utd->connections;
    struct sockaddr_in addr;
    unsigned int addrlen = sizeof(addr);
    char *msgbuf;
    int unused;

    if (recvfrom(input_fd, (char*) & unused, 4, MSG_PEEK,
		 (struct sockaddr *) &addr, &addrlen) != 4) {
        return;
    }
    while (ucd != NULL) {
	if (memcmp(&addr, &ucd->dest_addr, sizeof(addr)) == 0) {
	    utd->svc->trace_out(trans->cm, "UDP data available on existing connetion, IP addr %lx\n", 
			   ucd->udp_IP);
	    break;
	}
	ucd = ucd->next;
    }
    if (ucd == NULL) {
	CMConnection conn;
	attr_list conn_attr_list;
	ucd = create_udp_conn_data(utd->svc);

	conn_attr_list = create_attr_list();

	conn = utd->svc->connection_create(trans, ucd, conn_attr_list);
	ucd->dest_addr = addr;
	ucd->udp_IP = ntohl(addr.sin_addr.s_addr);
	ucd->udp_port = ntohs(addr.sin_port);
	ucd->utd = utd;
	ucd->conn = conn;
	ucd->attrs = conn_attr_list;
	ucd->next = NULL;
	add_connection(utd, ucd);
	add_attr(conn_attr_list, CM_UDP_ADDR, Attr_Int4,
		 (attr_value) (long)ucd->udp_IP);
	add_attr(conn_attr_list, CM_UDP_PORT, Attr_Int4,
		 (attr_value) (long)ucd->udp_port);

	utd->svc->trace_out(trans->cm, "UDP data available on new connetion, IP addr %lx\n", 
			   ucd->udp_IP);
    }

    ucd->read_buffer = utd->svc->get_data_buffer(trans->cm, MSGBUFSIZE + 4);

    msgbuf = &((char*)ucd->read_buffer->buffer)[0];
    if ((nbytes = recvfrom(input_fd, msgbuf, MSGBUFSIZE, 0,
			   (struct sockaddr *) &addr, &addrlen)) < 0) {
	perror("recvfrom");
	exit(1);
    }
    ucd->read_buf_len = nbytes;
    /* kick this upstairs */
    trans->data_available(vtrans, ucd->conn);
    utd->svc->return_data_buffer(trans->cm, ucd->read_buffer);
}

/* 
 * Initiate a connection to a udp group.
 */
extern CMConnection
libcmudp_LTX_initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{
    udp_conn_data_ptr udp_conn_data = create_udp_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;

    if (initiate_udp_conn(cm, svc, trans, attrs, udp_conn_data, conn_attr_list) != 1) {
	return NULL;
    }

    add_attr(conn_attr_list, CM_UDP_ADDR, Attr_Int4,
	     (attr_value) (long)udp_conn_data->udp_IP);
    add_attr(conn_attr_list, CM_UDP_PORT, Attr_Int4,
	     (attr_value) (long)udp_conn_data->udp_port);

    conn = svc->connection_create(trans, udp_conn_data, conn_attr_list);
    add_connection(udp_conn_data->utd, udp_conn_data);
    udp_conn_data->conn = conn;
    udp_conn_data->attrs = conn_attr_list;
    svc->connection_addref(conn);  /* one ref count went to select (and CM), 
				the other to the user */
    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 */
extern int
libcmudp_LTX_self_check(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{
    udp_transport_data_ptr utd = trans->trans_data;
    int host_addr;
    int int_port_num;
    char *host_name;
    char my_host_name[256];
    static int IP = 0;

    if (IP == 0) {
	IP = get_self_ip_addr(cm, svc);
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "CMself check UDP transport found no IP_HOST attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_UDP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_addr)) {
	svc->trace_out(cm, "CMself check UDP transport found no UDP_ADDR attribute");
	if (host_name == NULL) return 0;
	host_addr = 0;
    }
    if (!query_attr(attrs, CM_UDP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "CMself check UDP transport found no UDP_PORT attribute");
	return 0;
    }
    get_qual_hostname(cm, my_host_name, sizeof(my_host_name) - 1, svc, NULL, NULL);

    if (host_name && (strcmp(host_name, my_host_name) != 0)) {
	svc->trace_out(cm, "CMself check - Hostnames don't match");
	return 0;
    }
    if (host_addr && (IP != host_addr)) {
	svc->trace_out(cm, "CMself check - Host IP addrs don't match, %lx, %lx", IP, host_addr);
	return 0;
    }
    if (int_port_num != utd->self_port) {
	svc->trace_out(cm, "CMself check - Ports don't match");
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

extern int
libcmudp_LTX_connection_eq(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, udp_conn_data_ptr ucd)
{

    int int_port_num;
    int requested_IP = -1;
    char *host_name = NULL;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "UDP transport found no UDP_HOST attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "UDP transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_UDP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &int_port_num)) {
	svc->trace_out(cm, "Conn Eq CMUdp transport found no UDP_PORT attribute");
	return 0;
    }
    if (!query_attr(attrs, CM_UDP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &requested_IP)) {
	svc->trace_out(cm, "CMUdp transport found no UDP_ADDR attribute");
    }
    svc->trace_out(cm, "CMUdp Conn_eq comparing IP/ports %x/%d and %x/%d",
		   ucd->udp_IP, ucd->udp_port,
		   requested_IP, int_port_num);

    if (requested_IP == -1) {
	check_host(host_name, (void *) &requested_IP);
	svc->trace_out(cm, "IP translation for hostname %s is %x", host_name,
		       requested_IP);
    }

    if ((ucd->udp_IP == requested_IP) &&
	(ucd->udp_port == int_port_num)) {
	svc->trace_out(cm, "CMUdp Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "CMUdp Conn_eq returning FALSE");
    return 0;
}


extern attr_list
libcmudp_LTX_non_blocking_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_info)
{
    udp_transport_data_ptr utd = trans->trans_data;
    int int_port_num = 0;
    u_short port_num;
    attr_list listen_list;
    unsigned int nl;
    int one = 1;
    SOCKET socket_fd;
    struct sockaddr_in addr;
    int IP = get_self_ip_addr(cm, svc);

    if (listen_info != NULL &&
	(!query_attr(listen_info, CM_UDP_PORT, /* type pointer */ NULL,
		     (attr_value *) (intptr_t) &int_port_num))) {
	svc->trace_out(cm, "CMUDP transport found no UDP_PORT attribute");
	int_port_num = 0;
    } else {
	if (int_port_num > USHRT_MAX || int_port_num < 0) {
	    fprintf(stderr, "Requested port number %d is invalid\n", int_port_num);
	    return NULL;
	}
	svc->trace_out(cm, "CMUDP transport connect to port %d", int_port_num);
    }
    if ((socket_fd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	perror("socket");
	exit(1);
    }
    port_num = int_port_num;

    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);	/* N.B.: differs from *
						 * sender */
    addr.sin_port = htons(port_num);
#ifdef SO_REUSEPORT
    if (setsockopt(socket_fd, SOL_SOCKET, SO_REUSEPORT, (char *) &one, sizeof(one)) == -1) {
	perror("setsockopt reuseport");
    }
#else
    if (setsockopt(socket_fd, SOL_SOCKET, SO_REUSEADDR, (char *) &one, sizeof(one)) == -1) {
	perror("setsockopt reuseport");
    }
#endif
    /* bind to receive address */
    if (bind(socket_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
	perror("bind");
	exit(1);
    }
    nl = sizeof(struct sockaddr_in);
    if (getsockname(socket_fd, (struct sockaddr *) &addr, &nl) != 0)
	perror("getsockname");
    

    addr.sin_addr.s_addr = htonl(IP);
    listen_list = create_attr_list();
    add_attr(listen_list, CM_UDP_ADDR, Attr_Int4,
	     (attr_value) (long)IP);
    add_attr(listen_list, CM_UDP_PORT, Attr_Int4,
	     (attr_value) (long) ntohs(addr.sin_port));
    add_attr(listen_list, CM_TRANSPORT, Attr_String,
	     (attr_value) strdup("udp"));
    svc->trace_out(cm, "CMudp Adding libcmudp_data_available as action on fd %d", socket_fd);
    svc->fd_add_select(cm, socket_fd, libcmudp_data_available,
		       (void *) trans, (void *) (intptr_t)socket_fd);
    utd->socket_fd = socket_fd;
    utd->self_ip = IP;
    utd->self_port = ntohs(addr.sin_port);
    return listen_list;
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

/* 
 *  This function will not be used unless there is no read_to_buffer function
 *  in the transport.  It is an example, meant to be copied in transports 
 *  that are more efficient if they allocate their own buffer space.
 */
extern void *
libcmudp_LTX_read_block_func(CMtrans_services svc, udp_conn_data_ptr ucd, size_t *actual_len, size_t *offset_ptr)
{
    *actual_len = ucd->read_buf_len;
    *offset_ptr = 0;
    ucd->read_buf_len = 0;
    return ucd->read_buffer;
}

#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif

extern int
libcmudp_LTX_writev_func(CMtrans_services svc, udp_conn_data_ptr ucd, struct iovec *iov, size_t iovcnt, attr_list attrs)
{
    SOCKET fd = ucd->utd->socket_fd;
    if (ucd->utd->socket_fd == -1) {
	if ((ucd->utd->socket_fd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	    perror("socket");
	    exit(1);
	}
    }
    fd = ucd->utd->socket_fd;
    svc->trace_out(ucd->utd->cm, "CMUdp writev of %d vectors on fd %d",
		   iovcnt, fd);
#ifndef _MSC_VER
    struct sockaddr_in addr = ucd->dest_addr;
    struct msghdr msg;
    memset(&msg, 0, sizeof(msg));
    msg.msg_name = (void*)&addr;
    msg.msg_namelen = sizeof(addr);
    msg.msg_iov = &iov[0];
    msg.msg_iovlen = iovcnt;
    if (sendmsg(fd, &msg, 0) < 0) {
	perror("write sendmsg");
	exit(1);
    }
#else
    // no reimplementation for windows currently
#endif
    return (int)iovcnt;
}

#ifdef HAVE_WINDOWS_H
static int socket_global_init = 0;
/* Winsock init stuff, ask for ver 1.1 */
static WORD wVersionRequested = MAKEWORD(1, 1);
static WSADATA wsaData;
#endif

static void
free_udp_data(CManager cm, void *utdv)
{
    udp_transport_data_ptr utd = (udp_transport_data_ptr) utdv;
    CMtrans_services svc = utd->svc;
    free_attr_list(utd->characteristics);
    svc->free_func(utd);
}

extern void *
libcmudp_LTX_initialize(CManager cm, CMtrans_services svc)
{
    static int atom_init = 0;

    udp_transport_data_ptr udp_data;
    svc->trace_out(cm, "Initialize CMUdp transport");
#ifdef HAVE_WINDOWS_H
    if (socket_global_init == 0) {
	int nErrorStatus;
	/* initialize the winsock package */
	nErrorStatus = WSAStartup(wVersionRequested, &wsaData);
	if (nErrorStatus != 0) {
	    fprintf(stderr, "Could not initialize windows socket library!");
	    WSACleanup();
	    exit(-1);
	}
    }
#endif
    if (atom_init == 0) {
	CM_UDP_PORT = attr_atom_from_string("UDP_PORT");
	CM_UDP_ADDR = attr_atom_from_string("UDP_ADDR");
	CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	CM_TRANSPORT_RELIABLE = attr_atom_from_string("CM_TRANSPORT_RELIABLE");
	atom_init++;
    }
    udp_data = svc->malloc_func(sizeof(struct udp_transport_data));
    udp_data->cm = cm;
    udp_data->svc = svc;
    udp_data->socket_fd = -1;
    udp_data->self_ip = 0;
    udp_data->self_port = 0;
    udp_data->connections = NULL;
    udp_data->characteristics = create_attr_list();
    add_int_attr(udp_data->characteristics, CM_TRANSPORT_RELIABLE, 0);
    svc->add_shutdown_task(cm, free_udp_data, (void *) udp_data, FREE_TASK);
    return (void *) udp_data;
}

extern attr_list
libcmudp_LTX_get_transport_characteristics(transport_entry trans, CMtrans_services svc,
					   udp_transport_data_ptr utd)
{
    add_ref_attr_list(utd->characteristics);
    return utd->characteristics;
}
extern transport_entry
cmudp_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport = svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("udp");
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)libcmudp_LTX_initialize;
    transport->listen = (CMTransport_listen_func)libcmudp_LTX_non_blocking_listen;
    transport->initiate_conn = (CMTransport_conn_func)libcmudp_LTX_initiate_conn;
    transport->self_check = (CMTransport_self_check_func)libcmudp_LTX_self_check;
    transport->connection_eq = (CMTransport_connection_eq_func)libcmudp_LTX_connection_eq;
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)libcmudp_LTX_shutdown_conn;
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)NULL;
    transport->read_block_func = (CMTransport_read_block_func)libcmudp_LTX_read_block_func;;
    transport->writev_func = (CMTransport_writev_func)libcmudp_LTX_writev_func;
    transport->NBwritev_func = (CMTransport_writev_func)NULL;
    transport->set_write_notify = (CMTransport_set_write_notify_func)NULL;
    transport->get_transport_characteristics = (CMTransport_get_transport_characteristics) libcmudp_LTX_get_transport_characteristics;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}
