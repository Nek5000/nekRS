/***** Includes *****/
#include "config.h"
#include <sys/types.h>

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

#include <stdio.h>
#include <fcntl.h>
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
#include <algorithm>
#include <iostream>

#include <atl.h>
#include "evpath.h"
#include "cm_transport.h"
#include "udt.h"

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif

#if defined (__INTEL_COMPILER)
#pragma warning (disable: 869)
#pragma warning (disable: 310)
#pragma warning (disable: 1418)
#pragma warning (disable: 180)
#pragma warning (disable: 2259)
#pragma warning (disable: 177)
#endif

typedef struct func_list_item
{
    select_list_func func;
    void *arg1;
    void *arg2;
} FunctionListElement;

typedef struct socket_data_map_entry *socket_data_map_ptr;

typedef struct udt4_transport_data
{
    CManager cm;
    char *hostname;
    int listen_port;
    attr_list characteristics;
    CMtrans_services svc;
    int eid;
    int listen_sock;
    int conn_count;
    socket_data_map_ptr map;
} *udt4_transport_data_ptr;

typedef struct udt4_connection_data
{
    int remote_IP;
    int remote_contact_port;
    UDTSOCKET udtsock;
    udt4_transport_data_ptr utd;
    CMConnection conn;
    CMbuffer read_buffer;
    int read_buffer_len;
} *udt4_conn_data_ptr;

struct socket_data_map_entry
{
    UDTSOCKET udtsock;
    udt4_conn_data_ptr ucd;
};

#ifdef WSAEWOULDBLOCK
#define EWOULDBLOCK WSAEWOULDBLOCK
#define EAGAIN WSAEINPROGRESS
#define EINTR WSAEINTR
#define errno GetLastError()
#define read(fd, buf, len) recv(fd, buf, len, 0)
#define write(fd, buf, len) send(fd, buf, len, 0)
#endif

static atom_t CM_FD = -1;
static atom_t CM_THIS_CONN_PORT = -1;
static atom_t CM_PEER_CONN_PORT = -1;
static atom_t CM_PEER_IP = -1;
static atom_t CM_PEER_HOSTNAME = -1;
static atom_t CM_PEER_LISTEN_PORT = -1;
static atom_t CM_TRANSPORT = -1;
static atom_t CM_TRANSPORT_RELIABLE = -1;
static atom_t CM_IP_PORT = -1;
static atom_t CM_IP_HOSTNAME = -1;
static atom_t CM_IP_ADDR = -1;

#define TIMING_GUARD_START {     struct timeval t0,t1,diff; gettimeofday(&t0, NULL);
#define TIMING_GUARD_STOP gettimeofday(&t1, NULL);    timersub(&t1, &t0, &diff); if (diff.tv_sec > 0) fprintf(stderr, "TIME GUARD at %s:%d exceeded, time was was <%ld.%06ld> secs\n", __FILE__, __LINE__, (long)diff.tv_sec, (long)diff.tv_usec);}

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
	*((int *) sin_addr) = *((int *) &addr);
    } else {
	memcpy(sin_addr, host_addr->h_addr, host_addr->h_length);
    }
    return 1;
}

static udt4_conn_data_ptr
create_udt4_conn_data(CMtrans_services svc)
{
    udt4_conn_data_ptr ucd =
	(udt4_conn_data_ptr) svc->
	malloc_func(sizeof(struct udt4_connection_data));
    memset(ucd, 0, sizeof(struct udt4_connection_data));
    ucd->remote_contact_port = -1;
    ucd->udtsock = -1;
    return ucd;
}

static void udt4_service_epoll(void *void_trans, void *void_conn_sock);

static void
add_sock_map_data(udt4_transport_data_ptr udt, UDTSOCKET sock,
		  udt4_conn_data_ptr ucd)
{
    if (udt->conn_count == 0) {
	udt->map = (socket_data_map_entry *) malloc(sizeof(udt->map[0]));
    } else {
	udt->map =
	    (socket_data_map_entry *) realloc(udt->map,
					      (udt->conn_count +
					       1) * sizeof(udt->map[0]));
    }
    udt->map[udt->conn_count].udtsock = sock;
    udt->map[udt->conn_count].ucd = ucd;
    udt->conn_count++;
}

extern "C" void
libcmudt4_LTX_shutdown_conn(CMtrans_services svc, udt4_conn_data_ptr ucd)
{
    svc->connection_deref(ucd->conn);
    svc->trace_out(ucd->utd->cm, "UDT4, closing socket %x\n",
		   ucd->udtsock);
    UDT::epoll_remove_usock(ucd->utd->eid, ucd->udtsock);
    UDT::close(ucd->udtsock);
    free(ucd);
}


static int
is_private_192(int IP)
{
    return ((IP & 0xffff0000) == 0xC0A80000);	/* equal 192.168.x.x */
}

static int
is_private_182(int IP)
{
    return ((IP & 0xffff0000) == 0xB6100000);	/* equal 182.16.x.x */
}

static int
is_private_10(int IP)
{
    return ((IP & 0xff000000) == 0x0A000000);	/* equal 10.x.x.x */
}

extern "C" attr_list
libcmudt4_LTX_non_blocking_listen(CManager cm, CMtrans_services svc,
				  transport_entry trans,
				  attr_list listen_info);

static int
initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans,
	      attr_list attrs, udt4_conn_data_ptr ucd,
	      attr_list conn_attr_list)
{
    int sock;

    int int_port_num;
    u_short port_num;
    udt4_transport_data_ptr utd =
	(udt4_transport_data_ptr) trans->trans_data;
    char *host_name;
    int remote_IP = -1;
    static int host_ip = 0;
    union
    {
	struct sockaddr s;
	struct sockaddr_in s_I4;
	struct sockaddr_in6 s_l6;
    } sock_addr;

    if (utd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, utd->cm));
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *) (long) &host_name)) {
	svc->trace_out(cm, "udt4 transport found no IP_HOST attribute");
	host_name = NULL;
    } else {
	svc->trace_out(cm, "udt4 transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *) (long) &host_ip)) {
	svc->trace_out(cm, "udt4 transport found no IP_ADDR attribute");
	/* wasn't there */
	host_ip = 0;
    } else {
	svc->trace_out(cm, "udt4 transport connect to host_IP %lx",
		       host_ip);
    }
    if ((host_name == NULL) && (host_ip == 0))
	return -1;

    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
		    /* value pointer */ 
		    (attr_value *) (long) &int_port_num)) {
	svc->trace_out(cm, "udt4 transport found no IP_PORT attribute");
	return -1;
    } else {
	svc->trace_out(cm, "udt4 transport connect to port %d",
		       int_port_num);
    }
    port_num = int_port_num;

    /* we should already be listening, but if not, listen */
    if (utd->listen_port == -1) {
	attr_list l =
	    libcmudt4_LTX_non_blocking_listen(cm, svc, trans, NULL);
	if (l)
	    free_attr_list(l);
    }
    /* INET socket connection, host_name is the machine name */
    char ip_str[INET_ADDRSTRLEN];

    if ((sock = UDT::socket(AF_INET, SOCK_STREAM, 0)) == SOCKET_ERROR) {
	svc->trace_out(cm,
		       " UDT4 connect FAILURE --> Couldn't create socket");
	return -1;
    }
    ((struct sockaddr_in *) &sock_addr)->sin_family = AF_INET;
    if (host_name != NULL) {
	if (check_host(host_name, (void *) &sock_addr.s_I4.sin_addr) == 0) {
	    if (host_ip == 0) {
		svc->trace_out(cm,
			       "UDT4 connect FAILURE --> Host not found \"%s\", no IP addr supplied in contact list",
			       host_name);
	    } else {
		svc->trace_out(cm,
			       "CMUDT4 --> Host not found \"%s\", Using supplied IP addr %x",
			       host_name == NULL ? "(unknown)" : host_name,
			       host_ip);
		sock_addr.s_I4.sin_addr.s_addr = ntohl(host_ip);
	    }
	}
    } else {
	sock_addr.s_I4.sin_addr.s_addr = ntohl(host_ip);
    }
    sock_addr.s_I4.sin_port = htons(port_num);
    remote_IP = ntohl(sock_addr.s_I4.sin_addr.s_addr);
    if (is_private_192(remote_IP)) {
	svc->trace_out(cm,
		       "Target IP is on a private 192.168.x.x network");
    }
    if (is_private_182(remote_IP)) {
	svc->trace_out(cm, "Target IP is on a private 182.16.x.x network");
    }
    if (is_private_10(remote_IP)) {
	svc->trace_out(cm, "Target IP is on a private 10.x.x.x network");
    }
    inet_ntop(AF_INET, &sock_addr.s_I4.sin_addr, ip_str, INET_ADDRSTRLEN);
    svc->trace_out(cm,
		   "Attempting udt4 socket connection, host=\"%s\", IP = %s, port %d",
		   host_name == 0 ? "(unknown)" : host_name, ip_str,
		   ntohs(sock_addr.s_I4.sin_port));
    if (UDT::
	connect(sock, (struct sockaddr *) &sock_addr,
		sizeof(sock_addr.s_I4)) == SOCKET_ERROR) {
	printf("Errno was %d\n", errno);
	svc->trace_out(cm,
		       "UDT4 connect FAILURE --> Connect() to IP %s failed",
		       ip_str);
	UDT::close(sock);
    }
    int local_listen_port = htonl(utd->listen_port);
    if (UDT::send(sock, (char *) &local_listen_port, 4, 0) < 0) {
	svc->trace_out(cm, "UDT4 send of listen port failed\n");
	return -1;
    }

    svc->trace_out(cm, "--> Connection established");
    ucd->remote_IP = remote_IP;
    ucd->remote_contact_port = int_port_num;
    ucd->udtsock = sock;
    ucd->utd = utd;

    add_attr(conn_attr_list, CM_FD, Attr_Int4, (attr_value) (long) sock);
    add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	     (attr_value) (long) int_port_num);
    add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
	     (attr_value) (long) ucd->remote_IP);
    add_sock_map_data(utd, sock, ucd);
    int localFD = UDT::usock_getfd(sock);
    svc->trace_out(cm, "Adding UDT4 UDP socket %d to select list\n",
		   localFD);
    svc->fd_add_select(cm, localFD, udt4_service_epoll, (void *) trans,
		       (void *) utd);
    return sock;
}

/* 
 * Initiate a socket connection with another data exchange.  If port_num is -1,
 * establish a unix socket connection (name_str stores the file name of
 * the waiting socket).  Otherwise, establish an INET socket connection
 * (name_str stores the machine name).
 */
extern "C" CMConnection
libcmudt4_LTX_initiate_conn(CManager cm, CMtrans_services svc,
			    transport_entry trans, attr_list attrs)
{
    udt4_conn_data_ptr ucd = create_udt4_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;
    int sock;
    udt4_transport_data_ptr utd =
	(udt4_transport_data_ptr) trans->trans_data;

    if (utd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, utd->cm));
    }
    if ((sock =
	 initiate_conn(cm, svc, trans, attrs, ucd, conn_attr_list)) < 0)
	return NULL;

    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (long) ucd->remote_contact_port);
    conn = svc->connection_create(trans, ucd, conn_attr_list);
    ucd->utd = utd;
    ucd->conn = conn;

    svc->trace_out(cm,
		   "Cmudt4 Adding trans->data_available as action on fd %d",
		   sock);
    add_sock_map_data(utd, sock, ucd);
    UDT::epoll_add_usock(utd->eid, sock, NULL);

    free_attr_list(conn_attr_list);
    /* dump_sockinfo("initiate ", sock); */
    svc->connection_addref(conn);	/* one ref count went to select
					 * (and CM), the other to the
					 * user */
    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * For udt4, this involves checking to see if the host name is the 
 * same as ours and if the IP_PORT matches the one we are listening on.
 */
extern "C" int
libcmudt4_LTX_self_check(CManager cm, CMtrans_services svc,
			 transport_entry trans, attr_list attrs)
{

    udt4_transport_data_ptr utd =
	(udt4_transport_data_ptr) trans->trans_data;
    int host_addr;
    int int_port_num;
    char *host_name;
    char my_host_name[256];
    static int IP = 0;

    get_IP_config(my_host_name, sizeof(host_name), &IP, NULL, NULL, NULL,
		  NULL, svc->trace_out, (void *) cm);

    if (IP == 0) {
	if (IP == 0)
	    IP = INADDR_LOOPBACK;
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *) (long) &host_name)) {
	svc->trace_out(cm,
		       "CMself check udt4 transport found no IP_HOST attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *) (long) &host_addr)) {
	svc->trace_out(cm,
		       "CMself check udt4 transport found no IP_ADDR attribute");
	if (host_name == NULL)
	    return 0;
	host_addr = 0;
    }
    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
		    /* value pointer */ 
		    (attr_value *) (long) &int_port_num)) {
	svc->trace_out(cm,
		       "CMself check udt4 transport found no IP_PORT attribute");
	return 0;
    }
    if (host_name && (strcmp(host_name, my_host_name) != 0)) {
	svc->trace_out(cm, "CMself check - Hostnames don't match");
	return 0;
    }
    if (host_addr && (IP != host_addr)) {
	svc->trace_out(cm,
		       "CMself check - Host IP addrs don't match, %lx, %lx",
		       IP, host_addr);
	return 0;
    }
    if (int_port_num != utd->listen_port) {
	svc->trace_out(cm, "CMself check - Ports don't match, %d, %d",
		       int_port_num, utd->listen_port);
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

extern "C" int
libcmudt4_LTX_connection_eq(CManager cm, CMtrans_services svc,
			    transport_entry trans, attr_list attrs,
			    udt4_conn_data_ptr ucd)
{

    int int_port_num;
    int requested_IP = -1;
    char *host_name = NULL;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *) (long) &host_name)) {
	svc->trace_out(cm, "udt4 transport found no IP_HOST attribute");
    }
    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
		    /* value pointer */ 
		    (attr_value *) (long) &int_port_num)) {
	svc->trace_out(cm,
		       "Conn Eq udt4 transport found no IP_PORT attribute");
	return 0;
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
		    /* value pointer */ 
		    (attr_value *) (long) &requested_IP)) {
	svc->trace_out(cm, "udt4 transport found no IP_ADDR attribute");
    }
    if (requested_IP == -1) {
	check_host(host_name, (void *) &requested_IP);
	requested_IP = ntohl(requested_IP);
	svc->trace_out(cm, "IP translation for hostname %s is %x",
		       host_name, requested_IP);
    }

    svc->trace_out(cm, "Socket Conn_eq comparing IP/ports %x/%d and %x/%d",
		   ucd->remote_IP, ucd->remote_contact_port,
		   requested_IP, int_port_num);
    if ((ucd->remote_IP == requested_IP) &&
	(ucd->remote_contact_port == int_port_num)) {
	svc->trace_out(cm, "Socket Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "Socket Conn_eq returning FALSE");
    return 0;
}

static attr_list
build_listen_attrs(CManager cm, CMtrans_services svc,
		   udt4_transport_data_ptr sd, attr_list listen_info,
		   int int_port_num)
{
    char host_name[256];
    attr_list ret_list;
    int IP;
    int use_hostname = 0;

    svc->trace_out(cm, "CMUdt4 listen succeeded on port %d", int_port_num);
    get_IP_config(host_name, sizeof(host_name), &IP, NULL, NULL,
		  &use_hostname, listen_info, svc->trace_out, (void *) cm);

    ret_list = create_attr_list();

    if (sd) {
	sd->hostname = strdup(host_name);
	sd->listen_port = int_port_num;
    }
    if ((IP != 0) && !use_hostname) {
	add_attr(ret_list, CM_IP_ADDR, Attr_Int4, (attr_value) (long) IP);
    }
    if ((getenv("CMUseHostname") != NULL) || use_hostname) {
	add_attr(ret_list, CM_IP_HOSTNAME, Attr_String,
		 (attr_value) strdup(host_name));
    } else if (IP == 0) {
	add_int_attr(ret_list, CM_IP_ADDR, INADDR_LOOPBACK);
    }
    add_attr(ret_list, CM_IP_PORT, Attr_Int4,
	     (attr_value) (long) int_port_num);

    add_attr(ret_list, CM_TRANSPORT, Attr_String,
	     (attr_value) strdup("udt4"));
    return ret_list;
}

/* 
 * Create an IP socket for connection from other CMs
 */
extern "C" attr_list
libcmudt4_LTX_non_blocking_listen(CManager cm, CMtrans_services svc,
				  transport_entry trans,
				  attr_list listen_info)
{
    udt4_transport_data_ptr utd =
	(udt4_transport_data_ptr) trans->trans_data;
    int length;
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    int conn_sock;
    int attr_port_num = 0;
    u_short port_num = 0;
    int port_range_low, port_range_high;
    int use_hostname = 0;
    int int_port_num;
    int IP;
    char host_name[256];

    if (utd->listen_port != -1) {
	/* we're already listening */
	if (port_num == 0) {
	    /* not requesting a specific port, return what we have */
	    return build_listen_attrs(cm, svc, NULL, listen_info,
				      utd->listen_port);
	} else {
	    printf
		("CMlisten_specific() requesting a specific port follows other Udt4 operation which initiated listen at another port.  Only one listen allowed, second listen fails.\n");
	    return NULL;
	}
    }

    conn_sock = UDT::socket(AF_INET, SOCK_STREAM, 0);

    if (conn_sock == SOCKET_ERROR) {
	fprintf(stderr, "Cannot open INET socket\n");
	return NULL;
    }
    if (utd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, utd->cm));
    }
    /* 
     *  Check to see if a bind to a specific port was requested
     */
    if (listen_info != NULL
	&& !query_attr(listen_info, CM_IP_PORT,
		       NULL, (attr_value *) (long) &attr_port_num)) {
	port_num = 0;
    } else {
	if (attr_port_num > USHRT_MAX || attr_port_num < 0) {
	    fprintf(stderr, "Requested port number %d is invalid\n",
		    attr_port_num);
	    return NULL;
	}
	port_num = attr_port_num;
    }

    svc->trace_out(cm, "UDT4 begin listen, requested port %d",
		   attr_port_num);
    get_IP_config(host_name, sizeof(host_name), &IP, &port_range_low,
		  &port_range_high, &use_hostname, listen_info,
		  svc->trace_out, (void *) cm);

    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = htons(port_num);

    if (sock_addr.sin_port != 0) {
	/* specific port requested. set REUSEADDR, REUSEPORT because
	 * previous server might have died badly */
	if (UDT::
	    setsockopt(conn_sock, SOL_SOCKET, UDT_REUSEADDR,
		       (char *) &sock_opt_val,
		       sizeof(sock_opt_val)) != 0) {
	    fprintf(stderr,
		    "Failed to set REUSEADDR on INET socket before bind\n");
	    perror("setsockopt(SO_REUSEADDR) failed");
	    return NULL;
	}
	svc->trace_out(cm, "UDT4 trying to bind selected port %d",
		       port_num);
	if (UDT::
	    bind(conn_sock, (struct sockaddr *) &sock_addr,
		 sizeof sock_addr) == SOCKET_ERROR) {
	    fprintf(stderr, "Cannot bind INET socket\n");
	    return NULL;
	}
	/* begin listening for conns */
	if (UDT::listen(conn_sock, FD_SETSIZE)) {
	    fprintf(stderr, "listen failed %s\n",
		    UDT::getlasterror().getErrorMessage());
	    return NULL;
	}
    } else if (port_range_high == -1) {
	svc->trace_out(cm, "UDT4 trying to bind to any available port");
	sock_addr.sin_port = 0;
	if (UDT::
	    bind(conn_sock, (struct sockaddr *) &sock_addr,
		 sizeof sock_addr) == SOCKET_ERROR) {
	    fprintf(stderr, "Cannot bind INET socket\n");
	    return NULL;
	}
	/* begin listening for conns */
	if (UDT::listen(conn_sock, FD_SETSIZE)) {
	    fprintf(stderr, "listen failed %s\n",
		    UDT::getlasterror().getErrorMessage());
	    return NULL;
	}

    } else {
	long seedval = time(NULL) + getpid();
	/* port num is free.  Constrain to range to standards */
	int size = port_range_high - port_range_low;
	int tries = 30;
	int result = SOCKET_ERROR;
	srand48(seedval);
	while (tries > 0) {
	    int target = port_range_low + size * drand48();
	    sock_addr.sin_port = htons(target);
	    int_port_num = target;
	    svc->trace_out(cm, "UDT4 trying to bind port %d", target);
	    result = UDT::bind(conn_sock, (struct sockaddr *) &sock_addr,
			       sizeof sock_addr);
	    tries--;
	    /* begin listening for conns */
	    result = UDT::listen(conn_sock, FD_SETSIZE);
	    if (result != UDT::ERROR) {
		tries = 0;
	    } else {
		/* 
		 *  For some reason, UDT4 gives conflicting port errors
		 *  on listen rather than bind.  You can't bind a socket
		 *  more than once, so on error we have to close the old
		 *  and start with a new UDT::socket.
		 */
		UDT::close(conn_sock);
		conn_sock = UDT::socket(AF_INET, SOCK_STREAM, 0);
	    }

	    if (tries % 5 == 4) {
		/* try reseeding in case we're in sync with another
		 * process */
		srand48(time(NULL) + getpid());
	    }
	    if (tries == 20) {
		/* damn, tried a lot, increase the range (This might
		 * violate specified range) */
		size *= 10;
	    }
	    if (tries == 10) {
		/* damn, tried a lot more, increase the range (This might
		 * violate specified range) */
		size *= 10;
	    }
	}
	if (result == SOCKET_ERROR) {
	    fprintf(stderr, "Cannot bind INET socket\n");
	    return NULL;
	}
    }
    length = sizeof sock_addr;
    if (UDT::
	getsockname(conn_sock, (struct sockaddr *) &sock_addr,
		    &length) < 0) {
	fprintf(stderr, "Cannot get socket name\n");
	return NULL;
    }
    /* set the port num as one we can be contacted at */

    UDT::epoll_add_usock(utd->eid, conn_sock, NULL);
    int localFD = UDT::usock_getfd(conn_sock);
    svc->trace_out(cm, "Adding UDT4 UDP socket %d to select list\n",
		   localFD);
    svc->fd_add_select(cm, localFD, udt4_service_epoll, (void *) trans,
		       (void *) utd);
    utd->listen_sock = conn_sock;

    {
	attr_list ret_list;

	svc->trace_out(cm, "UDT4 listen succeeded on port %d, fd %d",
		       int_port_num, conn_sock);
	ret_list = create_attr_list();

	utd->hostname = strdup(host_name);
	utd->listen_port = int_port_num;
	add_attr(ret_list, CM_TRANSPORT, Attr_String,
		 (attr_value) strdup("udt4"));
	if ((IP != 0) && (!use_hostname)) {
	    add_attr(ret_list, CM_IP_ADDR, Attr_Int4,
		     (attr_value) (long) IP);
	}
	if ((getenv("Cmudt4UseHostname") != NULL) || use_hostname) {
	    add_attr(ret_list, CM_IP_HOSTNAME, Attr_String,
		     (attr_value) strdup(host_name));
	} else if (IP == 0) {
	    add_attr(ret_list, CM_IP_ADDR, Attr_Int4,
		     (attr_value) INADDR_LOOPBACK);
	}
	add_attr(ret_list, CM_IP_PORT, Attr_Int4,
		 (attr_value) (long) int_port_num);

	return ret_list;
    }
}

#ifdef NEED_IOVEC_DEFINE
struct iovec
{
    void *iov_base;
    long iov_len;
};

#endif


#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif


extern "C" int
libcmudt4_LTX_writev_func(CMtrans_services svc, udt4_conn_data_ptr ucd,
			  void *iovs, int iovcnt, attr_list attrs)
{
    int left = 0;
    struct iovec *iov = (struct iovec *) iovs;
    /* sum lengths */
    for (int i = 0; i < iovcnt; i++)
	left += iov[i].iov_len;

    svc->trace_out(ucd->utd->cm, "UDT4 writev of %d bytes on socket %d",
		   left, ucd->udtsock);
    /* 
     * stupid implementation of writev first
     *
     * we'll copy all the data into a single block, then send it with one send()
     */
    char *buffer = (char *) malloc(left), *tmp;
    tmp = buffer;
    for (int i = 0; i < iovcnt; i++) {
	memcpy(tmp, iov[i].iov_base, iov[i].iov_len);
	tmp += iov[i].iov_len;
    }
    tmp = buffer;
    while (left > 0) {
	int sent = UDT::send(ucd->udtsock, tmp, left, 0);

	if (sent < 0) {
	    svc->trace_out(ucd->utd->cm,
			   "UDT4 send failed, error was %s\n",
			   UDT::getlasterror().getErrorMessage());
	    return 0;
	}
	left -= sent;
	tmp += sent;
    }
    free(buffer);
    return iovcnt;
}

static void
free_udt4_data(CManager cm, void *utdv)
{
    udt4_transport_data_ptr utd = (udt4_transport_data_ptr) utdv;
    CMtrans_services svc = utd->svc;
    if (utd->hostname != NULL)
	svc->free_func(utd->hostname);
    svc->free_func(utd);
}

extern "C" int
libcmudt4_LTX_read_to_buffer_func(CMtrans_services svc,
				  udt4_conn_data_ptr ucd, void *buffer,
				  int requested_len, int non_blocking)
{
    int left, iget;
    int rcv_size, var_size = sizeof(rcv_size);
    svc->trace_out(ucd->utd->cm,
		   "CMUDT4 read of %d bytes on fd %d, non_block %d",
		   requested_len, ucd->udtsock, non_blocking);
    UDT::getsockopt(ucd->udtsock, 0, UDT_RCVDATA, &rcv_size, &var_size);
    if (non_blocking && (rcv_size == 0)) {
	svc->trace_out(ucd->utd->cm,
		       "CMUDT4 socket state is %d, rcv_size = %d, non_blocking read would block, returning 0\n",
		       UDT::getsockstate(ucd->udtsock), rcv_size);
	return 0;
    }

    if (UDT::ERROR ==
	(iget =
	 UDT::recv(ucd->udtsock, (char *) buffer, requested_len, 0))) {
	svc->trace_out(ucd->utd->cm, "CMUDT4 recv failed, error was %s\n",
		       UDT::getlasterror().getErrorMessage());;
    }
    if ((iget == -1) || (iget == 0)) {
	int lerrno = errno;
	if ((lerrno != EWOULDBLOCK) &&
	    (lerrno != EAGAIN) && (lerrno != EINTR)) {
	    /* serious error */
	    svc->trace_out(ucd->utd->cm,
			   "CMUDT4 iget was -1, errno is %d, returning 0 for read",
			   lerrno);
	    return -1;
	} else {
	    if (non_blocking) {
		svc->trace_out(ucd->utd->cm,
			       "CMUDT4 iget was -1, would block, errno is %d",
			       lerrno);
		return 0;
	    }
	    return -1;
	}
    }
    left = requested_len - iget;
    while (left > 0) {
	int lerrno;
	if (UDT::ERROR ==
	    (iget =
	     UDT::recv(ucd->udtsock,
		       ((char *) buffer) + requested_len - left, left,
		       0))) {
	    svc->trace_out(ucd->utd->cm, "recv failed, message was %s\n",
			   UDT::getlasterror().getErrorMessage());
	}
	lerrno = errno;
	if (iget == -1) {
	    if ((lerrno != EWOULDBLOCK) &&
		(lerrno != EAGAIN) && (lerrno != EINTR)) {
		/* serious error */
		svc->trace_out(ucd->utd->cm,
			       "Cmudt4 iget was -1, errno is %d, returning %d for read",
			       lerrno, requested_len - left);
		return (requested_len - left);
	    } else {
		iget = 0;
	    }
	} else if (iget == 0) {
	    svc->trace_out(ucd->utd->cm,
			   "Cmudt4 iget was 0, errno is %d, returning %d for read",
			   lerrno, requested_len - left);
	    return requested_len - left;	/* end of file */
	}
	left -= iget;
    }
    return requested_len;
}

static void
udt4_service_epoll(void *void_trans, void *not_used)
{
    transport_entry trans = (transport_entry) void_trans;
    udt4_transport_data_ptr utd =
	(udt4_transport_data_ptr) trans->trans_data;
    CManager cm = utd->cm;
    CMtrans_services svc = utd->svc;
    std::set < UDTSOCKET > readfds;
    if (!(CM_LOCKED(svc, utd->cm))) {
	printf("UTD4 service network, CManager not locked\n");
    }

    int ret = UDT::epoll_wait(utd->eid, &readfds, NULL, 50);
    svc->trace_out(utd->cm, "Epoll_wait returned %d\n", ret);
    for (std::set < UDTSOCKET >::iterator i = readfds.begin();
	 i != readfds.end(); ++i) {
	if (*i == utd->listen_sock) {
	    sockaddr_storage clientaddr;
	    attr_list conn_attr_list;
	    struct sockaddr_in remote;
	    int remote_len = sizeof(remote);
	    int addrlen = sizeof(clientaddr);
	    UDTSOCKET new_sock =
		UDT::accept(*i, (sockaddr *) & clientaddr, &addrlen);
	    udt4_conn_data_ptr ucd = create_udt4_conn_data(svc);
	    ucd->udtsock = new_sock;
	    ucd->utd = utd;
	    conn_attr_list = create_attr_list();
	    CMConnection conn =
		svc->connection_create(trans, ucd, conn_attr_list);
	    ucd->conn = conn;
	    UDT::getpeername(ucd->udtsock, (struct sockaddr *) &remote,
			     &remote_len);
	    add_int_attr(conn_attr_list, CM_PEER_IP,
			 ntohl(remote.sin_addr.s_addr));
	    ucd->remote_IP = ntohl(remote.sin_addr.s_addr);	/* remote_IP 
								 * is in
								 * host
								 * byte
								 * order */
	    if (UDT::
		recv(ucd->udtsock, (char *) &ucd->remote_contact_port, 4,
		     0) < 0) {
		svc->trace_out(utd->cm,
			       "Remote host dropped connection without data");
		return;
	    } else {
		ucd->remote_contact_port = ntohl(ucd->remote_contact_port);
	    }
	    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
		     (attr_value) (long) ucd->remote_contact_port);
	    struct in_addr addr;
	    addr.s_addr = htonl(ucd->remote_IP);
	    svc->trace_out(trans->cm,
			   "Remote host (IP %s) is listening at port %d\n",
			   inet_ntoa(addr), ucd->remote_contact_port);
	    free_attr_list(conn_attr_list);
	    add_sock_map_data(utd, new_sock, ucd);
	    UDT::epoll_add_usock(utd->eid, new_sock);
	} else {
	    // int rcv_size = 1024000;
	    // int var_size = sizeof(int), rs;
	    udt4_conn_data_ptr ucd = NULL;
	    for (int j = 0; j < utd->conn_count; j++) {
		if (utd->map[j].udtsock == *i) {
		    ucd = utd->map[j].ucd;
		}
	    }
	    if (ucd == NULL) {
		printf
		    ("Internal consistency error, conn data for sock %x not found\n",
		     *i);
	    }
	    svc->trace_out(cm, "Calling data available on connection %p\n",
			   ucd->conn);
	    (trans->data_available) (trans, ucd->conn);
	}
    }
}

extern "C" void *
libcmudt4_LTX_initialize(CManager cm, CMtrans_services svc,
			 transport_entry trans)
{
    static int atom_init = 0;

    (void) trans;
    udt4_transport_data_ptr udt4_data;
    svc->trace_out(cm, "Initialize UDP4 transport built in %s",
		   EVPATH_MODULE_BUILD_DIR);

    if (atom_init == 0) {
	CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
	CM_IP_PORT = attr_atom_from_string("IP_PORT");
	CM_IP_ADDR = attr_atom_from_string("IP_ADDR");
	CM_FD = attr_atom_from_string("CONNECTION_FILE_DESCRIPTOR");
	CM_THIS_CONN_PORT = attr_atom_from_string("THIS_CONN_PORT");
	CM_PEER_CONN_PORT = attr_atom_from_string("PEER_CONN_PORT");
	CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
	CM_PEER_IP = attr_atom_from_string("PEER_IP");
	CM_PEER_HOSTNAME = attr_atom_from_string("PEER_HOSTNAME");
	CM_PEER_LISTEN_PORT = attr_atom_from_string("PEER_LISTEN_PORT");
	CM_TRANSPORT_RELIABLE =
	    attr_atom_from_string("CM_TRANSPORT_RELIABLE");
	atom_init++;
    }
    UDT::startup();
    udt4_data =
	(udt4_transport_data_ptr) svc->
	malloc_func(sizeof(struct udt4_transport_data));
    memset(udt4_data, 0, sizeof(struct udt4_transport_data));
    udt4_data->cm = cm;
    udt4_data->hostname = NULL;
    udt4_data->listen_port = -1;
    udt4_data->svc = svc;
    udt4_data->characteristics = create_attr_list();
    add_int_attr(udt4_data->characteristics, CM_TRANSPORT_RELIABLE, 1);
    svc->add_shutdown_task(cm, free_udt4_data, (void *) udt4_data,
			   FREE_TASK);
    udt4_data->eid = UDT::epoll_create();
    int localFD = UDT::epoll_getfd(udt4_data->eid);
    svc->trace_out(cm, "Adding FD %d to select list", localFD);

    svc->fd_add_select(cm, localFD, udt4_service_epoll,
		       (void *) trans, (void *) udt4_data);
    // svc->add_periodic_task(cm, 1, 0, (CMPollFunc) udt4_poll_epoll,
    // (void *)trans);
    return (void *) udt4_data;
}

extern "C" attr_list
libcmudt4_LTX_get_transport_characteristics(transport_entry trans,
					    CMtrans_services svc,
					    void *vutd)
{
    struct udt4_transport_data *utd = (struct udt4_transport_data *) vutd;
    return utd->characteristics;
}

extern "C" void *
libcmudt4_LTX_read_block_func(CMtrans_services svc,
			      udt4_conn_data_ptr conn_data,
			      int *actual_len, int *offset_ptr)
{
    CMbuffer cb;

    if (conn_data->read_buffer_len == -1)
	return NULL;

    *actual_len = conn_data->read_buffer_len;
    *offset_ptr = 0;
    cb = conn_data->read_buffer;
    conn_data->read_buffer_len = 0;
    conn_data->read_buffer = NULL;
    return cb;
}

extern "C" transport_entry
cmudt4_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport =
	(transport_entry) svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("udt4");
    transport->cm = cm;
    transport->transport_init =
	(CMTransport_func) libcmudt4_LTX_initialize;
    transport->listen =
	(CMTransport_listen_func) libcmudt4_LTX_non_blocking_listen;
    transport->initiate_conn =
	(CMTransport_conn_func) libcmudt4_LTX_initiate_conn;
    transport->self_check =
	(CMTransport_self_check_func) libcmudt4_LTX_self_check;
    transport->connection_eq =
	(CMTransport_connection_eq_func) libcmudt4_LTX_connection_eq;
    transport->shutdown_conn =
	(CMTransport_shutdown_conn_func) libcmudt4_LTX_shutdown_conn;
    transport->read_block_func = (CMTransport_read_block_func) NULL;
    transport->read_to_buffer_func =
	(CMTransport_read_to_buffer_func)
	libcmudt4_LTX_read_to_buffer_func;
    transport->writev_func =
	(CMTransport_writev_func) libcmudt4_LTX_writev_func;
    transport->NBwritev_func = (CMTransport_writev_func) NULL;

    transport->set_write_notify = (CMTransport_set_write_notify_func) NULL;
    transport->get_transport_characteristics =
	(CMTransport_get_transport_characteristics)
	libcmudt4_LTX_get_transport_characteristics;
    if (transport->transport_init) {
	transport->trans_data =
	    transport->transport_init(cm, svc, transport);
    }
    return transport;
}
