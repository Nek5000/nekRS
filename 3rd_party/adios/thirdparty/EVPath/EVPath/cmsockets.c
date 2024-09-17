/***** Includes *****/
#include "config.h"
#include <sys/types.h>

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

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif

#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 2259)
#  pragma warning (disable: 177)
#endif

typedef struct func_list_item {
    select_list_func func;
    void *arg1;
    void *arg2;
} FunctionListElement;

typedef struct socket_client_data {
    CManager cm;
    char *hostname;
    int listen_count;
    SOCKET *listen_fds;
    int *listen_ports;
    attr_list characteristics;
    CMtrans_services svc;
} *socket_client_data_ptr;

typedef enum {Block, Non_Block} socket_block_state;

typedef struct socket_connection_data {
    int remote_IP;
    int remote_contact_port;
    SOCKET fd;
    socket_client_data_ptr sd;
    socket_block_state block_state;
    CMConnection conn;
} *socket_conn_data_ptr;

#ifdef _MSC_VER
#define read(fd, buf, len) recv(fd, buf, len, 0)
#define write(fd, buf, len) send(fd, buf, len, 0)
#define close(x) closesocket(x)
#define INST_ADDRSTRLEN 50
#endif

static atom_t CM_FD = -1;
static atom_t CM_THIS_CONN_PORT = -1;
static atom_t CM_PEER_CONN_PORT = -1;
static atom_t CM_PEER_IP = -1;
static atom_t CM_PEER_HOSTNAME = -1;
static atom_t CM_PEER_LISTEN_PORT = -1;
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
	if (inet_pton(PF_INET, hostname, &addr) == 0) {
//	if (inet_aton(hostname, &addr) == 0) {
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

static socket_conn_data_ptr 
create_socket_conn_data(CMtrans_services svc)
{
    socket_conn_data_ptr socket_conn_data =
    svc->malloc_func(sizeof(struct socket_connection_data));
    memset(socket_conn_data, 0, sizeof(struct socket_connection_data));
    socket_conn_data->remote_contact_port = -1;
    socket_conn_data->fd = 0;
    socket_conn_data->block_state = Block;
    return socket_conn_data;
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

#ifndef INET_ADDRSTRLEN
#define INET_ADDRSTRLEN 50
#endif

/* 
 * Accept socket connection
 */
static void
socket_accept_conn(void *void_trans, void *void_conn_sock)
{
    transport_entry trans = (transport_entry) void_trans;
    int conn_sock = (int) (intptr_t) void_conn_sock;
    socket_client_data_ptr sd = (socket_client_data_ptr) trans->trans_data;
    CMtrans_services svc = sd->svc;
    socket_conn_data_ptr socket_conn_data;
    SOCKET sock;
    struct sockaddr sock_addr;
    unsigned int sock_len = sizeof(sock_addr);
    int int_port_num;
    struct linger linger_val;
    int sock_opt_val = 1;

#ifdef TCP_NODELAY
    int delay_value = 1;
#endif
    CMConnection conn;
    attr_list conn_attr_list = NULL;

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    svc->trace_out(sd->cm, "Trying to accept something, socket %d\n", conn_sock);
    linger_val.l_onoff = 1;
    linger_val.l_linger = 60;
    if ((sock = accept(conn_sock, (struct sockaddr *) 0, (unsigned int *) 0)) == SOCKET_ERROR) {
	perror("Cannot accept socket connection");
	svc->fd_remove_select(sd->cm, conn_sock);
	fprintf(stderr, "failure in CMsockets  removing socket connection\n");
	return;
    }
    sock_opt_val = 1;
    setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, (char *) &sock_opt_val,
	       sizeof(sock_opt_val));
    if (setsockopt(sock, SOL_SOCKET, SO_LINGER, (char *) &linger_val,
		   sizeof(struct linger)) != 0) {
	perror("set SO_LINGER");
	return;
    }
#ifdef TCP_NODELAY
    setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));
#endif
    socket_conn_data = create_socket_conn_data(svc);
    socket_conn_data->sd = sd;
    socket_conn_data->fd = sock;
    conn_attr_list = create_attr_list();
    conn = svc->connection_create(trans, socket_conn_data, conn_attr_list);
    socket_conn_data->conn = conn;

    add_attr(conn_attr_list, CM_FD, Attr_Int4,
	     (attr_value) (intptr_t)sock);

    sock_len = sizeof(sock_addr);
    memset(&sock_addr, 0, sock_len);
    getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
    int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
    add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)int_port_num);

    memset(&sock_addr, 0, sizeof(sock_addr));
    sock_len = sizeof(sock_addr);
    if (getpeername(sock, &sock_addr, &sock_len) == 0) {
	int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
	add_attr(conn_attr_list, CM_PEER_CONN_PORT, Attr_Int4,
		 (attr_value) (intptr_t)int_port_num);
	socket_conn_data->remote_IP = ntohl(((struct sockaddr_in *) &sock_addr)->sin_addr.s_addr);
	add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
		 (attr_value) (intptr_t)socket_conn_data->remote_IP);
    }
    {
        char str[INET_ADDRSTRLEN];

        inet_ntop(AF_INET, &((struct sockaddr_in *) &sock_addr)->sin_addr, str, INET_ADDRSTRLEN);
	svc->trace_out(sd->cm, "Accepted TCP/IP socket connection from host at IP %s", 
		       str);
    }
    if (read(sock, (char *) &socket_conn_data->remote_contact_port, 4) != 4) {
	svc->trace_out(sd->cm, "Remote host dropped connection without data");
	return;
    }
    socket_conn_data->remote_contact_port =
	ntohs(socket_conn_data->remote_contact_port);
    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)socket_conn_data->remote_contact_port);
    svc->trace_out(sd->cm, "Remote host (IP %x) is listening at port %d\n",
		   socket_conn_data->remote_IP,
		   socket_conn_data->remote_contact_port);

/* dump_sockinfo("accept ", sock); */
    if (trans->data_available) {
        svc->fd_add_select(sd->cm, sock,
                           (void (*)(void *, void *)) trans->data_available,
                           (void *) trans, (void *) conn);
    }
    free_attr_list(conn_attr_list);
}

extern void
libcmsockets_LTX_shutdown_conn(CMtrans_services svc, socket_conn_data_ptr scd)
{
    svc->connection_deref(scd->conn);
    svc->fd_remove_select(scd->sd->cm, scd->fd);
    close(scd->fd);
    free(scd);
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

static SOCKET
initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, socket_conn_data_ptr socket_conn_data, attr_list conn_attr_list)
{
    SOCKET sock;

#ifdef TCP_NODELAY
    int delay_value = 1;
#endif
    struct linger linger_val;
    int sock_opt_val = 1;
    int int_port_num;
    u_short port_num;
    socket_client_data_ptr sd = (socket_client_data_ptr) trans->trans_data;
    char *host_name;
    int remote_IP = -1;
    static int host_ip = 0;
    unsigned int sock_len;
    union {
	struct sockaddr s;
	struct sockaddr_in s_I4;
	struct sockaddr_in6 s_l6;
    } sock_addr;

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "TCP/IP transport found no IP_HOST attribute");
	host_name = NULL;
    } else {
        svc->trace_out(cm, "TCP/IP transport connect to host %s", host_name);
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_ip)) {
	svc->trace_out(cm, "TCP/IP transport found no IP_ADDR attribute");
	/* wasn't there */
	host_ip = 0;
    } else {
        svc->trace_out(cm, "TCP/IP transport connect to host_IP %lx", host_ip);
    }
    if ((host_name == NULL) && (host_ip == 0))
	return -1;

    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "TCP/IP transport found no IP_PORT attribute");
	return -1;
    } else {
        svc->trace_out(cm, "TCP/IP transport connect to port %d", int_port_num);
    }
    port_num = int_port_num;
    linger_val.l_onoff = 1;
    linger_val.l_linger = 60;

    if (int_port_num == -1) {
#if defined(AF_UNIX) && !defined(HAVE_WINDOWS_H)
	/* unix socket connection, host_name is the file name */
	struct sockaddr_un sock_addru;
	if ((sock = socket(AF_UNIX, SOCK_STREAM, 0)) < 0) {
	    return -1;
	}
	sock_addru.sun_family = AF_UNIX;
	strcpy(sock_addru.sun_path, host_name);
	if (connect(sock, (struct sockaddr *) &sock_addru,
		    sizeof sock_addru) < 0) {
	    return -1;
	}
#else
	fprintf(stderr, "socket initiate_conn port_num parameter == -1 and unix sockets not available.\n");
	return -1;
#endif
    } else {
	/* INET socket connection, host_name is the machine name */
        char ip_str[INET_ADDRSTRLEN];

	if ((sock = socket(AF_INET, SOCK_STREAM, 0)) == SOCKET_ERROR) {
	    svc->trace_out(cm, " CMSocket connect FAILURE --> Couldn't create socket");
	    return -1;
	}
	((struct sockaddr_in *) &sock_addr)->sin_family = AF_INET;
	if (host_name != NULL) {
	    if (check_host(host_name, (void *) &sock_addr.s_I4.sin_addr) == 0) {
		if (host_ip == 0) {
		    svc->trace_out(cm, "CMSocket connect FAILURE --> Host not found \"%s\", no IP addr supplied in contact list", host_name);
		} else {
		    svc->trace_out(cm, "CMSOCKET --> Host not found \"%s\", Using supplied IP addr %x",
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
	    svc->trace_out(cm, "Target IP is on a private 192.168.x.x network");
	}
	if (is_private_182(remote_IP)) {
	    svc->trace_out(cm, "Target IP is on a private 182.16.x.x network");
	}
	if (is_private_10(remote_IP)) {
	    svc->trace_out(cm, "Target IP is on a private 10.x.x.x network");
	}
        inet_ntop(AF_INET, &sock_addr.s_I4.sin_addr, ip_str, INET_ADDRSTRLEN);
	svc->trace_out(cm, "Attempting TCP/IP socket connection, host=\"%s\", IP = %s, port %d",
		       host_name == 0 ? "(unknown)" : host_name, ip_str,
		       ntohs(sock_addr.s_I4.sin_port));
	if (connect(sock, (struct sockaddr *) &sock_addr,
		    sizeof(sock_addr.s_I4)) == SOCKET_ERROR) {
#ifdef WSAEWOULDBLOCK
	    int err = WSAGetLastError();
	    if (err != WSAEWOULDBLOCK || err != WSAEINPROGRESS) {
#endif
		svc->trace_out(cm, "CMSocket connect FAILURE --> Connect() to IP %s failed", ip_str);
		close(sock);
#ifdef WSAEWOULDBLOCK
	    }
#endif
	}
    }

    sock_opt_val = 1;
    setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, (char *) &sock_opt_val,
	       sizeof(sock_opt_val));
    setsockopt(sock, SOL_SOCKET, SO_LINGER, (char *) &linger_val,
	       sizeof(struct linger));

#ifdef TCP_NODELAY
    setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));
#endif

    {
	int local_listen_port = 0;
	if (sd->listen_count) {
	    local_listen_port = htons(sd->listen_ports[0]);
	}
	if (write(sock, (const char *) & local_listen_port, 4) != 4) {
	    svc->trace_out(cm, "Write failed\n");
	    return -1;
	}
    }
    svc->trace_out(cm, "--> Connection established");
    socket_conn_data->remote_IP = remote_IP;
    socket_conn_data->remote_contact_port = int_port_num;
    socket_conn_data->fd = sock;
    socket_conn_data->sd = sd;

    add_attr(conn_attr_list, CM_FD, Attr_Int4,
	     (attr_value) (intptr_t)sock);
    sock_len = sizeof(sock_addr);
    getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
    int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
    add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)int_port_num);
    add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
	     (attr_value) (intptr_t)socket_conn_data->remote_IP);
    return sock;
}

/* 
 * Initiate a socket connection with another data exchange.  If port_num is -1,
 * establish a unix socket connection (name_str stores the file name of
 * the waiting socket).  Otherwise, establish an INET socket connection
 * (name_str stores the machine name).
 */
extern CMConnection
libcmsockets_LTX_initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{
    socket_conn_data_ptr socket_conn_data = create_socket_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;
    SOCKET sock;
    socket_client_data_ptr sd = trans->trans_data;

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    if ((sock = initiate_conn(cm, svc, trans, attrs, socket_conn_data, conn_attr_list)) < 0)
	return NULL;

    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (intptr_t)socket_conn_data->remote_contact_port);
    conn = svc->connection_create(trans, socket_conn_data, conn_attr_list);
    socket_conn_data->conn = conn;

    svc->trace_out(cm, "CMSockets Adding trans->data_available as action on fd %d", sock);
    if (trans->data_available) {
        svc->fd_add_select(cm, sock, (select_list_func) trans->data_available,
                           (void *) trans, (void *) conn);
    }

    free_attr_list(conn_attr_list);
/* dump_sockinfo("initiate ", sock); */
    svc->connection_addref(conn);  /* one ref count went to select (and CM), 
				the other to the user */
    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * For sockets, this involves checking to see if the host name is the 
 * same as ours and if the IP_PORT matches the one we are listening on.
 */
extern int
libcmsockets_LTX_self_check(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{

    socket_client_data_ptr sd = trans->trans_data;
    int host_addr;
    int int_port_num;
    char *host_name;
    char my_host_name[256];
    static int IP = 0;

    get_IP_config(my_host_name, sizeof(host_name), &IP, NULL, NULL, NULL,
		  NULL, svc->trace_out, (void *)cm);

    if (IP == 0) {
	if (IP == 0) IP = INADDR_LOOPBACK;
    }
    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "CMself check TCP/IP transport found no IP_HOST attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_addr)) {
	svc->trace_out(cm, "CMself check TCP/IP transport found no IP_ADDR attribute");
	if (host_name == NULL) return 0;
	host_addr = 0;
    }
    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "CMself check TCP/IP transport found no IP_PORT attribute");
	return 0;
    }
    if (host_name && (strcmp(host_name, my_host_name) != 0)) {
	svc->trace_out(cm, "CMself check - Hostnames don't match");
	return 0;
    }
    if (host_addr && (IP != host_addr)) {
	svc->trace_out(cm, "CMself check - Host IP addrs don't match, %lx, %lx", IP, host_addr);
	return 0;
    }
    int port_match = 0;
    for (int i = 0; i < sd->listen_count; i++) {
	if (int_port_num == sd->listen_ports[i]) {
	    port_match = sd->listen_ports[i];
	}
    }
    if (!port_match) {
	svc->trace_out(cm, "CMself check - Ports don't match, %d, %d", int_port_num, port_match);
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

extern int
libcmsockets_LTX_connection_eq(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, socket_conn_data_ptr scd)
{

    int int_port_num;
    int requested_IP = -1;
    char *host_name = NULL;

    if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & host_name)) {
	svc->trace_out(cm, "TCP/IP transport found no IP_HOST attribute");
    }
    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & int_port_num)) {
	svc->trace_out(cm, "Conn Eq TCP/IP transport found no IP_PORT attribute");
	return 0;
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *)(intptr_t) & requested_IP)) {
	svc->trace_out(cm, "TCP/IP transport found no IP_ADDR attribute");
    }
    if (requested_IP == -1) {
	check_host(host_name, (void *) &requested_IP);
	requested_IP = ntohl(requested_IP);
	svc->trace_out(cm, "IP translation for hostname %s is %x", host_name,
		       requested_IP);
    }

    svc->trace_out(cm, "Socket Conn_eq comparing IP/ports %x/%d and %x/%d",
		   scd->remote_IP, scd->remote_contact_port,
		   requested_IP, int_port_num);
    if ((scd->remote_IP == requested_IP) &&
	(scd->remote_contact_port == int_port_num)) {
	svc->trace_out(cm, "Socket Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "Socket Conn_eq returning FALSE");
    return 0;
}

/* 
 * Create an IP socket for connection from other CMs
 */
extern attr_list
libcmsockets_LTX_non_blocking_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_info)
{
    socket_client_data_ptr sd = trans->trans_data;
    unsigned int length;
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    SOCKET conn_sock = 0;
    int attr_port_num = 0;
    u_short port_num = 0;
    int port_range_low, port_range_high;
    int use_hostname = 0;
    int IP;
    char host_name[256];

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    /* 
     *  Check to see if a bind to a specific port was requested
     */
    if (listen_info != NULL
	&& !query_attr(listen_info, CM_IP_PORT,
		       NULL, (attr_value *)(intptr_t) & attr_port_num)) {
	port_num = 0;
    } else {
	if (attr_port_num > USHRT_MAX || attr_port_num < 0) {
	    fprintf(stderr, "Requested port number %d is invalid\n", attr_port_num);
	    return NULL;
	}
	port_num = attr_port_num;
    }

    svc->trace_out(cm, "CMSocket begin listen, requested port %d", attr_port_num);
    get_IP_config(host_name, sizeof(host_name), &IP, &port_range_low, &port_range_high, 
		  &use_hostname, listen_info, svc->trace_out, (void *)cm);

    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = htons(port_num);

    for (int i = 0; i < sd->listen_count; i++) {
	if ((sd->listen_ports[i] == port_num) ||
	    ((sd->listen_ports[i] >= port_range_low) &&
	     (sd->listen_ports[i] <= port_range_high))) {
	    conn_sock = sd->listen_fds[i];
	}
    }
    if (!conn_sock) {
	conn_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (conn_sock == SOCKET_ERROR) {
	    fprintf(stderr, "Cannot open INET socket\n");
	    return NULL;
	}
	if (sock_addr.sin_port != 0) {
	    /* specific port requested. set REUSEADDR, REUSEPORT because previous server might have died badly */
	    if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
			   sizeof(sock_opt_val)) != 0) {
		fprintf(stderr, "Failed to set REUSEADDR on INET socket before bind\n");
		perror("setsockopt(SO_REUSEADDR) failed");
		return NULL;
	    }
#ifdef SO_REUSEPORT
	    sock_opt_val = 1;
	    if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEPORT, (const char*)&sock_opt_val, sizeof(sock_opt_val)) != 0) {
		fprintf(stderr, "Failed to set REUSEADDR on INET socket before bind\n");
		perror("setsockopt(SO_REUSEPORT) failed");
		return NULL;
	    }
#endif
	    svc->trace_out(cm, "CMSocket trying to bind selected port %d", port_num);
	    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
		     sizeof sock_addr) == SOCKET_ERROR) {
		fprintf(stderr, "Cannot bind INET socket\n");
		return NULL;
	    }
	} else if (port_range_high == -1) {
	    /* bind to any port, range unconstrained */
	    sock_addr.sin_port = 0;
	    svc->trace_out(cm, "CMSocket trying to bind to any available port");
	    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
		     sizeof sock_addr) == SOCKET_ERROR) {
		fprintf(stderr, "Cannot bind INET socket\n");
		return NULL;
	    }
	} else {
	    long seedval = (long) time(NULL) + getpid();
	    /* port num is free.  Constrain to range to standards */
	    int size = port_range_high - port_range_low;
	    int tries = 30;
	    int result = SOCKET_ERROR;
	    srand(seedval);
	    while (tries > 0) {
		int target = port_range_low + (rand() % size);
		sock_addr.sin_port = htons(target);
		svc->trace_out(cm, "CMSocket trying to bind port %d", target);
		result = bind(conn_sock, (struct sockaddr *) &sock_addr,
			      sizeof sock_addr);
		tries--;
		if (result != SOCKET_ERROR) tries = 0;
		if (tries%5 == 4) {
		    /* try reseeding in case we're in sync with another process */
		    srand((int)time(NULL) + (int)getpid());
		}
		if (tries == 20) {
		    /* damn, tried a lot, increase the range (This might violate specified range) */
		    size *= 10;
		}
		if (tries == 10) {
		    /* damn, tried a lot more, increase the range (This might violate specified range) */
		    size *= 10;
		}
	    }
	    if (result == SOCKET_ERROR) {
		fprintf(stderr, "Cannot bind INET socket\n");
		return NULL;
	    }
	}
	/* begin listening for conns and set the backlog */
	if (listen(conn_sock, FD_SETSIZE)) {
	    fprintf(stderr, "listen failed\n");
	    return NULL;
	}
	svc->trace_out(cm, "CMSockets Adding socket_accept_conn as action on fd %d", conn_sock);
	svc->fd_add_select(cm, conn_sock, (select_list_func)socket_accept_conn,
			   (void *) trans, (void *) (intptr_t)conn_sock);

	length = sizeof(sock_addr);
	if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
	    fprintf(stderr, "Cannot get socket name\n");
	    return NULL;
	}
	sd->listen_fds = realloc(sd->listen_fds,
				 sizeof(SOCKET)*(sd->listen_count+1));
	sd->listen_ports = realloc(sd->listen_ports,
				 sizeof(int)*(sd->listen_count+1));
	sd->listen_fds[sd->listen_count] = conn_sock;
	sd->listen_ports[sd->listen_count] = ntohs(sock_addr.sin_port);
	sd->listen_count++;
    } else {
	length = sizeof(sock_addr);
	if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
	    fprintf(stderr, "Cannot get socket name\n");
	    return NULL;
	}
	svc->trace_out(cm, "CMSockets reusing prior listen, fd %d, port %d\n", conn_sock, ntohs(sock_addr.sin_port));
    }
    /* set the port num as one we can be contacted at */
    {
	int int_port_num = ntohs(sock_addr.sin_port);
	attr_list ret_list;

	svc->trace_out(cm, "CMSocket listen succeeded on port %d, fd %d",
		       int_port_num, conn_sock);
	ret_list = create_attr_list();

	if (sd->hostname != NULL)
	    svc->free_func(sd->hostname);
	sd->hostname = strdup(host_name);
	if ((IP != 0) && (!use_hostname)) {
	    add_attr(ret_list, CM_IP_ADDR, Attr_Int4,
		     (attr_value) (intptr_t)IP);
	}
	if ((getenv("CMSocketsUseHostname") != NULL) || 
	    use_hostname) {
	    add_attr(ret_list, CM_IP_HOSTNAME, Attr_String,
		     (attr_value) strdup(host_name));
	} else if (IP == 0) {
	    add_attr(ret_list, CM_IP_ADDR, Attr_Int4, 
		     (attr_value)INADDR_LOOPBACK);
	}
	add_attr(ret_list, CM_IP_PORT, Attr_Int4,
		 (attr_value) (intptr_t)int_port_num);

	return ret_list;
    }
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

extern void
libcmsockets_LTX_set_write_notify(transport_entry trans, CMtrans_services svc, socket_conn_data_ptr scd, int enable)
{
    if (enable != 0) {
	svc->fd_write_select(trans->cm, scd->fd, (select_list_func) trans->write_possible,
			     (void *)trans, (void *) scd->conn);
    } else {
	/* remove entry */
	svc->fd_write_select(trans->cm, scd->fd, NULL, NULL, NULL);
    }	
}


static void
set_block_state(CMtrans_services svc, socket_conn_data_ptr scd,
		socket_block_state needed_block_state)
{
#ifndef _MSC_VER
    int fdflags = fcntl(scd->fd, F_GETFL, 0);
    if (fdflags == -1) {
	perror("getflags\n");
	return;
    }
    if ((needed_block_state == Block) && (scd->block_state == Non_Block)) {
	fdflags &= ~O_NONBLOCK;
	if (fcntl(scd->fd, F_SETFL, fdflags) == -1) 
	    perror("fcntl block");
	scd->block_state = Block;
	svc->trace_out(scd->sd->cm, "CMSocket switch fd %d to blocking",
		       scd->fd);
    } else if ((needed_block_state == Non_Block) && 
	       (scd->block_state == Block)) {
	fdflags |= O_NONBLOCK;
	if (fcntl(scd->fd, F_SETFL, fdflags) == -1) 
	    perror("fcntl nonblock");
	scd->block_state = Non_Block;
	svc->trace_out(scd->sd->cm, "CMSocket switch fd %d to nonblocking",
		       scd->fd);
    }
#else
#endif
}

extern ssize_t
libcmsockets_LTX_read_to_buffer_func(CMtrans_services svc, socket_conn_data_ptr scd, void *buffer, ssize_t requested_len, int non_blocking)
{
    ssize_t left, iget;
#ifndef _MSC_VER
    // GSE
    int fdflags = fcntl(scd->fd, F_GETFL, 0);
    if (fdflags == -1) {
	perror("getflags\n");
	return -1;
    }
#endif
    if (scd->block_state == Block) {
	svc->trace_out(scd->sd->cm, "CMSocket fd %d state block", scd->fd);
    } else {
	svc->trace_out(scd->sd->cm, "CMSocket fd %d state nonblock", scd->fd);
    }
    svc->trace_out(scd->sd->cm, "CMSocket read of %zd bytes on fd %d, non_block %d", requested_len,
		   scd->fd, non_blocking);
    if (non_blocking && (scd->block_state == Block)) {
	svc->trace_out(scd->sd->cm, "CMSocket switch to non-blocking fd %d",
		       scd->fd);
	set_block_state(svc, scd, Non_Block);
    }
    iget = read(scd->fd, (char *) buffer, (int)requested_len);
    if ((iget == -1) || (iget == 0)) {
	int lerrno = errno;
	if ((lerrno != EWOULDBLOCK) &&
	    (lerrno != EAGAIN) &&
	    (lerrno != EINTR)) {
	    /* serious error */
	    svc->trace_out(scd->sd->cm, "CMSocket iget was -1, errno is %d, returning 0 for read",
			   lerrno);
	    return -1;
	} else {
	    if (non_blocking) {
		svc->trace_out(scd->sd->cm, "CMSocket iget was -1, would block, errno is %d",
			   lerrno);
		return 0;
	    }
	    return -1;
	}
    }
    left = requested_len - iget;
    while (left > 0) {
	int lerrno;
	iget = read(scd->fd, (char *) buffer + requested_len - left,
		    (int)left);
	lerrno = errno;
	if (iget == -1) {
	    if ((lerrno != EWOULDBLOCK) &&
		(lerrno != EAGAIN) &&
		(lerrno != EINTR)) {
		/* serious error */
		svc->trace_out(scd->sd->cm, "CMSocket iget was -1, errno is %d, returning %d for read", 
			   lerrno, requested_len - left);
		return (requested_len - left);
	    } else {
		iget = 0;
		if (!non_blocking && (scd->block_state == Non_Block)) {
		    svc->trace_out(scd->sd->cm, "CMSocket switch to blocking fd %d",
				   scd->fd);
		    set_block_state(svc, scd, Block);
		}
	    }
	} else if (iget == 0) {
	    svc->trace_out(scd->sd->cm, "CMSocket iget was 0, errno is %d, returning %d for read", 
			   lerrno, requested_len - left);
	    return requested_len - left;	/* end of file */
	}
	left -= iget;
    }
    return requested_len;
}


#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif

#ifndef HAVE_WRITEV
static
ssize_t
writev(fd, iov, iovcnt)
int fd;
struct iovec *iov;
int iovcnt;
{
    ssize_t wrote = 0;
    int i;
    for (i = 0; i < iovcnt; i++) {
	size_t left = iov[i].iov_len;
	ssize_t iget = 0;

	while (left > 0) {
	    errno = 0;
	    size_t this_write = left;
	    char *this_base = ((char*)iov[i].iov_base) + iov[i].iov_len - left;
	    iget = write(fd, this_base, (int)this_write);
	    if (iget == -1) {
		int lerrno = errno;
		if ((lerrno != EWOULDBLOCK) &&
		    (lerrno != EAGAIN) &&
		    (lerrno != EINTR)) {
		    /* serious error */
		    return -1;
		} else {
		    iget = 0;
		}
	    }
	    left -= iget;
	}
	wrote += iov[i].iov_len;
    }
    return wrote;
}
#endif

int long_writev(CMtrans_services svc, socket_conn_data_ptr scd, void *iovs, int iovcnt)
{
    assert(0);   // for right now, don't try this
    return 0;
}

#ifndef MAX_RW_COUNT
// Not actually defined outside the kernel as far as I know  - GSE
#define MAX_RW_COUNT 0x7ffff000
#endif

extern ssize_t
libcmsockets_LTX_writev_func(CMtrans_services svc, socket_conn_data_ptr scd, void *iovs, int iovcnt, attr_list attrs)
{
    SOCKET fd = scd->fd;
    ssize_t left = 0;
    ssize_t iget = 0;
    ssize_t iovleft, i;
    iovleft = iovcnt;
    struct iovec * iov = (struct iovec*) iovs;
    /* sum lengths */
    for (i = 0; i < iovcnt; i++) {
	left += iov[i].iov_len;
    }

    svc->trace_out(scd->sd->cm, "CMSocket writev of %zd bytes on fd %d",
		   left, fd);
    if (left > MAX_RW_COUNT) {
	// more to write than unix lets us do in one call
	return long_writev(svc, scd, iovs, iovcnt);
    }
		    
    while (left > 0) {
	size_t write_count = iovleft;
	if (write_count > IOV_MAX)
	    write_count = IOV_MAX;
	iget = writev(fd, (struct iovec *) &iov[iovcnt - iovleft],
		      write_count);
	if (iget == -1) {
	    svc->trace_out(scd->sd->cm, "	writev failed, errno was %d", errno);
	    if ((errno != EWOULDBLOCK) && (errno != EAGAIN)) {
		/* serious error */
		return (iovcnt - iovleft);
	    } else {
		if (errno == EWOULDBLOCK) {
		    svc->trace_out(scd->sd->cm, "CMSocket writev blocked - switch to blocking fd %d",
				   scd->fd);
		    set_block_state(svc, scd, Block);
		}
		iget = 0;
	    }
	}
	if (iget == left) {
	    return iovcnt;
	}
	svc->trace_out(scd->sd->cm, "	writev partial success, %d bytes written", iget);
	left -= iget;
	while (iget > 0) {
	    iget -= iov[iovcnt - iovleft].iov_len;
	    iovleft--;
	}

	if (iget < 0) {
	    /* 
	     * Only part of the last block was written.  Modify IO 
	     * vector to indicate the remaining block to be written.
	     */
	    /* restore iovleft and iget to cover remaining block */
	    iovleft++;
	    iget += iov[iovcnt - iovleft].iov_len;

	    /* adjust count down and base up by number of bytes written */
	    iov[iovcnt - iovleft].iov_len -= iget;
	    iov[iovcnt - iovleft].iov_base =
		(char *) (iov[iovcnt - iovleft].iov_base) + iget;
	}
    }
    return iovcnt;
}

/* non blocking version */
extern ssize_t
libcmsockets_LTX_NBwritev_func(CMtrans_services svc, socket_conn_data_ptr scd, void *iovs, int iovcnt, attr_list attrs)
{
    SOCKET fd = scd->fd;
    ssize_t init_bytes, left = 0;
    ssize_t iget = 0;
    ssize_t iovleft, i;
    struct iovec * iov = (struct iovec*) iovs;
    iovleft = iovcnt;

    /* sum lengths */
    for (i = 0; i < iovcnt; i++)
	left += iov[i].iov_len;

    init_bytes = left;

    svc->trace_out(scd->sd->cm, "CMSocket Non-blocking writev of %zd bytes on fd %d",
		   left, fd);
    set_block_state(svc, scd, Non_Block);
    while (left > 0) {
	ssize_t write_count = iovleft;
	ssize_t this_write_bytes = 0;
	if (write_count > IOV_MAX)
	    write_count = IOV_MAX;
	for (i = 0; i < write_count; i++)
	    this_write_bytes += iov[i].iov_len;

	iget = writev(fd, (struct iovec *) &iov[iovcnt - iovleft],
		      write_count);
	if (iget == -1) {
	    svc->trace_out(scd->sd->cm, "CMSocket writev returned -1, errno %d",
		   errno);
	    if ((errno != EWOULDBLOCK) && (errno != EAGAIN)) {
		/* serious error */
		return -1;
	    } else {
		return init_bytes - left;
	    }
	}
	svc->trace_out(scd->sd->cm, "CMSocket writev returned %d", iget);
	left -= iget;
	if (iget != this_write_bytes) {
	    /* didn't write everything, the rest would block, return */
	    svc->trace_out(scd->sd->cm, "CMSocket blocked, return %d", 
			   init_bytes -left);
	    return init_bytes - left;
	}
	iovleft -= write_count;
    }
    return init_bytes - left;
}

static int socket_global_init = 0;

#ifdef HAVE_WINDOWS_H
/* Winsock init stuff, ask for ver 2.2 */
static WORD wVersionRequested = MAKEWORD(2, 2);
static WSADATA wsaData;
#endif

static void
free_socket_data(CManager cm, void *sdv)
{
    socket_client_data_ptr sd = (socket_client_data_ptr) sdv;
    CMtrans_services svc = sd->svc;
    if (sd->hostname != NULL)
	svc->free_func(sd->hostname);
    free_attr_list(sd->characteristics);
    for(int i = 0 ; i < sd->listen_count; i++) {
	close(sd->listen_fds[i]);
    }
    svc->free_func(sd->listen_fds);
    svc->free_func(sd->listen_ports);
    svc->free_func(sd);
}

extern void *
libcmsockets_LTX_initialize(CManager cm, CMtrans_services svc, transport_entry trans)
{
    static int atom_init = 0;

    (void)trans;
    socket_client_data_ptr socket_data;
    svc->trace_out(cm, "Initialize TCP/IP Socket transport built in %s",
		   EVPATH_MODULE_BUILD_DIR);
    if (socket_global_init == 0) {
#ifdef HAVE_WINDOWS_H
	int nErrorStatus;
	/* initialize the winsock package */
	nErrorStatus = WSAStartup(wVersionRequested, &wsaData);
	if (nErrorStatus != 0) {
	    fprintf(stderr, "Could not initialize windows socket library!");
	    WSACleanup();
	    exit(-1);
	}
#endif
	/* 
	 * ignore SIGPIPE's  (these pop up when ports die.  we catch the 
	 * failed writes) 
	 */
#ifdef SIGPIPE
	signal(SIGPIPE, SIG_IGN);
#endif
    }
    if (atom_init == 0) {
	CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
	CM_IP_PORT = attr_atom_from_string("IP_PORT");
	CM_IP_ADDR = attr_atom_from_string("IP_ADDR");
	CM_FD = attr_atom_from_string("CONNECTION_FILE_DESCRIPTOR");
	CM_THIS_CONN_PORT = attr_atom_from_string("THIS_CONN_PORT");
	CM_PEER_CONN_PORT = attr_atom_from_string("PEER_CONN_PORT");
	CM_PEER_IP = attr_atom_from_string("PEER_IP");
	CM_PEER_HOSTNAME = attr_atom_from_string("PEER_HOSTNAME");
	(void)CM_PEER_HOSTNAME;
	CM_PEER_LISTEN_PORT = attr_atom_from_string("PEER_LISTEN_PORT");
	CM_TRANSPORT_RELIABLE = attr_atom_from_string("CM_TRANSPORT_RELIABLE");
	atom_init++;
    }
    socket_data = svc->malloc_func(sizeof(struct socket_client_data));
    socket_data->cm = cm;
    socket_data->hostname = NULL;
    socket_data->svc = svc;
    socket_data->characteristics = create_attr_list();
    socket_data->listen_count = 0;
    socket_data->listen_fds = malloc(sizeof(SOCKET));
    socket_data->listen_ports = malloc(sizeof(int));
    add_int_attr(socket_data->characteristics, CM_TRANSPORT_RELIABLE, 1);
    svc->add_shutdown_task(cm, free_socket_data, (void *) socket_data, FREE_TASK);
    return (void *) socket_data;
}

extern attr_list
libcmsockets_LTX_get_transport_characteristics(transport_entry trans, CMtrans_services svc,
					       void* vsd)
{
    struct socket_client_data * sd = (struct socket_client_data *) vsd;
    add_ref_attr_list(sd->characteristics);
    return sd->characteristics;
}

extern transport_entry
cmsockets_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport = svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("sockets");
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)libcmsockets_LTX_initialize;
    transport->listen = (CMTransport_listen_func)libcmsockets_LTX_non_blocking_listen;
    transport->initiate_conn = (CMTransport_conn_func)libcmsockets_LTX_initiate_conn;
    transport->self_check = (CMTransport_self_check_func)libcmsockets_LTX_self_check;
    transport->connection_eq = (CMTransport_connection_eq_func)libcmsockets_LTX_connection_eq;
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)libcmsockets_LTX_shutdown_conn;
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)libcmsockets_LTX_read_to_buffer_func;
    transport->read_block_func = (CMTransport_read_block_func)NULL;
    transport->writev_func = (CMTransport_writev_func)libcmsockets_LTX_writev_func;
    transport->NBwritev_func = (CMTransport_writev_func)libcmsockets_LTX_NBwritev_func;
    
    transport->set_write_notify = (CMTransport_set_write_notify_func)    libcmsockets_LTX_set_write_notify;
    transport->get_transport_characteristics = (CMTransport_get_transport_characteristics) libcmsockets_LTX_get_transport_characteristics;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}
