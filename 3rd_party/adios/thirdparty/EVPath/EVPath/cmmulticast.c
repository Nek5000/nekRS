/***** Includes *****/
#include "config.h"
#include <sys/types.h>

#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <ws2ipdef.h>
#include <windows.h>
#define getpid()	_getpid()
#define close(x) closesocket(x)
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

typedef struct func_list_item {
    select_list_func func;
    void *arg1;
    void *arg2;
} FunctionListElement;

typedef struct multicast_transport_data {
    CManager cm;
    CMtrans_services svc;
} *multicast_transport_data_ptr;


#include <stdio.h>
#include "config.h"

#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 177)
#  pragma warning (disable: 2259)
#endif

#define MSGBUFSIZE 25600

typedef struct mcast_connection_data {
    int mcast_IP;
    int mcast_port;
    SOCKET input_fd;
    SOCKET output_fd;
    struct sockaddr_in output_addr;
    struct sockaddr_in my_addr;
    char read_buffer[MSGBUFSIZE];
    int read_buffer_len;
    int read_pointer;
    CMConnection conn;
    multicast_transport_data_ptr mtd;
} *mcast_conn_data_ptr;

#ifdef _MSC_VER
#define read(fd, buf, len) recv(fd, buf, len, 0)
#define write(fd, buf, len) send(fd, buf, len, 0)
#endif

static atom_t CM_MCAST_PORT = -1;
static atom_t CM_MCAST_ADDR = -1;
static atom_t CM_FD = -1;

static mcast_conn_data_ptr
create_mcast_conn_data(CMtrans_services svc)
{
    mcast_conn_data_ptr mcast_conn_data =
    svc->malloc_func(sizeof(struct mcast_connection_data));
    mcast_conn_data->mcast_port = -1;
    mcast_conn_data->input_fd = 0;
    mcast_conn_data->output_fd = 0;
    memset(&mcast_conn_data->my_addr, 0, sizeof(mcast_conn_data->my_addr));
    return mcast_conn_data;
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
libcmmulticast_LTX_shutdown_conn(CMtrans_services svc, mcast_conn_data_ptr mcd)
{
    svc->fd_remove_select(mcd->mtd->cm, mcd->input_fd);
    close(mcd->input_fd);
    close(mcd->output_fd);
    free(mcd);
}


#include "qual_hostname.c"

static SOCKET
initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, mcast_conn_data_ptr mcast_conn_data, attr_list conn_attr_list, int no_more_redirect)
{
    int one = 1;
    SOCKET input_fd, output_fd;
    int int_port_num;
    u_short port_num;
    multicast_transport_data_ptr mtd = (multicast_transport_data_ptr) trans->trans_data;
    static int mcast_ip = 0;
    struct sockaddr_in addr;
    struct sockaddr_in output_addr;
    struct ip_mreq mreq;
    char my_host_name[256];

    (void) no_more_redirect;
    if (!query_attr(attrs, CM_MCAST_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &mcast_ip)) {
	svc->trace_out(cm, "CMMulticast transport found no MCAST_ADDR attribute");
	/* wasn't there */
	mcast_ip = 0;
    } else {
	svc->trace_out(cm, "CMMulticast transport connect to mcast_IP %lx", mcast_ip);
    }
    if (mcast_ip == 0)
	return -1;

    if (!query_attr(attrs, CM_MCAST_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &int_port_num)) {
	svc->trace_out(cm, "CMMulticast transport found no MCAST_PORT attribute");
	return -1;
    } else {
	svc->trace_out(cm, "CMMulticast transport connect to port %d", int_port_num);
    }
    if ((input_fd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	perror("socket");
	exit(1);
    }
    port_num = int_port_num;

    get_qual_hostname(cm, my_host_name, sizeof(my_host_name) - 1, svc, NULL, NULL);

    memset(&addr, 0, sizeof(addr));
    addr.sin_family = (unsigned short) AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);	/* N.B.: differs from *
						 * sender */
    addr.sin_port = (unsigned short) htons(port_num);
#ifdef SO_REUSEPORT
    if (setsockopt(input_fd, SOL_SOCKET, SO_REUSEPORT, (char *) &one, sizeof(one)) == -1) {
	perror("setsockopt reuseport");
    }
#else
    if (setsockopt(input_fd, SOL_SOCKET, SO_REUSEADDR, (char *) &one, sizeof(one)) == -1) {
	perror("setsockopt reuseport");
    }
#endif
    /* bind to receive address */
    if (bind(input_fd, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
	perror("bind");
	exit(1);
    }
    /* use setsockopt() to request that the kernel join a multicast group */
    mreq.imr_multiaddr.s_addr = htonl(mcast_ip);
    mreq.imr_interface.s_addr = htonl(INADDR_ANY);
    if (setsockopt(input_fd, IPPROTO_IP, IP_ADD_MEMBERSHIP, (const char*) & mreq, sizeof(mreq)) < 0) {
	perror("setsockopt");
	exit(1);
    }
    /* create what looks like an ordinary UDP socket */
    if ((output_fd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
	perror("socket");
	exit(1);
    }
    /* set up destination address */
    memset(&addr, 0, sizeof(addr));
    memset(&output_addr, 0, sizeof(output_addr));
    output_addr.sin_family = AF_INET;
    output_addr.sin_addr.s_addr = htonl(mcast_ip);
    output_addr.sin_port = (unsigned short) htons(port_num);

    svc->trace_out(cm, "--> Connection established");

    mcast_conn_data->mcast_IP = mcast_ip;
    mcast_conn_data->mcast_port = int_port_num;
    mcast_conn_data->input_fd = input_fd;
    mcast_conn_data->output_fd = output_fd;
    mcast_conn_data->output_addr = output_addr;
    mcast_conn_data->mtd = mtd;

    add_attr(conn_attr_list, CM_FD, Attr_Int4,
	     (attr_value) (long)input_fd);
    return input_fd;
}

static void
libcmmulticast_data_available(void *vtrans, void *vmcd)
{
    transport_entry trans = vtrans;
    int nbytes;
    mcast_conn_data_ptr mcd = vmcd;
    struct sockaddr_in addr;
    unsigned int addrlen = sizeof(addr);

    char *msgbuf = &mcd->read_buffer[0];
    if ((nbytes = recvfrom(mcd->input_fd, msgbuf, MSGBUFSIZE, 0,
			   (struct sockaddr *) &addr, &addrlen)) < 0) {
	perror("recvfrom");
	exit(1);
    }
    if (mcd->my_addr.sin_port == 0) {
	unsigned int nl;
	int IP = get_self_ip_addr(NULL, mcd->mtd->svc);
	nl = sizeof(struct sockaddr_in);
	if (getsockname(mcd->output_fd, (struct sockaddr *) &mcd->my_addr, &nl) != 0)
	    perror("getsockname");


	mcd->my_addr.sin_addr.s_addr = htonl(IP);
    }
    if (memcmp(&addr, &mcd->my_addr, sizeof(addr)) == 0) {
	/* ignore stuff we've sent */
	return;
    }
    mcd->read_buffer_len = nbytes;
    mcd->read_pointer = 0;
    /* kick this upstairs */
    trans->data_available(vtrans, mcd->conn);
}

/* 
 * Initiate a connection to a multicast group.
 */
extern CMConnection
libcmmulticast_LTX_initiate_conn(CManager cm,CMtrans_services svc, transport_entry trans, attr_list attrs)
{
    mcast_conn_data_ptr mcast_conn_data = create_mcast_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;
    SOCKET sock;

    if ((sock = initiate_conn(cm, svc, trans, attrs, mcast_conn_data, conn_attr_list, 0)) < 0)
	return NULL;

    conn = svc->connection_create(trans, mcast_conn_data, conn_attr_list);
    mcast_conn_data->conn = conn;

    svc->trace_out(cm, "CMMulticast Adding libcmmulticast_data_available as action on fd %d", sock);
    svc->fd_add_select(cm, sock, (select_list_func) libcmmulticast_data_available,
		       (void *) trans, (void *) mcast_conn_data);

    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * This only makes sense for connection-based transports, not multicast.
 */
extern int
libcmmulticast_LTX_self_check(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{
    return 0;
}

extern int
libcmmulticast_LTX_connection_eq(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, mcast_conn_data_ptr mcd)
{

    int int_port_num;
    int requested_IP = -1;

    if (!query_attr(attrs, CM_MCAST_PORT, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &int_port_num)) {
	svc->trace_out(cm, "Conn Eq CMMulticast transport found no MCAST_PORT attribute");
	return 0;
    }
    if (!query_attr(attrs, CM_MCAST_ADDR, /* type pointer */ NULL,
    /* value pointer */ (attr_value *) (intptr_t) &requested_IP)) {
	svc->trace_out(cm, "CMMulticast transport found no MCAST_ADDR attribute");
    }
    svc->trace_out(cm, "CMMulticast Conn_eq comparing IP/ports %x/%d and %x/%d",
		   mcd->mcast_IP, mcd->mcast_port,
		   requested_IP, int_port_num);
    if ((mcd->mcast_IP == requested_IP) &&
	(mcd->mcast_port == int_port_num)) {
	svc->trace_out(cm, "CMMulticast Conn_eq returning TRUE");
	return 1;
    }
    svc->trace_out(cm, "CMMulticast Conn_eq returning FALSE");
    return 0;
}


extern attr_list
libcmmulticast_LTX_non_blocking_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_info)
{
    /* meaningless in muticast */
    return NULL;
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
libcmmulticast_LTX_read_func(CMtrans_services svc, mcast_conn_data_ptr mcd, int requested_len, int *actual_len)
{
    char *ret = &mcd->read_buffer[mcd->read_pointer];
    *actual_len = requested_len;
    mcd->read_pointer += requested_len;
    return ret;
}

#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif
#ifdef _MSC_VER
#define msghdr _WSAMSG
#include <ws2def.h>
#endif
extern int
libcmmulticast_LTX_writev_func(CMtrans_services svc, mcast_conn_data_ptr mcd, struct iovec *iov, int iovcnt, attr_list attrs)
{
    SOCKET fd = mcd->output_fd;
#ifndef _MSC_VER
    // no real equivalent on windows
    struct sockaddr_in addr = mcd->output_addr;
    struct msghdr msg;
    svc->trace_out(mcd->mtd->cm, "CMMcast writev of %d vectors on fd %d",
		   iovcnt, fd);
    memset(&msg, 0, sizeof(msg));
    msg.msg_name = (void*)&addr;
    msg.msg_namelen = sizeof(addr);
    msg.msg_iov = &iov[0];
    msg.msg_iovlen = iovcnt;
    if (sendmsg(fd, &msg, 0) < 0) {
	perror("write sendmsg");
	exit(1);
    }
#endif
    if (mcd->my_addr.sin_port == 0) {
	unsigned int nl;
	int IP = get_self_ip_addr(NULL, svc);
	nl = sizeof(struct sockaddr_in);
	if (getsockname(fd, (struct sockaddr *) &mcd->my_addr, &nl) != 0)
	    perror("getsockname");


	mcd->my_addr.sin_addr.s_addr = htonl(IP);
    }
    return iovcnt;
}

#ifdef HAVE_WINDOWS_H
static int socket_global_init = 0;
/* Winsock init stuff, ask for ver 1.1 */
static WORD wVersionRequested = MAKEWORD(1, 1);
static WSADATA wsaData;
#endif

static void
free_mcast_data(CManager cm, void *mtdv)
{
    multicast_transport_data_ptr mtd = (multicast_transport_data_ptr) mtdv;
    CMtrans_services svc = mtd->svc;
    svc->free_func(mtd);
}

extern void *
libcmmulticast_LTX_initialize(CManager cm, CMtrans_services svc)
{
    static int atom_init = 0;

    multicast_transport_data_ptr mcast_data;
    svc->trace_out(cm, "Initialize CMMulticast transport");
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
	CM_MCAST_PORT = attr_atom_from_string("MCAST_PORT");
	CM_MCAST_ADDR = attr_atom_from_string("MCAST_ADDR");
	CM_FD = attr_atom_from_string("CM_FD");
	atom_init++;
    }
    mcast_data = svc->malloc_func(sizeof(struct multicast_transport_data));
    mcast_data->cm = cm;
    mcast_data->svc = svc;
    svc->add_shutdown_task(cm, free_mcast_data, (void *) mcast_data, FREE_TASK);
    return (void *) mcast_data;
}
