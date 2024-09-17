/***** Includes *****/
#include "config.h"
#include <sys/types.h>

#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <windows.h>
#define getpid()    _getpid()
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
#include <netinet/tcp.h>
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
#include "cm_internal.h"

#include <infiniband/verbs.h>
#include <sys/queue.h>
#include <stdlib.h>

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif

#define LISTSIZE 1024
#define _WITH_IB_
#define PIGGYBACK 1025*100

#ifdef _WITH_IB_


#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 2259)
#  pragma warning (disable: 177)
#endif


//all the message types have a queue associated with them
static char *msg_string[] = {"Request", "Response", "Piggyback"};
enum {msg_request = 0, msg_response = 1, msg_piggyback = 2} msg_type;

struct request
{
    int magic;    
    uint32_t length;    
    uint64_t request_ID;
};


struct response
{
    uint64_t remote_addr;
    uint32_t rkey;
    uint32_t max_length;    
    uint64_t request_ID;
};

struct piggyback
{
    uint32_t total_length;
    char  body[4];
};


struct control_message
{
    int type;
    union {
	struct request req;
	struct response resp;
	struct piggyback pb;
    } u;
};

struct ibparam
{
	int lid;
	int psn;
	int qpn;
	int port;
	/*anything else? */
};

int page_size = 0;
int perftrace =0;


#define ptr_from_int64(p) (void *)(unsigned long)(p)
#define int64_from_ptr(p) (u_int64_t)(unsigned long)(p)


typedef struct ib_client_data {
	CManager cm;
	char *hostname;
	int listen_port;
	CMtrans_services svc;
	int lid;
	int qpn;
	int psn;
	int port;
	struct ibv_device *ibdev;
	struct ibv_context *context;
	struct ibv_comp_channel *send_channel;
	struct ibv_comp_channel *recv_channel;
	struct ibv_pd *pd;
	struct ibv_cq *recv_cq;
	struct ibv_cq *send_cq;
	struct ibv_srq *srq;
	int max_sge;    
} *ib_client_data_ptr;


typedef struct notification
{
	int done;
	struct ibv_mr *mr;    
	struct ibv_send_wr wr;
	struct ibv_sge sg;    
	struct ibv_recv_wr rwr;    
	struct ibv_recv_wr *badrwr;    
}notify;

typedef struct remote_info
{
    int mrlen;
    notify isDone;
    struct ibv_mr **mrlist;
    struct ibv_send_wr *wr; 
    CMcompletion_notify_func notify_func;
    void *notify_client_data;
}rinfo;


typedef struct ib_connection_data {
	char *remote_host;
	int remote_IP;
	int remote_contact_port;
	int fd;
	void *read_buffer;
	int read_buffer_len;
	struct tbuffer_ *tb;    
	ib_client_data_ptr sd;
	CMConnection conn;
	struct ibv_qp *dataqp;    
	notify isDone;
	int infocount;
	rinfo infolist[10];
	int max_imm_data;   
} *ib_conn_data_ptr;

typedef struct tbuffer_
{
	CMbuffer buf;
	struct ibv_mr *mr;
	ib_conn_data_ptr scd;
	uint64_t size;
	uint64_t offset;
	struct tbuffer_ *parent;
	int childcount;    
        int inuse;
	LIST_ENTRY(tbuffer_) entries;    
}tbuffer;

LIST_HEAD(tblist, tbuffer_) memlist;
LIST_HEAD(inuselist, tbuffer_) uselist;


char *ibv_wc_opcode_str[] = {
	"IBV_WC_SEND",
	"IBV_WC_RDMA_WRITE",
	"IBV_WC_RDMA_READ",
	"IBV_WC_COMP_SWAP",
	"IBV_WC_FETCH_ADD",
	"IBV_WC_BIND_MW",
	"IBV_WC_RECV_RDMA_WITH_IMM"
};

static struct ibv_qp * initqp(ib_conn_data_ptr ib_conn_data,
                              ib_client_data_ptr sd);
static int connectqp(ib_conn_data_ptr ib_conn_data,
                     ib_client_data_ptr sd,
                     struct ibparam lparam,
                     struct ibparam rparam);
static struct ibv_mr ** regblocks(ib_client_data_ptr sd,
                                  struct iovec *iovs, int iovcnt, int flags,
                                  int *mrlen);
static tbuffer *findMemory(ib_conn_data_ptr scd, ib_client_data_ptr sd, 
                           CMtrans_services svc, int req_size);


static struct ibv_send_wr * createwrlist(ib_conn_data_ptr conn, 
                                         struct ibv_mr **mrlist,
                                         struct iovec *iovlist,
                                         int mrlen, int *wrlen); 

static int waitoncq(ib_conn_data_ptr scd,
                    ib_client_data_ptr sd,
                    CMtrans_services svc, struct ibv_cq *cq);

static inline int msg_offset()
{
	return ((int) (((char *) (&(((struct control_message *)NULL)->u.pb.body[0]))) - ((char *) NULL)));
}


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
static atom_t CM_NETWORK_POSTFIX = -1;
static atom_t CM_IP_PORT = -1;
static atom_t CM_IP_HOSTNAME = -1;
static atom_t CM_IP_ADDR = -1;
static atom_t CM_TRANSPORT = -1;

static double getlocaltime()
{
	struct timeval t;
	double dt;
	gettimeofday(&t, NULL);
	dt = (double) t.tv_usec / 1e6 + t.tv_sec;
	return dt;
}

static void free_func(void *cbd)
{
	tbuffer *tb = (tbuffer*)cbd;
	tbuffer *temp;
    
	tb->inuse = 0;
	if(tb->childcount == 0)
	{
		//we can merge upwards
		temp = tb->parent;
		if(!temp)
		{
			//is the top level - at the top it means we don't need to merge or do anything other than drop offset to 0
			tb->offset = 0;     
		}
		else
		{
			//is a child of someone - we merge up
			temp->size += tb->size;
			temp->childcount --;
        
			//we can free the CMbuffer now
			free(tb->buf);
        
			//now remove tb from list and free it
			LIST_REMOVE(tb, entries);
			free(tb);       
		}
	}
    
    
    
}

inline static struct ibv_device * IB_getdevice(char *name)
{
	struct ibv_device **dev_list;
	struct ibv_device *ib_dev;

	dev_list = ibv_get_device_list(NULL);
	if(!dev_list || !*(dev_list))
	{
		fprintf(stderr, "%s %d:%s - Couldn't get IB device list\n",
		        __FILE__, __LINE__, __FUNCTION__);
		return NULL;
	}
    
	if(name)
	{
		for(; (ib_dev= *dev_list); ++dev_list)
		{
			if(!strcmp(ibv_get_device_name(ib_dev), name))
			{
				break;
			}
		}
		if(!ib_dev)
		{
			fprintf(stderr, "%s %d:%s - Couldn't get IB device of name %s\n",
			        __FILE__, __LINE__, __FUNCTION__, name);
		}
	}
	else
		ib_dev = *dev_list; //return very first device so obtained

	return ib_dev; //could be null
}


static inline uint16_t get_local_lid(struct ibv_context *context, int port)
{
	struct ibv_port_attr attr;

	if (ibv_query_port(context, port, &attr))
		return 0;

	return attr.lid;
}

static int
check_host(char *hostname,void *sin_addr)
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

static ib_conn_data_ptr 
create_ib_conn_data(CMtrans_services svc)
{
	ib_conn_data_ptr ib_conn_data = svc->malloc_func(sizeof(struct ib_connection_data));
	memset(ib_conn_data, 0, sizeof(struct ib_connection_data));
	ib_conn_data->remote_host = NULL;
	ib_conn_data->remote_contact_port = -1;
	ib_conn_data->fd = 0;
	ib_conn_data->read_buffer = NULL;
	ib_conn_data->read_buffer_len = 0;

	return ib_conn_data;
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

static int internal_read_stream(CMtrans_services svc,
                                ib_conn_data_ptr scd,
                                int length,
                                unsigned char * buffer)
{
	//this function is called whenever we want to read
	//from the fd
	
	int left;
	int iget = 0;
	int requested_len;
	int lerrno;

	left = length;
	requested_len = left;
	
	while(left > 0)
	{
		iget = read(scd->fd, (char*)buffer+requested_len - left,
		            left);
		lerrno = errno;
		if (iget == -1) {
			if ((lerrno != EWOULDBLOCK) &&
			    (lerrno != EAGAIN) &&
			    (lerrno != EINTR)) {
				/* serious error */
				svc->trace_out(scd->sd->cm, "CMIB iget was -1, errno is %d, returning %d for read", 
				               lerrno, requested_len - left);
				return iget;
			}
		}
		if (iget == 0) {
		  return -1;
		}
		left -= iget;
	}

	return iget;

}

static int internal_write_piggyback(CMtrans_services svc,
                                    ib_conn_data_ptr scd,
                                    int length, struct iovec *iov, int iovcnt)
{
	//this function is only called if length < piggyback
	struct control_message *msg;
	char *point;
	int offset = msg_offset();
	int i;
	
	if(length >= PIGGYBACK)
	{
		//should never happen
		return -1;
	}

	msg = malloc(offset + length);
	memset(msg, 0, offset+length);
	msg->type = msg_piggyback;
	msg->u.pb.total_length = length + offset;
	point = &msg->u.pb.body[0];
	svc->trace_out(scd->sd->cm, "CMIB sending piggyback msg of length %d,", 
		       length);

	
	for (i = 0; i < iovcnt; i++)
	{
		memcpy(point, iov[i].iov_base, iov[i].iov_len);
		point += iov[i].iov_len;
	}

	if(write(scd->fd, msg, msg->u.pb.total_length) != msg->u.pb.total_length)
	{
		//error in writing the whole block
		svc->trace_out(scd->sd->cm, "CMIB error in writing to %d", scd->fd);
		free(msg);
		return -1;
	}
	
	free(msg);
	return 0;
}


static int internal_write_response(CMtrans_services svc,
                                   ib_conn_data_ptr scd,
                                   tbuffer *tb,
                                   int length,
				   int64_t request_ID)
{
	struct control_message msg;
	int iget = 0;

	msg.type = msg_response;
	if(tb != NULL)
	{
		msg.u.resp.remote_addr = int64_from_ptr(tb->buf->buffer);
		msg.u.resp.rkey = tb->mr->rkey;
		msg.u.resp.max_length = length;
		msg.u.resp.request_ID = request_ID;
	}
	else
	{
		msg.u.resp.remote_addr = 0;
		msg.u.resp.rkey = 0;
		msg.u.resp.max_length = 0;
		msg.u.resp.request_ID = request_ID;
	}

	iget = write(scd->fd, &msg, sizeof(struct control_message));
	if(iget != sizeof(struct control_message))
	{
		svc->trace_out(scd->sd->cm, "CMIB error in writing response to %d", scd->fd);
		return -1;
	}

	return 0;
}

static int internal_write_request(CMtrans_services svc,
                                  ib_conn_data_ptr scd,
                                  int length,
    				  void *request_ID)
{
	struct control_message msg;
	int iget = 0;

	msg.type = msg_request;
	msg.u.req.magic = 0xdeadbeef;
	msg.u.req.length = length;
	msg.u.req.request_ID = int64_from_ptr(request_ID);

	svc->trace_out(scd->sd->cm, "Doing internal write request, writing %d bytes", sizeof(struct control_message));
	iget = write(scd->fd, &msg, sizeof(struct control_message));
	if(iget != sizeof(struct control_message))
	{
		svc->trace_out(scd->sd->cm, "CMIB error in writing request to %d", scd->fd);
		return -1;
	}

	return 0;
}
	
static double write_t = 0;
static double da_t = 0;

static void
free_sg_list(struct ibv_sge *sg_list, int len, int free_data_elements)
{
    int i;
    if (free_data_elements) {
	for(i=0; i<len; i++) {
	    free(ptr_from_int64(sg_list[i].addr));
	}
    }
    free(sg_list);     
}

static int handle_response(CMtrans_services svc,
                           ib_conn_data_ptr scd,
                           struct control_message *msg)
{

    //read back response
    int iget = 0;
    struct ibv_send_wr *bad_wr;    
    int retval = 0;
    int i = 0;
    struct ibv_wc wc;
    struct response *rep;
    rinfo *write_request;
    int free_data_elements;
    
    rep = &msg->u.resp;
    write_request = ptr_from_int64(rep->request_ID);
    
    free_data_elements = (write_request->notify_func == NULL);
    if (write_request == NULL) {
	fprintf(stderr, "Failed to get work request - aborting write\n");
	return -0x01000;
    }

    write_request->wr->wr.rdma.remote_addr = rep->remote_addr;    
    write_request->wr->wr.rdma.rkey = rep->rkey;

    
    if(write_request->wr == NULL)
    {
	fprintf(stderr, "failed to get work request - aborting write\n");
	return -0x01000;    
    }
    
    retval = ibv_req_notify_cq(scd->sd->send_cq, 0);
    if(retval)
    {
	scd->sd->svc->trace_out(scd->sd->cm, "CMib notification request failed");
	return -1;  
    }
    
    if(rep->remote_addr == 0)
    {
	//error on the recieving side just send a -1 as the isdone and return
	scd->isDone.done = -1;
	
	retval = ibv_post_send(scd->dataqp, &scd->isDone.wr, 
			       &bad_wr);
        
	if(retval)
	{
	    //we got an error - ideally we'll fall through and post an error on the connection socket
	    svc->trace_out(scd->sd->cm, "CMib unable to notify over ib\n");
	}
	
	for(i = 0; i < write_request->mrlen && write_request->mrlist; i ++)
	{
	    ibv_dereg_mr(write_request->mrlist[i]);	    
	}

	if(write_request->wr)
	{       
	    free_sg_list(write_request->wr->sg_list, write_request->mrlen, free_data_elements);
	    free(write_request->wr);
	}
	
	//no we wait for the send even to get anoutput
	retval = waitoncq(scd, scd->sd, svc, scd->sd->send_cq);
	if(retval)
	{
	    svc->trace_out(scd->sd->cm, "Error while waiting\n");
	    return -1;      
	}
        
        
	
    }
    else
    {    
	retval = ibv_post_send(scd->dataqp, write_request->wr, &bad_wr);
	if(retval)
	{
	    svc->trace_out(scd->sd->cm, "CMIB unable to post send %d\n", retval);
	    //we can get the error from the *bad_wr
	    return retval;  
	    
	}
	
	retval = waitoncq(scd, scd->sd, svc, scd->sd->send_cq);
	if(retval)
	{
	    svc->trace_out(scd->sd->cm, "Error while waiting\n");
	    return -1;      
	}
	
    }
    
    
    do
    {
	//empty the poll cq 1 by 1
	iget = ibv_poll_cq(scd->sd->send_cq, 1, &wc);
	if(iget > 0 && wc.status == IBV_WC_SUCCESS && wc.opcode == IBV_WC_RDMA_WRITE)
	{
	    //send completeled for RDMA write
	    //we can break out after derigstering the memory
	    
	    //now post a send
	    scd->isDone.done = 0;
	    
	    retval = ibv_post_send(scd->dataqp, &scd->isDone.wr, 
				   &bad_wr);
	    
	    if(retval)
	    {
		//we got an error - ideally we'll fall through and post an error on the connection socket
		svc->trace_out(scd->sd->cm, "CMib unable to notify over ib\n");
		break;      
	    }
	    
#ifdef DO_DEREG_ON_FINISH
	    for(i = 0; i < write_request->mrlen; i ++)
	    {
		ibv_dereg_mr(write_request->mrlist[i]);       
		write_request->mrlist[i] = NULL;
	    }
	    
	    free_sg_list(write_request->wr->sg_list, write_request->mrlen, free_data_elements);
	    free(write_request->mrlist);
	    free(write_request->wr);
	    write_request->wr = NULL;
#endif

	    //no we wait for the send even to get anoutput
	    retval = waitoncq(scd, scd->sd, svc, scd->sd->send_cq);
	    if(retval)
	    {
		svc->trace_out(scd->sd->cm, "Error while waiting\n");
		return -1;      
	    }
	    
	    
	}
	else if(iget > 0 && wc.status == IBV_WC_SUCCESS && wc.opcode == IBV_WC_SEND)
	{
	    //cool beans - send completed we can go on with our life
	    
	    svc->trace_out(scd->sd->cm, "notification for send done %p\n", 
			   ptr_from_int64(wc.wr_id));
	    break;      
	}
	else if(iget == 0)
	{
	    //cq is empty - we shouldn't even be here!
	    continue;       
	}
	else
	{
	    svc->trace_out(scd->sd->cm, "error in polling queue %X %d %d", 
			   wc.wr_id, wc.status,  wc.opcode);
	    
	    return -1;      
	}
	
    }while(1);
    if (write_request->notify_func) {
	(write_request->notify_func)( write_request->notify_client_data);
    }
    svc->wake_any_pending_write(scd->conn);
    return 0;   
    
}

static void handle_request(CMtrans_services svc,
                           ib_conn_data_ptr scd,
                           struct control_message *msg)
{
	//handling the request message

	//first read the request message from the socket
	struct ibv_wc wc;
	tbuffer *tb;
	int retval = 0;
	struct request *req;

	req = &msg->u.req;
	
	tb = findMemory(scd, scd->sd, svc, req->length);
	if(tb == NULL)
	{
		svc->trace_out(scd->sd->cm, "Failed to get memory\n");
		internal_write_response(svc, scd, NULL, 0, 0);
		goto wait_on_q; 
	}

	scd->tb = tb;

	svc->set_pending_write(scd->conn);
	internal_write_response(svc, scd, tb, req->length, req->request_ID);

wait_on_q:
	//now the sender will start the data transfer. When it finishes he will 
	//issue a notification after which the data is in memory
	retval = waitoncq(scd, scd->sd, svc, scd->sd->recv_cq);    
	if(retval)
	{
		svc->trace_out(scd->sd->cm, "Error while waiting\n");
		return;     
	}

	//now we start polling the cq

	do
	{
		svc->trace_out(scd->sd->cm, "CMIB poll cq -start");    
		retval = ibv_poll_cq(scd->sd->recv_cq, 1, &wc);
		svc->trace_out(scd->sd->cm, "CMIB poll cq -end");    
		if(retval > 0 && wc.status == IBV_WC_SUCCESS && wc.opcode == IBV_WC_RECV)
		{
			//cool beans - send completed we can go on with our life
        
			//issue a reccieve so we don't run out
			retval = ibv_post_recv(scd->dataqp, &scd->isDone.rwr, 
			                       &scd->isDone.badrwr);
			if(retval)
			{
				scd->sd->svc->trace_out(scd->sd->cm, "CMib unable to post recv %d\n", retval);
			}
    
			break;      
		}
		else
		{
			svc->trace_out(scd->sd->cm, "Error polling for write completion\n");
			break;      
		}

		svc->trace_out(scd->sd->cm, "CMIB poll cq looping");    
    
	}while(1);

	scd->read_buffer = tb->buf;
	svc->trace_out(scd->sd->cm, "FIrst 16 bytes of receive buffer (len %d) are %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x \n", req->length, ((unsigned char*)tb->buf->buffer)[0], ((unsigned char*)tb->buf->buffer)[1], ((unsigned char*)tb->buf->buffer)[2], ((unsigned char*)tb->buf->buffer)[3], ((unsigned char*)tb->buf->buffer)[4], ((unsigned char*)tb->buf->buffer)[5], ((unsigned char*)tb->buf->buffer)[6], ((unsigned char*)tb->buf->buffer)[7], ((unsigned char*)tb->buf->buffer)[8], ((unsigned char*)tb->buf->buffer)[9], ((unsigned char*)tb->buf->buffer)[10], ((unsigned char*)tb->buf->buffer)[11], ((unsigned char*)tb->buf->buffer)[12], ((unsigned char*)tb->buf->buffer)[13], ((unsigned char*)tb->buf->buffer)[14], ((unsigned char*)tb->buf->buffer)[15]);
	scd->read_buffer_len = req->length;
	svc->trace_out(scd->sd->cm, "CMIB handle_request completed");    
	svc->wake_any_pending_write(scd->conn);
}

int handle_piggyback(CMtrans_services svc,
                      ib_conn_data_ptr scd,
                      struct control_message *msg)
{
	// read rest and deliver up
	int offset = msg_offset();
	int iget = 0;
	
	CMbuffer cb = svc->get_data_buffer(scd->sd->cm, msg->u.pb.total_length - offset);
	memcpy(cb->buffer, &(msg->u.pb.body[0]), sizeof(struct control_message)-offset);
	int left = msg->u.pb.total_length - sizeof(struct control_message);
	char *buffer = (char*)cb->buffer + sizeof(struct control_message)- offset;
	iget = internal_read_stream(svc, scd, left, (unsigned char *)buffer);
	if(iget < 0)
	{
		svc->return_data_buffer(scd->sd->cm, cb);
		svc->connection_close(scd->conn);
		return -1;
	}
	svc->trace_out(scd->sd->cm, "CMIB received piggyback msg of length %d, added to read_buffer", 
		       scd->read_buffer_len);

	scd->read_buffer_len = msg->u.pb.total_length - offset;
	scd->read_buffer = cb;
	return 0;
   
}


void
CMIB_data_available(transport_entry trans, CMConnection conn)
{    
	int iget;
	ib_client_data_ptr sd = (ib_client_data_ptr) trans->trans_data;
	CMtrans_services svc = sd->svc;
	struct control_message msg;
	ib_conn_data_ptr scd;
    
	da_t = getlocaltime();
    
	scd = (ib_conn_data_ptr) svc->get_transport_data(conn);

	iget = internal_read_stream(svc, scd, sizeof(struct control_message),
	                            (unsigned char*)&msg);
	if(iget <= 0)
	{
	    svc->trace_out(scd->sd->cm, "CMIB iget is %d, dying\n", 
			   iget);

	    svc->connection_close(conn);
	    return;
	}

	svc->trace_out(scd->sd->cm, "CMIB data available type = %s(%d)", 
		       msg_string[msg.type], msg.type);

	switch(msg.type) {
	case msg_piggyback:
		iget = handle_piggyback(svc, scd, &msg);
		if(iget == 0)
		{
		    CMbuffer tmp = scd->read_buffer;
		    trans->data_available(trans, conn);
		    svc->return_data_buffer(trans->cm, tmp);
		}
		break;
	case msg_response:
		handle_response(svc, scd, &msg);
		break;
	case msg_request:
		handle_request(svc, scd, &msg);
		if(scd->isDone.done == 0)
		{
			trans->data_available(trans, conn);
		}
		else
		{
			svc->trace_out(scd->sd->cm, "CMib data available error in the protocol");
		}
		break;
	default:
		printf("Bad message type %d\n", msg.type);
	}

	//returning control to CM
	svc->trace_out(scd->sd->cm, "CMIB data_available returning");
}

/* 
 * Accept socket connection
 */
static void
ib_accept_conn(void *void_trans, void *void_conn_sock)
{
	transport_entry trans = (transport_entry) void_trans;
	int conn_sock = (int) (long) void_conn_sock;
	ib_client_data_ptr sd = (ib_client_data_ptr) trans->trans_data;
	CMtrans_services svc = sd->svc;
	ib_conn_data_ptr ib_conn_data;
	int sock;
	struct sockaddr sock_addr;
	unsigned int sock_len = sizeof(sock_addr);
	int int_port_num;
	struct linger linger_val;
	int sock_opt_val = 1;

	int delay_value = 1;
	CMConnection conn;
	attr_list conn_attr_list = NULL;
	struct ibparam param, remote_param;

	//ib stuff
	struct ibv_qp_init_attr  qp_init_attr;
    


	svc->trace_out(sd->cm, "Trying to accept something, socket %d\n", conn_sock);
	linger_val.l_onoff = 1;
	linger_val.l_linger = 60;
	if ((sock = accept(conn_sock, (struct sockaddr *) 0, (unsigned int *) 0)) == SOCKET_ERROR) {
		perror("Cannot accept socket connection");
		svc->fd_remove_select(sd->cm, conn_sock);
		fprintf(stderr, "failure in Cmib  removing socket connection\n");
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
	setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	           sizeof(delay_value));
	ib_conn_data = create_ib_conn_data(svc);
	ib_conn_data->sd = sd;
	ib_conn_data->fd = sock;

	//initialize the dataqp that will be used for all RC comms
	memset(&qp_init_attr, 0, sizeof(struct ibv_qp_init_attr));
	qp_init_attr.qp_context = sd->context;
	qp_init_attr.send_cq = sd->send_cq;
	qp_init_attr.recv_cq = sd->recv_cq;
	qp_init_attr.cap.max_recv_wr = LISTSIZE;
	qp_init_attr.cap.max_send_wr = LISTSIZE;
	qp_init_attr.cap.max_send_sge = 32;
	qp_init_attr.cap.max_recv_sge = 1;
	qp_init_attr.cap.max_inline_data = 32;
	qp_init_attr.qp_type = IBV_QPT_RC;
    

	ib_conn_data->dataqp = initqp(ib_conn_data, sd);    
	if(ib_conn_data->dataqp == NULL)
	{
		svc->trace_out(sd->cm, "CMIB can't create qp\n");
		return;
    
	}


	conn_attr_list = create_attr_list();
	conn = svc->connection_create(trans, ib_conn_data, conn_attr_list);
	ib_conn_data->conn = conn;

	add_attr(conn_attr_list, CM_FD, Attr_Int4,
	         (attr_value) (long)sock);

	sock_len = sizeof(sock_addr);
	memset(&sock_addr, 0, sock_len);
	getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
	int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
	add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	         (attr_value) (long)int_port_num);

	memset(&sock_addr, 0, sizeof(sock_addr));
	sock_len = sizeof(sock_addr);
	if (getpeername(sock, &sock_addr, &sock_len) == 0) {
		int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
		add_attr(conn_attr_list, CM_PEER_CONN_PORT, Attr_Int4,
		         (attr_value) (long)int_port_num);
		ib_conn_data->remote_IP = ntohl(((struct sockaddr_in *) &sock_addr)->sin_addr.s_addr);
		add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
		         (attr_value) (long)ib_conn_data->remote_IP);
		if (sock_addr.sa_family == AF_INET) {
			struct hostent *host;
			struct sockaddr_in *in_sock = (struct sockaddr_in *) &sock_addr;
			host = gethostbyaddr((char *) &in_sock->sin_addr,
			                     sizeof(struct in_addr), AF_INET);
			if (host != NULL) {
				ib_conn_data->remote_host = strdup(host->h_name);
				add_attr(conn_attr_list, CM_PEER_HOSTNAME, Attr_String,
				         (attr_value) strdup(host->h_name));
			}
		}
	}
	if (ib_conn_data->remote_host != NULL) {
		svc->trace_out(sd->cm, "Accepted CMIB socket connection from host \"%s\"",
		               ib_conn_data->remote_host);
	} else {
		svc->trace_out(sd->cm, "Accepted CMIB socket connection from UNKNOWN host");
	}

	//here we read the incoming remote contact port number. 
	//in IB we'll extend this to include ib connection parameters
	param.lid  = sd->lid;
	param.qpn  = ib_conn_data->dataqp->qp_num;
	param.port = sd->port;
	param.psn  = sd->psn;

    
	if (read(sock, (char *) &ib_conn_data->remote_contact_port, 4) != 4) {
		svc->trace_out(sd->cm, "Remote host dropped connection without data");
		return;
	}
	if (read(sock, (char *) &remote_param, sizeof(remote_param)) != sizeof(remote_param)) {
		svc->trace_out(sd->cm, "CMIB Remote host dropped connection without data");
		return;
	}
	if(write(sock, &param, sizeof(param)) != sizeof(param))
	{
		svc->trace_out(sd->cm, "CMIB remote side failed to send its parameters");
		return; 
	}
    
	if(connectqp(ib_conn_data, sd, param, remote_param))
	{
		svc->trace_out(sd->cm, "CMIB connectqp failed in accept connection");
		return; 
	}
    
    
	ib_conn_data->remote_contact_port =
		ntohs(ib_conn_data->remote_contact_port);
	add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	         (attr_value) (long)ib_conn_data->remote_contact_port);
	svc->trace_out(sd->cm, "Remote host (IP %x) is listening at port %d\n",
	               ib_conn_data->remote_IP,
	               ib_conn_data->remote_contact_port);

/* dump_sockinfo("accept ", sock); */
	svc->fd_add_select(sd->cm, sock,
	                   (void (*)(void *, void *)) CMIB_data_available,
	                   (void *) trans, (void *) conn);

	svc->trace_out(sd->cm, "Falling out of accept conn\n");
    free_attr_list(conn_attr_list);
}

extern void
libcmib_LTX_shutdown_conn(CMtrans_services svc, ib_conn_data_ptr scd)
{
	svc->trace_out(scd->sd->cm, "CMIB shutdown_conn, removing select %d\n",
	               scd->fd);
	svc->fd_remove_select(scd->sd->cm, scd->fd);
	close(scd->fd);
	//free(scd->remote_host);
	//free(scd->read_buffer);
	free(scd);
}


#include "qual_hostname.c"

static int
is_private_192(int IP)
{
	return ((IP & 0xffff0000) == 0xC0A80000);   /* equal 192.168.x.x */
}

static int
is_private_182(int IP)
{
	return ((IP & 0xffff0000) == 0xB6100000);   /* equal 182.16.x.x */
}

static int
is_private_10(int IP)
{
	return ((IP & 0xff000000) == 0x0A000000);   /* equal 10.x.x.x */
}

static int
initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, ib_conn_data_ptr ib_conn_data, attr_list conn_attr_list, int no_more_redirect)
{
	int sock;

	int delay_value = 1;
	struct linger linger_val;
	int sock_opt_val = 1;
	int int_port_num;
	u_short port_num;
	ib_client_data_ptr sd = (ib_client_data_ptr) trans->trans_data;
	char *host_name;
	int remote_IP = -1;
	static int host_ip = 0;
	unsigned int sock_len;
	union {
	    struct sockaddr s;
	    struct sockaddr_in s_I4;
	    struct sockaddr_in6 s_l6;
	} sock_addr;
	struct ibparam param, remote_param;

	//ib stuff

	int retval = 0;
    
    

	if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_name)) {
		svc->trace_out(cm, "CMIB transport found no IP_HOST attribute");
		host_name = NULL;
	} else {
		svc->trace_out(cm, "CMIB transport connect to host %s", host_name);
	}
	if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_ip)) {
		svc->trace_out(cm, "CMIB transport found no IP_ADDR attribute");
		/* wasn't there */
		host_ip = 0;
	} else {
		svc->trace_out(cm, "CMIB transport connect to host_IP %lx", host_ip);
	}
	if ((host_name == NULL) && (host_ip == 0))
		return -1;

	if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & int_port_num)) {
		svc->trace_out(cm, "CMIB transport found no IP_PORT attribute");
		return -1;
	} else {
		svc->trace_out(cm, "CMIB transport connect to port %d", int_port_num);
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
		char *network_string;
		if ((sock = socket(AF_INET, SOCK_STREAM, 0)) == SOCKET_ERROR) {
			svc->trace_out(cm, " CMIB connect FAILURE --> Couldn't create socket");
			return -1;
		}
		sock_addr.s.sa_family = AF_INET;
		if (((network_string = getenv("CM_NETWORK")) != NULL) &&
		    (host_name != NULL)) {
			int name_len = strlen(host_name) + 2 + strlen(network_string);
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
			if (check_host(new_host_name, (void *) &sock_addr.s_I4.sin_addr) == 0) {
				/* host has no NETWORK interface */
				if (check_host(host_name, (void *) &sock_addr.s_I4.sin_addr) == 0) {
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
				if (check_host(host_name, (void *) &sock_addr.s_I4.sin_addr) == 0) {
					if (host_ip == 0) {
						svc->trace_out(cm, "CMIB connect FAILURE --> Host not found \"%s\", no IP addr supplied in contact list", host_name);
					} else {
						svc->trace_out(cm, "CMIB --> Host not found \"%s\", Using supplied IP addr %x",
						               host_name == NULL ? "(unknown)" : host_name,
						               host_ip);
						((struct sockaddr_in *) &sock_addr)->sin_addr.s_addr = ntohl(host_ip);
					}
				}
			} else {
				sock_addr.s_I4.sin_addr.s_addr = ntohl(host_ip);
			}
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
		svc->trace_out(cm, "Attempting CMIB socket connection, host=\"%s\", IP = %s, port %d",
		               host_name == 0 ? "(unknown)" : host_name, 
		               inet_ntoa(sock_addr.s_I4.sin_addr),
		               int_port_num);
		if (connect(sock, (struct sockaddr *) &sock_addr,
		            sizeof sock_addr) == SOCKET_ERROR) {
#ifdef WSAEWOULDBLOCK
			int err = WSAGetLastError();
			if (err != WSAEWOULDBLOCK || err != WSAEINPROGRESS) {
#endif
				svc->trace_out(cm, "CMIB connect FAILURE --> Connect() to IP %s failed", inet_ntoa(sock_addr.s_I4.sin_addr));
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

	setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	           sizeof(delay_value));

	//initialize the dataqp that will be used for all RC comms

	ib_conn_data->dataqp = initqp(ib_conn_data, sd);
	if(ib_conn_data->dataqp == NULL)
	{
		svc->trace_out(sd->cm, "CMIB initqp failed in initiate_conn\n");
		return -1;  
	}
    
//here we write out the connection port to the other side. 
//for sockets thats all thats required. For IB we can use this to exchange information about the 
//IB parameters for the other side

//What does no_more_redirect check?
	if (!no_more_redirect) {
		int local_listen_port = htons(sd->listen_port);
		if (write(sock, &local_listen_port, 4) != 4) {
			return -1;
		}
    
	}

	param.lid  = sd->lid;
	param.qpn  = ib_conn_data->dataqp->qp_num;
	param.port = sd->port;
	param.psn  = sd->psn;
    
	svc->trace_out(sd->cm, "Writing param, size %d\n", sizeof(param));
	retval = write(sock, &param, sizeof(param));
	if(retval <= 0)
	{
		svc->trace_out(sd->cm, "CMIB write parameter to socket failed %d\n", retval);
		return retval;
    
	}
    
	svc->trace_out(sd->cm, "reading remote param, size %d\n", sizeof(param));
	retval = read(sock, &remote_param, sizeof(param));
	if(retval <= 0)
	{
		svc->trace_out(sd->cm, "CMIB write parameter to socket failed %d\n", retval);
		return retval;
	}    
    

	retval = connectqp(ib_conn_data, sd,
	                   param, remote_param);
	if(retval)
	{
		//svc->trace_out(sd->cm, "CMIB connectqp failed in initiate connection\n");
		return -1;
    
	}
    

	svc->trace_out(cm, "--> Connection established");
	ib_conn_data->remote_host = host_name == NULL ? NULL : strdup(host_name);
	ib_conn_data->remote_IP = remote_IP;
	ib_conn_data->remote_contact_port = int_port_num;
	ib_conn_data->fd = sock;
	ib_conn_data->sd = sd;

	//fixed sized right now, but we will change it to a list/table
	memset(ib_conn_data->infolist, 0, sizeof(rinfo)*10);
	ib_conn_data->infocount = 0;

    

	add_attr(conn_attr_list, CM_FD, Attr_Int4,
	         (attr_value) (long)sock);
	sock_len = sizeof(sock_addr);
	getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
	int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
	add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	         (attr_value) (long)int_port_num);
	add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
	         (attr_value) (long)ib_conn_data->remote_IP);
	if (getpeername(sock, &sock_addr.s, &sock_len) == 0) {
		int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
		add_attr(conn_attr_list, CM_PEER_CONN_PORT, Attr_Int4,
		         (attr_value) (long)int_port_num);
		if (sock_addr.s.sa_family == AF_INET) {
			struct hostent *host;
			struct sockaddr_in *in_sock = (struct sockaddr_in *) &sock_addr;
			host = gethostbyaddr((char *) &in_sock->sin_addr,
			                     sizeof(struct in_addr), AF_INET);
			if (host != NULL) {
				ib_conn_data->remote_host = strdup(host->h_name);
				add_attr(conn_attr_list, CM_PEER_HOSTNAME, Attr_String,
				         (attr_value) strdup(host->h_name));
			}
		}
	}

	svc->trace_out(sd->cm, "Falling out of init conn\n");
	return sock;
}

/* 
 * Initiate a socket connection with another data exchange.  If port_num is -1,
 * establish a unix socket connection (name_str stores the file name of
 * the waiting socket).  Otherwise, establish an INET socket connection
 * (name_str stores the machine name).
 */
extern CMConnection
libcmib_LTX_initiate_conn(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{
	ib_conn_data_ptr ib_conn_data = create_ib_conn_data(svc);
	attr_list conn_attr_list = create_attr_list();
	CMConnection conn;
	int sock;

	if ((sock = initiate_conn(cm, svc, trans, attrs, ib_conn_data, conn_attr_list, 0)) < 0)
		return NULL;

	add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	         (attr_value) (long)ib_conn_data->remote_contact_port);
	conn = svc->connection_create(trans, ib_conn_data, conn_attr_list);
	ib_conn_data->conn = conn;

	svc->trace_out(cm, "Cmib Adding trans->data_available as action on fd %d", sock);
	svc->fd_add_select(cm, sock, (select_list_func) CMIB_data_available,
	                   (void *) trans, (void *) conn);

/* dump_sockinfo("initiate ", sock); */
	return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * For sockets, this involves checking to see if the host name is the 
 * same as ours and if the IP_PORT matches the one we are listening on.
 */
extern int
libcmib_LTX_self_check(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{

	ib_client_data_ptr sd = trans->trans_data;
	int host_addr;
	int int_port_num;
	char *host_name;
	char my_host_name[256];
	static int IP = 0;

	if (IP == 0) {
		IP = get_self_ip_addr(cm, svc);
	}
	if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_name)) {
		svc->trace_out(cm, "CMself check CMIB transport found no IP_HOST attribute");
		host_name = NULL;
	}
	if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_addr)) {
		svc->trace_out(cm, "CMself check CMIB transport found no IP_ADDR attribute");
		if (host_name == NULL) return 0;
		host_addr = 0;
	}
	if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & int_port_num)) {
		svc->trace_out(cm, "CMself check CMIB transport found no IP_PORT attribute");
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
	if (int_port_num != sd->listen_port) {
		svc->trace_out(cm, "CMself check - Ports don't match, %d, %d", int_port_num, sd->listen_port);
		return 0;
	}
	svc->trace_out(cm, "CMself check returning TRUE");
	return 1;
}

extern int
libcmib_LTX_connection_eq(CManager cm, CMtrans_services svc, 
			  transport_entry trans, attr_list attrs,
			  ib_conn_data_ptr scd)
{

	int int_port_num;
	int requested_IP = -1;
	char *host_name = NULL;

	if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_name)) {
		svc->trace_out(cm, "CMIB transport found no IP_HOST attribute");
	}
	if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & int_port_num)) {
		svc->trace_out(cm, "Conn Eq CMIB transport found no IP_PORT attribute");
		return 0;
	}
	if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & requested_IP)) {
		svc->trace_out(cm, "CMIB transport found no IP_ADDR attribute");
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
libcmib_LTX_non_blocking_listen(CManager cm, CMtrans_services svc,
				transport_entry trans, attr_list listen_info)
{
	ib_client_data_ptr sd = trans->trans_data;
	unsigned int length;
	struct sockaddr_in sock_addr;
	int sock_opt_val = 1;
	int conn_sock;
	int attr_port_num = 0;
	u_short port_num = 0;
	char *network_string;

	conn_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (conn_sock == SOCKET_ERROR) {
		fprintf(stderr, "Cannot open INET socket\n");
		return NULL;
	}
	/* 
	 *  Check to see if a bind to a specific port was requested
	 */
	if (listen_info != NULL
	    && !query_attr(listen_info, CM_IP_PORT,
	                   NULL, (attr_value *)(long) & attr_port_num)) {
		port_num = 0;
	} else {
		if (attr_port_num > USHRT_MAX || attr_port_num < 0) {
			fprintf(stderr, "Requested port number %d is invalid\n", attr_port_num);
			return NULL;
		}
		port_num = attr_port_num;
	}

	svc->trace_out(cm, "CMIB begin listen, requested port %d", attr_port_num);
	sock_addr.sin_family = AF_INET;
	sock_addr.sin_addr.s_addr = INADDR_ANY;
	sock_addr.sin_port = htons(port_num);
	if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
	               sizeof(sock_opt_val)) != 0) {
		fprintf(stderr, "Failed to set 1REUSEADDR on INET socket\n");
		return NULL;
	}
	if (sock_addr.sin_port != 0) {
		/* specific port requested */
		svc->trace_out(cm, "CMIB trying to bind selected port %d", port_num);
		if (bind(conn_sock, (struct sockaddr *) &sock_addr,
		         sizeof sock_addr) == SOCKET_ERROR) {
			fprintf(stderr, "Cannot bind INET socket\n");
			return NULL;
		}
	} else {
		/* port num is free.  Constrain to range 26000 : 26100 */
		int low_bound = 26000;
		int high_bound = 26100;
		int size = high_bound - low_bound;
		int tries = 10;
		int result = SOCKET_ERROR;
		srand48(time(NULL)+getpid());
		while (tries > 0) {
			int target = low_bound + size * drand48();
			sock_addr.sin_port = htons(target);
			svc->trace_out(cm, "CMIB trying to bind port %d", target);
			result = bind(conn_sock, (struct sockaddr *) &sock_addr,
			              sizeof sock_addr);
			tries--;
			if (result != SOCKET_ERROR) tries = 0;
		}
		if (result == SOCKET_ERROR) {
			fprintf(stderr, "Cannot bind INET socket\n");
			return NULL;
		}
	}
	sock_opt_val = 1;
	if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
	               sizeof(sock_opt_val)) != 0) {
		perror("Failed to set 2REUSEADDR on INET socket");
		return NULL;
	}
	length = sizeof sock_addr;
	if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
		fprintf(stderr, "Cannot get socket name\n");
		return NULL;
	}
	/* begin listening for conns and set the backlog */
	if (listen(conn_sock, FD_SETSIZE)) {
		fprintf(stderr, "listen failed\n");
		return NULL;
	}
	/* set the port num as one we can be contacted at */

	svc->trace_out(cm, "Cmib Adding ib_accept_conn as action on fd %d", conn_sock);
	svc->fd_add_select(cm, conn_sock, ib_accept_conn,
	                   (void *) trans, (void *) (long)conn_sock);

	/* in the event the DE is shut down, close the socket */
	/* 
	 *  -- Don't do this...  Close() seems to hang on sockets after 
	 *  listen() for some reason.  I haven't found anywhere that defines 
	 *  this behavior, but it seems relatively uniform. 
	 */
	/* DExchange_add_close(de, close_socket_fd, (void*)conn_sock, NULL); */

	{
		char host_name[256];
		int int_port_num = ntohs(sock_addr.sin_port);
		attr_list ret_list;
		int IP = get_self_ip_addr(cm, svc);
		int network_added = 0;

		svc->trace_out(cm, "CMIB listen succeeded on port %d, fd %d",
		               int_port_num, conn_sock);
		ret_list = create_attr_list();
#if !NO_DYNAMIC_LINKING
		get_qual_hostname(cm, host_name, sizeof(host_name) - 1 , svc, listen_info, 
		                  &network_added);
#endif 

		sd->hostname = strdup(host_name);
		sd->listen_port = int_port_num;
		add_attr(ret_list, CM_TRANSPORT, Attr_String,
		         (attr_value) strdup("ib"));
		if ((IP != 0) && (getenv("CM_NETWORK") == NULL) &&
		    (!query_attr(listen_info, CM_NETWORK_POSTFIX, NULL,
		                 (attr_value *) (long)& network_string))) {
			add_attr(ret_list, CM_IP_ADDR, Attr_Int4,
			         (attr_value) (long)IP);
		}
		if ((getenv("CmibUseHostname") != NULL) || 
		    (getenv("CM_NETWORK") != NULL) ||
		    (query_attr(listen_info, CM_NETWORK_POSTFIX, NULL,
		                (attr_value *) (long)& network_string))) {
			add_attr(ret_list, CM_IP_HOSTNAME, Attr_String,
			         (attr_value) strdup(host_name));
			if (network_added) {
				if (query_attr(listen_info, CM_NETWORK_POSTFIX, NULL,
				               (attr_value *) (long)& network_string)) {
					add_attr(ret_list, CM_NETWORK_POSTFIX, Attr_String,
					         (attr_value) strdup(network_string));
				}
			}
		} else if (IP == 0) {
			add_attr(ret_list, CM_IP_ADDR, Attr_Int4, 
			         (attr_value)INADDR_LOOPBACK);
		}
		add_attr(ret_list, CM_IP_PORT, Attr_Int4,
		         (attr_value) (long)int_port_num);

		return ret_list;
	}
}

#if defined(HAVE_WINDOWS_H) && !defined(NEED_IOVEC_DEFINE)
#define NEED_IOVEC_DEFINE
#endif

#ifdef NEED_IOVEC_DEFINE
struct iovec {
	void *iov_base;
	long iov_len;
};

#endif

extern void
libcmib_LTX_set_write_notify(transport_entry trans, CMtrans_services svc,
			     ib_conn_data_ptr scd, int enable)
{
	if (enable != 0) {
		svc->fd_write_select(trans->cm, scd->fd, (select_list_func) trans->write_possible,
		                     (void *)trans, (void *) scd->conn);
	} else {
		/* remove entry */
		svc->fd_write_select(trans->cm, scd->fd, NULL, NULL, NULL);
	}   
}


extern CMbuffer
libcmib_LTX_read_block_func(CMtrans_services svc, ib_conn_data_ptr scd, int *len_ptr, int *offset_ptr)
{
	*len_ptr = scd->read_buffer_len;
	*offset_ptr = 0;
	if (scd->read_buffer) {
	    CMbuffer tmp = scd->read_buffer;
	    scd->read_buffer = NULL;
	    scd->read_buffer_len = 0;
	    return tmp;
	}
	if(scd->tb)
	    return scd->tb->buf;
	return NULL;  
}


#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif

static double reg_t = 0;

static double writev_t = 0;


extern int
libcmib_LTX_writev_complete_notify_func(CMtrans_services svc, 
					ib_conn_data_ptr scd,
					void *iovs,
					int iovcnt,
					attr_list attrs,
					CMcompletion_notify_func notify_func,
					void *notify_client_data)
{
    	int fd = scd->fd;
	int left = 0;
	int iget = 0;
	int i;
	struct iovec * iov = (struct iovec*) iovs;
	struct iovec * tmp_iov;
	int wrlen = 0;
	double start = 0, end = 0;
	struct control_message msg; 
	int can_reuse_mapping = 0;
	rinfo *last_write_request = &scd->infolist[scd->infocount];

	writev_t = getlocaltime();
    
	start = getlocaltime();
    
	for (i = 0; i < iovcnt; i++) {
	    left += iov[i].iov_len;
	}

	svc->trace_out(scd->sd->cm, "CMIB writev of %d bytes on fd %d",
	               left, fd);
	
	if (left < PIGGYBACK)
	{
		//total size is less than the piggyback size
		iget = internal_write_piggyback(svc, scd, left, iov, iovcnt);
		if (notify_func) {
		    (notify_func)(notify_client_data);
		}
		if(iget < 0)
		{
			svc->trace_out(scd->sd->cm, "CMIB error in writing piggyback");
			return -1;
		}
		if(iget == 0)
		{
			return iovcnt;
		}
		return -1;
	}

	svc->set_pending_write(scd->conn);

	if (notify_func) {
	    can_reuse_mapping = 1;
	    /* OK, we're not going to copy the data */
	    if (last_write_request->mrlen == iovcnt) {
		int i;
		for(i=0; i < last_write_request->mrlen; i++) {
		    if ((iov[i].iov_len != last_write_request->wr->sg_list[i].length) ||
			(int64_from_ptr(iov[i].iov_base) != last_write_request->wr->sg_list[i].addr)) {
			can_reuse_mapping = 0;
			svc->trace_out(scd->sd->cm, "CMIB already mapped data, doesn't match write, buf %d, %p vs. %p, %d vs. %d",
				       i, iov[i].iov_base, last_write_request->wr->sg_list[i].addr, iov[i].iov_len, last_write_request->wr->sg_list[i].length);
		    break;
		    }
		}
	    } else {
		svc->trace_out(scd->sd->cm, "CMIB either no already mapped data, or wrong buffer count");
		can_reuse_mapping = 0;
	    }
	} else {
	    svc->trace_out(scd->sd->cm, "CMIB User-owned data with no notify, so no reuse\n");
	}
#ifndef DO_DEREG_ON_FINISH
	if (last_write_request->wr && !(can_reuse_mapping)) {
	    int free_data_elements = (last_write_request->notify_func == NULL);
	    for(i = 0; i < last_write_request->mrlen; i ++)
	    {
		ibv_dereg_mr(last_write_request->mrlist[i]);       
		last_write_request->mrlist[i] = NULL;
	    }
	    
	    free_sg_list(last_write_request->wr->sg_list, last_write_request->mrlen, free_data_elements);
	    free(last_write_request->mrlist);
	    free(last_write_request->wr);
	    last_write_request->wr = NULL;
	}
#endif
	if (notify_func == NULL) {
	    /* 
	     * Semantics are that *MUST* be done with the data when we return,
	     * so, for now, copy all data 
	     */
	    tmp_iov = malloc(sizeof(tmp_iov[0]) * iovcnt);
	    for (i = 0; i < iovcnt; i++) {
		tmp_iov[i].iov_len = iov[i].iov_len;
		tmp_iov[i].iov_base = malloc(iov[i].iov_len);
		memcpy(tmp_iov[i].iov_base, iov[i].iov_base, iov[i].iov_len);
	    }
	} else {
	    /*
	     *  Cool.  The app doesn't need the data back right away.  
	     *  We can keep it and tell the upper levels when we're done.
	     *  Don't copy.  The reply message will trigger the notification.
	     */
	    tmp_iov = iov;
	}
	memset(&msg, 0, sizeof(msg));

    
	end = getlocaltime();
	write_t = end - start;
    
	start = getlocaltime();
	if (!can_reuse_mapping) {
	    scd->infolist[scd->infocount].mrlist = regblocks(scd->sd, tmp_iov, iovcnt,
							     IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_WRITE
							     |IBV_ACCESS_REMOTE_READ,
							     &scd->infolist[scd->infocount].mrlen);
	}
	scd->infolist[scd->infocount].notify_func = notify_func;
	scd->infolist[scd->infocount].notify_client_data = notify_client_data;
	end = getlocaltime();
    
	reg_t = end - start;
    
	if(scd->infolist[scd->infocount].mrlist == NULL)
	{
		svc->trace_out(scd->sd->cm, "CMIB writev error in registereing memory");
		return -1;
	}
	scd->infolist[scd->infocount].wr = createwrlist(scd,
	                                                scd->infolist[scd->infocount].mrlist, tmp_iov,
	                                                scd->infolist[scd->infocount].mrlen, &wrlen);
	
	svc->trace_out(scd->sd->cm, "FIrst 16 bytes of send buffer (len %d) are %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x\n", left, ((unsigned char*)iov[0].iov_base)[0], ((unsigned char*)iov[0].iov_base)[1], ((unsigned char*)iov[0].iov_base)[2], ((unsigned char*)iov[0].iov_base)[3], ((unsigned char*)iov[0].iov_base)[4], ((unsigned char*)iov[0].iov_base)[5], ((unsigned char*)iov[0].iov_base)[6], ((unsigned char*)iov[0].iov_base)[7], ((unsigned char*)iov[0].iov_base)[8], ((unsigned char*)iov[0].iov_base)[9], ((unsigned char*)iov[0].iov_base)[10], ((unsigned char*)iov[0].iov_base)[11], ((unsigned char*)iov[0].iov_base)[12], ((unsigned char*)iov[0].iov_base)[13], ((unsigned char*)iov[0].iov_base)[14], ((unsigned char*)iov[0].iov_base)[15]);
	iget = internal_write_request(svc, scd, left, &scd->infolist[scd->infocount]);
	if(iget < 0)
	{
		svc->trace_out(scd->sd->cm, "CMIB error in writing request");
		return -1;
	}

	if (notify_func == NULL) {
	    /* it was our tmp_iov */
	    free(tmp_iov);
	}

	return iovcnt;
}

extern int
libcmib_LTX_writev_func(CMtrans_services svc, ib_conn_data_ptr scd,
			void *iovs, int iovcnt, attr_list attrs)
{
    return libcmib_LTX_writev_complete_notify_func(svc, scd, iovs, iovcnt, 
						   attrs, NULL, NULL);
}

static int socket_global_init = 0;

static void
free_ib_data(CManager cm, void *sdv)
{
	ib_client_data_ptr sd = (ib_client_data_ptr) sdv;
	CMtrans_services svc = sd->svc;
	if (sd->hostname != NULL)
		svc->free_func(sd->hostname);
	svc->free_func(sd);
}

extern void *
libcmib_LTX_initialize(CManager cm, CMtrans_services svc)
{
	static int atom_init = 0;

	ib_client_data_ptr socket_data;
	svc->trace_out(cm, "Initialize CM IB transport built in %s\n",
	               EVPATH_MODULE_BUILD_DIR);
	page_size = sysconf(_SC_PAGE_SIZE);
	if (socket_global_init == 0) {
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
		CM_PEER_LISTEN_PORT = attr_atom_from_string("PEER_LISTEN_PORT");
		CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
		CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
		atom_init++;
	}
	socket_data = svc->malloc_func(sizeof(struct ib_client_data));
	socket_data->cm = cm;
	socket_data->hostname = NULL;
	socket_data->listen_port = -1;
	socket_data->svc = svc;
	socket_data->ibdev = IB_getdevice(NULL);
	socket_data->context = ibv_open_device(socket_data->ibdev);
	socket_data->port = 1; //need to somehow get proper port here
	socket_data->lid = get_local_lid(socket_data->context, socket_data->port);
	socket_data->pd = ibv_alloc_pd(socket_data->context);
	socket_data->send_channel = ibv_create_comp_channel(socket_data->context);
	socket_data->send_cq = ibv_create_cq(socket_data->context, 1024, 
	                                     (void*)socket_data, socket_data->send_channel, 0);
    

	socket_data->recv_channel = ibv_create_comp_channel(socket_data->context);
	socket_data->recv_cq = ibv_create_cq(socket_data->context, 1024, 
	                                     (void*)socket_data, socket_data->recv_channel, 0);
    

	//create srq
	struct ibv_srq_init_attr sqa;
    
	sqa.attr.max_wr = 64;
	sqa.attr.max_sge = 1;
	sqa.attr.srq_limit = 1;
	sqa.srq_context = (void*)socket_data;
    
    
    
	socket_data->srq = ibv_create_srq(socket_data->pd, 
	                                  &sqa);
	if(socket_data->srq == NULL)
	{
		svc->trace_out(cm, "unable to create srq\n");
	}    
    
	socket_data->psn = lrand48()%256;

	// //set up padding
	// socket_data->pad.mr = ibv_reg_mr(socket_data->pd, socket_data->pad.pad, 
	//                   sizeof(socket_data->pad.pad), IBV_ACCESS_LOCAL_WRITE);
    

	svc->add_shutdown_task(cm, free_ib_data, (void *) socket_data, FREE_TASK);

	//here we will add the first 4MB memory buffer
	LIST_INIT(&memlist);
	LIST_INIT(&uselist);

    
	int bsize = 4*1024*1024;
	void *buffer;
    
	if (posix_memalign(&buffer, page_size, bsize) != 0) {
		svc->trace_out(socket_data->cm, "Unable to register initial memory - this is bad!\n");
		return NULL;    
	}    
	tbuffer *tb = (tbuffer*)malloc(sizeof(tbuffer));
	CMbuffer cb = svc->create_data_and_link_buffer(socket_data->cm, buffer, bsize);
	cb->return_callback = free_func;
	cb->return_callback_data = (void*)tb;
    
	tb->buf = cb;
	tb->inuse = 1;
	tb->scd = NULL;
	tb->size = bsize;
	tb->offset = 0;
	tb->parent = NULL;
	tb->childcount = 0;    
	tb->mr = ibv_reg_mr(socket_data->pd, tb->buf->buffer, bsize,
	                    IBV_ACCESS_LOCAL_WRITE | 
	                    IBV_ACCESS_REMOTE_WRITE | 
	                    IBV_ACCESS_REMOTE_READ);
	if(!tb->mr)
	{
		svc->trace_out(socket_data->cm, "Unable to register initial memory - this is bad!\n");
		return NULL;    
	}    

	LIST_INSERT_HEAD(&memlist, tb, entries);
	perftrace = (getenv("CMIBTransportVerbose") != NULL);    

	return (void *) socket_data;
}

static struct ibv_qp * initqp(ib_conn_data_ptr ib_conn_data,
                              ib_client_data_ptr sd)                    
{
	struct ibv_qp_init_attr qp_init_attr;
	struct ibv_qp_attr qp_attr;
	struct ibv_qp *dataqp;
	int retval = 0;

	memset(&qp_init_attr, 0, sizeof(struct ibv_qp_init_attr));
	qp_init_attr.qp_context = sd->context;
	qp_init_attr.send_cq = sd->send_cq;
	qp_init_attr.recv_cq = sd->recv_cq;
	qp_init_attr.cap.max_recv_wr = LISTSIZE;
	qp_init_attr.cap.max_send_wr = LISTSIZE;
	qp_init_attr.cap.max_send_sge = 32;
	qp_init_attr.cap.max_recv_sge = 1;
	qp_init_attr.cap.max_inline_data = 0;
	qp_init_attr.qp_type = IBV_QPT_RC;
	qp_init_attr.srq = NULL;
//    qp_init_attr.sq_sig_all = 1;

    
 try_again:
	dataqp = ibv_create_qp(sd->pd, &qp_init_attr);
	if(dataqp == NULL) {
	    if (qp_init_attr.cap.max_send_sge > 1) {
	        sd->svc->trace_out(sd->cm, "CMIB ibv_create_qp failed with max_send_sge = %d, trying smaller", qp_init_attr.cap.max_send_sge);
		qp_init_attr.cap.max_send_sge >>= 1;
		goto try_again;
	    } else {
	        perror("ibv_create_qp");
		sd->svc->trace_out(sd->cm, "CMIB can't create qp\n");
		return NULL;
	    }
	}

	memset(&qp_attr, 0, sizeof(qp_attr));
	qp_attr.qp_state = IBV_QPS_INIT;
	qp_attr.pkey_index = 0;
	qp_attr.port_num = sd->port;
	qp_attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE;
	qp_attr.qkey = 0x11111111;

	retval = ibv_modify_qp(dataqp, &qp_attr, 
	                       IBV_QP_STATE |
	                       IBV_QP_PKEY_INDEX | 
	                       IBV_QP_PORT | 
	                       IBV_QP_ACCESS_FLAGS);
	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMIB unable to set qp to INIT %d\n", retval);
		ibv_destroy_qp(dataqp); 
		return NULL;
	}

	//register the notification memory block
	memset(&ib_conn_data->isDone, 0, sizeof(notify));
	ib_conn_data->isDone.mr = ibv_reg_mr(sd->pd, 
	                                     &ib_conn_data->isDone.done, 
	                                     sizeof(ib_conn_data->isDone.done), 
	                                     IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ);

	if(ib_conn_data->isDone.mr == NULL)
	{
		sd->svc->trace_out(sd->cm, "CMib unable to create notification mr\n");
		ibv_destroy_qp(dataqp);
		return NULL;
	}


	//register padding

	ib_conn_data->isDone.sg.addr = int64_from_ptr(&ib_conn_data->isDone.done);
	ib_conn_data->isDone.sg.length = sizeof(ib_conn_data->isDone.done);
	ib_conn_data->isDone.sg.lkey = ib_conn_data->isDone.mr->lkey;
    
	ib_conn_data->isDone.wr.wr_id = int64_from_ptr(ib_conn_data);
	ib_conn_data->isDone.wr.sg_list = &ib_conn_data->isDone.sg;
	ib_conn_data->isDone.wr.num_sge = 1;
	ib_conn_data->isDone.wr.opcode = IBV_WR_SEND;
	ib_conn_data->isDone.wr.send_flags = IBV_SEND_SIGNALED;
	ib_conn_data->isDone.wr.next = NULL;

	ib_conn_data->isDone.rwr.wr_id = int64_from_ptr(ib_conn_data);
	ib_conn_data->isDone.rwr.next = NULL;
	ib_conn_data->isDone.rwr.sg_list = &ib_conn_data->isDone.sg;
	ib_conn_data->isDone.rwr.num_sge = 1;


	//issue some 10 receives on the qp  so we can get 10 notifications 
	int i = 0;
	for(i = 0; i < 10; i++)
	{
    
		retval = ibv_post_recv(dataqp, &ib_conn_data->isDone.rwr, 
		                       &ib_conn_data->isDone.badrwr);
		if(retval)
		{
			sd->svc->trace_out(sd->cm, "CMib unable to post recv %d\n", retval);
			ibv_dereg_mr(ib_conn_data->isDone.mr);
			ibv_destroy_qp(dataqp);
        
			return NULL;        
		}
    
	}
    
    
    
	return dataqp;
}


static int connectqp(ib_conn_data_ptr ib_conn_data,
                     ib_client_data_ptr sd,
                     struct ibparam lparam,
                     struct ibparam rparam)
{
	struct ibv_qp_attr qp_attr;
    
	int retval = 0;
    
	if(ib_conn_data == NULL || ib_conn_data->dataqp == NULL)
		return -1;    

	memset(&qp_attr, 0, sizeof(struct ibv_qp_attr));

	qp_attr.qp_state = IBV_QPS_RTR;
	qp_attr.dest_qp_num = rparam.qpn;
	qp_attr.rq_psn = rparam.psn;
	qp_attr.sq_psn = lparam.psn;
	qp_attr.ah_attr.is_global = 0;
	qp_attr.ah_attr.dlid = rparam.lid;
	qp_attr.ah_attr.sl = 0;
	qp_attr.ah_attr.src_path_bits = 0;
	qp_attr.ah_attr.port_num = rparam.port; 
	qp_attr.path_mtu = IBV_MTU_1024;
	qp_attr.max_dest_rd_atomic = 16;
	qp_attr.min_rnr_timer = 24;
	qp_attr.timeout = 28;
	qp_attr.retry_cnt = 18;
	qp_attr.rnr_retry = 18;
	qp_attr.max_rd_atomic = 4;
    
    
	retval = ibv_modify_qp(ib_conn_data->dataqp, &qp_attr,
	                       IBV_QP_STATE |
	                       IBV_QP_AV |
	                       IBV_QP_PATH_MTU |
	                       IBV_QP_DEST_QPN |
	                       IBV_QP_RQ_PSN |  
	                       IBV_QP_MIN_RNR_TIMER | IBV_QP_MAX_DEST_RD_ATOMIC);
	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMIB unable to set qp to RTR %d\n", retval);
		return retval;  

	}
    

	//   qp_attr.cap.max_inline_data = 1;    
	// qp_attr.cap.max_send_wr = 1024;
	// qp_attr.cap.max_recv_wr = 1024;
	// qp_attr.cap.max_send_sge = 32;
	// qp_attr.cap.max_recv_sge = 1;
    
	qp_attr.qp_state = IBV_QPS_RTS;
	retval = ibv_modify_qp(ib_conn_data->dataqp, &qp_attr, IBV_QP_STATE|
	                       IBV_QP_TIMEOUT| 
	                       IBV_QP_RETRY_CNT|
	                       IBV_QP_RNR_RETRY|
	                       IBV_QP_SQ_PSN| IBV_QP_MAX_QP_RD_ATOMIC |
	                       IBV_QP_MAX_QP_RD_ATOMIC);

	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMIB unable to set qp to RTS %d\n", retval);
		return retval;  

	}

	retval = ibv_req_notify_cq(sd->send_cq, 0);
	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMib notification request failed\n");

		//cleaqnup
		return -1;  
	}


	retval = ibv_req_notify_cq(sd->recv_cq, 0);
	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMib notification request failed\n");

		//cleaqnup
		return -1;  
	}

	return 0;    
}


static struct ibv_mr ** regblocks(ib_client_data_ptr sd,
                                  struct iovec *iovs, int iovcnt, int flags, 
                                  int *mrlen)                 
{
	int i =0;
    
	struct ibv_mr **mrlist;

    
	mrlist = (struct ibv_mr**) malloc(sizeof(struct ibv_mr *) * iovcnt);
	if(mrlist == NULL)
	{
		//failed to allocate memory - big issue
		return NULL;    
	}
    
	for(i = 0; i < iovcnt; i++)
	{
    
		mrlist[i] = ibv_reg_mr(sd->pd, iovs[i].iov_base, 
		                       iovs[i].iov_len, 
		                       flags);
		if(mrlist[i] == NULL)
		{
			fprintf(stderr, "registeration failed \n");
			for(; i > 0; i--)
			{
				ibv_dereg_mr(mrlist[i-1]);      
			}
			free(mrlist);
			return NULL;        
		}   
	}
	*mrlen = iovcnt;
    
	return mrlist;    
}



static struct ibv_send_wr * createwrlist(ib_conn_data_ptr conn, 
                                         struct ibv_mr **mrlist,
                                         struct iovec *iovlist,
                                         int mrlen, int *wrlen)
{
	//create an array of work requests that can be posted for the transter
	struct ibv_qp_attr attr;
	struct ibv_qp_init_attr init_attr;
	struct ibv_sge *sge;
	struct ibv_send_wr *wr;
	int retval = 0;
	ib_client_data_ptr sd = conn->sd;
    
    
    
	memset(&attr, 0, sizeof(attr));
	memset(&init_attr, 0, sizeof(init_attr));
    
	//query to get qp params
	retval = ibv_query_qp(conn->dataqp, &attr, IBV_QP_CAP, &init_attr);
	if(retval)
	{
		sd->svc->trace_out(sd->cm, "CMIB unable to query initial state %d\n", retval);
		return NULL;    
	}
    
	conn->max_imm_data = attr.cap.max_inline_data;
    
    
	// fprintf(stderr, "wr = %d\tsge = %d\timm = %d %d\n",
	//      attr.cap.max_send_wr, attr.cap.max_send_sge, attr.cap.max_inline_data, 
	//      init_attr.cap.max_inline_data);
	// fprintf(stderr, "mrlen = %d\n", mrlen);
    

	if(mrlen > attr.cap.max_send_sge)
	{
		fprintf(stderr, "too many sge fall back to slow mode\n");
		//do the slow mode here
		//TODO still
	}
	else
		*wrlen = 1;
    

	sge = (struct ibv_sge*)malloc(sizeof(struct ibv_sge) * (mrlen));
	if(sge == NULL)
	{
		fprintf(stderr, "couldn't allocate memory\n");
		return NULL;
    
	}
    
	wr=(struct ibv_send_wr*)malloc(sizeof(struct ibv_send_wr)*(*wrlen));
	if(wr == NULL)
	{
		fprintf(stderr, "malloc failed for wr\n");
		free(sge);
		return NULL;    
	}
    
    
	wr->wr_id = int64_from_ptr(conn);
	wr->next = NULL;    
	wr->sg_list = sge;    
	wr->num_sge = mrlen;    
	wr->opcode = IBV_WR_RDMA_WRITE;    
	wr->send_flags = IBV_SEND_FENCE| IBV_SEND_SIGNALED ;    
	wr->imm_data = 0;    

	int i = 0;
    
	int len = 0;
	for(i = 0; i <mrlen; i++)
	{
		len += iovlist[i].iov_len;
	}
	for(i = 0; i <mrlen; i++) {
	    sge[i].addr = int64_from_ptr(iovlist[i].iov_base);
	    sge[i].length = iovlist[i].iov_len;
	    sge[i].lkey = mrlist[i]->lkey;  
	}
	// sge[mrlen].addr = int64_from_ptr(sd->pad.pad);
	// sge[mrlen].length = sizeof(sd->pad.pad);
	// sge[mrlen].lkey = sd->pad.mr->lkey;
    
    
	return wr;
}

static int waitoncq(ib_conn_data_ptr scd,
                    ib_client_data_ptr sd,
                    CMtrans_services svc, struct ibv_cq *cq)
{


	struct ibv_wc wc;    
	int retval  = 0;
	struct ibv_cq *ev_cq;

	memset(&wc, 0, sizeof(wc));    

	retval = ibv_req_notify_cq(cq, 0);
	if(retval)
	{
		svc->trace_out(scd->sd->cm, "CMib notification request failed\n");

		//cleaqnup
		return -1;  
	}
    
    
	retval = ibv_get_cq_event(cq->channel,
	                          &ev_cq, (void*)scd);
	if(retval)
	{
		svc->trace_out(sd->cm, "Failed to get event\n");
		//cleanup
		return -1;  
	}
    
	//ack the event
	ibv_ack_cq_events(cq, 1);

	//reequest notify on cq 

	retval = ibv_req_notify_cq(ev_cq, 0);
	if(retval)
	{
		svc->trace_out(sd->cm, "CMib notification request failed\n");   
		//cleanup
		return -1;  
	}
    
	return 0;    
}


static tbuffer *findMemory(ib_conn_data_ptr scd, ib_client_data_ptr sd, 
                           CMtrans_services svc, int req_size)
{
	tbuffer *temp = NULL, *prov= NULL;
    
	for(temp = memlist.lh_first;temp != NULL; temp = temp->entries.le_next)
	{
	    if(((temp->size - temp->offset) >= req_size) && (temp->inuse ==0))
		{
			//possible match
			if(!prov || (prov->size - prov->offset) >= (temp->size - temp->offset))
			{
				prov = temp;
			}
		}   
	}
	if(!prov)
	{
		//couldn't find matching memory
		//allocate new buffer
		tbuffer *tb = (tbuffer*)malloc(sizeof(tbuffer));
		void *buffer;
		if (posix_memalign(&buffer, page_size, req_size) != 0) {
			svc->trace_out(sd->cm, "Unable to register initial memory - this is bad!\n");
			return NULL;    
		}
    
		CMbuffer cb = svc->create_data_and_link_buffer(sd->cm, buffer, req_size);
		tb->buf = cb;
		tb->scd = scd;
		tb->size = req_size;
		tb->inuse = 1;
		tb->offset = 0;
		cb->return_callback = free_func;
		cb->return_callback_data = (void*)tb;
		tb->childcount = 0;
		tb->parent = NULL;  

		tb->mr = ibv_reg_mr(sd->pd, buffer, req_size,
		                    IBV_ACCESS_LOCAL_WRITE | 
		                    IBV_ACCESS_REMOTE_WRITE | 
		                    IBV_ACCESS_REMOTE_READ);
		if(!tb->mr)
		{
			svc->trace_out(sd->cm, "Unable to register initial memory - this is bad!\n");
			return NULL;    
		}   
		LIST_INSERT_HEAD(&memlist, tb, entries);    
		return tb;              
	}
	else
	{
		//found matching memory but we can't just use this block still 
		//because FFS will blow up if we don't make new CMbuffer
		//on the other hand we don't have to register atleast!
		tbuffer *tb = (tbuffer*)malloc(sizeof(tbuffer));
		//tb = new tbuffer, prov = old tbuffer
		tb->parent = prov;
		prov->childcount ++;
		// fprintf(stderr, "Original buffer = %p req_size = %d offset = %d\n",
		//  prov->buf, req_size, prov->offset);

		void *buffer = ptr_from_int64((int64_from_ptr(prov->buf->buffer) + prov->offset));
    
    
		uint64_t oldsize = prov->size;
    
		uint64_t ptr = int64_from_ptr(buffer);
		buffer = ptr_from_int64(((ptr+(8-1))& ~(8-1)));
		uint64_t newsize = int64_from_ptr(buffer) - int64_from_ptr(prov->buf->buffer);
		uint64_t bsize = prov->size - newsize;  
		// fprintf(stderr, "Resizing %p tbuff to %d from %d\n", prov, prov->size, 
		//  newsize);

		prov->size = newsize;
		tb->size = bsize;
		if(oldsize != (newsize + bsize))
		{
			fprintf(stderr, "lost some memory here\n");     
		}

		tb->buf = svc->create_data_and_link_buffer(sd->cm, buffer, tb->size);
		tb->buf->return_callback = free_func;
		tb->buf->return_callback_data = tb;
		tb->inuse = 1;
		tb->offset = req_size;

		tb->scd = scd;
		tb->childcount = 0;
		tb->mr = prov->mr;


		LIST_INSERT_HEAD(&memlist, tb, entries);	




		return tb;  
	}
    
    
}


extern transport_entry
cmib_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport = svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("ib");
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)libcmib_LTX_initialize;
    transport->listen = (CMTransport_listen_func)libcmib_LTX_non_blocking_listen;
    transport->initiate_conn = (CMConnection(*)())libcmib_LTX_initiate_conn;
    transport->self_check = (int(*)())libcmib_LTX_self_check;
    transport->connection_eq = (int(*)())libcmib_LTX_connection_eq;
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)libcmib_LTX_shutdown_conn;
    transport->read_block_func = (CMTransport_read_block_func)libcmib_LTX_read_block_func;
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)NULL;
    transport->writev_func = (CMTransport_writev_func)libcmib_LTX_writev_func;
    transport->writev_complete_notify_func = (CMTransport_writev_complete_notify_func)libcmib_LTX_writev_complete_notify_func;
    transport->get_transport_characteristics = NULL;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}

#endif

