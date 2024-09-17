/***** Includes *****/
#include "config.h"
#include <sys/types.h>

#include <inttypes.h>
#include <getopt.h>


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
#include <rdma/fabric.h>
#include <rdma/fi_endpoint.h>
#include <rdma/fi_rma.h>
#include <rdma/fi_cm.h>
#include <rdma/fi_errno.h>
#include <rdma/fi_eq.h>

#include <atl.h>
#include "evpath.h"
#include "cm_transport.h"
#include "cm_internal.h"
#include "ev_select.h"

#include <stdlib.h>

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif

#define _WITH_IB_
#define PIGGYBACK 1025*10

#ifdef _WITH_IB_


#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 2259)
#  pragma warning (disable: 177)
#endif


//   BEGIN from shared.h in fabtests

/* haven't tested with earlier than version 1.2 */
#define FT_FIVERSION FI_VERSION(1,2)

void cq_readerr(struct fid_cq *cq, char *cq_str);


#define FT_PRINTERR(call, retv) \
	do { fprintf(stderr, call "(): %d, %d (%s)\n", __LINE__, (int) retv, fi_strerror((int) -retv)); } while (0)

#define FT_ERR(fmt, ...) \
	do { fprintf(stderr, "%s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); } while (0)

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

//   END from shared.h in fabtests


//all the message types have a queue associated with them
static char *msg_string[] = {"Request", "Response", "Piggyback"};
enum {msg_request = 0, msg_response = 1, msg_piggyback = 2} msg_type;

struct remote_entry {
    uint64_t remote_addr;
    uint32_t length;
    uint64_t rkey;
};

struct request
{
    uint64_t length;
    uint32_t iovcnt;
    uint32_t piggyback_length;
    uint64_t request_ID;
    struct remote_entry  read_list[1];
};


struct response
{
    uint64_t max_length;    
    uint64_t request_ID;
};

struct piggyback
{
    uint32_t total_length;
    uint32_t padding;
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

#define ptr_from_int64(p) (void *)(unsigned long)(p)
#define int64_from_ptr(p) (u_int64_t)(unsigned long)(p)


struct msg_item;

typedef enum pull_status_t {PULL_UNSUBMITTED = 0, PULL_SUBMITTED, PULL_COMPLETED, PULL_ACCOUNTED} pull_status_t;

struct pull_record {
    struct remote_entry remote;
    char *dest;
    pull_status_t status;
    struct msg_item *parent;
};

struct mr_list_item;
struct fabric_connection_data;

typedef struct msg_item {
    struct fabric_connection_data *conn_data;
    char *buffer;
    long cur_buffer_size;
    long message_size;
    uint64_t request_ID;
    struct fid_mr *mr;
    int pull_count;
    struct pull_record *pulls;
    struct mr_list_item *mr_rec;
    struct msg_item *next;
} *msg_item_p;

struct fabric_client_data;

typedef struct mr_list_item {
    struct fabric_client_data *fabd;
    long buffer_length;
    void *buffer;
    struct fid_mr *mr;
    int in_use;
    struct mr_list_item *next;
} *mr_list;

typedef struct fabric_client_data {
    CManager cm;
    CMtrans_services svc;
    transport_entry trans;
    struct fi_info *hints;

    struct fid_fabric *fab;
    struct fid_pep *listen_ep;
    struct fid_domain *dom;
    struct fid_eq *cmeq;

    char *hostname;
    int listen_port;
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

    struct timeval pull_schedule_base;
    struct timeval pull_schedule_period;
    CMavail_period_ptr avail;

    mr_list existing_mr_list;
    
    int thread_init;
    msg_item_p pull_queue;
    msg_item_p completed_queue;
    pthread_mutex_t pull_queue_mutex;
    int thread_should_run;
    struct fid_wait *send_waitset;
    fd_set readset;
    int nfds;
    int wake_read_fd;
    int wake_write_fd;
    struct fabric_connection_data **fcd_array_by_sfd;
    int ncqs;
    struct fid **cq_array;
} *fabric_client_data_ptr;


typedef struct remote_info
{
    int iovcnt;
    struct fid_mr **mrlist;
    struct iovec *iov;
    CMcompletion_notify_func notify_func;
    void *notify_client_data;
}rinfo;

typedef struct fabric_connection_data {
    fabric_client_data_ptr fabd;
    struct fid_cq *rcq, *scq;
    struct fid_mr *read_mr;
    struct fid_mr *send_mr;
    struct fid_ep *conn_ep;
    size_t buffer_size;
    void *mapped_recv_buf;
    char *send_buf;
    CMbuffer read_buf;
    int max_credits;
    int read_buffer_len;
    int read_offset;

    char *remote_host;
    int remote_IP;
    int remote_contact_port;
    int fd;
    int sfd;   /* fd for the scq */
    CMConnection conn;
    struct ibv_qp *dataqp;    
    int infocount;
    rinfo last_write;
    int max_imm_data;   
} *fabric_conn_data_ptr;

static inline int msg_offset()
{
	return ((int) (((char *) (&(((struct control_message *)NULL)->u.pb.body[0]))) - ((char *) NULL)));
}

static int alloc_cm_res(fabric_client_data_ptr fabd);
static int alloc_ep_res(fabric_conn_data_ptr fcd, struct fi_info *fi);
static int bind_ep_res(fabric_conn_data_ptr fcd);
static void free_ep_res(fabric_conn_data_ptr fcd);
static void
hand_to_pull_thread(CMtrans_services svc, fabric_client_data_ptr fabd,
		    msg_item_p msg);


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
static atom_t CM_IP_INTERFACE = -1;
static atom_t CM_TRANSPORT = -1;

static int
check_host(hostname, sin_addr)
	char *hostname;
void *sin_addr;
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

static fabric_conn_data_ptr 
create_fabric_conn_data(CMtrans_services svc)
{
    fabric_conn_data_ptr fabric_conn_data = svc->malloc_func(sizeof(struct fabric_connection_data));
    memset(fabric_conn_data, 0, sizeof(struct fabric_connection_data));
    fabric_conn_data->remote_host = NULL;
    fabric_conn_data->remote_contact_port = -1;
    fabric_conn_data->fd = 0;
    fabric_conn_data->read_buf = NULL;
    fabric_conn_data->read_buffer_len = 0;
    
    return fabric_conn_data;
}


static void
handle_scq_completion(struct fi_cq_data_entry *comp)
{
    if (comp->flags & FI_READ) {
	struct pull_record *this_pull = 
	    (struct pull_record *)comp->op_context;
//	printf("Got an FI_READ completion\n");
	this_pull->status = PULL_COMPLETED;
    }
    if (comp->flags & FI_SEND) {
	*(int*)comp->op_context = 1;
    }
}

static int internal_write_piggyback(CMtrans_services svc,
                                    fabric_conn_data_ptr fcd,
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
    svc->trace_out(fcd->fabd->cm, "CMFABRIC sending piggyback msg of length %d,", 
		   length);

	
    for (i = 0; i < iovcnt; i++)
    {
	memcpy(point, iov[i].iov_base, iov[i].iov_len);
	point += iov[i].iov_len;
    }
    
    {
	
	memcpy(fcd->send_buf, msg, msg->u.pb.total_length);
	int ret, sent = 0;
	
	svc->trace_out(fcd->fabd->cm, "fi_send on conn_ep, length %d, send_buf %p\n", msg->u.pb.total_length, fcd->send_buf);
	ret = fi_send(fcd->conn_ep, fcd->send_buf, msg->u.pb.total_length, fi_mr_desc(fcd->send_mr), 0, &sent);
	if (ret) {
	    FT_PRINTERR("fi_send", ret);
	    return ret;
	}
	
	/* Read send queue */
	do {
	    struct fi_cq_data_entry comp;
	    ret = fi_cq_read(fcd->scq, &comp, 1);
	    if (ret < 0 && ret != -FI_EAGAIN) {
		FT_PRINTERR("fi_cq_read", ret);
		cq_readerr(fcd->scq, " in internal write piggyback");
		return ret;
	    }
	    if (ret == 1) handle_scq_completion(&comp);
	} while (!sent);
    }
    
    free(msg);
    return 0;
}


static int internal_write_response(CMtrans_services svc,
                                   fabric_conn_data_ptr fcd,
                                   int length,
				   int64_t request_ID)
{
    struct control_message msg;

    msg.type = msg_response;
    msg.u.resp.max_length = length;
    msg.u.resp.request_ID = request_ID;

    memcpy(fcd->send_buf, &msg, sizeof(msg));
    int ret, sent = 0;
	
    svc->trace_out(fcd->fabd->cm, "fi_send for write response\n");
    ret = fi_send(fcd->conn_ep, fcd->send_buf, sizeof(msg), fi_mr_desc(fcd->send_mr), 0, &sent);
    if (ret) {
	FT_PRINTERR("fi_send", ret);
	return ret;
    }
	
    /* Read send queue */
    do {
	struct fi_cq_data_entry comp;
	svc->trace_out(fcd->fabd->cm, "fi_cq_read for send completion in write response\n");

	ret = fi_cq_read(fcd->scq, &comp, 1);
	if (ret < 0 && ret != -FI_EAGAIN) {
	    FT_PRINTERR("fi_cq_read", ret);
	    cq_readerr(fcd->scq, " in internal write response");
	    return ret;
	}
	if (ret == 1) handle_scq_completion(&comp);
    } while (!sent);

    return 0;
}

static int internal_write_request(CMtrans_services svc,
                                  fabric_conn_data_ptr fcd,
                                  int length,
    				  rinfo *request_info)
{
	struct control_message *msg;
	int ret, i, piggyback_prefix_vectors;
	int size = sizeof(*msg) + request_info->iovcnt * sizeof(struct remote_entry);
	msg = malloc(size);

	msg->type = msg_request;
	msg->u.req.length = length;
	msg->u.req.iovcnt = request_info->iovcnt;
	msg->u.req.request_ID = int64_from_ptr(request_info);
	msg->u.req.piggyback_length = 0;

	piggyback_prefix_vectors = 0;
	for (i=0; i < request_info->iovcnt; i++) {
	    if (size + request_info->iov[i].iov_len < PIGGYBACK) {
		size += request_info->iov[i].iov_len;
		piggyback_prefix_vectors++;
	    } else {
		break;
	    }
	}
	msg->u.req.iovcnt -= piggyback_prefix_vectors;
	/* redo size */
	size = sizeof(*msg) + (request_info->iovcnt-piggyback_prefix_vectors) 
	    * sizeof(struct remote_entry);

	for (i=0; i < piggyback_prefix_vectors; i++) {
	    msg = realloc(msg, size + request_info->iov[i].iov_len);
	    memcpy(((char*)msg) + size, request_info->iov[i].iov_base, 
		   request_info->iov[i].iov_len);
	    size += request_info->iov[i].iov_len;
	    msg->u.req.piggyback_length += request_info->iov[i].iov_len;
	}
	/* handle remaining */
	for (; i < request_info->iovcnt; i++) {
	    int j = i - piggyback_prefix_vectors;
	    msg->u.req.read_list[j].remote_addr = int64_from_ptr(request_info->iov[i].iov_base);
	    msg->u.req.read_list[j].length = request_info->iov[i].iov_len;
	    msg->u.req.read_list[j].rkey = fi_mr_key(request_info->mrlist[i]);
	    svc->trace_out(fcd->fabd->cm, "Adding source buffer[%d] %lx, len %d, key %lx\n",
		   i, msg->u.req.read_list[j].remote_addr, msg->u.req.read_list[j].length,
		   msg->u.req.read_list[j].rkey);
	}

	svc->trace_out(fcd->fabd->cm, "Doing internal write request, writing %d bytes", size);

	memcpy(fcd->send_buf, msg, size);
	int sent = 0;
	ret = fi_send(fcd->conn_ep, fcd->send_buf, size, fi_mr_desc(fcd->send_mr), 0, &sent);
	if (ret) {
	    FT_PRINTERR("fi_send", ret);
	    return ret;
	}
	
	/* Read send queue */
	do {
	    struct fi_cq_data_entry comp;
	    ret = fi_cq_read(fcd->scq, &comp, 1);
	    if (ret < 0 && ret != -FI_EAGAIN) {
		FT_PRINTERR("fi_cq_read", ret);
		cq_readerr(fcd->scq, " in internal write request");
		return ret;
	    }
	    if (ret == 1) handle_scq_completion(&comp);
	} while (!sent);
	return 0;
}
	
static int handle_response(CMtrans_services svc,
                           fabric_conn_data_ptr fcd,
                           struct control_message *msg)
{

    //read back response
    struct response *rep;
    rinfo *write_request;
    
    rep = &msg->u.resp;
    write_request = ptr_from_int64(rep->request_ID);
    
    if (rep->max_length == -1) {
	fprintf(stderr, "WRITE FAILED, request %p\n", write_request);
	return -1;
    }

    if (write_request->notify_func) {
	(write_request->notify_func)( write_request->notify_client_data);
    }
    svc->wake_any_pending_write(fcd->conn);
    return 0;   
    
}

static mr_list
get_free_mr(fabric_client_data_ptr fabd, long length)
{
    mr_list list, return_item = NULL;
    pthread_mutex_lock(&fabd->pull_queue_mutex);
    list = fabd->existing_mr_list;
    while (list != NULL) {
	if ((!list->in_use) && (list->buffer_length >= length)) {
	    list->in_use = 1;
	    return_item = list;
	    break;
	}
	list = list->next;
    }
    pthread_mutex_unlock(&fabd->pull_queue_mutex);
    return return_item;
}

static void
add_mr_to_list(fabric_client_data_ptr fabd, mr_list mr_rec)
{
    pthread_mutex_lock(&fabd->pull_queue_mutex);
    mr_rec->fabd = fabd;
    mr_rec->next = fabd->existing_mr_list;
    fabd->existing_mr_list = mr_rec;
    pthread_mutex_unlock(&fabd->pull_queue_mutex);
}

static void 
free_func(void *mr_item_v)
{
    fabric_client_data_ptr fabd;
    mr_list list, mr_item = mr_item_v;
    fabd = mr_item->fabd;
    pthread_mutex_lock(&fabd->pull_queue_mutex);
    list = fabd->existing_mr_list;
    while (list != NULL) {
	if (list ==  mr_item_v) {
	    list->in_use = 0;
	    break;
	}
	list = list->next;
    }
    if (list == NULL) {
	fprintf(stderr, "libfabric MR list inconsistency in free_func\n");
    }
    pthread_mutex_unlock(&fabd->pull_queue_mutex);
}
	

static void
add_to_pull_queue(CMtrans_services svc, fabric_conn_data_ptr fcd, 
	     struct control_message *msg)
{
    int i;
    struct request *req = &msg->u.req;
    int header_size = sizeof(*msg) + req->iovcnt * sizeof(struct remote_entry);


    msg_item_p msg_pull_item = calloc(1, sizeof(struct msg_item));
    msg_pull_item->message_size = req->length;
    msg_pull_item->buffer = malloc(req->piggyback_length);
    msg_pull_item->cur_buffer_size = req->piggyback_length;
    /* copy portion of msg that was piggybacked */
    memcpy(msg_pull_item->buffer, (char*)msg + header_size, req->piggyback_length);
    msg_pull_item->pull_count = req->iovcnt;
    msg_pull_item->pulls = calloc(req->iovcnt, sizeof(struct pull_record));
    msg_pull_item->conn_data = fcd;
    msg_pull_item->request_ID = req->request_ID;
    for (i = 0; i < msg_pull_item->pull_count; i++) {
	msg_pull_item->pulls[i].parent = msg_pull_item;
	msg_pull_item->pulls[i].remote = req->read_list[i];
    }
    hand_to_pull_thread(svc, fcd->fabd, msg_pull_item);
}

static int
perform_pull(CMtrans_services svc, fabric_conn_data_ptr fcd, 
	     struct control_message *msg)
{
    int ret;
    int count;
    int i;
    char *ptr;
    int header_size;
    struct request *req;
    void *buffer;
    struct fid_mr *mr;
    struct mr_list_item *mr_rec;
    CMbuffer cb;
    req = &msg->u.req;
    
    svc->trace_out(fcd->fabd->cm, "In handle request, len is %ld, iovcnt %d, request_ID %lx, piggyback_length = %d\n",
	   req->length, req->iovcnt, req->request_ID, req->piggyback_length);
    mr_rec = get_free_mr(fcd->fabd, req->length);
    if (mr_rec == NULL) {
	mr_rec = calloc(1, sizeof(*mr_rec));
	buffer = malloc(req->length);
	mr_rec->buffer_length = req->length;
	mr_rec->buffer = buffer;
	mr_rec->in_use = 1;
	mr_rec->next = NULL;
	svc->trace_out(fcd->fabd->cm, "fi_mr_reg, buff %p, size %ld with attrs REMOTE_READ,REMOTE_WRITE, SEND, RECV\n", buffer, req->length);
	ret = fi_mr_reg(fcd->fabd->dom, buffer, req->length,
			FI_REMOTE_READ | FI_REMOTE_WRITE | FI_SEND | FI_RECV,
			0, 0, 0, &mr, NULL);
	if(ret) {
	    FT_PRINTERR("fi_mr_reg", ret);
	    svc->trace_out(fcd->fabd->cm, "Failed to get memory\n");
	    internal_write_response(svc, fcd, -1, req->request_ID);
	    return ret;
	}

	mr_rec->mr = mr;
	add_mr_to_list(fcd->fabd, mr_rec);
    } else {
	svc->trace_out(fcd->fabd->cm, "Reusing exiting MR buff %p, size %ld\n", mr_rec->buffer, mr_rec->buffer_length);
	buffer = mr_rec->buffer;
	mr = mr_rec->mr;
	mr_rec->in_use = 1;
    }
    cb = svc->create_data_and_link_buffer(fcd->fabd->cm, buffer, req->length);
    cb->return_callback = free_func;
    cb->return_callback_data = (void*)mr_rec;
    header_size = sizeof(*msg) + req->iovcnt * sizeof(struct remote_entry);
    memcpy(cb->buffer, (char*)msg + header_size, req->piggyback_length);

    ptr = cb->buffer + req->piggyback_length;
    count = req->iovcnt;
    for (i=0; i < req->iovcnt; i++) {
    	struct fi_cq_data_entry comp;
	struct pull_record this_pull;
	this_pull.status = PULL_SUBMITTED;
	int call_count = 0;
	svc->trace_out(fcd->fabd->cm, "fi_read, buffer %p, len %d, (keys.rkey %lx, keys.addr %lx)\n", ptr, req->read_list[i].length, req->read_list[i].rkey, req->read_list[i].remote_addr);
	ssize_t ret = fi_read(fcd->conn_ep, ptr, req->read_list[i].length, fi_mr_desc(mr),
			      0, req->read_list[i].remote_addr, req->read_list[i].rkey, (void*)&this_pull);
	ptr += req->read_list[i].length;
	do {
	    /* Read send queue */
	    ret = fi_cq_read(fcd->scq, &comp, 1);
	    if (ret < 0 && ret != -FI_EAGAIN) {
		FT_PRINTERR("fi_cq_read", ret);
		cq_readerr(fcd->scq, "scq");
	    }
	    call_count++;
	    if (ret == 1) handle_scq_completion(&comp);
	} while (this_pull.status != PULL_COMPLETED);

    	if (comp.flags & FI_READ) {
    	    count--;
	}
    }


    while (count != 0) {
    	struct fi_cq_data_entry comp;
	ret = fi_cq_read(fcd->scq, &comp, 1);
	if (ret == -FI_EAGAIN) {printf("Eagain\n"); continue;}
	if (ret < 0 && ret != -FI_EAGAIN) {
	    if (ret == -FI_EAVAIL) {
		cq_readerr(fcd->scq, "scq");
	    } else {
		FT_PRINTERR("fi_cq_read", ret);
	    }
	    continue;
	} else if (ret > 0) {
	    printf("Successful remote read, Completion size is %ld, data %p\n", comp.len, comp.buf);
	}
	if (ret == 1) handle_scq_completion(&comp);

    	if (comp.flags & FI_READ) {
    	    count--;
    	}
    }

    fcd->read_buf = cb;
    svc->trace_out(fcd->fabd->cm, "FIrst 16 bytes of receive buffer (len %d) are %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x \n", req->length, ((unsigned char*)cb->buffer)[0], ((unsigned char*)cb->buffer)[1], ((unsigned char*)cb->buffer)[2], ((unsigned char*)cb->buffer)[3], ((unsigned char*)cb->buffer)[4], ((unsigned char*)cb->buffer)[5], ((unsigned char*)cb->buffer)[6], ((unsigned char*)cb->buffer)[7], ((unsigned char*)cb->buffer)[8], ((unsigned char*)cb->buffer)[9], ((unsigned char*)cb->buffer)[10], ((unsigned char*)cb->buffer)[11], ((unsigned char*)cb->buffer)[12], ((unsigned char*)cb->buffer)[13], ((unsigned char*)cb->buffer)[14], ((unsigned char*)cb->buffer)[15]);
    fcd->read_buffer_len = req->length;
    svc->trace_out(fcd->fabd->cm, "CMFABRIC handle_request completed");    
    internal_write_response(svc, fcd, req->length, req->request_ID);
    return 0;
}


static void
wake_pull_thread(fabric_client_data_ptr fabd) 
{
    static char buffer = 'W';  /* doesn't matter what we write */
    if (!fabd->send_waitset) {
	if (write(fabd->wake_write_fd, &buffer, 1) != 1) {
	    if (fabd->thread_should_run)
		printf("Whoops, wake_pull_thread write failed\n");
	}
    } else {
	/* send a cq wake */
    }
}

static void
hand_to_pull_thread(CMtrans_services svc, fabric_client_data_ptr fabd,
		    msg_item_p msg)
{
    struct msg_item *queue = fabd->pull_queue;
    msg->next = NULL;
    pthread_mutex_lock(&fabd->pull_queue_mutex);
    if (fabd->pull_queue == NULL) {
	fabd->pull_queue = msg;
    } else {
	while (queue->next) queue = queue->next;
	queue->next = msg;
    }
    pthread_mutex_unlock(&fabd->pull_queue_mutex);
    /* wake the thread if he's sleeping */
    wake_pull_thread(fabd);
}

static void
return_completed_pull(CMtrans_services svc, fabric_client_data_ptr fabd,
		      msg_item_p msg)
{
    pthread_mutex_lock(&fabd->pull_queue_mutex);
    /* find msg in pull queue and take it out */
    if (fabd->pull_queue == msg) {
	fabd->pull_queue = msg->next;
	msg->next = NULL;
    } else {
	struct msg_item *queue = fabd->pull_queue;
	while (queue->next != msg) queue = queue->next;
	queue->next = msg->next;
	msg->next = NULL;
    }

    /* find add msg to completion queue */
    if (fabd->completed_queue == NULL) {
	fabd->completed_queue = msg;
    } else {
	struct msg_item *queue = fabd->completed_queue;
	while (queue->next) queue = queue->next;
	queue->next = msg;
    }
    svc->trace_out(fabd->cm, "Returning completed msg to the CM network thread\n");
    pthread_mutex_unlock(&fabd->pull_queue_mutex);
    /* 
     *  this should wake the CM server thread and we'll check the 
     *  pull completion queue in our handler 
     */
    fabd->svc->wake_comm_thread(fabd->cm);

}

static void
kill_thread(fabric_client_data_ptr fabd)
{
    fabd->thread_should_run = 0;
    if (fabd->thread_init) wake_pull_thread(fabd);
    fabd->thread_init = 0;
}

static void
update_period_base(fabric_client_data_ptr fabd, struct timeval *now)
{
    
    struct timeval next;
    /* add period to base until it's greater than now.  
     * Return base before the last add so that base is just less than now.
     */
    timeradd(&fabd->pull_schedule_base, &fabd->pull_schedule_period, 
	     &next);
    while(timercmp(&next, now, <)) {
	fabd->pull_schedule_base = next;
	timeradd(&fabd->pull_schedule_base, &fabd->pull_schedule_period, 
		 &next);
    }
}
	
static int
in_pull_period(fabric_client_data_ptr fabd, struct timeval *now)
{
    int i = 0;
    int period_count = 0;
    struct timeval cur_offset, zero = {0,0};
    timersub(now, &fabd->pull_schedule_base, &cur_offset);
    update_period_base(fabd, now);
    while (timercmp(&fabd->avail[period_count].duration, &zero, !=))  {
	period_count++;
    }
    for (i = 0; i < period_count; i++) {
	struct timeval end;
	timeradd(&fabd->avail[i].offset, &fabd->avail[i].duration, &end);
	if (timercmp(&cur_offset, &end, <) && 
	    timercmp(&cur_offset, &fabd->avail[i].offset, >)) return 1;
    }
    return 0;
}

static struct timeval
delay_until_wake(fabric_client_data_ptr fabd, struct timeval *now)
{
    /* calculate the delay until the next pull period opens */
    int i = 0;
    int period_count = 0;
    struct timeval cur_offset, zero = {0,0};
    timersub(now, &fabd->pull_schedule_base, &cur_offset);
    update_period_base(fabd, now);
    while (timercmp(&fabd->avail[period_count].duration, &zero, !=))  {
	period_count++;
    }
    /* the one we want is the next start time past now */
    for (i = 0; i < period_count; i++) {
	if (timercmp(&fabd->avail[i].offset, &cur_offset, >)) {
	    struct timeval ret;
	    timersub(&fabd->avail[i].offset, &cur_offset, &ret);
	    return ret;
	}
    }
    /* we must have been past the last offset, 
       return the first offset of the next period */
    struct timeval ret;
    timeradd(&fabd->pull_schedule_period, &fabd->avail[0].offset, &ret);
    timersub(&ret, &cur_offset, &ret);
    return ret;
}

static void
add_completion_fd(fabric_client_data_ptr fabd, fabric_conn_data_ptr fcd)
{
    if (fcd->sfd > fabd->nfds) {
	fabd->fcd_array_by_sfd = realloc(fabd->fcd_array_by_sfd,
					 sizeof(fabd->fcd_array_by_sfd[0]) * fcd->sfd+1);
	fabd->nfds = fcd->sfd;
    }
    fabd->fcd_array_by_sfd[fcd->sfd] = fcd;
    fabd->cq_array = realloc(fabd->cq_array, 
			     sizeof(fabd->cq_array[0]) * (fabd->ncqs+1));
    fabd->cq_array[fabd->ncqs] = &fcd->scq->fid;
    fabd->ncqs++;
    FD_SET(fcd->sfd, &fabd->readset);
}


static void
remove_completion_fd(fabric_client_data_ptr fabd, fabric_conn_data_ptr fcd)
{
    int i;
    for (i = 0; i < fabd->ncqs; i++) {
	if (fabd->cq_array[i] == &fcd->scq->fid) {
	    memcpy(&fabd->cq_array[i+1], &fabd->cq_array[i], 
		   sizeof(fabd->cq_array[0]) * (fabd->ncqs - i - 1));
	}
    }
    fabd->ncqs--;
    FD_CLR(fcd->sfd, &fabd->readset);
}

#include <poll.h>

static int
can_do_something(fabric_client_data_ptr fabd)
{
    CMtrans_services svc = fabd->svc;
    int did_something = 0;
    int i;
    msg_item_p pull = fabd->pull_queue;
    while (pull) {
	msg_item_p next = pull->next;
	fabric_conn_data_ptr fcd = pull->conn_data;
	svc->trace_out(fabd->cm, "Got a pull message with buffer count %d, next %p\n",
		       pull->pull_count, pull->next);
	if (!pull->mr_rec) {
	    int ret;
	    /* try to get a mr */
	    mr_list mr_rec = get_free_mr(fcd->fabd, pull->message_size);
	    char *buffer = NULL;
	    if (mr_rec == NULL) {
		if (0 /*!allow_register_new_mr(fcd, pull->message_size)*/) {
		    pull = pull->next;
		    continue;
		}
		mr_rec = calloc(1, sizeof(*mr_rec));
		buffer = malloc(pull->message_size);
		mr_rec->buffer_length = pull->message_size;
		mr_rec->buffer = buffer;
		mr_rec->in_use = 1;
		mr_rec->next = NULL;
		svc->trace_out(fcd->fabd->cm, "fi_mr_reg, buff %p, size %ld with attrs REMOTE_READ,REMOTE_WRITE, SEND, RECV\n", buffer, pull->message_size);
		ret = fi_mr_reg(fcd->fabd->dom, buffer, pull->message_size,
				FI_REMOTE_READ | FI_REMOTE_WRITE | FI_SEND | FI_RECV,
				0, 0, 0, &mr_rec->mr, NULL);
		if(ret) {
		    FT_PRINTERR("fi_mr_reg", ret);
		    svc->trace_out(fcd->fabd->cm, "Failed to get memory\n");
		    internal_write_response(svc, fcd, -1, pull->request_ID);
		    return ret;
		}
		pull->mr = mr_rec->mr;
		pull->mr_rec = mr_rec;
		add_mr_to_list(fcd->fabd, mr_rec);
	    } else {
		svc->trace_out(fcd->fabd->cm, "Reusing exiting MR buff %p, size %ld\n", mr_rec->buffer, mr_rec->buffer_length);
		mr_rec->in_use = 1;
		pull->mr_rec = mr_rec;
		buffer = mr_rec->buffer;
		pull->mr = mr_rec->mr;
	    }
	    /* piggyback part of msg is in old pull->buffer */
	    memcpy(buffer, pull->buffer, pull->cur_buffer_size);
	    free(pull->buffer);
	    pull->buffer = buffer;

	    /* setup local buffers, now that we have memory */
	    char *ptr = pull->buffer + pull->cur_buffer_size;
	    pull->cur_buffer_size = pull->message_size;
	    for ( i = 0; i < pull->pull_count; i++) {
		pull->pulls[i].dest = ptr;
		ptr += pull->pulls[i].remote.length;
	    }
	    pull->cur_buffer_size = mr_rec->buffer_length;
	    did_something = 1;
	}
	/* at this point we should have an MR and mr_rec */
	for ( i = 0; i < pull->pull_count; i++) {
	    if (pull->pulls[i].status == PULL_UNSUBMITTED) {

//		struct fi_cq_data_entry comp;
		struct remote_entry *remote = &pull->pulls[i].remote;
		struct pull_record *this_pull  = &pull->pulls[i];
		if (0 /* can't do more pulls for some reason */) {
		    continue;
		}
		
		svc->trace_out(fcd->fabd->cm, "fi_read, buffer %p, len %d, (keys.rkey %lx, keys.addr %lx)\n", this_pull->dest, remote->length, remote->rkey, remote->remote_addr);
//		struct pollfd poll_list[1];
//		poll_list[0].fd = fcd->sfd;
//		poll_list[0].events = POLLIN|POLLPRI;
//		int pret = poll(&poll_list[0], (unsigned long)1, 0);
//		printf("Before read poll list ret %d, revents = %x\n", pret, poll_list[0].revents);
		ssize_t ret = fi_read(fcd->conn_ep, this_pull->dest, 
				      remote->length, fi_mr_desc(pull->mr),
				      0, remote->remote_addr, remote->rkey, 
				      (void*)this_pull);
		if (ret) continue;
		this_pull->status = PULL_SUBMITTED;
		did_something = 1;
		
		add_completion_fd(fabd, fcd);
/* 		do { */
/* 		    /\* Read send queue *\/ */
/* 		    ret = fi_cq_read(fcd->scq, &comp, 1); */
/* 		    if (ret < 0 && ret != -FI_EAGAIN) { */
/* 			FT_PRINTERR("fi_cq_read", ret); */
/* 			cq_readerr(fcd->scq, "scq"); */
/* 		    } */
/* //		    pret = poll(&poll_list[0], (unsigned long)1, 0); */
/* //		    printf("AFTER read poll list ret %d, revents = %x\n", pret, poll_list[0].revents); */
/* 		} while (ret == -FI_EAGAIN); */
/* 		svc->trace_out(fcd->fabd->cm, "FI_READ Completion op_context is %p, flags %lx, len %ld, buf %p, data %lx\n", */
/* 		       comp.op_context, comp.flags, comp.len, comp.buf, comp.data); */
/* 		if (ret == 1) { */
/* 		    printf("Immediate completion\n"); */
/* 		    handle_scq_completion(&comp); */
/* 		} */
	    }
	}
	/* if we get to this point, all pulls must have been submitted. */
	if (pull->mr_rec) {
	    int all_done = 1;
	    for ( i = 0; i < pull->pull_count; i++) {
		if (pull->pulls[i].status == PULL_COMPLETED) {
		    pull->pulls[i].status = PULL_ACCOUNTED;
		}
		if (pull->pulls[i].status != PULL_ACCOUNTED) {
		    all_done = 0;
		}
	    }
	    if (all_done) {
		remove_completion_fd(fabd, fcd);
		return_completed_pull(svc, fcd->fabd, pull);
	    }
	}
	pull = next;
    }
    if (!did_something) svc->trace_out(fabd->cm, "PULL_THREAD Nothing to do\n");
    return did_something;
}

static void
handle_completions(fabric_client_data_ptr fabd, fd_set *readset)
{
    /* everything left should be scq completions */
    int i;
    if (readset) {
	for (i=0; i < fabd->nfds; i++) {
	    if (FD_ISSET(i, readset)) {
		int ret;
		struct fi_cq_data_entry comp;
		/* Read send queue */
		fabric_conn_data_ptr fcd = fabd->fcd_array_by_sfd[i];
		ret = fi_cq_read(fcd->scq, &comp, 1);
		if (ret < 0 && ret != -FI_EAGAIN) {
		    FT_PRINTERR("fi_cq_read", ret);
		    cq_readerr(fcd->scq, "scq");
		}
		if (ret == 1) {
		    handle_scq_completion(&comp);
		}
	    }
	}
    } else {
	for (i=0; i < fabd->ncqs; i++) {
	    int ret;
	    struct fi_cq_data_entry comp;
	    /* Read send queue */
	    ret = fi_cq_read((struct fid_cq *) fabd->cq_array[i], &comp, 1);
	    if (ret < 0 && ret != -FI_EAGAIN) {
		FT_PRINTERR("fi_cq_read", ret);
		cq_readerr((struct fid_cq *)fabd->cq_array[i], "scq");
	    }
	    if (ret == 1) {
		handle_scq_completion(&comp);
	    }
	}
    }
}
		

static void
pull_thread(fabric_client_data_ptr fabd)
{
//    thread_setup();  /* extract wait fd, other? */
//    CMtrans_services svc = fabd->svc;
    while(fabd->thread_should_run) {
	struct timeval now, delay;
	fd_set readset;
	gettimeofday(&now, NULL);
	if (in_pull_period(fabd, &now)) {
	    while (can_do_something(fabd));
        }
	/* calculate time to next wake */
	/* sleep until something happens, *or* the next pull period begins */ 
	delay = delay_until_wake(fabd, &now);
	readset = fabd->readset;
	if (fabd->ncqs == 0) {
	    /* we're just waiting on wake fd */
	    select(fabd->nfds+1, &readset, NULL, NULL, &delay);
	} else if (fi_trywait(fabd->fab, &fabd->cq_array[0], fabd->ncqs) == FI_SUCCESS) {
	    select(fabd->nfds+1, &readset, NULL, NULL, &delay);
	    if (FD_ISSET(fabd->wake_read_fd, &readset)) {
		/* read and discard wake byte */
		char buffer;
		if (read(fabd->wake_read_fd, &buffer, 1) != 1) {
		    perror("wake read failed\n");
		}
		FD_CLR(fabd->wake_read_fd, &readset);
	    }
	    handle_completions(fabd, &readset);
	} else {
	    /* test all CQs */
	    handle_completions(fabd, NULL);
	}
    }
//    thread_free_resources;
}

/*
 * handle a pull request message 
 */
static void handle_request(CMtrans_services svc,
                           fabric_conn_data_ptr fcd,
                           struct control_message *msg)
{
    //handling the request message
    if (fcd->fabd->avail) {
	add_to_pull_queue(svc, fcd, msg);
	wake_pull_thread(fcd->fabd);
    } else {
	/* immediate pull and handle */
	(void) perform_pull(svc, fcd, msg);
    }
}

void cq_readerr(struct fid_cq *cq, char *cq_str)
{ 
	struct fi_cq_err_entry cq_err;
	const char *err_str;
	int ret;

	ret = fi_cq_readerr(cq, &cq_err, 0);
	if (ret < 0)
		FT_PRINTERR("fi_cq_readerr", ret);

	err_str = fi_cq_strerror(cq, cq_err.prov_errno, cq_err.err_data, NULL, 0);
	fprintf(stderr, "%s %s (%d)\n", cq_str, err_str, cq_err.prov_errno);
}

void
CMFABRIC_data_available(transport_entry trans, CMConnection conn)
{    
	fabric_client_data_ptr fabd = (fabric_client_data_ptr) trans->trans_data;
	CMtrans_services svc = fabd->svc;
	struct control_message *msg;
	fabric_conn_data_ptr fcd;
	int ret, call_data_available;
	CMbuffer CMbuffer_to_return = NULL;
	struct fid **fids = malloc(sizeof(fids[0]));
	fcd = (fabric_conn_data_ptr) svc->get_transport_data(conn);
	fids[0] = &fcd->rcq->fid;
	fabd->trans = trans;

	svc->trace_out(fabd->cm, "At the beginning of CMFabric_data_available: ");
	ret = fi_trywait(fcd->fabd->fab, fids, 1);
	switch (ret) {
	case FI_SUCCESS:
	    svc->trace_out(fabd->cm, "Try wait on rcq returned FI_SUCCESS");
	    break;
	case -FI_EAGAIN:
	    svc->trace_out(fabd->cm, "Try wait on rcq returned FI_EAGAIN, read cq events");
	    break;
	default:
	    svc->trace_out(fabd->cm, "Try wait on rcq returned %d", ret);
	}
	{
	    struct fi_cq_data_entry comp;
		ret = fi_cq_read(fcd->rcq, &comp, 1);
		if (ret == -FI_EAGAIN) {
		    return;
		}
		if (ret < 0 && ret != -FI_EAGAIN) {
			if (ret == -FI_EAVAIL) {
				cq_readerr(fcd->rcq, "rcq");
			} else {
				FT_PRINTERR("fi_cq_read", ret);
				return;
			}
			return;
		} else if (ret > 0) {
		    //
		}
		svc->trace_out(fabd->cm, "FI_RECV Completion op_context is %p, flags %lx, len %ld, buf %p, data %lx\n",
		       comp.op_context, comp.flags, comp.len, comp.buf, comp.data);
	}

	fcd = (fabric_conn_data_ptr) svc->get_transport_data(conn);
	msg = (struct control_message *) fcd->mapped_recv_buf;
	svc->trace_out(fcd->fabd->cm, "CMFABRIC data available type = %s(%d)", 
		       msg_string[msg->type], msg->type);

	call_data_available = 0;
	CMbuffer_to_return = NULL;
	switch(msg->type) {
	case msg_piggyback: {
	    	int offset = msg_offset();
		fcd->read_buffer_len = msg->u.pb.total_length - offset;
		svc->trace_out(fcd->fabd->cm, "CMFABRIC received piggyback msg of length %d, added to read_buffer", 
			       fcd->read_buffer_len);
		
		fcd->read_buf = fcd->fabd->svc->get_data_buffer(trans->cm, fcd->read_buffer_len);
		memcpy(fcd->read_buf->buffer, &msg->u.pb.body[0], fcd->read_buffer_len);
		fcd->read_offset = 0;

		CMbuffer_to_return = fcd->read_buf;
		call_data_available = 1;
		break;
	}
	case msg_response:
		handle_response(svc, fcd, msg);
		break;
	case msg_request:
		handle_request(svc, fcd, msg);
		call_data_available = 1;
		CMbuffer_to_return = fcd->read_buf;
		break;
	default:
		printf("Bad message type %d\n", msg->type);
	}
	/* post the next receive before relinquishing control */
	ret = fi_recv(fcd->conn_ep, fcd->mapped_recv_buf, fcd->buffer_size, fi_mr_desc(fcd->read_mr), 0, fcd->mapped_recv_buf);
	if (ret)
		FT_PRINTERR("fi_recv", ret);

	if (call_data_available) {
	    trans->data_available(trans, conn);
	}
	if (CMbuffer_to_return) {
	    svc->return_data_buffer(trans->cm, CMbuffer_to_return);
	}
	svc->trace_out(fcd->fabd->cm, "CMFABRIC data_available returning");
}

static int server_connect(fabric_conn_data_ptr fcd);

/* 
 * Accept socket connection
 */
static void
fabric_accept_conn(void *void_trans, void *void_conn_sock)
{
    transport_entry trans = (transport_entry) void_trans;
    fabric_client_data_ptr fabd = (fabric_client_data_ptr) trans->trans_data;
    CMtrans_services svc = fabd->svc;
    fabric_conn_data_ptr fcd;
    int fd, ret;
    struct sockaddr sock_addr;
    unsigned int sock_len = sizeof(sock_addr);

    CMConnection conn;
    attr_list conn_attr_list = NULL;

    //ib stuff
    fcd = create_fabric_conn_data(svc);
    fcd->fabd = fabd;

    server_connect(fcd);
    //initialize the dataqp that will be used for all RC comms

    conn_attr_list = create_attr_list();
    conn = svc->connection_create(trans, fcd, conn_attr_list);
    fcd->conn = conn;

    sock_len = sizeof(sock_addr);
    memset(&sock_addr, 0, sock_len);
//    getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
//    int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
//    add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
//	     (attr_value) (long)int_port_num);

    memset(&sock_addr, 0, sizeof(sock_addr));
    sock_len = sizeof(sock_addr);
    /* if (getpeername(sock, &sock_addr, &sock_len) == 0) { */
    /* 	int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port); */
    /* 	add_attr(conn_attr_list, CM_PEER_CONN_PORT, Attr_Int4, */
    /* 		 (attr_value) (long)int_port_num); */
    /* 	fcd->remote_IP = ntohl(((struct sockaddr_in *) &sock_addr)->sin_addr.s_addr); */
    /* 	add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4, */
    /* 		 (attr_value) (long)fcd->remote_IP); */
    /* 	if (sock_addr.sa_family == AF_INET) { */
    /* 	    struct hostent *host; */
    /* 	    struct sockaddr_in *in_sock = (struct sockaddr_in *) &sock_addr; */
    /* 	    host = gethostbyaddr((char *) &in_sock->sin_addr, */
    /* 				 sizeof(struct in_addr), AF_INET); */
    /* 	    if (host != NULL) { */
    /* 		fcd->remote_host = strdup(host->h_name); */
    /* 		add_attr(conn_attr_list, CM_PEER_HOSTNAME, Attr_String, */
    /* 			 (attr_value) strdup(host->h_name)); */
    /* 	    } */
    /* 	} */
    /* } */
    if (fcd->remote_host != NULL) {
	svc->trace_out(fabd->cm, "Accepted CMFABRIC socket connection from host \"%s\"",
		       fcd->remote_host);
    } else {
	svc->trace_out(fabd->cm, "Accepted CMFABRIC socket connection from UNKNOWN host");
    }

    
    if ((ret = fi_control (&fcd->rcq->fid, FI_GETWAIT, (void *) &fd))) {
	FT_PRINTERR("fi_control(FI_GETWAIT)", ret);
    }
    add_attr(conn_attr_list, CM_FD, Attr_Int4,
	     (attr_value) (long)fd);

    svc->trace_out(fabd->cm, "Cmfabric Adding trans->data_available as action on fd %d", fd);
    svc->fd_add_select(fabd->cm, fd, (select_list_func) CMFABRIC_data_available,
		       (void *) trans, (void *) conn);

    svc->trace_out(fabd->cm, "Falling out of accept conn\n");
    free_attr_list(conn_attr_list);
    fcd->fd = fd;
    if ((ret = fi_control (&fcd->scq->fid, FI_GETWAIT, (void *) &fcd->sfd))) {
	FT_PRINTERR("fi_control(FI_GETWAIT)", ret);
    }
}

/* 
 * incoming event on CM eq
 */
static void
fabric_service_incoming(void *void_trans, void *void_eq)
{
    transport_entry trans = (transport_entry) void_trans;
    fabric_client_data_ptr fabd = (fabric_client_data_ptr) trans->trans_data;
    struct fi_eq_cm_entry entry;
    uint32_t event;
    ssize_t rd;

    rd = fi_eq_sread(fabd->cmeq, &event, &entry, sizeof entry, -1, FI_PEEK);
    if (rd != sizeof entry) {
	if (rd == -FI_EAVAIL) {
	    struct fi_eq_err_entry error = {0};
	    int rc = fi_eq_readerr(fabd->cmeq, &error, 0);
	    if (rc) {
		char buf[1024];
		fprintf(stderr, "error event: %s\n", fi_eq_strerror(fabd->cmeq, error.prov_errno,
      error.err_data, buf, 1024));
	    }
	} else {
	    FT_PRINTERR("fi_eq_sread", rd);
	}
	return;
    }
    
    if (event == FI_CONNREQ) {
	fabric_accept_conn(void_trans, void_eq);
    } else {
	rd = fi_eq_sread(fabd->cmeq, &event, &entry, sizeof entry, -1, 0);
	if (event == FI_SHUTDOWN){
	    fabd->svc->trace_out(fabd->cm, "CMFABRIC got a shutdown event for some conn, who knows which one?\n");
	} else {
	    printf("Unexpected event in service incoming,%s %d\n", fi_tostr(&event, FI_TYPE_EQ_EVENT), event);
	}
    }
}

extern void
libcmfabric_LTX_shutdown_conn(svc, fcd)
	CMtrans_services svc;
fabric_conn_data_ptr fcd;
{
	svc->trace_out(fcd->fabd->cm, "CMFABRIC shutdown_conn, removing select %d\n",
	               fcd->fd);
	svc->fd_remove_select(fcd->fabd->cm, fcd->fd);
	close(fcd->fd);
	//free(fcd->remote_host);
	//free(fcd->read_buffer);
	if (fcd->last_write.iov) free(fcd->last_write.iov);
	free(fcd);
}


static int client_connect(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs, fabric_conn_data_ptr fcd)
{
    fabric_client_data_ptr fabd = fcd->fabd;
    struct fi_eq_cm_entry entry;
    uint32_t event;
    struct fi_info *fi;
    ssize_t rd;
    int ret, int_port_num;
    struct in_addr dest_ip;
    char *host_name, *host_rep = NULL;
    char *dst_port = NULL;
    int i;

    /* Get fabric info */
    if (!get_int_attr(attrs, CM_IP_ADDR,(int*) & dest_ip.s_addr)) {
	svc->trace_out(cm, "CMFABRIC transport found no IP_ADDR attribute");
    } else {
	host_rep = malloc(16);
	dest_ip.s_addr = htonl(dest_ip.s_addr);
	sprintf(host_rep, "%s", inet_ntoa(dest_ip));

    }
    if (!get_int_attr(attrs, CM_IP_PORT, (int*) & int_port_num)) {
	svc->trace_out(cm, "CMFABRIC transport found no IP_PORT attribute");
    } else {
	dst_port = malloc(10);
	sprintf(dst_port, "%d", int_port_num);
    }
    svc->trace_out(fabd->cm, "Connecting to addr, %s, port %s\n", host_rep, dst_port);
    if (!get_string_attr(attrs, CM_IP_HOSTNAME, &host_name)) {
	svc->trace_out(cm, "CMFABRIC transport found no IP_HOSTNAME attribute");
    } else {
      host_rep = malloc(strlen(host_name));
      for (i = 0; i < (strlen(host_name)/2); i++) {
	sscanf(&host_name[i*2], "%2hhx", &host_rep[i]);
      }
      /* printf("name len is %d\n", (int)strlen(host_name)/2); */
      /* for(i = 0; i < strlen(host_name)/2; i++) { */
      /* 	printf("%02x", (unsigned char) host_rep[i]); */
      /* } */
      /* printf(" done\n"); */
    }
    ret = fi_getinfo(FT_FIVERSION, host_rep, dst_port, 0, fabd->hints, &fi);
    svc->trace_out(cm, "%s return value fi is %s\n", "client", fi_tostr(fi, FI_TYPE_INFO));
    if (ret) {
	FT_PRINTERR("fi_getinfo", ret);
	goto err0;
    }

    /* Open fabric */
    ret = fi_fabric(fi->fabric_attr, &fabd->fab, NULL);
    if (ret) {
	FT_PRINTERR("fi_fabric", ret);
	goto err1;
    }

    /* Open domain */
    ret = fi_domain(fabd->fab, fi, &fabd->dom, NULL);
    if (ret) {
	FT_PRINTERR("fi_domain", ret);
	goto err2;
    }
	
    ret = alloc_cm_res(fabd);
    if (ret)
	goto err4;

    ret = alloc_ep_res(fcd, fi);
    if (ret)
	goto err5;

    ret = bind_ep_res(fcd);
    if (ret)
	goto err6;

    /* Connect to server */
    ret = fi_connect(fcd->conn_ep, fi->dest_addr, NULL, 0);
    if (ret) {
	FT_PRINTERR("fi_connect", ret);
	goto err6;
    }
    
    /* Wait for the connection to be established */
    rd = fi_eq_sread(fabd->cmeq, &event, &entry, sizeof entry, -1, 0);
    if (rd != sizeof entry) {
	if (ret == -FI_EAVAIL) {
	    struct fi_eq_err_entry error = {0};
	    int rc = fi_eq_readerr(fabd->cmeq, &error, 0);
	    if (rc) {
		char buf[1024];
		fprintf(stderr, "error event: %s\n", fi_eq_strerror(fabd->cmeq, error.prov_errno,
      error.err_data, buf, 1024));
	    }
	} else {
	    FT_PRINTERR("fi_eq_sread", rd);
	}
	goto err6;
    }

    if (event != FI_CONNECTED || entry.fid != &fcd->conn_ep->fid) {
	FT_ERR("Unexpected CM event %d fid %p (ep %p)\n", event, entry.fid, fcd->conn_ep);
	ret = -FI_EOTHER;
	goto err6;
    }

    fi_freeinfo(fi);
    return 0;

err6:
    free_ep_res(fcd);
err5:
    fi_close(&fabd->cmeq->fid);
err4:
    fi_close(&fabd->dom->fid);
err2:
    fi_close(&fabd->fab->fid);
err1:
    fi_freeinfo(fi);
err0:
    return ret;
}

static int
initiate_conn(cm, svc, trans, attrs, fcd, conn_attr_list, no_more_redirect)
	CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
fabric_conn_data_ptr fcd;
attr_list conn_attr_list;
int no_more_redirect;
{
	int int_port_num;
	fabric_client_data_ptr fabd = (fabric_client_data_ptr) trans->trans_data;
	char *host_name;
	int remote_IP = -1;
	static int host_ip = 0;

	//fabric stuff

	if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_name)) {
		svc->trace_out(cm, "CMFABRIC transport found no IP_HOST attribute");
		host_name = NULL;
	} else {
		svc->trace_out(cm, "CMFABRIC transport connect to host %s", host_name);
	}
	if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_ip)) {
		svc->trace_out(cm, "CMFABRIC transport found no IP_ADDR attribute");
		/* wasn't there */
		host_ip = 0;
	} else {
		svc->trace_out(cm, "CMFABRIC transport connect to host_IP %lx", host_ip);
	}

	if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & int_port_num)) {
		svc->trace_out(cm, "CMFABRIC transport found no IP_PORT attribute");
//		return -1;
	} else {
		svc->trace_out(cm, "CMFABRIC transport connect to port %d", int_port_num);
	}

	client_connect(cm, svc, trans, attrs, fcd);


//here we write out the connection port to the other side. 
//for sockets thats all thats required. For IB we can use this to exchange information about the 
//IB parameters for the other side

	svc->trace_out(cm, "--> Connection established");
	fcd->remote_host = host_name == NULL ? NULL : strdup(host_name);
	fcd->remote_IP = remote_IP;
	fcd->remote_contact_port = int_port_num;
	fcd->fd = 0;
	fcd->fabd = fabd;

	memset(&fcd->last_write, 0, sizeof(rinfo));
	fcd->infocount = 0;

    

	add_attr(conn_attr_list, CM_THIS_CONN_PORT, Attr_Int4,
	         (attr_value) (long)int_port_num);
	add_attr(conn_attr_list, CM_PEER_IP, Attr_Int4,
	         (attr_value) (long)fcd->remote_IP);
/*	if (getpeername(sock, &sock_addr.s, &sock_len) == 0) {
		int_port_num = ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
		add_attr(conn_attr_list, CM_PEER_CONN_PORT, Attr_Int4,
		         (attr_value) (long)int_port_num);
		if (sock_addr.s.sa_family == AF_INET) {
			struct hostent *host;
			struct sockaddr_in *in_sock = (struct sockaddr_in *) &sock_addr;
			host = gethostbyaddr((char *) &in_sock->sin_addr,
			                     sizeof(struct in_addr), AF_INET);
			if (host != NULL) {
				fcd->remote_host = strdup(host->h_name);
				add_attr(conn_attr_list, CM_PEER_HOSTNAME, Attr_String,
				         (attr_value) strdup(host->h_name));
			}
		}
	}
*/
	svc->trace_out(fabd->cm, "Falling out of init conn\n");
	return 1;
}

/* 
 * Initiate a socket connection with another data exchange.  If port_num is -1,
 * establish a unix socket connection (name_str stores the file name of
 * the waiting socket).  Otherwise, establish an INET socket connection
 * (name_str stores the machine name).
 */
extern CMConnection
libcmfabric_LTX_initiate_conn(cm, svc, trans, attrs)
	CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
{
    fabric_conn_data_ptr fcd = create_fabric_conn_data(svc);
    attr_list conn_attr_list = create_attr_list();
    CMConnection conn;
    int fd, ret;

    fcd->fabd = trans->trans_data;

    if (initiate_conn(cm, svc, trans, attrs, fcd, conn_attr_list, 0) < 0)
	return NULL;

    add_attr(conn_attr_list, CM_PEER_LISTEN_PORT, Attr_Int4,
	     (attr_value) (long)fcd->remote_contact_port);
    conn = svc->connection_create(trans, fcd, conn_attr_list);
    fcd->conn = conn;

    if ((ret = fi_control (&fcd->rcq->fid, FI_GETWAIT, (void *) &fd))) {
	FT_PRINTERR("fi_control(FI_GETWAIT)", ret);
    }
    svc->trace_out(cm, "Cmfabric Adding trans->data_available as action on fd %d", fd);
    svc->fd_add_select(cm, fd, (select_list_func) CMFABRIC_data_available,
		       (void *) trans, (void *) conn);

    fcd->fd = fd;
    if ((ret = fi_control (&fcd->scq->fid, FI_GETWAIT, (void *) &fcd->sfd))) {
	FT_PRINTERR("fi_control(FI_GETWAIT)", ret);
    }
    return conn;
}

/* 
 * Check to see that if we were to attempt to initiate a connection as
 * indicated by the attribute list, would we be connecting to ourselves?
 * For sockets, this involves checking to see if the host name is the 
 * same as ours and if the IP_PORT matches the one we are listening on.
 */
extern int
libcmfabric_LTX_self_check(CManager cm, CMtrans_services svc, transport_entry trans, attr_list attrs)
{

    fabric_client_data_ptr fd = trans->trans_data;
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
		    /* value pointer */ (attr_value *)(long) & host_name)) {
	svc->trace_out(cm, "CMself check CMFABRIC transport found no IP_HOST attribute");
	host_name = NULL;
    }
    if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *)(long) & host_addr)) {
	svc->trace_out(cm, "CMself check CMFABRIC transport found no IP_ADDR attribute");
	if (host_name == NULL) return 0;
	host_addr = 0;
    }
    if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
		    /* value pointer */ (attr_value *)(long) & int_port_num)) {
	svc->trace_out(cm, "CMself check CMFABRIC transport found no IP_PORT attribute");
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
    if (int_port_num != fd->listen_port) {
	svc->trace_out(cm, "CMself check - Ports don't match, %d, %d", int_port_num, fd->listen_port);
	return 0;
    }
    svc->trace_out(cm, "CMself check returning TRUE");
    return 1;
}

extern int
libcmfabric_LTX_connection_eq(cm, svc, trans, attrs, fcd)
	CManager cm;
CMtrans_services svc;
transport_entry trans;
attr_list attrs;
fabric_conn_data_ptr fcd;
{

	int int_port_num;
	int requested_IP = -1;
	char *host_name = NULL;

	if (!query_attr(attrs, CM_IP_HOSTNAME, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & host_name)) {
		svc->trace_out(cm, "CMFABRIC transport found no IP_HOST attribute");
	}
	if (!query_attr(attrs, CM_IP_PORT, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & int_port_num)) {
		svc->trace_out(cm, "Conn Eq CMFABRIC transport found no IP_PORT attribute");
		return 0;
	}
	if (!query_attr(attrs, CM_IP_ADDR, /* type pointer */ NULL,
	                /* value pointer */ (attr_value *)(long) & requested_IP)) {
		svc->trace_out(cm, "CMFABRIC transport found no IP_ADDR attribute");
	}
	if (requested_IP == -1) {
		check_host(host_name, (void *) &requested_IP);
		requested_IP = ntohl(requested_IP);
		svc->trace_out(cm, "IP translation for hostname %s is %x", host_name,
		               requested_IP);
	}

	svc->trace_out(cm, "Socket Conn_eq comparing IP/ports %x/%d and %x/%d",
	               fcd->remote_IP, fcd->remote_contact_port,
	               requested_IP, int_port_num);
	if ((fcd->remote_IP == requested_IP) &&
	    (fcd->remote_contact_port == int_port_num)) {
		svc->trace_out(cm, "Socket Conn_eq returning TRUE");
		return 1;
	}
	svc->trace_out(cm, "Socket Conn_eq returning FALSE");
	return 0;
}


static void free_lres(fabric_client_data_ptr fd)
{
	fi_close(&fd->cmeq->fid);
}

static int alloc_cm_res(fabric_client_data_ptr fd)
{
	struct fi_eq_attr cm_attr;
	int ret;

	memset(&cm_attr, 0, sizeof cm_attr);
	cm_attr.wait_obj = FI_WAIT_FD;
	ret = fi_eq_open(fd->fab, &cm_attr, &fd->cmeq, NULL);
	if (ret)
		FT_PRINTERR("fi_eq_open", ret);

	return ret;
}

static void free_ep_res(fabric_conn_data_ptr fcd)
{
	fi_close(&fcd->conn_ep->fid);
	fi_close(&fcd->send_mr->fid);
	fi_close(&fcd->read_mr->fid);
	fi_close(&fcd->rcq->fid);
	fi_close(&fcd->scq->fid);
}

static int alloc_ep_res(fabric_conn_data_ptr fcd, struct fi_info *fi)
{
	fabric_client_data_ptr fabd = fcd->fabd;
	struct fi_cq_attr cq_attr;
	uint64_t access_mode;
	int ret;

	fcd->buffer_size = PIGGYBACK;
	
	fcd->read_buf = fabd->svc->get_data_buffer(fabd->cm, MAX(fcd->buffer_size, sizeof(uint64_t)));
	if (!fcd->read_buf) {
		perror("malloc");
		return -1;
	}
	fcd->max_credits = 512;
	memset(&cq_attr, 0, sizeof cq_attr);
	cq_attr.format = FI_CQ_FORMAT_DATA;
	cq_attr.size = fcd->max_credits << 1;
	if (fabd->send_waitset) {
	    cq_attr.wait_obj = FI_WAIT_SET;
	    cq_attr.wait_set = fabd->send_waitset;
	} else {
	    cq_attr.wait_obj = FI_WAIT_FD;
	}
	ret = fi_cq_open(fabd->dom, &cq_attr, &fcd->scq, NULL);
	if (ret) {
		FT_PRINTERR("fi_cq_open, on fcd->scq", ret);
		goto err1;
	}

	struct fi_cq_attr attrs;
	memset(&attrs, 0, sizeof(attrs));
	attrs.format = FI_CQ_FORMAT_DATA;
	attrs.wait_obj = FI_WAIT_FD;
	attrs.size = fcd->max_credits << 1;
	ret = fi_cq_open(fabd->dom, &cq_attr, &fcd->rcq, NULL);
	if (ret) {
		FT_PRINTERR("fi_cq_open", ret);
		goto err2;
	}
	
	access_mode = FI_REMOTE_READ;
	access_mode |= FI_RECV | FI_READ | FI_WRITE | FI_REMOTE_WRITE;
	fcd->send_buf = malloc(MAX(fcd->buffer_size, sizeof(uint64_t)));
	fcd->mapped_recv_buf = malloc(MAX(fcd->buffer_size, sizeof(uint64_t)));
	if (!fcd->send_buf) {
		perror("malloc");
		return -1;
	}
	ret = fi_mr_reg(fabd->dom, fcd->mapped_recv_buf, MAX(fcd->buffer_size, sizeof(uint64_t)), 
			access_mode, 0, 0, 0, &fcd->read_mr, NULL);
	if (ret) {
		FT_PRINTERR("fi_mr_reg", ret);
		goto err3;
	}

	access_mode = FI_REMOTE_WRITE | FI_WRITE;
	printf("fi_mr_reg length %lu, send_buf %p\n", MAX(fcd->buffer_size, sizeof(uint64_t)), fcd->send_buf);
	ret = fi_mr_reg(fabd->dom, fcd->send_buf, MAX(fcd->buffer_size, sizeof(uint64_t)), 
			access_mode, 0, 0, 0, &fcd->send_mr, NULL);
	if (ret) {
		FT_PRINTERR("fi_mr_reg", ret);
		goto err3;
	}

	ret = fi_endpoint(fabd->dom, fi, &fcd->conn_ep, NULL);
	if (ret) {
		FT_PRINTERR("fi_endpoint", ret);
		goto err4;
	}

	if (!fabd->cmeq) {
		ret = alloc_cm_res(fabd);
		if (ret)
			goto err4;
	}

	return 0;

err4:
	fi_close(&fcd->read_mr->fid);
	fi_close(&fcd->send_mr->fid);
err3:
	fi_close(&fcd->rcq->fid);
err2:
	fi_close(&fcd->scq->fid);
err1:
	free(fcd->send_buf);
	return ret;
}

static int bind_ep_res(fabric_conn_data_ptr fcd)
{
	int ret;

	ret = fi_ep_bind(fcd->conn_ep, &fcd->fabd->cmeq->fid, 0);
	if (ret) {
		FT_PRINTERR("fi_ep_bind", ret);
		return ret;
	}

	ret = fi_ep_bind(fcd->conn_ep, &fcd->scq->fid, FI_SEND);
	if (ret) {
		FT_PRINTERR("fi_ep_bind", ret);
		return ret;
	}

	ret = fi_ep_bind(fcd->conn_ep, &fcd->rcq->fid, FI_RECV);
	if (ret) {
		FT_PRINTERR("fi_ep_bind", ret);
		return ret;
	}

	ret = fi_enable(fcd->conn_ep);
	if (ret) {
		FT_PRINTERR("fi_enable", ret);
		return ret;
	}
	
	/* Post the first recv buffer */
	ret = fi_recv(fcd->conn_ep, fcd->mapped_recv_buf, fcd->buffer_size, fi_mr_desc(fcd->read_mr), 0, fcd->mapped_recv_buf);
	if (ret)
		FT_PRINTERR("fi_recv", ret);

	return ret;
}


typedef struct {
    /* OFI objects */
    /* int avtid; */
    /* struct fid_domain *domain; */
    /* struct fid_fabric *fabric; */
    /* struct fid_av     *av; */
    /* struct fid_ep     *ep; */
    /* struct fid_cq     *p2p_cq; */
    /* struct fid_cntr   *rma_ctr; */

    /* Queryable limits */
    uint64_t        max_buffered_send;
    uint64_t        max_buffered_write;
    uint64_t        max_send;
    uint64_t        max_write;
    uint64_t        max_short_send;
    uint64_t        max_mr_key_size;
    int             max_windows_bits;
    int             max_huge_rma_bits;
    int             max_huge_rmas;
    int             huge_rma_shift;
    int             context_shift;
    size_t          iov_limit;
    size_t          rma_iov_limit;

    /* /\* Mutexex and endpoints *\/ */
    /* MPIDI_OFI_cacheline_mutex_t mutexes[4]; */
    /* MPIDI_OFI_context_t         ctx[MPIDI_OFI_MAX_ENDPOINTS]; */

    /* /\* Window/RMA Globals *\/ */
    /* void                             *win_map; */
    /* uint64_t                          cntr; */
    /* MPIDI_OFI_atomic_valid_t  win_op_table[MPIDI_OFI_DT_SIZES][MPIDI_OFI_OP_SIZES]; */

    /* Active Message Globals */
    /* struct iovec                           am_iov[MPIDI_OFI_NUM_AM_BUFFERS]; */
    /* struct fi_msg                          am_msg[MPIDI_OFI_NUM_AM_BUFFERS]; */
    /* void                                  *am_bufs[MPIDI_OFI_NUM_AM_BUFFERS]; */
    /* MPIDI_OFI_am_repost_request_t  am_reqs[MPIDI_OFI_NUM_AM_BUFFERS]; */
    /* MPIDI_NM_am_target_handler_fn      am_handlers[MPIDI_OFI_MAX_AM_HANDLERS_TOTAL]; */
    /* MPIDI_NM_am_origin_handler_fn      am_send_cmpl_handlers[MPIDI_OFI_MAX_AM_HANDLERS_TOTAL]; */
    /* MPIU_buf_pool_t                       *am_buf_pool; */
    /* OPA_int_t                              am_inflight_inject_emus; */
    /* OPA_int_t                              am_inflight_rma_send_mrs; */

    /* Completion queue buffering */
    /* MPIDI_OFI_cq_buff_entry_t cq_buffered[MPIDI_OFI_NUM_CQ_BUFFERED]; */
    /* struct slist                      cq_buff_list; */
    /* int                               cq_buff_head; */
    /* int                               cq_buff_tail; */

    /* Process management and PMI globals */
    /* int    pname_set; */
    /* int    pname_len; */
    /* int    jobid; */
    /* char   addrname[FI_NAME_MAX]; */
    /* size_t addrnamelen; */
    /* char   kvsname[MPIDI_KVSAPPSTRLEN]; */
    /* char   pname[MPI_MAX_PROCESSOR_NAME]; */
    /* int    port_name_tag_mask[MPIR_MAX_CONTEXT_MASK]; */
} MPIDI_OFI_global_t;

MPIDI_OFI_global_t       MPIDI_Global;

static int server_listen(fabric_client_data_ptr fd, attr_list listen_info)
{
    struct fi_info *fi, *prov_use;
    CMtrans_services svc = fd->svc;
    int ret;
    int attr_port_num;
    char *port_str = NULL;

//    ret = fi_getinfo(FT_FIVERSION, fd->opts.src_addr, fd->opts.src_port, FI_SOURCE,
//		     fd->hints, &fi);
    if (listen_info != NULL
	&& !get_int_attr(listen_info, CM_IP_PORT, &attr_port_num)) {
	attr_port_num = 0;
    } else {
	if (attr_port_num > USHRT_MAX || attr_port_num < 0) {
	    fprintf(stderr, "Requested port number %d is invalid\n", attr_port_num);
	    return 1;
	}
	port_str = malloc(10);
	sprintf(port_str, "%d", attr_port_num);
    }

    ret = fi_getinfo(FT_FIVERSION, NULL, port_str, FI_SOURCE, fd->hints, &fi);

    if (((struct sockaddr_in *)fi->src_addr)->sin_addr.s_addr == htonl(INADDR_LOOPBACK) &&
	(strcmp(fi->fabric_attr->prov_name, "verbs") == 0)) {
	char host_name[256];
	fi_freeinfo(fi);
	
	if (listen_info) {
	    listen_info = attr_copy_list(listen_info);
	} else {
	    listen_info = create_attr_list();
	}
	set_string_attr(listen_info, CM_IP_INTERFACE, strdup("ib"));
	
	svc->trace_out(fd->cm, "CMFabric begin listen, requested port %d", attr_port_num);
	get_IP_config(host_name, sizeof(host_name), NULL, NULL, NULL, 
		      NULL, listen_info, svc->trace_out, (void *)fd->cm);
	ret = fi_getinfo(FT_FIVERSION, host_name, port_str, FI_SOURCE, fd->hints, &fi);
	svc->trace_out(fd->cm, "%s return value fi is %s\n", "server", fi_tostr(fi, FI_TYPE_INFO));
	fd->hostname = strdup(host_name);
	free_attr_list(listen_info);
    } else {
	svc->trace_out(fd->cm, "%s return value fi is %s\n", "server", fi_tostr(fi, FI_TYPE_INFO));
    }
	    
    prov_use = fi;
    MPIDI_Global.max_buffered_send  = prov_use->tx_attr->inject_size;
    MPIDI_Global.max_buffered_write = prov_use->tx_attr->inject_size;
    MPIDI_Global.max_send           = prov_use->ep_attr->max_msg_size;
    MPIDI_Global.max_write          = prov_use->ep_attr->max_msg_size;
    svc->trace_out(fd->cm, "Max send is %ld, max write is %ld\n", MPIDI_Global.max_send, MPIDI_Global.max_write);
    /* MPIDI_Global.iov_limit          = MIN(prov_use->tx_attr->iov_limit,MPIDI_OFI_IOV_MAX); */
    /* MPIDI_Global.rma_iov_limit      = MIN(prov_use->tx_attr->rma_iov_limit,MPIDI_OFI_IOV_MAX); */
    MPIDI_Global.max_mr_key_size    = prov_use->domain_attr->mr_key_size;

    /* if(MPIDI_Global.max_mr_key_size >= 8) { */
    /*     MPIDI_Global.max_windows_bits   = MPIDI_OFI_MAX_WINDOWS_BITS_64; */
    /*     MPIDI_Global.max_huge_rma_bits  = MPIDI_OFI_MAX_HUGE_RMA_BITS_64; */
    /*     MPIDI_Global.max_huge_rmas      = MPIDI_OFI_MAX_HUGE_RMAS_64; */
    /*     MPIDI_Global.huge_rma_shift     = MPIDI_OFI_HUGE_RMA_SHIFT_64; */
    /*     MPIDI_Global.context_shift      = MPIDI_OFI_CONTEXT_SHIFT_64; */
    /* } else if(MPIDI_Global.max_mr_key_size >= 4) { */
    /*     MPIDI_Global.max_windows_bits   = MPIDI_OFI_MAX_WINDOWS_BITS_32; */
    /*     MPIDI_Global.max_huge_rma_bits  = MPIDI_OFI_MAX_HUGE_RMA_BITS_32; */
    /*     MPIDI_Global.max_huge_rmas      = MPIDI_OFI_MAX_HUGE_RMAS_32; */
    /*     MPIDI_Global.huge_rma_shift     = MPIDI_OFI_HUGE_RMA_SHIFT_32; */
    /*     MPIDI_Global.context_shift      = MPIDI_OFI_CONTEXT_SHIFT_32; */
    /* } else if(MPIDI_Global.max_mr_key_size >= 2) { */
    /*     MPIDI_Global.max_windows_bits   = MPIDI_OFI_MAX_WINDOWS_BITS_16; */
    /*     MPIDI_Global.max_huge_rma_bits  = MPIDI_OFI_MAX_HUGE_RMA_BITS_16; */
    /*     MPIDI_Global.max_huge_rmas      = MPIDI_OFI_MAX_HUGE_RMAS_16; */
    /*     MPIDI_Global.huge_rma_shift     = MPIDI_OFI_HUGE_RMA_SHIFT_16; */
    /*     MPIDI_Global.context_shift      = MPIDI_OFI_CONTEXT_SHIFT_16; */
    /* } else { */
    /*     MPIR_ERR_SETFATALANDJUMP4(mpi_errno, */
    /*                               MPI_ERR_OTHER, */
    /*                               "**ofid_rma_init", */
    /*                               "**ofid_rma_init %s %d %s %s", */
    /*                               __SHORT_FILE__, */
    /*                               __LINE__, */
    /*                               FCNAME, */
    /*                               "Key space too small"); */
    /* } */
	if (ret) {
		FT_PRINTERR("fi_getinfo", ret);
		return ret;
	}

	ret = fi_fabric(fi->fabric_attr, &fd->fab, NULL);

	if (ret) {
		FT_PRINTERR("fi_fabric", ret);
		goto err0;
	}

	struct fi_wait_attr attrs;
	attrs.wait_obj = FI_WAIT_FD;
	attrs.flags = 0;
	ret = fi_wait_open(fd->fab, &attrs, &fd->send_waitset);
	if (ret) {
	    /* sigh */
	    fd->send_waitset = NULL;
	}

	ret = fi_passive_ep(fd->fab, fi, &fd->listen_ep, NULL);
	if (ret) {
		FT_PRINTERR("fi_passive_ep", ret);
		goto err1;
	}

	ret = alloc_cm_res(fd);
	if (ret)
		goto err2;

	ret = fi_pep_bind(fd->listen_ep, &fd->cmeq->fid, 0);
	if (ret) {
		FT_PRINTERR("fi_pep_bind", ret);
		goto err3;
	}


	ret = fi_listen(fd->listen_ep);
	if (ret) {
		FT_PRINTERR("fi_listen", ret);
		goto err3;
	}

	fi_freeinfo(fi);
	return 0;
err3:
	free_lres(fd);
err2:
	fi_close(&fd->listen_ep->fid);
err1:
	fi_close(&fd->fab->fid);
err0:
	fi_freeinfo(fi);
	return ret;
}

static int server_connect(fabric_conn_data_ptr fcd)
{
	struct fi_eq_cm_entry entry;
	uint32_t event;
	struct fi_info *info = NULL;
	ssize_t rd;
	int ret;
	fabric_client_data_ptr fabd = fcd->fabd;

	rd = fi_eq_sread(fabd->cmeq, &event, &entry, sizeof entry, -1, 0);
	if (rd != sizeof entry) {
	    if (rd == -FI_EAVAIL) {
		struct fi_eq_err_entry error = {0};
		int rc = fi_eq_readerr(fabd->cmeq, &error, 0);
		if (rc) {
		    char buf[1024];
		    fprintf(stderr, "error event: %s\n", fi_eq_strerror(fabd->cmeq, error.prov_errno,
									error.err_data, buf, 1024));
		}
	    } else {
		FT_PRINTERR("fi_eq_sread", rd);
	    }
	    return (int) rd;
	}

	info = entry.info;
	if (event != FI_CONNREQ) {
		fprintf(stderr, "Unexpected CM event %d\n", event);
		ret = -FI_EOTHER;
		goto err1;
	}

	ret = fi_domain(fabd->fab, info, &fabd->dom, NULL);
	if (ret) {
		FT_PRINTERR("fi_domain", ret);
		goto err1;
	}


	ret = fi_endpoint(fabd->dom, info, &fcd->conn_ep, NULL);
	if (ret) {
		FT_PRINTERR("fi_endpoint", -ret);
		goto err1;
	}

	ret = alloc_ep_res(fcd, info);
	if (ret)
		 goto err1;

	ret = bind_ep_res(fcd);
	if (ret)
		goto err3;

	ret = fi_accept(fcd->conn_ep, NULL, 0);
	if (ret) {
		FT_PRINTERR("fi_accept", ret);
		goto err3;
	}

	rd = fi_eq_sread(fabd->cmeq, &event, &entry, sizeof entry, -1, 0);
 	if (rd != sizeof entry) {
	    if (ret == -FI_EAVAIL) {
		struct fi_eq_err_entry error = {0};
		int rc = fi_eq_readerr(fabd->cmeq, &error, 0);
		if (rc) {
		    char buf[1024];
		    fprintf(stderr, "error event: %s\n", fi_eq_strerror(fabd->cmeq, error.prov_errno,
									error.err_data, buf, 1024));
		}
	    } else {
		FT_PRINTERR("fi_eq_sread", rd);
	    }
	    goto err3;
 	}

	if (event != FI_CONNECTED || entry.fid != &fcd->conn_ep->fid) {
		fprintf(stderr, "Unexpected CM event %d fid %p (ep %p)\n",
			event, entry.fid, fcd->conn_ep);
 		ret = -FI_EOTHER;
 		goto err3;
 	}
 
 	fi_freeinfo(info);
 	return 0;

err3:
	free_ep_res(fcd);
err1:
 	fi_reject(fabd->listen_ep, info->handle, NULL, 0);
 	fi_freeinfo(info);
 	return ret;
}

/* 
 * Create an passive endpoint to listen for connections
 */
extern attr_list
libcmfabric_LTX_non_blocking_listen(CManager cm, CMtrans_services svc, transport_entry trans, attr_list listen_info)
{
    fabric_client_data_ptr fd = trans->trans_data;
    int wait_sock = -1;
    int ret, IP, port_num;
    size_t addrlen;
    struct sockaddr_in local_addr;

    if (cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, cm));
    }
    int result = server_listen(fd, listen_info);
    if (result != 0) {
	fprintf(stderr, "Cannot bind INET socket\n");
	return NULL;
    }

    addrlen = sizeof(local_addr);
    ret = fi_getname(&fd->listen_ep->fid, (void*)&local_addr, &addrlen);
    IP = ntohl(local_addr.sin_addr.s_addr);
    port_num = ntohs(local_addr.sin_port);
    if (ret) {
	FT_PRINTERR("fi_getname", ret);
	return NULL;
    }

    ret = fi_control (&fd->cmeq->fid, FI_GETWAIT, (void *) &wait_sock);
    if (ret) {
	FT_PRINTERR("fi_control(FI_GETWAIT)", ret);
    } else {
	svc->trace_out(cm, "Cmfabric Adding fabric_service_incoming as action on fd %d", wait_sock);
	svc->fd_add_select(cm, wait_sock, fabric_service_incoming,
			   (void *) trans, (void *) fd->listen_ep);
    }
    {
	attr_list ret_list;
	
	svc->trace_out(cm, "CMFABRIC listen succeeded on port %d, fd %d",
		       port_num, wait_sock);
	ret_list = create_attr_list();
	
	fd->listen_port = port_num;
	add_attr(ret_list, CM_TRANSPORT, Attr_String,
		 (attr_value) strdup("fabric"));
	if ((getenv("CMFabricUseHostname") != NULL) || 
	    (getenv("CM_NETWORK") != NULL)) {
	    add_attr(ret_list, CM_IP_HOSTNAME, Attr_String,
		     (attr_value) strdup(fd->hostname));
	} else if (IP == 0) {
	    add_attr(ret_list, CM_IP_ADDR, Attr_Int4, 
		     (attr_value)INADDR_LOOPBACK);
	} else {
	    add_int_attr(ret_list, CM_IP_ADDR, (int)IP);
	}
	add_attr(ret_list, CM_IP_PORT, Attr_Int4,
		 (attr_value) (long)port_num);
	
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
libcmfabric_LTX_set_write_notify(trans, svc, fcd, enable)
	transport_entry trans;
CMtrans_services svc;
fabric_conn_data_ptr fcd;
int enable;
{
	if (enable != 0) {
		svc->fd_write_select(trans->cm, fcd->fd, (select_list_func) trans->write_possible,
		                     (void *)trans, (void *) fcd->conn);
	} else {
		/* remove entry */
		svc->fd_write_select(trans->cm, fcd->fd, NULL, NULL, NULL);
	}   
}


extern CMbuffer
libcmfabric_LTX_read_block_func(CMtrans_services svc, fabric_conn_data_ptr fcd, int *len_ptr, int *offset_ptr)
{
    *len_ptr = fcd->read_buffer_len;
    if (fcd->read_buf) {
	CMbuffer tmp = fcd->read_buf;
	fcd->read_buf = NULL;
	fcd->read_buffer_len = 0;
	if (offset_ptr) *offset_ptr = fcd->read_offset;
	return tmp;
    }
    return NULL;  
}


#ifndef IOV_MAX
/* this is not defined in some places where it should be.  Conservative. */
#define IOV_MAX 16
#endif

extern int
libcmfabric_LTX_writev_complete_notify_func(CMtrans_services svc, 
					fabric_conn_data_ptr fcd,
					void *iovs,
					int iovcnt,
					attr_list attrs,
					CMcompletion_notify_func notify_func,
					void *notify_client_data)
{
    	int fd = fcd->fd;
	int left = 0;
	int iget = 0;
	int i;
	int can_reuse_mapping = 0;
	struct iovec * iov = (struct iovec*) iovs, *tmp_iov;
	rinfo *last_write_request = &fcd->last_write;
	struct control_message msg; 
    
	for (i = 0; i < iovcnt; i++) {
	    left += iov[i].iov_len;
	}

	svc->trace_out(fcd->fabd->cm, "CMFABRIC writev of %d bytes on fd %d",
	               left, fd);
	
	if (left + msg_offset() < PIGGYBACK)
	{
		//total size is less than the piggyback size
		iget = internal_write_piggyback(svc, fcd, left, iov, iovcnt);
		if (notify_func) {
		    (notify_func)(notify_client_data);
		}
		if(iget < 0)
		{
			svc->trace_out(fcd->fabd->cm, "CMFABRIC error in writing piggyback");
			return -1;
		}
		if(iget == 0)
		{
			return iovcnt;
		}
		return -1;
	}

	svc->set_pending_write(fcd->conn);

	if (notify_func) {
	    can_reuse_mapping = 1;
	    /* OK, we're not going to copy the data */
	    if (last_write_request->iovcnt == iovcnt) {
		int i;
		for(i=0; i < last_write_request->iovcnt; i++) {
		    if ((iov[i].iov_len != last_write_request->iov[i].iov_len) ||
			(iov[i].iov_base != last_write_request->iov[i].iov_base)) {
			can_reuse_mapping = 0;
			svc->trace_out(fcd->fabd->cm, "CMFABRIC already mapped data, doesn't match write, buf %d, %p vs. %p, %d vs. %d",
				       i, iov[i].iov_base, last_write_request->iov[i].iov_base, iov[i].iov_len, last_write_request->iov[i].iov_len);
		    break;
		    }
		}
	    } else {
		svc->trace_out(fcd->fabd->cm, "CMFABRIC either no already mapped data, or wrong buffer count");
		can_reuse_mapping = 0;
	    }
	} else {
	    svc->trace_out(fcd->fabd->cm, "CMFABRIC User-owned data with no notify, so no reuse\n");
	}
#ifndef DO_DEREG_ON_FINISH
	if (last_write_request->iovcnt && !(can_reuse_mapping)) {
	    for(i = 0; i < last_write_request->iovcnt; i ++)
	    {
		fi_close(&last_write_request->mrlist[i]->fid);
		last_write_request->mrlist[i] = NULL;
	    }
	    
	    free(last_write_request->mrlist);
	    free(last_write_request->iov);
	    last_write_request->iov = NULL;
	    last_write_request->iov = NULL;
	    last_write_request->mrlist = NULL;
	    last_write_request->iovcnt = 0;
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

    
	svc->trace_out(fcd->fabd->cm, "Can reuse mapping is %d\n", can_reuse_mapping);
	if (!can_reuse_mapping) {
	    int i;
	    fcd->last_write.iovcnt = iovcnt;
	    if (fcd->last_write.iov) free(fcd->last_write.iov);
	    fcd->last_write.iov = malloc(sizeof(tmp_iov[0]) * iovcnt);
	    memcpy(fcd->last_write.iov, tmp_iov, sizeof(tmp_iov[0]) * iovcnt);

	    fcd->last_write.mrlist = malloc(iovcnt * sizeof(fcd->last_write.mrlist[0]));
	    for (i=0; i < iovcnt; i++) {
		int ret;
		svc->trace_out(fcd->fabd->cm, "fi_mr_reg %d, addr %p, len %ld, with attrs REMOTE_READ,REMOTE_WRITE, SEND, RECV)\n", i, tmp_iov[i].iov_base, tmp_iov[i].iov_len);
		ret = fi_mr_reg(fcd->fabd->dom, tmp_iov[i].iov_base, tmp_iov[i].iov_len, 
				FI_REMOTE_READ | FI_REMOTE_WRITE | FI_SEND | FI_RECV, 0, 0, 0, &fcd->last_write.mrlist[i], NULL);
		if (ret) {
		    FT_PRINTERR("fi_mr_reg of elements in bulk write", ret);
		    return -1;
		}

	    }
	}
	fcd->last_write.notify_func = notify_func;
	fcd->last_write.notify_client_data = notify_client_data;
    
	if(fcd->last_write.mrlist == NULL)
	{
		svc->trace_out(fcd->fabd->cm, "CMFABRIC writev error in registereing memory");
		return -1;
	}
	
	svc->trace_out(fcd->fabd->cm, "FIrst 16 bytes of send buffer (len %d) are %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x\n", left, ((unsigned char*)iov[0].iov_base)[0], ((unsigned char*)iov[0].iov_base)[1], ((unsigned char*)iov[0].iov_base)[2], ((unsigned char*)iov[0].iov_base)[3], ((unsigned char*)iov[0].iov_base)[4], ((unsigned char*)iov[0].iov_base)[5], ((unsigned char*)iov[0].iov_base)[6], ((unsigned char*)iov[0].iov_base)[7], ((unsigned char*)iov[0].iov_base)[8], ((unsigned char*)iov[0].iov_base)[9], ((unsigned char*)iov[0].iov_base)[10], ((unsigned char*)iov[0].iov_base)[11], ((unsigned char*)iov[0].iov_base)[12], ((unsigned char*)iov[0].iov_base)[13], ((unsigned char*)iov[0].iov_base)[14], ((unsigned char*)iov[0].iov_base)[15]);
	iget = internal_write_request(svc, fcd, left, &fcd->last_write);
	if(iget < 0)
	{
		svc->trace_out(fcd->fabd->cm, "CMFABRIC error in writing request");
		return -1;
	}

	return iovcnt;
}

extern int
libcmfabric_LTX_writev_func(svc, fcd, iovs, iovcnt, attrs)
CMtrans_services svc;
fabric_conn_data_ptr fcd;
void *iovs;
int iovcnt;
attr_list attrs;
{
    return libcmfabric_LTX_writev_complete_notify_func(svc, fcd, iovs, iovcnt, 
						   attrs, NULL, NULL);
}

static void
free_fabric_data(CManager cm, void *fdv)
{
    fabric_client_data_ptr fd = (fabric_client_data_ptr) fdv;
    CMtrans_services svc = fd->svc;
    kill_thread(fd);
    if (fd->hostname != NULL)
	svc->free_func(fd->hostname);
    svc->free_func(fd);
}

static void
check_completed_pull(CManager cm, fabric_client_data_ptr fabd)
{
    msg_item_p msg;
    fabric_conn_data_ptr fcd;
    CMtrans_services svc = fabd->svc;
    if (!fabd->completed_queue) return;

    pthread_mutex_lock(&fabd->pull_queue_mutex);
    msg = fabd->completed_queue;
    fabd->completed_queue = msg->next;
    pthread_mutex_unlock(&fabd->pull_queue_mutex);

    svc->acquire_CM_lock(fabd->cm, __FILE__, __LINE__);
    /* OK, msg is ours alone */
    fcd = msg->conn_data;
    free(msg->pulls);
    CMbuffer cb = svc->create_data_and_link_buffer(fcd->fabd->cm, 
						   msg->buffer, 
						   msg->message_size);
    cb->return_callback = free_func;
    cb->return_callback_data = (void*)msg->mr_rec;
    fcd->read_buf = cb;
    svc->trace_out(fcd->fabd->cm, "FIrst 16 bytes of receive buffer (len %d) are %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x %02x \n", msg->message_size, ((unsigned char*)cb->buffer)[0], ((unsigned char*)cb->buffer)[1], ((unsigned char*)cb->buffer)[2], ((unsigned char*)cb->buffer)[3], ((unsigned char*)cb->buffer)[4], ((unsigned char*)cb->buffer)[5], ((unsigned char*)cb->buffer)[6], ((unsigned char*)cb->buffer)[7], ((unsigned char*)cb->buffer)[8], ((unsigned char*)cb->buffer)[9], ((unsigned char*)cb->buffer)[10], ((unsigned char*)cb->buffer)[11], ((unsigned char*)cb->buffer)[12], ((unsigned char*)cb->buffer)[13], ((unsigned char*)cb->buffer)[14], ((unsigned char*)cb->buffer)[15]);
    fcd->read_buffer_len = msg->message_size;
    svc->trace_out(fcd->fabd->cm, "CMFABRIC handle_request completed");    
    internal_write_response(svc, fcd, msg->message_size, msg->request_ID);
    fabd->trans->data_available(fabd->trans, fcd->conn);
    svc->return_data_buffer(fabd->trans->cm, cb);
    svc->drop_CM_lock(fabd->cm, __FILE__, __LINE__);

}

extern
void libcmfabric_LTX_install_pull_schedule(CMtrans_services svc,
					   transport_entry trans,
					   struct timeval *base_time, 
					   struct timeval *period, 
					   CMavail_period_ptr avail)
{
    fabric_client_data_ptr fabd = trans->trans_data;
    fabd->pull_schedule_base = *base_time;
    fabd->pull_schedule_period = *period;
    CMavail_period_ptr tmp = fabd->avail;
    fabd->avail = avail;
    free(tmp);
    if (!fabd->thread_init) {
	svc->trace_out(fabd->cm, "Starting pull thread!\n");
	pthread_mutex_init(&fabd->pull_queue_mutex, NULL);
	fabd->thread_should_run = 1;
	pthread_t thr;
	if (!fabd->send_waitset) {
	    /* no waitset support */
	    int filedes[2];
	    if (pipe(filedes) != 0) {
		perror("Pipe for wake not created.  Wake mechanism inoperative.");
		return;
	    }
	    fabd->wake_read_fd = filedes[0];
	    fabd->wake_write_fd = filedes[1];
	    fabd->nfds = fabd->wake_read_fd;
	    FD_SET(fabd->wake_read_fd, &fabd->readset);
	    fabd->fcd_array_by_sfd = malloc(sizeof(char*));
	    fabd->cq_array = malloc(sizeof(char*));
	}
	svc->add_poll(fabd->cm, (CMPollFunc) check_completed_pull, fabd);
	pthread_create(&thr, NULL, (void*(*)(void*))&pull_thread, fabd);
	fabd->thread_init = 1;
    }
}

extern void *
libcmfabric_LTX_initialize(CManager cm, CMtrans_services svc)
{
    static int atom_init = 0;

    fabric_client_data_ptr fabd;
    svc->trace_out(cm, "Initialize CM fabric transport built in %s\n",
		   EVPATH_MODULE_BUILD_DIR);
    if (atom_init == 0) {
	CM_IP_HOSTNAME = attr_atom_from_string("IP_HOST");
	CM_IP_PORT = attr_atom_from_string("IP_PORT");
	CM_IP_ADDR = attr_atom_from_string("IP_ADDR");
	CM_IP_INTERFACE = attr_atom_from_string("IP_INTERFACE");
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
    fabd = svc->malloc_func(sizeof(struct fabric_client_data));
    fabd->svc = svc;
    memset(fabd, 0, sizeof(struct fabric_client_data));
    fabd->cm = cm;
    fabd->hostname = NULL;
    fabd->listen_port = -1;
    fabd->svc = svc;
    fabd->port = 1; //need to somehow get proper port here
    
    fabd->psn = lrand48()%256;

    fabd->hints = fi_allocinfo();
    
//	fabd->hints->cap		= FI_CONTEXT;
    fabd->hints->ep_attr->type	= FI_EP_MSG;
    fabd->hints->caps		= FI_MSG | FI_RMA;
//	fabd->hints->caps		= FI_MSG;
    fabd->hints->mode		= FI_CONTEXT | FI_LOCAL_MR | FI_RX_CQ_DATA;
    fabd->hints->addr_format	= FI_SOCKADDR;
//	fabd->hints->tx_attr->op_flags  = FI_DELIVERY_COMPLETE | FI_COMPLETION;
    fabd->hints->tx_attr->op_flags  = FI_COMPLETION;

    fabd->hints->domain_attr->mr_mode = FI_MR_BASIC;
    fabd->hints->domain_attr->threading        =  FI_THREAD_SAFE;
    fabd->hints->domain_attr->control_progress =  FI_PROGRESS_AUTO;
    fabd->hints->domain_attr->data_progress    =  FI_PROGRESS_AUTO;
    svc->add_shutdown_task(cm, free_fabric_data, (void *) fabd, FREE_TASK);
    
    fabd->wake_read_fd = -1;
    EVPATH_FD_ZERO(&fabd->readset);
    fabd->nfds = 0;
    return (void *) fabd;
}


extern transport_entry
cmfabric_add_static_transport(CManager cm, CMtrans_services svc)
{
    transport_entry transport;
    transport = svc->malloc_func(sizeof(struct _transport_item));
    memset(transport, 0, sizeof(*transport));
    transport->trans_name = strdup("fabric");
    transport->cm = cm;
    transport->transport_init = (CMTransport_func)libcmfabric_LTX_initialize;
    transport->listen = (CMTransport_listen_func)libcmfabric_LTX_non_blocking_listen;
    transport->initiate_conn = (CMConnection(*)())libcmfabric_LTX_initiate_conn;
    transport->self_check = (int(*)())libcmfabric_LTX_self_check;
    transport->connection_eq = (int(*)())libcmfabric_LTX_connection_eq;
    transport->shutdown_conn = (CMTransport_shutdown_conn_func)libcmfabric_LTX_shutdown_conn;
    transport->read_block_func = (CMTransport_read_block_func)libcmfabric_LTX_read_block_func;
    transport->read_to_buffer_func = (CMTransport_read_to_buffer_func)NULL;
    transport->writev_func = (CMTransport_writev_func)libcmfabric_LTX_writev_func;
    transport->writev_complete_notify_func = (CMTransport_writev_complete_notify_func)libcmfabric_LTX_writev_complete_notify_func;
    transport->install_pull_schedule_func = (CMTransport_install_pull_schedule)libcmfabric_LTX_install_pull_schedule;
    transport->get_transport_characteristics = NULL;
    if (transport->transport_init) {
	transport->trans_data = transport->transport_init(cm, svc, transport);
    }
    return transport;
}

#endif

