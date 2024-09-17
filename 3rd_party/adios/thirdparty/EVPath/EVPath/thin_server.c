#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifdef HAVE_WINSOCK2_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#else
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "evpath.h"
#include "cm_internal.h"

#if defined (__INTEL_COMPILER)
//  Allow int conversions
#  pragma warning (disable: 2259)
#endif
static void socket_accept_thin_client(void *cmv, void * sockv);
extern void CMget_qual_hostname(CManager cm, char *buf, int len);

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif
extern void
CMget_port_range(CManager cm, int *high_bound, int *low_bound);

extern int
EVthin_socket_listen(CManager cm,  char **hostname_p, int *port_p)
{

    unsigned int length;
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    SOCKET conn_sock;
    int int_port_num = 0;
    u_short port_num = 0;
    char host_name[256];
    int high_bound, low_bound;
    CMget_port_range(cm, &high_bound, &low_bound);

    conn_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (conn_sock == SOCKET_ERROR) {
	fprintf(stderr, "Cannot open INET socket\n");
	return 0;
    }
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = (unsigned short) (htons(port_num));
    if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	fprintf(stderr, "Failed to set 1REUSEADDR on INET socket\n");
	return 0;
    }
    {
	long seedval = (long)time(NULL) + (long)getpid();
	/* port num is free.  Constrain to range to standards */
	int size = high_bound - low_bound;
	int tries = 30;
	int result = SOCKET_ERROR;
	srand48(seedval);
	while (tries > 0) {
	    int target = low_bound + (int)(size * drand48());
	    sock_addr.sin_port = htons(target);
	    CMtrace_out(cm, CMConnectionVerbose, "CMSocket trying to bind port %d", target);
	    result = bind(conn_sock, (struct sockaddr *) &sock_addr,
			  sizeof sock_addr);
	    int_port_num = target;
	    tries--;
	    if (result != SOCKET_ERROR) tries = 0;
	    if (tries%5 == 4) {
		/* try reseeding in case we're in sync with another process */
		srand48(time(NULL) + getpid());
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
    }
    /*    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
	     sizeof sock_addr) == SOCKET_ERROR) {
	fprintf(stderr, "Cannot bind INET socket\n");
	return 0;
	}*/
    sock_opt_val = 1;
    if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	perror("Failed to set 2REUSEADDR on INET socket");
	return 0;
    }
    length = sizeof sock_addr;
    if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
	fprintf(stderr, "Cannot get socket name\n");
	return 0;
    }
    /* begin listening for conns and set the backlog */
    if (listen(conn_sock, FD_SETSIZE)) {
	fprintf(stderr, "listen failed\n");
	return 0;
    }
    /* set the port num as one we can be contacted at */
    
    CM_fd_add_select(cm, conn_sock, socket_accept_thin_client,
		    (void *) cm, (void *) (intptr_t)conn_sock);
    

    //    int_port_num = (unsigned short)(ntohs(sock_addr.sin_port));

    CMget_qual_hostname(cm, host_name, sizeof(host_name));
	
    *hostname_p = strdup(host_name);
    *port_p = int_port_num;
    return 1;
}


typedef struct thin_conn {
    FFSFile ffsfile;
    SOCKET fd;
    int target_stone;
    int format_count;
    FMStructDescList *format_list;
    int max_src_list;
    EVsource *src_list;
} *thin_conn_data;

static void thin_free_func(void *event_data, void *client_data)
{
    (void) client_data;
    free(event_data);
}

static void
thin_data_available(void *cmv, void * conn_datav)
{
    thin_conn_data cd = conn_datav;
    CManager cm = (CManager) cmv;
    int i;

    CManager_unlock(cm);
    switch(FFSnext_record_type(cd->ffsfile)) {
    case FFSindex:
	break;
    case FFSend:
    case FFSerror:
	close_FFSfile(cd->ffsfile);
	free_FFSfile(cd->ffsfile);
	for (i=0; i < cd->format_count; i++) {
	    int j = 0;
	    if (cd->format_list[i] == NULL) continue;
	    while((cd->format_list[i])[j].format_name != NULL) {
		int k = 0;
		free((void*)cd->format_list[i][j].format_name);
		while (cd->format_list[i][j].field_list[k].field_name != NULL) {
		    free((char*)cd->format_list[i][j].field_list[k].field_name);
		    free((char*)cd->format_list[i][j].field_list[k].field_type);
		    k++;
		}
		free(cd->format_list[i][j].field_list);
		j++;
	    }
	    free(cd->format_list[i]);
	}
	free(cd->format_list);
	for (i=0; i <= cd->max_src_list; i++) {
	    if (cd->src_list[i] != NULL) {
		EVfree_source(cd->src_list[i]);
	    }
	}
	free(cd->src_list);
	CM_fd_remove_select(cm, cd->fd);
	free(cd);
	break;
    case FFSformat: {
	FFSTypeHandle next_format = FFSread_format(cd->ffsfile);
	FMStructDescList formats = 
	    get_localized_formats(FMFormat_of_original(next_format));
	FFSTypeHandle target = 
	    FFSset_fixed_target(FFSContext_of_file(cd->ffsfile), 
				formats);
	int format_num = FMformat_index(FMFormat_of_original(target));
	if (cd->format_list == NULL) {
	    cd->format_list = malloc(sizeof(cd->format_list[0]));
	    cd->format_count = 1;
	}
	if (cd->format_count < format_num) {
	    cd->format_list = 
		realloc(cd->format_list, 
			(format_num + 1) * sizeof(cd->format_list[0]));
	    memset(cd->format_list + cd->format_count, 0,
		   sizeof(cd->format_list[0]) * (format_num -cd->format_count +1));
	    cd->format_count = format_num + 1;
	}
	cd->format_list[format_num] = formats;
	break;
    }
    case FFSdata: {
	FFSTypeHandle next_format = FFSnext_type_handle(cd->ffsfile);
	size_t len = FFSnext_data_length(cd->ffsfile);
	int format_num = FMformat_index(FMFormat_of_original(next_format));
	void *data = malloc(len);
	FFSread(cd->ffsfile, data);
	if (cd->max_src_list < format_num) {
	    cd->src_list = realloc(cd->src_list, 
				   (format_num+1) * sizeof(cd->src_list[0]));
	    memset(&cd->src_list[cd->max_src_list], 0,
		   (format_num - cd->max_src_list + 1) * sizeof(cd->src_list[0]));
	    cd->max_src_list = format_num;
	}
	if (cd->src_list[format_num] == NULL) {
	    cd->src_list[format_num] = 
		EVcreate_submit_handle_free(cm, cd->target_stone, 
					    cd->format_list[format_num],
					    thin_free_func, cd);
	}
	EVsubmit(cd->src_list[format_num], data, NULL);
	break;
    }
    case FFScomment: {
	char *comment = FFSread_comment(cd->ffsfile);
	if (strncmp(comment, "Stone ", 6) == 0) {
	    int tmp_stone;
	    if (sscanf(comment, "Stone %d", &tmp_stone) == 1) {
		cd->target_stone = tmp_stone;
	    }
	}
	break;
    }
    }
    CManager_lock(cm);
}

static void
socket_accept_thin_client(void *cmv, void * sockv)
{
    CManager cm = (CManager) cmv;
    SOCKET conn_sock = (int) (intptr_t)sockv;
    SOCKET sock;
    struct sockaddr sock_addr;
    unsigned int sock_len = sizeof(sock_addr);
    int int_port_num;
    struct linger linger_val;
    int sock_opt_val = 1;
    thin_conn_data cd;

#ifdef TCP_NODELAY
    int delay_value = 1;
#endif

    linger_val.l_onoff = 1;
    linger_val.l_linger = 60;
    if ((sock = accept(conn_sock, (struct sockaddr *) 0, (unsigned int *) 0)) == SOCKET_ERROR) {
	perror("Cannot accept socket connection");
	CM_fd_remove_select(cm, conn_sock);
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

    sock_len = sizeof(sock_addr);
    memset(&sock_addr, 0, sock_len);
    getsockname(sock, (struct sockaddr *) &sock_addr, &sock_len);
    int_port_num = (unsigned short) ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);

    memset(&sock_addr, 0, sizeof(sock_addr));
    sock_len = sizeof(sock_addr);
    if (getpeername(sock, &sock_addr,
		    &sock_len) == 0) {
	int_port_num = (unsigned short) ntohs(((struct sockaddr_in *) &sock_addr)->sin_port);
    }

    (void) int_port_num;
    cd = malloc(sizeof(*cd));
    memset(cd, 0, sizeof(*cd));
    cd->ffsfile = open_FFSfd((void*)(intptr_t)sock, "r");
    cd->fd = sock;
    cd->src_list = malloc(sizeof(cd->src_list[0]));
    cd->src_list[0] = NULL;
    INT_CM_fd_add_select(cm, sock,
		     (void (*)(void *, void *)) thin_data_available,
		       (void *) cm, (void *) cd);
}
