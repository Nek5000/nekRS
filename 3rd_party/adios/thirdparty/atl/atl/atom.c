#include "config.h"

#undef NDEBUG
#  include <assert.h>
#  include <string.h>
#  include <stdlib.h>
#  ifdef HAVE_MALLOC_H
#    include <malloc.h>
#  endif
#  include <stdio.h>
#  include <stdint.h>
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#  include <errno.h>
#  include <sys/types.h>
#  ifndef HAVE_WINDOWS_H
#    include <netdb.h>
#    include <sys/socket.h>
#    include <netinet/in.h>
#    include <netinet/tcp.h>
#    include <arpa/inet.h>
#  else
#    include <winsock2.h>
#    include <ws2tcpip.h>
#    pragma comment(lib, "Ws2_32.lib")
#  endif
#  include <fcntl.h>

#ifdef _MSC_VER
    #define strdup _strdup
    #include <io.h>
#pragma warning(disable: 4996)
#endif

#include "atl.h"
#include "atom_internal.h"

#define MAXDATASIZE 100
#include "tclHash.h"

#ifndef HAVE_WINDOWS_H
typedef int SOCKET;
#else
#define write(fs, addr, len) send(fs, addr, len, 0)
#define read(fs, addr, len) recv(fs, addr, len, 0)
#define close(fs) closesocket(fs);
#endif

/* opaque type for atom server handle */
typedef struct _atom_server {
    int sockfd;
    SOCKET tcp_fd;
    int use_tcp;
    int no_server;
    struct hostent *he;
    struct sockaddr_in their_addr;
    int flags;
    char *server_id;
    Tcl_HashTable string_hash_table;
    Tcl_HashTable value_hash_table;
} atom_server_struct;

static char *atom_server_host = NULL;
static int establish_server_connection(atom_server as, int do_fallback);

#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif

#ifndef O_NONBLOCK
#define O_NONBLOCK 0x80
#endif
static void
set_blocking(atom_server as, int block)
{
    if (block && ((as->flags & O_NONBLOCK) == 0)) {
	return;			/* already blocking */
    }
    if (!block && ((as->flags & O_NONBLOCK) == O_NONBLOCK)) {
	return;			/* already non-blocking */
    }
    if (block) {
	as->flags &= (~O_NONBLOCK);
    } else {
	as->flags |= O_NONBLOCK;
    }
#ifndef HAVE_WINDOWS_H
    if (fcntl(as->sockfd, F_SETFL, as->flags) < 0) {
	perror("fcntl");
	exit(1);
    }
    if (as->tcp_fd > 0) {
	if (fcntl(as->tcp_fd, F_SETFL, as->flags) < 0) {
	    perror("TCP_FD fcntl");
	}
    }
#else
    u_long block_val = !block;
    if (ioctlsocket(as->sockfd, FIONBIO, &block_val) != 0) {
	perror("ioctlsocket");
	exit(1);
    }
    ioctlsocket(as->tcp_fd, FIONBIO, &block_val);
#endif

}

static void
handle_unexpected_msg(atom_server as, char *msg)
{
    switch (msg[0]) {
    case 'E':{
	    Tcl_HashEntry *entry = NULL;
	    char *str;
	    int atom;
	    atom = strtol(&msg[1], &str, 10);
	    str++;
	    entry = Tcl_FindHashEntry(&as->string_hash_table, str);
	    if (entry != NULL) {
		send_get_atom_msg_ptr atom_entry =
		(send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
		if ((atom_entry != NULL) && (atom_entry->atom != atom)) {
		    printf("Warning:  Atom use inconsistency.\n");
		    printf("\tThis program associates the string \"%s\" with atom value %d, %x, '%c%c%c%c'\n",
			   str, atom_entry->atom, atom_entry->atom,
			   ((char*)&atom_entry->atom)[0],
			   ((char*)&atom_entry->atom)[1],
			   ((char*)&atom_entry->atom)[2],
			   ((char*)&atom_entry->atom)[3]);
		    printf("\tOther programs use the atom value %d, %x, '%c%c%c%c'\n", atom,
			   atom, ((char*)&atom)[0], ((char*)&atom)[1], 
			((char*)&atom)[2], ((char*)&atom)[3]
			);
		}
	    }
	    entry = Tcl_FindHashEntry(&as->value_hash_table, (char *) (int64_t)atom);
	    if (entry != NULL) {
		send_get_atom_msg_ptr atom_entry =
		(send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
		if ((atom_entry != NULL) &&
		    (strcmp(atom_entry->atom_string, str) != 0)) {
		    printf("Warning:  Atom use inconsistency.\n");
		    printf("\tThis program associates the string \"%s\" with atom value %d, %x, '%c%c%c%c'\n",
			   atom_entry->atom_string, atom_entry->atom,
			   atom_entry->atom, 
			   ((char*)&atom_entry->atom)[0],
			   ((char*)&atom_entry->atom)[1],
			   ((char*)&atom_entry->atom)[2],
			   ((char*)&atom_entry->atom)[3]);
		    printf("\tOther programs associate the string \"%s\" with that value\n", str);
		}
		printf("Atom cache inconsistency, tried to associate value %d %x, '%c%c%c%c' with string \"%s\"\n	Previous association was string \"%s\"\n",
		       atom, atom, ((char*)&atom)[0], ((char*)&atom)[1], 
		       ((char*)&atom)[2], ((char*)&atom)[3], str, atom_entry->atom_string);
	    }
	    break;
	}
    default:
	printf("Warning: Got an unexpected message \"%s\"\n", msg);
    }
}

static
int
enter_atom_into_cache(atom_server as, send_get_atom_msg_ptr msg)
{
    int new;
    char *str;
    send_get_atom_msg_ptr stored;
    Tcl_HashEntry *entry = NULL;

    if ((msg->atom_string == NULL) || (msg->atom == -1))
	return 0;
    str = strdup(msg->atom_string);
    stored = (send_get_atom_msg_ptr) malloc(sizeof(send_get_atom_msg));
    stored->atom_string = str;
    stored->atom = msg->atom;

    /* enter into string hash table */
    entry = Tcl_CreateHashEntry(&as->string_hash_table, str, &new);
    if (!new) {
	/* already inserted by someone else */
	free(stored);
	free(str);
	return 0;
    }
    Tcl_SetHashValue(entry, stored);
    /* enter into value hash table */
    entry = Tcl_CreateHashEntry(&as->value_hash_table,
				(char *) (int64_t) stored->atom, &new);
    if (!new) {
	printf("Serious internal error in atom cache.  Duplicate value hash entry.\n");
	exit(1);
    }
    Tcl_SetHashValue(entry, stored);
    return 1;
}

void
set_string_and_atom(atom_server as, char *str, atom_t atom)
{
    send_get_atom_msg tmp_value;
    Tcl_HashEntry *entry = NULL, *entry2 = NULL;
    long numbytes, len;
    unsigned char buf[MAXDATASIZE];
    socklen_t addr_len = sizeof(struct sockaddr);
    int new;

    entry = Tcl_FindHashEntry(&as->string_hash_table, str);
    if (entry != NULL) {
	send_get_atom_msg_ptr atom_entry =
	(send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
	if ((atom_entry != NULL) && (atom_entry->atom != atom)) {
	    printf("Atom cache inconsistency, tried to associate string \"%s\" with value %d, %x, '%c%c%c%c'\n	Previous association was value %d, %x, '%c%c%c%c'\n", 

		   str, atom, atom, ((char*)&atom)[0], ((char*)&atom)[1], 
		   ((char*)&atom)[2], ((char*)&atom)[3], atom_entry->atom,
		   atom_entry->atom, ((char*)&atom_entry->atom)[0],
		   ((char*)&atom_entry->atom)[1],
		   ((char*)&atom_entry->atom)[2],
		   ((char*)&atom_entry->atom)[3]);
	    return;
	}
    }
    entry2 = Tcl_FindHashEntry(&as->value_hash_table, (char *) (int64_t) atom);
   if (entry2 != NULL) {
	send_get_atom_msg_ptr atom_entry =
	(send_get_atom_msg_ptr) Tcl_GetHashValue(entry2);
	if ((atom_entry != NULL) &&
	    (strcmp(atom_entry->atom_string, str) != 0)) {
	    printf("Atom cache inconsistency, tried to associate value %d, %x, '%c%c%c%c' with string \"%s\"\n	Previous association was string \"%s\"\n",
		   atom, atom, ((char*)&atom)[0], ((char*)&atom)[1], 
		   ((char*)&atom)[2], ((char*)&atom)[3], str, 
		   atom_entry->atom_string);
	    return;
	}
    }
    tmp_value.atom = atom;
    tmp_value.atom_string = str;
    new = enter_atom_into_cache(as, &tmp_value);
    if (as->no_server) return;
    if (!new) return;
    sprintf((char *)&buf[1], "A%d %s", atom, str);
    len = (long) strlen((char*)&buf[1]);
    if (as->use_tcp) {
	set_blocking(as, 1);
	buf[0] = (unsigned char) len;
	if (establish_server_connection(as, 1) == 0) return;
	
	if ((numbytes = write(as->tcp_fd, (char*)buf, len+1)) != len +1) {
	    close(as->tcp_fd);
	    return;
	}
	set_blocking(as, 0);
	if (read(as->tcp_fd, (char*)buf, 1) != 1) {
	    return;
	}
	if (read(as->tcp_fd, (char*)&buf[1], buf[0]) != buf[0]) {
	    return;
	}
	buf[buf[0]+1] = 0;
	handle_unexpected_msg(as, (char*) &buf[1]);
    } else {
	if (as->their_addr.sin_addr.s_addr == 0) return;
	set_blocking(as, 0);	/* set server fd nonblocking */
	if ((numbytes = sendto(as->sockfd, (char*)&buf[1], len, 0,
			       (struct sockaddr *) &(as->their_addr), sizeof(struct sockaddr))) == -1) {
	    /* don't try that again... */
	    as->their_addr.sin_addr.s_addr = 0;
	    return;
	}
	if ((numbytes = recvfrom(as->sockfd, (char*)&buf[1], MAXDATASIZE - 1, 0,
				 (struct sockaddr *) &(as->their_addr), &addr_len)) != -1) {
	    /* actually got a message back ! */
	    buf[numbytes+1] = 0;
	    handle_unexpected_msg(as, (char*) &buf[1]);
	}
    }
}

static int atom_server_verbose = 1;

static int
fill_hostaddr(void *addr, char *hostname)
{
    struct hostent *host_addr;
    
    host_addr = gethostbyname(hostname);
    if (host_addr == NULL) {
	int address;
	address = inet_addr(hostname);
	if (address == -1) {
	    /* 
	     *  not translatable as a hostname or 
	     * as a dot-style string IP address
	     */
	    return 0;
	}
	assert(sizeof(int) == sizeof(struct in_addr));
	*((int*)addr) = (int)address;
    } else {
	memcpy(addr, host_addr->h_addr, host_addr->h_length);
    }
    return 1;
}

static int
establish_server_connection(atom_server as, int do_fallback)
{
    SOCKET sock;
    int delay_value = 1;
    char ping_char = 0;

    if (atom_server_verbose == -1) {
	if (getenv("ATOM_SERVER_VERBOSE") == NULL) {
	    atom_server_verbose = 0;
	} else {
	    atom_server_verbose = 1;
	}
    }
    if (as->tcp_fd == -2) return 0;
    if ((as->tcp_fd == -1) || 
	(write(as->tcp_fd, &ping_char, 1) != 1)) {
	/* reestablish connection, name_str is the machine name */
	struct sockaddr_in sock_addr;

	fprintf(stderr, "Establish server connection, write failed, creating socket\n");
	if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
	    fprintf(stderr, "Failed to create socket for ATL atom server connection.  Not enough File Descriptors?\n");
	    return 0;
	}
	
	sock_addr.sin_family = AF_INET;
		
	fprintf(stderr, "Establish server connection, fill_host_addr\n");
	if (fill_hostaddr(&sock_addr.sin_addr, atom_server_host) == 0) {
	    fprintf(stderr, "Unknown Host \"%s\" specified as ATL atom server.\n",
		    atom_server_host);
	    as->tcp_fd = -2;
	    return 0;
	}
	sock_addr.sin_port = htons(TCP_PORT);

	if (atom_server_verbose) {
	    printf("Trying connection to atom server on %s ...  ",
		   atom_server_host);
	}
	if (connect(sock, (struct sockaddr *) &sock_addr,
		    sizeof sock_addr) < 0) {

	    if (atom_server_verbose) {
		printf("failed\n");
	    }
	    if (!do_fallback) return 0;

	    /* fallback */
	    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
		fprintf(stderr, "Failed to create socket for ATL atom server connection.  Not enough File Descriptors?\n");
	       return 0;
	    }
	    atom_server_host = "atomhost.cercs.gatech.edu";
	    sock_addr.sin_family = AF_INET;
	    if (fill_hostaddr(&sock_addr.sin_addr, atom_server_host) == 0){
		fprintf(stderr, "Unknown Host \"%s\" specified as ATL atom server.\n",
			atom_server_host);
		as->tcp_fd = -2;
		return 0;
	    }
	    sock_addr.sin_port = htons(TCP_PORT);
	    if (atom_server_verbose) {
		printf("Trying fallback connection to atom server on %s ...  ",
		       atom_server_host);
	    }
	    if (connect(sock, (struct sockaddr *) &sock_addr,
			sizeof sock_addr) < 0) {
		fprintf(stderr, "Failed to connect to primary or fallback atom servers.\n");
		as->tcp_fd = -2;
	        return 0;
	    }
	}
	if (atom_server_verbose) {
	    printf("succeeded\n");
	}
	setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
		   sizeof(delay_value));
	as->tcp_fd = sock;
	/* 
	 * ignore SIGPIPE's  (these pop up when ports die.  we catch the 
	 * failed writes) 
	 */
#ifdef SIGPIPE
	signal(SIGPIPE, SIG_IGN);
#endif
    }
    return 1;
}

extern atom_t ATLget_hash(const char *str);

extern
 atom_t
atom_from_string(atom_server as, char *str)
{
    atom_t atom;

    atom = ATLget_hash(str);

    set_string_and_atom(as, str, atom);

    return atom;
}

extern
char *
string_from_atom(atom_server as, atom_t atom)
{
    send_get_atom_msg tmp_rec;
    send_get_atom_msg_ptr stored;
    Tcl_HashEntry *entry = NULL;
    int numbytes;
    char buf[MAXDATASIZE];

    entry = Tcl_FindHashEntry(&as->value_hash_table, (char *) (int64_t) atom);

    if (entry == NULL) {
	sprintf(&buf[1], "N%d", atom);
	if (establish_server_connection(as, 1) == 0) return NULL;
	buf[0] = (char) strlen(&buf[1]);
	if (write(as->tcp_fd, buf, buf[0]+1) != buf[0] + 1) {
	    perror("write");
	    return NULL;
	}
	set_blocking(as, 1);	/* set server fd blocking */
	buf[1] = 0;
	while (buf[1] != 'S') {
	    if ((numbytes = read(as->tcp_fd, buf, 1)) == -1) {
		perror("read");
		return NULL;
	    }
	    if ((numbytes = read(as->tcp_fd, &buf[1], buf[0])) != buf[0]) {
		perror("read2");
		return NULL;
	    }
	    buf[numbytes+1] = 0;
	    if (buf[1] != 'S')
		handle_unexpected_msg(as, &buf[1]);
	}

	if (buf[2] == 0) {
	    return NULL;
	}
	tmp_rec.atom_string = &buf[2];
	tmp_rec.atom = atom;

	(void) enter_atom_into_cache(as, &tmp_rec);
	stored = &tmp_rec;
    } else {
	stored = (send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
    }
    if (stored->atom_string != NULL) {
	return strdup(stored->atom_string);
    } else {
	return NULL;
    }
}

extern
char *
get_server_id(atom_server as)
{
    return as->server_id;
}

#ifdef HAVE_WINDOWS_H
/* Winsock init stuff, ask for ver 1.1 */
static WORD wVersionRequested = MAKEWORD(1, 1);
static WSADATA wsaData;

static void
nt_socket_init_func()
{
    int nErrorStatus;
    nErrorStatus = WSAStartup(wVersionRequested, &wsaData);
    if (nErrorStatus != 0) {
	fprintf(stderr, "Could not initialize windows socket library!");
	WSACleanup();
	exit(-1);
    }
}
#else
static void nt_socket_init_func(){}
#endif


static char *in_use_values[] = {
"CM_BW_MEASURED_COF",
"CM_BW_MEASURED_VALUE",
"CM_BW_MEASURE_INTERVAL",
"CM_BW_MEASURE_SIZE",
"CM_BW_MEASURE_SIZEINC",
"CM_BW_MEASURE_TASK",
"CM_CONN_BLOCKING",
"CM_ENET_ADDR",
"CM_ENET_HOST",
"CM_ENET_PORT",
"CM_EVENT_SIZE",
"CM_FD",
"CM_INCOMING_CONNECTION",
"CM_NETWORK_POSTFIX",
"CM_NNTI_TRANSPORT",
"CM_REBWM_REPT",
"CM_REBWM_RLEN",
"CM_REG_BW_REPEAT_CNT",
"CM_REG_BW_RUN_LEN",
"CM_SHM_MAX_PAYLOAD",
"CM_SHM_NUM_SLOTS",
"CM_TRANSPORT",
"CM_TRANSPORT_RELIABLE",
"CM_TRANS_MEGABITS_SEC",
"CM_TRANS_TEST_DURATION_SECS",
"CM_TRANS_TEST_RECEIVED_COUNT",
"CM_TRANS_TEST_REPEAT",
"CM_TRANS_TEST_REUSE_WRITE_BUFFER",
"CM_TRANS_TEST_SIZE",
"CM_TRANS_TEST_TAKEN_CORRUPT",
"CM_TRANS_TEST_TAKE_RECEIVE_BUFFER",
"CM_TRANS_TEST_VECS",
"CM_TRANS_TEST_VERBOSE",
"CMdemo_test_atom",
"CONNECTION_FILE_DESCRIPTOR",
"DoReconfig",
"ECHO_EVENT_NETWORK",
"ECHO_EVENT_TRANSPORT",
"ECHO_USE_EVENT_TRANSPORT",
"ECho_attr_test_atom",
"EV_BACKPRESSURE_HIGH",
"EV_BACKPRESSURE_LOW",
"EV_EVENT_COUNT",
"EV_EVENT_LSUM",
"EventCount",
"IP_ADDR",
"IP_HOST",
"IP_PORT",
"MCAST_ADDR",
"MCAST_PORT",
"NNTI_ADDR",
"NNTI_ENET_CONTROL",
"NNTI_IMMEDIATE_PULL_WAIT",
"NNTI_PARAMS",
"NNTI_PORT",
"NNTI_SHM",
"PEER_CONN_PORT",
"PEER_HOSTNAME",
"PEER_IP",
"PEER_LISTEN_PORT",
"SSL_PORT",
"THIS_CONN_PORT",
"UDP_ADDR",
"UDP_PORT",
"application_reconfiguration_atom",
"fp_dst_condition",
"fp_dst_rank",
"fp_flush_id",
"fp_size",
"fp_starttime",
"hop_count_atom",
"index_atom",
"iteration",
"level",
"mpisize",
"test_value",
NULL};

static void
preload_in_use_atoms(atom_server as)
{
    int i=0;
    while (in_use_values[i] != NULL) {
	(void) atom_from_string(as, in_use_values[i++]);
    }
}

void
free_atom_server(atom_server as)
{
  Tcl_HashSearch search;
  Tcl_HashEntry * entry = Tcl_FirstHashEntry(&as->string_hash_table, &search);
  while (entry) {
    send_get_atom_msg_ptr stored;
    stored = (send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
    free(stored->atom_string);
    free(stored);
    entry = Tcl_NextHashEntry(&search);
  }
  Tcl_DeleteHashTable(&as->string_hash_table);
  Tcl_DeleteHashTable(&as->value_hash_table);
  free(as);
}

atom_server
init_atom_server(atom_cache_type cache_style)
{
    atom_server as = (atom_server) malloc(sizeof(atom_server_struct));

    nt_socket_init_func();
    if (atom_server_host == NULL) {	/* environment override */
	atom_server_host = getenv("ATOM_SERVER_HOST");
    }
    if (atom_server_host == NULL) {
	atom_server_host = ATOM_SERVER_HOST;	/* from configure */
    }
    as->server_id = atom_server_host;
    as->tcp_fd = -1;
    as->use_tcp = (getenv("ATL_USE_TCP") != NULL);
    as->no_server = 1;

    Tcl_InitHashTable(&as->string_hash_table, TCL_STRING_KEYS);
    Tcl_InitHashTable(&as->value_hash_table, TCL_ONE_WORD_KEYS);

    if ((as->he = gethostbyname(atom_server_host)) == NULL) {
	as->he = NULL;
	as->their_addr.sin_addr.s_addr = 0;
    } else {
	as->their_addr.sin_addr = *((struct in_addr *) as->he->h_addr);
    }
    if ((as->sockfd = (int) socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
	perror("socket");
	exit(1);
    }
#ifndef HAVE_WINDOWS_H
    as->flags = fcntl(as->sockfd, F_GETFL);
#else
    as->flags = 0;
#endif
    as->their_addr.sin_family = AF_INET;
    as->their_addr.sin_port = htons(UDP_PORT);
    memset(&(as->their_addr.sin_zero), '\0', 8);

    preload_in_use_atoms(as);
    as->no_server = 0;
    return as;
}
