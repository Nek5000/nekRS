#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <errno.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#ifndef HAVE_WINDOWS_H
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <sys/time.h>
#else
#include <windows.h>
#include <winsock.h>
#ifndef EWOULDBLOCK
#define EWOULDBLOCK WSAEWOULDBLOCK
#endif
#endif
#include <fcntl.h>
#include <tclHash.h>
typedef int atom_t;
#include <stdint.h>
#include "atom_internal.h"
#ifdef _MSC_VER
    #define strdup _strdup
    #include <io.h>
#pragma warning(disable: 4996)
#endif

#define MAXBUFLEN 100

extern Tcl_HashTable *stringhash;
extern Tcl_HashTable *valuehash;
extern char *atom_to_string(Tcl_HashTable *, Tcl_HashTable *, int);
extern int string_to_atom(Tcl_HashTable *, Tcl_HashTable *, char *);
extern int set_string_and_atom(Tcl_HashTable *, Tcl_HashTable *, char *, int);
static void close_client(int client);
static void server_init(int udp_socket, int tcp_socket);
static int poll_and_handle();


Tcl_HashTable *stringhash;
Tcl_HashTable *valuehash;

int verbose = 0;

#define LOG printf

char *
atom_to_string(Tcl_HashTable * string_hash_table, Tcl_HashTable * value_hash_table, int value)
{
    Tcl_HashEntry *entry = NULL;
    send_get_atom_msg_ptr value_string;

    if (verbose)
	printf("Doing a atom_to_string\n");

    entry = Tcl_FindHashEntry(value_hash_table, (char *) (int64_t)value);

    if (entry) {
	value_string = (send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
	return value_string->atom_string;
    }
    return 0;

}

int
string_to_atom(Tcl_HashTable * string_hash_table, Tcl_HashTable * value_hash_table, char *a)
{
    Tcl_HashEntry *entry = NULL;
    send_get_atom_msg_ptr return_msg;

    if (verbose)
	printf("Doing a string_to_atom\n");

    entry = Tcl_FindHashEntry(string_hash_table, a);

    if (entry) {
	return_msg = (send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
	return return_msg->atom;
    }
    return -1;	/* the string was not in the db */

}

int
set_string_and_atom(Tcl_HashTable * string_hash_table, Tcl_HashTable * value_hash_table, char *a, int set_atom)
{

    send_get_atom_msg_ptr return_msg;
    Tcl_HashEntry *entry = NULL;
    int new;
    send_get_atom_msg_ptr stored;

    stored = (send_get_atom_msg_ptr) malloc(sizeof(send_get_atom_msg));
    stored->atom_string = strdup(a);

    stored->atom = set_atom;
    entry = Tcl_CreateHashEntry(string_hash_table, a, &new);
    Tcl_SetHashValue(entry, stored);
    entry = Tcl_CreateHashEntry(value_hash_table, (char *) (int64_t)stored->atom,
				&new);
    Tcl_SetHashValue(entry, stored);
    return_msg = (send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
    return return_msg->atom;	/* Success */

}

void
Initialize(void)
{
    stringhash = (Tcl_HashTable *) malloc(sizeof(Tcl_HashTable));
    Tcl_InitHashTable(stringhash, TCL_STRING_KEYS);
    valuehash = (Tcl_HashTable *) malloc(sizeof(Tcl_HashTable));
    Tcl_InitHashTable(valuehash, TCL_ONE_WORD_KEYS);
}


extern int
establish_server_connection()
{
    int sock;
#ifdef MODULE
    struct socket *socket;	
#endif
    struct sockaddr_in sock_addr;


    if ((sock = (int) socket(AF_INET, SOCK_STREAM, 0)) < 0) {
	fprintf(stderr, "Failed to create socket for ATL atom server connection.  Not enough File Descriptors?\n");
	return 0;
    }
	
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_port = htons(TCP_PORT);
    sock_addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    
    if (connect(sock, (struct sockaddr *) &sock_addr, sizeof sock_addr) < 0) {
	return 0;
    }
    close(sock);
    return 1;
}

int
main(int argc, char **argv)
{
    int sockfd, tcpfd;
    int quiet = 0;
    int run = 1;
    int no_fork = 0;
    struct sockaddr_in my_addr;	/* my address information */
    struct sockaddr_in sock_addr;	/* connector's address information */
    int sock_opt_val = 1;

    int i;

    for (i = 1; i < argc; i++) {
	if (strcmp(argv[i], "-no_fork") == 0) {
	    no_fork++;
	} else if (strcmp(argv[i], "-quiet") == 0) {
	    quiet = 1;
	    verbose = 0;
	} else if (strcmp(argv[i], "-verbose") == 0) {
	    quiet = 0;
	    verbose = 1;
	} else {
	    fprintf(stderr, "Unknown argument \"%s\"\n", argv[i]);
	    exit(1);
	}
    }

    if ((sockfd = (int) socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
	perror("socket");
	exit(1);
    }
    my_addr.sin_family = AF_INET;	/* host byte order */
    my_addr.sin_port = htons(UDP_PORT);	/* short, network byte order */
    my_addr.sin_addr.s_addr = INADDR_ANY; /* automatically fill with my IP */
    memset(&(my_addr.sin_zero), '\0', 8);	/* zero the rest of the *
						 * struct */
    {
	if (establish_server_connection()) {
	    if (!quiet)
		printf("\n\tAtom server already running\n");
	    exit(0);

	}
    }

    if (!quiet)
	printf("\n\tNo Atom server found - ");
#ifdef HAVE_FORK
    if (!no_fork) {
	if (!quiet)
	    printf("Forking server to background\n");
	if (fork() != 0) {
	    /* I'm the parent, return now */
	    exit(0);
	}
    } else {
	if (!quiet)
	    printf("Running...\n");
    }
#endif

    if (bind(sockfd, (struct sockaddr *) &my_addr,
	     sizeof(struct sockaddr)) == -1) {
	perror("bind");
	exit(1);
    }

    tcpfd = (int) socket(AF_INET, SOCK_STREAM, 0);
    if (tcpfd < 0) {
	fprintf(stderr, "Cannot open INET socket %d\n", tcpfd);
	return -1;
    }
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = htons((unsigned short) TCP_PORT);
    sock_opt_val = 1;
    if (setsockopt((long) tcpfd, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	fprintf(stderr, "Failed to set 1REUSEADDR on INET socket\n");
	return -1;
    }
    if (bind(tcpfd, (struct sockaddr *) &sock_addr,
	     sizeof sock_addr) < 0) {
	fprintf(stderr, "Cannot bind INET socket\n");
	return -1;
    }
    sock_opt_val = 1;
    if (setsockopt(tcpfd, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	perror("Failed to set 2REUSEADDR on INET socket");
	return -1;
    }
    /* begin listening for conns and set the backlog */
    if (listen(tcpfd, FD_SETSIZE)) {
	fprintf(stderr, "listen failed\n");
	return -1;
    }
    run = 1;
    Initialize();
    server_init(sockfd, tcpfd);

    while (run) {
	poll_and_handle();
    }
    /* NOTREACHED */
    close(sockfd);
    return 0;
}

static void
process_data(char* buf, char *response)
{
    switch (buf[0]) {
    case 'P':
	/* Ping */
	if (verbose)
	    printf("Sending %c\n", 'R');
	response[0] = 'R';
	response[1] = 0;
	return;
    case 'N':{
	/* request translation of numeric value to a string */
	char *str;
	int atom = 0;
	if (sscanf(&buf[1], "%d", &atom) != 1)
	    break;
	str = atom_to_string(stringhash, valuehash, atom);
	response[0] = 'S';
	if (str) {
	    strncpy(&response[1], str, MAXBUFLEN-1);
	} else {
	    response[1] = 0;
	}
	if (verbose)
	    printf("Sending %s\n", response);
	return;
    }
    case 'S': {
	/* request translation of string to a numeric value */
	int value = string_to_atom(stringhash, valuehash, &buf[1]);

	sprintf(response, "N%d", value);
	if (verbose)
	    printf("Sending %s\n", response);
	return;
    }
    case 'A':{
	/* create an association between a string and a value */
	int atom = 0;
	char *str;
	Tcl_HashEntry *entry = NULL;

	atom = strtol(&buf[1], &str, 10);
	str++;
	entry = Tcl_FindHashEntry(stringhash, str);
	if (entry != NULL) {
	    send_get_atom_msg_ptr atom_entry =
		(send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
	    if ((atom_entry != NULL) && (atom_entry->atom != atom)) {
		if (verbose)
		    printf("Atom cache inconsistency, tried to associate string \"%s\" with value %d\n	Previous association was value %d\n",
			   str, atom, atom_entry->atom);
		sprintf(response, "E%d %s", atom_entry->atom,
			atom_entry->atom_string);
		if (verbose)
		    printf("Sending %s\n", response);
		
		return;
	    }
	}
	entry = Tcl_FindHashEntry(valuehash, (char *) (int64_t) atom);
	if (entry != NULL) {
	    send_get_atom_msg_ptr atom_entry =
		(send_get_atom_msg_ptr) Tcl_GetHashValue(entry);
	    if ((atom_entry != NULL) &&
		(strcmp(atom_entry->atom_string, str) != 0)) {
		if (verbose)
		    printf("Atom cache inconsistency, tried to associate value %d with string \"%s\"\n	Previous association was string \"%s\"\n",
			   atom, str, atom_entry->atom_string);
		sprintf(response, "E%d %s", atom_entry->atom,
			atom_entry->atom_string);
		if (verbose)
		    printf("Sending %s\n", response);

		return;
	    }
	}
	set_string_and_atom(stringhash, valuehash, str, atom);
    }
    }
}

static void
handle_udp_data(int sockfd)
{
    int numbytes;
    socklen_t addr_len = sizeof(struct sockaddr);
    char buf[MAXBUFLEN];
    char response[MAXBUFLEN];
    struct sockaddr_in their_addr;	/* connector's address information */
    if ((numbytes = recvfrom(sockfd, buf, MAXBUFLEN - 1, 0,
			     (struct sockaddr *) &their_addr, &addr_len)) == -1) {
	perror("recvfrom");
	exit(1);
    }
    buf[numbytes] = 0;

    if (verbose)
	printf("UDP recd. %s\n", buf);

    response[0] = 0;
    process_data(buf, response);

    if (response[0] != 0) {
	int len = (int)strlen(response);
	if ((numbytes = sendto(sockfd, &response[0], len, 0,
			       (struct sockaddr *) &their_addr,
			       sizeof(struct sockaddr))) == -1) {
	    perror("sendto");
	    exit(1);
	}
    }
}

static void
handle_tcp_data(int sockfd)
{
    int numbytes;
    unsigned char len;
    char buf[MAXBUFLEN];
    char response[MAXBUFLEN];
    if ((numbytes = read(sockfd, &len, 1)) != 1) {
	if ((numbytes == -1) && (errno != EBADF)){
	    perror("read");
	}
	close_client(sockfd);
    }
    if ((numbytes = read(sockfd, &buf, len)) != len) {
	if ((numbytes == -1) && (errno != EBADF)){
	    perror("read");
	}
	close_client(sockfd);
    }
    buf[len] = 0;

    if (verbose)
	printf("TCP recd. %s\n", buf);

    response[1] = 0;
    process_data(buf, &response[1]);

    if (response[1] != 0) {
	unsigned len = (unsigned) strlen(&response[1]);
	response[0] = len;
	if ((numbytes = write(sockfd, &response[0], len +1) != len + 1)) {
	    perror("write - handle_tcp_data");
	    close_client(sockfd);
	}
    }
}

typedef struct _ASClient {
    time_t created;
    time_t timestamp;
} *ASClient;

#define Max(i,j) ((i<j) ? j : i)

static fd_set server_fdset;
static ASClient clients = NULL;
static int conn_sock_inet = -1;
static int udp_sock = -1;

static void
close_client(int client)
{
    if (verbose) LOG("Closing client %d\n", client);
    if (client == conn_sock_inet) return;
    if (client == udp_sock) return;
    clients[client].created = -1;
    close(client);
    FD_CLR(client, &server_fdset);
}

static void
server_init(int udp_socket, int tcp_socket)
{
    int max_fd = FD_SETSIZE;	/* returns process fd table size */
    int i;
    FD_ZERO(&server_fdset);
    udp_sock = udp_socket;
    conn_sock_inet = tcp_socket;
    FD_SET(udp_socket, &server_fdset);
    FD_SET(tcp_socket, &server_fdset);
    clients = (ASClient)malloc(sizeof(struct _ASClient)*max_fd);
    for(i=0; i < max_fd; i++) {
	clients[i].created = -1;
    }
}

#define CONN_TIMEOUT_INTERVAL 3600

static void
timeout_old_conns()
{
    static time_t last_timeout_time = 0;
    time_t now = time(NULL);
    int i;

    if (last_timeout_time == 0) {
	last_timeout_time = time(NULL);
	return;
    }
    if (now - last_timeout_time < CONN_TIMEOUT_INTERVAL) {
	return;
    }
    last_timeout_time = now;
    for (i = 0; i < FD_SETSIZE; i++) {
	if (clients[i].created != -1) {
	    if (now - clients[i].timestamp > CONN_TIMEOUT_INTERVAL) {
		LOG("Timeout -> fd %d", i);
		close_client(i);
	    }
	}
    }
}

static void
timeout_oldest_conn()
{
    time_t oldest_time = 0;
    int oldest_index = 0;
    int i;

    for (i = 0; i < FD_SETSIZE; i++) {
	if (clients[i].created != -1) {
	    if (oldest_time == 0) {
		oldest_time = clients[i].timestamp;
		oldest_index = i;
	    } else if (oldest_time > clients[i].timestamp) {
		oldest_time = clients[i].timestamp;
		oldest_index = i;
	    }
	}
    }
    if (oldest_time != 0) {
	LOG("Closing oldest -> fd %d", oldest_index);
	close_client(oldest_index);
    }
}

extern void
accept_conn_sock(int conn_sock)
{
    int sock;
    int delay_value = 1;

    struct linger linger_val;

    linger_val.l_onoff = 1;
    linger_val.l_linger = 60;

    if ((sock = (int) accept(conn_sock, (struct sockaddr *) 0, (socklen_t *) 0)) < 0) {
	LOG("accept failed");

	timeout_oldest_conn();

	if ((sock = (int) accept(conn_sock, (struct sockaddr *) 0, (socklen_t *) 0)) < 0) {
	    LOG("accept failed twice");
	    fprintf(stderr, "Cannot accept socket connection\n");
	    return;
	}
    }
    if (verbose) LOG("\nAccepting fd %d\n", sock);
    if (setsockopt(sock, SOL_SOCKET, SO_LINGER, (char *) &linger_val,
		   sizeof(struct linger)) != 0) {
	perror("set SO_LINGER");
	return;
    }
    setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));

    clients[sock].timestamp = clients[sock].created = time(NULL);

    FD_SET(sock, &server_fdset);
    if (sock > FD_SETSIZE) {
	fprintf(stderr, "Internal Error, stupid WINSOCK large FD bug.\n");
	fprintf(stderr, "Increase FD_SETSIZE.  Item not added to fdset.\n");
    }

/*    out_domain_rejection(sock);*/
}

static int
poll_and_handle()
{
    int i, res;
    fd_set rd_set = server_fdset;
    struct timeval timeout;

#if !defined(__FD_SET) && !defined(HAVE_WINDOWS_H) && 0
    LOG("  selecting, rd_set is %lx,%lx,%lx,%lx,%lx,%lx", rd_set.fds_bits[0],
	rd_set.fds_bits[1], rd_set.fds_bits[2], rd_set.fds_bits[3],
	rd_set.fds_bits[4], rd_set.fds_bits[5]);
#endif
    /* scan existing conn for incoming data & do something */
    timeout.tv_usec = 0;
    timeout.tv_sec = CONN_TIMEOUT_INTERVAL + 5;

    res = select(FD_SETSIZE, &rd_set, (fd_set *) NULL,
		 (fd_set *) NULL, &timeout);
    if (res == -1) {
	if (errno == EBADF) {
	    int i;
	    int found_one = 0;
	    for (i = 0; i < FD_SETSIZE; i++) {
		if (FD_ISSET(i, &rd_set) && (i != (long) conn_sock_inet)) {
		    fd_set test_set;
		    timeout.tv_usec = 0;
		    timeout.tv_sec = 0;
		    FD_ZERO(&test_set);
		    FD_SET(i, &test_set);
		    errno = 0;
		    select(FD_SETSIZE, &test_set, (fd_set *) NULL,
			   (fd_set *) NULL, &timeout);
		    if (errno == EBADF) {
			fprintf(stderr, "Select failed, fd %d\n is bad.  Removing from select list.\n",
				i);
			LOG("REmoving bad fd %d", i);
			clients[i].created = -1;
		    } else if ((i !=  conn_sock_inet) && 
			       (i !=  udp_sock)) {
			close_client(i);
		    }
		    found_one++;
		}
	    }
	    if (found_one == 0) {
		fprintf(stderr, "Bad file descriptor in select.  Failed to localize.  Exiting.\n");
		exit(1);
	    }
	} else if (errno == EINTR) {
	    return 0;
	}
	fprintf(stderr, "select failed, errno %d\n", errno);
	return -1;
    } else if (res == 0) {
	timeout_old_conns();
	return 0;
    }
    if (conn_sock_inet >= 0 && FD_ISSET(conn_sock_inet, &rd_set)) {
	FD_CLR(conn_sock_inet, &rd_set);
	accept_conn_sock(conn_sock_inet);
    }
    if (udp_sock >= 0 && FD_ISSET(udp_sock, &rd_set)) {
	FD_CLR(udp_sock, &rd_set);
	handle_udp_data(udp_sock);
    }
    for (i = 0; i < FD_SETSIZE; i++) {	/* check each possible read fd */
	if ((clients[i].created == -1) || (!FD_ISSET(i, &rd_set))) continue;
	clients[i].timestamp = time(NULL);
	handle_tcp_data(i);

    }
    timeout_old_conns();
    return res;
}

