#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#undef FD_SETSIZE
#define FD_SETSIZE 1024
#include <winsock.h>
#define getpid()	_getpid()
#define getcwd(a,b)	_getcwd(a,b)
#define __ANSI_CPP__
#else
#include <unistd.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#include <sys/socket.h>
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#endif
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#include <signal.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <signal.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif

#include "fm.h"
#include "fm_internal.h"
#include "io_interface.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include <assert.h>
#include <pthread.h>

typedef struct IOFormatRep {
    server_ID_type server_ID;
    format_rep server_format_rep;
} *IOFormatRep;

typedef struct format_list {
    IOFormatRep format;
    time_t last_reference;
    struct format_list *next;
} format_list;

#define MAX_CLIENT_FD 256

typedef struct _FSclient {
    int port;			/* port's port number */
    char *hostname;		/* port's host name */
    char *usock_name;		/* port's unix socket name */
    void *fd;
    int byte_reversal;
    format_server fs;
    int provisional;
    int version;
    long created;
    int input_bytes;
    int output_bytes;
    int formats_registered;
    int formats_fetched;
    unsigned char *key;
    int key_len;
    pthread_t handler_thread;
} *FSClient;

typedef struct _format_server {
    int port;			/* server port as announced to world */
    char *hostname;		/* server's hostname */
    int format_server_identifier; /* pseudo-random generated id */
    time_t start_time;    
    int stdout_verbose;
    void *conn_sock_inet;	/* public sock for IP socket connection */
    void *conn_sock2_inet;	/* 2nd public sock for IP connection */
    char *usock_name;		/* file name for unix socket connection */
    int portCount;		/* number of current ports */
    FSClient *ports;		/* array of ports index by socket fd */
    time_t *timestamp;		/* array of timestamps */
    fd_set fdset;		/* bitmap of those conns */
    char *data_buffer;		/* buffer for data */
    int buffer_size;
    FMContext proxy_context_to_master;
    int pending_unregistered;
    format_list *lists[1];	/* first dimension float format */
    int fork_threads;
    pthread_mutex_t lock;
    pthread_mutex_t log_lock;
} _fsserver;

static FSClient
 format_server_accept_conn_sock(format_server fs, void *conn_sock);
static void read_formats_from_file(format_server fs);
static format_server format_server_create();
static void format_server_handle_data(format_server fs, FSClient fsc);
static int format_server_listen(format_server fs, int port_num);
static int format_server_poll_and_handle(format_server fs);
static void timeout_old_conns(format_server server_fs);
static void timeout_oldest_conn(format_server server_fs);
static void send_stats(FSClient fsc);
static FMContext stats_context = NULL;
static IOFormatRep fetch_format(format_server fs, FSClient fsc, 
				      const unsigned char *format_id,
				      int format_id_length);
static int try_master_connect(format_server fs);

#ifndef SELECT_DEFINED
extern int select(int width, fd_set * readfds, fd_set * writefds,
			fd_set * exceptfds, struct timeval * timeout);
#endif
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <time.h>
extern pid_t getpid();
#ifndef HAVE_GETDOMAINNAME
extern int getdomainname(char *name, size_t len);
#endif
#ifndef PATH_MAX
#define PATH_MAX _MAX_PATH
#endif
static char format_server_log[PATH_MAX + 1];
static char FS_restart_file[PATH_MAX + 1];

static void dump_stats_to_log(format_server fs);

static FILE *log = (void *) 0;
#define _GNU_SOURCE
#include <unistd.h>

static int no_log = 0;
void
ffs_set_no_log(int in_no_log)
{
    no_log = in_no_log;
}

extern void
LOG(format_server fs, char *format,...)
{
    va_list ap;
    struct stat stat_buf;
    static int log_count = 0;
    time_t rawtime;
    struct tm *info;
    char buffer[80];

    if (no_log) return;
    pthread_mutex_lock(&fs->log_lock);
    if ((log == (void *) -1) || (log == (void *) 0)) {
	int restart = 0;
	if (log == (void *) 0) {
	    restart++;
	    if (stat("/users/c/chaos/bin/format_server_report", &stat_buf) == 0) {
	        if(system("/users/c/chaos/bin/format_server_report restart") == -1) {
		    printf("RESTART FAILED\n");
		}
	    }
	}
	log = fopen(format_server_log, "a");

	if (log == NULL) {
	    log = (void *) -1;
	    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	    pthread_mutex_unlock(&fs->log_lock);
	    return;
	}
	if (restart) {
	    fprintf(log, "\n\nRestarting format server \n");
	}
    }
    time( &rawtime );

    info = localtime( &rawtime );

    strftime(buffer,80,"%x - %I:%M%p", info);

    fprintf(log, "%s - ", buffer);

#ifdef STDC_HEADERS
    va_start(ap, format);
#else
    va_start(ap);
#endif
    vfprintf(log, format, ap);
    va_end(ap);
    fprintf(log, "\n");
    fflush(log);
    if (fs->stdout_verbose) {
#ifdef STDC_HEADERS
	va_start(ap, format);
#else
	va_start(ap);
#endif
	vfprintf(stdout, format, ap);
	va_end(ap);
	fprintf(stdout, "\n");
    }
    stat(format_server_log, &stat_buf);

    if (log_count++ >= 1000) {
	log_count = 0;
	if (stat_buf.st_size > 1024000) {
	    dump_stats_to_log(fs);
	    fflush(log);
	    fclose(log);
	    if (stat("/users/c/chaos/bin/format_server_report", &stat_buf) == 0) {
	        if(system("/users/c/chaos/bin/format_server_report log") == -1) {
		    printf("log failed\n");
	        }
	    }
	    unlink(format_server_log);
	    
	    log = (void *) -1;
	}
    }
    pthread_mutex_unlock(&fs->log_lock);
}

static int
LogAtomicRead(format_server fs, void *fd, void *buffer, long length)
{
    char *junk_result_str = NULL;
    int junk_errno;
    int ret = ffs_server_read_func(fd, buffer, length, &junk_errno,
				  &junk_result_str);

    if (getenv("BAD_CLIENT") && (drand48() < 0.0001)) sleep(600);
    if (ret != length) {
	LOG(fs, "Read failed for fd %d, return %d, length %d, errno %d\n",
	    fd, ret, length, junk_errno);
    }
    return ret;
}

#include "self_ip_addr.c"

static int my_pid = 0;
static INT4 my_IP_addr = 0;

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 64
#endif

void
general_format_server(int port, int do_restart, int verbose, int do_proxy)
{
    format_server fs = format_server_create();
    fs->stdout_verbose = verbose;
    if (do_proxy) {
	fs->proxy_context_to_master = create_FMcontext();
    }
    if (format_server_listen(fs, port) != -1) {
	if (do_restart) read_formats_from_file(fs);
	LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
	pthread_mutex_lock(&fs->lock);
	LOG(fs, "  - got the mutex at Line %d ", __LINE__);
	while (1) {
	    format_server_poll_and_handle(fs);
	}
#ifdef __NVCOMPILER
#pragma diag_suppress 111
#endif
	LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	pthread_mutex_unlock(&fs->lock);
#ifdef __NVCOMPILER
#pragma diag_default 111
#endif
    }
    return;
}

static void FSClient_close(FSClient fsc);
static void FSClient_force_close(FSClient fsc);
#define CONN_TIMEOUT_INTERVAL 300

static int test_count = 0;
static int registration_count = 0;
static int get_count = 0;
static int select_handle_count = 0;
static int total_client_count = 0;
static int rejected_client_count = 0;

static int
format_server_poll_and_handle(format_server fs)
{
    int i, res;
    fd_set rd_set = fs->fdset;
    struct timeval timeout;

#if !defined(__FD_SET) && !defined(HAVE_WINDOWS_H)
    LOG(fs, "  selecting, rd_set is %lx,%lx,%lx,%lx,%lx,%lx", rd_set.fds_bits[0],
	rd_set.fds_bits[1], rd_set.fds_bits[2], rd_set.fds_bits[3],
	rd_set.fds_bits[4], rd_set.fds_bits[5]);
#endif
    /* scan existing conn for incoming data & do something */
    timeout.tv_usec = 0;
    timeout.tv_sec = CONN_TIMEOUT_INTERVAL + 5;

    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
    pthread_mutex_unlock(&fs->lock);
    res = select(FD_SETSIZE, &rd_set, (fd_set *) NULL,
		 (fd_set *) NULL, &timeout);
    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
    pthread_mutex_lock(&fs->lock);
    LOG(fs, "  - got the mutex at Line %d ", __LINE__);
    select_handle_count++;
    if (res == -1) {
	if (errno == EBADF) {
	    int i;
	    int found_one = 0;
	    for (i = 0; i < FD_SETSIZE; i++) {
		if (FD_ISSET(i, &rd_set) 
		    && (i != (long) fs->conn_sock_inet)
		    && (i != (long) fs->conn_sock2_inet)) {
		    fd_set test_set;
		    timeout.tv_usec = 0;
		    timeout.tv_sec = 0;
#ifdef __NVCOMPILER
#pragma diag_suppress 550
#endif
		    FD_ZERO(&test_set);
#ifdef __NVCOMPILER
#pragma diag_default 550
#endif
		    FD_SET(i, &test_set);
		    errno = 0;
		    select(FD_SETSIZE, &test_set, (fd_set *) NULL,
			   (fd_set *) NULL, &timeout);
		    if (errno == EBADF) {
			fprintf(stderr, "Select failed, fd %d\n is bad.  Removing from select list.\n",
				i);
			LOG(fs, "REmoving bad fd %d", i);
			if (fs->ports[i]) {
			    FSClient_force_close(fs->ports[i]);
			} else if (i != (long) fs->conn_sock_inet) {
			    FD_CLR((unsigned long) i, &fs->fdset);
			    ffs_close_func((void *) (long) i);
			}
			found_one++;
		    }
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
	timeout_old_conns(fs);
	if (fs->pending_unregistered) try_master_connect(fs);
	return 0;
    }
    if ((long) fs->conn_sock_inet >= 0 &&
	FD_ISSET((unsigned long) fs->conn_sock_inet, &rd_set)) {
	FD_CLR((unsigned long) fs->conn_sock_inet, &rd_set);
	format_server_accept_conn_sock(fs, fs->conn_sock_inet);
    }
    if ((long) fs->conn_sock2_inet >= 0 &&
	FD_ISSET((unsigned long) fs->conn_sock2_inet, &rd_set)) {
	FD_CLR((unsigned long) fs->conn_sock2_inet, &rd_set);
	format_server_accept_conn_sock(fs, fs->conn_sock2_inet);
    }
    for (i = 0; i < FD_SETSIZE; i++) {	/* check each possible read fd */
	FSClient fsc = fs->ports[i];
	if (fsc == NULL || !FD_ISSET(i, &rd_set))
	    continue;
	fs->timestamp[i] = time(NULL);
	format_server_handle_data(fs, fsc);

    }
    timeout_old_conns(fs);
    if (fs->pending_unregistered) try_master_connect(fs);
    return res;
}

static int
format_eq(IOFormatRep form1, IOFormatRep form2)
{
    int i;
    int all_zero = 1;
    int all_eq = 1;
    if (form1->server_ID.length != 0) {
	if (form1->server_ID.length != form2->server_ID.length)
	    return 0;
	for (i = 0; i < form1->server_ID.length; i++) {
	    if (form1->server_ID.value[i] != 0) {
		all_zero = 0;
	    }
	    if (form1->server_ID.value[i] !=
		form2->server_ID.value[i]) {
		all_eq = 0;
	    }
	}
    }
    if (all_eq && !all_zero)
	return 1;
    if ((form1->server_format_rep != 0) && (form2->server_format_rep != 0)) {
	/* both have server reps.  Equal if server reps are equal */
	if (ntohs(form1->server_format_rep->format_rep_length) !=
	    ntohs(form2->server_format_rep->format_rep_length))
	    return 0;
	return (memcmp(form1->server_format_rep, form2->server_format_rep,
		       ntohs(form1->server_format_rep->format_rep_length)) == 0);
    }
    /* if the format IDs were not equal or the reps were not equal, then they aren't equal */
    return 0;
}

static int server_format_count = 0;

static void
byte_swap(char *data, int size)
{
    int i;
    assert((size % 2) == 0);
    for (i = 0; i < size / 2; i++) {
	char tmp = data[i];
	data[i] = data[size - i - 1];
	data[size - i - 1] = tmp;
    }
}

static void
dump_formats_to_file(format_server fs)
{
    long fd = open(FS_restart_file, O_WRONLY | O_CREAT | O_TRUNC, 0600);

    if (fd == -1) {
	perror("opening dump file");
	return;
    }
    LOG(fs, "Attempting to Dump format to file %s \n", FS_restart_file);
    format_list *list = fs->lists[0];
    while (list != NULL) {
	IOFormatRep ioformat = list->format;
	format_rep rep = ioformat->server_format_rep;
	UINT2 id_len = ioformat->server_ID.length;
	UINT2 rep_len = ntohs(rep->format_rep_length);
	void *id = ioformat->server_ID.value;
	char *junk_result_str;
	int junk_errno;
	if (ffs_server_write_func((void *) fd, &id_len, sizeof(id_len),
				 &junk_errno, &junk_result_str)
	    != sizeof(id_len)) {
	    fprintf(stderr, "error writing dump id len\n");
	    return;
	}
	if (ffs_server_write_func((void *) fd, id, id_len, &junk_errno,
				 &junk_result_str) != id_len) {
	    fprintf(stderr, "error writing dump id \n");
	    return;
	}
	if (ffs_server_write_func((void *) fd, rep, rep_len, &junk_errno,
				 &junk_result_str) != rep_len) {
	    fprintf(stderr, "error writing dump format rep\n");
	    return;
	}
	list = list->next;
    }
    close(fd);

    LOG(fs, "Successfully Dumped format to file %s \n", FS_restart_file);

}

static void
read_formats_from_file(format_server fs)
{
    long fd = open(FS_restart_file, O_RDONLY, 0600);
    format_list *last;
    format_rep rep;

    last = NULL;

    LOG(fs, "Attempting to Read format from file %s \n", FS_restart_file);
    while (1) {
	IOFormatRep ioformat;
	format_list *new;
	UINT2 id_len, rep_len;
	unsigned char id[128];		/* too big */
	char buffer[256];
	unsigned char *tmp_id;

	if (LogAtomicRead(fs, (void *) fd, (void*)&id_len, sizeof(id_len)) != sizeof(id_len)) {
	    close(fd);
	    LOG(fs, "Successfully Read %d formats from file %s \n", server_format_count, FS_restart_file);
	    return;
	}
	if (id_len > sizeof(id)) {
	    close(fd);
	    LOG(fs, "ID size too large\n");
	    LOG(fs, "Read %d formats from file %s \n", server_format_count, FS_restart_file);
	    return;
	}
	if (LogAtomicRead(fs, (void *) fd, &id[0], (long)id_len) != id_len) {
	    close(fd);
	    LOG(fs, "Abnormal restart termination after id_len read, truncated file?\n");
	    LOG(fs, "Read %d formats from file %s \n", server_format_count, FS_restart_file);
	    return;
	}
	if (LogAtomicRead(fs, (void *) fd, (void*)&rep_len, sizeof(rep_len)) != sizeof(rep_len)) {
	    close(fd);
	    LOG(fs, "Abnormal restart termination after id read, truncated file?\n");
	    LOG(fs, "Read %d formats from file %s \n", server_format_count, FS_restart_file);
	    return;
	}
	rep = malloc(rep_len);
	rep->format_rep_length = htons(rep_len);
	if (LogAtomicRead(fs, (void *) fd, (void*)(((char *) rep) + sizeof(rep_len)),
			     rep_len - sizeof(rep_len)) !=
	    (rep_len - sizeof(rep_len))) {
	    close(fd);
	    LOG(fs, "Abnormal restart termination after rep_len read, truncated file?\n");
	    LOG(fs, "Read %d formats from file %s \n", server_format_count, FS_restart_file);
	    return;
	}
	stringify_server_ID(id, buffer, sizeof(buffer));
	LOG(fs, "Read format %s from file.\n", buffer);
	tmp_id = malloc(id_len);
	memcpy(tmp_id, id, id_len);
	ioformat = fetch_format(fs, NULL, tmp_id, id_len);
	if (ioformat == NULL) {
	    ioformat = malloc(sizeof(*ioformat));
	    ioformat->server_format_rep = rep;
	    ioformat->server_ID.length = id_len;
	    ioformat->server_ID.value = (char*) tmp_id;
	    new = (format_list *) malloc(sizeof(format_list));
	    new->format = ioformat;
	    new->next = NULL;
	    if (last == NULL) {
		fs->lists[0] = new;
	    } else {
		last->next = new;
	    }
	    last = new;
	    server_format_count++;
	} else {
	    /* ioformat already known */
	    free(rep);
	}
    }
}


static void
register_all_to_master(format_server fs);

static int
try_master_connect(format_server fs)
{
    static time_t last_try = 0;
    struct timeval now;
    int try_connection = 0;
    if (fs->proxy_context_to_master == NULL) return 0;
    if (fs->proxy_context_to_master->server_fd != (void*)-1) {
	/* assume connection still good */
	if (fs->stdout_verbose) printf("Master connection still good\n");
	return 1;
    }
    gettimeofday(&now, NULL);
    if (last_try == 0) try_connection++;
    if ((now.tv_sec - last_try) > 300) try_connection++;
    if (!try_connection) {
	if (fs->stdout_verbose) printf("Master connection dead, too soon to retry\n");
	return 0;
    }
    last_try = now.tv_sec;
    if (fs->stdout_verbose) printf("Master connection dead, try reconnect\n");
    establish_server_connection(fs->proxy_context_to_master, never_local);
    if (fs->proxy_context_to_master->server_fd != (void*)-1) {
	int s = (int)(long)fs->proxy_context_to_master->server_fd;
	int optval;
	socklen_t optlen = sizeof(optval);
	if (fs->stdout_verbose) printf("Master connection dead, reconnect successful!\n");
	/* Set the option active */
	optval = 1;
	optlen = sizeof(optval);
	if(setsockopt(s, SOL_SOCKET, SO_KEEPALIVE, &optval, optlen) < 0) {
	    perror("setsockopt()");
	    close(s);
	    exit(EXIT_FAILURE);
	}
	register_all_to_master(fs);
	return 1;
    } 
    if (fs->stdout_verbose) printf("Master connection dead, reconnect failed\n");
    return 0;
}

static void
register_format_to_master(format_server fs, IOFormatRep ioformat)
{
    
    if (fs->stdout_verbose) printf("Trying to register to master\n");
    if (!try_master_connect(fs)) {
	fs->pending_unregistered++;
	return;
    } else {
	struct {
	    char reg[2];
	    unsigned short len;
	} tmp = {{'P', 2}, 0 };	/* format push, version 2 */
	int ret;
	void* server_fd = fs->proxy_context_to_master->server_fd;
	int errno;
	char *errstr;

	format_rep rep = ioformat->server_format_rep;

	tmp.len = ntohs(ioformat->server_ID.length);
	ret = ffs_server_write_func(server_fd, &tmp, sizeof(tmp),
				   &errno, &errstr);
	if (ret != sizeof(tmp)) {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Write failed1\n");
	    return;
	}
	ret = ffs_server_write_func(server_fd, 
				   ioformat->server_ID.value,
				   ioformat->server_ID.length,
				   &errno, &errstr);
	if (ret != ioformat->server_ID.length) {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Write failed2\n");
	    return;
	}

	if (fs->stdout_verbose) {
	    printf("Pushing ");
	    print_server_ID( (unsigned char *) ioformat->server_ID.value);
	    printf("Writing %d bytes of format rep\n", ntohs(rep->format_rep_length));
	}
	ret = ffs_server_write_func(server_fd, &rep->format_rep_length, 2, &errno, &errstr);
	if (ret != 2)  {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Write failed3\n");
	    return;
	}
	ret = ffs_server_write_func(server_fd, (char*) rep,
				   ntohs(rep->format_rep_length),
				   &errno, &errstr);
	if (ret != ntohs(rep->format_rep_length)) {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Write failed\n");
	    return;
	}

	if (fs->stdout_verbose) printf("Done\n");
    }
}

static void
register_all_to_master(format_server fs)
{
    if (fs->stdout_verbose) printf("Trying to register all master\n");
    if (!try_master_connect(fs)) return;
    format_list *list = fs->lists[0];
    while (list != NULL) {
	register_format_to_master(fs, list->format);
	list = list->next;
    }
    fs->pending_unregistered = 0;
}

static IOFormatRep
get_format_from_master(format_server fs, IOFormatRep ioformat)
{
    char get[2] = {'g', 8};		/* format get, size */
    char block_version;
    UINT2 length;
    char return_char = 0;
    format_rep rep;
    void* server_fd;
    int errno, id_size;
    char *errstr;
    int ret;

    if (fs->stdout_verbose) printf("Trying to get from master\n");
    if (!try_master_connect(fs)) return NULL;

    server_fd = fs->proxy_context_to_master->server_fd;
    id_size = get[1] = ioformat->server_ID.length;
    if ((ret = ffs_server_write_func(server_fd, &get[0], 2, &errno, &errstr)) != 2) {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Write failed2.1, ret is %d\n", ret);
	return NULL;
    }
    if (ffs_server_write_func(server_fd, ioformat->server_ID.value, id_size, &errno, &errstr) != id_size) {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Write failed3\n");
	return NULL;
    }
    if (LogAtomicRead(fs, server_fd, &return_char, (long)1) != 1) {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Read failed4\n");
	return NULL;
    }
    if (return_char == 'P') {
        if (LogAtomicRead(fs, server_fd, &return_char, (long)1) != 1) {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Read failed4\n");
	    return NULL;
	}
    }
    if (return_char != 'f') {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Bad character\n");
	return NULL;
    }
    if (LogAtomicRead(fs, server_fd, &block_version, (long)1) != 1) {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Bad character\n");
	return NULL;
    }
    if (block_version != 1) {
	fprintf(stderr, "Unknown version \"%d\"in block registration\n", block_version);
	return NULL;
    }
    if (LogAtomicRead(fs, server_fd, (void*)&length, sizeof(length)) !=
	sizeof(length)) {
	close((int)(long)server_fd);
	fs->proxy_context_to_master->server_fd = (void*)-1;
	LOG(fs, "Read failed\n");
	return NULL;
    }
    length = ntohs(length);
    if (length == 0) {
	return NULL;
    } else {
	rep = malloc(length);
	rep->format_rep_length = htons((short)length);
	if (LogAtomicRead(fs, server_fd, (void*)(((char *) rep) + sizeof(length)),
			     length - sizeof(length)) != (length - sizeof(length))) {
	    close((int)(long)server_fd);
	    fs->proxy_context_to_master->server_fd = (void*)-1;
	    LOG(fs, "Read failed\n");
	    return NULL;
	}
	ioformat->server_format_rep = rep;
	return ioformat;
    }
#ifdef __NVCOMPILER
#pragma diag_suppress 111
#endif
    return NULL;
#ifdef __NVCOMPILER
#pragma diag_default 111
#endif
}


static IOFormatRep
find_format(format_server fs, FSClient fsc, IOFormatRep ioformat, int new_format_mode, int requested_id_version)
{
    format_list *list = fs->lists[0];
    format_list *last = NULL, *new = NULL;
    /* always, search list for this hosts byte order */
    if (ioformat->server_ID.value) {
	if ( version_of_format_ID(ioformat->server_ID.value) !=
	     requested_id_version) {
	    LOG(fs, "Mismatched versions in request\n");
	}
    }

    while (list != NULL) {
	int server_ID_version =
	    version_of_format_ID(list->format->server_ID.value);
	if (format_eq(ioformat, list->format) &&
	    (requested_id_version == server_ID_version)) {
	    if (ioformat->server_ID.value) {
		if (list->format->server_ID.length != ioformat->server_ID.length) {
		    LOG(fs, "Version 2 IDs differ in length\n");
		}
		if (memcmp(list->format->server_ID.value, ioformat->server_ID.value,
			   list->format->server_ID.length) != 0) {
		    LOG(fs, "Version 2 ID values differ\n");
		}
	    }
	    if (ioformat->server_ID.value != NULL) free(ioformat->server_ID.value);
	    if (ioformat->server_format_rep != NULL) free(ioformat->server_format_rep);
	    free(ioformat);
	    return list->format;
	}
	last = list;
	list = list->next;
    }
    if (new_format_mode) {
	INT2 port = fs->port;
	char *tmp;
	/*
	 * Didn't find it in the list, create a blank.
	 * This is the point at which we know we are registering a 
	 * previously-unregistered format.
	 */
	new = (format_list *) malloc(sizeof(format_list));
	new->format = ioformat;
	switch (requested_id_version) {
	case 0:
	    tmp = (char *) ioformat;
	    if (!WORDS_BIGENDIAN)
		byte_swap(tmp, (int)sizeof(tmp));

	    tmp += (0x7 & my_pid);	/* add bottom 3 bits of pid to
					 * format id */
	    ioformat->server_ID.length = 8;
	    ioformat->server_ID.value = malloc(8);
	    memcpy(ioformat->server_ID.value, &tmp, sizeof(ioformat));
	    memset(ioformat->server_ID.value + sizeof(ioformat), 0,
		   ioformat->server_ID.length - sizeof(ioformat));
	    if (sizeof(ioformat) == 8) {
		/* make sure the top byte is zero */
		*((char *) ioformat->server_ID.value) = 0;
	    }
	    break;
	case 1:
	    ioformat->server_ID.length = 10;
	    ioformat->server_ID.value = malloc(10);
	    ((version_1_format_ID *) ioformat->server_ID.value)->version = 1;
	    ((version_1_format_ID *) ioformat->server_ID.value)->salt =
		(my_pid & 0xff);
	    ((version_1_format_ID *) ioformat->server_ID.value)->port = port;
	    ((version_1_format_ID *) ioformat->server_ID.value)->IP_addr =
		my_IP_addr;
	    ((version_1_format_ID *) ioformat->server_ID.value)->format_identifier = server_format_count;
	    break;
	case 2:
	    if (!ioformat->server_ID.value) {
		generate_format3_server_ID(&ioformat->server_ID, ioformat->server_format_rep);
	    } else {
		server_ID_type tmp;
		generate_format3_server_ID(&tmp, ioformat->server_format_rep);
		if (tmp.length != ioformat->server_ID.length) {
		    LOG(fs, "Version 2 IDs differ in length\n");
		}
		if (memcmp(tmp.value, ioformat->server_ID.value,
			   tmp.length) != 0) {
		    LOG(fs, "Version 2 ID values differ\n");
		}
		free(tmp.value);
	    }
	    break;
	}
	new->next = NULL;
	server_format_count++;
	if (last == NULL) {
	    fs->lists[0] = new;
	} else {
	    last->next = new;
	}
	register_format_to_master(fs, ioformat);
	return ioformat;
    } else {
	if (fs->proxy_context_to_master != NULL) {
	    /* we're a proxy, perhaps the master knows this format */
	    IOFormatRep tmp = get_format_from_master(fs, ioformat);
	    if (tmp == NULL) return NULL;

	    /* we got a format, free the local version */
	    if (ioformat->server_ID.value != NULL) free(ioformat->server_ID.value);
	    if (ioformat->server_format_rep != NULL) free(ioformat->server_format_rep);
	    free(ioformat);
	    return tmp;
	}
	return NULL;
    }
}

/* 
 * Close an on-line connection.
 */
static void
FSClient_close(FSClient fsc)
{
    format_server fs = fsc->fs;
    int fd = (int) (long) fsc->fd;

    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
    pthread_mutex_lock(&fs->lock);
    LOG(fs, "  - got the mutex at Line %d ", __LINE__);
    LOG(fs, "Closing client fd %d - %s", fd, fsc->hostname ?
	fsc->hostname : "NULL");
    if (fd != 0) {
	fs->ports[fd] = NULL;
	fs->timestamp[fd] = 0;
	FD_CLR((unsigned long) fd, &fs->fdset);
	fs->portCount--;
    }
    if (fsc->hostname) {
	free(fsc->hostname);
	fsc->hostname = NULL;
    }
    free(fsc);
    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
    pthread_mutex_unlock(&fs->lock);
    /* moved close outside of the mutex because sometimes it blocks */
    ffs_close_func((void *) (long) fd);
    if (fs->fork_threads) pthread_exit(NULL);
}

static void
FSClient_force_close(FSClient fsc)
{
    format_server fs = fsc->fs;
    int fd = (int) (long) fsc->fd;
    LOG(fs, "Forcibly closing client fd %d- %s", fd, fsc->hostname ?
	fsc->hostname : "NULL");
    if (fd != 0) {
	fs->ports[fd] = NULL;
	fs->timestamp[fd] = 0;
	FD_CLR((unsigned long) fd, &fs->fdset);
	shutdown(fd, SHUT_RDWR);
//	pthread_kill(fsc->handler_thread, 13);
	pthread_join(fsc->handler_thread, NULL);
	ffs_close_func((void *) (long) fd);
	fs->portCount--;
//	fsc->fd = 0;
    }
//    if (fsc->hostname) free(fsc->hostname);
//    free(fsc);
}

static IOFormatRep
fetch_format(format_server fs, FSClient fsc, const unsigned char *format_id,
	     int format_id_length)
{
    IOFormatRep ioformat;
    ioformat = malloc(sizeof(*ioformat));
    ioformat->server_format_rep = NULL;
    ioformat->server_ID.length = format_id_length;
    ioformat->server_ID.value = (char *) format_id;
    if (fs->stdout_verbose) {
	print_server_ID( (unsigned char *) ioformat->server_ID.value);
	printf("\n");
    }
    ioformat = find_format(fs, fsc, ioformat,
			   0 /* not new format mode */ ,
			   version_of_format_ID(ioformat->server_ID.value));
    return ioformat;
}


typedef struct host_info {
    struct in_addr ip_addr;
    time_t intro_time;
} *host_info_p;

static host_info_p hostlist = NULL;
static int host_count = 0;


extern int
get_internal_format_server_identifier(format_server fs)
{
    return fs->format_server_identifier;
}

#define FILE_INT INT4

static int
put_serverAtomicInt(void *fd, FILE_INT *file_int_ptr, FMContext fmc)
{
#if SIZEOF_INT == 4
    int tmp_value = *file_int_ptr;
    int junk_errno;
    char *junk_result_str;
    if (ffs_server_write_func(fd, &tmp_value, 4, &junk_errno, &junk_result_str) != 4) {
	printf("SERVER WRITE FAILED, ERRNO = %d\n", junk_errno);
	return 0;
    }
#else
    Baaad shit;
#endif
    return 1;
}

static int
get_serverAtomicInt(void *fd, FILE_INT *file_int_ptr, int byte_reversal)
{
#if SIZEOF_INT == 4
    int tmp_value;
    int junk_errno;
    char *junk_result_str;
    if (ffs_server_read_func(fd, &tmp_value, 4, &junk_errno, &junk_result_str) != 4) {
	printf("SERVER READ FAILED, ERRNO = %d\n", junk_errno);

	return 0;
    }
#else
    Baaad shit;
#endif
    if (byte_reversal)
	byte_swap((char *) &tmp_value, 4);
    *file_int_ptr = tmp_value;
    return 1;
}

extern int
unix_timeout_read_func(void *conn, void *buffer, int length, 
		       int *errno_p, char **result_p);

static int
get_serverTimeoutInt(void *fd, FILE_INT *file_int_ptr, int byte_reversal)
{
#if SIZEOF_INT == 4
    int tmp_value;
    int junk_errno;
    char *junk_result_str;
    if (unix_timeout_read_func(fd, &tmp_value, 4, &junk_errno, &junk_result_str) != 4)
	return 0;
#else
    Baaad shit;
#endif
    if (byte_reversal)
	byte_swap((char *) &tmp_value, 4);
    *file_int_ptr = tmp_value;
    return 1;
}

static void
server_read_header(FSClient fsc)
{
    FILE_INT magic;
    FILE_INT float_format;
    FILE_INT pid = getpid();

    void *fd = fsc->fd;
    int byte_reversal = 0;
    int version = -1;

    get_serverTimeoutInt(fd, &magic, byte_reversal);

    switch (magic) {
    case MAGIC_NUMBER:
	break;
    case REVERSE_MAGIC_NUMBER:
	byte_reversal = 1;
	break;
    case MAGIC_NUMBER + 1:
	version = 1;
	break;
    case (FILE_INT) REVERSE_MAGIC_NUMBER + 0x1000000:
	version = 1;
	byte_reversal = 1;
	break;
    case MAGIC_NUMBER + 2:
	version = 2;
	break;
    case (FILE_INT) REVERSE_MAGIC_NUMBER + 0x2000000:
	version = 2;
	byte_reversal = 1;
	break;
    case MAGIC_NUMBER + 3:
	version = 3;
	break;
    case (FILE_INT) REVERSE_MAGIC_NUMBER + 0x3000000:
	version = 3;
	byte_reversal = 1;
	break;
    case MAGIC_NUMBER + 4:
	version = 4;
	break;
    case (FILE_INT) REVERSE_MAGIC_NUMBER + 0x4000000:
	version = 4;
	byte_reversal = 1;
	break;
    default:
	close((int)(long)fd);
	return;
    }
    if (version <= 2) {
	get_serverAtomicInt(fd, &float_format, byte_reversal);
    } else {
	FILE_INT key_len;
	unsigned char *key;
	get_serverAtomicInt(fd, &key_len, byte_reversal);
	if (key_len > 1024) {
	    close((int)(long)fd);
	    return;
	}
	fsc->key = NULL;
	fsc->key_len = key_len;
	if (key_len > 0) {
	    key = malloc(key_len);
	    LogAtomicRead(fsc->fs, fd, key, (long)key_len);
	    fsc->key = key;
	    fsc->key_len = key_len;
	}
    }
    magic = MAGIC_NUMBER;
    put_serverAtomicInt(fd, &magic, (FMContext) NULL);
    if (version >= 2) {
	put_serverAtomicInt(fd, &pid, (FMContext) NULL);
    }
    if (version >= 3) {
	int format_server_identifier = get_internal_format_server_identifier(fsc->fs);
	put_serverAtomicInt(fd, &format_server_identifier, (FMContext) NULL);
    }
    fsc->byte_reversal = byte_reversal;
    fsc->version = version;
}

static void
format_server_handle_data(format_server fs, FSClient fsc)
{
    char next_action;
    int input_bytes = 0;
    void* fd = fsc->fd;
    LOG(fs, "    waiting for read on fd %d - %s", fsc->fd, fsc->hostname);
    if (LogAtomicRead(fs, fd, &next_action, (long)1) != 1) {
	LOG(fs, "Close on atomic read fail");
	FSClient_close(fsc);
	return;
    }
    input_bytes++;
    switch (next_action) {
    case 'T':
	LOG(fs, "    Testing ");
	test_count++;
	break;
    case 'C':
	LOG(fs, "    Close ");
	FSClient_close(fsc);
	break;
    case 'f': {
	    IOFormatRep ioformat;
	    char block_version;
	    UINT2 length;

	    format_rep rep;
	    if (fs->stdout_verbose) {
		printf("New style block format registration\n");
	    }
	    LOG(fs, "    Registration ");
	    if (LogAtomicRead(fs, fd, &block_version, (long)1) != 1) {
		LOG(fs, "Close on atomic read fail 2");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes++;
	    if (block_version > 2) {
		fprintf(stderr, "Unknown version in block registration\n");
		LOG(fs, "Close on unknown version");
		FSClient_close(fsc);
		return;
	    }
	    if (LogAtomicRead(fs, fd, (void*)&length, sizeof(length)) !=
		sizeof(length)) {
		LOG(fs, "Close on atomic read fail 3");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes += sizeof(length);
	    length = ntohs(length);
	    rep = malloc(length);
	    rep->format_rep_length = htons((short)length);
	    if (LogAtomicRead(fs, fd, (void*)(((char *) rep) + sizeof(length)),
				 length - sizeof(length)) !=
		(length - sizeof(length))) {
		LOG(fs, "Close on atomic read fail 4");
		FSClient_close(fsc);
		return;
	    }
    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
	    pthread_mutex_lock(&fs->lock);
	LOG(fs, "  - got the mutex at Line %d ", __LINE__);
	    fsc->input_bytes += input_bytes + (length - sizeof(length));
	    ioformat = malloc(sizeof(*ioformat));
	    ioformat->server_format_rep = rep;
	    ioformat->server_ID.length = 0;
	    ioformat->server_ID.value = NULL;
	    ioformat = find_format(fs, fsc, ioformat,
				   1 /* new format mode */ , block_version);
	    if (fs->stdout_verbose) {
		printf("Returning -> ");
		print_server_ID( (unsigned char *) ioformat->server_ID.value);
		printf("\n");
	    }
    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	    pthread_mutex_unlock(&fs->lock);
	    {
		char ret[2];
		ret[0] = 'I';
		/*		if (fsc->provisional != 0) {
		    ret[0] = 'P';
		    } */
		ret[1] = ioformat->server_ID.length;
		if (ffs_server_write_func(fd, &ret[0], 2, NULL, NULL) != 2) {
		    LOG(fs, "Close on atomic write fail");
		    FSClient_close(fsc);
		    return;
		}
		if (ffs_server_write_func(fd, ioformat->server_ID.value,
				  ioformat->server_ID.length, NULL,
			     NULL) != ioformat->server_ID.length) {
		    LOG(fs, "Close on atomic write fail 2");
		    FSClient_close(fsc);
		    return;
		}
	    }
	    registration_count++;
	    fsc->formats_registered++;
	    break;
	}
    case 'P': {
	    char block_version;
	    UINT2 length;

	    format_rep rep;
	    IOFormatRep ioformat;
	    char *format_ID;
	    int format_ID_len;
	    int input_bytes = 0;
	    if (fs->stdout_verbose) {
		printf("Push block format registration\n");
	    }
	    LOG(fs, "    Push Registration ");
	    if (LogAtomicRead(fs, fd, &block_version, (long)1) != 1) {
		LOG(fs, "Close on atomic read fail P1");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes++;
	    if (block_version <= 1) {
		fprintf(stderr, "Unknown version in block registration\n");
		LOG(fs, "Close on unknown version block reg 1");
		FSClient_close(fsc);
		return;
	    }
	    /* read pushed format ID */
	    if (LogAtomicRead(fs, fd, (void*)&length, sizeof(length)) !=
		sizeof(length)) {
		LOG(fs, "Close on atomic read fail P2");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes += sizeof(length);
	    format_ID_len = length = ntohs(length);
	    format_ID = malloc(format_ID_len);
	    if (LogAtomicRead(fs, fd, ((char *) format_ID), (long)format_ID_len)
		!= format_ID_len) {
		LOG(fs, "Close on atomic read fail P3");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes += (length - sizeof(length));

	    /* read body */
	    if (LogAtomicRead(fs, fd, (void*)&length, sizeof(length)) !=
		sizeof(length)) {
		LOG(fs, "Close on atomic read fail P4");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes += sizeof(length);
	    length = ntohs(length);
	    rep = malloc(length);
	    rep->format_rep_length = htons(length);
	    if (LogAtomicRead(fs, fd, ((char *) rep), (long)length) != length) {
		LOG(fs, "Close on atomic read fail P5");
		FSClient_close(fsc);
		return;
	    }

	    
    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
	    pthread_mutex_lock(&fs->lock);
	LOG(fs, "  - got the mutex at Line %d ", __LINE__);
	    fsc->input_bytes += input_bytes + (length - sizeof(length));
	    if (fs->stdout_verbose) {
		printf("Got Pushed -> ");
		print_server_ID( (unsigned char *) format_ID);
		printf("\n");
	    }
	    /* GSE create format rep */
	    ioformat = malloc(sizeof(*ioformat));
	    ioformat->server_format_rep = rep;
	    ioformat->server_ID.length = format_ID_len;
	    ioformat->server_ID.value = format_ID;
	    ioformat = find_format(fs, fsc, ioformat,
				   1 /* new format mode */ , block_version);


	    registration_count++;
	    fsc->formats_registered++;
    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	    pthread_mutex_unlock(&fs->lock);
	    break;
	}
    case 'g': {
	    IOFormatRep ioformat;
	    char format_id_length;
	    int input_bytes = 0;
	    struct {
		char reg[2];
		unsigned short len;
	    } tmp = { { 'f', 1 }, 0 };		/* format reg, version 1 */
	    int ret;
	    int errno;
	    char *errstr;
	    unsigned char *format_id_value;
	    LOG(fs, "    new style get format");
	    if (fs->stdout_verbose)
		printf("new style Get Format   -> ");
	    if (LogAtomicRead(fs, fd, &format_id_length, (long)1) != 1) {
		LOG(fs, "Close on atomic read fail g1");
		FSClient_close(fsc);
		return;
	    }
	    input_bytes++;
	    format_id_value = malloc(format_id_length);
	    if (LogAtomicRead(fs, fd, format_id_value,
			      (long)format_id_length) != format_id_length) {
		LOG(fs, "Close on atomic read fail g2");
		FSClient_close(fsc);
		return;
	    }
    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
	    pthread_mutex_lock(&fs->lock);
	LOG(fs, "  - got the mutex at Line %d ", __LINE__);
	    fsc->input_bytes += input_bytes + format_id_length;
	    ioformat = fetch_format(fs, fsc, format_id_value, 
				    format_id_length);
	    if (fs->stdout_verbose) {
		printf("Returning -> ");
		if (ioformat == NULL) {
		    printf("	NULL\n");
		}
	    }
	    if ((ioformat != NULL) &&
		(ioformat->server_format_rep != NULL)) {
		tmp.len = ioformat->server_format_rep->format_rep_length;
	    } else {
		tmp.len = 0;
	    }
	    fsc->formats_fetched++;
    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	    pthread_mutex_unlock(&fs->lock);
	    ret = ffs_server_write_func(fd, &tmp, sizeof(tmp),
				       &errno, &errstr);
	    if (ret != sizeof(tmp)) {
		LOG(fs, "Close on atomic read fail g3");
		FSClient_close(fsc);
		break;
	    }
	    if ((ioformat != NULL) &&
		(ioformat->server_format_rep != NULL)) {
		format_rep rep = ioformat->server_format_rep;
		ret = ffs_server_write_func(fd,
					   (char *) rep + sizeof(rep->format_rep_length),
					   ntohs(rep->format_rep_length) - sizeof(rep->format_rep_length),
					   &errno, &errstr);
		if (ret != ntohs(rep->format_rep_length) - sizeof(rep->format_rep_length)) {
		    LOG(fs, "Close on atomic read fail g4");
		    FSClient_close(fsc);
		    break;
		}
	    }
	    get_count++;
	    break;
	}
    case 'D':
	LOG(fs, "	Dump known format IDs");
	{
	    int junk_errno;
	    char *junk_str;
	    int out_count = 0, i;
	    int byte_order;
	    int tmp;
	    server_ID_type *out_list = malloc(sizeof(out_list[0]));
    LOG(fs, "Doing Mutex Lock at Line %d ", __LINE__);
	    pthread_mutex_lock(&fs->lock);
	LOG(fs, "  - got the mutex at Line %d ", __LINE__);
	    for (byte_order = 0; byte_order <= 1; byte_order++) {
		format_list *list = fs->lists[0];
		while (list != NULL) {
		    out_list = realloc(out_list, sizeof(out_list[0]) *
				       (out_count + 1));
		    out_list[out_count] = list->format->server_ID;
		    out_count++;
		    list = list->next;
		}
	    }
	    tmp = htonl(out_count);
    LOG(fs, "Doing Mutex unLock at Line %d ", __LINE__);
	    pthread_mutex_unlock(&fs->lock);
	    ffs_server_write_func(fd, &tmp, 4, &junk_errno,
				 &junk_str);
	    for (i = 0; i < out_count; i++) {
		tmp = htonl(out_list[i].length);
		ffs_server_write_func(fd, &tmp, 4,
				     &junk_errno, &junk_str);
		ffs_server_write_func(fd, out_list[i].value,
				     out_list[i].length,
				     &junk_errno, &junk_str);
	    }
	    free(out_list);
	}
	break;
    case 'd':
	LOG(fs, "	Dump known hosts");
	{
	    int junk_errno;
	    char *junk_str;
	    int out_count = host_count, i;
	    int tmp;

	    tmp = htonl(out_count);
	    ffs_server_write_func(fd, &tmp, 4, &junk_errno,
				 &junk_str);
	    for (i = 0; i < out_count; i++) {
		char *time_str;
		ffs_server_write_func(fd, &hostlist[i].ip_addr, 4,
				     &junk_errno, &junk_str);
		time_str = (char *) ctime(&hostlist[i].intro_time);
		tmp = htonl(strlen(time_str) + 1);
		ffs_server_write_func(fd, &tmp, 4, &junk_errno,
				     &junk_str);
		ffs_server_write_func(fd, time_str, strlen(time_str) + 1,
				     &junk_errno, &junk_str);
	    }
	}
	break;
    case 'c':
	LOG(fs, "	Command interface");
	{
	    struct sockaddr_in in_sock;
	    struct sockaddr *sock_addr = (struct sockaddr *) &in_sock;
	    int sock_len = sizeof(in_sock);
	    int junk_errno;
	    char *junk_str;
	    int local_connection = 1;

	    memset(&in_sock, 0, sizeof(in_sock));
	    if (getpeername((int) (long) fd, sock_addr, (socklen_t*)&sock_len) != 0) {
		if (fs->stdout_verbose) {
		    printf("Rejecting command connection from non peer host\n");
		}
		LOG(fs, "Close on non peer host");
		FSClient_close(fsc);
		break;
	    }
	    if (my_IP_addr != htonl(in_sock.sin_addr.s_addr)) {
		local_connection = 0;
	    }
	    if (LogAtomicRead(fs, fd, &next_action, (long)1) != 1) {
		LOG(fs, "Close on atomic read fail c1");
		FSClient_close(fsc);
		return;
	    }
	    switch (next_action) {
	    case 'N':
		if (fs->stdout_verbose) {
		    printf("NO command action\n");
		}
		break;
	    /* case 'X': */
	    /* 	if (fs->stdout_verbose) { */
	    /* 	    printf("EXIT command action\n"); */
	    /* 	} */
	    /* 	exit(0); */
	    /* 	break; */
	    case 'D':
		if (fs->stdout_verbose) {
		    printf("dump formats to file\n");
		}
		dump_formats_to_file(fs);
		break;
	    case 'R':
		if (!local_connection) {
		    if (fs->stdout_verbose) {
			printf("Rejecting command connection not from this host\n");
		    }
		    LOG(fs, "Close on command reject 1");
		    FSClient_close(fsc);
		    break;
		}
		if (fs->stdout_verbose) {
		    printf("read formats from file\n");
		}
		read_formats_from_file(fs);
		break;
	    case 'S':
		if (fs->stdout_verbose) {
		    printf("send statistics to client\n");
		}
		send_stats(fsc);
		break;
	    case 'p':
		if (fs->stdout_verbose) {
		    printf("ping from client\n");
		}
		/* send it back */
		ffs_server_write_func(fd, &next_action, 1, &junk_errno,
				     &junk_str);
		break;
	    case 'P':
		if (fs->stdout_verbose) {
		    printf("ping/request for ID from client\n");
		}
		/* send it back */
		ffs_server_write_func(fd, &fs->format_server_identifier, 4, &junk_errno,
				     &junk_str);
		break;
	    }
	}
	break;
    default:
	LOG(fs, "Unknown request from fd %d- %s", fd, fsc->hostname);
	printf("Unknown request %c, %x\n", next_action, (unsigned char) next_action);
	FSClient_close(fsc);
	break;
    }
}


/* 
 * Create and initialize a connection data structure.
 */
static FSClient
FSClient_create(format_server fs)
{
    FSClient fsc = (FSClient) malloc(sizeof(*fsc));
    fsc->fd = NULL;
    fsc->hostname = NULL;
    fsc->port = 0;
    fsc->provisional = 0;
    fsc->fs = fs;
    fsc->version = 0;
    fsc->created = time(NULL);
    fsc->input_bytes = 0;
    fsc->output_bytes = 0;
    fsc->formats_registered = 0;
    fsc->formats_fetched = 0;
    fsc->key = NULL;
    fsc->key_len = 0;
    return fsc;
}

static void
get_qual_hostname(char *buf, int len)
{
    struct hostent *host = NULL;

    char *network_string = getenv("CM_NETWORK");
    char *hostname_string = getenv("CERCS_HOSTNAME");
    if (hostname_string != NULL) {
	strncpy(buf, hostname_string, len);
	return;
    }
    buf[0] = 0;
    gethostname(buf, len);
    buf[len - 1] = '\0';
    if (memchr(buf, '.', strlen(buf)) == NULL) {
	/* no dots, probably not fully qualified */
#ifdef HAVE_GETDOMAINNAME
	int end = strlen(buf);
	buf[end] = '.';
	if (getdomainname((&buf[end]) + 1, len - strlen(buf)) == -1) {
	    buf[end+1] = 0;
	}
	if (buf[end + 1] == 0) {
	    struct hostent *tmp;
	    buf[end] = 0;
	    /* getdomainname was useless, hope that gethostbyname helps */
	    tmp = gethostbyname(buf);
	    if (tmp) {
		strncpy(buf, tmp->h_name, len);
	    }
	}
#else
	{
	    /* no getdomainname, hope that gethostbyname will help */
	    struct hostent *tmp;
	    tmp = gethostbyname(buf);
	    if (tmp) {
		strncpy(buf, tmp->h_name, len);
	    }
	}
#endif
	buf[len - 1] = '\0';
    }
    if (memchr(buf, '.', strlen(buf)) == NULL) {
	/* useless hostname if it's not fully qualified */
	buf[0] = 0;
    }
    if ((buf[0] != 0) && ((host = gethostbyname(buf)) == NULL)) {
	/* useless hostname if we can't translate it */
	buf[0] = 0;
    }
    if (host != NULL) {
	char **p;
	int good_addr = 0;
	for (p = host->h_addr_list; *p != 0; p++) {
	    struct in_addr *in = *(struct in_addr **) p;
	    if (ntohl(in->s_addr) != INADDR_LOOPBACK) {
		good_addr++;
	    }
	}
	if (good_addr == 0) {
	    /* 
	     * even a fully qualifiedhostname that doesn't get us a valid
	     * IP addr is useless
	     */
	    buf[0] = 0;
	}
    }
    if (buf[0] == 0) {
	/* bloody hell, what do you have to do? */
	struct in_addr IP;
	IP.s_addr = htonl(get_self_ip_addr());
	host = gethostbyaddr((char *) &IP, sizeof(IP), AF_INET);
	if (host != NULL) {
	    strncpy(buf, host->h_name, len);
	}
    }
    if (network_string != NULL) {
	int name_len = strlen(buf) + 2 + strlen(network_string);
	char *new_name_str = malloc(name_len);
	char *first_dot = strchr(buf, '.');

	/* stick the CM_NETWORK value in there */
	memset(new_name_str, 0, name_len);
	*first_dot = 0;
	first_dot++;
	sprintf(new_name_str, "%s%s.%s", buf, network_string, first_dot);
	if (gethostbyname(new_name_str) != NULL) {
	    /* host has no NETWORK interface */
	    strcpy(buf, new_name_str);
	}
	free(new_name_str);
    }

    if (((host = gethostbyname(buf)) == NULL) ||
	(memchr(buf, '.', strlen(buf)) == NULL)) {
	/* just use the bloody IP address */
	struct in_addr IP;
	IP.s_addr = htonl(get_self_ip_addr());
	if (IP.s_addr != 0) {
	    struct in_addr ip;
	    char *tmp;
	    ip.s_addr = htonl(get_self_ip_addr());
	    tmp = inet_ntoa(ip);
	    strncpy(buf, tmp, len);
	} else {
	    static int warn_once = 0;
	    if (warn_once == 0) {
		warn_once++;
		printf("Attempts to establish your fully qualified hostname, or indeed any\nuseful network name, have failed horribly.  Sorry.\n");
	    }
	    strncpy(buf, "Unknown", len);
	}
    }
}

/* 
 * Create and initialize a data exchange.
 */
static format_server
format_server_create()
{
    int j;
    struct hostent *host;
    char buf[MAXHOSTNAMELEN];
    format_server fs = (format_server) malloc(sizeof(*fs));
    char *format_server_file_dir = NULL;
    int max_fd = FD_SETSIZE;	/* returns process fd table size */
    long seed;
    if (fs == NULL)
	return NULL;

/*    if (ffs_sockets_init_func != NULL) ffs_sockets_init_func();*/

    my_pid = (int) getpid();
    
    gethostname(buf, MAXHOSTNAMELEN);
    host = gethostbyname(buf);
    if (host != NULL) {
	my_IP_addr = htonl(*((int *) host->h_addr_list[0]));
    } else {
	fprintf(stderr, "Warning, you have no useful hostname.  (I.E. the hostname\n");
	fprintf(stderr, " provided by gethostbyname() cannot be translated to an IP address).\n");
	fprintf(stderr, " Using the loopback address 127.0.0.1 as IP.\n");

	my_IP_addr = INADDR_LOOPBACK;
    }
    seed = time(NULL);
    seed ^= my_IP_addr;
    srand48(seed);
    fs->format_server_identifier = lrand48();
    while ((fs->format_server_identifier == 0) || (fs->format_server_identifier == -1)) {
	/* neither 0 nor -1 are acceptable identifiers, exclude them. */
	fs->format_server_identifier = lrand48();
    }
    if (getenv("FORMAT_SERVER_PWD") != 0) {
	format_server_file_dir = malloc(PATH_MAX + 1);
	format_server_file_dir = getcwd(format_server_file_dir, PATH_MAX);
    } else {
	format_server_file_dir = strdup("/tmp");
    }

    sprintf(format_server_log, "%s/format_server_log",
	    format_server_file_dir);
    sprintf(FS_restart_file, "%s/FS_restart_file",
	    format_server_file_dir);

    /* initialize data structs */
    fs->port = -1;
    fs->hostname = NULL;
    fs->conn_sock_inet = (void *) -1;
    fs->conn_sock2_inet = (void *) -1;
    fs->usock_name = NULL;
    fs->start_time = time(NULL);
    fs->stdout_verbose = 0;
    pthread_mutex_init(&fs->lock, NULL);
    pthread_mutex_init(&fs->log_lock, NULL);
    for (j = 0; j < 1; j++) {
	fs->lists[j] = NULL;
    }

    fs->portCount = 0;
    fs->ports = (FSClient *) malloc(sizeof(FSClient) * max_fd);
    memset((char *) fs->ports, 0, sizeof(FSClient) * max_fd);

    fs->timestamp = (time_t *) malloc(sizeof(time_t) * max_fd);
    memset((char *) fs->timestamp, 0, sizeof(FSClient) * max_fd);

#ifdef __NVCOMPILER
#pragma diag_suppress 550
#endif
    FD_ZERO(&fs->fdset);
#ifdef __NVCOMPILER
#pragma diag_default 550
#endif
    fs->data_buffer = (char *) malloc(1);
    fs->buffer_size = 1;
    fs->proxy_context_to_master = NULL;
    fs->pending_unregistered = 0;
    fs->fork_threads = 1;
    /* 
     * ignore SIGPIPE's  (these pop up when ports die.  we catch the
     * failed writes)
     */
#ifdef SIGPIPE
    signal(SIGPIPE, SIG_IGN);
#endif

    return fs;
}


static int
server_listen(int port)
{
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    int length;
    int conn_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (conn_sock < 0) {

	fprintf(stderr, "Cannot open INET socket %d\n", conn_sock);
	return -1;
    }
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = htons((unsigned short) port);
    if (setsockopt((long) conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	fprintf(stderr, "Failed to set 1REUSEADDR on INET socket\n");
	return -1;
    }
    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
	     sizeof sock_addr) < 0) {
	fprintf(stderr, "Cannot bind INET socket, port %d, errno %d\n", port, errno);
	return -1;
    }
    sock_opt_val = 1;
    if (setsockopt(conn_sock, SOL_SOCKET, SO_REUSEADDR, (char *) &sock_opt_val,
		   sizeof(sock_opt_val)) != 0) {
	perror("Failed to set 2REUSEADDR on INET socket");
	return -1;
    }
    length = sizeof sock_addr;
    if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, (socklen_t*)&length) < 0) {
	fprintf(stderr, "Cannot get socket name\n");
	return -1;
    }

    /* begin listening for conns and set the backlog */
    if (listen(conn_sock, FD_SETSIZE)) {
	fprintf(stderr, "listen failed\n");
	return -1;
    }
    return conn_sock;
}

/* 
 * Create an IP socket for connection from other data exchanges.
 */
static int
format_server_listen(format_server fs, int port_num)
{
    char buf[256];
    int conn_sock;

    if ((long) fs->conn_sock_inet >= 0) {
	fprintf(stderr, "INET connection socket already open %lx\n",
		(long) fs->conn_sock_inet);
	return -1;
    }
    conn_sock = server_listen(port_num);
    fs->port = port_num;

    /* record the server's hostname */
    get_qual_hostname(buf, sizeof(buf));
    fs->hostname = (char *) strdup(buf);

    fs->conn_sock_inet = (void *) (long) conn_sock;
    if (conn_sock > FD_SETSIZE) {
	fprintf(stderr, "Internal Error, stupid WINSOCK large FD bug.\n");
	fprintf(stderr, "Increase FD_SETSIZE.  Item not added to fdset.\n");
    }
    FD_SET(conn_sock, &fs->fdset);

    if (getenv("FORMAT_SERVER_PORT2") != NULL) {
	char *port_string = getenv("FORMAT_SERVER_PORT2");
	int tmp_port, conn_sock2;
	if (sscanf(port_string, "%d", &tmp_port) != 1) {
	    printf("FORMAT_SERVER_PORT2 spec \"%s\" not understood.\n", 
		   port_string);
	    return 0;
	}
	conn_sock2 = server_listen(tmp_port);
	fs->conn_sock2_inet = (void*)(long)conn_sock2;
	FD_SET(conn_sock2, &fs->fdset);
    }

    return 0;
}

static void out_domain_rejection(int fd, FSClient fsc);
static int verify_magic_cookie (FSClient fsc);

static void *
client_handler_thread(void * fscv)
{
    FSClient fsc = (FSClient)fscv;
    server_read_header(fsc);
    while(1) {
	format_server_handle_data(fsc->fs, fsc);
    }
}

/* 
 * Accept socket connection from other data exchs.  If conn_sock is -1,
 * check if there is connection request at conn_sock_inet or
 * conn_sock_unix.
 */
extern FSClient
format_server_accept_conn_sock(format_server fs, void *conn_sock)
{
    int sock;
    FSClient fsc;
    int delay_value = 1;

    struct linger linger_val;

    linger_val.l_onoff = 1;
    linger_val.l_linger = 60;
    if (fs == NULL)
	return NULL;
    if ((long) conn_sock == -1) {
	fd_set fds;
	struct timeval timeout;
#ifdef __NVCOMPILER
#pragma diag_suppress 550
#endif
	FD_ZERO(&fds);
#ifdef __NVCOMPILER
#pragma diag_default 550
#endif
	LOG(fs, "minus 1 variation");
	if ((long) fs->conn_sock_inet >= 0) {
	    FD_SET((unsigned long) fs->conn_sock_inet, &fds);
	    if ((long) fs->conn_sock_inet > FD_SETSIZE) {
		fprintf(stderr, "Internal Error, stupid WINSOCK large FD bug.\n");
		fprintf(stderr, "Increase FD_SETSIZE.  Item not added to fdset.\n");
	    }
	}
	timeout.tv_sec = 0L;
	timeout.tv_usec = 0L;

	if (select(FD_SETSIZE, &fds, (fd_set *) NULL, (fd_set *) NULL,
		   &timeout) > 0) {
	    if (FD_ISSET((int) (long) fs->conn_sock_inet, &fds)) {
		format_server_accept_conn_sock(fs, fs->conn_sock_inet);
	    }
	}
	return NULL;
    }
    if ((sock = accept((unsigned int) (long) conn_sock, (struct sockaddr *) 0, (socklen_t*)(int *) 0)) < 0) {
	LOG(fs, "accept failed");

	timeout_oldest_conn(fs);

	if ((sock = accept((unsigned int) (long) conn_sock, (struct sockaddr *) 0, (socklen_t*)(int *) 0)) < 0) {
	    LOG(fs, "accept failed twice");
	    fprintf(stderr, "Cannot accept socket connection\n");
	    return NULL;
	}
    }
    LOG(fs, "Accepting fd %d", sock);
    if (setsockopt(sock, SOL_SOCKET, SO_LINGER, (char *) &linger_val,
		   sizeof(struct linger)) != 0) {
	perror("set SO_LINGER");
	return NULL;
    }
    setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));

    fsc = FSClient_create(fs);
    fsc->fd = (void *) (long) sock;

    fs->ports[sock] = fsc;
    fs->timestamp[sock] = time(NULL);
    fs->portCount++;

    if (fs->stdout_verbose) {
	printf("FSxchange_accept_conn_sock successful\n");
    }
    if (fsc->hostname == NULL)
	fsc->hostname = strdup("none");
    if (fs->stdout_verbose) {
	printf("fsc->fd = %d, port = %d, hostname = %s\n",
	       (int) (long) fsc->fd, fsc->port,
	       fsc->hostname);
    }
    total_client_count++;
    if (!fs->fork_threads) {
	FD_SET(sock, &fs->fdset);
	if (sock > FD_SETSIZE) {
	    fprintf(stderr, "Internal Error, stupid WINSOCK large FD bug.\n");
	    fprintf(stderr, "Increase FD_SETSIZE.  Item not added to fdset.\n");
	}
	server_read_header(fsc);
	if (!verify_magic_cookie(fsc)) {
	    out_domain_rejection(sock, fsc);
	} else {
	    if (fs->stdout_verbose) {
		printf("client verified by key authentication\n");
	    }
	}
    } else {
	pthread_create(&fsc->handler_thread, NULL, client_handler_thread, fsc);
	pthread_detach(fsc->handler_thread);
    }
    return fsc;
}

static time_t
get_time_for_host(struct in_addr ip_addr, char *hostname)
{
    int i;
    if (hostlist == NULL) {
	hostlist = malloc(sizeof(struct host_info));
    }
    for (i = 0; i < host_count; i++) {
	if (memcmp(&ip_addr, &hostlist[i].ip_addr, sizeof(struct in_addr)) == 0) {
	    return hostlist[i].intro_time;
	}
    }
    hostlist = realloc(hostlist, sizeof(struct host_info) * (host_count + 1));
    memcpy(&hostlist[host_count].ip_addr, &ip_addr, sizeof(struct in_addr));
    hostlist[host_count].intro_time = time(NULL);
    host_count++;
    return hostlist[host_count - 1].intro_time;
}

#ifndef FORMAT_SERVICE_DOMAIN
#define FORMAT_SERVICE_DOMAIN ""
#endif
static char *postfix_string = FORMAT_SERVICE_DOMAIN;
static char **postfix_list = NULL;

/* 3 days of free use of our daemons for every host */
#define GRACE_PERIOD_SEC 60 * 60 * 24 * 3

static void
out_domain_rejection(int fd, FSClient fsc)
{
    struct sockaddr sock_addr;
    int sock_len = sizeof(sock_addr);
    time_t intro_time;
    int i;
    int postfix_count = 1;

    if (postfix_list == NULL) {
	char *tmp = NULL;
	if ((tmp = getenv("FORMAT_SERVICE_DOMAIN")) == NULL) {
	  tmp = postfix_string = FORMAT_SERVICE_DOMAIN;
	}
	if (tmp[0] == '"') tmp++;
	if (postfix_string[strlen(postfix_string)-1] == '"') postfix_string[strlen(postfix_string)-1] = 0;
	postfix_string = tmp;
	while ((tmp = strchr(tmp, ':')) != NULL) {
	    postfix_count++;
	    tmp++;
	}
	postfix_list = malloc(sizeof(char *) * (postfix_count + 1));
	postfix_list[0] = strdup(postfix_string);
	for (i = 1; i < postfix_count; i++) {
	    postfix_list[i] = strchr(postfix_list[i - 1], ':');
	    if (postfix_list[i] != NULL) {
		*postfix_list[i] = 0;	/* kill : */
		postfix_list[i]++;
	    }
	}
	postfix_list[postfix_count] = NULL;
    }
	    
    if (postfix_list[0][0] == 0) {
	/* null postfix list */
	return;
    }
    memset(&sock_addr, 0, sock_len);
    if (getpeername(fd, &sock_addr, (socklen_t*)&sock_len) == 0) {
	struct hostent *host;
	if (sock_addr.sa_family == AF_INET) {
	    char *hostname;
	    struct sockaddr_in *in_sock = (struct sockaddr_in *) &sock_addr;
	    host = gethostbyaddr((char *) &in_sock->sin_addr,
				 sizeof(struct in_addr), AF_INET);
	    if (host != NULL) {
		char **tmp_list = postfix_list;
		while (*tmp_list != NULL) {
		    char *postfix = *(tmp_list++);
		    char **alias;
		    int postlen = strlen(postfix);
		    int hostlen = strlen(host->h_name);
		    /* 
		     * Check to see if our FORMAT_SERVICE_DOMAIN is a postfix
		     * of the hostname.  If not, we'll reject the host as being
		     * outside of our service domain.
		     */
		    if (hostlen > postlen) {
			if (strcasecmp(&host->h_name[hostlen - postlen], postfix) == 0) {
			    /* good host */
			    if (fsc->hostname) free(fsc->hostname);
			    fsc->hostname = strdup(host->h_name);
			    LOG(fsc->fs, "    Connection from %s", host->h_name);
			    return;
			}
		    } else {
			if (strcasecmp(host->h_name, "localhost") == 0) {
			    if (fsc->hostname) free(fsc->hostname);
			    fsc->hostname = strdup(host->h_name);
			    LOG(fsc->fs, "    Connection from %s", host->h_name);
			    return;
			}
		    }
		    for (alias = host->h_aliases; *alias != 0; alias++) {
			hostlen = strlen(*alias);
			if (hostlen > postlen) {
			    if (strcasecmp(alias[hostlen - postlen], postfix) == 0) {
				/* good host */
				if (fsc->hostname) free(fsc->hostname);
				fsc->hostname = strdup(host->h_name);
				LOG(fsc->fs, "    Connection from %s", host->h_name);
				return;
			    }
			}
		    }
		}
		hostname = host->h_name;
		if (fsc->hostname) free(fsc->hostname);
		fsc->hostname = strdup(hostname);
	    } else {
		/* no reverse DNS ??? */
		hostname = (char *) inet_ntoa(in_sock->sin_addr);
		if (fsc->hostname) free(fsc->hostname);
		fsc->hostname = strdup(hostname);
	    }

	    intro_time = get_time_for_host(in_sock->sin_addr, hostname);
	    if ((time(NULL) - intro_time) < GRACE_PERIOD_SEC) {
		if (fsc->version > 0) {
		    LOG(fsc->fs, "accepting provisional host, %s", hostname);
		    fsc->provisional = 1;
		} else {
		    LOG(fsc->fs, "rejecting old provisional host, %s", hostname);
		    close(fd);
		}
		rejected_client_count++;
		return;
	    }
	    printf("Rejecting out-domain host %s, aliases: ", hostname);
	    if (host != NULL) {
		char **alias;
		for (alias = host->h_aliases; *alias != 0; alias++) {
		    printf(" %s", *alias);
		}
	    }
	    printf("\n");
	}
    }
    /* must not be a good host, close it down ... */
    LOG(fsc->fs, "Rejecting unknown host");
    rejected_client_count++;
    close(fd);
}


static int
verify_magic_cookie (FSClient fsc)
{
    unsigned char enc_cookie[1024];
    int my_key_len;
    
    if (fsc->key_len == 0) return 0;

    my_key_len = FFS_gen_authentication(&enc_cookie[0]);

    if (my_key_len == 0) return 0;
    if (fsc->key_len != my_key_len) {
	LOG(fsc->fs, "Incoming client key length is %d, server key length is %d, rejected",
	    fsc->key_len, my_key_len);
	FSClient_close(fsc);
	return 0;
    }
    
    if (memcmp (fsc->key, enc_cookie, my_key_len) == 0) {
	return 1;
    } else {
	LOG(fsc->fs, "Incoming client key differs from server key");
	if (my_key_len >= 16) {
	    LOG(fsc->fs, " server encoded begins %02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
		enc_cookie[0],enc_cookie[1],enc_cookie[2],enc_cookie[3],
		enc_cookie[4],enc_cookie[5],enc_cookie[6],enc_cookie[7],
		enc_cookie[8],enc_cookie[9],enc_cookie[10],enc_cookie[11],
		enc_cookie[12],enc_cookie[13],enc_cookie[14],enc_cookie[15]);
	    LOG(fsc->fs, " client encoded begins %02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
		fsc->key[0],fsc->key[1],fsc->key[2],fsc->key[3],
		fsc->key[4],fsc->key[5],fsc->key[6],fsc->key[7],
		fsc->key[8],fsc->key[9],fsc->key[10],fsc->key[11],
		fsc->key[12],fsc->key[13],fsc->key[14],fsc->key[15]);
	}
	FSClient_close (fsc);
	rejected_client_count++;
	return 0;
    }
}

static void
timeout_old_conns(format_server server_fs)
{
    static int last_timeout_time = 0;
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
	if (server_fs->ports[i] != NULL) {
	    if (now - server_fs->timestamp[i] > CONN_TIMEOUT_INTERVAL) {
		LOG(server_fs, "Timeout -> fd %d", i);
		FSClient_force_close(server_fs->ports[i]);
	    }
	}
    }
}

static void
timeout_oldest_conn(format_server server_fs)
{
    int oldest_time = 0;
    int oldest_index = 0;
    int i;

    for (i = 0; i < FD_SETSIZE; i++) {
	if (server_fs->ports[i] != NULL) {
	    if (oldest_time == 0) {
		oldest_time = server_fs->timestamp[i];
		oldest_index = i;
	    } else if (oldest_time > server_fs->timestamp[i]) {
		oldest_time = server_fs->timestamp[i];
		oldest_index = i;
	    }
	}
    }
    if (oldest_time != 0) {
	LOG(server_fs, "Closing oldest -> fd %d", i);
	FSClient_force_close(server_fs->ports[oldest_index]);
    }
}

static
void 
dump_stats_to_log(format_server server_fs)
{
    int i, fsc_count = 0;
    time_t now = time(NULL);
    fprintf(log, "\n\n");
    fprintf(log, "Currently registered formats  -- %d\n", server_format_count);
    fprintf(log, "Data handle count -- %d\n", select_handle_count);
    fprintf(log, "Test count -- %d\n", test_count);
    fprintf(log, "Registration count -- %d\n", registration_count);
    fprintf(log, "Format Get count -- %d\n", get_count);
    fprintf(log, "\n format_server portCount = %d\n", server_fs->portCount);
    for (i = 0; i < FD_SETSIZE; i++) {
	if (server_fs->ports[i] != NULL) {
	    fsc_count++;
	    fprintf(log, "	fd %d - last use %s                                  age %ld\n", i,
		    ctime(&server_fs->timestamp[i]),
		    (long) now - server_fs->timestamp[i]);
	}
    }
    fprintf(log, " format_server client count = %d\n", fsc_count);
}

static FMFormat
register_server_format(format_server fs, FMContext iocontext, FMStructDescList str_list)
{
    FMFormat ioformat = register_data_format(iocontext, str_list);
    IOFormatRep iofr = malloc(sizeof(*iofr));
    int junk;
    iofr->server_format_rep = (format_rep) get_server_rep_FMformat(ioformat, &junk);
    iofr->server_ID.length = 0;
    iofr->server_ID.value = NULL;
    iofr = find_format(fs, (FSClient)NULL, iofr, /* new format */1, 
			   /* version */ 2);
    add_format_to_iofile(iocontext, ioformat,
			     iofr->server_ID.length,
			     iofr->server_ID.value, -1);
    return ioformat;
}

typedef struct client_stats {
    char *hostname;
    int byte_order_different;
    int marked_for_provisional_use;
    int FS_protocol_version;
    char *connection_initiation_time;
    int bytes_from_client;
    int bytes_to_client;
    int formats_registered;
    int formats_fetched;
} *client_stats;

typedef struct FS_stats {
    char *hostname;
    char *up_since;
    int format_count;
    int test_count;
    int registration_count;
    int fetch_count;
    int total_client_count;
    int rejected_client_count;
    int current_client_count;
    client_stats clients;
}*FS_stats;

FMField client_field_list[] = {
    {"hostname", "string", sizeof(char*), FMOffset(client_stats, hostname)},
    {"byte_order_different", "boolean", sizeof(int), 
     FMOffset(client_stats, byte_order_different)},
    {"marked_for_provisional_use", "boolean", sizeof(int), 
     FMOffset(client_stats, marked_for_provisional_use)},
    {"format_server_protocol_version", "integer", sizeof(int),
     FMOffset(client_stats, FS_protocol_version)},
    {"connection_initiation_time", "string", sizeof(char*),
     FMOffset(client_stats, connection_initiation_time)},
    {"bytes_from_client", "integer", sizeof(int),
     FMOffset(client_stats, bytes_from_client)},
    {"bytes_to_client", "integer", sizeof(int),
     FMOffset(client_stats, bytes_to_client)},
    {"formats_registered", "integer", sizeof(int),
     FMOffset(client_stats, formats_registered)},
    {"formats_fetched", "integer", sizeof(int),
     FMOffset(client_stats, formats_fetched)},
    {NULL, NULL, 0, 0}};

static char *client_xml_markup = 
"Stats for client : \n\
	Client Hostname  : <FFS:data field_name=hostname>\n\
	Byte order different : <FFS:data field_name=byte_order_different>\n\
	Provisional Use : <FFS:data field_name=marked_for_provisional_use>\n\
	Server Protocol Version  : <FFS:data field_name=format_server_protocol_version>\n\
	Connected Since : <FFS:data field_name=connection_initiation_time>\n\
	Bytes received from client : <FFS:data field_name=bytes_from_client>\n\
	Bytes sent to client : <FFS:data field_name=bytes_to_client>\n\
	Number of Formats Registered : <FFS:data field_name=formats_registered>\n\
	Number of Formats Fetched : <FFS:data field_name=formats_fetched>\n";

FMField stats_field_list[] = {
    {"server_hostname", "string", sizeof(char*), 
     FMOffset(FS_stats, hostname)},
    {"up_since", "string", sizeof(char*), FMOffset(FS_stats, up_since)},
    {"format_count", "integer", sizeof(int), 
     FMOffset(FS_stats, format_count)},
    {"test_count", "integer", sizeof(int), 
     FMOffset(FS_stats, test_count)},
    {"registration_count", "integer", sizeof(int), 
     FMOffset(FS_stats, registration_count)},
    {"fetch_count", "integer", sizeof(int), 
     FMOffset(FS_stats, fetch_count)},
    {"total_client_count", "integer", sizeof(int), 
     FMOffset(FS_stats, total_client_count)},
    {"rejected_client_count", "integer", sizeof(int), 
     FMOffset(FS_stats, rejected_client_count)},
    {"current_client_count", "integer", sizeof(int), 
     FMOffset(FS_stats, current_client_count)},
    {"clients", "client_stats[current_client_count]", 
     sizeof(struct client_stats), FMOffset(FS_stats, clients)},
    {NULL, NULL, 0, 0}};

static char *stats_xml_markup = 
"Statistics for FFS format server running on <FFS:data field_name=server_hostname>\n\
    Server Up Since : <FFS:data field_name=up_since>\n\
    Known Format Count : <FFS:data field_name=format_count>\n\
    Test Char Count : <FFS:data field_name=test_count>\n\
    Format Registration Count : <FFS:data field_name=registration_count>\n\
    Format Fetch Count : <FFS:data field_name=fetch_count>\n\
    Client Count (since initiation) : <FFS:data field_name=total_client_count>\n\
    Clients Rejected : <FFS:data field_name=rejected_client_count>\n\
    Current Client Count : <FFS:data field_name=current_client_count>\n\
    <FFS:data field_name=clients>";

static void
send_stats(FSClient fsc)
{
    static FMFormat stats_format = NULL;
    char *stats_block = NULL;
    int stats_block_len = 0, tmp, rep_len, id_len;
    struct FS_stats stats;
    int junk_errno;
    char *junk_str, *start;
    int i, fsc_count = 0;
    char *server_rep;
    char *server_ID;

    if (stats_context == NULL) {
	FMStructDescRec str_list[3];
	FMOptInfo client_opt_info[2];
	FMOptInfo stats_opt_info[2];
	client_opt_info[0].info_type = XML_OPT_INFO;
	client_opt_info[0].info_len = strlen(client_xml_markup);
	client_opt_info[0].info_block = client_xml_markup;
	client_opt_info[1].info_len = 0;
	client_opt_info[1].info_block = NULL;
	stats_opt_info[0].info_type = XML_OPT_INFO;
	stats_opt_info[0].info_len = strlen(stats_xml_markup);
	stats_opt_info[0].info_block = stats_xml_markup;
	stats_opt_info[1].info_len = 0;
	stats_opt_info[1].info_block = NULL;
	str_list[0].format_name = "format_server stats";
	str_list[0].field_list = stats_field_list;
	str_list[0].struct_size = sizeof(struct FS_stats);
	str_list[0].opt_info = stats_opt_info;
	str_list[1].format_name = "client_stats";
	str_list[1].field_list = client_field_list;
	str_list[1].struct_size = sizeof(struct client_stats);
	str_list[1].opt_info = client_opt_info;;
	str_list[2].format_name = NULL;

	stats_context = create_local_FMcontext();
	
	stats_format = register_server_format(fsc->fs, 
					      (FMContext) stats_context, 
					      str_list);
    }
    start = strdup(ctime(&fsc->fs->start_time));
    start[strlen(start) - 1] = 0;  /* stomp on <CR> */
    stats.hostname = fsc->fs->hostname;
    stats.up_since = start;
    stats.format_count = server_format_count;
    stats.test_count = test_count;
    stats.registration_count = registration_count;
    stats.fetch_count = get_count;
    stats.total_client_count = total_client_count;
    stats.rejected_client_count = rejected_client_count;
    for (i = 0; i < FD_SETSIZE; i++) {
	if (fsc->fs->ports[i] != NULL) {
	    fsc_count++;
	}
    }
    stats.current_client_count = fsc_count;
    stats.clients = malloc(sizeof(struct client_stats) * fsc_count);
    fsc_count = 0;
    for (i = 0; i < FD_SETSIZE; i++) {
	if (fsc->fs->ports[i] != NULL) {
	    char *client_start;
	    stats.clients[fsc_count].hostname = 
		fsc->fs->ports[i]->hostname;
	    stats.clients[fsc_count].byte_order_different = 
		fsc->fs->ports[i]->byte_reversal;
	    stats.clients[fsc_count].marked_for_provisional_use = 
		fsc->fs->ports[i]->provisional;
	    stats.clients[fsc_count].FS_protocol_version = 
		fsc->fs->ports[i]->version;
	    client_start = 
		strdup(ctime((time_t*)&fsc->fs->ports[i]->created));
	    client_start[strlen(start) - 1] = 0;  /* stomp on <CR> */
	    stats.clients[fsc_count].connection_initiation_time = 
		client_start;
	    stats.clients[fsc_count].bytes_from_client =
		fsc->fs->ports[i]->input_bytes;
	    stats.clients[fsc_count].bytes_to_client =
		fsc->fs->ports[i]->output_bytes;
	    stats.clients[fsc_count].formats_registered =
		fsc->fs->ports[i]->formats_registered;
	    stats.clients[fsc_count].formats_fetched =
		fsc->fs->ports[i]->formats_fetched;
	    fsc_count++;
	}
    }
/*    stats_block = encode_IOcontext_buffer(stats_context, stats_format,
      &stats, &stats_block_len);*/

    server_rep = 
	get_server_rep_FMformat(stats_format, &rep_len);
    server_ID = get_server_ID_FMformat(stats_format, &id_len);

    for (i = 0; i < fsc_count; i++) {
	free(stats.clients[i].connection_initiation_time);
    }
    free(start);
    tmp = htonl(id_len);
    ffs_server_write_func(fsc->fd, &tmp, 4, &junk_errno,
			 &junk_str);
    ffs_server_write_func(fsc->fd, server_ID, id_len,
			 &junk_errno, &junk_str);
    tmp = htonl(rep_len);
    ffs_server_write_func(fsc->fd, &tmp, 4, &junk_errno,
			 &junk_str);
    ffs_server_write_func(fsc->fd, server_rep, rep_len,
			 &junk_errno, &junk_str);
    tmp = htonl(stats_block_len);
    ffs_server_write_func(fsc->fd, &tmp, 4, &junk_errno,
			 &junk_str);
    ffs_server_write_func(fsc->fd, stats_block, stats_block_len,
			 &junk_errno, &junk_str);
}
