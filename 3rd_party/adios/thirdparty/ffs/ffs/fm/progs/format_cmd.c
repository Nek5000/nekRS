
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include <string.h>
#include "fm.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define sleep(x) Sleep(1000*x)
#else
#include <signal.h>
#include <netinet/in.h>
#endif
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif

#include "fm_internal.h"

extern int (*establish_server_connection_ptr)(FMContext fmc, action_t action);

static int
serverAtomicWrite(void* fd, void *buffer, int length)
{
     char *junk_result_str;
     int junk_errno;
     return ffs_server_write_func(fd, buffer, length, &junk_errno,
				&junk_result_str);
}

static void
usage_exit()
{
    printf("Usage:  format_cmd {checkpoint,restart,stats,ping}\n");
    exit(1);
}

static char *format_server_host = "\"not set\"";

static void
die_with_error(int sig)
{
    fprintf(stderr, "Format server %s not responding.  Timeout.\n", 
	    format_server_host);
    exit(1);
}

int
main(int argc, char **argv)
{
    FMContext context = create_FMcontext();
    char format_command_char = 'c';
    char format_dump_char = 'D';
    char format_read_char = 'R';
    char format_stats_char = 'S';
    char format_ping_char = 'p';
    char format_action_char;

/*    if (os_sockets_init_func != NULL) os_sockets_init_func();*/

    if (argc == 2) {
	if (strncmp(argv[1], "checkpoint", strlen(argv[1])) == 0) {
	    format_action_char = format_dump_char;
	} else if (strncmp(argv[1], "restart", strlen(argv[1])) == 0) {
	    format_action_char = format_read_char;
	} else if (strncmp(argv[1], "stats", strlen(argv[1])) == 0) {
	    format_action_char = format_stats_char;
	} else if (strncmp(argv[1], "ping", strlen(argv[1])) == 0) {
	    format_action_char = format_ping_char;
#ifndef HAVE_WINDOWS_H
	    signal(SIGALRM, die_with_error);
	    alarm(10);
#endif
	} else {
	    usage_exit();
	}
    } else {
	usage_exit();
    }

    /* 
     *  this mirrors code in server_acts.c, reproduced here to avoid
     *  exporting information
     */
    format_server_host = getenv("FORMAT_SERVER_HOST");
    if (format_server_host == NULL) {
	format_server_host = FORMAT_SERVER_HOST;	/* from configure */
    }

    if (establish_server_connection_ptr(context, host_only) == 0) {
	printf("Failed to contact format server on host %s\n",
	       format_server_host);
	return 0;
    }

    if (serverAtomicWrite(context->server_fd, &format_command_char, 1) != 1) {
	fprintf(stderr, "command write failed\n");
	exit(1);
    }

    if (serverAtomicWrite(context->server_fd, &format_action_char, 1) != 1) {
	fprintf(stderr, "action write failed\n");
	exit(1);
    }

    if (format_action_char == 'S') {
	/* get stats */
	int len = 0;
	char *buffer;
	int id_len = 0;
	char *id;
	int rep_len = 0;
	char *server_rep;
	FMFormat format;
	if (serverAtomicRead(context->server_fd, &id_len, 4) != 4) {
	    fprintf(stderr, "length read failed...  Command rejected by server.\n");
	    exit(1);
	}
	id_len = ntohl(id_len);
	id = malloc(id_len);
	if (serverAtomicRead(context->server_fd, id, id_len) != len) {
	    fprintf(stderr, "stats block read failed\n");
	}
	printf("Read %d bytes\n", id_len);
	if (serverAtomicRead(context->server_fd, &rep_len, 4) != 4) {
	    fprintf(stderr, "length read failed...  Command rejected by server.\n");
	    exit(1);
	}
	rep_len = ntohl(rep_len);
	server_rep = malloc(rep_len);
	if (serverAtomicRead(context->server_fd, server_rep, rep_len) != len) {
	    fprintf(stderr, "stats block read failed\n");
	}
	printf("Read %d bytes\n", rep_len);

	load_external_format_FMcontext(context, id, rep_len, server_rep);
	if (serverAtomicRead(context->server_fd, &len, 4) != 4) {
	    fprintf(stderr, "length read failed...  Command rejected by server.\n");
	    exit(1);
	}
	len = ntohl(len);
	buffer = malloc(len);
	if (serverAtomicRead(context->server_fd, buffer, len) != len) {
	    fprintf(stderr, "stats block read failed\n");
	}
	printf("Read %d bytes\n", len);
	format = FMformat_from_ID(context, buffer);
	FMdump_encoded_data(format, buffer, -1);
    } else if (format_action_char == 'p') {
	/* ping */
	char response;
	if (serverAtomicRead(context->server_fd, &response, 1) != 1) {
	    fprintf(stderr, "Ping response read failed.\n");
	    exit(1);
	}
	if (response != 'p') {
	    printf("Format server %s failed ping response.\n",
		   format_server_host);
	    exit(1);
	}
	printf("Format server %s is responsive.\n", format_server_host);
    }
    return 0;
}	
