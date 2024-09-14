
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include <string.h>
#include "fm.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define sleep(x) Sleep(1000*x)
#else
#include <netinet/in.h>
#include <arpa/inet.h>
extern int sleep();
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

int
main(int argc, char **argv)
{
    FMContext context = create_FMcontext();
    char format_dump_char = 'D';

/*    if (os_sockets_init_func != NULL) os_sockets_init_func();*/

    if (argc > 1) format_dump_char = 'd';
    if (establish_server_connection_ptr(context, host_and_fallback) == 0) {
	printf("Failed to contact format server\n");
	return 0;
    }

    if (serverAtomicWrite(context->server_fd, &format_dump_char, 1) != 1) {
	fprintf(stderr, "dump write failed\n");
	exit(1);
    }
    if (format_dump_char == 'D') {
	/* dumping formats */
	int format_count, i;
	server_ID_type *format_list;

	if (serverAtomicRead(context->server_fd, &format_count, 4) != 4) {
	    fprintf(stderr, "dump count read failed\n");
	    exit(1);
	}
	format_count = ntohl(format_count);
	printf("There are %d formats registered\n", format_count);
	format_list = malloc(sizeof(format_list[0]) * format_count);
	for (i=0; i < format_count; i++) {
	    if (serverAtomicRead(context->server_fd, &format_list[i].length, 
				 4) != 4) {
		fprintf(stderr, "dump list read failed\n");
		exit(1);
	    }
	    format_list[i].length = ntohl(format_list[i].length);
	    format_list[i].value = malloc(format_list[i].length);
	    if (serverAtomicRead(context->server_fd, format_list[i].value,
				 format_list[i].length)
		!= format_list[i].length) {
		fprintf(stderr, "dump list read failed\n");
		exit(1);
	    }
	}
	for (i=0; i < format_count; i++) {
	    FMFormat tmp_format = FMformat_from_ID(context, format_list[i].value);
	    printf("\n\nFormat %d (server size %d) is:\n", i,
		   tmp_format->server_format_rep->format_rep_length);
	    if ( tmp_format != NULL ) {
		dump_FMFormat(tmp_format);
	    } else {
		printf("Get format failed!\n");
	    }
	}
    } else {
	/* dumping hosts */
	int host_count, i;
	if (serverAtomicRead(context->server_fd, &host_count, 4) != 4) {
	    fprintf(stderr, "dump count read failed\n");
	    exit(1);
	}
	host_count = ntohl(host_count);
	printf("There are %d foreign hosts which have used this server\n", host_count);
	for (i=0; i < host_count; i++) {
	    int ip_addr;
	    int str_len;
	    char time_str[1024];
	    char *hostname;
	    struct hostent *host;

	    if (serverAtomicRead(context->server_fd, &ip_addr, 4) != 4) {
		fprintf(stderr, "dump ip read failed\n");
		exit(1);
	    }
	    if (serverAtomicRead(context->server_fd, &str_len, 4) != 4) {
		fprintf(stderr, "dump str len read failed\n");
		exit(1);
	    }
	    str_len = ntohl(str_len);
	    if (str_len > sizeof(time_str)) {
		fprintf(stderr, "time string too long\n");
		exit(1);
	    }
	    if (serverAtomicRead(context->server_fd, &time_str[0], str_len) !=
		str_len) {
		fprintf(stderr, "dump time string read failed\n");
		exit(1);
	    }
	    host = gethostbyaddr((char*)&ip_addr, sizeof(ip_addr), AF_INET);
	    if (host != NULL) {
		hostname = host->h_name;
	    } else {
		struct in_addr IP;
		IP.s_addr = ip_addr;
		hostname = (char*)inet_ntoa(IP);
	    }
	    printf("Host %s, active since %s\n", hostname, time_str);
	}
    }
    return 0;
}	
