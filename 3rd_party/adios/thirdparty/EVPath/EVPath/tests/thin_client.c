#include "../config.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifndef HAVE_WINDOWS_H
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#else
#include <WinSock2.h>
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define close(x) closesocket(x)
#endif
#include <string.h>
#include "ffs.h"

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;

static FMField simple_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {"short_field", "integer",
     sizeof(short), FMOffset(simple_rec_ptr, short_field)},
    {"long_field", "integer",
     sizeof(long), FMOffset(simple_rec_ptr, long_field)},
    {"double_field", "float",
     sizeof(double), FMOffset(simple_rec_ptr, double_field)},
    {"char_field", "char",
     sizeof(char), FMOffset(simple_rec_ptr, char_field)},
    {"scan_sum", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, scan_sum)},
    {NULL, NULL, 0, 0}
};

FMStructDescRec thin_formats[] = {
    {"thin_message", simple_field_list, sizeof(simple_rec), NULL},
    {NULL, NULL, 0, NULL}};

#ifndef _MSC_VER
#define SOCKET int
#endif
static SOCKET do_connection(char* host, int port);
static void generate_record (simple_rec_ptr event);

int
main(int argc, char **argv)
{
    char *remote_host;
    int remote_port;
    int stone;
    FFSFile out_connection;
    FMFormat ioformat;
    FMContext fmc;
    simple_rec rec;
    SOCKET conn;
    char *comment = strdup("Stone xxxxx");

    if ((argc != 4) || (sscanf(argv[2], "%d", &remote_port) != 1)
			|| (sscanf(argv[3], "%d", &stone) != 1)) {
	printf("Usage \"thin_client remote_host remote_port\"\n");
	exit(1);
    }
    remote_host = argv[1];

    conn = do_connection(remote_host, remote_port);

    if (conn == -1) {
	printf("Connection to %s:%d failed\n", remote_host, remote_port);
	exit(1);
    }
    out_connection = open_FFSfd((void*)(intptr_t)conn, "w");

    sprintf(comment, "Stone %d", stone);
    write_comment_FFSfile(out_connection, comment);
    fmc = FMContext_of_file(out_connection);
    ioformat = register_data_format(fmc, thin_formats);
    generate_record(&rec);
    write_FFSfile(out_connection, ioformat, &rec);
    close_FFSfile(out_connection);
    free_FFSfile(out_connection);
    free(comment);
    return 0;
}

SOCKET
do_connection(char * remote_host, int port)
{
    struct hostent *host_addr;
    struct sockaddr_in sin;

    memset((char*)&sin, 0, sizeof(sin));

    SOCKET conn = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);

    if (conn == -1) return -1;

    host_addr = gethostbyname(remote_host);
    if (!host_addr) {
	sin.sin_addr.s_addr = inet_addr(remote_host);
	if(host_addr == NULL) return -1;
    } else {
	memcpy((char*)&sin.sin_addr, host_addr->h_addr, host_addr->h_length);
    }
    sin.sin_port = (unsigned short) htons(port);
    sin.sin_family = AF_INET;
    if (connect(conn, (struct sockaddr *) &sin,
		sizeof sin) == -1) {
#ifdef WSAEWOULDBLOCK
	int err = WSAGetLastError();
	if (err != WSAEWOULDBLOCK || err != WSAEINPROGRESS) {
#endif
	    close(conn);
	    return -1;
#ifdef WSAEWOULDBLOCK
	}
#endif
    }
    return conn;
}

static
void 
generate_record(simple_rec_ptr event)
{
    long sum = 0;
    memset(event, 0, sizeof(*event));
    event->integer_field = (int) lrand48() % 100;
    sum += event->integer_field % 100;
    event->short_field = ((short) lrand48());
    sum += event->short_field % 100;
    event->long_field = ((long) lrand48());
    sum += event->long_field % 100;

    event->double_field = drand48();
    sum += ((int) (event->double_field * 100.0)) % 100;
    event->char_field = lrand48() % 128;
    sum += event->char_field;
    sum = sum % 100;
    event->scan_sum = (int) sum;
}

