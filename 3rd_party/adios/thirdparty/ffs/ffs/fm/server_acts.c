
#ifdef MODULE
#ifndef __KERNEL__
#define __KERNEL__
#endif
#endif

#ifndef MODULE
#include "config.h"
#include <fcntl.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <errno.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <time.h>
#ifdef HAVE_SYS_SOCKET_H
#include <sys/socket.h>
#endif
#ifdef HAVE_SOCKLIB_H
#include "sockLib.h"
#include "hostLib.h"
#endif
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifndef HAVE_WINSOCK_H
#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#endif
#include <net/if.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#else
#include <winsock.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#include "stdio.h"
#include "assert.h"
#include "fm.h"
#include "fm_internal.h"
#include "io_interface.h"

#else
#include "config.h"
#include <linux/types.h>
#include <linux/fcntl.h>
#include <linux/slab.h>
#include <linux/string.h>
#include <linux/errno.h>
#include <linux/unistd.h>
#include <linux/time.h>
#include <linux/net.h>
#include <linux/mm.h>
#include <linux/socket.h>
#include <linux/un.h>
#include <linux/in.h>
#include <linux/uio.h>
#include <linux/fs.h>
#include <linux/artemis.h>
#include <net/sock.h>
#include <asm/uaccess.h>
#include "kernel/kpbio.h"
#include "kernel/library.h"
#include "assert.h"

#define printf         printk
#define perror         printk
#define fprintf(fmt,.. printk (args) 
#define malloc         (void *)DAllocMM
#define free(a)        DFreeMM((addrs_t)a)
#define realloc(a,b)   (void *)DReallocMM((addrs_t)a,b)
#define exit           sys_exit
#define qsort          _quicksort
#define getpid()       current->pid
#define signal(a,b)    
#define gethostname    sys_gethostname
#define gethostbyname  lookup_name
#endif

typedef enum {Raw, Authenticated} Server_Protocol;
#define CIPHER_BUF_SIZE 1024
#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

static int
fill_hostaddr(void *addr, char *hostname, Server_Protocol *protocol)
{
    char *colon = strchr(hostname, ':');
    struct hostent *host_addr;

    *protocol = Raw;
    if (colon != NULL) {
	if (strncasecmp(hostname, "auth:", 5) == 0) {
	    *protocol = Authenticated;
	    hostname = hostname + 5;
	} else {
	    fprintf(stderr, "Unknown protocol specified with FORMAT_SERVER_HOST, \"%s\"\n",
		    hostname);
	    hostname = colon + 1;
	}
    }
#ifdef MODULE
    struct hostent *host_addr;

    if ((host_addr = lookup_name(hostname)) == NULL) {
	return 0;
    }
    memcpy(addr, host_addr->h_addr, host_addr->h_length);
    DFreeMM((addrs_t)host_addr->h_name);
    DFreeMM((addrs_t)host_addr);
    return 1;
#else
    host_addr = gethostbyname(hostname);
    if (host_addr == NULL) {
	int _addr = inet_addr(hostname);
	if (_addr == -1) {
	    /* 
	     *  not translatable as a hostname or 
	     * as a dot-style string IP address
	     */
	    return 0;
	}
	assert(sizeof(int) == sizeof(struct in_addr));
	*((int*)addr) = (int)_addr;
    } else {
	memcpy(addr, host_addr->h_addr, host_addr->h_length);
    }
    return 1;
#endif
}

static int format_server_verbose = -1;
static int ffs_format_server_port = DEFAULT_FS_PORT;
static char *format_server_host = NULL;

extern void
set_format_server(char *hostname, int port)
{
    format_server_host = strdup(hostname);
    ffs_format_server_port = port;
}


extern int
server_write_header(FMContext fmc, int enc_len, unsigned char *enc_buffer);

extern int
establish_server_connection(FMContext iofile, action_t action)
{
    SOCKET sock;
    int ret;
    Server_Protocol protocol = Raw;
    int conn_is_dead = 0;
    int connected = 0;
#ifndef MODULE
    int delay_value = 1;
#else
    struct socket *socket;	
#endif
    if (ffs_sockets_init_func) ffs_sockets_init_func();

    if (iofile->server_fd != (void*)-1) {
#ifndef MODULE
	fd_set rd_set;
	struct timeval timeout;
	int ret;

#ifdef __NVCOMPILER
#pragma diag_suppress 550
#endif
	FD_ZERO(&rd_set);
#ifdef __NVCOMPILER
#pragma diag_default 550
#endif
	timeout.tv_sec = 0;
	timeout.tv_usec = 0;
	FD_SET( (int)(intptr_t)iofile->server_fd, &rd_set);
	ret = select(FD_SETSIZE, &rd_set, (fd_set*)NULL, 
		     (fd_set *) NULL, &timeout);
	if (ret == -1) {
	    if (format_server_verbose)
	      printf("Dead connection, Select return is %d, server fd is %p, errno is %d\n", ret, (void*)(intptr_t)iofile->server_fd, errno);
	    conn_is_dead = 1;
	}
#else
	int junk_errno;
	char *junk_str;
	int ret = os_server_write_func(iofile->server_fd, &test_byte, 1, 
				       &junk_errno, &junk_str);
	if (ret != 1) conn_is_dead = 1;
#endif
    }
    if ((iofile->server_fd == (void*)-1) || (conn_is_dead)) {
	/* reestablish connection, name_str is the machine name */
	struct sockaddr_in sock_addr;
	char *tmp_server_host;

	if (format_server_host == NULL) {
	    char *format_server_port = NULL;
	    format_server_host = getenv("FORMAT_SERVER_HOST");
	    format_server_port = getenv("FORMAT_SERVER_PORT");
	    if (format_server_port != NULL) {
		int tmp_port;
		if (sscanf(format_server_port, "%d", &tmp_port) != 1) {
		    printf("FORMAT_SERVER_PORT spec \"%s\" not understood.\n",
			   format_server_port);
		} else {
		    ffs_format_server_port = tmp_port;
		}
	    }
	}

	if (format_server_verbose == -1) {
	    if (getenv("FORMAT_SERVER_VERBOSE") == NULL) {
		format_server_verbose = 0;
	    } else {
		format_server_verbose = 1;
	    }
	}
	if (format_server_host == NULL) {
	    format_server_host = FORMAT_SERVER_HOST;	/* from configure */
	}
	if (format_server_verbose && conn_is_dead) {
	    printf("detected dead link to format server, restarting\n");
	}
	if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
	    fprintf(stderr, "Failed to create socket for FFS format server connection.  Not enough File Descriptors?\n");
	    return 0;
	}
	
	sock_addr.sin_family = AF_INET;
		
	tmp_server_host = format_server_host;
	if (action == local_only) {
	    tmp_server_host = "localhost";
	}
	if ((action != never_local) || (strcmp(tmp_server_host, "localhost") != 0)) {
	    if (fill_hostaddr(&sock_addr.sin_addr, tmp_server_host, 
			      &protocol) == 0) {
		if (iofile->self_server_fallback == 0) {
		    fprintf(stderr, "Unknown Host \"%s\" specified as FFS format server.\n",
			    tmp_server_host);
		}
		return 0;
	    }
	    if (format_server_verbose == 1) {
		printf("Trying connection to format server on %s ...  ",
		       tmp_server_host);
	    }
	    sock_addr.sin_port = htons(ffs_format_server_port);
	    
	    if (connect(sock, (struct sockaddr *) &sock_addr,
			sizeof sock_addr) < 0) {

		if (format_server_verbose) {
		    printf("failed\n");
		}
	    } else {
		connected++;
	    }

	    if (!connected) {
		if ((action == local_only) || (action == host_only)) return 0;
		/* fallback */

		if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
		    fprintf(stderr, "Failed to create socket for FFS format server connection.  Not enough File Descriptors?\n");
		    return 0;
		}
		
		format_server_host = "formathost.cercs.gatech.edu";
		sock_addr.sin_family = AF_INET;
		if (fill_hostaddr(&sock_addr.sin_addr, format_server_host, 
				  &protocol) == 0) {
		    if (iofile->self_server_fallback == 0) {
			fprintf(stderr, "Unknown Host \"%s\" specified as FFS format server.\n",
				format_server_host);
		    }
		    return 0;
		}
		sock_addr.sin_port = htons(ffs_format_server_port);
		if (format_server_verbose == 1) {
		    printf("Trying fallback connection to format server on %s ...  ",
			   format_server_host);
		}
		if (connect(sock, (struct sockaddr *) &sock_addr,
			    sizeof sock_addr) < 0) {
		    if (iofile->self_server_fallback == 0) {
			fprintf(stderr, "Failed to connect to primary or fallback format servers.\n");
		    }
		    return 0;
		}
	    }
	}
	if (format_server_verbose == 1) {
	    printf("succeeded\n");
	}
#ifndef MODULE
	setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
		   sizeof(delay_value));
	delay_value = 1;
	setsockopt(sock, SOL_SOCKET, SO_KEEPALIVE, (const char *)& delay_value, sizeof(delay_value));
#endif
	iofile->server_fd = (void*) (intptr_t) sock;
	
	if (protocol == Authenticated) {
	    unsigned char outbuf[CIPHER_BUF_SIZE];
	    int enc_len = FFS_gen_authentication(&outbuf[0]);
	    ret = server_write_header(iofile, enc_len, &outbuf[0]);
	    if (format_server_verbose) {
		printf(" encoded key length %d, begins %02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
		       enc_len, outbuf[0],outbuf[1],outbuf[2],outbuf[3],
		       outbuf[4],outbuf[5],outbuf[6],outbuf[7],
		       outbuf[8],outbuf[9],outbuf[10],outbuf[11],
		       outbuf[12],outbuf[13],outbuf[14],outbuf[15]);
	    }
	} else {
	    ret = server_write_header(iofile, 0, NULL);
	}	    
	
	/* 
	 * ignore SIGPIPE's  (these pop up when ports die.  we catch the 
	 * failed writes) 
	 */
#ifdef SIGPIPE
	signal(SIGPIPE, SIG_IGN);
#endif
	if (ret != 1) {
	    /* 
	     * got a server, but it was unacceptable in some way, 
	     * didn't give us a magic number or had restarted since 
	     *  last time. 
	     */
	    return 0;
	}
    
    }
    return 1;
}
	

#ifdef USE_ENCRYPTED_AUTHENTICATOR
#include <gcrypt.h>
int
FFS_gen_authentication (unsigned char *outbuf)
{
    unsigned char key[16];

    unsigned char cookie[CIPHER_BUF_SIZE];
    int got_key = 0;
    static char enc_cookie[CIPHER_BUF_SIZE];
    char in_buf[CIPHER_BUF_SIZE];
    static int enc_len = 0;
    static int first_time = 1;
    char *auth_file_name;
    struct stat statbuf;
    FILE *auth;
    int good_cookie = 0;
  
    /* initialize the cipher handle */
    gcry_cipher_hd_t handle;
    gcry_error_t err;
    gcry_check_version( "1.0.0" );
    err = gcry_cipher_open( &handle, GCRY_CIPHER_BLOWFISH, GCRY_CIPHER_MODE_CBC, 0 );

    if (gcry_err_code( err ) != GPG_ERR_NO_ERROR) {
      /*  some kind of error  */
      fprintf (stderr, "Error in setting up cipher context.");
      return 0;
    }

    if (!first_time) {
	memcpy(outbuf, &enc_cookie[0], enc_len);
	return enc_len;
    }
    first_time = 0;
    auth_file_name = getenv("FFS_AUTHENTICATION_FILE");
    if (auth_file_name == NULL) return 0;

    if (stat(auth_file_name, &statbuf) == -1) {
	fprintf(stderr, "failed to stat() FFS_AUTHENTICATION_FILE \"%s\"  ",
		auth_file_name);
	perror("error:");
	return 0;
    }
    auth = fopen(auth_file_name,"r");
    if (auth == NULL) {
	fprintf(stderr, "failed to open FFS_AUTHENTICATION_FILE \"%s\"  ",
		auth_file_name);
	perror("error:");
	return 0;
    }
    while (fgets(&in_buf[0], sizeof(in_buf) -1, auth) != NULL) {
	char *start = &in_buf[0];
	char *cr;
	in_buf[sizeof(in_buf) -1] = 0;
	cr = start + strlen(start) - 1;
	if (*cr == '\n') *cr = 0;
	start += strspn(start, " \t");
	if (strlen(start) == 0) continue;
	if (*start == '#') continue;
	if (strncmp(in_buf, "COOKIE=", 7) == 0) {
	    /* got the cookie line */
	    char *val = &in_buf[7];
	    val += strspn(val, " \t");
	    if (val[0] == '"') {
		char *end = strrchr(++val, '"');
		if (end) *end = 0;
	    }
	    /* cookie = malloc(strlen(val) + 1); */
	    strcpy((char*)cookie, val);
	    good_cookie++;
	} else if (strncmp(in_buf, "KEY=", 4) == 0) {
	    /* got the cookie line */
	    char *val = &in_buf[4];
	    int i;
	    val += strspn(val, " \t");
	    if (strspn(val, "0123456789abcdef") != 32) {
		fprintf(stderr, "Incorrect charcter count in FFS auth key, expected 32 hex characters, got \"%s\"\n", val);
		continue;
	    }
	    for (i=0 ; i < 16; i++) {
		unsigned int key_item;
		sscanf(&val[2*i], "%2x", &key_item);
		key[i] = key_item;
	    }
	    got_key++;
	} else {
	    printf("unknown line in FFS_AUTH_FILE : \"%s\"\n", &in_buf[0]);
	}
    }
    if (!good_cookie) {
	fprintf(stderr, "Acceptable COOKIE entry not found in FFS_AUTHORIZATION_FILE\n");
	return 0;
    }
    if (got_key == 0) {
	fprintf(stderr, "Acceptable KEY entry not found in FFS_AUTHORIZATION_FILE\n");
	return 0;
    }
    /*
     * for now encrypt the cookie
     */
    gcry_cipher_setkey( handle, &key[0], 128 );
    /* not strictly necessary to use an IV */
    /*    gcry_cipher_setiv( handle, &iv[0], 64 );  */

    err = gcry_cipher_encrypt( handle, outbuf, CIPHER_BUF_SIZE, cookie, CIPHER_BUF_SIZE );

    if (gcry_err_code( err ) != GPG_ERR_NO_ERROR) {
      /*  some kind of error  */
      fprintf (stderr, "Error in encrypting cookie.");
      return 0;
    }

    memcpy(&enc_cookie[0], outbuf, CIPHER_BUF_SIZE);

    gcry_cipher_close( handle );

    return CIPHER_BUF_SIZE;
}
#else
int 
FFS_gen_authentication (unsigned char *outbuf)
{ 
    return 0;
}
#endif


extern int (*establish_server_connection_ptr)(FMContext iofile, action_t action);
extern int ffs_fixed_format_ids;
extern INT4 FFS_self_server_IP_addr;  /* declared in io_formats.c */

#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 64
#endif

extern
FMContext
create_FMcontext()
{
    FMContext iofile;

    if (establish_server_connection_ptr == NULL) {
	/* fill in function pointer used in fm_formats.c */
	establish_server_connection_ptr = establish_server_connection;
    }
    iofile = new_FMContext();
    return iofile;
}
