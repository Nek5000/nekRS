#include "config.h"
#include <sys/types.h>
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#define FD_SETSIZE 1024
#include <windows.h>
#include <stdio.h>
#include <fcntl.h>
#include <io.h>
#include "ffs.h"
#include "io_interface.h"
#include "ffs_internal.h"

static int
nt_file_read_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    int left = length;
    DWORD iget;
    BOOL bResult;
    bResult = ReadFile(conn, (char *) buffer, length, &iget, NULL);

    if (iget == 0) {
	if (result_p) *result_p = "End of file";
	if (errno_p) *errno_p = 0;
	return 0;		/* end of file */
    } else {
	if (!bResult) {
	    return -1;
	} else {
		if(errno_p) *errno_p = 0;
	}
    }

    left = length - iget;
    while (left > 0) {
	bResult = ReadFile((HANDLE) conn, (char *) buffer + length - left,
			   left, &iget, NULL);
	if (iget == 0) {
	    if (result_p) *result_p = "End of file";
	    if (errno_p) *errno_p = 0;
	    return length - left;	/* end of file */
	} else {
	    if (!bResult) {
		    return (length - left);
		} else {
		    if (errno_p) *errno_p = 0;
		}
	}
	left -= iget;
    }
    return length;
}

static int
nt_socket_read_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    int left = length;
    int iget;
    iget = recv((unsigned int) (intptr_t)conn, (char *) buffer + length - left, left, 0);
    if (iget == 0) {
	if (result_p) *result_p = NULL;
	if (errno_p) *errno_p = 0;
	return 0;		/* No more socket data */
    } else if (iget == SOCKET_ERROR) {
	DWORD tmp = WSAGetLastError();
	if (errno_p) *errno_p = tmp;
	if ((tmp != WSAEWOULDBLOCK) &&
	    (tmp != WSAEINPROGRESS) &&
		(tmp != WSAECONNRESET)  &&
	    (tmp != WSAEINTR)) {
	    /* serious error */
	    fprintf(stderr, "WINSOCK ERROR during receive, %i on socket %p\n",
		    (int)tmp, conn);
	    return -1;
	} else {
		if (tmp == WSAECONNRESET)
			return -1;
	    if (errno_p) *errno_p = 0;
	    iget = 0;
	}
    }
    left = length - iget;
    while (left > 0) {
	iget = recv((unsigned int)(intptr_t) conn, (char *) buffer + length - left,
		    left, 0);
	if (iget == 0) {
	    if (result_p) *result_p = NULL;
	    if (errno_p) *errno_p = 0;
	    return length - left;	/* no more socket data */
	} else {
	    if (iget == SOCKET_ERROR) {
		DWORD tmp = WSAGetLastError();
		if (errno_p) *errno_p = tmp;
		if ((tmp != WSAEWOULDBLOCK) &&
		    (tmp != WSAEINPROGRESS) &&
			(tmp != WSAECONNRESET)  &&
		    (tmp != WSAEINTR)) {

		    /* serious error */
		    fprintf(stderr, "WINSOCK ERROR during receive2, %i on socket %p\n",
			    (int) tmp, conn);
		    return (length - left);
		} else {
			if (tmp == WSAECONNRESET)
				return -1;
		    if(errno_p) *errno_p = 0;
		    iget = 0;
		}
	    }
	}
	left -= iget;
    }

    return length;
}


static int
nt_file_write_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    int left = length;
    int iget = 0;
    BOOL bResult;

    while (left > 0) {
	bResult = WriteFile((HANDLE) conn, (char *) buffer + length - left, 
			    left, (unsigned long *)&iget, NULL);
	if (!bResult) {
	    DWORD tmp = GetLastError();
	    if ((tmp != WSAEWOULDBLOCK) &&
		(tmp != WSAEINPROGRESS) &&
		(tmp != WSAEINTR)) {
		/* serious error */
		return (length - left);
	    } else {
		if(errno_p) *errno_p = 0;
		iget = 0;
	    }
	}
	left -= iget;
    }
    return length;
}


static int
nt_socket_write_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    int left = length;
    int iget = 0;

    while (left > 0) {
	iget = send((unsigned int) (intptr_t)conn, (char *) buffer + length - left,
		    left, 0);
	if (iget == SOCKET_ERROR) {
	    DWORD tmp = GetLastError();
	    if (errno_p) *errno_p = tmp;
	    if ((tmp != WSAEWOULDBLOCK) &&
		(tmp != WSAEINPROGRESS) &&
		(tmp != WSAEINTR)) {
		/* serious error */
		return (length - left);
	    } else {
		if (errno_p) *errno_p = 0;
		iget = 0;
	    }
	}
	left -= iget;
    }
    return length;
}


static int
nt_close_func(conn)
void *conn;
{
    DWORD status;
    /* make sure handle exists before we close it. *otherwise -- an * * *
     * access error occurs */
    if (GetHandleInformation(conn, &status)) {
	CloseHandle((HANDLE) conn);
	return 1;
    }
    return 0;
}

static void *
nt_file_open_func(const char *path, const char *flag_str, int *input, int *output)
{
    int readfile = 0;
    int writefile = 0;
    void *file;
    long tmp_flags = (long)(intptr_t)flag_str;
    if (input) *input = 0;
    if (output) *output = 0;

    tmp_flags &= ~(O_TRUNC);
    tmp_flags &= ~(O_CREAT);

    if ((O_RDONLY == tmp_flags) ||
	(O_WRONLY == tmp_flags)) {
	 /* must be old style call */
	if (input) *input = (O_RDONLY == (long) (intptr_t)flag_str);
	if (output) *output = (O_WRONLY & (long) (intptr_t)flag_str);
    } else {
	if (strcmp(flag_str, "r") == 0) {
	    if (input) *input = TRUE;
	    readfile = 1;
	} else if (strcmp(flag_str, "w") == 0) {
	    if (output) *output = TRUE;
	    writefile = 1;
	 } else if (strcmp(flag_str, "a") == 0) {
	     if (output) *output = 1;
	     if (input) *input = 1;
	     readfile = writefile = 1;
	} else {
	    fprintf(stderr, "Open flags value not understood for file \"%s\"\n",
		    path);
	    return NULL;
	}
    }

    if (readfile) {
	file = CreateFile(path, GENERIC_READ, FILE_SHARE_READ,
		      NULL, OPEN_EXISTING, FILE_ATTRIBUTE_ARCHIVE, NULL);

    } else {
	file = CreateFile(path, GENERIC_WRITE, FILE_SHARE_READ,
		      NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_ARCHIVE, NULL);
    }
    if (file == INVALID_HANDLE_VALUE) {
	return NULL;
    } else {
	return file;
    }
}

static int
nt_file_lseek_func (void *file, size_t pos, int origin)
{
    return SetFilePointer((HANDLE)file, (long)pos, 0, origin);
}


static int
nt_socket_readv_func(conn, iov, icount, errno_p, result_p)
void *conn;
struct iovec *iov;
int icount;
int *errno_p;
char **result_p;
{

    int i = 0;
    for (; i < icount; i++) {
	if (nt_socket_read_func(conn, (void*)iov[i].iov_base, iov[i].iov_len,
				errno_p, result_p) != iov[i].iov_len) {
	    return i;
	}
    }
    return icount;
}


static int
null_file_readv_func(conn, iov, icount, errno_p, result_p)
void *conn;
struct iovec *iov;
int icount;
int *errno_p;
char **result_p;
{

    int i = 0;
    for (; i < icount; i++) {
	if (nt_file_read_func(conn, (void*)iov[i].iov_base, iov[i].iov_len, errno_p,
			      result_p) != iov[i].iov_len) {
	    return i;
	}
    }
    return icount;
}

static int
null_file_writev_func(conn, iov, icount, errno_p, result_p)
void* conn;
struct iovec* iov;
int icount;
int* errno_p;
char** result_p;
{

    int i = 0;
    for (; i < icount; i++) {
	if (nt_file_write_func(conn, (void*)iov[i].iov_base, iov[i].iov_len, errno_p,
	    result_p) != iov[i].iov_len) {
	    return i;
	}
    }
    return icount;
}

/* Winsock init stuff, ask for ver 1.1 */
static WORD wVersionRequested = MAKEWORD(2, 2);
static WSADATA wsaData;

static void
nt_socket_init_func()
{
    static int once = 0;
    if (once) return;
    once = 1;
    int nErrorStatus;
    nErrorStatus = WSAStartup(wVersionRequested, &wsaData);
    if (nErrorStatus != 0) {
	fprintf(stderr, "Could not initialize windows socket library!");
	WSACleanup();
	exit(-1);
    }
}


static int
nt_poll_func(conn)
void *conn;
{
    int fd = (int) (intptr_t) conn;
    struct timeval time;
    fd_set read_fds;
    int ret_val;

    time.tv_sec = time.tv_usec = 0;
    FD_ZERO(&read_fds);
    FD_SET(fd, &read_fds);
    if (fd > FD_SETSIZE) {
	fprintf(stderr, "Internal Error, stupid WINSOCK large FD bug.\n");
	fprintf(stderr, "Increase FD_SETSIZE.  Item not added to fdset.\n");
    }
    ret_val = select(FD_SETSIZE, &read_fds, NULL, NULL, &time);
    return (ret_val > 0);
}

IOinterface_func ffs_file_read_func = (IOinterface_func)nt_file_read_func;
IOinterface_func ffs_file_write_func = (IOinterface_func)nt_file_write_func;
IOinterface_funcv ffs_file_readv_func = (IOinterface_funcv)null_file_readv_func;
IOinterface_funcv ffs_file_writev_func = (IOinterface_funcv)null_file_writev_func;
IOinterface_lseek ffs_file_lseek_func = (IOinterface_lseek)nt_file_lseek_func;


IOinterface_func ffs_read_func = (IOinterface_func)nt_socket_read_func;
IOinterface_func ffs_write_func = (IOinterface_func)nt_socket_write_func;
IOinterface_funcv ffs_readv_func = (IOinterface_funcv)nt_socket_readv_func;
IOinterface_funcv ffs_writev_func = NULL;
int ffs_max_iov = 1;


IOinterface_open ffs_file_open_func = (IOinterface_open)nt_file_open_func;
IOinterface_close ffs_close_func = (IOinterface_close) nt_close_func;
IOinterface_poll  ffs_poll_func = (IOinterface_poll)nt_poll_func;
IOinterface_func ffs_server_read_func = (IOinterface_func)nt_socket_read_func;
IOinterface_func ffs_server_write_func = (IOinterface_func)nt_socket_write_func;
IOinterface_init ffs_sockets_init_func = (IOinterface_init)nt_socket_init_func;
