#include "config.h"
#include <stdlib.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <string.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#ifdef HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#include <stdio.h>
#include "fm.h"
#include "fm_internal.h"
#include "assert.h"


static int
null_read_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    return 0;
    
}

static int
null_readv_func(conn, iov, icount, errno_p, result_p)
void *conn;
struct iovec *iov;
int icount;
int *errno_p;
char **result_p;
{
    return 0;
}

static int
null_write_func(conn, buffer, length, errno_p, result_p)
void *conn;
void *buffer;
int length;
int *errno_p;
char **result_p;
{
    return 0;
}

static int
null_writev_func(conn, iov, iovcnt, errno_p, result_p)
void *conn;
struct iovec *iov;
int iovcnt;
int *errno_p;
char **result_p;
{
    return 0;
}

static int
null_close_func(conn)
void *conn;
{
    return 0;
    
}

static void *
null_file_open_func(path, flag_str, input, output)
const char *path;
const char *flag_str;
int *input;
int *output;
{
    return NULL;
    
}
    
static int
null_poll_func(conn)
void *conn;
{
    return 0;
    
}

IOinterface_func ffs_file_read_func = (IOinterface_func)null_read_func;
IOinterface_func ffs_file_write_func = (IOinterface_func)null_write_func;
/*IOinterface_funcv ffs_file_readv_func = unix_readv_func;
  IOinterface_funcv ffs_file_writev_func = unix_writev_func;*/

IOinterface_func ffs_read_func = (IOinterface_func) null_read_func;
IOinterface_func ffs_write_func = (IOinterface_func)null_write_func;
/*IOinterface_funcv ffs_readv_func = unix_readv_func;
  IOinterface_funcv ffs_writev_func = unix_writev_func;*/
#ifndef IOV_MAX
#define IOV_MAX 16
#endif
int ffs_max_iov = IOV_MAX;
IOinterface_close ffs_close_func = (IOinterface_close)null_close_func;
IOinterface_poll  ffs_poll_func = (IOinterface_poll)null_poll_func;
IOinterface_open ffs_file_open_func = (IOinterface_open)null_file_open_func;
IOinterface_func ffs_server_read_func = (IOinterface_func)null_read_func;
IOinterface_func ffs_server_write_func = (IOinterface_func)null_write_func;
IOinterface_init ffs_sockets_init_func = NULL;
