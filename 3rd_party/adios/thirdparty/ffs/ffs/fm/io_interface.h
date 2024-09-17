#if defined(_MSC_VER) && !defined(_STRUCT_IOVEC)
#define _STRUCT_IOVEC
struct	iovec {
    const void *iov_base;
    size_t	iov_len;
};
#endif

#ifndef _MSC_VER
#define SOCKET int
#endif

#ifndef FM_INTERNAL_H
typedef int (*IOinterface_func)(void *conn, void *buffer, size_t length,
				      int *errno_p, char **result_p);

typedef int (*IOinterface_funcv)(void *conn, struct iovec *iov, 
				       int icount, int *errno_p, 
				       char **result_p);

typedef int (*IOinterface_close)(void *conn);

typedef int (*IOinterface_poll)(void *conn);

typedef void *(*IOinterface_open)(const char *path, const char *flag_str, int *input, int *output);
typedef void (*IOinterface_init)(void );
typedef int (*IOinterface_lseek)(void* conn, size_t pos, int cmd);
#endif

extern IOinterface_func ffs_file_read_func;
extern IOinterface_func ffs_file_write_func;
extern IOinterface_funcv ffs_file_readv_func;
extern IOinterface_funcv ffs_file_writev_func;
extern IOinterface_lseek ffs_file_lseek_func;

extern IOinterface_func ffs_read_func;
extern IOinterface_func ffs_write_func;
extern IOinterface_funcv ffs_readv_func;
extern IOinterface_funcv ffs_writev_func;
extern int ffs_max_iov;
extern IOinterface_close ffs_close_func;
extern IOinterface_open ffs_file_open_func;
extern IOinterface_init ffs_sockets_init_func;

#ifdef __FFS_H__
extern void
set_interface_FFSFile(FFSFile f, IOinterface_func write_func, 
		      IOinterface_func read_func, 
		      IOinterface_funcv writev_func,
		      IOinterface_funcv readv_func, int max_iov,
		      IOinterface_close close_func);
extern void *
get_conn_FFSFile(FFSFile f);

extern void
set_conn_FFSFile(FFSFile f, void *conn);

extern FFSFile
open_created_FFSFile(FFSFile f, char *flags);

extern void
set_socket_interface_FFSFile(FFSFile f);

extern void
set_file_interface_FFSFile(FFSFile f);

#endif
