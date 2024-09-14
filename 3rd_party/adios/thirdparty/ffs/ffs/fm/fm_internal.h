#ifndef FM_INTERNAL_H
#define FM_INTERNAL_H

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#define strncpy(dest, src, size) strcpy_s(dest, size, src)
#define snprintf sprintf_s
#endif

#define MAGIC_NUMBER 0x4356ffa9	/* random magic */
#define REVERSE_MAGIC_NUMBER 0xa9ff5643		/* byte reversed random
						 * magic */
#define MAGIC_FLOAT 0.0078125	/* random float */

extern FMdata_type FMstr_to_data_type(const char *str);
extern FMdata_type FMarray_str_to_data_type(const char *str, 
						long *element_count_ptr);
extern const char *data_type_to_str(FMdata_type data_type);
extern int field_offset_compar(const void *a, const void *b);
extern void dump_FMFormat(FMFormat ioformat);
extern FMContext new_FMContext();

typedef enum {
    FMType_pointer, FMType_array, FMType_string, FMType_subformat, FMType_simple
}FMTypeEnum;
	
typedef struct FMTypeDesc {
    struct FMTypeDesc *next;
    FMTypeEnum type;
    FMdata_type data_type;
    int pointer_recursive;
    int field_index;
    int static_size;
    int control_field_index;
} FMTypeDesc;

typedef struct FMDimen {
    int static_size;
    int control_field_index;
} FMDimen;

typedef struct _FMVarInfoStruct {
    int string;
    int var_array;
    int byte_vector;
    FMdata_type data_type;
    int dimen_count;
    FMDimen *dimens;
    FMTypeDesc type_desc;
} FMVarInfoStruct, *FMVarInfoList;

typedef struct _xml_output_info *xml_output_info;

struct _format_wire_format_1;
typedef struct _format_wire_format_1 *format_rep;

typedef struct _server_ID_struct {
    int length;
    char *value;
} server_ID_type;

typedef struct _FMFormatBody {
    int ref_count;
    FMContext context;
    char *format_name;
    int format_index;
    server_ID_type server_ID;
    int record_length;
    int byte_reversal;
    FMfloat_format float_format;
    int pointer_size;
    int IOversion;
    int field_count;
    int variant;
    int recursive;
    int alignment;
    int column_major_arrays;
    FMStructDescList master_struct_list;
    FMStructDescList orig_struct_list;	/* we don't own this memory */
    FMFormat superformat;
    FMFormat *subformats;
    FMFieldList field_list;
    FMVarInfoList var_list;
    FMFormat *field_subformats;
    FMOptInfo *opt_info;
    xml_output_info xml_out;
    format_rep server_format_rep;    
    void *ffs_info;
    void (*free_ffs_info)(void*);
} FMFormatBody;

extern long
FMget_array_element_count(FMFormat f, FMVarInfoList var, char *data, 
			  int encode);

extern FMTypeDesc*
gen_FMTypeDesc(FMFieldList fl, int field, const char *typ);

typedef struct _FMgetFieldStruct {
    size_t offset;
    int size;
    FMdata_type data_type;
    unsigned char byte_swap;
    unsigned char src_float_format;
    unsigned char target_float_format;
} FMgetFieldStruct;

extern FMfloat_format ffs_reverse_float_formats[];
extern char *base_data_type(const char *str);

struct _format_server;

typedef struct _FMContextStruct {
    int ref_count;
    int reg_format_count;
    int byte_reversal;
    int native_pointer_size;
    FMfloat_format native_float_format;
    int native_column_major_arrays;  /* False for C, True for Fortran */
    int errno_val;
    char *result;

    FMContext master_context;
    int self_server;
    int self_server_fallback;
    void *server_client_data;
    void *server_fd;
    int server_pid;
    int format_server_identifier;
    int server_byte_reversal;

    int format_list_size;
    FMFormat *format_list;
} FMContextStruct;

typedef struct _format_server *format_server;

#if SIZEOF_INT == 4
#define INT4 int
#define UINT4 unsigned int
#endif

#define INT2 short
#define UINT2 unsigned short

struct _field_wire_format_0 {  /* 16 bytes total */
    UINT2 field_name_offset;
    UINT2 field_type_offset;
    INT4 field_size;
    INT4 field_offset;
};

struct _field_wire_format_1 {  /* 16 bytes total */
    UINT4 field_name_offset;
    UINT4 field_type_offset;
    INT4 field_size;
    INT4 field_offset;
};

struct _opt_info_wire_format {  /* 12 bytes total */
    INT4 info_type;
    INT4 info_len;
    INT4 info_offset;
};

struct _format_wire_format_0 {
/*byte 0*/    UINT2 format_rep_length;  /* transmitted in net byte order */
/*byte 2*/    unsigned char record_byte_order;
/*byte 3*/    unsigned char server_rep_version;
/*byte 4*/    unsigned char subformat_count;
/*byte 5*/    unsigned char recursive_flag;
/*byte 6*/    unsigned char unused1_in_format_0;
/*byte 7*/    unsigned char unused2_in_format_0;
};

struct _format_wire_format_1 {
/*byte 0*/    UINT2 format_rep_length;  /* transmitted in net byte order */
/*byte 2*/    unsigned char record_byte_order;
/*byte 3*/    unsigned char server_rep_version;
/*byte 4*/    unsigned char subformat_count;
/*byte 5*/    unsigned char recursive_flag;
/*byte 6*/    UINT2 top_bytes_format_rep_length;  /* transmitted in net byte order */
};

struct _subformat_wire_format_0 {  /* 20 bytes for base */
/*byte 0*/    UINT2 subformat_rep_length;  /* transmitted in net byte order */
/*byte 2*/    unsigned char server_rep_version;
/*byte 3*/    unsigned char record_byte_order;
/*byte 4*/    unsigned char pointer_size;
/*byte 5*/    unsigned char header_size;
/*byte 6*/    UINT2 name_offset;	/* native host byte order */
/*byte 8*/    UINT2 field_count;
/*byte 10*/   UINT2 floating_point_rep;
/*byte 12*/   INT4 record_length;
/*byte 16*/   UINT2 opt_info_offset;
/*byte 18*/   unsigned char column_major_arrays;  /* false for C, true for Fortran */
/*byte 19*/   unsigned char alignment;
};

struct _subformat_wire_format_1 {  /* 28 bytes for base */
/*byte 0*/    UINT2 subformat_rep_length;  /* transmitted in net byte order */
/*byte 2*/    unsigned char server_rep_version;
/*byte 3*/    unsigned char record_byte_order;
/*byte 4*/    UINT4 name_offset;	/* native host byte order */
/*byte 8*/    UINT4 field_count;
/*byte 12*/   INT4 record_length;
/*byte 16*/    unsigned char pointer_size;
/*byte 17*/    unsigned char header_size;
/*byte 18*/   UINT2 floating_point_rep;
/*byte 20*/   UINT2 opt_info_offset;
/*byte 22*/   unsigned char column_major_arrays;  /* false for C, true for Fortran */
/*byte 23*/   unsigned char alignment;
/*byte 24*/   UINT2 top_bytes_subformat_rep_length;  /* transmitted in net byte order */
/*byte 26*/   UINT2 top_bytes_opt_info_offset;
};

struct _format_wire_format {
    union {
	struct _format_wire_format_0 f0;
	struct _format_wire_format_1 f1;
    }f;
};

struct _subformat_wire_format {
    union {
	struct _subformat_wire_format_0 f0;
	struct _subformat_wire_format_1 f1;
    }f;
};

typedef struct {
    unsigned char version;
    unsigned char salt;
    unsigned INT2 port;
    unsigned INT4 IP_addr;
    unsigned INT2 format_identifier;
} version_1_format_ID;

typedef struct {
    unsigned char version;
    unsigned char unused;
    unsigned INT2 rep_len;
    unsigned INT4 hash1;
    unsigned INT4 hash2;
} version_2_format_ID;

typedef struct {
    unsigned char version;
    unsigned char top_byte_rep_len;
    unsigned INT2 rep_len;
    unsigned INT4 hash1;
    unsigned INT4 hash2;
} version_3_format_ID;

#define DEFAULT_FS_PORT 5347

#if SIZEOF_LONG_DOUBLE != 0 && SIZEOF_LONG_DOUBLE != SIZEOF_DOUBLE
#define MAX_FLOAT_TYPE long double
#else
#define MAX_FLOAT_TYPE double
#define MAX_FLOAT_GET get_FMdouble
#endif
#if SIZEOF_LONG_LONG != 0
#define MAX_INTEGER_TYPE long long
#else
#define MAX_INTEGER_TYPE long
#endif

#define MAX_UNSIGNED_TYPE unsigned MAX_INTEGER_TYPE

typedef int (*IOinterface_func)(void *conn, void *buffer, size_t length,
				      int *errno_p, char **result_p);

#if !defined(HAVE_IOVEC_DEFINE) && !defined(_STRUCT_IOVEC) && !(defined(_BITS_UIO_H))
#define _STRUCT_IOVEC
struct	iovec {
    const void *iov_base;
    size_t iov_len;
};
#else
#ifdef _MSC_VER
#include "winsock.h"
#else
#include <sys/socket.h>
#endif
#endif

typedef int (*IOinterface_funcv)(void *conn, struct iovec *iov, 
				 int icount, int *errno_p, 
				 char **result_p);

typedef int (*IOinterface_close)(void *conn);

typedef int (*IOinterface_poll)(void *conn);

typedef void *(*IOinterface_open)(const char *path,
					const char *flag_str, 
					int *input, int *output);
typedef int (*IOinterface_lseek)(void* conn, size_t pos, int cmd);
typedef void (*IOinterface_init)(void );

extern IOinterface_func ffs_server_read_func;
extern IOinterface_func ffs_server_write_func;
extern IOinterface_init ffs_sockets_init_func;

extern int version_of_format_ID(void *server_ID);
extern int FFS_gen_authentication (unsigned char *outbuf);
extern int serverAtomicRead(void *fd, void *buffer, int length);
extern void stringify_server_ID(unsigned char *ID, char *buffer, int len);
extern void 
generate_format3_server_ID(server_ID_type *server_ID,
			   struct _format_wire_format_1 *server_format_rep);
extern void
add_format_to_iofile(FMContext fmc, FMFormat ioformat, int id_size, 
		     void *id_buffer, int index);
typedef enum {
    local_only, never_local, host_only, host_and_fallback
} action_t;

extern int establish_server_connection(FMContext iofile, action_t action);
extern void general_format_server(int port, int do_restart, int verbose, int do_proxy);
extern void dump_FMFormat(FMFormat ioformat);
extern int format_server_restarted(FMContext context);
extern int FMhas_XML_info(FMFormat format);
extern int get_internal_format_server_identifier(format_server fs);

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#define malloc(x) ffs_malloc(x)
#define realloc(ptr, s) ffs_realloc(ptr, s)
void* ffs_malloc(size_t s);
void* ffs_realloc(void* ptr, size_t s);
#endif
