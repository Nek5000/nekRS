#ifndef FFS_INTERNAL_H
#include "../fm/fm.h"
#include "../fm/fm_internal.h"

#define MAGIC_NUMBER 0x4356ffa9	/* random magic */
#define REVERSE_MAGIC_NUMBER 0xa9ff5643		/* byte reversed random
						 * magic */
#define MAGIC_FLOAT 0.0078125
	/* random float */
#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define DEFAULT_FS_PORT 5347

extern char *architecture;

#ifdef __LIBCATAMOUNT__

#define htons(X) bswap_16 (X)
#define htonl(X) bswap_32 (X)
#define ntohs(X) bswap_16 (X)
#define ntohl(X) bswap_32 (X)
extern short bswap_16(short s);
extern int bswap_32(int l);
#endif

#if defined(_MSC_VER) && !defined(strdup)
#define strdup _strdup
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#define strncpy(to, from, len) strcpy_s(to, len, from)
#endif

struct _FFSBuffer {
    void *tmp_buffer;
    ssize_t tmp_buffer_size;
    ssize_t tmp_buffer_in_use_size;
};

typedef struct _internal_iovec {
    void *iov_base;
    size_t iov_offset;
    size_t iov_len;
} internal_iovec;

struct _IOContextStruct {
    FMContext fmc;

    struct _FFSBuffer tmp;
    int handle_list_size;
    FFSTypeHandle *handle_list;
};

struct _IOConversionStruct;

typedef struct _IOConversionStruct *IOConversionPtr;

typedef enum {
    not_checked, conversion_set, none_available
} Conversion_status;

struct _FFSTypeHandle {
    FFSContext context;
    int format_id;
    FFSTypeHandle *subformats;
    FFSTypeHandle *field_subformats;
    IOConversionPtr conversion;
    Conversion_status status;
    FFSTypeHandle conversion_target;
    FMFormat body;
    int is_fixed_target;
};

struct _IOgetFieldStruct {
    int offset;
    size_t size;
    FMdata_type data_type;
    unsigned char byte_swap;
    unsigned char src_float_format;
    unsigned char target_float_format;
};

typedef enum {
    none_required, direct_to_mem, buffer_and_convert, copy_dynamic_portion
} IOconversion_type;

typedef enum {
    no_row_column_swap, swap_source_column_major, swap_source_row_major
} row_column_swap_type;

typedef struct _IOconvFieldStruct {
    struct _FMgetFieldStruct src_field;
    FMVarInfoStruct *iovar;
    ssize_t dest_offset;
    int dest_size;
    void *default_value;
    row_column_swap_type rc_swap;
    IOConversionPtr subconversion;
} IOconvFieldStruct;

typedef void (*conv_routine) (void *src, void *dest, void *cs);

typedef struct _IOConversionStruct {
    IOconversion_type conversion_type;
    int notify_of_format_change;
    FMFieldList native_field_list;
    int conv_count;
    size_t base_size_delta;	/* native size - file record length */
    double max_var_expansion;
    int target_pointer_size;
    FFSContext context;
    FFSTypeHandle ioformat;
    void *free_data;
    void (*free_func)(void*);
    char *conv_pkg;
    conv_routine conv_func;
    conv_routine conv_func4;
    conv_routine conv_func2;
    conv_routine conv_func1;
    int required_alignment;
    size_t string_offset_size;
    int converted_strings;
    IOconvFieldStruct conversions[1];
} IOConversionStruct;

typedef struct _FMFormat_Comp_result{
	int diff1;
	int diff2;
} FMFormat_Comp_result;

#if SIZEOF_INT == 4
#define FILE_INT int
#else
Whoops...  Unhandled.
#endif

extern void *io_malloc(int size);
extern void io_free(void *ptr);
extern void *io_realloc(void *ptr, int size);
extern char *io_strdup(const char *str);
extern FMdata_type str_to_data_type(const char *str);
extern FMdata_type array_str_to_data_type(const char *str, 
						long *element_count_ptr);
extern const char *data_type_to_str(FMdata_type data_type);
extern int field_offset_compar(const void *a, const void *b);
extern void
FFSconvert_record(IOConversionPtr conv, void *src, void *dest, 
		  void *final_string_base, void *src_string_base);

extern void
str_to_val(char *val_str, FMdata_type dtype, int size, void **default_value);

#ifdef NOT_DEF
extern int AtomicRead(void *fd, void *buffer, int length, IOFile iofile);
extern int AtomicWrite(void *fd, void *buffer, int length, IOFile iofile);
extern int get_AtomicInt(IOFile iofile, FILE_INT * buffer);
extern int put_AtomicInt(IOFile iofile, FILE_INT * buffer);
extern void reset_read_ahead(IOFile iofile);
extern void expand_IOFile(IOFile iofile);
extern int write_IOformat(IOFile iofile, IOFormat ioformat);
#endif

extern void IOconvert_record(IOConversionPtr conv, void *src,
				   void *dest, void *final_string_area,
				   void* src_string_area);
extern int sdump_value(char *str, const char*field_type, int field_size,
			     int field_offset, FFSContext iofile, void *data,
			     void *string_base, int byte_reversal, 
			     int float_format, int encode);
extern int field_name_compar(const void *a, const void *b);
extern int count_FMfield(FMFieldList field_list);
extern char *get_IOstring_base(FMFieldPtr iofield, void *data, void *string_base);
extern void *get_IOaddr(FMFieldPtr iofield, void *data, void *string_base, int encode);
extern char *base_data_type(const char *data_type);
extern void dump_IOConversion(IOConversionPtr conv_ptr);
extern void dump_IOConversion_as_XML(IOConversionPtr conv_ptr);
extern int get_var_array_control(const char *str, FMFieldList fields);
extern int validate_data_type(const char *str);
extern int min_align_type(FMdata_type typ, int size);
extern void init_float_formats();
extern int ffs_getreg(void* s, int *reg_p, int type, int var_tmp);
extern int ffs_putreg(void* s, int reg, int type);
int ffs_localb(void* s, int size);
extern conv_routine generate_conversion(IOConversionPtr conv, 
					int src_alignment, int dest_alignment);
extern void ffs_internal_convert_field(FMFieldPtr src_spec, void *src,
				       FMdata_type dest_type, 
				       int dest_size, void *dest);
extern int field_is_flat(FMFormat f, FMTypeDesc *t);

extern void 
set_IOconversion_for_format(FFSContext iofile, FMFormat file_ioformat,
			    FMFieldList native_field_list,
				  int native_struct_size);

extern
char *
make_tmp_buffer(FFSBuffer buf, int64_t size);

extern size_t
FFS_decode_length_format(FFSContext context, FFSTypeHandle ioformat, 
			 size_t record_length);
extern void
FFS_determine_conversion(FFSContext c, FFSTypeHandle format);

extern int
FFSheader_size(FFSTypeHandle ioformat);

#define Max(i,j) ((i<j) ? j : i)

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
#define FFS_NO_LEAF_COPY 0x1

#endif
