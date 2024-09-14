
#include "config.h"

#ifndef MODULE
#include <stdio.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <ctype.h>
#include <sys/types.h>
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#define HAVE_IOVEC_DEFINE
#endif
#else
/* kernel build */
#include "kernel/pbio_kernel.h"
#include "kernel/kpbio.h"
#include "kernel/library.h"
#endif

#include "fm.h"
#include "fm_internal.h"
#include "assert.h"

#if !defined(tolower) & !defined(_MSC_VER)
extern int tolower(int);
#endif

extern FMfloat_format fm_my_float_format;
static size_t get_offset(void *, int, int);
static MAX_UNSIGNED_TYPE get_big_unsigned(FMFieldPtr field, void *data);
static MAX_FLOAT_TYPE get_big_float(FMFieldPtr field, void *data);

static int get_long_warn = 0;

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

static MAX_INTEGER_TYPE
get_big_int(FMFieldPtr field, void *data)
{
    if (field->data_type == integer_type) {
	if (field->size == sizeof(char)) {
	    char tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(char));
	    return (long) tmp;
	} else if (field->size == sizeof(short)) {
	    short tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(short));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int) sizeof(short));
	    return (long) tmp;
	} else if (field->size == sizeof(int)) {
	    int tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(int));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int) sizeof(int));
	    return (long) tmp;
	} else if (field->size == sizeof(size_t)) {
	    size_t tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(size_t));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(size_t));
	    return tmp;
	} else if (field->size == 2 * sizeof(long)) {
	    long tmp;
	    size_t low_bytes_offset = field->offset;
	    if (WORDS_BIGENDIAN) {
		if (!field->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (field->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    }
	    memcpy(&tmp, (char *) data + low_bytes_offset, sizeof(long));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(long));
	    return tmp;
	} else {
	    if (!get_long_warn) {
		fprintf(stderr, "Get Long failed!  Size problems.  File int size is %d.\n", field->size);
		get_long_warn++;
	    }
	    return -1;
	}
    } else if (field->data_type == unsigned_type) {
	MAX_UNSIGNED_TYPE tmp = get_big_unsigned(field, data);
	return (MAX_UNSIGNED_TYPE) tmp;
    } else if (field->data_type == float_type) {
	MAX_FLOAT_TYPE tmp = get_big_float(field, data);
#ifndef METICULOUS_FLOATS_AND_LONGS
	return (MAX_INTEGER_TYPE) (long) (double) tmp;
#else
	return (MAX_INTEGER_TYPE) tmp;
#endif
    } else {
	fprintf(stderr, "Get IOlong failed on invalid data type!\n");
	exit(1);
    }
    /* NOTREACHED */
    return 0;
}

static MAX_UNSIGNED_TYPE
get_big_unsigned(FMFieldPtr field, void *data)
{
    if ((field->data_type == unsigned_type) || 
	(field->data_type == enumeration_type) || 
	(field->data_type == boolean_type)) {
	if (field->size == sizeof(char)) {
	    unsigned char tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(char));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (field->size == sizeof(short)) {
	    unsigned short tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(short));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(short));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (field->size == sizeof(int)) {
	    unsigned int tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(int));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(int));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (field->size == sizeof(size_t)) {
	    size_t tmp;
	    memcpy(&tmp, (char *) data + field->offset, sizeof(size_t));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(size_t));
	    return tmp;
	} else if (field->size == 2 * sizeof(long)) {
	    unsigned long tmp;
	    size_t low_bytes_offset = field->offset;
	    if (WORDS_BIGENDIAN) {
		if (!field->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (field->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    }
	    memcpy(&tmp, (char *) data + low_bytes_offset, sizeof(long));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(long));
	    return tmp;
	} else {
	    if (!get_long_warn) {
		fprintf(stderr, "Get Long failed!  Size problems.  File int size is %d.\n", field->size);
		get_long_warn++;
	    }
	    return 0;
	}
    } else if (field->data_type == integer_type) {
	MAX_INTEGER_TYPE tmp = get_big_int(field, data);
	return (MAX_UNSIGNED_TYPE) tmp;
    } else if (field->data_type == float_type) {
	MAX_FLOAT_TYPE tmp = get_big_float(field, data);
#ifndef METICULOUS_FLOATS_AND_LONGS
	return (MAX_UNSIGNED_TYPE) (long) (double) tmp;
#else
	return (MAX_UNSIGNED_TYPE) tmp;
#endif
    } else {
	fprintf(stderr, "Get IOulong failed on invalid data type!\n");
	exit(1);
    }
    /* NOTREACHED */
    return 0;
}

extern FMfloat_format fm_reverse_float_formats[];
#define CONV(a,b) ((a<<16) + b)
static void
float_conversion(unsigned char*value, int size, FMfloat_format src_format,
		 FMfloat_format dest_format)
{
    int tmp;
    if (src_format == dest_format) return;
    if (src_format == fm_reverse_float_formats[dest_format]) {
	byte_swap((char *) value, size);
	return;
    }
    switch (CONV(src_format, dest_format)) {
    case CONV(Format_IEEE_754_bigendian, Format_IEEE_754_mixedendian):
    case CONV(Format_IEEE_754_mixedendian, Format_IEEE_754_bigendian):
	byte_swap((char*)&value[0], 4);
	byte_swap((char*)&value[4], 4);
	break;
    case CONV(Format_IEEE_754_littleendian, Format_IEEE_754_mixedendian):
    case CONV(Format_IEEE_754_mixedendian, Format_IEEE_754_littleendian):
	tmp = *(int*)&value[0];
	*(int*)&value[0] = *(int*)&value[4];
	*(int*)&value[4] = tmp;
	break;
    default:
	printf("unanticipated float conversion \n");
    }
}

static MAX_FLOAT_TYPE
get_big_float(FMFieldPtr field, void *data)
{
    if (field->data_type == float_type) {
	if (field->size == sizeof(float)) {
	    float tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + field->offset), sizeof(float));
	    if (field->byte_swap)
		byte_swap((char *) &tmp, (int)sizeof(float));
	    tmp2 = tmp;
	    return tmp2;
	} else if (field->size == sizeof(double)) {
	    double tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + field->offset), sizeof(double));
	    float_conversion((unsigned char *)&tmp, sizeof(double), 
			     (FMfloat_format)field->src_float_format,
			     (FMfloat_format)field->target_float_format);
	    tmp2 = tmp;
	    return tmp2;
#if SIZEOF_LONG_DOUBLE != 0
	} else if (field->size == sizeof(long double)) {
	    long double tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + field->offset),
		   sizeof(long double));
	    if (field->byte_swap)
		byte_swap((char *) &tmp,
			  (int)sizeof(long double));
	    tmp2 = tmp;
	    return tmp2;
#endif
	} else {
	    static int get_double_warn = 0;
	    if (!get_double_warn) {
		fprintf(stderr, "Get Double failed!  Size problems.  File double size is %d.\n", field->size);
		get_double_warn++;
	    }
	    return 0.0;
	}
    } else if (field->data_type == integer_type) {
	MAX_INTEGER_TYPE tmp = get_big_int(field, data);
#ifndef METICULOUS_FLOATS_AND_LONGS
	/* 
	 * A concession to inter-compiler interoperability...  We
	 * shouldn't need the (double)(long) casts here.  If we don't
	 * use them, AND we're on a machine which doesn't support "long 
	 * long" and "long double" in native code, gcc generates calls 
	 * to libgcc and the resulting library can't be linked with
	 * anything but GCC.  This can be a problem for installed
	 * libraries.  Using the casts avoids that particular problem,
	 * at a cost of data loss in the case of someone converting a
	 * "long long" that doesn't fit in a long to a floating point
	 * value. 
	 */
	return (MAX_FLOAT_TYPE) (double) (long) tmp;
#else
	return (MAX_FLOAT_TYPE) tmp;
#endif
    } else if (field->data_type == unsigned_type) {
	MAX_UNSIGNED_TYPE tmp = get_big_unsigned(field, data);
#ifndef METICULOUS_FLOATS_AND_LONGS
	return (MAX_FLOAT_TYPE) (double) (long) tmp;
#else
	return (MAX_FLOAT_TYPE) tmp;
#endif
    } else {
	fprintf(stderr, "Get Double failed on invalid data type!\n");
	exit(1);
    }
    /* NOTREACHED */
    return 0;
}

extern float
get_FMfloat(FMFieldPtr field, void *data)
{
    float tmp;
    tmp = (float)get_big_float(field, data);
    return tmp;
}

extern double
get_FMdouble(FMFieldPtr field, void *data)
{
    double tmp;
    tmp = get_big_float(field, data);

    return tmp;
}

#if SIZEOF_LONG_DOUBLE != 0
extern long double
get_FMlong_double(FMFieldPtr field, void *data)
{
    long double tmp;
    tmp = get_big_float(field, data);
    return tmp;
}
#endif

extern short
get_FMshort(FMFieldPtr field, void *data)
{
    short tmp = (short)get_big_int(field, data);
    return tmp;
}

extern unsigned short
get_FMushort(FMFieldPtr field, void *data)
{
    unsigned short tmp = (unsigned short) get_big_unsigned(field, data);
    return tmp;
}

extern int
get_FMint(FMFieldPtr field, void *data)
{
    int tmp = (int) get_big_int(field, data);
    return tmp;
}

extern unsigned int
get_FMuint(FMFieldPtr field, void *data)
{
    unsigned int tmp = (unsigned int) get_big_unsigned(field, data);
    return tmp;
}


extern size_t
get_FMlong(FMFieldPtr field, void *data)
{
    ssize_t tmp = get_big_int(field, data);
    return tmp;
}

extern size_t
get_FMulong(FMFieldPtr field, void *data)
{
    size_t tmp = get_big_unsigned(field, data);
    return tmp;
}

#if SIZEOF_LONG_LONG != 0
extern long long
get_FMlong_long(FMFieldPtr field, void *data)
{
    long long tmp= get_big_int(field, data);
    return tmp;
}

extern unsigned long long
get_FMulong_long(FMFieldPtr field, void *data)
{
    unsigned long long tmp = get_big_unsigned(field, data);
    return tmp;
}

#endif

extern void
get_FMlong8(FMFieldPtr field, void *data, unsigned long *low_long, long *high_long)
{
    *low_long = 0;
    if (high_long)
	*high_long = 0;
    if (field->data_type == integer_type) {
	if (field->size == 2 * sizeof(long)) {
	    size_t low_bytes_offset = field->offset;
	    size_t high_bytes_offset = field->offset;
	    FMgetFieldStruct tmp_field;  /*OK */
	    tmp_field = *field;
	    if (WORDS_BIGENDIAN) {
		if (field->byte_swap) {
		    high_bytes_offset += sizeof(long);
		} else {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (field->byte_swap) {
		    low_bytes_offset += sizeof(long);
		} else {
		    high_bytes_offset += sizeof(long);
		}
	    }
	    tmp_field.size = sizeof(long);
	    tmp_field.offset = low_bytes_offset;
	    *low_long = (unsigned long) get_FMulong(&tmp_field, data);
	    if (high_long) {
		tmp_field = *field;
		tmp_field.size = sizeof(long);
		tmp_field.offset = high_bytes_offset;
		*high_long = (long) get_FMlong(&tmp_field, data);
	    }
	} else {
	    *low_long = (long) get_FMlong(field, data);
	}
    } else if (field->data_type == float_type) {
	MAX_FLOAT_TYPE tmp;
	tmp = get_big_float(field, data);
	*low_long = (long) tmp;
    } else {
	fprintf(stderr, "Get IOlong8 failed on invalid data type!\n");
	exit(1);
    }
}

extern int
get_FMulong8(FMFieldPtr field, void *data, unsigned long *low_long, unsigned long *high_long)
{
    *low_long = 0;
    if (high_long)
	*high_long = 0;
    if (field->data_type == unsigned_type) {
	if (field->size == 2 * sizeof(long)) {
	    size_t low_bytes_offset = field->offset;
	    size_t high_bytes_offset = field->offset;
	    FMgetFieldStruct tmp_field;  /*OK */
	    tmp_field = *field;
	    if (WORDS_BIGENDIAN) {
		if (field->byte_swap) {
		    high_bytes_offset += sizeof(long);
		} else {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (field->byte_swap) {
		    high_bytes_offset += sizeof(long);
		} else {
		    low_bytes_offset += sizeof(long);
		}
	    }
	    tmp_field.size = sizeof(unsigned long);
	    tmp_field.offset = low_bytes_offset;
	    *low_long = (unsigned long) get_FMulong(&tmp_field, data);
	    if (high_long) {
		tmp_field = *field;
		tmp_field.size = sizeof(unsigned long);
		tmp_field.offset = high_bytes_offset;
		*high_long = (unsigned long) get_FMulong(&tmp_field, data);
	    }
	    return 0;
	} else {
	    *low_long = (unsigned long) get_FMulong(field, data);
	    return 0;
	}
    } else if (field->data_type == integer_type) {
	/* this is a bit to ugly to handle at the moment.. */
	assert(0);
	return 0;
    } else if (field->data_type == float_type) {
	MAX_FLOAT_TYPE tmp;
	tmp = get_big_float(field, data);
	*low_long = (long) tmp;
	return 0;
    } else {
	fprintf(stderr, "Get IOlong8 failed on invalid data type!\n");
	exit(1);
    }
    /* NOTREACHED */
    return 0;
}

extern char
get_FMchar(FMFieldPtr field, void *data)
{
    if (field->data_type == char_type) {
	return (char) *((char *) data + field->offset);
    } else {
	fprintf(stderr, "Get Char failed on invalid data type!\n");
	exit(1);
    }
    /* NOTREACHED */
    return 0;
}

extern int
get_FMenum(FMFieldPtr field, void *data)
{
    FMgetFieldStruct tmp_field;  /*OK */
    tmp_field = *field;
    tmp_field.data_type = integer_type;
    return get_FMint(&tmp_field, data);
}

extern char *
get_FMstring(FMFieldPtr field, void *data)
{
    return get_FMstring_base(field, data, data);
}

extern char *
get_FMstring_base(FMFieldPtr field, void *data, void *string_base)
{
    uintptr_t offset = get_offset((void *) ((char *) data + field->offset),
			     field->size, field->byte_swap);
    if (offset == 0) {
	return NULL;
    } else if (offset > (uintptr_t) data) {		/* probably *
							 * converted *
							 * string */
	return (char *) offset;
    } else {
	return (char *) string_base + offset;
    }
}

extern void *
get_FMaddr(FMFieldPtr field, void *data, void *string_base, int encode)
{
    uintptr_t offset = get_offset((void *) ((char *) data + field->offset),
			     field->size, field->byte_swap);
    if (offset == 0) {
	return NULL;
    } else if (!encode) {
	return (char *) offset;
    } else {
	return (char *) string_base + offset;
    }
}

static size_t
get_offset(void *data, int size, int swap)
{
    FMgetFieldStruct field;  /*OK */
    field.offset = 0;
    field.size = size;
    field.data_type = integer_type;
    field.byte_swap = swap;
    if (size == sizeof(int)) {
	return get_FMlong(&field, data);
    } else {
	field.offset = size - sizeof(size_t);
	field.size = sizeof(size_t);
	return get_FMlong(&field, data);
    }
}

extern int
FM_field_type_eq(const char *str1, const char *str2)
{
    FMdata_type t1, t2;
    long t1_count, t2_count;

    t1 = FMarray_str_to_data_type(str1, &t1_count);
    t2 = FMarray_str_to_data_type(str2, &t2_count);

    if ((t1_count == -1) && (t2_count == -1)) {
	/* variant array */
	char *tmp_str1 = base_data_type(str1);
	char *tmp_str2 = base_data_type(str2);
	
	char *colon1 = strchr(tmp_str1, ':');
	char *colon2 = strchr(tmp_str2, ':');
	char *lparen1 = strchr(str1, '[');
	char *lparen2 = strchr(str2, '[');
	intptr_t count1 = 0;
	intptr_t count2 = 0;

	if (colon1 != NULL) {
	    count1 = colon1 - tmp_str1;
	} else {
	    count1 = strlen(tmp_str1);
	}
	if (colon2 != NULL) {
	    count2 = colon2 - tmp_str2;
	} else {
	    count2 = strlen(tmp_str2);
	}
	/*compare base type */
	if (strncmp(tmp_str1, tmp_str2,(count1>count2)?count1:count2) != 0) {
	    /* base types differ */
	    return 0;
	}
	free(tmp_str1);
	free(tmp_str2);
	if ((lparen1 == NULL) || (lparen2 == NULL)) return -1;
	return (strcmp(lparen1, lparen2) == 0);
    }
    return ((t1 == t2) && (t1_count == t2_count));
}

extern int
IOget_array_size_dimen(const char *str, FMFieldList fields, int dimen, int *control_field)
{
    char *left_paren, *end;
    long static_size;

    *control_field = -1;
    if ((left_paren = strchr(str, '[')) == NULL) {
	return 0;
    }	
    while (dimen != 0) {
	left_paren = strchr(left_paren + 1, '[');
	if (left_paren == NULL) return 0;
	dimen--;
    }
    static_size = strtol(left_paren + 1, &end, 0);
    if (left_paren + 1 == end) {
	/* dynamic element */
	char field_name[1024];
	int count = 0;
	int i = 0;
	while (((left_paren+1)[count] != ']') &&
	       ((left_paren+1)[count] != 0)) {
	    field_name[count] = (left_paren+1)[count];
	    count++;
	}
	field_name[count] = 0;
	while (fields[i].field_name != NULL) {
	    if (strcmp(field_name, fields[i].field_name) == 0) {
		if (FMstr_to_data_type(fields[i].field_type) ==
		    integer_type) {
		    *control_field = i;
		    return -1;
		} else {
		    fprintf(stderr, "Variable length control field \"%s\" not of integer type.\n", field_name);
		    return 0;
		}
	    }
	    i++;
	}
	fprintf(stderr, "Array dimension \"%s\" in type spec\"%s\" not recognized.\n",
		field_name, str);
	fprintf(stderr, "Dimension must be a field name (for dynamic arrays) or a positive integer.\n");
	fprintf(stderr, "To use a #define'd value for the dimension, use the IOArrayDecl() macro.\n");
	return -1;
    }
    if (*end != ']') {
	fprintf(stderr, "Malformed array dimension, unexpected character '%c' in type spec \"%s\"\n",
		*end, str);
	fprintf(stderr, "Dimension must be a field name (for dynamic arrays) or a positive integer.\n");
	fprintf(stderr, "To use a #define'd value for the dimension, use the IOArrayDecl() macro.\n");
	return -1;
    }
    if (static_size <= 0) {
	fprintf(stderr, "Non-positive array dimension %ld in type spec \"%s\"\n",
		static_size, str);
	fprintf(stderr, "Dimension must be a field name (for dynamic arrays) or a positive integer.\n");
	fprintf(stderr, "To use a #define'd value for the dimension, use the IOArrayDecl() macro.\n");
	return -1;
    }
    return static_size;
}

extern FMFieldPtr
get_FMfieldPtrFromList(FMFieldList field_list, const char *fieldname)
{
    int index;
    FMFieldPtr ret_val;
    FMdata_type data_type;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL)
	return NULL;

    data_type = FMstr_to_data_type(field_list[index].field_type);
    if (data_type == unknown_type) {
	fprintf(stderr, "Unknown field type for field %s\n",
		field_list[index].field_name);
	return NULL;
    }
    ret_val = (FMFieldPtr) malloc(sizeof(*ret_val));
    ret_val->src_float_format = Format_Unknown;
    ret_val->target_float_format = fm_my_float_format;
    ret_val->offset = field_list[index].field_offset;
    ret_val->size = field_list[index].field_size;
    ret_val->data_type = data_type;
    ret_val->byte_swap = 0;
    return ret_val;
}

extern void *
get_FMfieldAddr_by_name(FMFieldList field_list, const char *fieldname, void *data)
{
    int index;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL)
	return NULL;

    return ((char*)data) + field_list[index].field_offset;
}

extern int
get_FMfieldInt_by_name(FMFieldList field_list, const char *fieldname, void *data)
{
    int index;
    FMgetFieldStruct descr;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL) {
	printf("Field %s not found\n", fieldname);
	return 0;
    }

    descr.offset = field_list[index].field_offset;
    descr.size = field_list[index].field_size;
    descr.data_type = integer_type;
    descr.byte_swap = 0;
    descr.src_float_format = descr.target_float_format = Format_IEEE_754_bigendian;
    return get_FMint(&descr, data);
}

extern size_t
get_FMfieldLong_by_name(FMFieldList field_list, const char *fieldname, void *data)
{
    int index;
    FMgetFieldStruct descr;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL) {
	printf("Field %s not found\n", fieldname);
	return 0;
    }

    descr.offset = field_list[index].field_offset;
    descr.size = field_list[index].field_size;
    descr.data_type = integer_type;
    descr.byte_swap = 0;
    descr.src_float_format = descr.target_float_format = Format_IEEE_754_bigendian;
    return get_FMlong(&descr, data);
}

extern void *
get_FMPtrField_by_name(FMFieldList field_list, const char *fieldname, void *data, int encode)
{
    int index;
    struct _FMgetFieldStruct field;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL)
	return NULL;

    field.offset = field_list[index].field_offset;
    field.size = sizeof(void*);
    field.byte_swap = 0;

    size_t offset = get_offset((void *) ((char *) data + field.offset),
			     field.size, field.byte_swap);
    if (offset == 0) {
	return NULL;
    } else if (!encode) {
	return (char *) offset;
    } else {
	return (char *) data + offset;
    }
}

extern int
set_FMPtrField_by_name(FMFieldList field_list, const char *fieldname, void *data, void *ptr_value)
{
    int index;

    for (index = 0; field_list[index].field_name != NULL; index++) {
	if (strcmp(field_list[index].field_name, fieldname) == 0) {
	    break;
	}
    }
    if (field_list[index].field_name == NULL)
	return 0;

    
    *((void**)(((char*)data) + field_list[index].field_offset)) = ptr_value;
    return 1;
}

extern FMFieldPtr
get_FMfieldPtr(FMFormat format, const char *fieldname)
{
    int index;
    FMFieldPtr ret_val;
    FMdata_type data_type;
    long junk;

    if (format == NULL)
	return NULL;

    for (index = 0; index < format->field_count; index++) {
	if (strcmp(format->field_list[index].field_name,
		   fieldname) == 0) {
	    break;
	}
    }
    if (index >= format->field_count)
	return NULL;

    data_type = FMarray_str_to_data_type(format->field_list[index].field_type,
				       &junk);
    if (data_type == unknown_type) {
	fprintf(stderr, "Unknown field type for field %s\n",
		format->field_list[index].field_name);
	return NULL;
    }
    ret_val = (FMFieldPtr) malloc(sizeof(*ret_val));
    ret_val->offset = format->field_list[index].field_offset;
    ret_val->size = format->field_list[index].field_size;
    ret_val->data_type = data_type;
    ret_val->byte_swap = format->byte_reversal;
    ret_val->src_float_format = format->float_format;
    ret_val->target_float_format = fm_my_float_format;
    return ret_val;
}

extern FMFieldPtr
get_local_FieldPtr(FMFormat format, const char *fieldname)
{
    FMFieldPtr ret_val = get_FMfieldPtr(format, fieldname);

    /* 
     * get_local_FieldPtr() differs from get_FMfieldPtr() only in
     * that the byte_swap value is always false.  This is because
     * local Field ptrs are to be used with the records returned
     * from local conversions.  A local conversion already has
     * performed the byte_swapping. 
     */
    if (ret_val)
	ret_val->byte_swap = 0;
    return ret_val;
}
