
#include "config.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_STRING_H
#include <string.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <limits.h>
#include <string.h>
#include "assert.h"
#include "fm.h"
#include "fm_internal.h"

#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif

int FMdumpVerbose = 0;
extern FMfloat_format fm_my_float_format;

typedef struct addr_list {
    void *addr;
    size_t offset;
} addr_list_entry;

typedef struct dump_state {
    int encoded;
    size_t output_len;
    ssize_t output_limit;
    int use_XML;
    int indent;
    char *offset_base;
    int use_file_out;
    char *output_string;
    int realloc_string;
    FILE *out;
    /* 
     * entries below are to prevent multiple visits in recursively-defined types
     */
    int malloc_addr_size;
    int addr_list_is_stack;
    int addr_list_cnt;
    addr_list_entry *addr_list;
}*dstate;

static void free_addr_list(dstate s);
static ssize_t search_addr_list(dstate s, void *addr);
static void add_to_addr_list(dstate s, void *addr, size_t offset);
static int dump_subfields(void *base, FMFormat f, dstate s, size_t data_offset);

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

extern int
dump_output(dstate s, size_t length_estimate, char *format, ...)
{
    char buf[1024];
    char *tmp = &buf[0];
    va_list ap;
    int free_tmp = 0;
    int use_buf = 1;
    if ((s->output_limit != -1) && 
	(s->output_len + length_estimate > (size_t)s->output_limit)) {
	return 0;
    }
    if (s->output_string != NULL) {
	if (s->realloc_string) {
	    s->output_string = realloc(s->output_string, s->output_len + length_estimate + 1);
	}
	use_buf = 0;
	tmp = s->output_string + s->output_len;
    }
    if (use_buf && (length_estimate > 1024)) {
	tmp = malloc(length_estimate + 1);
	free_tmp++;
    }
#ifdef STDC_HEADERS
    va_start(ap, format);
#else
    va_start(ap);
#endif
    vsnprintf(tmp, sizeof(buf), format, ap);
    va_end(ap);
    s->output_len += strlen(tmp);
    if (s->use_file_out) {
	fprintf(s->out, "%s", tmp);
    }
    if (free_tmp) free(tmp);
    return 1;
}

#if SIZEOF_LONG != 8
#ifndef WORDS_BIGENDIAN
static int words_bigendian = -1;

static int
set_bigendian () {
  /* Are we little or big endian?  From Harbison&Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  words_bigendian = (u.c[sizeof (long) - 1] == 1);
  return words_bigendian;
}
#define WORDS_BIGENDIAN ((words_bigendian == -1) ? set_bigendian() : words_bigendian)
#endif
#endif

static unsigned long
quick_get_ulong(FMFieldPtr iofield, void *data)
{
    data = (void *) ((char *) data + iofield->offset);
    /* only used when field type is an integer and aligned by its size */
    switch (iofield->size) {
    case 1:
	return (unsigned long) (*((unsigned char *) data));
    case 2:
	return (unsigned long) (*((unsigned short *) data));
    case 4:
	return (unsigned long) (*((unsigned int *) data));
    case 8:
#if SIZEOF_LONG == 8
	if ((((long) data) & 0x0f) == 0) {
	    /* properly aligned */
	    return (unsigned long) (*((unsigned long *) data));
	} else {
	    union {
		unsigned long tmp;
		int tmpi[2];
	    } u;
	    u.tmpi[0] = ((int *) data)[0];
	    u.tmpi[1] = ((int *) data)[1];
	    return u.tmp;
	}
#else
	/* must be fetching 4 bytes of the 8 available */
	if (WORDS_BIGENDIAN)
	    return (*(((unsigned long *) data) +1));
	else 
	    return (*((unsigned long *) data));
#endif
    }
    return 0;
}

static void *
quick_get_pointer(FMFieldPtr field, void *data)
{
    union {
	void *p;
	unsigned long tmp;
	int tmpi[2];
    } u;
    unsigned short s;
    unsigned int i;
    data = (void *) ((char *) data + field->offset);
    /* only used when field type is an integer and aligned by its size */
    switch (field->size) {
    case 1:
	u.tmp = (unsigned long) (*((unsigned char *) data));
	break;
    case 2:
	s = (*((unsigned short *) data));
	if (field->byte_swap) byte_swap((char*)&s, 2);
	u.tmp = (unsigned long) s;
	break;
    case 4:
	i = (*((unsigned int *) data));
	if (field->byte_swap)  byte_swap((char *) &i, 4);
	u.tmp = (unsigned long) i;
	break;
    case 8:
#if SIZEOF_LONG == 8
	if ((((long) data) & 0x0f) == 0) {
	    /* properly aligned */
	    u.tmp = (unsigned long) (*((unsigned long *) data));
	} else {
	    u.tmpi[0] = ((int *) data)[0];
	    u.tmpi[1] = ((int *) data)[1];
	    if (field->byte_swap)  byte_swap((char *) &u.p, 8);
	    return u.p;
	}
#else
	/* must be fetching 4 bytes of the 8 available */
	if (WORDS_BIGENDIAN)
	    u.tmp = (*(((unsigned long *) data) +1));
	else 
	    u.tmp = (*((unsigned long *) data));

#endif
	break;
    }
    return u.p;
}

void
init_dump_state(dstate state)
{
    state->output_len = 0;
    state->output_limit = -1;
    state->use_XML = 0;
    state->addr_list_cnt = 0;
    state->indent = 0;
    state->use_file_out = 1;
    state->offset_base = NULL;
    state->output_string = NULL;
    state->realloc_string = 0;
    state->out = stdout;
}

static int
internal_dump_data(FMFormat format, void *data, dstate state);

extern int
FMdump_data(FMFormat format, void *data, int character_limit)
{
    return FMfdump_data(stdout, format, data, character_limit);
}

extern int
FMfdump_data(void *out, FMFormat format, void *data, int character_limit)
{
    int ret;
    struct dump_state state;
    init_dump_state(&state);
    state.encoded = 0;
    state.output_limit = character_limit;
    state.out = (FILE*)out;
    ret = internal_dump_data(format, data, &state);
    return ret;
}

#define STACK_ARRAY_SIZE 100

static int
internal_dump_data(FMFormat format, void *data, dstate state)
{
    int ret;
    addr_list_entry stack_addr_list[STACK_ARRAY_SIZE];
    state->addr_list_is_stack = 1;
    state->addr_list_cnt = 0;
    state->addr_list = stack_addr_list;

    if (format->recursive) {
	state->addr_list[state->addr_list_cnt].addr = data;
	state->addr_list[state->addr_list_cnt].offset = 0;
	state->addr_list_cnt++;
    }

    ret = dump_subfields(data, format, state, 0);

    free_addr_list(state);
    return ret;
}

static int
dump_subfield(void *base, FMFormat f, dstate s, size_t data_offset, void* parent_base, FMTypeDesc *t);

#define FALSE 0
#define TRUE 1

#ifdef NOTDEF
static int
field_is_flat(FMFormat f, FMTypeDesc *t)
{
    switch (t->type) {
    case FMType_pointer:
	return FALSE;
    case FMType_array:
	return field_is_flat(f, t->next);
    case FMType_string:
	return FALSE;
    case FMType_subformat:
	return !f->field_subformats[t->field_index]->variant;
    case FMType_simple:
	return TRUE;
    default:
	assert(FALSE);
    }
    /* notreached */
    return FALSE;
}
#endif

static void
do_indent(dstate s, int indent)
{
    int i;
    for (i=0; i < indent; i++) {
	dump_output(s, 2, "  ");
    }
}

static void
start_field(dstate s, FMFieldList f, FMTypeDesc *t)
{
    if ((s->indent != -1) && (t->type != FMType_simple) && (t->type != FMType_pointer) && (t->type != FMType_string)) {
	do_indent(s, s->indent++);
    }
    if (s->use_XML) {
	dump_output(s, strlen(f->field_name) + 2, "<%s>", f->field_name);
    } else {
	dump_output(s, strlen(f->field_name) + 3, "%s = ", f->field_name);
    }
    if ((t->type != FMType_simple) && (t->type != FMType_pointer) && (t->type != FMType_string)) 
	dump_output(s, 1, "\n");
}

static void
stop_field(dstate s, FMFieldList f, FMTypeDesc *t)
{
    if ((s->indent != -1) && (t->type != FMType_simple) && (t->type != FMType_pointer) && (t->type != FMType_string)) {
	do_indent(s, --(s->indent));
    }
    if (s->use_XML) {
	dump_output(s, strlen(f->field_name) + 3, "</%s>", f->field_name);
    } else {
	dump_output(s, 1, ",");
    }
    if ((t->type != FMType_simple) && (t->type != FMType_pointer) && (t->type != FMType_string)) 
	dump_output(s, 1, "\n");
}

static int
dump_subfields(void *base, FMFormat f, dstate s, size_t data_offset)
{
    int i;

    for (i = 0; i < f->field_count; i++) {
	size_t subfield_offset = data_offset + f->field_list[i].field_offset;
	FMFieldList fmfield = &f->field_list[i];
	int ret;
	start_field(s, fmfield, &f->var_list[i].type_desc);

	ret = dump_subfield(base, f, s, subfield_offset, (char*)base + data_offset, 
			      &f->var_list[i].type_desc);
	stop_field(s, fmfield, &f->var_list[i].type_desc);
	if (ret != 1) return 0;
	if ((s->output_limit != -1) && (s->output_len >= (size_t)s->output_limit)) return 0;
    }
    return 1;
}

static int
determine_dump_size(FMFormat f, void *data, void* parent_base, FMTypeDesc *t)
{
    switch (t->type) {
    case FMType_pointer:
    case FMType_string:
	return f->pointer_size;
    case FMType_array: {
	int size = 1;
	while (t->type == FMType_array) {
	    if (t->static_size == 0) {
		struct _FMgetFieldStruct src_spec;
		int field = t->control_field_index;
		memset(&src_spec, 0, sizeof(src_spec));
		src_spec.size = f->field_list[field].field_size;
		src_spec.offset = f->field_list[field].field_offset;
		int tmp = quick_get_ulong(&src_spec, parent_base);
		size = size * tmp;
	    } else {
		size *= t->static_size;
	    }
	    t = t->next;
	}
	size *= determine_dump_size(f, data, parent_base, t);
	return size;
    }
    case FMType_subformat:
	return f->field_subformats[t->field_index]->record_length;
    case FMType_simple:
	return f->field_list[t->field_index].field_size;
    }
    return -1;
}

/* big enough for most single values */
#define MAX_VALUE_SIZE 512


static int
sdump_value(dstate s, const char *field_type, int field_size, size_t field_offset, FMFormat top_format,
	    void *data, void *string_base, int byte_reversal, int float_format, int encode)
{
    FMgetFieldStruct descr;  /*OK */
    long junk;
    descr.offset = field_offset;
    descr.size = field_size;
    descr.data_type = FMarray_str_to_data_type(field_type, &junk);
    descr.byte_swap = byte_reversal;
    descr.src_float_format = float_format ;
    descr.target_float_format = fm_my_float_format;

    if (descr.data_type == integer_type) {
	if (field_size <= sizeof(long)) {
	    size_t tmp = get_FMlong(&descr, data);
	    dump_output(s, 25, "%ld ", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long;
	    long high_long;
	    get_FMlong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		dump_output(s, 12, "%ld ", low_long);
	    } else {
		dump_output(s, 20, "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    dump_output(s, 10, "+long int+ ");
	} else {
	    dump_output(s, 20, "+int size %d+ ", field_size);
	}
    } else if (descr.data_type == unsigned_type) {
	if (field_size <= sizeof(unsigned long)) {
	    ssize_t tmp = get_FMulong(&descr, data);
	    dump_output(s, 25, "%lu ", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long, high_long;
	    get_FMulong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		dump_output(s, 12, "%lu ", low_long);
	    } else {
		dump_output(s, 20, "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    dump_output(s, 12, "+ulong int+ ");
	} else {
	    dump_output(s, 20, "+uint size %u+ ", field_size);
	}
    } else if (descr.data_type == enumeration_type) {
	size_t tmp = get_FMulong(&descr, data);
	dump_output(s, 25, "%lu ", tmp);
    } else if (descr.data_type == boolean_type) {
	size_t tmp = get_FMulong(&descr, data);
	dump_output(s, 25, "%lu ", tmp);
	if (tmp == 0) {
	    dump_output(s, 5, "false ");
	} else {
	    dump_output(s, 4, "true ");
	}
    } else if (descr.data_type == float_type) {
	if (field_size == sizeof(float)) {
	    float tmp = get_FMfloat(&descr, data);
	    dump_output(s, 20, "%g ", tmp);
	} else if (field_size == sizeof(double)) {
	    double tmp = get_FMdouble(&descr, data);
	    dump_output(s, 20, "%g ", tmp);
#if SIZEOF_LONG_DOUBLE != 0 && SIZEOF_LONG_DOUBLE != SIZEOF_DOUBLE
	} else if (field_size == sizeof(long double)) {
	    long double tmp;
	    memcpy(&tmp, (float *) ((char *) data + field_offset),
		   sizeof(double));
	    dump_output(s, 30, "%Lg ", tmp);
#endif
	} else {
	    if (field_size < sizeof(float)) {
		dump_output(s, 12, "+tiny float+ ");
	    } else if (field_size > sizeof(double)) {
		dump_output(s, 12, "+big float+ ");
	    } else {
		dump_output(s, 20, "+float size %u+ ", field_size);
	    }
	}
    } else if (descr.data_type == char_type) {
	dump_output(s, 1, "%c ", *(char *) ((char *) data + field_offset));
    } else if (descr.data_type == string_type) {
	char *tmp_str = (char *) get_FMaddr(&descr, data, string_base, encode);
	if (tmp_str == 0) {
	    dump_output(s, 4, "\" \"");
	} else {
	    dump_output(s, strlen(tmp_str) + 3, "\"%s \"", tmp_str);
	}
    } else if (strcmp(field_type, "R3vector") == 0) {
	dump_output(s, 60, "(%g, %g, %g) ", *(double *) ((char *) data + field_offset),
	     *(double *) ((char *) data + field_offset + sizeof(double)),
	*(double *) ((char *) data + field_offset + 2 * sizeof(double)));
    } else {
	return 0;
    }
    return 1;
}

static int
dump_subfield(void*base, FMFormat f, dstate s, size_t data_offset, void* parent_base, FMTypeDesc *t)
{

    if ((s->output_limit != -1) && (s->output_len > (size_t)s->output_limit)) return 0;

    switch (t->type) {
    case FMType_pointer:
    {
	struct _FMgetFieldStruct src_spec;
	size_t new_offset;
	char *ptr_value;
	memset(&src_spec, 0, sizeof(src_spec));
	src_spec.size = f->pointer_size;
	src_spec.byte_swap = f->byte_reversal;
	ptr_value = quick_get_pointer(&src_spec, (char*)base + data_offset);
	if (ptr_value == NULL) {
	    dump_output(s, 4, "NULL");
	    return 1;
	} else {
	    dump_output(s, 20, "%p ", ptr_value);
	}
	if (s->encoded) {
#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 810)
#endif
	  new_offset =  (intptr_t) ptr_value;
#if defined (__INTEL_COMPILER)
#  pragma warning (default: 810)
#endif
	    ptr_value = (intptr_t)ptr_value + s->offset_base;
	} else {
	    new_offset = (intptr_t)ptr_value - (intptr_t)s->offset_base;
	}
	if (f->recursive) {
	    ssize_t previous_offset = search_addr_list(s, ptr_value);
	    if (previous_offset != -1) {
		/* already visited this */
		return 1;
	    }
	}

	if (f->recursive) {
	    add_to_addr_list(s, ptr_value, new_offset);
	}
	return dump_subfield(ptr_value, f, s, 0, parent_base, t->next);
    }
    case FMType_string:
    {
	struct _FMgetFieldStruct src_spec;
	char *ptr_value;
	memset(&src_spec, 0, sizeof(src_spec));
	src_spec.size = f->pointer_size;
	src_spec.byte_swap = f->byte_reversal;
	ptr_value = quick_get_pointer(&src_spec, (char*)base + data_offset);
	if (ptr_value == NULL) {
	    dump_output(s, 5, "NULL ");
	} else {
	    if (s->encoded) {
		ptr_value = (intptr_t)ptr_value + s->offset_base;
	    }
	    dump_output(s, strlen(ptr_value) + 2, "\"%s\"", ptr_value);
	}
	break;
    }
    case FMType_array:
    {
	int elements = 1, i;
	int element_size;
	FMTypeDesc *next = t;
	while (next->type == FMType_array) {
	    if (next->static_size == 0) {
		struct _FMgetFieldStruct src_spec;
		int field = next->control_field_index;
		memset(&src_spec, 0, sizeof(src_spec));
		src_spec.size = f->field_list[field].field_size;
		src_spec.offset = f->field_list[field].field_offset;
		int tmp = quick_get_ulong(&src_spec, parent_base);
		elements = elements * tmp;
	    } else {
		elements = elements * next->static_size;
	    }
	    next = next->next;
	}
	element_size = determine_dump_size(f, base, parent_base, next);
	for (i = 0; i < elements ; i++) {
	    size_t element_offset = data_offset + i * element_size;
	    if (!dump_subfield(base, f, s, element_offset, parent_base, next)) return 0;
	}
	break;
    }
    case FMType_subformat:
    {
	int field_index = t->field_index;
	FMFormat subformat = f->field_subformats[field_index];
	int ret;
	FMFieldList fmfield = &f->field_list[field_index];
	start_field(s, fmfield, &f->var_list[field_index].type_desc);
	ret = dump_subfields(base, subformat, s, data_offset);
	stop_field(s, fmfield, &f->var_list[field_index].type_desc);
	return ret;
    }
    case FMType_simple: {
	FMFieldList fmfield = &f->field_list[t->field_index];
	size_t field_offset = data_offset;
	int field_size = fmfield->field_size;
	const char *field_type = fmfield->field_type;
	int byte_reversal = f->byte_reversal;
	int float_format = f->float_format;
	sdump_value(s, field_type, field_size, field_offset, (FMFormat)NULL, base,
		    base, byte_reversal, float_format, s->encoded);
	break;
    }
    default:
	assert(0);
    }
    return 1;
}

static void
add_to_addr_list(dstate s, void *addr, size_t offset)
{
    if (s->addr_list_is_stack) {
	if (s->addr_list_cnt == STACK_ARRAY_SIZE) {
	    /* reached the max size for stack-based addr list */
	    addr_list_entry *malloc_list;
	    s->addr_list_is_stack = 0;
	    s->malloc_addr_size = 2*STACK_ARRAY_SIZE;
	    malloc_list = malloc(2*STACK_ARRAY_SIZE * sizeof(addr_list_entry));
	    memcpy(malloc_list, s->addr_list, STACK_ARRAY_SIZE * sizeof(addr_list_entry));
	    s->addr_list = malloc_list;
	}
    } else {
       /* malloc'd addr_list */
       if (s->addr_list_cnt == s->malloc_addr_size) {
	   s->malloc_addr_size *= 2;
	   s->addr_list = 
	       realloc(s->addr_list, sizeof(addr_list_entry)*s->malloc_addr_size);
       }
    }
    s->addr_list[s->addr_list_cnt].addr = addr;
    s->addr_list[s->addr_list_cnt].offset = offset;
    s->addr_list_cnt++;
}

static ssize_t
search_addr_list(dstate s, void *addr)
{
    int i;
    ssize_t previous_offset = -1;
    for (i=0; i < s->addr_list_cnt; i++) {
	if (s->addr_list[i].addr == addr) {
	    previous_offset = s->addr_list[i].offset;
	}
    }
    return previous_offset;
}

static void
free_addr_list(dstate s)
{
    if (s->addr_list_is_stack == 0) {
	free(s->addr_list);
	s->addr_list = NULL;
    }
}

extern int
dump_raw_FMrecord(FMContext fmc, FMFormat format, void *data)
{
    struct dump_state state;
    init_dump_state(&state);
    state.output_limit = -1;
    state.encoded = 1;
    state.offset_base = data;
    if (FMdumpVerbose)
	dump_output(&state, strlen(format->format_name) + 15, "Record type %s :", format->format_name);
    internal_dump_data(format, data, &state);
    dump_output(&state, 1, "\n");
    return 0;
}

extern char *
dump_raw_FMrecord_to_string(FMContext fmc, FMFormat format, void *data)
{
    struct dump_state state;
    init_dump_state(&state);
    state.output_limit = -1;
    state.encoded = 1;
    state.offset_base = data;
    state.use_file_out = 0;
    state.output_string = malloc(1);
    state.realloc_string++;
    if (FMdumpVerbose)
	dump_output(&state, strlen(format->format_name) + 15, "Record type %s :", format->format_name);
    (void)internal_dump_data(format, data, &state);
    dump_output(&state, 1, "\n");
    return state.output_string;
}

extern int
FMfdump_encoded_data(void *out, FMFormat format, void *data, int character_limit)
{
    int ret;
    int header_size = format->server_ID.length;
    if (format->variant) {
	if (format->IOversion <= 3) 
	    header_size += sizeof(INT4);
	else 
	    header_size += 8;
    }
    header_size += (8 - header_size) & 0x7;
    data = (char*)data + header_size;

    struct dump_state state;
    init_dump_state(&state);
    state.output_limit = character_limit;
    state.encoded = 1;
    state.offset_base = data;
    state.out = (FILE*)out;
    if (FMdumpVerbose)
	dump_output(&state, strlen(format->format_name) + 15, "Record type %s :", format->format_name);
    ret = internal_dump_data(format, data, &state);
    dump_output(&state, 1, "\n");
    return ret;
}

extern int
FMdump_encoded_data(FMFormat format, void *data, int character_limit)
{
    return FMfdump_encoded_data((void*)stdout, format, data, character_limit);
}

void *
FMheader_skip(FMContext c, void *data)
{
    FMFormat format = FMformat_from_ID(c, data);
    int header_size = format->server_ID.length;

    if (format->variant) {
	if (format->IOversion <= 3) 
	    header_size += sizeof(INT4);
	else 
	    header_size += 8;
    }
    header_size += (8 - header_size) & 0x7;
    data = (char*)data + header_size;
    return data;
}

extern void
FMfdump_encoded_XML(void *out, FMContext c, void *data, int limit)
{
    FMFormat format = FMformat_from_ID(c, data);
    int header_size = format->server_ID.length;
    struct dump_state state;

    if (format->variant) {
	if (format->IOversion <= 3) 
	    header_size += sizeof(INT4);
	else 
	    header_size += 8;
    }
    header_size += (8 - header_size) & 0x7;
    data = (char*)data + header_size;

    
    if (FMhas_XML_info(format)) {
	FMdump_XML(format, data, 1);
	return;
    }
    init_dump_state(&state);
    state.output_limit = -1;
    state.encoded = 1;
    state.offset_base = data;
    state.use_XML = 1;
    state.out = (FILE*)out;
    dump_output(&state,	strlen(format->format_name) + 3, "<%s>\n", format->format_name);
    internal_dump_data(format, data, &state);
    dump_output(&state, strlen(format->format_name) + 4, "</%s>\n", format->format_name);
}

extern void
FMdump_encoded_XML(FMContext c, void *data, int limit)
{
    FMfdump_encoded_XML(stdout, c, data, limit);
}

extern void
dump_unencoded_FMrecord_as_XML(FMContext fmc, FMFormat format, void *data)
{
    struct dump_state state;
    if (FMhas_XML_info(format)) {
	FMdump_XML(format, data, 0);
	return;
    }
    init_dump_state(&state);
    state.output_limit = -1;
    state.encoded = 0;
    state.use_XML = 1;
    dump_output(&state, strlen(format->format_name) + 4, "<%s>\n", format->format_name);
    internal_dump_data(format, data, &state);
    dump_output(&state, strlen(format->format_name) + 4, "</%s>\n", format->format_name);
}

