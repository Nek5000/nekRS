
#include "config.h"

#include "assert.h"
#include "ffs.h"
#include "cod.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "fm_internal.h"
#include "ffs_internal.h"
#include "ffs_marshal.h"

static void *
quick_get_pointer(FMFieldPtr iofield, void *data);
static size_t
quick_get_ulong(FMFieldPtr iofield, void *data);
void
quick_put_ulong(FMFieldPtr iofield, size_t value, void *data);

static ssize_t add_to_tmp_buffer(FFSBuffer buf, size_t size);
static int64_t
final_variant_size_for_record(int64_t input_record_len, IOConversionPtr conv);

static void byte_swap(char *data, int size);

typedef struct addr_list {
    void *addr;
    ssize_t offset;
} addr_list_entry;

typedef struct encode_state {
    int copy_all;
    int64_t output_len;
    int iovec_is_stack;
    int iovcnt;
    internal_iovec *iovec;
    int malloc_vec_size;
    int no_leaf_copy;

    int addr_list_is_stack;
    int addr_list_cnt;
    addr_list_entry *addr_list;
    ssize_t malloc_addr_size;
    ssize_t saved_offset_difference;
    void *orig_data;
}*estate;

static void free_addr_list(estate s);

void
setup_header(FFSBuffer buf, FMFormat f, estate s) 
{
    int header_size = f->server_ID.length;
    ssize_t tmp_data;	/* offset for header */
    int align_pad;

    if (f->variant) {
	header_size += 8;	/* length info */
    }
    align_pad = (8 - header_size) & 0x7;	/* align to 8 */
    header_size += align_pad;

    tmp_data = add_to_tmp_buffer(buf, header_size);

    memcpy((char *) buf->tmp_buffer + tmp_data, f->server_ID.value,
	   f->server_ID.length);

    memset((char*)buf->tmp_buffer + tmp_data + f->server_ID.length,
	   0, header_size - f->server_ID.length);

    /* fill in an IOV field for the header */
    s->iovec[0].iov_len = header_size;
    s->iovec[0].iov_offset = tmp_data;
    s->iovec[0].iov_base = NULL;	/* offset is in tmp_data */
    s->iovcnt++;
    s->output_len = header_size;
}


static char zeroes[8] = {0, 0, 0, 0, 0, 0, 0, 0};
#define STACK_ARRAY_SIZE 100

static void
ensure_writev_room(estate s, int add_count)
{
    /* 
     * we need to track the length of the strings we're writing
     * assume a small number (STACK_ARRAY_SIZE) but expand
     * dynamically if we have more...  
     */
    if (!(s->malloc_vec_size) &&
	(s->iovcnt >= STACK_ARRAY_SIZE - add_count)) {
	int j;
	internal_iovec *tmpl = (internal_iovec *)
	    malloc(sizeof(internal_iovec) * (2 * STACK_ARRAY_SIZE + add_count));
	s->malloc_vec_size = (2 * STACK_ARRAY_SIZE + add_count);
	for (j = 0; j < s->iovcnt; j++) {
	    tmpl[j].iov_len = (s->iovec)[j].iov_len;
	    tmpl[j].iov_base = (s->iovec)[j].iov_base;
	    tmpl[j].iov_offset = (s->iovec)[j].iov_offset;
	}
	s->iovec = tmpl;
	s->iovec_is_stack = 0;
    } else if ((s->malloc_vec_size != 0) &&
	       (s->iovcnt >= s->malloc_vec_size - add_count)) {
	(s->malloc_vec_size) *= 2;
	s->iovec = realloc(s->iovec, sizeof(internal_iovec) *
				(s->malloc_vec_size));
    }
}


size_t
allocate_tmp_space(estate s, FFSBuffer buf, size_t length, int req_alignment, size_t *tmp_data_loc)
{
    int pad = (req_alignment - s->output_len) & (req_alignment -1);  /*  only works if req_align is power of two */
    size_t tmp_data, msg_offset;
    switch (req_alignment) {
    case 1: case 2: case 4: case 8: case 16: break;
    default:
	assert(0);
    }
    ensure_writev_room(s, 2);
    tmp_data = add_to_tmp_buffer(buf, length + pad);
    if (tmp_data == (size_t)-1) return (size_t) -1;
    if (pad != 0) {
	if (s->iovec[s->iovcnt-1].iov_base == NULL) {
	    /* last was tmp too */
	    memset((char*)buf->tmp_buffer + tmp_data, 0, pad);	/* zero pad */
	    tmp_data += pad;
	    s->iovec[s->iovcnt-1].iov_len += pad;
	} else {
	    s->iovec[s->iovcnt].iov_len = pad;
	    s->iovec[s->iovcnt].iov_offset = 0;
	    s->iovec[s->iovcnt].iov_base = zeroes;
	    s->iovcnt++;
	}
    }
    msg_offset = s->output_len + pad;
    if (tmp_data_loc) *tmp_data_loc = tmp_data;
    s->output_len += length + pad;
    return msg_offset;
}

size_t
copy_data_to_tmp(estate s, FFSBuffer buf, void *data, size_t length, int req_alignment, size_t *tmp_data_loc)
{
    size_t tmp_data;
    size_t msg_offset = allocate_tmp_space(s, buf, length, req_alignment, &tmp_data);
    if (length != 0) {
	memcpy((char*)buf->tmp_buffer + tmp_data, data, length);
	s->iovec[s->iovcnt].iov_len = length;
	s->iovec[s->iovcnt].iov_offset = tmp_data;
	s->iovec[s->iovcnt].iov_base = NULL;
	s->iovcnt++;
    }
    if (tmp_data_loc) *tmp_data_loc = tmp_data;
    return msg_offset;
}

/*
 *  same as copy_data_to_tmp, except doesn't actually copy the data.  This is used for FFS_NO_LEAF_COPY.
 */
size_t
reserve_space_for_data_in_tmp(estate s, FFSBuffer buf, void *data, size_t length, int req_alignment, size_t *tmp_data_loc)
{
    size_t tmp_data;
    size_t msg_offset = allocate_tmp_space(s, buf, length, req_alignment, &tmp_data);
    if (length != 0) {
	s->iovec[s->iovcnt].iov_len = length;
	s->iovec[s->iovcnt].iov_offset = tmp_data;
	s->iovec[s->iovcnt].iov_base = NULL;
	s->iovcnt++;
    }
    if (tmp_data_loc) *tmp_data_loc = tmp_data;
    return msg_offset;
}

size_t
add_data_iovec(estate s, FFSBuffer buf, void *data, size_t length, int req_alignment)
{
    size_t msg_offset;
    int pad = (req_alignment - s->output_len) & (req_alignment -1);  /*  only works if req_align is power of two */
    switch (req_alignment) {
    case 1: case 2: case 4: case 8: case 16: break;
    default:
	assert(0);
    }
    ensure_writev_room(s, 2);
    if (pad != 0) {
	s->iovec[s->iovcnt].iov_len = pad;
	s->iovec[s->iovcnt].iov_offset = 0;
	s->iovec[s->iovcnt].iov_base = zeroes;
	s->iovcnt++;
	s->output_len += pad;
    }
    if (length != 0) {
	s->iovec[s->iovcnt].iov_len = length;
	s->iovec[s->iovcnt].iov_offset = 0;
	s->iovec[s->iovcnt].iov_base = data;
	s->iovcnt++;
    }
    msg_offset = s->output_len;
    s->output_len += length;
    return msg_offset;
}

void
init_encode_state(estate state)
{
    state->output_len = 0;
    state->iovcnt = 0;
    state->malloc_vec_size = 0;
    state->addr_list_cnt = 0;
    state->no_leaf_copy = 0;
}

static int
handle_subfields(FFSBuffer buf, FMFormat f, estate s, size_t data_offset);

static char *
FFSencode_internal(FFSBuffer b, FMFormat fmformat, void *data, size_t *buf_size, int flags)
{
    internal_iovec stack_iov_array[STACK_ARRAY_SIZE];
    addr_list_entry stack_addr_list[STACK_ARRAY_SIZE];
    struct encode_state state;
    init_encode_state(&state);
    size_t base_offset = 0;
    size_t header_size;

    state.iovec_is_stack = 1;
    state.iovec = stack_iov_array;
    state.addr_list_is_stack = 1;
    state.addr_list_cnt = 0;
    state.addr_list = stack_addr_list;
    state.copy_all = 1;
    state.saved_offset_difference = 0;
    state.orig_data = data;
    if ((flags & FFS_NO_LEAF_COPY) == FFS_NO_LEAF_COPY) {
	state.no_leaf_copy = 1;
    }

    make_tmp_buffer(b, 0);

    /* setup header information */
    setup_header(b, fmformat, &state);

    header_size = state.output_len;
    state.saved_offset_difference = header_size;

    if (fmformat->variant || state.copy_all) {
	base_offset = copy_data_to_tmp(&state, b, data, 
				       fmformat->record_length, 1, NULL);
	if (base_offset == (size_t) -1) return NULL;
    }

    if (!fmformat->variant) {
	*buf_size = state.output_len;
	return b->tmp_buffer;
    }

    if (fmformat->recursive) {
	state.addr_list[state.addr_list_cnt].addr = data;
	state.addr_list[state.addr_list_cnt].offset = base_offset;
	state.addr_list_cnt++;
    }

    copy_data_to_tmp(&state, b, NULL, 0, 8, NULL);   /* force 64-bit align for var */

    if (!handle_subfields(b, fmformat, &state, base_offset)) return NULL;

    {
	/* fill in actual length of data marshalled after header */
	char *tmp_data = b->tmp_buffer;
	uint64_t record_len = state.output_len - header_size;
	tmp_data += fmformat->server_ID.length;
	memcpy(tmp_data, &record_len, 8);
    }
    free_addr_list(&state);
    *buf_size = state.output_len;
    if (!state.iovec_is_stack) {
	free(state.iovec);
    }
    return b->tmp_buffer;
}

char *
FFSencode(FFSBuffer b, FMFormat fmformat, void *data, size_t *buf_size)
{
    return FFSencode_internal(b, fmformat, data, buf_size, /*flags*/ 0);
}

char *
FFSencode_no_leaf_copy(FFSBuffer b, FMFormat fmformat, void *data, size_t *buf_size)
{
    return FFSencode_internal(b, fmformat, data, buf_size, /*flags*/ FFS_NO_LEAF_COPY);
}

static FFSEncodeVector
fixup_output_vector(FFSBuffer b, estate s)
{
    int size = (s->iovcnt + 5) * sizeof(struct FFSEncodeVec);
    ssize_t tmp_vec_offset = add_to_tmp_buffer(b, size);
    FFSEncodeVector ret;
    /* buffer will not be realloc'd after this */
    int i;

    if (tmp_vec_offset == -1) return NULL;

    tmp_vec_offset += sizeof(ret[0]) - 1;
    tmp_vec_offset -= tmp_vec_offset % sizeof(ret[0]);
    ret = (FFSEncodeVector)(((char*)b->tmp_buffer) + tmp_vec_offset);
    /* leave two blanks, see notes in ffs_file.c */
    /* the upshot is that we use these if we need to add headers */
    ret += 3;
    for (i=0; i < s->iovcnt; i++) {
	ret[i].iov_len = s->iovec[i].iov_len;
	if (s->iovec[i].iov_base != NULL) {
	    ret[i].iov_base = s->iovec[i].iov_base;
	} else {
	    ret[i].iov_base = ((char*)b->tmp_buffer) + s->iovec[i].iov_offset;
	}
    }
    ret[s->iovcnt].iov_base = NULL;
    ret[s->iovcnt].iov_len = 0;
    if (!s->iovec_is_stack) {
	free(s->iovec);
	s->iovec = NULL;
    }
    if (!s->addr_list_is_stack && (s->addr_list != NULL)) {
	free(s->addr_list);
	s->addr_list = NULL;
    }
    return ret;
}

static void
add_to_addr_list(estate s, void *addr, size_t offset)
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
search_addr_list(estate s, void *addr)
{
    int i;
    ssize_t previous_offset = (ssize_t) -1;
    for (i=0; i < s->addr_list_cnt; i++) {
	if (s->addr_list[i].addr == addr) {
	    previous_offset = s->addr_list[i].offset;
	}
    }
    return previous_offset;
}

static void
free_addr_list(estate s)
{
    if (s->addr_list_is_stack == 0) {
	free(s->addr_list);
	s->addr_list = NULL;
    }
}

FFSEncodeVector
FFSencode_vector(FFSBuffer b, FMFormat fmformat, void *data)
{
    internal_iovec stack_iov_array[STACK_ARRAY_SIZE];
    addr_list_entry stack_addr_list[STACK_ARRAY_SIZE];
    struct encode_state state;
    init_encode_state(&state);
    ssize_t base_offset = 0;
    ssize_t header_size;

    state.iovec_is_stack = 1;
    state.iovec = stack_iov_array;
    state.addr_list_is_stack = 1;
    state.addr_list = stack_addr_list;
    state.copy_all = 0;
    state.saved_offset_difference = 0;
    state.orig_data = data;

    make_tmp_buffer(b, 0);

    /* setup header information */
    setup_header(b, fmformat, &state);

    header_size = state.output_len;
    state.saved_offset_difference = header_size;

    if (fmformat->variant || state.copy_all) {
	base_offset = copy_data_to_tmp(&state, b, data, 
				       fmformat->record_length, 1, NULL);
	if (base_offset == -1) return NULL;
    } else if (!state.copy_all) {
	base_offset = add_data_iovec(&state, b, data, 
				     fmformat->record_length, 1);
    }
    if (!fmformat->variant) {
	return fixup_output_vector(b, &state);
    }

    if (fmformat->recursive) {
	state.addr_list[state.addr_list_cnt].addr = data;
	state.addr_list[state.addr_list_cnt].offset = base_offset;
	state.addr_list_cnt++;
    }

    copy_data_to_tmp(&state, b, NULL, 0, 8, NULL);   /* force 64-bit align for var */

    if (!handle_subfields(b, fmformat, &state, base_offset)) return NULL;
    {
	/* fill in actual length of data marshalled after header */
	char *tmp_data = b->tmp_buffer;
	int64_t record_len = state.output_len - header_size;
	tmp_data += fmformat->server_ID.length;
	memcpy(tmp_data, &record_len, 8);
    }
    free_addr_list(&state);
    return fixup_output_vector(b, &state);
}

static int
handle_subfield(FFSBuffer buf, FMFormat f, estate s, size_t data_offset, size_t parent_offset, FMTypeDesc *t);

extern int
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

static int
handle_subfields(FFSBuffer buf, FMFormat f, estate s, size_t data_offset)
{
    int i;
    /* if base is not variant (I.E. doesn't contain addresses), return;*/
    if (!f->variant) return 1;

    for (i = 0; i < f->field_count; i++) {
	size_t subfield_offset = data_offset + f->field_list[i].field_offset;
	int ret;
	if (field_is_flat(f, &f->var_list[i].type_desc)) continue;
	ret = handle_subfield(buf, f, s, subfield_offset, data_offset, 
			      &f->var_list[i].type_desc);
	if (ret != 1) return 0;
    }
    return 1;
}

static void
set_dynamic_array_size(FMFormat f, FFSBuffer buf, size_t parent_offset, FMTypeDesc *t, int new_value)
{
    struct _FMgetFieldStruct src_spec;
    int field = t->control_field_index;
    if ((t->type != FMType_array) || (t->static_size != 0)) {
	printf("Set array size failed!\n");
	return;
    }
    memset(&src_spec, 0, sizeof(src_spec));
    src_spec.size = f->field_list[field].field_size;
    src_spec.offset = f->field_list[field].field_offset;
    quick_put_ulong(&src_spec, new_value,
		    (char*)buf->tmp_buffer + parent_offset);
}


static size_t
determine_size(FMFormat f, FFSBuffer buf, size_t parent_offset, FMTypeDesc *t)
{
    switch (t->type) {
    case FMType_pointer:
    case FMType_string:
	return f->pointer_size;
    case FMType_array: {
	size_t size = 1;
	while (t->type == FMType_array) {
	    if (t->static_size == 0) {
		struct _FMgetFieldStruct src_spec;
		int field = t->control_field_index;
		memset(&src_spec, 0, sizeof(src_spec));
		src_spec.size = f->field_list[field].field_size;
		src_spec.offset = f->field_list[field].field_offset;
		size_t tmp = quick_get_ulong(&src_spec, (char*)buf->tmp_buffer + parent_offset);
		size = size * tmp;
	    } else {
		size *= t->static_size;
	    }
	    t = t->next;
	}
	size *= determine_size(f, buf, parent_offset, t);
	return size;
    }
    case FMType_subformat:
	return f->field_subformats[t->field_index]->record_length;
    case FMType_simple:
	return f->field_list[t->field_index].field_size;
    }
    assert(FALSE);
    return (size_t)-1;
}

static int
handle_subfield(FFSBuffer buf, FMFormat f, estate s, size_t data_offset, size_t parent_offset, FMTypeDesc *t)
{

    switch (t->type) {
    case FMType_pointer:
    {
	struct _FMgetFieldStruct src_spec;
	size_t size, new_offset, tmp_data_loc;
	char *ptr_value;
	field_marshal_info marshal_info;

	memset(&src_spec, 0, sizeof(src_spec));
	src_spec.size = f->pointer_size;
	ptr_value = quick_get_pointer(&src_spec, (char*)buf->tmp_buffer + data_offset);
	if (ptr_value == NULL) return 1;

	/* customized marshalling */
	if ((marshal_info = get_marshal_info(f, t)) != NULL) {
	    if (marshal_info->type == FFSDropField) {
		if (marshal_info->drop_row_func(s->orig_data)) {
		    /* drop the value */
		    quick_put_ulong(&src_spec, 0,
				    (char*)buf->tmp_buffer + data_offset);
		    return 1;
		}
	    }
	}

	if (f->recursive) {
	    ssize_t previous_offset = search_addr_list(s, ptr_value);
	    if (previous_offset != (ssize_t)-1) {
		quick_put_ulong(&src_spec, previous_offset,
				(char*)buf->tmp_buffer + data_offset);
		return 1;
	    }
	}
	size = determine_size(f, buf, parent_offset, t->next);
	if (size == 0)  {
	    quick_put_ulong(&src_spec, 0,
			    (char*)buf->tmp_buffer + data_offset);
	    return 1;
	}
	if ((marshal_info == NULL) || (marshal_info->type != FFSSubsampleArrayField)) {
	    if (!s->copy_all && field_is_flat(f, t->next)) {
		/* leave data where it sits */
		new_offset = add_data_iovec(s, buf, ptr_value, size, 8);
	    } else {
		if (s->no_leaf_copy) {
		    /* this path only used for FFSencode_no_leaf_copy() */
		    new_offset = reserve_space_for_data_in_tmp(s, buf, ptr_value, size, 8, &tmp_data_loc);
		} else {
		    new_offset = copy_data_to_tmp(s, buf, ptr_value, size, 8, &tmp_data_loc);
		}
	    }
	} else {
	    /* can't leave data where it sits */
	    /* 1. allocate temporary space for the destination array
	       2. call function to fill destination with part of source
	       3. return value is number of elements copied
	       4. adjust size of temporary to reflect actual data copied
	    */
	    size_t element_size = determine_size(f, buf, parent_offset, t->next->next);
	    size_t element_count = size / element_size;
	    cod_exec_context ec = marshal_info->ec;
	    struct subsample_marshal_data smd;
	    new_offset = allocate_tmp_space(s, buf, size, 8, &tmp_data_loc);
	    smd.element_count = (int)element_count;
	    smd.element_size = element_size;
	    smd.src_ptr = ptr_value;
	    smd.dst_ptr = (char*)buf->tmp_buffer + new_offset;
	    smd.marshalled_count = 0;
#ifdef DO_DCG
	    cod_assoc_client_data(ec,  0x534d4450, (intptr_t)&smd);
#endif
	    marshal_info->subsample_array_func(ec, s->orig_data, (int)element_count);
	    /* fixup size */
	    set_dynamic_array_size(f, buf, parent_offset, t->next, 
				   smd.marshalled_count);
	    if (size != 0) {
		s->iovec[s->iovcnt].iov_len = smd.marshalled_count * element_size;
		s->iovec[s->iovcnt].iov_offset = tmp_data_loc;
		s->iovec[s->iovcnt].iov_base = NULL;
		s->iovcnt++;
	    }
	    s->output_len -= (size - (smd.marshalled_count *element_size));
	    buf->tmp_buffer_in_use_size -= (size - (smd.marshalled_count *element_size));
	}
	quick_put_ulong(&src_spec, new_offset - s->saved_offset_difference, 
			(char*)buf->tmp_buffer + data_offset);
	if (f->recursive) {
	    add_to_addr_list(s, ptr_value, new_offset - s->saved_offset_difference);
	}
	if (field_is_flat(f, t->next)) return 1;
	return handle_subfield(buf, f, s, tmp_data_loc, parent_offset, t->next);
    }
    case FMType_string:
    {
	struct _FMgetFieldStruct src_spec;
	char *ptr_value;
	size_t size, str_offset;
	memset(&src_spec, 0, sizeof(src_spec));
	src_spec.size = f->pointer_size;
	ptr_value = quick_get_pointer(&src_spec, (char*)buf->tmp_buffer + data_offset);
	if (ptr_value == NULL) return 1;
	size = strlen(ptr_value) + 1;
	if (!s->copy_all) {
	    /* leave data where it sits */
	    str_offset = add_data_iovec(s, buf, ptr_value, size, 1);
	} else {
	    str_offset = copy_data_to_tmp(s, buf, ptr_value, size, 1, NULL);
	    if (str_offset == (size_t) -1) return 0;
	}
	quick_put_ulong(&src_spec, str_offset - s->saved_offset_difference,
			(char*)buf->tmp_buffer + data_offset);
	break;
    }
    case FMType_array:
    {
	size_t elements = 1, i;
	size_t element_size;
	FMTypeDesc *next = t;
	while (next->type == FMType_array) {
	    if (next->static_size == 0) {
		struct _FMgetFieldStruct src_spec;
		int field = next->control_field_index;
		memset(&src_spec, 0, sizeof(src_spec));
		src_spec.size = f->field_list[field].field_size;
		src_spec.offset = f->field_list[field].field_offset;
		size_t tmp = quick_get_ulong(&src_spec, (char*)buf->tmp_buffer + parent_offset);
		elements = elements * tmp;
	    } else {
		elements = elements * next->static_size;
	    }
	    next = next->next;
	}
	element_size = determine_size(f, buf, parent_offset, next);
	if (field_is_flat(f, next)) return 1;
	for (i = 0; i < elements ; i++) {
	    size_t element_offset = data_offset + i * element_size;
	    if (!handle_subfield(buf, f, s, element_offset, parent_offset, next)) return 0;
	}
	break;
    }
    case FMType_subformat:
    {
	int field_index = t->field_index;
	FMFormat subformat = f->field_subformats[field_index];
	return handle_subfields(buf, subformat, s, data_offset);
    }
    default:
	assert(0);
    }
    return 1;
}

extern FFSEncodeVector
copy_vector_to_FFSBuffer(FFSBuffer buf, FFSEncodeVector vec)
{
    ssize_t vec_offset = (char *) vec - (char *)buf->tmp_buffer;

    if (((char*)vec < (char*)buf->tmp_buffer)
        || ((char*)vec >= (char*)buf->tmp_buffer + buf->tmp_buffer_size)) {
        size_t i, remainder;
        for (i = 0; vec[i].iov_base; ++i);
        vec_offset = add_to_tmp_buffer(buf, (i + 2) * sizeof(*vec));
	remainder = vec_offset % sizeof(*vec);
	if (remainder != 0) {
	    vec_offset += sizeof(*vec) - remainder;
	}
        memcpy((char *) buf->tmp_buffer + vec_offset, vec, (i + 1) * sizeof(*vec));
    }

    return (FFSEncodeVector) ((char *) buf->tmp_buffer + vec_offset);
}

extern FFSEncodeVector
copy_all_to_FFSBuffer(FFSBuffer buf, FFSEncodeVector vec)
{
    int i = 0;
    intptr_t vec_offset = (intptr_t) vec - (intptr_t)buf->tmp_buffer;
    /* 
     * vec and some of the buffers may be in the memory managed by the
     * FFSBuffer.  The goal here to is put *everything* into the FFSBuffer.
     */
    int vec_count = 0;
    while (vec[vec_count].iov_base != NULL) {
      vec_count++;
    }
    assert(((uintptr_t)vec >= (uintptr_t)buf->tmp_buffer) && 
	   ((uintptr_t)vec < (uintptr_t)buf->tmp_buffer + buf->tmp_buffer_size));
    {
	size_t* already_in = malloc(sizeof(already_in[0]) * vec_count);
	while (vec[i].iov_base != NULL) {
	    if (((char*)vec[i].iov_base >= (char*)buf->tmp_buffer) &&
		((char*)vec[i].iov_base < (char*)buf->tmp_buffer + buf->tmp_buffer_size)) {
	        /* 
		 * remap pointers into temp so that they're offsets (must do
		 * this before we realloc the temp
		 */ 
	        vec[i].iov_base = (void*)((char *)vec[i].iov_base - (char *)buf->tmp_buffer + 1);
		already_in[i] = 1;
	    } else {
		already_in[i] = 0;
	    }
	    i++;
	}

	i = 0;
	while (((FFSEncodeVector)((intptr_t)buf->tmp_buffer + vec_offset))[i].iov_base !=
	       NULL) {
	    FFSEncodeVector v = (void*)((intptr_t) buf->tmp_buffer + vec_offset);
	    if (already_in[i] == 0) {
	        /* if this is an external buffer, copy it */
	        ssize_t offset = add_to_tmp_buffer(buf, v[i].iov_len);
		char * data = (char *) buf->tmp_buffer + offset;
		/* add_to_tmp_buffer() might have remapped vector */
		v = (void*)((char *) buf->tmp_buffer + vec_offset);
		memcpy(data, v[i].iov_base, v[i].iov_len);
		v[i].iov_base = (void*)(intptr_t)(offset + 1);
	    }
	    i++;
	}
	free(already_in);
    }
    /* reallocation done now */
    vec = (void*)((intptr_t)buf->tmp_buffer + vec_offset);
    i = 0;
    while (vec[i].iov_base != NULL) {
	if (((intptr_t)vec[i].iov_base > 0) &&
	    ((uintptr_t)vec[i].iov_base <= (uintptr_t)buf->tmp_buffer_size)) {
	    /* 
	     * remap pointers into temp so that they're addresses
	     */ 
	    vec[i].iov_base = (void*)((intptr_t)vec[i].iov_base + (char *)buf->tmp_buffer - 1);
	}
	i++;
    }
    return vec;
}

void
FFSfree_conversion(IOConversionPtr conv);

extern
void
free_FFSTypeHandle(FFSTypeHandle f)
{
    int i = 0;
    if (f->conversion) FFSfree_conversion(f->conversion);
    while(f->subformats && f->subformats[i]) {
	free_FFSTypeHandle(f->subformats[i]);
	f->subformats[i] = NULL;
	i++;
    }
    free(f->subformats);
    free(f->field_subformats);
    free(f);
}

extern
FFSTypeHandle
FFSTypeHandle_by_index(FFSContext c, int index)
{
    FFSTypeHandle handle = NULL;
    if (c->handle_list == NULL) {
	c->handle_list = malloc(sizeof(c->handle_list[0]) * (index+1));
	memset(c->handle_list, 0, sizeof(c->handle_list[0]) * (index+1));
	c->handle_list_size = index + 1;
    } else if (c->handle_list_size <= index) {
	c->handle_list = realloc(c->handle_list,
				 (index+1)*sizeof(c->handle_list[0]));
	memset(&c->handle_list[c->handle_list_size], 0, 
	       sizeof(c->handle_list[0]) * ((index+1) - c->handle_list_size));
	c->handle_list_size = index + 1;
    }
    if (c->handle_list[index] == NULL) {
	FMFormat fmf = FMformat_by_index(c->fmc, index);
	if (fmf == NULL) return NULL;
	c->handle_list[index] = malloc(sizeof(struct _FFSTypeHandle));
	handle = c->handle_list[index];
	handle->context = c;
	handle->format_id = index;
	handle->conversion = NULL;
	handle->status = not_checked;
	handle->body = FMformat_by_index(c->fmc, index);
	handle->is_fixed_target = 0;
	if ((fmf->subformats && (fmf->subformats[0] != NULL)) || fmf->recursive) {
	    int i, k, subformat_count = 0;
	    while (fmf->subformats[subformat_count] != NULL) subformat_count++;
	    handle->subformats = 
		malloc(sizeof(FFSTypeHandle)*(subformat_count+1));
	    for (i = 0 ; i < subformat_count ; i++) {
		handle->subformats[i] = malloc(sizeof(*handle->subformats[i]));
		handle->subformats[i]->context = c;
		handle->subformats[i]->format_id = -1;
		handle->subformats[i]->conversion = NULL;
		handle->subformats[i]->status = not_checked;
		handle->subformats[i]->body = fmf->subformats[i];
		handle->subformats[i]->subformats = NULL;
	    }
	    handle->subformats[subformat_count] = NULL;
	    handle->field_subformats = 
		malloc(fmf->field_count * sizeof(FFSTypeHandle));
	    memset(handle->field_subformats, 0, 
		   fmf->field_count * sizeof(FFSTypeHandle));
	    for (i = 0; i < fmf->field_count; i++) {
		if (fmf->field_subformats[i] != NULL) {
		    int j;
		    for (j = 0; j < subformat_count; j++) {
			if (fmf->field_subformats[i] == handle->subformats[j]->body) {
			    handle->field_subformats[i] = handle->subformats[j];
			}
		    }
		    if (fmf->field_subformats[i] == fmf) {
			/* recursive */
			handle->field_subformats[i] = handle;
		    }
		} else {
		    handle->field_subformats[i] = NULL;
		}
	    }
	    for (k = 0; k < subformat_count ; k++) {
		FFSTypeHandle sf = handle->subformats[k];
		FMFormat sfmf = fmf->subformats[k];
		sf->field_subformats = 
		    malloc(sfmf->field_count * sizeof(FFSTypeHandle));
		memset(sf->field_subformats, 0, 
		       sfmf->field_count * sizeof(FFSTypeHandle));
		for (i = 0; i < sfmf->field_count; i++) {
		    if (sfmf->field_subformats[i] != NULL) {
			int j;
			for (j = 0; j < subformat_count; j++) {
			    if (sfmf->field_subformats[i] == handle->subformats[j]->body) {
				sf->field_subformats[i] = handle->subformats[j];
			    }
			}
		    }
		}
	    }
	} else {
	    handle->subformats = NULL;
	    handle->field_subformats = NULL;
	}
    }
    return c->handle_list[index];
}

extern
FFSTypeHandle
FFSTypeHandle_from_encode(FFSContext c, char *buffer) 
{
    int index;
    FFSTypeHandle handle;
    /* first element in encoded buffer is format ID */
    FMFormat fmf = FMformat_from_ID(c->fmc, buffer);
    if (fmf == NULL) return NULL;
    index = fmf->format_index;
    handle = FFSTypeHandle_by_index(c, index);
    return handle;
}

extern
FFSContext
create_FFSContext_FM(FMContext fmc)
{
    FFSContext c;
    c = (FFSContext) malloc((size_t) sizeof(*c));
    init_float_formats();
    if (fmc == NULL) {
        fmc = create_FMcontext();
    } else {
        add_ref_FMcontext(fmc);
    }
    c->fmc = fmc;
    c->handle_list_size = 0;
    c->handle_list = NULL;
    c->tmp.tmp_buffer = NULL;
    c->tmp.tmp_buffer_size = 0;
    c->tmp.tmp_buffer_in_use_size = 0;

    return c;
}

extern FMContext
FMContext_from_FFS(FFSContext c)
{
    return c->fmc;
}

extern FMFormat
FMFormat_of_original(FFSTypeHandle h)
{
    return h->body;
}

extern void
free_FFSContext(FFSContext c)
{
    int i;
    free(c->tmp.tmp_buffer);
    for (i = 0; i < c->handle_list_size; i++) {
	if (c->handle_list[i]) free_FFSTypeHandle(c->handle_list[i]);
    }
    free(c->handle_list);
    free_FMcontext(c->fmc);
    free(c);
}

static void
reset_prior_conversions(FFSContext c)
{
    int i;
    for (i=0; i < c->handle_list_size; i++) {
	if (c->handle_list[i] != NULL) {
	    c->handle_list[i]->status = not_checked;
	}
    }
}

extern FFSTypeHandle
FFSset_fixed_target(FFSContext c, FMStructDescList struct_list)
{
    FMFormat fmf = register_data_format(c->fmc, struct_list);
    int index;
    FFSTypeHandle handle;
    index = fmf->format_index;
    handle = FFSTypeHandle_by_index(c, index);
    handle->is_fixed_target = 1;
    /* any new target may invalidate prior conversion decisions */
    reset_prior_conversions(c);
    return handle;
}

extern FFSTypeHandle
FFSset_simple_target(FFSContext c, char *format_name, FMFieldList field_list, size_t struct_size)
{
    FMStructDescRec struct_list[2];
    struct_list[0].format_name = format_name;
    struct_list[0].field_list = field_list;
    struct_list[0].struct_size = (int) struct_size;
    struct_list[0].opt_info = NULL;
    struct_list[1].format_name = NULL;
    struct_list[1].field_list = NULL;
    return FFSset_fixed_target(c, &struct_list[0]);
}

/* 
 * inputs:
 *  1.  maybe --> where data is now, base & variant
 *  2.  where we want data : base
 *  3.  should variant be buffered or go with base?
 *  5.  input_record_len
 *  
 * outputs:
 *  1.  where data should be read into: base & variant
 *              - or -
 *      if and to where data should be moved: base & variant
 *  2.  source address for conversion
 *  3.  dest address for conversion
 *  4.  src string address
 *  5.  final string address
 *  6.  base_string_offset
 */
typedef struct _conversion_action {
    char *cur_base;		/* where record is now (NULL if not in
				 * mem) */
    char *cur_variant;		/* where variant part is now (NULL if not
				 * in mem) */
    char *final_base;		/* where record should go (NULL if tmp
				 * buffer) */
    int variant_with_base;	/* true if variant should end up contig
				 * w/base */

/* outputs */
/* iovector read_iov[2]; int    read_icount; */
    void *src_address;
    void *dest_address;
    void *src_string_address;
    void *final_string_address;
    size_t base_string_offset;
} conversion_action, *conversion_action_ptr;

static int
in_place_base_conversion_possible(IOConversionPtr conv)
{
    switch (conv->conversion_type) {
    case buffer_and_convert:
	/* 
	 * Buffer_and_convert differs from Copy_Strings in that it 
	 * means no size-changing conversions need be done on the
	 * variant part of the record.
	 */
	return FALSE;		/* no in-place conversion */
    case copy_dynamic_portion:
	/* 
	 * Copy_Strings differs from Buffer_and_convert in that it 
	 * means size-changing conversions are required for the
	 * variant part of the record.  
	 */
	return FALSE;		/* no in-place conversion */
    case none_required:
    case direct_to_mem:
	return TRUE;
    default:
	assert(FALSE);
    }
    /* NOTREACHED */
    return FALSE;
}

static int
in_place_variant_conversion_possible(IOConversionPtr conv)
{
    switch (conv->conversion_type) {
    case copy_dynamic_portion:
	/* 
	 * Copy_Strings differs from Buffer_and_convert in that it 
	 * means size-changing conversions are required for the
	 * variant part of the record.  
	 */
	return FALSE;		/* no in-place conversion */
    case buffer_and_convert:
	/* 
	 * Buffer_and_convert differs from Copy_Strings in that it 
	 * means no size-changing conversions need be done on the
	 * variant part of the record.
	 */
	return TRUE;		/* yes, variant in-place conversion is
				 * possible */
    case none_required:
    case direct_to_mem:
	return TRUE;
    default:
	assert(FALSE);
    }
    /* NOTREACHED */
    return FALSE;
}

extern int
decode_in_place_possible(FFSTypeHandle format)
{
    return FFSdecode_in_place_possible(format);
}

extern int
FFSdecode_in_place_possible(FFSTypeHandle format)
{
    if (format->conversion != NULL) {
	if (format->body->variant) {
	    return in_place_base_conversion_possible(format->conversion) &&
		in_place_variant_conversion_possible(format->conversion);
	} else {
	    return in_place_base_conversion_possible(format->conversion);
	}
    } else {
	return 0;
    }
}

#define expand_size_to_align(size) ((((size) & 0x7) == 0) ? (size) : (((size) + 8) & (size_t) -8))

static int
set_conversion_params(FFSTypeHandle ioformat, int64_t input_record_len, IOConversionPtr conv, conversion_action_ptr params)
{
    FFSContext c = ioformat->context;
    size_t final_base_size;
    size_t src_base_size;
    int64_t possible_converted_variant_size;
    int64_t orig_variant_size;

    ssize_t dest_offset;
    void *dest_address;
    ssize_t src_offset;
    void *src_address;
    ssize_t src_string_offset;
    void *src_string_address;
    ssize_t final_string_offset;
    void *final_string_address;

    final_base_size = (size_t) expand_size_to_align((size_t)(ioformat->body->record_length + conv->base_size_delta));
    src_base_size = expand_size_to_align(ioformat->body->record_length);
    possible_converted_variant_size =
	final_variant_size_for_record(input_record_len, conv);
    orig_variant_size = input_record_len - expand_size_to_align(ioformat->body->record_length);

    make_tmp_buffer(&c->tmp, 0);
    /* set base dest values */
    if (params->final_base == NULL) {
	/* need memory for at least the base record in temp area */
	int64_t buffer_required = Max(final_base_size, src_base_size);
	dest_offset = add_to_tmp_buffer(&c->tmp, buffer_required);
	dest_address = NULL;
	if (dest_offset == -1) return 0;
    } else {
	dest_address = params->final_base;
	dest_offset = 0;
    }

    /* set base src values */
    if (in_place_base_conversion_possible(conv)) {
	/* 
	 * we can convert in place, so just use the already determined 
	 * destination for the src values
	 */
	src_offset = dest_offset;
	src_address = dest_address;
    } else {
	/* no converting in place, so source must be different than base */
	if ((params->cur_base == NULL) ||
	    (params->cur_base == (char *) dest_address)) {
	    /* 
	     * either the source is not in memory or it is the same as 
	     * where we want the record to end up.  Need temporary space.
	     */
	    int64_t source_base_size = expand_size_to_align(ioformat->body->record_length);
	    src_offset = add_to_tmp_buffer(&c->tmp, source_base_size);
	    src_address = NULL;
	    if (src_offset == -1) return 0;
	} else {
	    src_offset = 0;
	    src_address = params->cur_base;
	}
    }

    /* set final string address */
    if ((params->final_base == NULL) || (!params->variant_with_base)) {
	/* 
	 * either they didn't specified a specific base address, or they said
	 * to put the variant in temporary memory.  Either way, we need
	 * temporary memory for the final variant part.
	 *
	 * NOTE: the record length may have provided padding essential to the
	 * alignment of variant arrays in the string area.  Get a little extra
	 * memory so we can modify the base if necessary.
	 */
	int align_pad = (8 - ioformat->body->record_length) & 0x7;
	int64_t buffer_required = Max(possible_converted_variant_size + align_pad,
				  orig_variant_size + align_pad);
	buffer_required = expand_size_to_align(buffer_required);
	final_string_offset = add_to_tmp_buffer(&c->tmp, buffer_required);
	final_string_address = NULL;
	if (final_string_offset == -1) return 0;
    } else {
	final_string_offset = 0;
	final_string_address =
	    params->final_base + (size_t)expand_size_to_align((size_t)(ioformat->body->record_length +
								       conv->base_size_delta));
    }

    /* set variant src values */
    if (in_place_variant_conversion_possible(conv)) {
	/* 
	 * we can convert in place, so just use the already determined 
	 * destination for the src values
	 */
	src_string_offset = final_string_offset;
	src_string_address = final_string_address;
    } else {
	/* no converting in place, so source must be different than base */
	if ((params->cur_variant == NULL) ||
	    (params->cur_base == (char *) dest_address)) {
	    /* 
	     * either the source is not in memory or it is the same as 
	     * where we want the record to end up.  Need temporary space.
	     */
	    int64_t source_variant_size =	/* plus possible alignment of 8 */
		input_record_len - ioformat->body->record_length + 8;
	    src_string_offset = add_to_tmp_buffer(&c->tmp, source_variant_size);
	    src_string_address = NULL;
	    if (src_string_offset == -1) return 0;
	} else {
	    src_string_offset = 0;
	    src_string_address = params->cur_variant;
	}
    }
    if (dest_address == NULL) {
	params->dest_address = (char *) c->tmp.tmp_buffer + dest_offset;
    } else {
	params->dest_address = dest_address;
    }
    if (src_address == NULL) {
	params->src_address = (char *) c->tmp.tmp_buffer + src_offset;
    } else {
	params->src_address = src_address;
    }
    if (final_string_address == NULL) {
	params->final_string_address = (char *) c->tmp.tmp_buffer + final_string_offset;
    } else {
	params->final_string_address = final_string_address;
    }
    if (src_string_address == NULL) {
	params->src_string_address = (char *) c->tmp.tmp_buffer + src_string_offset;
    } else {
	params->src_string_address = src_string_address;
    }
    return 1;
}

static int64_t
final_variant_size_for_record(int64_t input_record_len, IOConversionPtr conv)
{
    return (int) ((input_record_len - conv->ioformat->body->record_length)
		  * conv->max_var_expansion);
}

extern size_t
FFS_decode_length_format(FFSContext context, FFSTypeHandle ioformat, 
			 size_t record_length)
{
    IOConversionPtr conv;
    int64_t variant_part, final_base_size, src_base_size;

    if (ioformat == NULL)
	return (size_t)-1;
    conv = ioformat->conversion;
    if (!ioformat->conversion)
	return record_length;
    variant_part = final_variant_size_for_record(record_length,
						 ioformat->conversion);
    final_base_size = expand_size_to_align((int64_t)(ioformat->body->record_length +
						     conv->base_size_delta));
    src_base_size = expand_size_to_align(ioformat->body->record_length);
    return variant_part + Max(final_base_size, src_base_size);
}

extern size_t
FFS_est_decode_length(FFSContext context, char *src, size_t record_length)
{
    FFSTypeHandle ioformat = FFSTypeHandle_from_encode(context, src);
    return FFS_decode_length_format(context, ioformat, record_length);
}

extern int
FFSheader_size(FFSTypeHandle ioformat)
{
    int header_size = ioformat->body->server_ID.length;
    int align_pad;
    if (ioformat->body->variant) {
	if (ioformat->body->IOversion <= 3) {
	    header_size += 4;
	} else {
	    header_size += 8;  // length info
	}
    }
    align_pad = (8 - header_size) & 0x7;
    return header_size + align_pad;
}

extern
int
FFShas_conversion(FFSTypeHandle ioformat)
{
    return (ioformat->conversion != NULL);
}

static int
FFSinternal_decode(FFSTypeHandle ioformat, char *src, void *dest, int to_buffer)
{
    FFSContext iofile = ioformat->context;
    IOConversionPtr conv;
    int64_t input_record_len;
    IOConversionStruct null_conv;
    int align_pad;
    int header_size = FFSheader_size(ioformat);
    int data_align_pad;
    conversion_action params;

    conv = ioformat->conversion;

    if (conv == NULL) {
	null_conv.conversion_type = none_required;
	null_conv.max_var_expansion = 1.0;
	null_conv.conv_count = 0;
	null_conv.base_size_delta = 0;
	null_conv.context = iofile;
	null_conv.ioformat = ioformat;
	null_conv.native_field_list = NULL;
	conv = &null_conv;
    }
    if (iofile != conv->context) {
	fprintf(stderr, "IOFile and conversion mismatch\n");
	return -1;
    }
    if (ioformat == NULL)
	return 0;

    if (ioformat->body->variant) {
	if (ioformat->body->IOversion <= 3) {
	    FILE_INT record_len;
	    int len_align_pad = (4 - ioformat->body->server_ID.length) & 3;
	    FILE_INT *len_ptr = (FILE_INT *) (src + ioformat->body->server_ID.length +
					      len_align_pad);
	    memcpy(&record_len, len_ptr, 4);
	    if (ioformat->body->byte_reversal)
		byte_swap((char *) &record_len, 4);
	    input_record_len = record_len;
	} else {
	    uint64_t record_len;
	    uint64_t *len_ptr = (uint64_t *) (src + ioformat->body->server_ID.length);
	    memcpy(&record_len, len_ptr, 8);
	    if (ioformat->body->byte_reversal)
		byte_swap((char *) &record_len, 8);
	    input_record_len = record_len;
	}
    } else {
	input_record_len = ioformat->body->record_length;
    }
    align_pad = (8 - header_size) & 0x7;
    header_size += align_pad;
    data_align_pad = (8 - ioformat->body->record_length) & 0x7;

    params.cur_base = src + header_size;
    params.cur_variant = params.cur_base + ioformat->body->record_length +
	data_align_pad;
    params.final_base = dest;
    params.variant_with_base = to_buffer;
    if (!set_conversion_params(ioformat, input_record_len, conv, &params))
	return 0;

    /* update stats 
    iofile->stats.decode_bytes += input_record_len;
    iofile->stats.decode_msg_count++;*/

    if (params.src_address != params.cur_base) {
	memcpy(params.src_address, params.cur_base, ioformat->body->record_length);
    }
    if (params.src_string_address != params.cur_variant) {
	if (input_record_len - ioformat->body->record_length - data_align_pad > 0) {
	    memcpy(params.src_string_address, params.cur_variant,
		   input_record_len - ioformat->body->record_length - data_align_pad);
	}
    }
    if (conv->conversion_type != none_required) {
	FFSconvert_record(conv, params.src_address, params.dest_address,
			 params.final_string_address,
			 params.src_string_address);
    }
    return 1;
}

static int
check_conversion(FFSTypeHandle ioformat)
{
    if (ioformat->conversion == NULL) {
	if (ioformat->status == not_checked) {
	    FFS_determine_conversion(ioformat->context, ioformat);
	}
	if (ioformat->status == none_available) {
	    fprintf(stderr, "FFS Warning:  Attempting to decode when no conversion has been set.  \n  Record is of type \"%s\", ioformat 0x%p.\n  No data returned.\n",
		    ioformat->body->format_name, ioformat);
	    ioformat->status = none_available;
	    return 0;
	}
    }
    return 1;
}

extern int
FFSdecode(FFSContext iocontext, char *src, char *dest)
{
    FFSTypeHandle ioformat;
    ioformat = FFSTypeHandle_from_encode(iocontext, src);
    if (ioformat == NULL) {
	return 0;
    }
    if (!check_conversion(ioformat)) {
	return 0;
    }
    return FFSinternal_decode(ioformat, src, dest, 0);
}

extern int
FFSdecode_in_place(FFSContext iocontext, char *src, void **dest_ptr)
{
    FFSTypeHandle ioformat = FFSTypeHandle_from_encode(iocontext, src);
    int header_size;
    int ret;

    if (ioformat == NULL) {
	return 0;
    }
    if (!check_conversion(ioformat)) {
	*dest_ptr = NULL;
	return 0;
    }
    header_size = FFSheader_size(ioformat);
    ret = FFSinternal_decode(ioformat, src, src + header_size, 1);
    *dest_ptr = src + header_size;
    return ret;
}

extern int
FFSdecode_to_buffer(FFSContext iocontext, char *src, void *dest)
{
    FFSTypeHandle ioformat;
    ioformat = FFSTypeHandle_from_encode(iocontext, src);
    if (ioformat == NULL) {
	return 0;
    }
    if (!check_conversion(ioformat)) {
	return 0;
    }
    return FFSinternal_decode(ioformat, src, dest, 1);
}

extern FFSTypeHandle
FFS_target_from_encode(FFSContext c, char *data)
{
    FFSTypeHandle f;
    /* first element in encoded buffer is format ID */
    f = FFSTypeHandle_from_encode(c, data);
    if (f == NULL) return NULL;
    if (f->status == not_checked) {
	FFS_determine_conversion(c, f);
    }
    if (f->status == none_available) return NULL;
    if (f->status == conversion_set) return f->conversion_target;
    assert(0);
    return NULL;
}

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

char *
FFSTypeHandle_name(FFSTypeHandle f)
{
    return f->body->format_name;
}

#define TMP_BUFFER_INIT_SIZE 1024

FFSBuffer
create_FFSBuffer()
{
    FFSBuffer buf = malloc(sizeof(struct _FFSBuffer));
    buf->tmp_buffer = NULL;
    buf->tmp_buffer_size = 0;
    buf->tmp_buffer_in_use_size = 0;
    return buf;
}

FFSBuffer
create_fixed_FFSBuffer(char *buffer, size_t size)
{
    FFSBuffer buf = malloc(sizeof(struct _FFSBuffer));
    buf->tmp_buffer = buffer;
    buf->tmp_buffer_size = -(ssize_t)size;
    buf->tmp_buffer_in_use_size = 0;
    return buf;
}

void
free_FFSBuffer(FFSBuffer buf)
{
    if ((buf->tmp_buffer_size > 0) && buf->tmp_buffer)
	free(buf->tmp_buffer);
    free(buf);
}

extern
char *
make_tmp_buffer(FFSBuffer buf, int64_t size)
{
    if (buf->tmp_buffer_size < 0) {
	/* fixed size buffer */
	if (size > (-buf->tmp_buffer_size)) return NULL;
	return buf->tmp_buffer;
    }
    if (buf->tmp_buffer_size == 0) {
	int64_t tmp_size = Max(size, TMP_BUFFER_INIT_SIZE);
	buf->tmp_buffer = malloc(tmp_size);
	buf->tmp_buffer_size = tmp_size;
    }
    if (size > buf->tmp_buffer_size) {
	buf->tmp_buffer = realloc(buf->tmp_buffer, size);
	if (buf->tmp_buffer) {
	    buf->tmp_buffer_size = size;
	} else {
	    buf->tmp_buffer_size = 0;
	}
    }
    buf->tmp_buffer_in_use_size = size;
    return buf->tmp_buffer;
}

static
ssize_t
add_to_tmp_buffer(FFSBuffer buf, size_t size)
{
    ssize_t old_size = buf->tmp_buffer_in_use_size;
    size += old_size;

    if (buf->tmp_buffer_size < 0) {
	/* fixed size buffer */
	if (size > (size_t)(-buf->tmp_buffer_size)) return -1;
    } else {
	if (buf->tmp_buffer_size == 0) {
	    int64_t tmp_size = Max(size, TMP_BUFFER_INIT_SIZE);
	    buf->tmp_buffer = malloc(tmp_size);
	}
	if (size > (size_t)buf->tmp_buffer_size) {
	    buf->tmp_buffer = realloc(buf->tmp_buffer, size);
	    buf->tmp_buffer_size = size;
	}
    }
    if (!buf->tmp_buffer) {
	buf->tmp_buffer_size = 0;
	return -1;
    } else {
	buf->tmp_buffer_in_use_size = size;
    }

    return old_size;
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

static size_t
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
quick_get_pointer(FMFieldPtr iofield, void *data)
{
    union {
	void *p;
	unsigned long tmp;
	int tmpi[2];
    } u;
    data = (void *) ((char *) data + iofield->offset);
    /* only used when field type is an integer and aligned by its size */
    switch (iofield->size) {
    case 1:
	u.tmp = (unsigned long) (*((unsigned char *) data));
	break;
    case 2:
	u.tmp = (unsigned long) (*((unsigned short *) data));
	break;
    case 4:
    {
	unsigned int tmpi;
	memcpy(&tmpi, data, 4);
	u.tmp = (unsigned long) tmpi;
	break;
    }
    case 8:
#if SIZEOF_SIZE_T == 8
	memcpy(&u.p, data, 8);
#else
	{
	  int tmp;
	  /* must be fetching 4 bytes of the 8 available */
	  if (WORDS_BIGENDIAN)
	    memcpy(&tmp, ((char*)data) + 4, 4);
	  else
	    memcpy(&tmp, data, 4);
	  u.tmp = (unsigned long) tmp;
	}
#endif
	break;
    }
    return u.p;
}

void
quick_put_ulong(FMFieldPtr iofield, size_t value, void *data)
{
    data = (void *) ((char *) data + iofield->offset);
    /* only used when field type is an integer and aligned by its size */
    switch (iofield->size) {
    case 1:
	*((unsigned char *) data) = (unsigned char) value;
	break;
    case 2:
	*((unsigned short *) data) = (unsigned short) value;
	break;
    case 4:{
	unsigned int tmp = (unsigned int) value;
	memcpy(data, &tmp, 4);
	break;
	}
    case 8:{
	size_t tmp = value;
	memcpy(data, &tmp, 8);
	break;
	}
    }
}

