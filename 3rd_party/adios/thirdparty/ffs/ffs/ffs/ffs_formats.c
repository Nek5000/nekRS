#include "config.h"

#include "assert.h"
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
#include <sys/types.h>
#include <ctype.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#ifdef DO_DCG
#include "dill.h"
#else 
#define dill_stream void*
#define dill_create_stream() NULL
/*stuff*/
enum {
    DILL_C,    /* char */
    DILL_UC,   /* unsigned char */
    DILL_S,    /* short */
    DILL_US,   /* unsigned short */
    DILL_I,    /* int */
    DILL_U,    /* unsigned */
    DILL_L,    /* long */
    DILL_UL,   /* unsigned long */
    DILL_P,    /* pointer */
    DILL_F,    /* floating */
    DILL_D,    /* double */
    DILL_V,    /* void */
    DILL_B,    /* block structure */
    DILL_ERR   /* no type */
};
static int
dill_type_align(dill_stream s, int t)
{
  switch (t) {
  case DILL_C: case DILL_UC: return 1;
  case DILL_S: case DILL_US: return 2;
  case DILL_I: case DILL_U: return 4;
  case DILL_L: case DILL_UL: return sizeof(long);
  case DILL_P: return sizeof(char*);
  case DILL_F: return sizeof(float);
  case DILL_D: return sizeof(double);
  default:
    return 0;
  }
}
#endif
#include "ffs.h"
#include "fm_internal.h"
#include "assert.h"

static int
local_size(int field_index, FMFormat f, FMTypeDesc *type, 
	   FMStructDescList subformats, int immediate) 
{
    switch(type->type) {
    case FMType_string:
	return sizeof(char*);
    case FMType_pointer:
	/* pointer is sizeof(char*), but we want to return size of item */
	if (immediate) return (sizeof(char*));
	return local_size(field_index, f, type->next, subformats, 0);
    case FMType_array:
	/* here, we always return the size of the element */
	return local_size(field_index, f, type->next, subformats, immediate);
    case FMType_subformat: {
	char *subformat_name = f->field_subformats[field_index]->format_name;
	while (subformats && (subformats->format_name != NULL)) {
	    if (strcmp(subformat_name, subformats->format_name) == 0) {
		return subformats->struct_size;
	    }
	    subformats++;
	}
	assert(0);
    }
    case FMType_simple: {
	int given_size = f->field_list[field_index].field_size;
	switch (type->data_type) {
	case integer_type:
	case unsigned_type:
	case boolean_type:
	case enumeration_type:
	    if ((given_size == 0) || (given_size > sizeof(long))) {
		given_size = sizeof(long);
	    }
	    break;
	case float_type:
	    if ((given_size == 0) || (given_size > SIZEOF_DOUBLE)) {
		given_size = sizeof(double);
	    }
	    break;
	case char_type:
	    if (given_size == 0) {
		given_size = 1;
	    }
	    break;
	case string_type:
	    given_size = sizeof(char*);
	    break;
	default:
	    assert(0);
	}
	return given_size;
    }
    }
    return 0;
}

static int
min_align_size(int size)
{
    int align_size = 8;		/* conservative on current machines */
    switch (size) {
    case 7:
    case 6:
    case 5:
    case 4:
	align_size = 4;
	break;
    case 3:
	align_size = 2;
	break;
    case 2:
    case 1:
	align_size = size;
	break;
    }
    return align_size;
}

static int
min_align_type(FMdata_type typ, int size)
{
#ifndef DO_DCG
    return min_align_size(size);
#else
    static dill_stream s = NULL;
    if (s == NULL) s = dill_create_raw_stream();
    switch (typ) {
    case float_type:
	if (size == dill_type_size(s, DILL_D)) return dill_type_align(s, DILL_D);
	if (size == dill_type_size(s, DILL_F)) return dill_type_align(s, DILL_F);
	/* punt */
	return min_align_size(size);
    case integer_type: case char_type: case string_type:
	if (size == dill_type_size(s, DILL_C)) return dill_type_align(s, DILL_C);
	if (size == dill_type_size(s, DILL_S)) return dill_type_align(s, DILL_S);
	if (size == dill_type_size(s, DILL_I)) return dill_type_align(s, DILL_I);
	if (size == dill_type_size(s, DILL_L)) return dill_type_align(s, DILL_L);
	/* punt */
	return min_align_size(size);
    case unsigned_type: case enumeration_type: case boolean_type:
	if (size == dill_type_size(s, DILL_UC)) return dill_type_align(s, DILL_UC);
	if (size == dill_type_size(s, DILL_US)) return dill_type_align(s, DILL_US);
	if (size == dill_type_size(s, DILL_U)) return dill_type_align(s, DILL_U);
	if (size == dill_type_size(s, DILL_UL)) return dill_type_align(s, DILL_UL);
	/* punt */
	return dill_type_align(s, DILL_B);
    default:
	return dill_type_align(s, DILL_B);
    }
#endif
}

static int
align_field(int local_size, int cur_offset, FMTypeDesc *type, dill_stream s, 
	    int *max_align_p)
{
    int new_offset, align_req;
    switch(type->type) {
    case FMType_array:
	if (type->static_size == 0) {
	    /* variant array */
	    if (s!= NULL) {
		align_req = dill_type_align(s, DILL_P);
	    } else {
		align_req = sizeof(char*);
	    }
	} else {
	    /* alignment of the elements */
	    return align_field(local_size, cur_offset, type->next, s, 
			       max_align_p);
	}
	break;
    case FMType_string:
    case FMType_pointer:
	if (s!= NULL) {
	    align_req = dill_type_align(s, DILL_P);
	} else {
	    align_req = sizeof(char*);
	}
	break;    
    case FMType_subformat:
        /* realistically, should be max align of subfields, too hard */
        if (s != NULL) {
	    align_req = 0;
	    int i;
	    for (i = DILL_C; i <DILL_V; i++) {
	        if (align_req < dill_type_align(s, i)) {
		    align_req = dill_type_align(s, i);
		}
	    }
	} else {
	    align_req = sizeof(long);   /* just a guess */
	}
	break;
    case FMType_simple: 
        align_req = min_align_type(type->data_type, local_size);
	break;
    default:
        align_req = -1;
	break;
    }
    if (align_req > *max_align_p) {
	*max_align_p = align_req;
    }
    new_offset = (cur_offset + align_req - 1) & (-align_req);
    return new_offset;
}

static void
generate_localized_subformat(FMFormat f, FMStructDescList subformats, 
			     dill_stream s)
{
    FMFieldList fl = copy_field_list(f->field_list);
    int i = 0;
    int last_field_end = 0;
    int max_align = 0;
    while (fl[i].field_name != NULL) {
	int align_size = local_size(i, f, &f->var_list[i].type_desc, 
				    subformats+1, 1);

	fl[i].field_size = local_size(i, f, &f->var_list[i].type_desc, 
				      subformats+1, 0);
	if (fl[i].field_offset < last_field_end) {
	    fl[i].field_offset = last_field_end;
	}
	fl[i].field_offset = align_field(align_size, fl[i].field_offset, 
					 &f->var_list[i].type_desc, s, &max_align);
	last_field_end = fl[i].field_offset + align_size;
	i++;
    }
    subformats->format_name = strdup(f->format_name);
    subformats->field_list = fl;
    subformats->struct_size = (last_field_end + max_align - 1) & (-max_align);
    subformats->opt_info = NULL;
}

extern FMStructDescList
get_localized_formats_dill(FMFormat f, dill_stream s)
{
    FMStructDescList ret_list;
    int subformat_count = 0;
    while(f->subformats && (f->subformats[subformat_count] != NULL)) {
	subformat_count++;
    }
    ret_list = malloc(sizeof(ret_list[0]) * (subformat_count+2));
    int i = 0;
    ret_list[subformat_count+1].format_name = NULL;
    ret_list[subformat_count+1].field_list = NULL;
    ret_list[subformat_count+1].struct_size = 0;
    ret_list[subformat_count+1].opt_info = 0;
    for (i=0; i < subformat_count; i++) {
	generate_localized_subformat(f->subformats[subformat_count-i-1], &ret_list[subformat_count - i], s);
    }
    generate_localized_subformat(f, ret_list, s);
    return ret_list;
}

extern FMStructDescList
get_localized_formats(FMFormat f)
{
    static dill_stream s = NULL;
    if (s == NULL) {
	s = dill_create_stream();
    }
    return get_localized_formats_dill(f, s);
}

