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
#include <stdint.h>

#ifdef HAVE_DILL_H
#include "dill.h"
#define static_ctx c 
#define VCALL7V(subr, argstr, arg1, arg2, arg3, arg4, arg5, arg6, arg7) dill_scallv(c, (void*)subr, argstr, arg1, arg2, arg3, arg4, arg5, arg6, arg7)
#define VCALL6V(subr, argstr, arg1, arg2, arg3, arg4, arg5, arg6) dill_scallv(c, (void*)subr, argstr, arg1, arg2, arg3, arg4, arg5, arg6)
#define VCALL5V(subr, argstr, arg1, arg2, arg3, arg4, arg5) dill_scallv(c, (void*)subr, argstr, arg1, arg2, arg3, arg4, arg5)
#define VCALL4V(subr, argstr, arg1, arg2, arg3, arg4) dill_scallv(c, (void*)subr, argstr, arg1, arg2, arg3, arg4)
#define VCALL4P(subr, argstr, arg1, arg2, arg3, arg4) dill_scallp(c, (void*)subr, argstr, arg1, arg2, arg3, arg4)
#define VCALL3P(subr, argstr, arg1, arg2, arg3, arg4) dill_scallp(c, (void*)subr, argstr, arg1, arg2, arg3, arg4)
#define VCALL3V(subr, argstr, arg1, arg2, arg3) dill_scallv(c, (void*)subr, argstr, arg1, arg2, arg3)
#define VCALL2V(subr, argstr, arg1, arg2) dill_scallv(c, (void*)subr, argstr, arg1, arg2)
#define TYPE_ALIGN(c, t) dill_type_align(c, t)
#define _vrr(x) x
#endif
#include "ffs.h"
#include "fm_internal.h"
#include "ffs_internal.h"
#if defined(HAVE_DILL_H)
#include "ffs_gen.h"
#endif
#include "assert.h"

#ifdef __has_feature
#if __has_feature(thread_sanitizer)
#define NO_SANITIZE_THREAD __attribute__((no_sanitize("thread")))
#endif
#endif

#ifndef NO_SANITIZE_THREAD
#define NO_SANITIZE_THREAD
#endif

static MAX_INTEGER_TYPE get_big_int(FMFieldPtr iofield, void *data);
static MAX_FLOAT_TYPE get_big_float(FMFieldPtr iofield, void *data);
static MAX_UNSIGNED_TYPE get_big_unsigned(FMFieldPtr iofield, void *data);
static void byte_swap(char *data, int size);

static int get_double_warn = 0;
static int get_long_warn = 0;

FMfloat_format ffs_my_float_format = Format_Unknown;
/* 
 * ffs_reverse_float_formats identifies for each format what, 
 * if any, format is its byte-swapped reverse.
*/
FMfloat_format ffs_reverse_float_formats[] = {
    Format_Unknown, /* no format complements unknown */
    Format_IEEE_754_littleendian, /* littleendian complements bigendian */
    Format_IEEE_754_bigendian, /* bigendian complements littleendian */
    Format_Unknown /* no exact opposite for mixed-endian (ARM) */
};

static unsigned char IEEE_754_4_bigendian[] = 
  {0x3c, 0x00, 0x00, 0x00};
static unsigned char IEEE_754_4_littleendian[] = 
  {0x00, 0x00, 0x00, 0x3c};
static unsigned char IEEE_754_8_bigendian[] = 
  {0x3f, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
static unsigned char IEEE_754_8_littleendian[] = 
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f};
static unsigned char IEEE_754_8_mixedendian[] = 
  {0x00, 0x00, 0x80, 0x3f, 0x00, 0x00, 0x00, 0x00};
static unsigned char IEEE_754_16_bigendian[] = 
  {0x3f, 0xf8, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
static unsigned char IEEE_754_16_littleendian[] = 
  {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf8, 0x3f};
static unsigned char IEEE_754_16_mixedendian[] = 
  {0x00, 0x00, 0xf8, 0x3f, 0x00, 0x00, 0x00, 0x00, 
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

static char *float_format_str[] = {
    "Unknown float format",
    "IEEE 754 float (bigendian)", 
    "IEEE 754 float (littleendian)",
    "IEEE 754 float (mixedendian)"};

static int IO_shut_up = 0;

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

static FMfloat_format
FFSinfer_float_format(char *float_magic, int object_len)
{
    switch (object_len) {
    case 4:
	if (memcmp(float_magic, &IEEE_754_4_bigendian[0], 4) == 0) {
	    return Format_IEEE_754_bigendian;
	} else if (memcmp(float_magic, &IEEE_754_4_littleendian[0], 4) == 0) {
	    return Format_IEEE_754_littleendian;
	}
	break;
    case 8:
	if (memcmp(float_magic, &IEEE_754_8_bigendian[0], 8) == 0) {
	    return Format_IEEE_754_bigendian;
	} else if (memcmp(float_magic, &IEEE_754_8_littleendian[0], 8) == 0) {
	    return Format_IEEE_754_littleendian;
	} else if (memcmp(float_magic, &IEEE_754_8_mixedendian[0], 8) == 0) {
	    return Format_IEEE_754_mixedendian;
	}
	break;
    case 16:
	if (memcmp(float_magic, &IEEE_754_16_bigendian[0], 16) == 0) {
	    return Format_IEEE_754_bigendian;
	} else if (memcmp(float_magic, &IEEE_754_16_littleendian[0], 16) ==0){
	    return Format_IEEE_754_littleendian;
	} else if (memcmp(float_magic, &IEEE_754_16_mixedendian[0], 16) == 0){
	    return Format_IEEE_754_mixedendian;
	}
	break;
    }
    return Format_Unknown;
}

extern void
init_float_formats()
{
    static int done = 0;
    if (!done) {
	double d = MAGIC_FLOAT;
	ffs_my_float_format = FFSinfer_float_format((char*)&d, sizeof(d));
	switch (ffs_my_float_format) {
	case Format_IEEE_754_bigendian:
	case Format_IEEE_754_littleendian:
	case Format_IEEE_754_mixedendian:
	    break;
	case Format_Unknown:
	    fprintf(stderr, "Warning, unknown local floating point format\n");
	    break;
	}
	done++;
    }
}
	

void
FFSfree_conversion(IOConversionPtr conv)
{
    int i;
    for (i = 0; i < conv->conv_count; i++) {
	if ((conv->conversions[i].subconversion) && 
	    (conv->conversions[i].subconversion != conv)) {
	    FFSfree_conversion(conv->conversions[i].subconversion);
	}
	if (conv->conversions[i].default_value) {
	    free(conv->conversions[i].default_value);
	}
    }
    if (conv->native_field_list) {
	for (i = 0; conv->native_field_list[i].field_name != NULL; i++) {
	    free((char*)conv->native_field_list[i].field_name);
	    free((char*)conv->native_field_list[i].field_type);
	}
	free(conv->native_field_list);
    }
    if (conv->free_func != NULL) {
	conv->free_func(conv->free_data);
    }
    free(conv);
}

/*
 * Copy fieldname from io_field to out_field_name and 
 * removes the default value from it.
 * default value is then stored in default_val after converting it to 
 * appropriate format.
 * */

static void
field_name_strip_get_default(const FMField *io_field, char *out_field_name, void **default_val)
{
    char *s, *s1;
    char *base_type = base_data_type(io_field->field_type);
    FMdata_type data_type = FMstr_to_data_type(base_type);
    strncpy(out_field_name, io_field->field_name, 128);
    s = strchr(out_field_name, '(');
    *default_val = NULL;
    free(base_type);
    if(s) {
	*s = '\0';
	s++;
	s1 = strchr(s, ')');
	if(s1) {
	    *s1 = '\0';
	}
	if((int)(s1-s)> 0){
	    str_to_val(s, data_type, io_field->field_size, default_val);
	}
    }
}
 
static void
create_default_conversion(FMField iofield, void *default_value, 
			  IOConversionPtr *conv_ptr_ptr, int conv_index)
{
    IOConversionPtr	conv_ptr =
	(IOConversionPtr) realloc(*conv_ptr_ptr, sizeof(IOConversionStruct) +
				  conv_index * sizeof(IOconvFieldStruct));
    *conv_ptr_ptr = conv_ptr;
    memset(&conv_ptr->conversions[conv_index].src_field, 0, 
	   sizeof(conv_ptr->conversions[conv_index].src_field));
    conv_ptr->conversions[conv_index].subconversion = NULL;
    conv_ptr->conversions[conv_index].iovar = NULL;
    conv_ptr->conversions[conv_index].dest_size = iofield.field_size;
    conv_ptr->conversions[conv_index].dest_offset = iofield.field_offset;
    conv_ptr->conversions[conv_index].default_value = default_value;
    conv_ptr->conversions[conv_index].rc_swap = no_row_column_swap;
}

/*
 * differs from find_field in fm in that it includes cur_field in search shortcut
 */
static long
find_field_for_conv(char *field_name, FMFieldList fields, int cur_field, void *search_help)
{
    int i;
    if (cur_field > 10) {
      /* search close first */
      for (i = cur_field; i > cur_field - 10; i--) {
	if (strcmp(field_name, fields[i].field_name) == 0) {
	    return i;
	}
      }
      for (i = cur_field + 1; i < cur_field + 10; i++) {
	if (fields[i].field_name == NULL) break;
	if (strcmp(field_name, fields[i].field_name) == 0) {
	    return i;
	}
      }
    }
    i = 0;
    while (fields[i].field_name != NULL) {
	if (strcmp(field_name, fields[i].field_name) == 0) {
	    return i;
	}
	i++;
    }
    return -1;
}

static
IOConversionPtr
create_conversion(FFSTypeHandle src_ioformat, FMFieldList target_field_list, int target_struct_size, int pointer_size,
		  int byte_reversal, FMfloat_format target_fp_format, IOconversion_type initial_conversion,
		  int target_column_major, int string_offset_size, int converted_strings, FMStructDescList target_list)
{
    int target_field_count = count_FMfield(target_field_list);
    FMFieldList nfl_sort = copy_field_list(target_field_list);
    FMFieldList input_field_list = src_ioformat->body->field_list;
    FMVarInfoList input_var_list = src_ioformat->body->var_list;
    IOconversion_type conv = initial_conversion;
    int input_index, conv_index = 0, i = 0;
    FMfloat_format src_float_format = src_ioformat->body->float_format;
    IOConversionPtr conv_ptr =
	(IOConversionPtr) malloc(sizeof(IOConversionStruct) +
				 target_field_count * sizeof(IOconvFieldStruct));
    int column_row_swap_necessary = (target_column_major != src_ioformat->body->column_major_arrays);
    
    if (target_fp_format == Format_Unknown) target_fp_format = ffs_my_float_format;

    conv_ptr->notify_of_format_change = 0;
    conv_ptr->context = src_ioformat->context;
    conv_ptr->ioformat = src_ioformat;
    conv_ptr->base_size_delta = target_struct_size -
	src_ioformat->body->record_length;
    conv_ptr->max_var_expansion = 1.0;
    conv_ptr->conv_count = 0;
    qsort(nfl_sort, target_field_count, sizeof(nfl_sort[0]),
	  field_offset_compar);
    conv_ptr->native_field_list = nfl_sort;
    conv_ptr->target_pointer_size = pointer_size;
    conv_ptr->required_alignment = 8; /* placeholder */
    conv_ptr->free_data = NULL;
    conv_ptr->free_func = NULL;
    conv_ptr->conv_func = NULL;
    conv_ptr->conv_func4 = NULL;
    conv_ptr->conv_func2 = NULL;
    conv_ptr->conv_func1 = NULL;
    conv_ptr->conv_pkg = NULL;
    conv_ptr->string_offset_size = string_offset_size;
    conv_ptr->converted_strings = converted_strings;

    /* 
     * We assume that the fields listed in the target_field_list are
     * those that the user is interested in.  Skipping input fields
     * is OK... 
     */

    if (src_ioformat->body->record_length > target_struct_size) {
	/* if input record is larger than target record, must do buffered */
	switch (conv) {
	case none_required:
	case direct_to_mem:
	    conv = buffer_and_convert;
	    break;
	case buffer_and_convert:
	case copy_dynamic_portion:
	    break;
	default:
	    assert(FALSE);
	}
    }
    /* try for no conversion/direct_to_memory first.. */
  restart:
    if (conv_index > 0) {
	/* really a restart, NULL a couple of things */
	int i = 0;
	for (i=0; i < conv_index; i++) {
	    conv_ptr->conversions[i].iovar = NULL;
	    if (conv_ptr->conversions[i].default_value) {
		free(conv_ptr->conversions[i].default_value);
	    }
	    if (conv_ptr->conversions[i].subconversion) {
		FFSfree_conversion(conv_ptr->conversions[i].subconversion);
	    }
	}
    }
    conv_ptr->conversion_type = conv;
    input_index = conv_index = 0;
    for (i = 0; i < target_field_count; i++) {
	FMField input_field;
	FMdata_type in_data_type, target_data_type;
	long in_elements, target_elements;
	void *default_val = NULL;
	char *tmp_field_name = NULL;
	char *search_name;
	int multi_dimen_array = 0;

	/* 
	 * all fields in the target list must appear at the same offset
	 * as in the input list.
	 */
	input_index = 0;
	
	if (strchr(nfl_sort[i].field_name, '(') == NULL) {
	    /* no default value */
	    search_name = (char *) nfl_sort[i].field_name;
	} else {
	    tmp_field_name = malloc(strlen(nfl_sort[i].field_name)); /* certainly big enough */
	    field_name_strip_get_default(&nfl_sort[i], tmp_field_name, &default_val);
	    search_name = tmp_field_name;
	}
	input_index = find_field_for_conv(search_name, input_field_list, i, NULL);
	if (input_index == -1) {
	    if(default_val){
	        if ((conv == buffer_and_convert) || 
		    (conv == copy_dynamic_portion)) {
		    input_index = -1; /* Basically invalidating input_index
				   Indication for using default_val */
		    break;
		} else {
		    if (default_val) {
		        free(default_val);
			default_val = NULL;
		    }
		    conv = buffer_and_convert;
		    if (tmp_field_name) free(tmp_field_name);
		    goto restart;
		}
	    }
	    fprintf(stderr,
		    "Requested field %s missing from input format\n",
		    nfl_sort[i].field_name);
	    FFSfree_conversion(conv_ptr);
	    if (tmp_field_name) free(tmp_field_name);
	    return NULL;
	}
	if (tmp_field_name) free(tmp_field_name);
	if(input_index == -1){
	    create_default_conversion(nfl_sort[i], default_val, &conv_ptr, 
				      conv_index);
	    conv_ptr->conversion_type = conv;
	    conv_ptr->conv_count = ++conv_index;
	    continue;
	} else {
	    if (default_val) {
		free(default_val);
		default_val = NULL;
	    }
	}
	input_field = input_field_list[input_index];
	in_data_type = FMarray_str_to_data_type(input_field.field_type,
					      &in_elements);
	if (in_elements != 1) {
	    char *first_bracket = memchr(input_field.field_type, '[', strlen(input_field.field_type));
	    if (first_bracket != NULL) {
		char *sec_bracket = memchr(first_bracket+1, '[', strlen(first_bracket));
		if (sec_bracket) multi_dimen_array++;
	    }
	}
	    
	if (in_elements == -1) {
	    in_elements = 1;	/* var array */
	}
	target_data_type = FMarray_str_to_data_type(nfl_sort[i].field_type,
						  &target_elements);
	if (target_elements == -1) {
	    target_elements = 1;	/* var array */
	}
	switch (conv) {
	case none_required:
	case direct_to_mem:
	    if (nfl_sort[i].field_offset !=
		input_field_list[input_index].field_offset) {

		/* 
		 * we were planning direct to memory transfer, but found a 
		 * field that we couldn't do that way.   start over
		 * planning buffer and convert. 
		 */
		conv = buffer_and_convert;
		goto restart;
	    }
	    if ((in_elements == target_elements) &&
		(in_data_type == target_data_type) &&
		(in_data_type != unknown_type) &&
		(!byte_reversal) &&
		((in_data_type != float_type) ||
		 (target_fp_format == src_float_format)) &&
		field_is_flat(src_ioformat->body, 
			      &input_var_list[input_index].type_desc) &&
		(!column_row_swap_necessary || !multi_dimen_array) &&
		(nfl_sort[i].field_size == input_field.field_size)) {
		/* nothing to do for this field */
		continue;
	    }
	    if ((in_elements != 1) &&
		(nfl_sort[i].field_size != input_field.field_size)) {
		/* Can't do direct to memory with array field interleaving 
		 * 
		 */
		conv = buffer_and_convert;
		goto restart;
	    }
	    if (column_row_swap_necessary && multi_dimen_array) {
		/* Can't transpose an array in place
		 * 
		 */
		conv = copy_dynamic_portion;
		goto restart;
	    }

	    if ((in_data_type == unknown_type) && 
		!input_var_list[input_index].var_array &&
		(input_var_list[input_index].type_desc.type != FMType_pointer)) {
		FFSTypeHandle format = src_ioformat->field_subformats[input_index];
		if ((format != NULL) && (format->conversion != NULL)) {
		    switch (format->conversion->conversion_type) {
		    case copy_dynamic_portion:
			assert(conv == copy_dynamic_portion);
			break;
		    case buffer_and_convert:
			conv = buffer_and_convert;
			goto restart;
		    case direct_to_mem:
			conv = direct_to_mem;
			break;
		    case none_required:
			continue;
		    }
		}
	    }
	    /* falling through */
	case buffer_and_convert:
	    if (input_var_list[input_index].var_array ||
	    	(input_var_list[input_index].type_desc.type == FMType_pointer)) {
	    	if (nfl_sort[i].field_size != input_field.field_size) {
	    	    /* argh.  Must buffer variant part too */
	    	    conv = copy_dynamic_portion;
	    	    goto restart;
	    	}
	    }
	case copy_dynamic_portion:
	    if (input_var_list[input_index].var_array ||
		(input_var_list[input_index].type_desc.type == FMType_pointer)) {
		/* 
		 * expansion value includes padding for possibly having
		 * to re-align the value to the proper boundary.
		 * (This is generally a gross overestimate, but it's safe.)
		 */
		int local_size = nfl_sort[i].field_size +
		min_align_size(nfl_sort[i].field_size);
		double expansion = ((double) local_size) / input_field.field_size;
		if (expansion > conv_ptr->max_var_expansion) {
		    conv_ptr->max_var_expansion = expansion;
		}
	    }
	    break;
	default:
	    assert(FALSE);
	}
	if ((in_elements != target_elements) || 
	    (in_data_type != target_data_type)) {

	    fprintf(stderr,
		    "Requested field %s base type %s \n   differs from source type %s\n\n",
		    nfl_sort[i].field_name, nfl_sort[i].field_type,
		    input_field.field_type);
	    FFSfree_conversion(conv_ptr);
	    return NULL;
	}
	if (conv == none_required) {
	    conv = direct_to_mem;
	}
	conv_ptr->conversion_type = conv;
	memset(&conv_ptr->conversions[conv_index].src_field, 0, 
	       sizeof(conv_ptr->conversions[conv_index].src_field));
	conv_ptr->conversions[conv_index].src_field.byte_swap = byte_reversal;
	conv_ptr->conversions[conv_index].src_field.src_float_format =
	    src_float_format;
	conv_ptr->conversions[conv_index].src_field.target_float_format =
	    target_fp_format;
	conv_ptr->conversions[conv_index].subconversion = NULL;
	conv_ptr->conversions[conv_index].iovar = NULL;
	conv_ptr->conversions[conv_index].rc_swap = no_row_column_swap;
	if (column_row_swap_necessary & multi_dimen_array) {
	    if (src_ioformat->body->column_major_arrays) {
		conv_ptr->conversions[conv_index].rc_swap = swap_source_column_major;
	    } else {
		conv_ptr->conversions[conv_index].rc_swap = swap_source_row_major;
	    }
	}
	in_data_type = FMarray_str_to_data_type(input_field.field_type, 
					      &in_elements);
	conv_ptr->conversions[conv_index].iovar =
	    &input_var_list[input_index];
	if (in_data_type == unknown_type) {
	    FFSTypeHandle format = src_ioformat->field_subformats[input_index];
	    if (format != NULL) {
		int format_index = 0;
		char *base_type = base_data_type(input_field.field_type);

		while(target_list[format_index].format_name && 
		      (strcmp(base_type, target_list[format_index].format_name) != 0)) {
		    format_index++;
		}
		free(base_type);
	        if (format == src_ioformat) {
		    conv_ptr->conversions[conv_index].subconversion = conv_ptr;
		} else if (target_list[format_index].format_name != NULL) {
		    IOConversionPtr subconv;
		    int target_struct_size = target_list[format_index].struct_size;
		    subconv =
		      create_conversion(format,
					target_list[format_index].field_list,
					target_struct_size, pointer_size,
					byte_reversal, target_fp_format,
					conv, target_column_major,
					string_offset_size,
					converted_strings, target_list);
		    if (conv_ptr->conversions[conv_index].subconversion) {
			printf("Subconversion already had value!\n");
		    }
		    conv_ptr->conversions[conv_index].subconversion = subconv;
		}
	    } else {
		
		fprintf(stderr, "Unknown field type for field %s ->\"%s\", format %p\n",
			input_field.field_name,
			src_ioformat->body->field_list[input_index].field_type,
			src_ioformat);
		FFSfree_conversion(conv_ptr);
		return NULL;
	    }
	}
	if (input_var_list[input_index].var_array) {
	    in_elements = 1;
	}
	conv_ptr->conversions[conv_index].src_field.data_type = in_data_type;
	conv_ptr->conversions[conv_index].src_field.offset =
	    input_field.field_offset;
	conv_ptr->conversions[conv_index].src_field.size =
	    input_field.field_size;
	conv_ptr->conversions[conv_index].dest_size =
	    nfl_sort[i].field_size;
	conv_ptr->conversions[conv_index].dest_offset =
	    nfl_sort[i].field_offset;
	conv_ptr->conversions[conv_index].default_value = NULL;
	conv_ptr->conv_count = ++conv_index;
	if (default_val) {
	    free(default_val);
	    default_val = NULL;
	}
    }
    conv_ptr->conv_func = generate_conversion(conv_ptr, src_ioformat->body->alignment, 4);
    switch(conv_ptr->required_alignment) {
    case 1:
	conv_ptr->conv_func1 = conv_ptr->conv_func;
	/* falling */
    case 2:
	conv_ptr->conv_func2 = conv_ptr->conv_func;
	/* falling */
    case 4:
	conv_ptr->conv_func4 = conv_ptr->conv_func;
	/* falling */
    case 8:
    case 16:
    case 32:
    case 64:
	/* really can't imagine these currently, but... */
    case 0:
	/* zero is no conversion, so no requirement */
	break;
    default:
	fprintf(stderr, "Funky alignment, %d, for conversion %s\n",
		conv_ptr->required_alignment, 
		conv_ptr->ioformat->body->format_name);
    }
    return conv_ptr;
}

extern
void
set_general_IOconversion_for_format(FFSContext iofile, FFSTypeHandle file_ioformat, FMFieldList native_field_list,
				    int native_struct_size, int pointer_size, FMStructDescList target_list)
{
    IOConversionPtr conv_ptr;
    IOconversion_type conv_type = none_required;
    int string_offset_size;
    int data_align_pad = (8 - file_ioformat->body->record_length) & 0x7;

    if (file_ioformat->body->byte_reversal)
	conv_type = direct_to_mem;

    string_offset_size = file_ioformat->body->record_length + data_align_pad;

    conv_ptr = create_conversion(file_ioformat, native_field_list,
				 native_struct_size, pointer_size,
				 file_ioformat->body->byte_reversal, 
				 ffs_my_float_format, conv_type,
/*				 iofile->native_column_major_arrays*/ 0,
				 string_offset_size, FALSE, target_list);

    if (conv_ptr == NULL) {
	fprintf(stderr, "Set_IOconversion failed for format name %s\n",
		file_ioformat->body->format_name);
	return;
    }
    conv_ptr->context = iofile;
    if (file_ioformat->conversion != NULL) {
	FFSfree_conversion(file_ioformat->conversion);
    }
    file_ioformat->conversion = conv_ptr;
}

static int
set_conversion_from_list(FFSContext iocontext, FFSTypeHandle ioformat, FMStructDescList struct_list)
{
    int i = 0;
    while(struct_list[i].format_name != NULL) {
	if (strcmp(struct_list[i].format_name, ioformat->body->format_name) == 0) {
	    break;
	}
	i++;
    }
    if (struct_list[i].format_name == NULL) {
	printf("Local structure description for type \"%s\" not found in IOStructDescList\n",
	       ioformat->body->format_name);
	return 0;
    }

    set_general_IOconversion_for_format(iocontext, ioformat,
					struct_list[i].field_list, 
					struct_list[i].struct_size,
					(int) sizeof(char *), struct_list);
    return 1;
}

#ifdef HAVE_DILL_H
static void
link_conversion_sets(FFSTypeHandle ioformat)
{
    dill_extern_entry* externs;
    int i, subformat_count = 0;
    while (ioformat->subformats && ioformat->subformats[subformat_count])
	subformat_count++;
    externs = malloc(sizeof(externs[0]) * (subformat_count + 2));
    
    for(i=0; i < subformat_count; i++) {
	FFSTypeHandle th = ioformat->subformats[i];
	externs[i].extern_name = FFSTypeHandle_name(th);
	if (th->conversion->conv_pkg != NULL) {
	    externs[i].extern_value = dill_package_entry(th->conversion->conv_pkg);
	} else {
	    externs[i].extern_value = (void *)th->conversion->conv_func;
	}
    }
    externs[subformat_count].extern_name = FFSTypeHandle_name(ioformat);
    if (ioformat->conversion->conv_pkg != NULL) {
	externs[subformat_count].extern_value = 
	    dill_package_entry(ioformat->conversion->conv_pkg);
    } else {
	externs[subformat_count].extern_value = (void*)ioformat->conversion->conv_func;
    }
    externs[subformat_count+1].extern_name = NULL;
    externs[subformat_count+1].extern_value = NULL;
}
#else
static void link_conversion_sets(FFSTypeHandle ioformat){}
#endif

extern
void
establish_conversion(FFSContext iocontext, FFSTypeHandle ioformat, FMStructDescList struct_list)
{
    int use_package = 0;

    if (!(set_conversion_from_list(iocontext, ioformat, struct_list))) {
	return;
    }
    if (use_package) {
	link_conversion_sets(ioformat);
    }
}

extern void
ffs_internal_convert_field(FMFieldPtr src_spec, void *src, FMdata_type dest_type, int dest_size, void *dest)
{
    int float_OK = 1;
    if (dest_type == float_type) {
	FMfloat_format reverse_src = 
	    ffs_reverse_float_formats[src_spec->src_float_format];

	if (!src_spec->byte_swap && 
	    src_spec->src_float_format != src_spec->target_float_format) {
	    float_OK = 0;
	} else if (src_spec->byte_swap &&
		   reverse_src != src_spec->target_float_format) {
	    float_OK = 0;
	}
    }		   
    /* quick check to see if it's just a copy... */
    if ((dest_type != string_type) && (dest_type == src_spec->data_type) &&
	(dest_size == src_spec->size) && (float_OK == 1)) {
	if (src_spec->byte_swap) {
	    int i;
	    char *destc = (char *) dest;
	    char *srcc = (char *) src + src_spec->offset;
	    for (i = 0; i < (dest_size >> 1); i++) {
		char tmp = srcc[dest_size - i - 1];
		destc[dest_size - i - 1] = srcc[i];
		destc[i] = tmp;
	    }
	    if ((dest_size & 0x1) != 0) {
		destc[dest_size >> 1] = srcc[dest_size >> 1];
	    }
	} else {
	    char *srcc = (char *) src + src_spec->offset;
	    if (dest != srcc) memcpy(dest, srcc, dest_size);
	}
	return;
    }
    switch (dest_type) {
    case integer_type:
	{
	    MAX_INTEGER_TYPE tmp = get_big_int(src_spec, src);
	    if (dest_size == sizeof(char)) {
		char *dest_field = (char *) dest;
		*dest_field = (char) tmp;
	    } else if (dest_size == sizeof(short)) {
		short *dest_field = (short *) dest;
		*dest_field = (short) tmp;
	    } else if (dest_size == sizeof(int)) {
		int *dest_field = (int *) dest;
		*dest_field = (int) tmp;
	    } else if (dest_size == sizeof(long)) {
		long* dest_field = (long*)dest;
		*dest_field = (long)tmp;
	    } else if (dest_size == sizeof(size_t)) {
		size_t* dest_field = (size_t*)dest;
		*dest_field = (long)tmp;
#if SIZEOF_LONG_LONG != 0
	    } else if (dest_size == sizeof(long long)) {
		long long lltmp = tmp;
		memcpy(dest, &lltmp, sizeof(long long));
#endif
	    }
	    break;
	}
    case char_type:
	{
	    char tmp = (char)get_big_int(src_spec, src);
	    char *dest_field = (char *) dest;
	    *dest_field = tmp;
	    break;
	}
    case boolean_type:
    case enumeration_type:
	{
	    int tmp = (int) get_big_int(src_spec, src);
	    int *dest_field = (int *) dest;
	    *dest_field = tmp;
	    break;
	}
    case unsigned_type:
	{
	    MAX_UNSIGNED_TYPE tmp = get_big_unsigned(src_spec, src);
	    if (dest_size == sizeof(unsigned char)) {
		unsigned char *dest_field = (unsigned char *) dest;
		*dest_field = (unsigned char) tmp;
	    } else if (dest_size == sizeof(unsigned short)) {
		unsigned short *dest_field = (unsigned short *) dest;
		*dest_field = (unsigned short) tmp;
	    } else if (dest_size == sizeof(int)) {
		unsigned int* dest_field = (unsigned int*)dest;
		*dest_field = (unsigned int)tmp;
	    } else if (dest_size == sizeof(size_t)) {
		size_t* dest_field = (size_t*)dest;
		*dest_field = tmp;
#if SIZEOF_LONG_LONG != 0
	    } else if (dest_size == sizeof(long long)) {
		unsigned long long *dest_field = (unsigned long long *) dest;
		*dest_field = tmp;
#endif
	    }
	    break;
	}
    case float_type:
	{
	    MAX_FLOAT_TYPE tmp = get_big_float(src_spec, src);
	    if (dest_size == sizeof(float)) {
		float *dest_field = (float *) dest;
		*dest_field = (float) tmp;
	    } else if (dest_size == sizeof(double)) {
		double *dest_field = (double *) dest;
		*dest_field = tmp;
#if SIZEOF_LONG_DOUBLE != 0
	    } else if (dest_size == sizeof(long double)) {
		long double *dest_field = (long double *) dest;
		*dest_field = tmp;
#endif
	    }
	    break;
	}
    default:
	assert(FALSE);
    }
}

extern long
get_static_array_element_count(FMVarInfoList var)
{
    int i;
    long count = 1;
    if (var == NULL) return 1;
    for (i = 0; i < var->dimen_count; i++) {
	if (var->dimens[i].static_size != 0) {
	    count = count * var->dimens[i].static_size;
	} else {
	    return -1;
	}
    }
    return count;
}

static int _debug_code_generation = -1;

static inline int NO_SANITIZE_THREAD debug_code_generation(void)
{
    return _debug_code_generation;
}

typedef struct conv_status {
    void *src;
    void *dest;
    void *src_pointer_base;
    void *dest_pointer_base;
    ssize_t   src_offset_adjust;
    ssize_t   dest_offset_adjust;
    size_t   cur_offset;
    size_t   *control_value;
    int target_pointer_size;
    int src_pointer_size;
    int register_args;
    IOConversionPtr global_conv;
} *ConvStatus;

typedef struct run_time_conv_status {
    void *src_pointer_base;
    void *dest_pointer_base;
    ssize_t   dest_offset_adjust;
} *RTConvStatus;

static void
internal_convert_record(IOConversionPtr conv, ConvStatus conv_status, 
			void *src, void *dest, int data_already_copied);

static void
print_IOConversion(IOConversionPtr conv_ptr, int indent)
{
    int i;
    int ind;
    if (indent == 0) {
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("IOConversion base type is ");
	if (conv_ptr == NULL) {
	    printf("NULL\n");
	    return;
	}
	switch (conv_ptr->conversion_type) {
	case none_required:
	    printf("None_Required\n");
	    break;
	case direct_to_mem:
	    printf("Direct_to_Memory\n");
	    break;
	case buffer_and_convert:
	    printf("Buffer_and_Convert\n");
	    break;
        case copy_dynamic_portion:
	    printf("Copy_Dynamic_Portion\n");
	    break;
	default:
	    assert(FALSE);
	    break;
	}
    }
    for (ind = 0; ind < indent; ind++)
	printf("    ");
    printf(" base_size_delta=%zd, max_var_exp=%g, target_pointer_size=%d, string_offset=%zd, converted_strings=%d\n",
	   conv_ptr->base_size_delta, conv_ptr->max_var_expansion,
	   conv_ptr->target_pointer_size, conv_ptr->string_offset_size,
	   conv_ptr->converted_strings);
    printf(" conversion_function= %p, required_align=%d\n",
	   conv_ptr->conv_func, conv_ptr->required_alignment);
    for (ind = 0; ind < indent; ind++)
	printf("    ");
    printf("  There are %d conversions registered:\n", conv_ptr->conv_count);
    for (i = 0; i < conv_ptr->conv_count; i++) {
	FMVarInfoStruct *iovar = conv_ptr->conversions[i].iovar;
	FMFieldPtr src_field = &conv_ptr->conversions[i].src_field;
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("  Conversion %d:\n", i);
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("    Base type : %s", data_type_to_str(src_field->data_type));
	if (conv_ptr->conversions[i].iovar != NULL) {
	    int j;
	    for (j = 0; j < iovar->dimen_count; j++) {
		if (iovar->dimens[j].static_size != 0) {
		    printf("[%d]", iovar->dimens[j].static_size);
		} else {
		    /* variant array */
		    int field = iovar->dimens[j].control_field_index;
		    FMFormat f = conv_ptr->ioformat->body;
		    int offset = f->field_list[field].field_offset;
		    int size = f->field_list[field].field_size;
		    printf("[ size at offset %d, %dbytes ]", offset, size);
		}
	    }
	}
	
	if (conv_ptr->conversions[i].rc_swap == swap_source_column_major) {
	    printf(" row/column swap required (SRC column-major) - ");
	}
	if (conv_ptr->conversions[i].rc_swap == swap_source_row_major) {
	    printf(" row/column swap required (SRC row-major) - ");
	}
	if (src_field->byte_swap) {
	    printf(" byte order reversal required\n");
	} else {
	    printf("\n");
	}
	if ((src_field->data_type == float_type) &&
	    (src_field->src_float_format != src_field->target_float_format)) {
	    printf("conversion from %s to %s required\n", 
		   float_format_str[(int)src_field->src_float_format],
		   float_format_str[(int)src_field->target_float_format]);
	}
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	if (conv_ptr->conversions[i].default_value == NULL) {
	    printf("    Src offset : %zd    size %d\n", src_field->offset,
		   src_field->size);
	} else {
	    int j;
	    printf("    Default value : 0x");
	    for (j=0; j< conv_ptr->conversions[i].dest_size; j++) {
		printf("%02x", ((unsigned char*)conv_ptr->conversions[i].default_value)[j]);
	    }
	    printf("\n");
	}
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("    Dst offset : %zd    size %d\n",
	       conv_ptr->conversions[i].dest_offset,
	       conv_ptr->conversions[i].dest_size);
	if (conv_ptr->conversions[i].subconversion) {
	    for (ind = 0; ind < indent; ind++)
		printf("    ");
	    if (conv_ptr->conversions[i].subconversion == conv_ptr) {
		printf("    Subconversion is recursive\n");
	    } else {
		printf("    Subconversion as follows:\n");
		print_IOConversion(conv_ptr->conversions[i].subconversion,
				   indent + 1);
	    }
	}
    }
}

static void
print_IOConversion_as_XML(IOConversionPtr conv_ptr, int indent)
{
    int i;
    int ind;
    if (indent == 0) {
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("<IOConversion baseType=\"");
	if (conv_ptr == NULL) {
	    printf("NULL\" />");
	    return;
	}
	switch (conv_ptr->conversion_type) {
	case none_required:
	    printf("None_Required");
	    break;
	case direct_to_mem:
	    printf("Direct_to_Memory");
	    break;
	case buffer_and_convert:
	    printf("Buffer_and_Convert");
	    break;
	case copy_dynamic_portion:
	    printf("Copy_Strings");
	    break;
	default:
	    assert(FALSE);
	    break;
	}
	printf("\">\n");
    }
    for (ind = 0; ind < indent; ind++)
	printf("    ");
    printf("<baseSizeDelta>%zd</baseSizeDelta>\n", conv_ptr->base_size_delta);
    printf("<maxVarExpansion>%g</maxVarExpansion>\n", conv_ptr->max_var_expansion);
    printf("<targetPointerSize>%d</targetPointerSize>\n", conv_ptr->target_pointer_size);
    printf("<stringOffsetSize>%zd</stringOffsetSize>\n", conv_ptr->string_offset_size);
    printf("<convertedStrings>%d</convertedStrings>\n", conv_ptr->converted_strings);
    for (ind = 0; ind < indent; ind++)
	printf("    ");
    for (i = 0; i < conv_ptr->conv_count; i++) {
	FMFieldPtr src_field = &conv_ptr->conversions[i].src_field;
	FMVarInfoStruct *iovar = conv_ptr->conversions[i].iovar;

	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("<registeredConversion>\n");
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("<baseType>%s</baseType>\n",
	 data_type_to_str(src_field->data_type));
	printf("<controlField>\n");
	if (conv_ptr->conversions[i].iovar != NULL) {
	    int j;
	    for (j = 0; j < iovar->dimen_count; j++) {
		if (iovar->dimens[j].static_size != 0) {
		    printf("<arrayDimension>%d</arrayDimension>", 
			   iovar->dimens[j].static_size);
		} else {
		    int field = iovar->dimens[j].control_field_index;
		    FMFormat f = conv_ptr->ioformat->body;
		    int offset = f->field_list[field].field_offset;
		    int size = f->field_list[field].field_size;
		    /* variant array */
		    printf("<offset>%d</offset><size units=\"bytes\">%d</size>\n",
			   offset, size);
		}
	    }
	}
	printf("</controlField>\n");
	if (src_field->byte_swap) {
	    printf("<byteReversal />\n");
	} else {
	    printf("\n");
	}
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("<sourceOffset>%zd</sourceOffset><sourceSize>%d</sourceSize>\n",
	       src_field->offset,
	       src_field->size);
	for (ind = 0; ind < indent; ind++)
	    printf("    ");
	printf("<destOffset>%zd</destOffset><destSize>%d</destSize>\n",
	       conv_ptr->conversions[i].dest_offset,
	       conv_ptr->conversions[i].dest_size);
	if (conv_ptr->conversions[i].subconversion) {
	    for (ind = 0; ind < indent; ind++)
		printf("    ");
	    if (conv_ptr->conversions[i].subconversion == conv_ptr) {
		printf("    Subconversion is recursive\n");
	    } else {
		print_IOConversion_as_XML(conv_ptr->conversions[i].subconversion,
					  indent + 1);
	    }
	}
    }
    printf("</IOConversion>\n");
}


extern void
dump_IOConversion(IOConversionPtr conv_ptr)
{
    print_IOConversion(conv_ptr, 0);
}

extern void
dump_IOConversion_as_XML(IOConversionPtr conv_ptr)
{
    print_IOConversion_as_XML(conv_ptr, 0);
}

void
FFSconvert_record(IOConversionPtr conv, void *src, void *dest, void *final_string_base, void *src_string_base)
{
    struct conv_status cs;
    if (src_string_base == NULL) {
	src_string_base = final_string_base;
    }
    if (conv->conv_func) {
	struct run_time_conv_status rtcs;
	if (debug_code_generation()) {
	    int i;
	    int limit = 30;
	    int *tmp = (int *) (((char *) src_string_base) -
				(((intptr_t) src_string_base) % 4));
	    printf("record of type \"%s\", contents :\n", 
		   conv->ioformat->body->format_name);
	    if (limit * sizeof(int) > conv->ioformat->body->record_length)
		limit = conv->ioformat->body->record_length / sizeof(int);
	    for (i = 0; i < limit; i += 4) {
		printf("%p: %8x %8x %8x %8x\n", ((char *) src) + (i * 4),
		       ((int *) src)[i], ((int *) src)[i + 1],
		       ((int *) src)[i + 2], ((int *) src)[i + 3]);
	    }
	    if (src_string_base != NULL) {
		printf("string contents :\n");
		limit = 10;
/*		if (conv->ioformat->body->variant) {
		    FILE_INT record_len;
		    int len_align_pad = (4 - conv->ioformat->body->server_ID.length) & 3;
		    FILE_INT *len_ptr = (FILE_INT *) (src + conv->ioformat->body->server_ID.length +
						      len_align_pad);
		    memcpy(&record_len, len_ptr, sizeof(FILE_INT));
		    if (conv->ioformat->body->byte_reversal)
			byte_swap((char *) &record_len, 4);
		    record_len -= conv->ioformat->body->record_length;
		    if (limit * sizeof(int) > record_len) {
			limit = record_len / sizeof(int);
		    }
		}
*/		for (i = 0; i < limit; i += 4) {
		    printf("%p: %8x %8x %8x %8x\n", ((char *) tmp) + (i * 4),
			   ((int *) tmp)[i],
			   ((int *) tmp)[i + 1],
			   ((int *) tmp)[i + 2],
			   ((int *) tmp)[i + 3]);
		}
	    }
	}
	rtcs.src_pointer_base = src_string_base;
	rtcs.dest_pointer_base = final_string_base;
	rtcs.dest_offset_adjust = -(ssize_t) conv->string_offset_size;
	conv->conv_func(src, dest, &rtcs);
	return;
    } else {
	cs.src_pointer_base = src_string_base;
	cs.dest_pointer_base = final_string_base;
	cs.src_offset_adjust = -(ssize_t)conv->string_offset_size;
	cs.dest_offset_adjust = -(ssize_t)conv->string_offset_size;
	cs.cur_offset = 0;
	cs.control_value = NULL;
	cs.target_pointer_size = conv->target_pointer_size;
	cs.src_pointer_size = conv->ioformat->body->pointer_size;
	cs.global_conv = conv;
	internal_convert_record(conv, &cs, src, dest, 1);
    }
}

void
transpose_array(size_t *dimens, char *src, char *dest, int source_column_major,
		FMdata_type dest_type, int dest_size, 
		FMFieldPtr src_spec)
{
    int dimen_count = 0;
    size_t *index;
    size_t i;
    ssize_t cur_index;
    size_t jump = 1;

    while (dimens[dimen_count] != 0) dimen_count++;
    struct _FMgetFieldStruct tmp_spec = *src_spec;

    if (dimen_count <= 1) return;
    index = malloc(sizeof(index[0]) * dimen_count);
    for(i = 0; i< dimen_count; i++) {
	index[i] = 0;
    }
    cur_index = 0;
    jump = 1;
    for (i = 0; i < dimen_count-1; i++) {
	jump = (jump * dimens[i]);
    }
    while(index[0] < dimens[0]) {
	if (cur_index == (dimen_count-1)) {
	    size_t col_index_base = 0;
	    size_t row_index_base = 0;
	    int i;
	    void * dest_field;
	    if (dimen_count >= 2) {
		col_index_base = index[dimen_count-1];
		for (i = dimen_count-1; i >= 0; i--) {
		    col_index_base = (col_index_base * dimens[i] + index[i]);
		}
		row_index_base = index[0];
		for (i = 1; i < dimen_count; i++) {
		    row_index_base = (row_index_base * dimens[i] + index[i]);
		}
	    }
	    if (source_column_major) {
		dest_field = ((char*)dest) + dest_size * row_index_base;
		tmp_spec.offset = tmp_spec.size * col_index_base;
	    } else {
		dest_field = ((char*)dest) + dest_size * col_index_base;
		tmp_spec.offset = tmp_spec.size * row_index_base;
	    }
	    for(i=0; i < dimens[cur_index]; i++) {
		if (dest_type != unknown_type) {
		    /* simple (native) field or variant array */
		    if (dest_type != string_type) {
/*GSE var char blocks here*/
			ffs_internal_convert_field(&tmp_spec, src,
						   dest_type, dest_size,
						   dest_field);
		    } else {
			printf("Bad type in transpose\n");
			free(index);
			return;
		    }
		} else {
		    printf("Bad type in transpose\n");
		    free(index);
		    return;
		}
		if (source_column_major) {
		    dest_field = (char*)dest_field + dest_size;
		    tmp_spec.offset += tmp_spec.size * jump;
		} else {
		    dest_field = (char*)dest_field + (dest_size * jump);
		    tmp_spec.offset += tmp_spec.size;
		}
		col_index_base += jump;
		row_index_base++;
	    }
	    cur_index--;
	    index[cur_index]++;
	} else {
	    if (index[cur_index] == dimens[cur_index]) {
		index[cur_index] = 0;
		cur_index--;
		if (cur_index == -1) {
		    free(index);
		    return;
		}
		index[cur_index]++;
	    } else {
		cur_index++;
	    }
	}
    }
    free(index);
}

static size_t
get_offset_for_addr(char *src_field_addr, ConvStatus conv_status, 
		    IOconvFieldStruct *conv)
{
    struct _FMgetFieldStruct tmp_src_field;  /* OK */
    MAX_INTEGER_TYPE tmp_int;
    tmp_src_field = conv->src_field;
    tmp_src_field.size = conv_status->src_pointer_size;
    tmp_src_field.data_type = integer_type;
    tmp_src_field.offset = 0;

    tmp_int = get_big_int(&tmp_src_field, src_field_addr);
    return tmp_int;
}
    
static void
convert_address_field(char *src_field_addr, char **output_source_ptr,
		      char *dest_field_addr, char **output_dest_ptr,
		      ConvStatus conv_status, IOconvFieldStruct *conv,
		      int required_alignment)
{
    int dest_size = conv_status->target_pointer_size;
    char *output_dest = NULL;
    char *output_source = NULL;
    size_t offset = get_offset_for_addr(src_field_addr, conv_status,
				     conv);
    if (offset != 0) {
	int align_tmp = 0;
	output_source = (char*)conv_status->src_pointer_base + offset + 
	  conv_status->src_offset_adjust;
	/* handle possibly different string base */
	output_dest = (char*)conv_status->dest_pointer_base + offset + 
	    conv_status->dest_offset_adjust;
	
	if ((align_tmp = (((uintptr_t)output_dest) % required_alignment)) != 0) {
	    output_dest += (required_alignment - align_tmp);
	    conv_status->dest_offset_adjust += (required_alignment - align_tmp);
	}
    }

    *output_source_ptr = output_source;
    *output_dest_ptr = output_dest;
    if (dest_size == sizeof(char *)) {
	(*(char**)dest_field_addr) = output_dest;   /* set field in dest */
    }
}


static void
new_convert_address_field(size_t offset, char **output_source_ptr,
		      char *dest_field_addr, char **output_dest_ptr,
		      ConvStatus conv_status, int required_alignment)
{
    int dest_size = conv_status->target_pointer_size;
    char *output_dest = NULL;
    char *output_source = NULL;

    if (offset != 0) {
	int align_tmp = 0;
	output_source = (char*)conv_status->src_pointer_base + offset + 
	    conv_status->src_offset_adjust;
	/* handle possibly different string base */
	output_dest = (char*)conv_status->dest_pointer_base + offset + 
	    conv_status->dest_offset_adjust;
	
	if ((align_tmp = (((uintptr_t)output_dest) % required_alignment)) != 0) {
	    output_dest += (required_alignment - align_tmp);
	    conv_status->dest_offset_adjust += (required_alignment - align_tmp);
	}
    }

    *output_source_ptr = output_source;
    *output_dest_ptr = output_dest;
    if (dest_size == sizeof(char *)) {
	(*(char**)dest_field_addr) = output_dest;   /* set field in dest */
    }
}

static ssize_t
decode_size_delta(ConvStatus conv_status, 
		  IOconvFieldStruct *conv, FMTypeDesc *type_desc)
{
    switch(type_desc->type) {
    case FMType_string:
    case FMType_pointer:
	return (conv_status->target_pointer_size - conv_status->src_pointer_size);
    case FMType_subformat:
	return conv->subconversion->base_size_delta;
    case FMType_simple:
	return(conv->dest_size - conv->src_field.size);
    default:
    case FMType_array:
	assert(0);
    }
    assert(0);
    return 0;
}

static int
item_size(ConvStatus conv_status, 
	   IOconvFieldStruct *conv, FMTypeDesc *type_desc)
{
    switch(type_desc->type) {
    case FMType_string:
    case FMType_pointer:
	return conv_status->src_pointer_size;
    case FMType_subformat:
	return conv->subconversion->ioformat->body->record_length;
    case FMType_simple:
	return conv->src_field.size;
    case FMType_array:
    default:
	assert(0);
    }
    return 0;
}

static void
new_convert_field(char *src_field_addr, char *dest_field_addr, 
		  ConvStatus conv_status, 
		  IOconvFieldStruct *conv, FMTypeDesc *type_desc,
		  int data_already_copied)
{
    switch(type_desc->type) {
    case FMType_pointer: {
	char *new_src, *new_dest;
	size_t offset = get_offset_for_addr(src_field_addr, conv_status, 
					 conv);
	new_convert_address_field(offset, &new_src, dest_field_addr, 
			      &new_dest, conv_status, 8);
	if (new_dest == NULL) break;
	if (offset != 0) {
	    /* if the offset is smaller than where we were working, we went backwards because of recursion */
	    if (offset <= conv_status->cur_offset) {
		break;
	    }
	}
	conv_status->cur_offset = offset;
	/* at this point...  We should tweak conv_status->dest_offset_adjust 
	 * IFF the src and dest sizes of the thing pointed to by this 
	 * pointer are different.   Also, memcpy src to dest if necessary.
	 */
	new_convert_field(new_src, new_dest, conv_status, conv,
			  type_desc->next, 0);
	break;
    }
    case FMType_string:{
	char *new_src, *new_dest;
	convert_address_field(src_field_addr, &new_src, dest_field_addr, 
			      &new_dest, conv_status, conv, 1);
	if (new_src != new_dest) { 
	    strcpy(new_dest, new_src);
	}
	break;
    }
    case FMType_subformat: {
	IOConversionPtr subtype_conv = conv->subconversion;
	struct conv_status cs;
	
	cs.src_pointer_base = conv_status->src_pointer_base;
	cs.dest_pointer_base = conv_status->dest_pointer_base;
	cs.src_offset_adjust = conv_status->src_offset_adjust;
	cs.dest_offset_adjust = conv_status->dest_offset_adjust;
	cs.cur_offset = conv_status->cur_offset;
	cs.control_value = NULL;
	cs.target_pointer_size = subtype_conv->target_pointer_size;
	cs.src_pointer_size = subtype_conv->ioformat->body->pointer_size;

	cs.global_conv = subtype_conv;
	internal_convert_record(subtype_conv, &cs, src_field_addr, 
				dest_field_addr, data_already_copied);
	conv_status->dest_offset_adjust = cs.dest_offset_adjust;
	conv_status->cur_offset = cs.cur_offset;
	break;
    }
    case FMType_simple: {
	struct _FMgetFieldStruct tmp_spec = conv->src_field;
	tmp_spec.offset = 0;
	ffs_internal_convert_field(&tmp_spec, src_field_addr,
				   tmp_spec.data_type, conv->dest_size,
				   dest_field_addr);
	break;
    }
    case FMType_array: {
	size_t elements = 1, i;
	char *new_src = src_field_addr;
	char *new_dest = dest_field_addr;
	FMTypeDesc *next = type_desc;
	while (next->type == FMType_array) {
	    if (next->static_size != 0) {
		elements *= next->static_size;
	    } else {
		elements *= conv_status->control_value[next->control_field_index];
		
	    }
	    next = next->next;
	}
	if (conv->rc_swap == no_row_column_swap) {
	    if (!data_already_copied) {
		ssize_t base_delta = decode_size_delta(conv_status, conv, next);
		size_t total_delta = base_delta * elements;
		conv_status->dest_offset_adjust += total_delta;
		if (conv_status->global_conv->conversion_type == copy_dynamic_portion) {
		    int base_size = item_size(conv_status, conv, next);
		    memcpy(dest_field_addr, src_field_addr, base_size * elements);
		}
		data_already_copied = 1;
	    }

	    int array_of_pointers = 0;
	    if (conv->iovar->type_desc.type == FMType_pointer) {
		array_of_pointers = (conv->iovar->type_desc.next->next->type == FMType_pointer);
	    }
	    /* changes here */
	    if (data_already_copied &&
		!conv->src_field.byte_swap &&
		(conv->src_field.src_float_format == conv->src_field.target_float_format) &&
		(conv->src_field.size == conv->dest_size) &&
		(conv->subconversion == NULL) &&
		(next->type != FMType_pointer) &&
		(next->type != FMType_string) &&
		((conv->src_field.data_type != string_type))) {
		/* Nothing but data movement required */
		    break;
	    }
	    for (i=0; i< elements ; i++) {
		new_convert_field(new_src, new_dest, conv_status, conv,
				  next, data_already_copied);
		if (array_of_pointers) {
		    new_src += conv_status->src_pointer_size;
		    new_dest += conv_status->target_pointer_size;
		} else {
		    new_src += conv->src_field.size;
		    new_dest += conv->dest_size;
		}
	    }
	} else {
	    /* row/column swap time */
	    int source_column_major = 
		(conv->rc_swap == swap_source_column_major);
	    int dimen_count = conv->iovar->dimen_count;
	    size_t *dimens = malloc(sizeof(dimens[0]) * (dimen_count + 1));
	    int i = 0;
	    FMdata_type dest_type = conv->src_field.data_type;
	    int dest_size = conv->dest_size;
	    void *dest_base = (void *) new_dest;
	    struct _FMgetFieldStruct tmp_spec = conv->src_field;
	    tmp_spec.offset = 0;

	    FMTypeDesc *next = type_desc;
	    while (next->type == FMType_array) {
		if (next->static_size != 0) {
		    dimens[i] = next->static_size;
		} else {
		    dimens[i] = conv_status->control_value[next->control_field_index];
		}
		i++;
		next = next->next;
	    }
	    dimens[dimen_count] = 0;
	    transpose_array(dimens, (char *) new_src, dest_base,
			    source_column_major, dest_type, dest_size,
			    &tmp_spec);
	    free(dimens);
	}
    }
	break;
    }
}

static void
internal_convert_record(IOConversionPtr conv, ConvStatus conv_status, void *src,
			void *dest, int data_already_copied)
{
    int i;
    size_t *control_value = NULL;
    if (conv->conversion_type == none_required) return;
    for (i = 0; i < conv->conv_count; i++) {
	FMTypeDesc *next = &conv->conversions[i].iovar->type_desc;
	while (next != NULL) {
	    if ((next->type == FMType_array) && 
		(next->static_size == 0)) {
		size_t elements;
		FMFormat f = conv->ioformat->body;
		int field = next->control_field_index;
		struct _FMgetFieldStruct tmp_src_spec;
		memset(&tmp_src_spec, 0, sizeof(tmp_src_spec));
		tmp_src_spec.size = f->field_list[field].field_size;
		tmp_src_spec.offset = f->field_list[field].field_offset;
		tmp_src_spec.data_type = integer_type;
		tmp_src_spec.byte_swap = conv->ioformat->body->byte_reversal;
		elements = get_big_int(&tmp_src_spec, src);
		if (control_value == NULL) {
		    int j;
		    control_value = (size_t *) malloc(sizeof(control_value[0]) * f->field_count);
		    for (j = 0; j < f->field_count; j++)
			control_value[j] = 0;
		    conv_status->control_value = control_value;
		}
		control_value[next->control_field_index] = elements;
	    }
	    next = next->next;
	}
    }
    for (i = 0; i < conv->conv_count; i++) {
	FMFieldPtr src_spec = &conv->conversions[i].src_field;
	FMTypeDesc *type_desc = &conv->conversions[i].iovar->type_desc;
	int elements = get_static_array_element_count(conv->conversions[i].iovar);
	int byte_swap = conv->conversions[i].src_field.byte_swap;
	char *src_field_addr = (char*)src + src_spec->offset;
	char *dest_field_addr = (char*)dest + conv->conversions[i].dest_offset;

	if (conv->conversions[i].src_field.size == 1) byte_swap = 0;
	if (conv->conversions[i].default_value){
	    void *dest_field = (void *) (conv->conversions[i].dest_offset + 
					 (char *) dest);
	    memcpy(dest_field, conv->conversions[i].default_value, 
		   conv->conversions[i].dest_size);
	} else if (!byte_swap &&
		   (src_spec->size == conv->conversions[i].dest_size) &&
		   (conv->conversions[i].subconversion == NULL) &&
		   (type_desc->type != FMType_pointer) &&
		   (type_desc->type != FMType_string) &&
		   (elements != -1 /* var array */) &&
		   (conv->conversions[i].rc_swap == no_row_column_swap) &&
		   ((src_spec->data_type != float_type) || 
		    (src_spec->src_float_format == src_spec->target_float_format)) &&
		   ((src_spec->data_type != string_type) || 
		    ((conv_status->dest_pointer_base) == NULL))) {
	    /* data movement is all that is required */
	    memcpy(dest_field_addr, src_field_addr,
		   src_spec->size * elements);
	    continue;
	}
	new_convert_field(src_field_addr, dest_field_addr, conv_status, 
			  &conv->conversions[i], type_desc, 1);
    }
    if (control_value != NULL) {
	free(control_value);
    }
}


static MAX_INTEGER_TYPE
get_big_int(FMFieldPtr iofield, void *data)
{
    if (iofield->data_type == integer_type) {
	if (iofield->size == sizeof(char)) {
	    char tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(char));
	    return (long) tmp;
	} else if (iofield->size == sizeof(short)) {
	    short tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(short));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(short));
	    return (long) tmp;
	} else if (iofield->size == sizeof(int)) {
	    int tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(int));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(int));
	    return (long) tmp;
	} else if (iofield->size == sizeof(long)) {
	    long tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(long));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(long));
	    return tmp;
	} else if (iofield->size == 2 * sizeof(long)) {
	    long tmp;
	    size_t low_bytes_offset = iofield->offset;
	    if (WORDS_BIGENDIAN) {
		if (!iofield->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (iofield->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    }
	    memcpy(&tmp, (char *) data + low_bytes_offset, sizeof(long));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(long));
	    return tmp;
	} else {
	    if (!IO_shut_up && !get_long_warn) {
		fprintf(stderr, "Get Long failed!  Size problems.  File int size is %d.\n", iofield->size);
		get_long_warn++;
	    }
	    return -1;
	}
    } else if (iofield->data_type == unsigned_type) {
	MAX_UNSIGNED_TYPE tmp = get_big_unsigned(iofield, data);
	return (MAX_UNSIGNED_TYPE) tmp;
    } else if (iofield->data_type == float_type) {
	MAX_FLOAT_TYPE tmp = get_big_float(iofield, data);
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
get_big_unsigned(FMFieldPtr iofield, void *data)
{
    if ((iofield->data_type == unsigned_type) || 
	(iofield->data_type == enumeration_type) || 
	(iofield->data_type == boolean_type)) {
	if (iofield->size == sizeof(char)) {
	    unsigned char tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(char));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (iofield->size == sizeof(short)) {
	    unsigned short tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(short));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(short));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (iofield->size == sizeof(int)) {
	    unsigned int tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(int));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(int));
	    return (MAX_UNSIGNED_TYPE) tmp;
	} else if (iofield->size == sizeof(long)) {
	    unsigned long tmp;
	    memcpy(&tmp, (char *) data + iofield->offset, sizeof(long));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(long));
	    return tmp;
	} else if (iofield->size == 2 * sizeof(long)) {
	    unsigned long tmp;
	    size_t low_bytes_offset = iofield->offset;
	    if (WORDS_BIGENDIAN) {
		if (!iofield->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    } else {
		if (iofield->byte_swap) {
		    low_bytes_offset += sizeof(long);
		}
	    }
	    memcpy(&tmp, (char *) data + low_bytes_offset, sizeof(long));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(long));
	    return tmp;
	} else {
	    if (!IO_shut_up && !get_long_warn) {
		fprintf(stderr, "Get Long failed!  Size problems.  File int size is %d.\n", iofield->size);
		get_long_warn++;
	    }
	    return 0;
	}
    } else if (iofield->data_type == integer_type) {
	MAX_INTEGER_TYPE tmp = get_big_int(iofield, data);
	return (MAX_UNSIGNED_TYPE) tmp;
    } else if (iofield->data_type == float_type) {
	MAX_FLOAT_TYPE tmp = get_big_float(iofield, data);
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

#define CONV(a,b) ((a<<16) + b)
static void
float_conversion(unsigned char*value, int size, FMfloat_format src_format,
		 FMfloat_format dest_format)
{
    int tmp;
    if (src_format == dest_format) return;
    if (src_format == ffs_reverse_float_formats[dest_format]) {
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
get_big_float(FMFieldPtr iofield, void *data)
{
    if (iofield->data_type == float_type) {
	if (iofield->size == sizeof(float)) {
	    float tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + iofield->offset), sizeof(float));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp, sizeof(float));
	    tmp2 = tmp;
	    return tmp2;
	} else if (iofield->size == sizeof(double)) {
	    double tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + iofield->offset), sizeof(double));
	    float_conversion((unsigned char *)&tmp, sizeof(double), 
			     (FMfloat_format)iofield->src_float_format,
			     (FMfloat_format)iofield->target_float_format);
	    tmp2 = tmp;
	    return tmp2;
#if SIZEOF_LONG_DOUBLE != 0
	} else if (iofield->size == sizeof(long double)) {
	    long double tmp;
	    MAX_FLOAT_TYPE tmp2;
	    memcpy(&tmp, ((char *) data + iofield->offset),
		   sizeof(long double));
	    if (iofield->byte_swap)
		byte_swap((char *) &tmp,
			  sizeof(long double));
	    tmp2 = tmp;
	    return tmp2;
#endif
	} else {
	    if (!IO_shut_up && !get_double_warn) {
		fprintf(stderr, "Get Double failed!  Size problems.  File double size is %d.\n", iofield->size);
		get_double_warn++;
	    }
	    return 0.0;
	}
    } else if (iofield->data_type == integer_type) {
	MAX_INTEGER_TYPE tmp = get_big_int(iofield, data);
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
    } else if (iofield->data_type == unsigned_type) {
	MAX_UNSIGNED_TYPE tmp = get_big_unsigned(iofield, data);
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

#if !defined(HAVE_DILL_H)
typedef void *dill_stream;
int ffs_putreg(void* s, int reg, int type) {return 0;}
int ffs_getreg(dill_stream s, int *reg_p, int type, int var_tmp){return 0;}
int ffs_localb(void*s, int size) {return 0;}
extern
 conv_routine
generate_conversion(IOConversionPtr conv, int src_alignment, int dest_alignment)
{
    return NULL;
}
#else

#undef max
#define max(x,y) (x<y?y:x)

static int
drisc_type(struct _FMgetFieldStruct *field)
{
    switch(field->data_type) {
    case integer_type:
    case unsigned_type:
	switch(field->size) {
	case 1:
	    return DILL_C;
	case 2:
	    return DILL_S;
	case 4:
	    return DILL_I;
	case 8:
	    return DILL_L;
	}
	return DILL_I;
    case float_type:
	if (field->size == SIZEOF_DOUBLE) {
	    return DILL_D;
	} else if (field->size == SIZEOF_FLOAT) {
	    return DILL_F;
	} else {
	    return DILL_I;
	}
    case char_type:
	return DILL_C;
    case string_type:
	return DILL_P;
    case boolean_type:
    case enumeration_type:
	return DILL_I;
    default:
	return DILL_I;
    }
}

static int
conv_required_alignment(dill_stream c, IOConversionPtr conv)
{
    if (conv->conv_count == 0) return 0;
    return conv->ioformat->body->alignment;
}

extern void
new_generate_conversion_code(dill_stream c, ConvStatus conv_status, IOConversionPtr conv, dill_reg *args, int assume_align, int register_args);

static int _ffs_conversion_generation = -1;
static int _generation_verbose = -1;

static inline int NO_SANITIZE_THREAD ffs_conversion_generation(void)
{
    return _ffs_conversion_generation;
}

static inline int NO_SANITIZE_THREAD generation_verbose(void)
{
    return _generation_verbose;
}

#define gen_fatal(str) do {fprintf(stderr, "%s\n", str); exit(1);} while (0)

int ffs_local(dill_stream s, int type)
{
#ifdef RAW
    return dill_local(s, type);
#else
    return dill_getvblock(s, dill_type_size(s, type));
#endif
}

int ffs_localb(void *vs, int size)
{
    dill_stream s = (dill_stream)vs;
#ifdef RAW
    return dill_localb(s, size);
#else
    return dill_getvblock(s, size);
#endif
}

int ffs_getreg(void *vs, int *reg_p, int type, int var_tmp)
{
    dill_stream s = (dill_stream)vs;
#ifdef RAW
    return dill_raw_getreg(s, reg_p, type, var_tmp);
#else
    int reg = dill_getreg(s, type);
    *reg_p = reg;
#endif
    return 1;
}

int ffs_putreg(void *vs, int reg, int type)
{
#ifdef RAW
    dill_stream s = (dill_stream)vs;
    return dill_raw_putreg(s, reg, type);
#endif
    return 1;
}


static void NO_SANITIZE_THREAD read_generation_environment_variables(void)
{
    if (_ffs_conversion_generation == -1) {
	char *gen_string = getenv("FFS_CONVERSION_GENERATION");
	int tmp = FFS_CONVERSION_GENERATION_DEFAULT;
	if (gen_string != NULL) {
#ifdef MODULE
	    tmp = strtol(gen_string, NULL, 10);
	    if ((tmp == LONG_MIN) || (tmp == LONG_MAX) || (tmp == 0)) {
#else
	    if (sscanf(gen_string, "%d", &tmp) != 1) {
#endif
		if (*gen_string == 0) {
		    /* empty string, just turn on generation */
		    tmp = 1;
		} else {
		    printf("Unable to parse FFS_CONVERSION_GENERATION environment variable \"%s\".\n", gen_string);
		}
	    }
	}
	_ffs_conversion_generation = tmp;
    }
    if (_debug_code_generation == -1)
	_debug_code_generation =
	    (getenv("FFS_CONVERSION_DEBUG") != NULL);
    if (_generation_verbose == -1)
	_generation_verbose =
	    (getenv("FFS_CONVERSION_VERBOSE") != NULL);
}


extern
conv_routine
generate_conversion(IOConversionPtr conv, int src_alignment, int dest_alignment)
{
    dill_stream c = NULL;
    void (*conversion_routine)();
    dill_reg args[6];
    dill_reg tmp_regs[10];
    char *format_name = conv->ioformat->body->format_name;
    int count = 0, register_args = 1;
    struct conv_status cs;
    read_generation_environment_variables();
    if (!ffs_conversion_generation())
	return NULL;
    if (generation_verbose()) {
	printf("For format %s ===================\n", format_name);
	dump_IOConversion(conv);
    }

    {
#ifdef RAW
	/* 
	 *  Determine whether or not some arguments should be left on the
	 *  stack.  The issue is that some architectures (x86) don't have
	 *  enough registers to follow the usual vcode convention of
	 *  immediately moving all arguments into registers for the
	 *  convenience of the generated subroutine.  If we let it do that,
	 *  we don't have enough temporary registers for use here.  So, we
	 *  first find out how many we have by doing getreg() until it fails
	 *  (or a max of 10 times).  If we don't have at least 8 registers
	 *  (four for args and four for temporary use), then we only want 
	 *  the first two arguments of the conversion code (the source and
	 *  destination address) to be kept in registers.  Trick vcode into
	 *  doing this by making all but two of the registers unavailable.
	 *  
	 *  IF WE DETERMINE WE HAVE ENOUGH REGISTERS, 'register_args' is set
	 *  to TRUE and the 'args' array contains the usual register
	 *  numbers.  If we don't have enough registers, 'register_args' is
	 *  false and args[2] and args[3] contain the stack offsets of the
	 *  third and fourth arguments.
	 */
	c = dill_create_raw_stream();
	count = 0;
	for (; count < sizeof(tmp_regs)/sizeof(tmp_regs[0]); count++) 
	    tmp_regs[count] = -1;
	count = 0;
	for (; count < sizeof(tmp_regs)/sizeof(tmp_regs[0]); count++) {
	    if (dill_raw_getreg(c, &tmp_regs[count], DILL_I, DILL_VAR) == 0) {
		break;
	    }
	}
	if (count <= 8) {
	    int i;
	    register_args = 0;
	    for (i= 2; i < count; i++) {  /* Make all but 2 unavail*/
		dill_raw_unavailreg(c, DILL_I, tmp_regs[i]);
	    }
	}
#else
	c = dill_create_stream();
	register_args = 1;
#endif
    }
    if (register_args) {
	/* Normal, lots of registers, case */

	dill_start_proc(c, "convert", DILL_I, "%p%p%p");
	args[0] = dill_param_reg(c, 0);
	args[1] = dill_param_reg(c, 1);
	args[2] = dill_param_reg(c, 2);

    } else {
	/* very few registers case */

	dill_parameter_type dr_params[3];	/* drisc param info */
	int i;
	for (i=0; i < 3; i++) {
	    dill_param_alloc(c, i, DILL_P, (dill_reg*)&args[i]);
	    dill_param_struct_alloc(c, i, DILL_P, &(dr_params[i]));
	}
	dill_start_simple_proc(c, "convert", DILL_I);

	/* store argument stack offsets in args[2]*/
	args[2] = dr_params[2].offset;
	for (count=2; count<sizeof(tmp_regs)/sizeof(tmp_regs[0]); count++) {
	    if (tmp_regs[count] != -1) {
		/* make the other registers available again */
		dill_raw_availreg(c, DILL_I, tmp_regs[count]);
	    }
	}
    }
    if (debug_code_generation()) {
	if (register_args) {
	    dill_reg src_pointer_base, dest_pointer_base, dest_offset_adjust;	    dill_scallv(c, (void*)printf, "printf", "%P%P%p%p%p",
		     "convert for %s called with src= %p, dest %p, rt_conv_status =%p\n",
		     format_name, args[0], args[1], args[2]);
	    ffs_getreg(c, &src_pointer_base, DILL_P, DILL_VAR);
	    dill_ldpi(c, src_pointer_base, args[2], FMOffset(RTConvStatus,src_pointer_base));
	    ffs_getreg(c, &dest_pointer_base, DILL_P, DILL_VAR);
	    dill_ldpi(c, dest_pointer_base, args[2], FMOffset(RTConvStatus,dest_pointer_base));
	    ffs_getreg(c, &dest_offset_adjust, DILL_I, DILL_VAR);
	    dill_ldpi(c, dest_offset_adjust, args[2], FMOffset(RTConvStatus,dest_offset_adjust));
	    dill_scallv(c, (void*)printf, "printf", "%P%p%p%p",
		     "rt_conv_status is src_pointer_base= %p, dest_pointer_base=%p, dest_offset_adjust =%zd\n",
		     src_pointer_base, dest_pointer_base, dest_offset_adjust);
	} else {
#ifdef HAVE_DILL_H	    
	    dill_reg v_at;
	    if (dill_raw_getreg(c, &v_at, DILL_I, DILL_TEMP) == 0) {
		gen_fatal("Out of regs 1\n");
	    }
#endif
	    dill_scallv(c, (void*)printf, "printf", "%P%P%p%p",
		     "convert for %s called with src= %p, dest %p\n",
		     format_name, args[0], args[1]);
	    dill_ldpi(c, v_at, dill_lp(c), args[2]);
	    dill_scallv(c, (void*)printf, "printf", "%P%p",
		     "               rt_conv_status %p\n",
		     v_at);
#ifdef HAVE_DILL_H	    
	    dill_raw_putreg(c, v_at, DILL_I);
#endif
	}
    }
    conv->required_alignment = conv_required_alignment(c, conv);
    if (register_args) {
	dill_reg tmp;
	int mask;
	if (ffs_getreg(c, &tmp, DILL_I, DILL_VAR) == 0) {
	    printf("out of regs for mod\n");
	}
	switch(conv->required_alignment) {
	case 2:
	    mask = 0x1;
	    break;
	case 4:
	    mask = 0x3;
	    break;
	case 8:
	    mask = 0x7;
	    break;
	default:
	    mask = 0;
	}
	if (mask != 0) {
	    int zero_target = dill_alloc_label(c, NULL);
	    dill_anduli(c, tmp, args[0], mask);
	    dill_beqli(c, tmp, 0, zero_target);
	    dill_scallv(c, (void*)printf, "printf", "%P%P%p%I",
		    "convert for %s called with bad align src= %p, align is %d\n",
		    format_name, args[0], conv->required_alignment);
	    dill_mark_label(c, zero_target);
	}
    }
    
    cs.src_pointer_base = 0;
    cs.dest_pointer_base = 0;
    cs.src_offset_adjust = -(ssize_t)conv->string_offset_size;
    cs.dest_offset_adjust = -(ssize_t)conv->string_offset_size;
    cs.cur_offset = 0;
    cs.control_value = NULL;
    cs.target_pointer_size = conv->target_pointer_size;
    cs.src_pointer_size = conv->ioformat->body->pointer_size;
    cs.register_args = register_args;
    cs.global_conv = conv;
    conv->conv_pkg = NULL;
    new_generate_conversion_code(c, &cs, conv, args, src_alignment, register_args);
    dill_retp(c, args[2]);
    if (conv->conv_pkg == (char*)-1) {
	int pkg_len;
	conv->conv_pkg = dill_finalize_package(c, &pkg_len);
	conv->free_data = conv->conv_pkg;
	conv->free_func = (void(*)(void*))&free;
	dill_free_stream(c);
	return NULL;
    } else {
	dill_exec_handle conversion_handle;

	conversion_handle = dill_finalize(c);
	dill_free_handle(conversion_handle);
	conversion_handle = dill_get_handle(c);
	conv->free_data = conversion_handle;
	conv->free_func = (void(*)(void*))&dill_free_handle;
	conversion_routine = (void(*)()) dill_get_fp(conversion_handle);
    }
    if (generation_verbose()) {
	dill_dump(c);
    }
    dill_free_stream(c);
    return (conv_routine) conversion_routine;
}
/*#define REG_DEBUG(x) printf x ;*/
#define REG_DEBUG(x)

static void
gen_mem_float_conv(dill_stream c, struct _FMgetFieldStruct src, int src_addr, 
		   size_t src_offset, int assume_align,
		   int dest_reg, size_t dest_offset,
		   int dest_size, int dst_aligned)
{
    FMfloat_format src_format = (FMfloat_format) src.src_float_format;
    FMfloat_format dst_format = (FMfloat_format) src.target_float_format;
    int src_drisc_type = drisc_type(&src);

    if (src_format == dst_format) {
	dill_reg tmp;
	switch (src_drisc_type) {
	case DILL_D:
	    if (assume_align >= TYPE_ALIGN(c, DILL_D)) {
		ffs_getreg(c, &tmp, DILL_D, DILL_TEMP);
		dill_lddi(c, tmp, src_addr, src_offset);
		dill_stdi(c, tmp, dest_reg, dest_offset);
		ffs_putreg(c, tmp, DILL_D);
		return;
	    }
	    break;
	case DILL_F:
	    if (assume_align >= TYPE_ALIGN(c, DILL_D)) {
		ffs_getreg(c, &tmp, DILL_F, DILL_TEMP);
		dill_ldfi(c, tmp, src_addr, src_offset);
		dill_stfi(c, tmp, dest_reg, dest_offset);
		ffs_putreg(c, tmp, DILL_F);
		return;
	    }
	    break;
	}
	gen_memcpy(c, src_addr, src_offset, dest_reg, dest_offset, 0,
		   dest_size);
	return;
    }
    if (src_format == ffs_reverse_float_formats[dst_format]) {
	switch(dest_size) {
	case sizeof(short): {
	    dill_reg tmp;
	    ffs_getreg(c, &tmp, DILL_S, DILL_TEMP);
	    dill_ldbssi(c, tmp, src_addr, src_offset);
	    dill_stsi(c, tmp, dest_reg, dest_offset);
	    ffs_putreg(c, tmp, DILL_S);
	    break;
	}
	case sizeof(int): {
	    dill_reg tmp;
	    ffs_getreg(c, &tmp, DILL_I, DILL_TEMP);
	    dill_ldbsii(c, tmp, src_addr, src_offset);
	    dill_stii(c, tmp, dest_reg, dest_offset);
	    ffs_putreg(c, tmp, DILL_I);
	    break;
	}
#if SIZEOF_LONG == 8
	case sizeof(long): 
	    if (((src_offset & 0x7) == 0) && (assume_align >= sizeof(long))) {
		dill_reg tmp;
		ffs_getreg(c, &tmp, DILL_L, DILL_TEMP);
		dill_ldbsli(c, tmp, src_addr, src_offset);
		dill_stli(c, tmp, dest_reg, dest_offset);
		ffs_putreg(c, tmp, DILL_L);
		break;
	    } else {
		dill_reg tmp, tmp2;
		int i;
		ffs_getreg(c, &tmp, DILL_I, DILL_TEMP);
		ffs_getreg(c, &tmp2, DILL_I, DILL_TEMP);
	REG_DEBUG(("Getting reg %d for float conv\n", tmp));
	REG_DEBUG(("Getting reg %d for float conv\n", tmp2));
		for (i = 0; i < (dest_size >> 1); i += sizeof(int)) {
		    int near_offset = i*sizeof(int);
		    int far_offset = dest_size - (i+1)*sizeof(int);
		    dill_ldbsii(c, tmp, src_addr, src_offset + near_offset);
		    dill_ldbsii(c, tmp2, src_addr, src_offset + far_offset);
		    dill_stii(c, tmp, dest_reg, dest_offset + far_offset);
		    dill_stii(c, tmp2, dest_reg, dest_offset + near_offset);
		}
	REG_DEBUG(("putting reg %d for float conv\n", tmp));
	REG_DEBUG(("putting reg %d for float conv\n", tmp2));
		ffs_putreg(c, tmp, DILL_I);
		ffs_putreg(c, tmp2, DILL_I);
		break;
	    }
#endif
	default: 
	{
	    dill_reg tmp, tmp2;
	    int i;
	    ffs_getreg(c, &tmp, DILL_L, DILL_TEMP);
	    ffs_getreg(c, &tmp2, DILL_L, DILL_TEMP);
	REG_DEBUG(("Getting reg %d for float conv\n", tmp));
	REG_DEBUG(("Getting reg %d for float conv\n", tmp2));
	    for (i = 0; i < (dest_size >> 1); i += sizeof(long)) {
		int near_offset = i*sizeof(int);
		int far_offset = dest_size - (i+1)*sizeof(int);
		dill_ldbsli(c, tmp, src_addr, src_offset + near_offset);
		dill_ldbsli(c, tmp2, src_addr, src_offset + far_offset);
		dill_stli(c, tmp, dest_reg, dest_offset + far_offset);
		dill_stli(c, tmp2, dest_reg, dest_offset + near_offset);
	    }
	REG_DEBUG(("Putting reg %d for float conv\n", tmp));
	REG_DEBUG(("Putting reg %d for float conv\n", tmp2));
	    ffs_putreg(c, tmp, DILL_L);
	    ffs_putreg(c, tmp2, DILL_L);
	}
	}
	return;
    }
    switch (CONV(src_format, dst_format)) {
    case CONV(Format_IEEE_754_bigendian, Format_IEEE_754_mixedendian):
    case CONV(Format_IEEE_754_mixedendian, Format_IEEE_754_bigendian):
	{
	    /* byte swap in place */
	    dill_reg tmp;
	    int i;
	    ffs_getreg(c, &tmp, DILL_I, DILL_TEMP);
	    for (i = 0; i < dest_size; i += sizeof(int)) {
		dill_ldbsii(c, tmp, src_addr, src_offset + i);
		dill_stii(c, tmp, dest_reg, dest_offset + i);
	    }
	    ffs_putreg(c, tmp, DILL_I);
	}
	break;
    case CONV(Format_IEEE_754_littleendian, Format_IEEE_754_mixedendian):
    case CONV(Format_IEEE_754_mixedendian, Format_IEEE_754_littleendian):
	{
	    /* swap words, no byteswapping */
	    dill_reg tmp, tmp2;
	    int i;
	    ffs_getreg(c, &tmp, DILL_I, DILL_TEMP);
	    ffs_getreg(c, &tmp2, DILL_I, DILL_TEMP);
	    for (i = 0; i < (dest_size >> 1); i += sizeof(int)) {
		int near_offset = i;
		int far_offset = dest_size - (i+sizeof(int));
		dill_ldii(c, tmp, src_addr, src_offset + near_offset);
		dill_ldii(c, tmp2, src_addr, src_offset + far_offset);
		dill_stii(c, tmp, dest_reg, dest_offset + far_offset);
		dill_stii(c, tmp2, dest_reg, dest_offset + near_offset);
	    }
	    ffs_putreg(c, tmp, DILL_I);
	    ffs_putreg(c, tmp2, DILL_I);
	}
	break;
    default:
	printf("unanticipated float conversion \n");
    }
}
	    
static void
gen_simple_field_conv(dill_stream c, struct _FMgetFieldStruct tmp_spec, int assume_align,
		      dill_reg src_addr, size_t src_offset, int dest_size,
		      FMdata_type dest_type, dill_reg dest_addr, size_t dest_offset)
{
    /* simple conversion */
    iogen_oprnd src_oprnd;
    struct _FMgetFieldStruct dest_spec;
    int src_drisc_type;
    int dst_drisc_type;
    int src_required_align;
    int dst_required_align;
    int src_is_aligned = 1;
    int dst_is_aligned = 1;

    dest_spec.data_type = dest_type;
    dest_spec.size = dest_size;
    src_drisc_type = drisc_type(&tmp_spec);
    dst_drisc_type = drisc_type(&dest_spec);
    src_required_align = TYPE_ALIGN(c, src_drisc_type);
    dst_required_align = TYPE_ALIGN(c, dst_drisc_type);
    if ((assume_align < src_required_align) || 
	((src_offset % src_required_align) != 0)) {
	src_is_aligned = 0;
    }
    if ((assume_align < dst_required_align) || 
	((dest_offset % dst_required_align) != 0)) {
	dst_is_aligned = 0;
    }

    if (tmp_spec.data_type != float_type) {
	assert(dest_type != string_type);

	src_oprnd = gen_operand(src_addr, src_offset, tmp_spec.size,
				tmp_spec.data_type,
				src_is_aligned,
				tmp_spec.byte_swap);

	if (src_oprnd.address) {
	    gen_load(c, &src_oprnd);
	}
	if (src_oprnd.data_type != dest_type) {
	    iogen_oprnd tmp_oprnd;
	    tmp_oprnd = gen_type_conversion(c, src_oprnd, dest_type);
	    free_oprnd(c, src_oprnd);
	    src_oprnd = tmp_oprnd;
	}
	if (src_oprnd.size != dest_size) {
	    iogen_oprnd tmp_oprnd;
	    tmp_oprnd = gen_size_conversion(c, src_oprnd, dest_size);
	    free_oprnd(c, src_oprnd);
	    src_oprnd = tmp_oprnd;
	}
	gen_store(c, src_oprnd, dest_addr, dest_offset,
		  dest_size, dest_type, dst_is_aligned);
	free_oprnd(c, src_oprnd);
    } else {
	if (tmp_spec.size == dest_size) {
	    gen_mem_float_conv(c, tmp_spec, src_addr, src_offset, assume_align,
			       dest_addr, dest_offset, dest_size, 
			       dst_is_aligned);
	    return;
	}
	if ((dst_drisc_type != DILL_I) && (src_drisc_type != DILL_I)) {
	    /* both float sizes supported */
	    int float_conv_offset = ffs_local(c, dst_drisc_type);
	    dill_reg float_reg;
	    gen_mem_float_conv(c, tmp_spec, src_addr, src_offset, assume_align,
			       dill_lp(c), float_conv_offset, dest_size, 1);
	    ffs_getreg(c, &float_reg, DILL_D, DILL_TEMP);
	    switch(dst_drisc_type) {
	    case DILL_D:
		dill_lddi(c, float_reg, dill_lp(c), float_conv_offset);
		dill_cvd2f(c, float_reg, float_reg);
		dill_stfi(c, float_reg, dest_addr, dest_offset);
		break;
	    case DILL_F:
		dill_ldfi(c, float_reg, dill_lp(c), float_conv_offset);
		dill_cvf2d(c, float_reg, float_reg);
		dill_stdi(c, float_reg, dest_addr, dest_offset);
		break;	    
	    }
	    ffs_putreg(c, float_reg, DILL_D);
	    return;
	}
	printf("must do call to conversion subroutine\n");
    }
}

static void
gen_convert_address_field(dill_stream c, struct _FMgetFieldStruct tmp_spec, int assume_align,
			  dill_reg src_addr, size_t src_offset, size_t dest_size, dill_reg dest_addr,
			  size_t dest_offset, size_t string_offset_size, dill_reg rt_conv_status,
			  size_t base_size_delta, dill_reg *string_src_reg, dill_reg *string_dest_reg,
			  int register_args, int null_target)
{
    iogen_oprnd src_oprnd;
    int src_drisc_type;
    int src_required_align;
    int src_is_aligned = 1;
    dill_reg tmp_dest_reg;

    src_drisc_type = drisc_type(&tmp_spec);
    src_required_align = TYPE_ALIGN(c, src_drisc_type);
    if ((assume_align < src_required_align) || 
	((src_offset % src_required_align) != 0)) {
	src_is_aligned = 0;
    }


    /* handle addresses */
    src_oprnd = gen_fetch(c, src_addr, src_offset, tmp_spec.size,
			  integer_type, src_is_aligned,
			  tmp_spec.byte_swap);

    *string_dest_reg = src_oprnd.vc_reg;
    /* src_oprnd now holds the offset value */
    if (dest_size >= sizeof(char *)) {

	if (src_oprnd.size != dest_size) {
	    /* make it the right size to operate on */
	    iogen_oprnd tmp_oprnd;
	    tmp_oprnd = gen_size_conversion(c, src_oprnd, sizeof(long));
	    free_oprnd(c, src_oprnd);
	    src_oprnd = tmp_oprnd;
	    *string_dest_reg = src_oprnd.vc_reg;
	}
#ifdef VERBOSE
	dill_scallv(c, (void*)printf, "printf", "%P%i%p%I",
		    "Fetched msgptr %d from address %p offset %d\n", src_oprnd.vc_reg,
		    src_addr, src_offset);
#endif

	/* generate : if it's zero, leave it zero  branch away */
	dill_beqli(c, src_oprnd.vc_reg, 0, null_target);

	/* else, sub the string_offset_size */
	dill_addli(c, src_oprnd.vc_reg, src_oprnd.vc_reg,
		string_offset_size);
 
	/* Moving to here to more efficiently use registers */   
	if (!ffs_getreg(c, string_src_reg, DILL_P, DILL_TEMP))
	  gen_fatal("gen field convert out of registers C\n");
	REG_DEBUG(("Getting reg %d for string src reg\n", *string_src_reg));
	if (!ffs_getreg(c, &tmp_dest_reg, DILL_P, DILL_TEMP))
	  gen_fatal("gen field convert out of registers D\n");
	REG_DEBUG(("Getting reg %d for string dest reg\n", tmp_dest_reg));

	/* calculate the address of this in the source */
	if (register_args) {
	    dill_ldpi(c, *string_src_reg, rt_conv_status, FMOffset(RTConvStatus,src_pointer_base));
	    dill_addl(c, *string_src_reg, src_oprnd.vc_reg, *string_src_reg);
	} else {
	    dill_ldpi(c, *string_src_reg, dill_lp(c), rt_conv_status);
	    dill_ldpi(c, *string_src_reg, *string_src_reg, FMOffset(RTConvStatus,src_pointer_base));
	    dill_addl(c, *string_src_reg, src_oprnd.vc_reg, *string_src_reg);
	}
	    
	/* and the address in the destination */
	if (register_args) {
	    dill_ldpi(c, tmp_dest_reg, rt_conv_status, FMOffset(RTConvStatus, dest_pointer_base));
	    dill_addl(c, src_oprnd.vc_reg, src_oprnd.vc_reg, tmp_dest_reg);
	} else {
	    dill_ldpi(c, tmp_dest_reg, dill_lp(c), rt_conv_status);
	    dill_ldpi(c, tmp_dest_reg, tmp_dest_reg, FMOffset(RTConvStatus, dest_pointer_base));
	    dill_addl(c, src_oprnd.vc_reg, src_oprnd.vc_reg, tmp_dest_reg);
	}
/*	dill_mark_label(c, null_target);*/

	dill_movp(c, *string_dest_reg, src_oprnd.vc_reg);
	if (dest_size > sizeof(char *)) {
	    iogen_oprnd tmp_oprnd;
	    printf("Doing gen size conversion\n");
	    tmp_oprnd = gen_size_conversion(c, src_oprnd,
					    (int)dest_size);
	    free_oprnd(c, src_oprnd);
	    src_oprnd = tmp_oprnd;
	}
	REG_DEBUG(("Regs %d and %d are src and dest \n",
		   _vrr(*string_src_reg), _vrr(*string_dest_reg)));
	free_oprnd(c, src_oprnd);
    } else {
	assert(FALSE);
    }
}

static void
generate_convert_field(dill_stream c, ConvStatus conv_status, dill_reg src_addr,
		       size_t src_offset, dill_reg dest_addr, size_t dest_offset,
		       dill_reg rt_conv_status, IOconvFieldStruct *conv, FMTypeDesc *type_desc,
		       int data_already_copied)
{
    switch(type_desc->type) {
    case FMType_pointer: {
	dill_reg actual_src_reg, actual_dest_reg;
	struct _FMgetFieldStruct tmp_spec = conv->src_field;
	int null_target = dill_alloc_label(c, NULL);
	tmp_spec.offset = 0;
	tmp_spec.size = conv_status->src_pointer_size;
	
	gen_convert_address_field(c, tmp_spec, 0, src_addr, src_offset,
				  conv_status->target_pointer_size, 
				  dest_addr, dest_offset, 
				  conv_status->src_offset_adjust, 
				  rt_conv_status,
				  decode_size_delta(conv_status, conv, type_desc),
				  &actual_src_reg, &actual_dest_reg, 
				  conv_status->register_args, null_target);
	{
#ifdef RAW
	    int src_storage = ffs_local(c, DILL_P);
	    int dest_storage = ffs_local(c, DILL_P);
	    /* save values of src_addr and dest_addr */
	    dill_stpi(c, src_addr, dill_lp(c), src_storage);
	    dill_stpi(c, dest_addr, dill_lp(c), dest_storage);
#else
	    dill_reg tmp_src_reg, tmp_dest_reg;
	    ffs_getreg(c, &tmp_dest_reg, DILL_P, DILL_VAR);
	    ffs_getreg(c, &tmp_src_reg, DILL_P, DILL_VAR);
#endif
	    dill_movp(c, tmp_src_reg, actual_src_reg);
	    dill_movp(c, tmp_dest_reg, actual_dest_reg);
	    REG_DEBUG(("Putting %d and %d for new src & dest\n", 
		       actual_src_reg, actual_dest_reg));
	    ffs_putreg(c, actual_src_reg, DILL_P);
	    ffs_putreg(c, actual_dest_reg, DILL_P);
	    generate_convert_field(c, conv_status, tmp_src_reg, 0, 
				   tmp_dest_reg, 0, rt_conv_status,
				   conv, type_desc->next, conv_status->global_conv->conversion_type == direct_to_mem);
	    
#ifdef RAW
	    dill_ldpi(c, src_addr, dill_lp(c), src_storage);
	    dill_ldpi(c, dest_addr, dill_lp(c), dest_storage);
#endif
	}
	dill_mark_label(c, null_target);
	
#ifdef VERBOSE
	dill_scallv(c, (void*)printf, "printf", "%P%p",
		    "storing POINTER value %p\n",actual_dest_reg);
#endif
	dill_stpi(c, actual_dest_reg, dest_addr, dest_offset);

	break;
    }
    case FMType_string:{
	dill_reg actual_src_reg, actual_dest_reg;
	struct _FMgetFieldStruct tmp_spec = conv->src_field;
	tmp_spec.offset = 0;
	int null_target = dill_alloc_label(c, NULL);
	gen_convert_address_field(c, tmp_spec, 0, src_addr, src_offset,
				  conv->dest_size, dest_addr, dest_offset, 
				  conv_status->src_offset_adjust, 
				  rt_conv_status,
				  decode_size_delta(conv_status, conv, type_desc),
				  &actual_src_reg, &actual_dest_reg, 
				  conv_status->register_args, null_target);
	/* generate strcpy */
	if (!data_already_copied) {
#ifdef VERBOSE
	    dill_scallv(c, (void*)printf, "printf", "%P%p%p",
			"Calling Strcpy with args %p, %p\n",actual_dest_reg, actual_src_reg);

#endif
	    dill_scallv(c, (void*)strcpy, "strcpy", "%p%p", actual_dest_reg, actual_src_reg);

	}
	dill_mark_label(c, null_target);
	
	dill_stpi(c, actual_dest_reg, dest_addr, dest_offset);
	REG_DEBUG(("Putting %d and %d for new src & dest\n", 
		   actual_src_reg, actual_dest_reg));
	ffs_putreg(c, actual_src_reg, DILL_P);
	ffs_putreg(c, actual_dest_reg, DILL_P);
	break;
    }
    case FMType_subformat: {
	FFSTypeHandle subformat = conv->subconversion->ioformat;
	char *name = FFSTypeHandle_name(subformat);
	if (conv->subconversion->conv_func == NULL) {
	    /* we're not linking to an address that's valid, 
	       must fill it in later */
	    conv_status->global_conv->conv_pkg = (char *) -1;
	}
	    
	if (conv_status->register_args) {
	    dill_reg new_src, new_dest;
	    if (!ffs_getreg(c, &new_src, DILL_P, DILL_TEMP) ||
		!ffs_getreg(c, &new_dest, DILL_P, DILL_TEMP))
		gen_fatal("temp vals in subcall\n");
	    REG_DEBUG(("Getting %d and %d for new src & dest\n", 
		       new_src, new_dest));
	    dill_addpi(c, new_src, src_addr, src_offset);
	    dill_addpi(c, new_dest, dest_addr, dest_offset);
	    (void) dill_scallp(c, (void*)conv->subconversion->conv_func, name, "%p%p%p", new_src,
			      new_dest, rt_conv_status);
	    REG_DEBUG(("Putting %d and %d for new src & dest\n", 
		       new_src, new_dest));
	    if (debug_code_generation()) {
/*	        VCALL2V(printf, "%P%p",
	    "After subroutine call, new src_string_base is %p\n", src_string_base);*/
	    }
	    ffs_putreg(c, new_src, DILL_P);
	    ffs_putreg(c, new_dest, DILL_P);
	} else {
	    dill_reg reg_rt_conv_status;
	    dill_reg new_src, new_dest;
	    if (!ffs_getreg(c, &new_src, DILL_P, DILL_TEMP) ||
		!ffs_getreg(c, &new_dest, DILL_P, DILL_TEMP))
		gen_fatal("temp vals in subcall\n");
#ifdef RAW
	    int src_storage = ffs_local(c, DILL_P);
	    int dest_storage = ffs_local(c, DILL_P);

	    /* save values of src_addr and dest_addr */
	    dill_stpi(c, src_addr, dill_lp(c), src_storage);
	    dill_stpi(c, dest_addr, dill_lp(c), dest_storage);
#endif

	    if (!ffs_getreg(c, &reg_rt_conv_status, DILL_P, DILL_TEMP))
		gen_fatal("temp string vals in subcall\n");
	    REG_DEBUG(("Getting %d and %d for reg_rt_conv_status\n", reg_rt_conv_status));
	    dill_addpi(c, new_src, src_addr, src_offset);
	    dill_addpi(c, new_dest, dest_addr, dest_offset);
	    dill_ldpi(c, reg_rt_conv_status, dill_lp(c), rt_conv_status);
	    dill_scallp(c, (void*)conv->subconversion->conv_func, name, "%p%p%p", new_src,
			new_dest, reg_rt_conv_status);
	    REG_DEBUG(("Putting %d reg_rt_conv_status\n", reg_rt_conv_status));
	    ffs_putreg(c, reg_rt_conv_status, DILL_P);
	    /* restore values of src_addr and dest_addr */
#ifdef RAW
	    dill_ldpi(c, src_addr, dill_lp(c), src_storage);
	    dill_ldpi(c, dest_addr, dill_lp(c), dest_storage);
#endif
	}
	break;
    }
    case FMType_simple: {
	struct _FMgetFieldStruct tmp_spec = conv->src_field;
	tmp_spec.offset = 0;
	gen_simple_field_conv(c, tmp_spec, 
			      conv_status->global_conv->required_alignment, 
			      src_addr, src_offset, 
			      conv->dest_size, tmp_spec.data_type, 
			      dest_addr, dest_offset);
	break;
    }
    case FMType_array:
	if (conv->rc_swap == no_row_column_swap) {
	    int static_elements = 1;
	    struct _FMgetFieldStruct tmp_spec = conv->src_field;
	    FMTypeDesc *next = type_desc;
	    while (next->type == FMType_array) {
		if (next->static_size != 0) {
		    static_elements *= next->static_size;
		}
		next = next->next;
	    }
	
	    dill_reg loop_var;
	    int loop_var_type;

	    int loop_head, loop_end;
	    if (data_already_copied &&
		!conv->src_field.byte_swap &&
		(conv->src_field.src_float_format == conv->src_field.target_float_format) &&
		(conv->src_field.size == conv->dest_size) &&
		(conv->subconversion == NULL) &&
		(next->type != FMType_pointer) &&
		(next->type != FMType_string) &&
		((conv->src_field.data_type != string_type))) {
		/* Nothing but data movement required */
		    break;
	    }
#ifdef RAW
	    int src_storage = ffs_local(c, DILL_P);
	    int dest_storage = ffs_local(c, DILL_P);
	    int loop_storage = 0;

	    /* save values of src_addr and dest_addr */
	    dill_stpi(c, src_addr, dill_lp(c), src_storage);
	    dill_stpi(c, dest_addr, dill_lp(c), dest_storage);
#else
	    {
		dill_reg tmp_src, tmp_dest;
		ffs_getreg(c, &tmp_src, DILL_P,DILL_TEMP);
		ffs_getreg(c, &tmp_dest, DILL_P,DILL_TEMP);
		dill_movp(c, tmp_src, src_addr);
		dill_movp(c, tmp_dest, dest_addr);
		src_addr = tmp_src;
		dest_addr = tmp_dest;
	    }
#endif

	    if (((conv->subconversion == NULL)) &&
		!debug_code_generation()) {
		if (!ffs_getreg(c, &loop_var, DILL_I, DILL_TEMP))
		    gen_fatal("gen field convert out of registers BB \n");
		loop_var_type = DILL_TEMP;
	    } else {
		/* may call a subconversion in here, use VARs */
		if (!ffs_getreg(c, &loop_var, DILL_I, DILL_VAR))
		    gen_fatal("gen field convert out of registers CC\n");
		loop_var_type = DILL_VAR;
	    }
	    REG_DEBUG(("Getting1 %d as loop_var\n", _vrr(loop_var)));
	    dill_addpi(c, src_addr, src_addr, src_offset);
	    dill_addpi(c, dest_addr, dest_addr, dest_offset);

	    /* gen conversion loop */
	    loop_head = dill_alloc_label(c, NULL);
	    loop_end = dill_alloc_label(c, NULL);
	    dill_seti(c, loop_var, static_elements);
	    next = type_desc;
	    while (next->type == FMType_array) {
		if (next->static_size == 0) {
		    dill_reg addr_reg = (dill_reg)(intptr_t)conv_status->control_value;
		    dill_reg val;
		    int field = next->control_field_index;
		    ffs_getreg(c, &val, DILL_I, DILL_TEMP);
		    dill_ldii(c, val, addr_reg, field*sizeof(int));
		    dill_muli(c, loop_var, loop_var, val);
		}
		next = next->next;
	    }
	    if (!data_already_copied) {
		ssize_t base_delta = decode_size_delta(conv_status, conv, next);
		if (conv_status->global_conv->conversion_type == copy_dynamic_portion) {
		    int base_size = item_size(conv_status, conv, next);
		    dill_reg size;
		    ffs_getreg(c, &size, DILL_I, DILL_TEMP);
		    dill_mulii(c, size, loop_var, base_size);
#ifdef VERBOSE
		    dill_scallv(c, (void*)printf, "printf", "%P%p%p",
				"Calling Memcpy with args %p, %p, %d\n",dest_addr, src_addr, size);
#endif
		    dill_scallv(c, (void*)memcpy, "memcpy", "%p%p%i",
				dest_addr, src_addr, size);
		}
		if (base_delta != 0) {
		    dill_reg dest_src_ptr, delta;
		    ffs_getreg(c, &delta, DILL_I, DILL_TEMP);
		    ffs_getreg(c, &dest_src_ptr, DILL_I, DILL_TEMP);
		    dill_mulii(c, delta, loop_var, base_delta);
#ifdef VERBOSE
		    dill_scallv(c, (void*)printf, "printf", "%P%p",
				"Adjusting dest_pointer_base by %p\n", delta);
#endif
		    dill_ldpi(c, dest_src_ptr, rt_conv_status, FMOffset(RTConvStatus,dest_pointer_base));
		    dill_addp(c, dest_src_ptr, dest_src_ptr, delta);
		    dill_stpi(c, dest_src_ptr, rt_conv_status, FMOffset(RTConvStatus,dest_pointer_base));
		}
	    }

	    if (debug_code_generation()) {
		dill_scallv(c, (void*)printf, "printf", "%P%P%p",
			    "format %s, field Initial loopvar = %x\n", conv_status->global_conv->ioformat->body->format_name, loop_var);
	    }
	    
	    dill_bleii(c, loop_var, 0, loop_end);
	    dill_mark_label(c, loop_head);
#if defined(NOT) & defined(RAW)
	    if (!register_args) {
		/* store away loop var and free the reg */
		loop_storage = ffs_local(c, DILL_I);
		dill_stii(c, loop_var, dill_lp(c), loop_storage);
		REG_DEBUG(("Putting %d as loop_var\n", _vrr(loop_var)));
		ffs_putreg(c, loop_var, DILL_I);
	    }
#endif
	    generate_convert_field(c, conv_status, src_addr, 0, dest_addr, 0, 
				   rt_conv_status, conv, next, 
				   data_already_copied);
	    
	    (void)loop_var_type;  /* avoid warning */
#if defined(NOT) & defined(RAW)
	    /* generate end of loop */
	    if (!register_args) {
		/* store away loop var and free the reg */
		ffs_getreg(c, &loop_var, DILL_I, loop_var_type);
		REG_DEBUG(("Getting %d as loop_var\n", _vrr(loop_var)));
		dill_ldii(c, loop_var, dill_lp(c), loop_storage);
	    }
#endif
	    dill_subii(c, loop_var, loop_var, 1);
	    int array_of_pointers = 0;
	    if (conv->iovar->type_desc.type == FMType_pointer) {
		array_of_pointers = (conv->iovar->type_desc.next->next->type == FMType_pointer);
	    }
	    if (array_of_pointers) {
		dill_addpi(c, src_addr, src_addr, conv_status->src_pointer_size);
		dill_addpi(c, dest_addr, dest_addr, conv_status->target_pointer_size);
	    } else {
		dill_addpi(c, src_addr, src_addr, tmp_spec.size);
		dill_addpi(c, dest_addr, dest_addr, conv->dest_size);
	    }		
	    if (debug_code_generation()) {
		dill_scallv(c, (void*)printf, "printf", "%P%p%p%p",
			    "loopvar = %x, src %x, dest %x\n", loop_var,
			    src_addr, dest_addr);
	    }
	    dill_bgtii(c, loop_var, 0, loop_head);
	    dill_mark_label(c, loop_end);
	    ffs_putreg(c, loop_var, DILL_I);
	    REG_DEBUG(("Putting %d as loop_var\n", _vrr(loop_var)));
	
	    /* restore values of src_addr and dest_addr */
#ifdef RAW
	    dill_ldpi(c, src_addr, dill_lp(c), src_storage);
	    dill_ldpi(c, dest_addr, dill_lp(c), dest_storage);
#endif
	} else {
	    /* generate a call to transpose */
	    struct _FMgetFieldStruct tmp_spec = conv->src_field;
	    dill_reg dimens, dimen_reg, spec, spec_reg, tmp;
	    int source_column_major = 
		(conv->rc_swap == swap_source_column_major);
	    int dimen_count = conv->iovar->dimen_count;
	    FMdata_type dest_type = conv->src_field.data_type;

	    ffs_getreg(c, &dimen_reg, DILL_P, DILL_TEMP);
	    ffs_getreg(c, &tmp, DILL_I, DILL_TEMP);
	    ffs_getreg(c, &spec_reg, DILL_P, DILL_TEMP);
	    dimens = ffs_localb(c, dimen_count * sizeof(int));
	    spec = ffs_localb(c, sizeof(struct _FMgetFieldStruct));
	    dill_virtual_lea(c, dimen_reg, dimens);
	    dill_virtual_lea(c, spec_reg, spec);
	    FMTypeDesc *next = type_desc;
	    int i = 0;
	    while (next->type == FMType_array) {
		if (next->static_size != 0) {
		    dill_seti(c, tmp, next->static_size);
		} else {
		    dill_reg addr_reg = (dill_reg)(intptr_t)conv_status->control_value;
		    int field = next->control_field_index;
		    dill_ldii(c, tmp, addr_reg, field*sizeof(int));
		}
		dill_stii(c, tmp, dimen_reg, i * sizeof(int));
		i++;
		next = next->next;
	    }
	    dill_seti(c, tmp, 0);
	    dill_stii(c, tmp, dimen_reg, dimen_count * sizeof(int));
	    dill_stii(c, tmp, spec_reg, FMOffset(struct _IOgetFieldStruct *, offset));
	    for (i=4 ; i < sizeof(struct _IOgetFieldStruct); i+= 4) {
		dill_seti(c, tmp, *((int*)((char*)&tmp_spec + i)));
		dill_stii(c, tmp, spec_reg, i);
	    }
	    {
		dill_reg tmp_src, tmp_dest;
		ffs_getreg(c, &tmp_src, DILL_P,DILL_TEMP);
		ffs_getreg(c, &tmp_dest, DILL_P,DILL_TEMP);
		dill_addpi(c, tmp_src, src_addr, src_offset);
		dill_addpi(c, tmp_dest, dest_addr, dest_offset);
		src_addr = tmp_src;
		dest_addr = tmp_dest;
	    }
	    dill_scallv(c, (void*)transpose_array, "transpose_array",
			"%p%p%p%I%I%I%p", dimen_reg, src_addr, dest_addr,
			source_column_major, dest_type, conv->dest_size, spec_reg);
	}
	break;
    }
}

extern void
new_generate_conversion_code(dill_stream c, ConvStatus conv_status, IOConversionPtr conv, dill_reg *args,
			     int assume_align, int register_args)
{
    int i;
    dill_reg src_addr = args[0];
    dill_reg dest_addr = args[1];
    dill_reg rt_conv_status = args[2];
    int control_base = -1;
    dill_reg addr_reg = -1;

    for (i = 0; i < conv->conv_count; i++) {
	FMTypeDesc *next = &conv->conversions[i].iovar->type_desc;
	while (next != NULL) {
	    if ((next->type == FMType_array) && (next->static_size == 0)) {
		int src_is_aligned = 1;
		int src_drisc_type;
		int src_required_align;
		FMFormat f = conv->ioformat->body;
		int field = next->control_field_index;
		struct _FMgetFieldStruct tmp_src_spec;
		iogen_oprnd src_oprnd;

		memset(&tmp_src_spec, 0, sizeof(tmp_src_spec));
		tmp_src_spec.size = f->field_list[field].field_size;
		tmp_src_spec.offset = f->field_list[field].field_offset;
		tmp_src_spec.data_type = integer_type;
		tmp_src_spec.byte_swap = conv->ioformat->body->byte_reversal;
		src_drisc_type = drisc_type(&tmp_src_spec);
		src_required_align = TYPE_ALIGN(c, src_drisc_type);
		if ((assume_align < src_required_align) || 
		    (((tmp_src_spec.offset) % src_required_align) == 0)) {
		    src_is_aligned = 0;
		}
		src_oprnd = gen_fetch(c, src_addr, 
				      tmp_src_spec.offset,
				      tmp_src_spec.size,
				      tmp_src_spec.data_type,
				      src_is_aligned, tmp_src_spec.byte_swap);
		if (src_oprnd.size != sizeof(int)) {
		    iogen_oprnd tmp_oprnd;
		    tmp_oprnd = gen_size_conversion(c, src_oprnd, sizeof(int));
		    free_oprnd(c, src_oprnd);
		    src_oprnd = tmp_oprnd;
		}

		if (control_base == -1) {
		    control_base = ffs_localb(c, sizeof(int) * f->field_count);
#ifdef RAW
		}
		gen_store(c, src_oprnd, dill_lp(c), control_base + field * sizeof(int),
			  sizeof(int), integer_type, TRUE /* aligned */ );
#else
		    addr_reg = dill_getreg(c, DILL_P);
		    dill_virtual_lea(c, addr_reg, control_base);
		    conv_status->control_value = (size_t*)(intptr_t)addr_reg;
	        }
	    assert(addr_reg != -1);
	        gen_store(c, src_oprnd, addr_reg, field*sizeof(int),
			  sizeof(int), integer_type, TRUE /* aligned */);
#endif
		free_oprnd(c, src_oprnd);
	    }
	    next = next->next;
	}
    }
    for (i = 0; i < conv->conv_count; i++) {
	FMFieldPtr src_spec = &conv->conversions[i].src_field;
	FMTypeDesc *type_desc = &conv->conversions[i].iovar->type_desc;
	int byte_swap = conv->conversions[i].src_field.byte_swap;
	size_t dest_offset = conv->conversions[i].dest_offset;
	size_t src_offset = src_spec->offset;
	int elements = 
	    get_static_array_element_count(conv->conversions[i].iovar);
	if (conv->conversions[i].src_field.size == 1) byte_swap = 0;
	if (conv->conversions[i].default_value) {
	    iogen_oprnd src_oprnd;
	    int dst_is_aligned = (assume_align >= sizeof(long)) &
		(((conv->conversions[i].dest_offset) % 8) == 0);
	    src_oprnd = gen_set(c, conv->conversions[i].dest_size, 
				conv->conversions[i].default_value);
	    gen_store(c, src_oprnd, dest_addr,
		      conv->conversions[i].dest_offset,
		      conv->conversions[i].dest_size,
		      src_oprnd.data_type, dst_is_aligned);
	    free_oprnd(c, src_oprnd);
	} else if (!byte_swap &&
		   (src_spec->src_float_format == src_spec->target_float_format) &&
		   (src_spec->size == conv->conversions[i].dest_size) &&
		   (conv->conversions[i].subconversion == NULL) &&
		   (type_desc->type != FMType_pointer) &&
		   (type_desc->type != FMType_string) &&
		   (elements != -1) &&
		   (conv->conversions[i].rc_swap == no_row_column_swap) &&
		   ((src_spec->data_type != string_type))) {
	    /* data movement is all that is required */
	    int total_size = src_spec->size * elements;

	    if (total_size <= sizeof(long)) {
		iogen_oprnd src_oprnd;
		int src_is_aligned = (assume_align >= sizeof(long)) &
		    (((src_spec->offset) % 8) == 0);
		int dst_is_aligned = (assume_align >= sizeof(long)) &
		    (((conv->conversions[i].dest_offset) % 8) == 0);
		src_oprnd = gen_fetch(c, src_addr, 
				      src_spec->offset,
				      total_size, src_spec->data_type,
				      src_is_aligned, 0);
		gen_store(c, src_oprnd, dest_addr,
			  conv->conversions[i].dest_offset,
			  conv->conversions[i].dest_size,
			  src_spec->data_type, dst_is_aligned);
		free_oprnd(c, src_oprnd);
	    } else {
		gen_memcpy(c, src_addr, src_spec->offset, 
			   dest_addr,
			   conv->conversions[i].dest_offset,
			   0, total_size);
	    }
	    continue;
	}
	generate_convert_field(c, conv_status, src_addr, src_offset, dest_addr, dest_offset,
			       rt_conv_status, &conv->conversions[i], type_desc, 
			       conv->conversion_type != copy_dynamic_portion);
    }
}

#endif
