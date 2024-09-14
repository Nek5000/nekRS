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
#include <sys/types.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#include <limits.h>

#include "assert.h"
#include "fm.h"
#include "fm_internal.h"

struct _xml_output_info {
    char *prestring;  /* output this before this field */
    int field_num;
    FMFieldPtr field_ptr;
    FMdata_type base_type;
    FMVarInfoStruct *iovar;	/* for dynamic arrays, array size */
    char *element_prestring;
    char *element_poststring;
    char *poststring;  /* output this after this field */
};

typedef struct dstring {
    char *string;
    int length;
    int max_alloc;
} *ffs_dstring;

static int get_int_attr(int *int_ptr, char *name, char *start, char *end);
static int get_str_attr(char **str_ptr, char *name, char *start, char *end);
static void internal_record_to_XML_string(FMFormat format, void *data,
					  void *string_base, ffs_dstring ds,
					  int encoded);
static void
generic_field_to_XML(FMFormat format, int field,
			   void *data, void *string_base, int encode, 
			   ffs_dstring ds);


static int
is_tag(char *tag, char *xml_stream)
{
    int len = (int) strlen(tag);
    if (strncmp(xml_stream, tag, len) != 0) return 0;
    if (isalnum((int)xml_stream[len])) return 0;  /* not end of stream tag */
    if (xml_stream[len] == '_') return 0;
    if (xml_stream[len] == '-') return 0;
    return 1;
}

extern void
free_XML_output_info(xml_output_info info)
{
    int i = 0;
    while(info[i].field_ptr == NULL) {
	if (info[i].prestring != NULL) free(info[i].prestring);
	if (info[i].poststring != NULL) free(info[i].poststring);
	if (info[i].element_prestring != NULL) 
	    free(info[i].element_prestring);
	if (info[i].element_poststring != NULL) 
	    free(info[i].element_poststring);
    }
    free(info);
}



static void
expand_dstring(ffs_dstring ds, int min)
{
    int len = ds->max_alloc;
    int add = ds->max_alloc >> 3;   /* add an eighth */
    if (add < min) add = min;
    if (add < 128) add = 128;
    ds->string = realloc(ds->string, len + add);
    ds->max_alloc = len + add;
}

static void
dcatstr(ffs_dstring ds, char *str)
{
    int len = (int) strlen(str);
    if (ds->length + len + 1 > ds->max_alloc) {
	expand_dstring(ds, len + 1);
    }
    strncat(&ds->string[ds->length], str, ds->max_alloc);
    ds->length += len;
}

static ffs_dstring
new_dstring()
{
    ffs_dstring dstr = malloc(sizeof(struct dstring));
    dstr->string = malloc(64);
    dstr->string[0] = 0;
    dstr->length = 0;
    dstr->max_alloc = 64;
    return dstr;
}

static ffs_dstring
record_to_XML_string(FMFormat format, void *data, int encoded)
{
    ffs_dstring dstr = new_dstring();
    internal_record_to_XML_string(format, data, data, dstr, encoded);
    return dstr;
}

static ffs_dstring
raw_to_XML_string(FMContext c, FMFormat format, void *data, int encoded)
{
    int index;
    ffs_dstring ds;
    if (FMhas_XML_info(format)) {
	return record_to_XML_string(format, data, encoded);
    }
    ds = new_dstring();
    dcatstr(ds, "<");
    dcatstr(ds, format->format_name);
    dcatstr(ds, ">\n");
    for (index = 0; index < format->field_count; index++) {
	generic_field_to_XML(format, index, data, data, 
			     encoded, ds);
    }
    dcatstr(ds, "</");
    dcatstr(ds, format->format_name);
    dcatstr(ds, ">\n");
    return ds;
}

extern void
dump_raw_IOrecord_as_XML(FMContext c, FMFormat format, void *data)
{
    ffs_dstring dstr = raw_to_XML_string(c, format, data, 1);
    printf("%s", dstr->string);
}

extern char*
FMunencoded_to_XML_string(FMContext fmcontext, FMFormat format, void *data)
{
    ffs_dstring dstr = raw_to_XML_string((FMContext)fmcontext, format, data, 0);
    char *str = dstr->string;
    free(dstr);
    return str;
}

extern void
internal_dump_XML_record(FMFormat format, void *data, void *string_base, int encoded)
{
    struct dstring ds;
    ds.string = malloc(64);
    ds.string[0] = 0;
    ds.length = 0;
    ds.max_alloc = 64;
    internal_record_to_XML_string(format, data, string_base, &ds, encoded);
    printf("%s", ds.string);
    free(ds.string);
}

static long
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

static void
internal_record_to_XML_string(FMFormat format, void *data, void *string_base, ffs_dstring ds, int encoded) 

{
    int i = -1;
    xml_output_info info;
    if (!FMhas_XML_info(format)) {
	FMdump_XML(format, data, encoded);
	return;
    }
    info = format->xml_out;
    while(info[++i].field_ptr != NULL) {
	int j;
	int array_dimen = get_static_array_element_count(info[i].iovar);
	void *base = (char*)data;
	FMFieldPtr field_ptr = info[i].field_ptr;
	struct _FMgetFieldStruct dummy_ptr;
	if (info[i].prestring != NULL) {
	    dcatstr(ds, info[i].prestring);
/*	    printf("%s", info[i].prestring);*/
	}
	if (array_dimen == -1) {
	    array_dimen = FMget_array_element_count(format, info[i].iovar, base, 1);
	    dummy_ptr = *field_ptr;
	    dummy_ptr.size = format->pointer_size;
	    dummy_ptr.data_type = integer_type;
	    base = get_FMaddr(&dummy_ptr, base, string_base, encoded);
	    dummy_ptr = *field_ptr;
	    field_ptr = &dummy_ptr;
	    dummy_ptr.offset = 0;
	}
	for(j=0; j<array_dimen; j++) {
	    if (info[i].element_prestring != NULL) {
		dcatstr(ds, info[i].element_prestring);
/*		printf("%s", info[i].element_prestring);*/
	    }
	    switch(info[i].base_type) {
	    case integer_type:
	    {
		char tmp[64];
		ssize_t l;
		l = get_FMlong(field_ptr, base);
		snprintf(tmp, sizeof(tmp), "%zd", l);
		dcatstr(ds, tmp);
		break;
	    }
	    case unsigned_type:
	    case enumeration_type:
	    {
		char tmp[64];
		size_t l;
		l = get_FMulong(field_ptr, base);
		snprintf(&tmp[0], sizeof(tmp), "%zu", l);
		dcatstr(ds, tmp);
		break;
	    }
	    case boolean_type:
	    {
		size_t l;
		l = get_FMulong(field_ptr, base);
		if (l) 
		    dcatstr(ds, "true");
		else
		    dcatstr(ds, "false");
		break;
	    }
	    case float_type:
	    {
		char tmp[64];
		double d;
		d = get_FMdouble(field_ptr, base);
		snprintf(tmp, sizeof(tmp), "%g", d);
		dcatstr(ds, tmp);
		break;
	    }
	    case char_type:
	    {
		char tmp[2];
		tmp[0] = get_FMchar(field_ptr, base);
		tmp[1] = 0;
		dcatstr(ds, tmp);
		break;
	    }
	    case string_type:
	    {
		char *s;
		s = get_FMaddr(field_ptr, base, string_base, encoded);
		if (s) dcatstr(ds, s);
		break;
	    }
	    case unknown_type:
	    {
		int field = info[i].field_num;
		size_t offset = field_ptr->offset;
		FMFormat sub = format->field_subformats[field];
		if (sub != NULL) {
		    internal_record_to_XML_string(sub, (char*)base+offset,
						  string_base, ds, encoded);
		} else {
		    dcatstr(ds, "<INCONSISTENCY!  SUBFORMAT NOT FOUND FOR FIELD TYPE \"");
		    dcatstr(ds, (char*)format->field_list[field].field_type);
		    dcatstr(ds, "\">\n");
		}
		break;
	    }
	    }  /* end switch */
	    if (info[i].element_poststring != NULL) {
		dcatstr(ds, info[i].element_poststring);
	    }
	    if ((info[i].element_prestring == NULL) && 
		(info[i].element_poststring == NULL)
		&& (info[i].base_type != unknown_type)) {
		dcatstr(ds, " ");
	    }

	    base =  (char *) base + field_ptr->size;
	}
	if (info[i].poststring != NULL) {
	    dcatstr(ds, info[i].poststring);
	}
    }
}

extern void
FMdump_XML(FMFormat format, void *data, int encoded)
{
    internal_dump_XML_record(format, data, data, encoded);
}

static int
get_field_num(int* field_num_ptr, FMFormat format, char *left, char *right)
{
    char *field_name = NULL;
    int field_num;
    if (!get_str_attr(&field_name, "field_name", left, right)) {
	if (!get_int_attr(&field_num, "field_id", left, right)) {
	    printf("Neither Field_name nor Field_id attribute found\n");
	    return 0;
	}
    } else {
	field_num = 0;
	while(strcmp(format->field_list[field_num].field_name,
		     field_name) != 0) {
	    field_num++;
	    if (field_num >= format->field_count) {
		printf("Field name %s not found in format\n", field_name);
		free(field_name);
		return 0;
	    }
	}
	free(field_name);
    }
    if (field_num >= format->field_count) {
	printf("Field number %d too big\n", field_num);
	return 0;
    }
    *field_num_ptr = field_num;
    return 1;
}

extern void *
get_optinfo_FMFormat(FMFormat format, int info_type, size_t *len_p)
{
    int i = 0;
    if (format->opt_info == NULL) return NULL;
    while(format->opt_info[i].info_type != 0) {
	if (format->opt_info[i].info_type == info_type) {
	    *len_p = format->opt_info[i].info_len;
	    return (void*)format->opt_info[i].info_block;
	}
	i++;
    }
    return NULL;
}

extern FMfloat_format fm_my_float_format;

static xml_output_info
build_XML_output_info(FMFormat format)
{
    ssize_t i = 0, last_mark;
    char *XML_opt_data = NULL;
    size_t XML_data_len = 0;
    xml_output_info info = NULL;
    size_t info_count = 0;

    XML_opt_data = get_optinfo_FMFormat(format, 0x584D4C20, &XML_data_len);

    if (XML_opt_data == NULL) return NULL;
    last_mark = 0;
    i = -1;
    while (++i < (intptr_t) XML_data_len) {
	if (XML_opt_data[i] != '<') {
	    continue;
	}
	if (is_tag("/FFS:array", &XML_opt_data[i+1])) {
	    /* array terminatfmn */
	    ssize_t post_len = i - last_mark + 1;
	    info[info_count].element_poststring = malloc(post_len);
	    strncpy(info[info_count].element_poststring, 
		    &XML_opt_data[last_mark], post_len - 1);
	    info[info_count].element_poststring[post_len -1] = 0;
	    info_count++;
	    i += 13;
	    last_mark = i;
	    continue;
	}
	if (is_tag("FFS:data", &XML_opt_data[i+1])) {
	    size_t pre_len = i - last_mark;
	    FMdata_type data_type;
	    int field_num;
	    long junk;
	    char *right_end = strchr(&XML_opt_data[i+9], '>');
	    char *prestring = malloc(pre_len + 1);
	    FMFieldPtr field_ptr = NULL;
	    
	    strncpy(prestring, &XML_opt_data[last_mark], pre_len);
	    prestring[pre_len] = 0;
	    if (!get_field_num(&field_num, format, XML_opt_data + i + 9,
			       right_end)) {
		i = right_end - XML_opt_data;
		continue;
	    }
	    data_type = FMarray_str_to_data_type(format->field_list[field_num].field_type,
					       &junk);
	    field_ptr = malloc(sizeof(struct _FMgetFieldStruct));
	    field_ptr->offset = format->field_list[field_num].field_offset;
	    field_ptr->size = format->field_list[field_num].field_size;
	    field_ptr->data_type = data_type;
	    field_ptr->byte_swap = format->byte_reversal;
	    field_ptr->src_float_format = format->float_format;
	    field_ptr->target_float_format = fm_my_float_format;
	    if (info_count == 0) {
		info = malloc(sizeof(*info) * 2);
	    } else {
		info = realloc(info, sizeof(*info) * (info_count + 2));
	    }
	    info[info_count+1].field_ptr = NULL;
	    info[info_count].prestring = prestring;
	    info[info_count].field_ptr = field_ptr;
	    info[info_count].base_type = data_type;
	    info[info_count].poststring = NULL;
	    info[info_count].field_num = field_num;
	    info[info_count].iovar = &format->var_list[field_num];
	    info[info_count].element_prestring = NULL;
	    info[info_count].element_poststring = NULL;
	    info_count++;
	    i = right_end - XML_opt_data;
	    last_mark = i + 1;
	    continue;
	}
	if (is_tag("FFS:array", &XML_opt_data[i+1])) {
	    /* array terminatfmn */
	    size_t pre_len = i - last_mark;
	    char *right_end = strchr(&XML_opt_data[i+10], '>');
	    char *prestring = malloc(pre_len + 1);
	    
	    strncpy(prestring, &XML_opt_data[last_mark], pre_len);
	    prestring[pre_len] = 0;
	    if (info_count == 0) {
		info = malloc(sizeof(*info) * 2);
	    } else {
		info = realloc(info, sizeof(*info) * (info_count + 2));
	    }
	    info[info_count+1].field_ptr = NULL;
	    info[info_count].prestring = prestring;
	    info[info_count].poststring = NULL;
	    i = right_end - XML_opt_data;
	    last_mark = i + 1;
	    continue;
	}
	if (is_tag("FFS:array_data_mark", &XML_opt_data[i+1])) {
	    /* array terminatfmn */
	    size_t pre_len = i - last_mark;
	    int field_num;
	    char *right_end = strchr(&XML_opt_data[i+9], '>');
	    char *prestring = malloc(pre_len + 1);
	    long elements;
	    FMdata_type data_type;
	    FMFieldPtr field_ptr = NULL;
	    
	    strncpy(prestring, &XML_opt_data[last_mark], pre_len);
	    prestring[pre_len] = 0;
	    if (!get_field_num(&field_num, format, XML_opt_data + i + 9, 
			       right_end)) {
		i = right_end - XML_opt_data;
		continue;
	    }
	    info[info_count].element_prestring = prestring;
	    data_type = FMarray_str_to_data_type(format->field_list[field_num].field_type,
					       &elements);
	    info[info_count].iovar = &format->var_list[field_num];
	    info[info_count].base_type = data_type;
	    info[info_count].field_num = field_num;
	    field_ptr = malloc(sizeof(struct _FMgetFieldStruct));
	    field_ptr->offset = format->field_list[field_num].field_offset;
	    field_ptr->size = format->field_list[field_num].field_size;
	    field_ptr->data_type = data_type;
	    field_ptr->byte_swap = format->byte_reversal;
	    field_ptr->src_float_format = format->float_format;
	    field_ptr->target_float_format = fm_my_float_format;
	    info[info_count].field_ptr = field_ptr;
	    i = right_end - XML_opt_data;
	    last_mark = i + 1;
	    continue;
	}
    }
    if (info_count != 0) {
	size_t post_len = XML_data_len - last_mark;
	char *poststring = malloc(post_len + 1);
	strncpy(poststring, &XML_opt_data[last_mark], post_len);
	poststring[post_len] = 0;
	info[info_count-1].poststring = poststring;
	info[info_count].field_ptr = NULL;
/*	    for(i=0; i < info_count; i++) {
	    printf("XML field info:\n");
	    printf("    prestring = \"%s\"\n", info[i].prestring);
	    printf("    field num = \"%d\"\n", info[i].field_num);
	    printf("    field name = \"%s\"\n",
	       format->field_list[info[i].field_num].field_name);
	    printf("    poststring = \"%s\"\n", info[i].poststring == NULL ?
	        "(NULL)":info[i].poststring);
	}*/
	return info;
    } else {
	return NULL;
    }
}

static int
get_int_attr(int *int_ptr, char *name, char *start, char *end)
{
    int name_len = (int)strlen(name);
    start -= 1;
    while(++start < (end - name_len)) {
	if (*start == *name) {
	    /* first char eq, check further */
	    if (is_tag(name, start)) {
		char *point = &start[name_len];
		char *tmp;
		while (isspace((int)*point)) point++;
		if (point > end) continue;
		if (*point != '=') continue;
		point++;
		*int_ptr = strtol(point, &tmp, 10);
		if (tmp != point) {
		    return 1;
		} else {
		    return 0;
		}
	    }
	}
    }
    return 0;
}

static int
get_str_attr(char **str_ptr, char *name, char *start, char *end)
{
    int name_len = (int)strlen(name);
    start -= 1;
    while(++start < (end - name_len)) {
	if (*start == *name) {
	    /* first char eq, check further */
	    if (is_tag(name, start)) {
		char *point = &start[name_len];
		char *attr = NULL;
		while (isspace((int)*point)) point++;
		if (point > end) continue;
		if (*point != '=') continue;
		point++;
		while (isspace((int)*point)) point++;
		if (point > end) continue;
		if (*point == '"') {
		    /* begin quote */
		    char *end_quote = strchr(point+1, '"');
		    point++;
		    if ((end_quote == 0) || (end_quote > end)) continue;
		    attr = malloc(end_quote - point + 1);
		    strncpy(attr, point, end_quote - point);
		    attr[end_quote - point] = 0;
		    *str_ptr = attr;
		    return 1;
		} else {
		    char *next_space = point;
		    while(isalnum((int)*next_space) || (*next_space=='_')) 
			next_space++;
		    attr = malloc(next_space - point+1);
		    strncpy(attr, point, next_space - point);
		    attr[next_space - point] = 0;
		    *str_ptr = attr;
		    return 1;
		}
	    }
	}
    }
    return 0;
}

extern int
FMhas_XML_info(FMFormat format)
{
    if (format->xml_out == (void*) -1) return 0;
    if (format->xml_out == NULL) {
	/* search opt info */
	format->xml_out = build_XML_output_info(format);
	if (format->xml_out == NULL) {
	    format->xml_out = (void*)-1;
	    return 0;
	}
    }
    return 1;
}

extern void
dump_FMfield_as_XML(FMContext c, FMFormat format, int field, void *data,
		    void *string_base, int encode, int verbose)
{
    ffs_dstring ds = new_dstring();
    generic_field_to_XML(format, field, data, string_base, encode,
			 ds);
    printf("%s", ds->string);
    free(ds->string);
    free(ds);
}

extern int
add_single_value_as_XML(const char *field_type, int field_size, int field_offset, void *data,
			void *string_base, int byte_reversal, int float_format, int encode,
			ffs_dstring ds)
{
    FMgetFieldStruct descr;  /* OK */
    char str[64];
    descr.offset = field_offset;
    descr.size = field_size;
    descr.data_type = FMstr_to_data_type(field_type);
    descr.byte_swap = byte_reversal;
    descr.src_float_format = float_format;
    descr.target_float_format = fm_my_float_format;

    str[0] = 0;
    if (descr.data_type == integer_type) {
	if (field_size <= sizeof(long)) {
	    size_t tmp = get_FMlong(&descr, data);
	    snprintf(str, sizeof(str), "%zd", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long;
	    long high_long;
	    get_FMlong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		snprintf(str, sizeof(str), "%ld ", low_long);
	    } else {
		snprintf(str, sizeof(str), "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    snprintf(str, sizeof(str), "<scalar type=\"long\" />");
	} else {
	    snprintf(str, sizeof(str), "<scalar type=\"int\" size=\"%d\" />", field_size);
	}
    } else if (descr.data_type == unsigned_type) {
	if (field_size <= sizeof(unsigned long)) {
	    size_t tmp = get_FMulong(&descr, data);
	    snprintf(str, sizeof(str), "%zu ", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long, high_long;
	    get_FMulong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		snprintf(str, sizeof(str), "%lu ", low_long);
	    } else {
		snprintf(str, sizeof(str), "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    snprintf(str, sizeof(str), "<scalar type=\"unsignedLong\" />");
	} else {
	    snprintf(str, sizeof(str), "<scalar type=\"unsignedInt\" size=\"%d\" />", field_size);
	}
    } else if (descr.data_type == enumeration_type) {
	snprintf(str, sizeof(str), "%u ", *(int *) ((char *) data + field_offset));
    } else if (descr.data_type == boolean_type) {
	if (*(int *) ((char *) data + field_offset) == 0) {
	    strncpy(str, "false ", sizeof(str));
	} else {
	    strncpy(str, "true ", sizeof(str));
	}
    } else if (descr.data_type == float_type) {
	if (field_size == sizeof(float)) {
	    float tmp = get_FMfloat(&descr, data);
	    snprintf(str, sizeof(str), "%g ", tmp);
	} else if (field_size == sizeof(double)) {
	    double tmp = get_FMdouble(&descr, data);
	    snprintf(str, sizeof(str), "%g ", tmp);
#if SIZEOF_LONG_DOUBLE != 0 && SIZEOF_LONG_DOUBLE != SIZEOF_DOUBLE
	} else if (field_size == sizeof(long double)) {
	    long double tmp;
	    memcpy(&tmp, (float *) ((char *) data + field_offset),
		   sizeof(double));
	    snprintf(str, sizeof(str), "%Lg ", tmp);
#endif
	} else {
	    if (field_size < sizeof(float)) {
		snprintf(str, sizeof(str), "<scalar type=\"small-float\" />");
	    } else if (field_size > sizeof(double)) {
		snprintf(str, sizeof(str), "<scalar type=\"big-float\" />");
	    } else {
		snprintf(str, sizeof(str), "<scalar type=\"float\" size=\"%u\" />", field_size);
	    }
	}
    } else if (descr.data_type == char_type) {
	snprintf(str, sizeof(str), "%c ", *(char *) ((char *) data + field_offset));
    } else if (descr.data_type == string_type) {
	char *tmp_str = (char *) get_FMaddr(&descr, data, string_base, encode);
	if (tmp_str != 0) {
	    dcatstr(ds, tmp_str);
	    return 1;
	}
    } else {
	return 0;
    }
    dcatstr(ds, str);
    return 1;
}

static void
add_value_as_XML(const char *field_type, int field_size, int field_offset, FMFormat format, int field,
		 void *data, void *string_base, int byte_reversal, int float_format, int encode,
		 int in_array, ffs_dstring ds)
{
    if (add_single_value_as_XML(field_type, field_size, field_offset, 
				data, string_base, byte_reversal, 
				float_format, encode, ds)) {
	/* nothing */
    } else {
	char *typ = base_data_type(field_type);
	if ((format = format->field_subformats[field]) != NULL) {
	    int index;
	    if (in_array) {
		dcatstr(ds, "<");
		dcatstr(ds, format->format_name);
		dcatstr(ds, ">\n");
	    } else {
		dcatstr(ds, "\n");
	    }
	    for (index = 0; index < format->field_count; index++) {
		generic_field_to_XML(format, index,
				     (void *) ((char *) data + field_offset),
				     string_base, encode, ds);
	    }
	    if (in_array) {
		dcatstr(ds, "</");
		dcatstr(ds, format->format_name);
		dcatstr(ds, ">\n");
	    }
	} else {
	    printf("<scalar type=\"unknown\" />\n");
	}
	free(typ);
    }
}

static void
generic_field_to_XML(FMFormat format, int field, void *data, void *string_base,
		     int encode, ffs_dstring ds)
{
    FMFieldList fmfield = &format->field_list[field];
    FMVarInfoList iovar = &format->var_list[field];
    int field_offset = fmfield->field_offset;
    int field_size = fmfield->field_size;
    const char *field_type = fmfield->field_type;
    int byte_reversal = format->byte_reversal;
    int float_format = format->float_format;
    char *left_paren = NULL;

    if (encode == 0) byte_reversal = 0;

    dcatstr(ds, "<");
    dcatstr(ds, (char*)fmfield->field_name);
    dcatstr(ds, ">");
    if ((left_paren = strchr(field_type, '[')) == NULL) {
	add_value_as_XML(field_type, field_size, field_offset, format, field, data,
			 string_base, byte_reversal, float_format, 
			 encode, 0, ds);
    } else if (strchr(left_paren + 1, '[') == NULL) {
	/* single dimen array */
	long dimension = 0;
	char sub_type[64];
	int sub_field_size;
	int offset = fmfield->field_offset;
	dcatstr(ds, "\n");
	*left_paren = 0;
	strncpy(sub_type, field_type, sizeof(sub_type));
	*left_paren = '[';
	dimension = strtol(left_paren + 1, NULL, 10);

	if ((dimension == LONG_MIN) || (dimension == LONG_MAX) ||
	    (dimension == 0)) {
	    if (!iovar->var_array) {
		fprintf(stderr, "Couldn't parse array size in \"%s\"\n",
			field_type);
		return;
	    } else {
		FMgetFieldStruct descr;  /* OK */
		size_t tmp_offset;

		dimension = FMget_array_element_count(format, iovar, data, encode);

		descr.offset = fmfield->field_offset;
		descr.size = format->pointer_size;
		descr.data_type = integer_type;
		descr.byte_swap = byte_reversal;
		tmp_offset = get_FMlong(&descr, data);
		if (encode) {
		    data = (char *) string_base + tmp_offset;
		} else {
		    data = (void *) tmp_offset;
		}
		offset = 0;
	    }
	    sub_field_size = fmfield->field_size;
	} else {
	    /* normal internal array */
	    sub_field_size = fmfield->field_size;
	}
	for (; dimension > 0; dimension--) {
	    add_value_as_XML(sub_type, sub_field_size, offset, format, field, data,
			     string_base, byte_reversal, float_format,
			     encode, 1, ds);
	    offset += sub_field_size;
	}
    } else {
	/* double dimen array */
	long dimension1 = 0;
	long dimension2 = 0;
	char *temp_ptr;
	char sub_type[64];
	int sub_field_size, offset = fmfield->field_offset;
	dcatstr(ds, "\n");
	*left_paren = 0;
	strncpy(sub_type, field_type, sizeof(sub_type));
	dimension1 = strtol(left_paren + 1, &temp_ptr, 10);
	dimension2 = strtol(temp_ptr + 2, &temp_ptr, 10);

	if ((dimension2 == LONG_MIN) || (dimension2 == LONG_MAX) ||
	    (dimension2 == 0) || (dimension1 == 0)) {
	    *left_paren = '[';
	    fprintf(stderr, "Couldn't parse array size in \"%s\"\n",
		    field_type);
	    return;
	}
	*left_paren = '[';
	sub_field_size = fmfield->field_size;
	for (; dimension2 > 0; dimension2--) {
	    int i = 0;
	    for (; i < dimension1; i++) {
		add_value_as_XML(sub_type, sub_field_size, offset, format, field, 
				 data, string_base, byte_reversal, 
				 float_format, encode, 1, ds); 
		offset += sub_field_size;
	    }
	}
    }

    dcatstr(ds, "</");
    dcatstr(ds, (char*)fmfield->field_name);
    dcatstr(ds, ">\n");
}

