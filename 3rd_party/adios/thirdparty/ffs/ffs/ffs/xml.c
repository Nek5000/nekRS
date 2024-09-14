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
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#include "unix_defs.h"
#include <limits.h>

#include "assert.h"
#include "ffs.h"
#include "ffs_internal.h"

struct _xml_output_info {
    char *prestring;  /* output this before this field */
    int field_num;
    IOFieldPtr field_ptr;
    IOdata_type base_type;
    IOVarInfoStruct *iovar;	/* for dynamic arrays, array size */
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
static void internal_record_to_XML_string(IOFormat ioformat, void *data,
					  void *string_base, ffs_dstring ds,
					  int encoded);
static void
generic_field_to_XML ARGS((IOFile iofile, IOFormat ioformat, int field,
			   void *data, void *string_base, int encode, 
			   ffs_dstring ds));


static int
is_tag(char *tag, char *xml_stream)
{
    int len = strlen(tag);
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
    int len = strlen(str);
    if (ds->length + len + 1 > ds->max_alloc) {
	expand_dstring(ds, len + 1);
    }
    strcat(&ds->string[ds->length], str);
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
record_to_XML_string(IOFormat ioformat, void *data, int encoded)
{
    ffs_dstring dstr = new_dstring();
    internal_record_to_XML_string(ioformat, data, data, dstr, encoded);
    return dstr;
}

static ffs_dstring
raw_to_XML_string(iofile, ioformat, data, encoded)
IOFile iofile;
IOFormat ioformat;
void *data;
int encoded;
{
    int index;
    ffs_dstring ds;
    if (IOhas_XML_info(ioformat)) {
	return record_to_XML_string(ioformat, data, encoded);
    }
    ds = new_dstring();
    dcatstr(ds, "<");
    dcatstr(ds, ioformat->body->format_name);
    dcatstr(ds, ">\n");
    for (index = 0; index < ioformat->body->field_count; index++) {
	generic_field_to_XML(iofile, ioformat, index, data, data, 
			     encoded, ds);
    }
    dcatstr(ds, "</");
    dcatstr(ds, ioformat->body->format_name);
    dcatstr(ds, ">\n");
    return ds;
}

extern void
dump_raw_IOrecord_as_XML(iofile, ioformat, data)
IOFile iofile;
IOFormat ioformat;
void *data;
{
    ffs_dstring dstr = raw_to_XML_string(iofile, ioformat, data, 1);
    printf("%s", dstr->string);
}

extern char*
IOencoded_to_XML_string(IOContext iocontext, void *data)
{
    IOFormat ioformat = get_format_IOcontext(iocontext, data);
    int header_size;
    ffs_dstring dstr = NULL;
    char *str = NULL;

    if (ioformat == NULL) return NULL;
    header_size = header_size_IOformat(ioformat);
    dstr = raw_to_XML_string((IOFile)iocontext, ioformat, 
			     (char*)data + header_size, 1);
    str = dstr->string;
    free(dstr);
    return str;
}

extern char*
IOunencoded_to_XML_string(IOContext iocontext, IOFormat ioformat, void *data)
{
    ffs_dstring dstr = raw_to_XML_string((IOFile)iocontext, ioformat, data, 0);
    char *str = dstr->string;
    free(dstr);
    return str;
}

extern void
dump_encoded_as_XML(IOContext context, void *data)
{
    IOFormat ioformat = get_format_IOcontext(context, data);
    int header_size;
    if (ioformat == NULL) {
	printf("Encoded record has no valid IOFormat\n");
	return;
    }
    header_size = header_size_IOformat(ioformat);
    dump_XML_record(ioformat, (char*)data + header_size, 1);
}

extern void
internal_dump_XML_record(IOFormat ioformat, void *data, void *string_base, int encoded)
{
    struct dstring ds;
    ds.string = malloc(64);
    ds.string[0] = 0;
    ds.length = 0;
    ds.max_alloc = 64;
    internal_record_to_XML_string(ioformat, data, string_base, &ds, encoded);
    printf("%s", ds.string);
    free(ds.string);
}

static void
internal_record_to_XML_string(IOFormat ioformat, void *data, void *string_base, ffs_dstring ds, int encoded) 

{
    int i = -1;
    xml_output_info info;
    if (!IOhas_XML_info(ioformat)) {
	dump_raw_IOrecord_as_XML(ioformat->context, ioformat, data);
	return;
    }
    info = ioformat->body->xml_out;
    while(info[++i].field_ptr != NULL) {
	int j;
	int array_dimen = get_static_array_element_count(info[i].iovar);
	void *base = (char*)data;
	IOFieldPtr field_ptr = info[i].field_ptr;
	struct _IOgetFieldStruct dummy_ptr;
	if (info[i].prestring != NULL) {
	    dcatstr(ds, info[i].prestring);
/*	    printf("%s", info[i].prestring);*/
	}
	if (array_dimen == -1) {
	    array_dimen = get_array_element_count(info[i].iovar, base, 1);
	    dummy_ptr = *field_ptr;
	    dummy_ptr.size = ioformat->body->pointer_size;
	    dummy_ptr.data_type = integer_type;
	    base = get_IOaddr(&dummy_ptr, base, string_base, encoded);
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
		long l;
		l = get_IOlong(field_ptr, base);
		sprintf(tmp, "%ld", l);
		dcatstr(ds, tmp);
		break;
	    }
	    case unsigned_type:
	    case enumeration_type:
	    {
		char tmp[64];
		unsigned long l;
		l = get_IOulong(field_ptr, base);
		sprintf(tmp, "%lu", l);
		dcatstr(ds, tmp);
		break;
	    }
	    case boolean_type:
	    {
		unsigned long l;
		l = get_IOulong(field_ptr, base);
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
		d = get_IOdouble(field_ptr, base);
		sprintf(tmp, "%g", d);
		dcatstr(ds, tmp);
		break;
	    }
	    case char_type:
	    {
		char tmp[2];
		tmp[0] = get_IOchar(field_ptr, base);
		tmp[1] = 0;
		dcatstr(ds, tmp);
		break;
	    }
	    case string_type:
	    {
		char *s;
		s = get_IOaddr(field_ptr, base, string_base, encoded);
		if (s) dcatstr(ds, s);
		break;
	    }
	    case unknown_type:
	    {
		int field = info[i].field_num;
		int offset = field_ptr->offset;
		IOFormat sub = ioformat->field_subformats[field];
		if (sub != NULL) {
		    internal_record_to_XML_string(sub, (char*)base+offset,
						  string_base, ds, encoded);
		} else {
		    dcatstr(ds, "<INCONSISTENCY!  SUBFORMAT NOT FOUND FOR FIELD TYPE \"");
		    dcatstr(ds, (char*)ioformat->body->field_list[field].field_type);
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
dump_XML_record(IOFormat ioformat, void *data, int encoded)
{
    internal_dump_XML_record(ioformat, data, data, encoded);
}

static int
get_field_num(int* field_num_ptr, IOFormat ioformat, char *left, char *right)
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
	while(strcmp(ioformat->body->field_list[field_num].field_name,
		     field_name) != 0) {
	    field_num++;
	    if (field_num >= ioformat->body->field_count) {
		printf("Field name %s not found in format\n", field_name);
		free(field_name);
		return 0;
	    }
	}
	free(field_name);
    }
    if (field_num >= ioformat->body->field_count) {
	printf("Field number %d too big\n", field_num);
	return 0;
    }
    *field_num_ptr = field_num;
    return 1;
}

extern void *
get_optinfo_IOFormat(IOFormat ioformat, int info_type, int *len_p)
{
    int i = 0;
    if (ioformat->body->opt_info == NULL) return NULL;
    while(ioformat->body->opt_info[i].info_type != 0) {
	if (ioformat->body->opt_info[i].info_type == info_type) {
	    *len_p = ioformat->body->opt_info[i].info_len;
	    return ioformat->body->opt_info[i].info_block;
	}
	i++;
    }
    return NULL;
}

static xml_output_info
build_XML_output_info(IOFormat ioformat)
{
    int i = 0, last_mark;
    char *XML_opt_data = NULL;
    int XML_data_len = 0;
    xml_output_info info = NULL;
    int info_count = 0;

    XML_opt_data = get_optinfo_IOFormat(ioformat, 0x584D4C20, &XML_data_len);

    if (XML_opt_data == NULL) return NULL;
    last_mark = 0;
    i = -1;
    while (++i < XML_data_len) {
	if (XML_opt_data[i] != '<') {
	    continue;
	}
	if (is_tag("/FFS:array", &XML_opt_data[i+1])) {
	    /* array termination */
	    int post_len = i - last_mark + 1;
	    info[info_count].element_poststring = malloc(post_len);
	    strncpy(info[info_count].element_poststring,
		    &XML_opt_data[last_mark], post_len -1);
	    info[info_count].element_poststring[post_len -1] = 0;
	    info_count++;
	    i += 13;
	    last_mark = i;
	    continue;
	}
	if (is_tag("FFS:data", &XML_opt_data[i+1])) {
	    int pre_len = i - last_mark;
	    IOdata_type data_type;
	    int field_num;
	    long junk;
	    char *right_end = strchr(&XML_opt_data[i+9], '>');
	    char *prestring = malloc(pre_len + 1);
	    IOFieldPtr field_ptr = NULL;
	    
	    strncpy(prestring, &XML_opt_data[last_mark], pre_len);
	    prestring[pre_len] = 0;
	    if (!get_field_num(&field_num, ioformat, XML_opt_data + i + 9,
			       right_end)) {
		i = right_end - XML_opt_data;
		continue;
	    }
	    data_type = array_str_to_data_type(ioformat->body->field_list[field_num].field_type,
					       &junk);
	    field_ptr = io_malloc(sizeof(struct _IOgetFieldStruct));
	    field_ptr->offset = ioformat->body->field_list[field_num].field_offset;
	    field_ptr->size = ioformat->body->field_list[field_num].field_size;
	    field_ptr->data_type = data_type;
	    field_ptr->byte_swap = ioformat->body->byte_reversal;
	    field_ptr->src_float_format = ioformat->body->float_format;
	    field_ptr->target_float_format = ffs_my_float_format;
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
	    info[info_count].iovar = &ioformat->body->var_list[field_num];
	    info[info_count].element_prestring = NULL;
	    info[info_count].element_poststring = NULL;
	    info_count++;
	    i = right_end - XML_opt_data;
	    last_mark = i + 1;
	    continue;
	}
	if (is_tag("FFS:array", &XML_opt_data[i+1])) {
	    /* array termination */
	    int pre_len = i - last_mark;
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
	    /* array termination */
	    int pre_len = i - last_mark;
	    int field_num;
	    char *right_end = strchr(&XML_opt_data[i+9], '>');
	    char *prestring = malloc(pre_len + 1);
	    long elements;
	    IOdata_type data_type;
	    IOFieldPtr field_ptr = NULL;
	    
	    strncpy(prestring, &XML_opt_data[last_mark], pre_len);
	    prestring[pre_len] = 0;
	    if (!get_field_num(&field_num, ioformat, XML_opt_data + i + 9, 
			       right_end)) {
		i = right_end - XML_opt_data;
		continue;
	    }
	    info[info_count].element_prestring = prestring;
	    data_type = array_str_to_data_type(ioformat->body->field_list[field_num].field_type,
					       &elements);
	    info[info_count].iovar = &ioformat->body->var_list[field_num];
	    info[info_count].base_type = data_type;
	    info[info_count].field_num = field_num;
	    field_ptr = io_malloc(sizeof(struct _IOgetFieldStruct));
	    field_ptr->offset = ioformat->body->field_list[field_num].field_offset;
	    field_ptr->size = ioformat->body->field_list[field_num].field_size;
	    field_ptr->data_type = data_type;
	    field_ptr->byte_swap = ioformat->body->byte_reversal;
	    field_ptr->src_float_format = ioformat->body->float_format;
	    field_ptr->target_float_format = ffs_my_float_format;
	    info[info_count].field_ptr = field_ptr;
	    i = right_end - XML_opt_data;
	    last_mark = i + 1;
	    continue;
	}
    }
    if (info_count != 0) {
	int post_len = XML_data_len - last_mark;
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
	       ioformat->body->field_list[info[i].field_num].field_name);
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
    int name_len = strlen(name);
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
    int name_len = strlen(name);
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
IOhas_XML_info(IOFormat ioformat)
{
    if (ioformat->body->xml_out == (void*) -1) return 0;
    if (ioformat->body->xml_out == NULL) {
	/* search opt info */
	ioformat->body->xml_out = build_XML_output_info(ioformat);
	if (ioformat->body->xml_out == NULL) {
	    ioformat->body->xml_out = (void*)-1;
	    return 0;
	}
    }
    return 1;
}

extern void
dump_IOfield_as_XML(iofile, ioformat, field, data, string_base, encode, verbose)
IOFile iofile;
IOFormat ioformat;
int field;
void *data;
void *string_base;
int encode;
int verbose;
{
    ffs_dstring ds = new_dstring();
    generic_field_to_XML(iofile, ioformat, field, data, string_base, encode,
			 ds);
    printf("%s", ds->string);
    free(ds->string);
    free(ds);
}

extern int
add_single_value_as_XML(field_type, field_size, field_offset, iofile, 
			data, string_base, byte_reversal, float_format,
			encode, ds)
char *field_type;
int field_size;
int field_offset;
IOFile iofile;
void *data;
void *string_base;
int byte_reversal;
int float_format;
int encode;
ffs_dstring ds;
{
    IOgetFieldStruct descr;  /* OK */
    char str[64];
    descr.offset = field_offset;
    descr.size = field_size;
    descr.data_type = str_to_data_type(field_type);
    descr.byte_swap = byte_reversal;
    descr.src_float_format = float_format;
    descr.target_float_format = ffs_my_float_format;

    str[0] = 0;
    if (descr.data_type == integer_type) {
	if (field_size <= sizeof(long)) {
	    long tmp = get_IOlong(&descr, data);
	    sprintf(str, "%ld", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long;
	    long high_long;
	    get_IOlong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		sprintf(str, "%ld ", low_long);
	    } else {
		sprintf(str, "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    sprintf(str, "<scalar type=\"long\" />");
	} else {
	    sprintf(str, "<scalar type=\"int\" size=\"%d\" />", field_size);
	}
    } else if (descr.data_type == unsigned_type) {
	if (field_size <= sizeof(unsigned long)) {
	    unsigned long tmp = get_IOulong(&descr, data);
	    sprintf(str, "%lu ", tmp);
	} else if (field_size == 2 * sizeof(long) && field_size == 8) {
	    unsigned long low_long, high_long;
	    get_IOulong8(&descr, data, &low_long, &high_long);
	    if (high_long == 0) {
		sprintf(str, "%lu ", low_long);
	    } else {
		sprintf(str, "0x%lx%08lx ", high_long, low_long);
	    }
	} else if (field_size > sizeof(long)) {
	    sprintf(str, "<scalar type=\"unsignedLong\" />");
	} else {
	    sprintf(str, "<scalar type=\"unsignedInt\" size=\"%d\" />", field_size);
	}
    } else if (descr.data_type == enumeration_type) {
	sprintf(str, "%u ", *(int *) ((char *) data + field_offset));
    } else if (descr.data_type == boolean_type) {
	if (*(int *) ((char *) data + field_offset) == 0) {
	    strcpy(str, "false ");
	} else {
	    strcpy(str, "true ");
	}
    } else if (descr.data_type == float_type) {
	if (field_size == sizeof(float)) {
	    float tmp = get_IOfloat(&descr, data);
	    sprintf(str, "%g ", tmp);
	} else if (field_size == sizeof(double)) {
	    double tmp = get_IOdouble(&descr, data);
	    sprintf(str, "%g ", tmp);
#if SIZEOF_LONG_DOUBLE != 0 && SIZEOF_LONG_DOUBLE != SIZEOF_DOUBLE
	} else if (field_size == sizeof(long double)) {
	    long double tmp;
	    memcpy(&tmp, (float *) ((char *) data + field_offset),
		   sizeof(double));
	    sprintf(str, "%Lg ", tmp);
#endif
	} else {
	    if (field_size < sizeof(float)) {
		sprintf(str, "<scalar type=\"small-float\" />");
	    } else if (field_size > sizeof(double)) {
		sprintf(str, "<scalar type=\"big-float\" />");
	    } else {
		sprintf(str, "<scalar type=\"float\" size=\"%u\" />", field_size);
	    }
	}
    } else if (descr.data_type == char_type) {
	sprintf(str, "%c ", *(char *) ((char *) data + field_offset));
    } else if (descr.data_type == string_type) {
	char *tmp_str = (char *) get_IOaddr(&descr, data, string_base, encode);
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
add_value_as_XML(field_type, field_size, field_offset, iofile, data, 
		  string_base, byte_reversal, float_format, encode, 
		 in_array, ds)
char *field_type;
int field_size;
int field_offset;
IOFile iofile;
void *data;
void *string_base;
int byte_reversal;
int float_format;
int encode;
int in_array;
ffs_dstring ds;
{
    IOFormat ioformat;
    if (add_single_value_as_XML(field_type, field_size, field_offset, 
				iofile, data, string_base, byte_reversal, 
				float_format, encode, ds)) {
	/* nothing */
    } else {
	char *typ = base_data_type(field_type);
	if ((ioformat = get_IOformat_by_name(iofile, typ)) != NULL) {
	    int index;
	    if (in_array) {
		dcatstr(ds, "<");
		dcatstr(ds, ioformat->body->format_name);
		dcatstr(ds, ">\n");
	    } else {
		dcatstr(ds, "\n");
	    }
	    for (index = 0; index < ioformat->body->field_count; index++) {
		generic_field_to_XML(iofile, ioformat, index,
				     (void *) ((char *) data + field_offset),
				     string_base, encode, ds);
	    }
	    if (in_array) {
		dcatstr(ds, "</");
		dcatstr(ds, ioformat->body->format_name);
		dcatstr(ds, ">\n");
	    }
	} else {
	    printf("<scalar type=\"unknown\" />\n");
	}
	free(typ);
    }
}

static void
generic_field_to_XML(iofile, ioformat, field, data, string_base, encode, ds)
IOFile iofile;
IOFormat ioformat;
int field;
void *data;
void *string_base;
int encode;
ffs_dstring ds;
{
    IOFieldList iofield = &ioformat->body->field_list[field];
    IOVarInfoList iovar = &ioformat->body->var_list[field];
    int field_offset = iofield->field_offset;
    int field_size = iofield->field_size;
    const char *field_type = iofield->field_type;
    int byte_reversal = ioformat->body->byte_reversal;
    int float_format = ioformat->body->float_format;
    char *left_paren = NULL;

    if (encode == 0) byte_reversal = 0;

    dcatstr(ds, "<");
    dcatstr(ds, (char*)iofield->field_name);
    dcatstr(ds, ">");
    if ((left_paren = strchr(field_type, '[')) == NULL) {
	add_value_as_XML(field_type, field_size, field_offset, iofile, data,
			 string_base, byte_reversal, float_format, 
			 encode, 0, ds);
    } else if (strchr(left_paren + 1, '[') == NULL) {
	/* single dimen array */
	long dimension = 0;
	char sub_type[64];
	int sub_field_size;
	int offset = iofield->field_offset;
	dcatstr(ds, "\n");
	*left_paren = 0;
	strcpy(sub_type, field_type);
	*left_paren = '[';
	dimension = strtol(left_paren + 1, NULL, 10);

	if ((dimension == LONG_MIN) || (dimension == LONG_MAX) ||
	    (dimension == 0)) {
	    if (iovar->var_array != TRUE) {
		fprintf(stderr, "Couldn't parse array size in \"%s\"\n",
			field_type);
		return;
	    } else {
		IOgetFieldStruct descr;  /* OK */
		long tmp_offset;
		descr.offset = iofield->field_offset;
		descr.size = ioformat->body->pointer_size;
		descr.data_type = integer_type;
		descr.byte_swap = byte_reversal;

		dimension = get_IOlong(iovar->dimens[0].control_field, data);

		tmp_offset = get_IOlong(&descr, data);
		if (encode) {
		    data = (char *) string_base + tmp_offset;
		} else {
		    data = (void *) tmp_offset;
		}
		offset = 0;
	    }
	    sub_field_size = iofield->field_size;
	} else {
	    /* normal internal array */
	    if (ioformat->body->IOversion <= 3) {
		sub_field_size = iofield->field_size / dimension;
	    } else {
		sub_field_size = iofield->field_size;
	    }
	}
	for (; dimension > 0; dimension--) {
	    add_value_as_XML(sub_type, sub_field_size, offset, iofile, data,
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
	int sub_field_size, offset = iofield->field_offset;
	dcatstr(ds, "\n");
	*left_paren = 0;
	strcpy(sub_type, field_type);
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
	if (ioformat->body->IOversion <= 3) {
	    sub_field_size = iofield->field_size / (dimension1 * dimension2);
	} else {
	    sub_field_size = iofield->field_size;
	}
	for (; dimension2 > 0; dimension2--) {
	    int i = 0;
	    for (; i < dimension1; i++) {
		add_value_as_XML(sub_type, sub_field_size, offset, iofile, 
				 data, string_base, byte_reversal, 
				 float_format, encode, 1, ds); 
		offset += sub_field_size;
	    }
	}
    }

    dcatstr(ds, "</");
    dcatstr(ds, (char*)iofield->field_name);
    dcatstr(ds, ">\n");
}

