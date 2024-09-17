#include "config.h"
#include "ffs.h"
#include <sys/types.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include "ffs_internal.h"

typedef struct struct_element {
    IOFormat format;
    void *base;
    int decrement;
} *struct_stack_p;

typedef struct _IOSaxParseInfo {
    int cur_depth;
    int stack_depth;
    IOContext c;
    struct_stack_p struct_stack;
    IOgetFieldStruct next_field;
    int verbose;
    IOSaxEndItemCallback end_item_callback;
    void *orig_base;
} *parser_info;

static void get_data_from_string(parser_info info, const char *st, int len);
extern void IOSaxEndElement(void *userData, const char *name);

void IOSaxStartElement(void *userData, const char *name, const char **atts)
{
    parser_info info = userData;
    IOFormat this_field_format = NULL;
    int size = 0, offset = 0, is_var_array = 0;
    IOdata_type base_data_type = unknown_type;
    int i;

    info->next_field.size = info->next_field.offset = 0;
    if ((info->cur_depth) == info->stack_depth) {
	/* make sure the stack is deep enough */
	info->struct_stack = realloc(info->struct_stack, 
				     sizeof(*info->struct_stack) * 
				     (info->cur_depth+3));
	info->struct_stack[info->cur_depth].format = NULL;
	info->struct_stack[info->cur_depth].base = NULL;
	info->struct_stack[info->cur_depth].decrement = 1;
	info->stack_depth = info->cur_depth + 1;
    }
    if (info->verbose) {
	int i;
	for(i = 0 ; i< info->cur_depth; i++) printf("  ");
	printf("<%s>\n", name);
    }
    if (info->cur_depth == 0) {
	this_field_format = get_IOformat_by_name_IOcontext(info->c, name);
	if (this_field_format == NULL) {
	    printf("Outermost format %s not found\n", name);
	}
    } else {
	/* check containing format */
	IOFormat f = info->struct_stack[info->cur_depth-1].format;
	for (i=0; i< f->body->field_count ; i++) {
	    if (strcmp(f->body->field_list[i].field_name, name) == 0) {
		IOFieldList fl = &f->body->field_list[i];
		IOVarInfoList var = &f->body->var_list[i];
		base_data_type = var->data_type;
		size = fl->field_size;
		offset = fl->field_offset;
		is_var_array = var->var_array;
		if (base_data_type == unknown_type) {
		    /* moving into a subformat */
		    this_field_format = 
			info->struct_stack[info->cur_depth].format = 
			f->field_subformats[i];
		    info->struct_stack[info->cur_depth].base =
			(char*)info->struct_stack[info->cur_depth-1].base 
			+ offset;
		}
		if (is_var_array) {
		    void *containing_base = 
			info->struct_stack[info->cur_depth-1].base;
		    long array_dimen = 
			FMget_array_element_count(f, var, containing_base, 0);
		    IOgetFieldStruct descr;
		    void *pointer;
		    descr.offset = fl->field_offset;
		    descr.size = f->body->pointer_size;
		    descr.data_type = integer_type;
		    descr.byte_swap = 0;
		    /* get the pointer */
		    pointer = (void*)get_IOlong(&descr, containing_base);
		    if (array_dimen == 0) {
			pointer = malloc(size);
			array_dimen = 1;
		    } else {
			array_dimen++;
			pointer = realloc(pointer, size * array_dimen);
		    }
		    put_IOlong(&descr, (long)pointer, containing_base);
		    put_IOlong(var->dimens[0].control_field, array_dimen,
			       containing_base);
		    info->struct_stack[info->cur_depth].base = 
			((char*)pointer) + ((array_dimen-1) * fl->field_size);
		}
		break;
	    }
	}
	if (size == 0) {
	    printf("Field not found %s depth %d\n", name, info->cur_depth);
	}
	info->next_field.offset = offset;
	info->next_field.size = size;
	info->next_field.data_type = base_data_type;
	info->next_field.byte_swap = 0;
    }
    info->struct_stack[info->cur_depth].format = this_field_format;
    if (this_field_format != NULL) {
	/* this is a structure we're starting, init some stuff */
	memset(info->struct_stack[info->cur_depth].base, 0, 
	       this_field_format->body->record_length);
#if 0
	for (i=0; i< this_field_format->body->field_count ; i++) {
	    if (this_field_format->body->var_list[i].var_array) {
		/* 
		 * got a variant array here, null the pointer and 
		 * zero the count
		 */
		IOFieldList fl = &this_field_format->body->field_list[i];
		IOVarInfoList var = &this_field_format->body->var_list[i];
		IOgetFieldStruct descr;
		descr.offset = fl->field_offset;
		descr.size = this_field_format->body->pointer_size;
		descr.data_type = integer_type;
		descr.byte_swap = 0;
		/* zero the pointer */
		put_IOlong(&descr, 0, 
			   info->struct_stack[info->cur_depth].base);
		/* zero the count/length value */
		put_IOlong(var->control_field, 0, 
			   info->struct_stack[info->cur_depth].base);
	    }
	}
#endif
    }
    info->cur_depth += 1;
    while ((atts != NULL) && (atts[0] != NULL)) {
	char *new_name = malloc(strlen(name) + strlen(atts[0]) + 2);
	strcpy(new_name, name);
	strcat(new_name, "_");
	strcat(new_name, atts[0]);
	IOSaxStartElement(info, new_name, NULL);
	get_data_from_string(info, atts[1], strlen(atts[1]));
	IOSaxEndElement(info, new_name);
	free(new_name);
	atts += 2;
    }
    if (is_var_array && (base_data_type != unknown_type)) {
	info->cur_depth += 1;
	info->struct_stack[info->cur_depth].decrement = 2;
	info->next_field.offset = 0;
    } else {
	info->struct_stack[info->cur_depth].decrement = 1;
    }
}

extern void
IOset_parse_callback(IOSaxParseInfo info, IOSaxEndItemCallback callback)
{
    info->end_item_callback = callback;
}

extern void 
IOSaxEndElement(void *userData, const char *name)
{
    parser_info info = userData;
    info->cur_depth -= info->struct_stack[info->cur_depth].decrement;
    info->struct_stack[info->cur_depth].format = NULL;
    info->struct_stack[info->cur_depth].decrement = 1;
    info->next_field.size = 0;
    if (info->verbose) {
	int i;
	for(i = 0 ; i< info->cur_depth; i++) printf("  ");
	printf("</%s>\n", name);
    }
    if ((info->cur_depth == 0) && (info->end_item_callback != NULL)) {
	(info->end_item_callback)(name, info->orig_base);
	info->struct_stack[0].format = 0;
	info->struct_stack[0].base = info->orig_base;
	info->struct_stack[1].format = 0;
	info->struct_stack[1].base = info->orig_base;
    }
}

static void
get_data_from_string(parser_info info, const char *st, int len)
{
    switch(info->next_field.data_type) {
    case integer_type:
    { 
	long value = 0;
/*	printf(" would get int from %s\n", st);*/
	sscanf(st, "%ld", &value);
/*	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOlong(&info->next_field, value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case unsigned_type:
    {
	unsigned long value = 0;
/*	printf(" would get uint from %s\n", st);*/
	sscanf(st, "%lu", &value);
/*	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOlong(&info->next_field, value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case float_type:
    {
	double value = 0.0;
/*	printf(" would get float from %s\n", st);*/
	sscanf(st, "%lg", &value);
/*	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOdouble(&info->next_field, value,
		     info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case char_type:
    {
	unsigned char value = *st;
/*	printf(" would get char %c\n", *st);
	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOchar(&info->next_field, value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case string_type:
    {
	char *value = NULL;
/*	printf(" would get string from %s\n", st);*/
	if (len != 0) {
	    value = malloc(len + 1);
	    strcpy(value, st);
	}
	put_IOlong(&info->next_field, (long)value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case enumeration_type:
    {
	unsigned long value = 0;
/*	printf("would get enum from %s\n", st);*/
	sscanf(st, "%lu", &value);
/*	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOlong(&info->next_field, value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case boolean_type:
    { 
	int value;
	if (strcmp(st, "false") == 0) {
	    value = 0;
	} else {
	    value = 1;
	}
/*	printf("would get boolean from %s, %d\n", st, value);
	printf("put to addr %lx, offset %lx\n", 
	       info->struct_stack[info->cur_depth -2].base,
	       info->next_field.offset);*/
	put_IOenum(&info->next_field, value,
		   info->struct_stack[info->cur_depth -2].base);
	break;
    }
    case unknown_type:
/*	printf("UNKNOWN %s\n", st);*/
	break;
    }
}

void IOSaxDataHandler(void *userData, const char *s, int len)
{
    parser_info info = userData;
    char tmp[256];
    char *st;
    if (info->next_field.size == 0) return;

    if (len < 255) {
	strncpy(tmp, s, len);
	tmp[len] = 0;
	st = &tmp[0];
    } else {
	st = malloc(len + 1);
	strncpy(st, s, len);
    }
    get_data_from_string(info, st, len);
    if (st != &tmp[0]) free(st);
}

extern IOSaxParseInfo
IOcreate_SaxParseInfo(IOContext c, void *init_base)
{
    IOSaxParseInfo ret = malloc(sizeof(struct _IOSaxParseInfo));
    ret->struct_stack = malloc(sizeof(*ret->struct_stack));
    ret->struct_stack[0].format = 0;
    ret->struct_stack[0].base = init_base;
    ret->c = c;
    ret->stack_depth = 1;
    ret->cur_depth = 0;
    ret->verbose = (getenv("SaxVerbose") != NULL);
    ret->end_item_callback = (IOSaxEndItemCallback) NULL;
    ret->orig_base = init_base;
    return ret;
}
