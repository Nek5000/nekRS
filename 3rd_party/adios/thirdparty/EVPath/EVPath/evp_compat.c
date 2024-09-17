
#include "config.h"
#include <memory.h>
#undef NDEBUG
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "evpath.h"
#include "cm_internal.h"

struct _CMformat_list {
    /*! the name to be associated with this structure */
    char *format_name;
    /*! the PBIO-style list of fields within this structure */
    FMFieldList field_list;
};
typedef struct _CMformat_list CMFormatRec;
typedef CMFormatRec *CMFormatList;

extern int
IOget_array_size_dimen(const char *str, FMFieldList fields, int dimen,
		       int *control_field);

static int
is_var_array_field(FMFieldList field_list, int field)
{
    int done = 0;
    int ret = 0;
    int dimen_count = 0;
    int control_val;
    while (!done) {
	int static_size = IOget_array_size_dimen(field_list[field].field_type,
						 field_list, dimen_count, 
						 &control_val);
	dimen_count++;
	if (static_size == 0) {
	    done++;
	    continue;
	}
	if ((static_size == -1) && (control_val == -1)) {
	    /* failed validation, errors already delivered */
	    return -1;
	}
	if (control_val != -1) {
	    /* dynamic array */
	    ret = 1;
	}
    }
    return ret;
}

#define Max(i,j) ((i<j) ? j : i)

extern FMdata_type FMarray_str_to_data_type(const char *str, 
					    long *element_count_ptr);
static
int
struct_size_field_list(FMFieldList list, int pointer_size)
{
    int i = 0;
    int struct_size = 0;
    while (list[i].field_name != NULL) {
	int field_size = 0;
	if (is_var_array_field(list, i) == 1) {
	    /* variant array, real_field_size is ioformat->pointer_size */
	    field_size = pointer_size;
	} else {
	    long elements;
	    FMarray_str_to_data_type(list[i].field_type, &elements);
	    field_size = list[i].field_size * elements;
	}
	assert(field_size > 0);
	struct_size = Max(struct_size,
			  (list[i].field_offset + field_size));
/*	printf("i=%d field_name=%s field_type=%s struct_size= %d, offset=%d size=%d\n", i, list[i].field_name, list[i].field_type, struct_size, list[i].field_offset, field_size);*/
	i++;
    }
    return struct_size;
}

extern CMFormat
old_CMregister_format(CManager cm, char *format_name,
		      FMFieldList field_list, CMFormatList subformat_list)
{
    FMStructDescList structs;
    int sub_count = 0, i;
    if (subformat_list && (subformat_list[sub_count].format_name != NULL)) sub_count++;
    structs = malloc(sizeof(structs[0]) * (sub_count + 2));
    structs[0].format_name = format_name;
    structs[0].field_list = field_list;
    structs[0].struct_size = struct_size_field_list(field_list, (int)sizeof(char*));
    structs[0].opt_info = NULL;
    for (i = 0; i < sub_count; i++) {
	structs[i+1].format_name = subformat_list[i].format_name;
	structs[i+1].field_list = subformat_list[i].field_list;
	structs[i+1].struct_size = struct_size_field_list(structs[i+1].field_list, (int)sizeof(char*));
	structs[i+1].opt_info = NULL;
    }
    structs[sub_count+1].format_name = NULL;
    structs[sub_count+1].field_list = NULL;
    return CMregister_format(cm, structs);
}

extern EVaction
old_EVassoc_terminal_action(CManager cm, EVstone stone, CMFormatList format_list, 
			EVSimpleHandlerFunc handler, void* client_data)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return EVassoc_terminal_action(cm, stone, structs, handler, client_data);
}

extern EVaction
old_EVassoc_filter_action(CManager cm, EVstone stone, CMFormatList format_list, 
			EVSimpleHandlerFunc handler, EVstone out_stone, void* client_data)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return EVassoc_filter_action(cm, stone, structs, handler, out_stone, client_data);
}

extern EVsource
old_EVcreate_submit_handle(CManager cm, EVstone stone, CMFormatList format_list)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return EVcreate_submit_handle(cm, stone, structs);
}

extern EVsource
old_EVcreate_submit_handle_free(CManager cm, EVstone stone, CMFormatList format_list,
			    EVFreeFunction free_func, void *client_data)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return EVcreate_submit_handle_free(cm, stone, structs, free_func, client_data);
}

extern char *
old_create_filter_action_spec(CMFormatList format_list, char *function)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return create_filter_action_spec(structs, function);
}

extern char *
old_create_router_action_spec(CMFormatList format_list, char *function)
{
    FMStructDescList structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;
    return create_router_action_spec(structs, function);
}

extern char *
old_create_transform_action_spec(CMFormatList format_list, CMFormatList out_format_list, char *function)
{
    FMStructDescList structs, out_structs;
    int count = 0, i;
    while (format_list && (format_list[count].format_name != NULL)) count++;
    structs = malloc(sizeof(structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	structs[i].format_name = format_list[i].format_name;
	structs[i].field_list = format_list[i].field_list;
	structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	structs[i].opt_info = NULL;
    }
    structs[count].format_name = NULL;
    structs[count].field_list = NULL;

    count = 0;
    while (out_format_list && (out_format_list[count].format_name != NULL)) count++;
    out_structs = malloc(sizeof(out_structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	out_structs[i].format_name = out_format_list[i].format_name;
	out_structs[i].field_list = out_format_list[i].field_list;
	out_structs[i].struct_size = struct_size_field_list(out_structs[i].field_list, (int)sizeof(char*));
	out_structs[i].opt_info = NULL;
    }
    out_structs[count].format_name = NULL;
    out_structs[count].field_list = NULL;
    return create_transform_action_spec(structs, out_structs, function);
}

extern char*
old_create_multityped_action_spec(CMFormatList *input_format_lists, CMFormatList out_format_list, char *function)
{
    FMStructDescList structs, out_structs, *struct_list;
    int struct_count = 0, j;
    int count, i;
    while (input_format_lists[struct_count] != NULL) struct_count++;
    struct_list = malloc(sizeof(struct_list[0]) * (struct_count+1));
    for (j = 0; j < struct_count ; j++) {
	CMFormatList format_list = input_format_lists[j];
	count = 0;
	while (format_list && (format_list[count].format_name != NULL)) count++;
	structs = malloc(sizeof(structs[0]) * (count + 1));
	for (i = 0; i < count; i++) {
	    structs[i].format_name = format_list[i].format_name;
	    structs[i].field_list = format_list[i].field_list;
	    structs[i].struct_size = struct_size_field_list(structs[i].field_list, (int)sizeof(char*));
	    structs[i].opt_info = NULL;
	}
	structs[count].format_name = NULL;
	structs[count].field_list = NULL;
	struct_list[j] = structs;
    }

    count = 0;
    while (out_format_list && (out_format_list[count].format_name != NULL)) count++;
    out_structs = malloc(sizeof(out_structs[0]) * (count + 1));
    for (i = 0; i < count; i++) {
	out_structs[i].format_name = out_format_list[i].format_name;
	out_structs[i].field_list = out_format_list[i].field_list;
	out_structs[i].struct_size = struct_size_field_list(out_structs[i].field_list, (int)sizeof(char*));
	out_structs[i].opt_info = NULL;
    }
    out_structs[count].format_name = NULL;
    out_structs[count].field_list = NULL;
    return create_multityped_action_spec(struct_list, function);
}

