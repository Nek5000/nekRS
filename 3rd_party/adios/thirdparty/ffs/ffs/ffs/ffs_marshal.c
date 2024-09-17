#include "config.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <ffs.h>
#include "cod.h"
#include "fm_internal.h"
#include "ffs_internal.h"
#include "../cod/cod_internal.h"
#include "ffs_marshal.h"

extern field_marshal_info
get_marshal_info(FMFormat f, FMTypeDesc *t)
{
    format_marshal_info info = (format_marshal_info) f->ffs_info;
    int i;
    if (info == NULL) return NULL;
    for (i=0; i < info->count; i++) {
	if (info->field_info[i].t == t) {
	    return &info->field_info[i];
	}
    }
    return NULL;
}

#ifdef DO_DCG
extern sm_ref
cod_build_type_node(const char *name, FMFieldList field_list);
extern sm_ref
cod_build_param_node(const char *id, sm_ref typ, int param_num);

static void
add_param(cod_parse_context parse_context, char *name, int param_num,
	  FMFormat format)
{
    FMStructDescList list = format_list_of_FMFormat(format);
    int i = 1;
    sm_ref type, param;
    while (list[i].format_name != NULL) {
	FMFieldList fl = list[i].field_list;
	/* step through input formats */
	cod_add_simple_struct_type(list[i].format_name, fl, parse_context);
	i++;
    }
    type = cod_build_type_node(list[0].format_name, list[0].field_list);
    cod_add_decl_to_parse_context(list[0].format_name, type, parse_context);

    param = cod_build_param_node(name, type, param_num);

    cod_add_decl_to_parse_context(name, param, parse_context);
}

field_marshal_info
add_marshal_info(FMFormat f)
{
    format_marshal_info info = (format_marshal_info) f->ffs_info;
    if (info == NULL) {
	f->ffs_info = info = malloc(sizeof(struct marshal_info));
	info->count = 1;
	info->field_info = malloc(sizeof(info->field_info[0]));
    } else {
	info->count++;
	info->field_info = realloc(info->field_info,
				   sizeof(info->field_info[0]) * info->count);
    }
    return &info->field_info[info->count - 1];
}

extern void
install_drop_code(FMFormat f, char *field, char*code_str)
{
    cod_code code;
    int (*func)(void*);
    field_marshal_info marshal_info;
    cod_parse_context parse_context = new_cod_parse_context();
    int i, field_num = -1;

    static char extern_string[] = "\
		int printf(string format, ...);\n\
		void *malloc(int size);\n\
		void free(void *pointer);";

    static cod_extern_entry externs[] = {
	{"printf", (void *) 0},
	{"malloc", (void*) 0},
	{"free", (void*) 0},
	{(void *) 0, (void *) 0}
    };

    /*
     * some compilers think it isn't a static initialization to put this
     * in the structure above, so do it explicitly.
     */
    externs[0].extern_value = (void *) (intptr_t) printf;
    externs[1].extern_value = (void *) (intptr_t) malloc;
    externs[2].extern_value = (void *) (intptr_t) free;


    for (i=0; i< f->field_count; i++) {
	if (strcmp(f->field_list[i].field_name, field) == 0) field_num = i;
    }
    if (field_num == -1) {
	printf("field \"%s\" not found in install drop code\n", field);
	return;
    }
    add_param(parse_context, "input", 0, f);
    cod_assoc_externs(parse_context, externs);
    cod_parse_for_context(extern_string, parse_context);

    code = cod_code_gen(code_str, parse_context);
    cod_free_parse_context(parse_context);
    if (code == NULL) {
	printf("Compilation failed, field \"%s\" in install drop code \n", field);
	return;
    }
    func = (int(*)(void *))code->func;
    marshal_info = add_marshal_info(f);
    marshal_info->t = &f->var_list[field_num].type_desc;
    marshal_info->type = FFSDropField;
    marshal_info->drop_row_func = func;
}


static void
copy_array_element(cod_exec_context ec, int element)
{
    struct subsample_marshal_data *smd = (void*)cod_get_client_data(ec, 0x534d4450);
    char *element_src, *element_dest;
    if ((element < 0) || (element >= smd->element_count)) return;
    if (smd->marshalled_count == smd->element_count) {
	printf("Already Marshalled %d elements of this array!\n", 
	       smd->element_count);
	return;
    }
    element_src = (char*)smd->src_ptr + (element * smd->element_size);
    element_dest = (char*)smd->dst_ptr + (smd->marshalled_count * smd->element_size);
    memcpy(element_dest, element_src, smd->element_size);
    smd->marshalled_count++;
}

extern void
install_subsample_code(FMFormat f, char *field, char*code_str)
{
    cod_code code;
    int (*func)(cod_exec_context, void *, int);
    field_marshal_info marshal_info;
    cod_parse_context parse_context = new_cod_parse_context();
    int i, field_num = -1;

    static char extern_string[] = "\
		int printf(string format, ...);\n\
		void *malloc(int size);\n\
		void memcpy(void* dest, void* src, int size);\n\
        void FFSMarshalArrayElement(cod_exec_context ec, int element);";

    static cod_extern_entry externs[] = {
	{"printf", (void *) 0},
	{"malloc", (void*) 0},
	{"memcpy", (void*) 0},
	{"FFSMarshalArrayElement", (void*) 0},
	{(void *) 0, (void *) 0}
    };

    /*
     * some compilers think it isn't a static initialization to put this
     * in the structure above, so do it explicitly.
     */
    externs[0].extern_value = (void *) (intptr_t) printf;
    externs[1].extern_value = (void *) (intptr_t) malloc;
    externs[2].extern_value = (void *) (intptr_t) memcpy;
    externs[3].extern_value = (void *) (intptr_t) copy_array_element;


    for (i=0; i< f->field_count; i++) {
	if (strcmp(f->field_list[i].field_name, field) == 0) field_num = i;
    }
    if (field_num == -1) {
	printf("field \"%s\" not found in install subsample code\n", field);
	return;
    }
    cod_add_param("ec", "cod_exec_context", 0, parse_context);
    add_param(parse_context, "input", 1, f);
    cod_add_param("element_count", "int", 2, parse_context);

    cod_set_return_type("void", parse_context);
    cod_assoc_externs(parse_context, externs);
    cod_parse_for_context(extern_string, parse_context);

    code = cod_code_gen(code_str, parse_context);
    cod_free_parse_context(parse_context);
    if (code == NULL) {
	printf("Compilation failed, field \"%s\" in install subsample code \n", field);
	return;
    }
    func = (int(*)(cod_exec_context, void *, int))code->func;
    marshal_info = add_marshal_info(f);
    marshal_info->ec = cod_create_exec_context(code);
    marshal_info->t = &f->var_list[field_num].type_desc;
    marshal_info->type = FFSSubsampleArrayField;
    marshal_info->subsample_array_func = func;
}


#else
extern void
install_drop_code(FMFormat f, char *field, char*code_str){}
extern void
install_subsample_code(FMFormat f, char *field, char*code_str){}
#endif
