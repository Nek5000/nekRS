#include "../config.h"
#include <assert.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <stdio.h>
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include "fm.h"

#include "test_funcs.h"

static int verbose = 0;

void
dump_format_list(FMStructDescList l)
{
    int i = 0;
    while (l[i].format_name != NULL) {
	int j = 0;
	FMFieldList fields = l[i].field_list;
	printf("Format name \"%s\", struct size %d\n",
	       l[i].format_name, l[i].struct_size);
	while(fields[j].field_name != NULL) {
	    printf("\t\t%s \"%s\"; size = %d; offset = %d\n",
		   fields[j].field_name,
		   fields[j].field_type,
		   fields[j].field_size,
		   fields[j].field_offset);
	    j++;
	}
	i++;
    }
}

extern void free_field_list(FMFieldList list);

void
do_test(FMStructDescList list)
{
    int i = 0;
    int failed = verbose;
    while (list[i].format_name != NULL) i++;
    FMStructDescList tmp = malloc(sizeof(tmp[0]) * (i+1));

    i = 0;
    while (list[i].format_name != NULL) {
	int j = 0;
	tmp[i].format_name = list[i].format_name;
	tmp[i].field_list = copy_field_list(list[i].field_list);
	tmp[i].struct_size = -1;
	tmp[i].opt_info = NULL;
	while(tmp[i].field_list[j].field_name != NULL) {
	    const char *typ = tmp[i].field_list[j].field_type;
	    tmp[i].field_list[j].field_offset = -1;
	    if (strncmp(typ, "integer", 7) == 0) {
		if (tmp[i].field_list[j].field_size == sizeof(int)) {
		    tmp[i].field_list[j].field_size = -1;
		}
	    } else if (strncmp(typ, "unsigned integer", 16) == 0) {
		if (tmp[i].field_list[j].field_size == sizeof(unsigned)) {
		    tmp[i].field_list[j].field_size = -1;
		}
	    } else if (strncmp(typ, "double", 6) == 0) {
		if (tmp[i].field_list[j].field_size == sizeof(double)) {
		    tmp[i].field_list[j].field_size = -1;
		}
	    } else if (strncmp(typ, "float", 5) == 0) {
		if (tmp[i].field_list[j].field_size == sizeof(double)) {
		    tmp[i].field_list[j].field_size = -1;
		}
	    } else if (strncmp(typ, "char", 4) == 0) {
	    } else {
		tmp[i].field_list[j].field_size = -1;
	    }
	    j++;
	}
	i++;
    }
    tmp[i].format_name = NULL;
    tmp[i].field_list = NULL;
    tmp[i].struct_size = 0;
    tmp[i].opt_info = NULL;
    

    FMlocalize_structs(tmp);
    i = 0;
    while (tmp[i].format_name != NULL) {
	int j = 0;
	FMFieldList fields = tmp[i].field_list;
	while(fields[j].field_name != NULL) {
	    if (fields[j].field_offset != list[i].field_list[j].field_offset) {
		printf("Failed alignment of field %s in format %s, got %d, expected %d\n",
		       fields[j].field_name, list[i].format_name, 
		       fields[j].field_offset,
		       list[i].field_list[j].field_offset);
		failed++;
	    }
	    if (fields[j].field_size != list[i].field_list[j].field_size) {
		printf("Failed size of field %s in format %s, got %d, expected %d\n",
		       fields[j].field_name, list[i].format_name, 
		       fields[j].field_size,
		       list[i].field_list[j].field_size);
		failed++;
	    }
	    j++;
	}
	if (tmp[i].struct_size != list[i].struct_size) {
	    printf("Struct sizes for \"%s\" differ, got %d, expected %d\n",
		   tmp[i].format_name, tmp[i].struct_size, 
		   list[i].struct_size);
	    failed++;
	}
	i++;
    }
    if (failed) {
	i = 0;
	printf("Results were:\n");
	dump_format_list(tmp);
	printf("\nExpected : \n");
	dump_format_list(list);
    }
    i = 0;
    while (tmp[i].format_name != NULL) {
	free_field_list(tmp[i].field_list);
	i++;
    }
    free(tmp);
}


int
main(int argc, char **argv)
{
    if (argc > 1) {
	if (strcmp(argv[1], "-v") == 0) {
	    verbose++;
	}
    }
    
    do_test(first_format_list);
    do_test(string_format_list);
    do_test(structured_format_list);
    do_test(two_string_format_list);
    do_test(fourth_format_list);
    do_test(later_format_list);
    do_test(nested_format_list);
    do_test(embedded_format_list);
    do_test(variant_format_list);
    do_test(string_array_format_list);
    do_test(derive_format_list);
    do_test(multi_array_format_list);
    do_test(triangle_format_list);
    do_test(add_action_format_list);
    do_test(node_format_list);
    return 0;
}
