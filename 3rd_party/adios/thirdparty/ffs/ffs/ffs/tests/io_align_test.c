#include "config.h"
#include <assert.h>
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include "ffs.h"

#include "test_funcs.h"

static int verbose = 0;

extern void 
local_align_field_list(FMFieldList field_list, int pointer_size);

void
do_test(FMContext c, char *name, FMFieldList fields)
{
    int i = 0;
    FMFieldList tmp = copy_field_list(fields);

    local_align_field_list(tmp, sizeof(char*));
    while(fields[i].field_name != NULL) {
	if (fields[i].field_offset != tmp[i].field_offset) {
	    printf("Failed alignment of field %s in format %s, got %d, expected %d\n",
		   fields[i].field_name, name, tmp[i].field_offset,
		   fields[i].field_offset);
	}
	i++;
    }
}


int
main(argc, argv)
int argc;
char **argv;
{
    FMContext c;
    do_test(c, "first format", field_list);
    do_test(c, "string format", field_list2);
    do_test(c, "two string format", field_list3);
    do_test(c, "internal array format", field_list4);
    do_test(c, "embedded", embedded_field_list);
    do_test(c, "structured array format", field_list5);
    do_test(c, "variant array format", field_list6);
    do_test(c, "later format", later_field_list);
    do_test(c, "nested format", nested_field_list);
    do_test(c, "EventVecElem", event_vec_elem_fields);
    do_test(c, "EventV", field_list9);
    do_test(c, "string_array", string_array_field_list);
    return 0;
}
