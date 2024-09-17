#include "../../ffs.h"
typedef struct _first_rec {
    int         i;
    long        j;
    double      d;
    char        c;
} first_rec, *first_rec_ptr;

static FMField field_list[] = {
    {"i", "integer", sizeof(int), FMOffset(first_rec_ptr, i)},
    {"j", "integer", sizeof(long), FMOffset(first_rec_ptr, j)},
    {"d", "float",   sizeof(double), FMOffset(first_rec_ptr, d)},
    {"c", "integer", sizeof(char), FMOffset(first_rec_ptr, c)},
    {NULL, NULL, 0, 0},
};

typedef struct _second_rec {
    int         *items;
    long        count;
    char *      name;
} second_rec, *second_rec_ptr;

static FMField field_list2[] = {
    {"items", "integer[count]", sizeof(int), FMOffset(second_rec_ptr, items)},
    {"count", "integer", sizeof(long), FMOffset(second_rec_ptr, count)},
    {"name", "string",   sizeof(char*), FMOffset(second_rec_ptr, name)},
    {NULL, NULL, 0, 0},
};

int
main(int argc, char** argv)
{
    FFSFile iofile = open_FFSfile("test_output", "r");
    FFSTypeHandle first_rec_handle, second_rec_handle, next_handle;
    FFSContext context = FFSContext_of_file(iofile);

    first_rec_handle = FFSset_simple_target(context, "first format", field_list, sizeof(first_rec));
    second_rec_handle = FFSset_simple_target(context, "second format", field_list2, sizeof(second_rec));
    while ((next_handle = FFSnext_type_handle(iofile)) != NULL) {
	if (next_handle == first_rec_handle) {
	    first_rec rec1;
	    FFSread(iofile, &rec1);
	    printf("Read first_rec : i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.c);
	} else 	if (next_handle == second_rec_handle) {
	    second_rec rec2;
	    int i;
	    FFSread(iofile, &rec2);
	    printf("Read items= ");
	    for (i=0; i<rec2.count; i++) printf("%d ", rec2.items[i]);
	    printf(", count = %d, name = %s\n", rec2.count, rec2.name);
	} else {
	    printf("Unknown record type file\n");
	}
    }
    close_FFSfile(iofile);
}
