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

int
main(int argc, char** argv)
{
    FFSFile iofile = open_FFSfile("test_output", "r");
    FFSTypeHandle first_rec_handle;
    FFSContext context = FFSContext_of_file(iofile);
    first_rec rec1;

    first_rec_handle = FFSset_simple_target(context, "first format", field_list, sizeof(first_rec));
    FFSread(iofile, &rec1);
    close_FFSfile(iofile);
    printf("Read i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.j);
}
