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

int main(int argc, char **argv)
{
    FFSFile file = open_FFSfile("test_output", "w");
    FMContext fmc = FMContext_of_file(file);
    FMFormat rec_ioformat;
    first_rec rec1;

    srand48(time());
    rec1.i = lrand48()%100;
    rec1.j = lrand48()%100;
    rec1.d = drand48();
    rec1.j = lrand48()%26+'A';
    rec_ioformat = register_simple_format(fmc, "first format", field_list, sizeof(first_rec));
    write_FFSfile(file, rec_ioformat, &rec1);
    close_FFSfile(file);
    printf("Wrote i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.j);
}
