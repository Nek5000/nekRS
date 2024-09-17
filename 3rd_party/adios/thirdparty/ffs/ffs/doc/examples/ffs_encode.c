#include <fcntl.h>
#include <stdlib.h>
#include "ffs.h"

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
    FMContext fmc = create_FMcontext();
    FFSBuffer buf = create_FFSBuffer();
    FMFormat rec_ioformat;
    first_rec rec1;
    char *output;
    int fd, output_size;

    srand48(time());
    rec1.i = lrand48()%100;
    rec1.j = lrand48()%100;
    rec1.d = drand48();
    rec1.j = lrand48()%26+'A';
    rec_ioformat = register_simple_format(fmc, "first format", field_list, sizeof(first_rec));
    output = FFSencode(buf, rec_ioformat, &rec1, &output_size);

    /* write the encoded data */
    fd = open("enc_file", O_WRONLY|O_CREAT|O_TRUNC, 0777);

    printf("Wrote i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.j);
    write(fd, output, output_size);
    close(fd);
}

