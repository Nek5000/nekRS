#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
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

int main()     /* receiving program */
{
    FFSTypeHandle first_rec_handle;
    FFSContext context = create_FFSContext();
    first_rec rec1;
    int fd, size;
    char *encoded_buffer;
    struct stat stat_buf;

    /* "receive" encoded record over a file */
    fd = open("enc_file", O_RDONLY, 0777);
    fstat(fd, &stat_buf);
    encoded_buffer = malloc(stat_buf.st_size);
    size = read(fd, encoded_buffer, stat_buf.st_size);

    first_rec_handle = FFSset_simple_target(context, "first format", field_list, sizeof(first_rec));
    FFS_target_from_encode(context, encoded_buffer);
    FFSdecode(context, encoded_buffer, (void*)&rec1);
    printf("Read i=%d, j=%ld, d=%g, j=%c\n", rec1.i, rec1.j, rec1.d, rec1.j);
}
