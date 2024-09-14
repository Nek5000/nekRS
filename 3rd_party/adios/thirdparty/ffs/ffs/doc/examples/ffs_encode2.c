#include <fcntl.h>
#include <stdlib.h>
#include "ffs.h"

typedef struct _second_rec {
    int         *items;
    long        count;
    char *      name;
} second_rec, *second_rec_ptr;

static FMField second_field_list[] = {
    {"items", "integer[count]", sizeof(int), FMOffset(second_rec_ptr, items)},
    {"count", "integer", sizeof(long), FMOffset(second_rec_ptr, count)},
    {"name", "string",   sizeof(char*), FMOffset(second_rec_ptr, name)},
    {NULL, NULL, 0, 0},
};

int main(int argc, char **argv)
{
    FMContext fmc = create_FMcontext();
    FFSBuffer buf = create_FFSBuffer();
    FMFormat rec_ioformat;
    second_rec rec2;
    FFSEncodeVector outvec;
    int fd, outvec_size;
    char str[32];
    int array[10];
    int j;

    srand48(time());
    sprintf(str, "random string %ld", lrand48()%100);
    rec2.name = &str[0];
    rec2.count = lrand48()%10;
    rec2.items = &array[0];
    for(j=0; j<rec2.count; j++) rec2.items[j] = lrand48()%5+2;
    rec_ioformat = register_simple_format(fmc, "second format", second_field_list, sizeof(second_rec));
    outvec = FFSencode_vector(buf, rec_ioformat, &rec2);
    for(outvec_size=0;outvec[outvec_size].iov_len!=0;outvec_size++);

    printf("Wrote name=%s, count=%ld, items = [", rec2.name, rec2.count);
    for(j=0; j<rec2.count; j++) printf("%d, ", rec2.items[j]);
    printf("]\n");

    /* write the encoded data */
    fd = open("enc_file2", O_WRONLY|O_CREAT|O_TRUNC, 0777);
    writev(fd, outvec, outvec_size);
    close(fd);
}

