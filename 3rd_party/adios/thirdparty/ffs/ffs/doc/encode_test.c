#include <fcntl.h>
#include <stdlib.h>
#include "ffs.h"

typedef struct _dyn_rec {
    char	*string;
    long	icount;
    double	*double_array;
} dyn_rec, *dyn_rec_ptr;

IOField dyn_field_list[] = {
    {"string field", "string", sizeof(char *), 
      IOOffset(dyn_rec_ptr, string)},
    {"icount", "integer", sizeof(long), 
      IOOffset(dyn_rec_ptr, icount)},
    {"double_array", "float[icount]", sizeof(double), 
      IOOffset(dyn_rec_ptr, double_array)},
    { NULL, NULL, 0, 0}
};

int main(int argc, char **argv)
{
    FMContext fmc = create_FMContext();
    FFSBuffer buf = create_FFSBuffer();
    FMFormat rec_ioformat;
    first_rec rec1;
    char *output;
    int fd, output_size;

    rec_ioformat = register_simple_format(fmc, "first format", field_list, sizeof(first_rec));
    output = FFSencode(buf, rec_ioformat, &rec1, &output_size);

    /* write the encoded data */
    fd = open("enc_file", O_WRONLY|O_CREAT|O_TRUNC, 0777);
    write(fd, output, output_size);
    close(fd);
}

