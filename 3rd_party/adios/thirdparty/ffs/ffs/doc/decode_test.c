#include <fcntl.h>
#include <stdlib.h>
#include "io.h"

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

void main()     /* receiving program */
{
    IOContext dest_context = create_IOcontext();
    IOFormat dyn_rec_ioformat;
    dyn_rec rec;
    int fd, i;
    char encoded_buffer[2048];  /* hopefully big enough */

    /* "receive" encoded record over a file */
    fd = open("enc_file", O_RDONLY, 0777);
    read(fd, encoded_buffer, sizeof(encoded_buffer));

    dyn_rec_ioformat = get_format_IOcontext(dest_context, encoded_buffer);
    set_conversion_IOcontext(dest_context, dyn_rec_ioformat,
			     dyn_field_list, sizeof(dyn_rec));
    decode_IOcontext(dest_context, encoded_buffer, &rec);
    printf("string is %s\n", rec.string);
    printf("icount is %d\n", rec.icount);
    for(i=0; i< rec.icount; i++)
	printf("element %d is %g\n", i, rec.double_array[i]);
}
