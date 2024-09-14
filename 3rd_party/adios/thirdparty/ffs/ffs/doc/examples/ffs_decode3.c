#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
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

int main()     /* receiving program */
{
    FFSTypeHandle second_rec_handle;
    FFSContext context = create_FFSContext();
    second_rec *rec2;
    int fd, encode_size, decode_size, j;
    char *encoded_buffer;
    struct stat stat_buf;

    /* "receive" encoded record over a file */
    fd = open("enc_file", O_RDONLY, 0777);
    fstat(fd, &stat_buf);
    encoded_buffer = malloc(stat_buf.st_size);
    encode_size = read(fd, encoded_buffer, sizeof(encoded_buffer));

    second_rec_handle = FFSset_simple_target(context, "second format", second_field_list, sizeof(second_rec));
    FFS_target_from_encode(context, encoded_buffer);
    if (decode_in_place_possible(FFSTypeHandle_from_encode(context, encoded_buffer))) {
	FFSdecode_in_place(context, encoded_buffer, (void**)&rec2);
    } else {
	decode_size = FFS_est_decode_length(context, encoded_buffer, encode_size);
	rec2 = malloc(decode_size);
	FFSdecode(context, encoded_buffer, (void*)rec2);
    }
    printf("Read name=%s, count=%ld, items = [", rec2->name, rec2->count);
    for(j=0; j<rec2->count; j++) printf("%d, ", rec2->items[j]);
    printf("]\n");
}
