// This is a copied program from ffs_write.c

#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "ffs.h"

typedef struct _second_rec {
    int len;
    double *array;
    char char_field;
} second_rec, *second_rec_ptr;

FMField field_list2[] = {
    {"len", "integer", 
       sizeof(int), FMOffset(second_rec_ptr, len)},
    {"array", "float[len]",
       sizeof(double), FMOffset(second_rec_ptr, array)},
    {"char field", "char",
       sizeof(char), FMOffset(second_rec_ptr, char_field)},
    { NULL, NULL, 0, 0}
};

extern void
install_drop_code(FMFormat f, char *field, char*code_str);

int
main(int argc, char **argv)
{
    FMContext src_context;
    FFSFile ffsfile;
    FMFormat second_rec_ioformat;

    second_rec rec2;

    FMStructDescRec str_list[5];
    char *output_file = "marshal_output";
    int i;

    ffsfile = open_FFSfile(output_file, "w");

    src_context = create_local_FMcontext();

    str_list[0].format_name = "string format";
    str_list[0].field_list = field_list2;
    str_list[0].struct_size = sizeof(second_rec);
    str_list[0].opt_info = NULL;
    str_list[1].format_name = NULL;
    second_rec_ioformat = register_data_format(src_context, str_list);

    memset((char *) &rec2, 0, sizeof(rec2));
    rec2.len = 10;
    rec2.array = malloc(sizeof(rec2.array[0]) * rec2.len);
    for (i=0; i < rec2.len; i++ ) {
	rec2.array[i] = i * 2.0 + 0.5;
    }
    rec2.char_field = 'A';
    rec2.array[0] = 0.5;
    install_drop_code(second_rec_ioformat, "array", "{if (input.array[0] == 0.5) return 1; return 0;}");
    if (!write_FFSfile(ffsfile, second_rec_ioformat, &rec2))
	printf("write failed\n");

    rec2.array[0] = 2.0 + .5;

    if (!write_FFSfile(ffsfile, second_rec_ioformat, &rec2))
	printf("write failed\n");

    close_FFSfile(ffsfile);
    free_FFSfile(ffsfile);

    ffsfile = open_FFSfile(output_file, "r");
    FFSset_fixed_target(FFSContext_of_file(ffsfile), str_list);
    memset(&rec2, 0, sizeof(rec2));
    FFSread(ffsfile, &rec2);

    if (rec2.array != NULL) {
        printf("Expected first rec array dropped\n");
	return 1;
    }
    FFSread(ffsfile, &rec2);
    if (rec2.array == NULL) {
        printf("First rec array dropped when it shouldn't have been\n");
	return 1;
    }
    if(rec2.array[0] != 2.5) {
        printf("First rec array element 0 wrong.\n");
	return 1;
    }
    close_FFSfile(ffsfile);
    free_FFSfile(ffsfile);
    unlink(output_file);
    return 0;
}
