// This is a copied program from ffs_write.c

#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "ffs.h"

extern void
install_subsample_code(FMFormat f, char *field, char*code_str);

typedef struct _second_rec {
    int len, len2;
    int *type_array;
    double *value_array;
    char char_field;
} second_rec, *second_rec_ptr;

FMField field_list2[] = {
    {"len", "integer", 
       sizeof(int), FMOffset(second_rec_ptr, len)},
    {"len2", "integer", 
       sizeof(int), FMOffset(second_rec_ptr, len2)},
    {"type_array", "integer[len]",
       sizeof(int), FMOffset(second_rec_ptr, type_array)},
    {"value_array", "float[len2]",
       sizeof(double), FMOffset(second_rec_ptr, value_array)},
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
    rec2.len2 = 10;
    rec2.type_array = malloc(sizeof(rec2.type_array[0]) * rec2.len);
    rec2.value_array = malloc(sizeof(rec2.value_array[0]) * rec2.len2);
    for (i=0; i < rec2.len; i++) {
	rec2.type_array[i] = 1;
	if ((i % 2) == 0) rec2.type_array[i] = 2;
	if ((i % 3) == 0) rec2.type_array[i] = 3;
	rec2.value_array[i] = 2*i + 10;
    }
    rec2.char_field = 'A';
    install_subsample_code(second_rec_ioformat, "type_array", "{int i; for (i=0;i<element_count;i++) { if (input.type_array[i] == 2) { FFSMarshalArrayElement(i);}}}");
    install_subsample_code(second_rec_ioformat, "value_array", "{int i; for (i=0;i<element_count;i++) { if (input.type_array[i] == 2) { FFSMarshalArrayElement(i);}}}");
    if (!write_FFSfile(ffsfile, second_rec_ioformat, &rec2))
	printf("write failed\n");

    for (i=0; i < rec2.len; i++) {
	rec2.type_array[i] = 1;
	if ((i % 3) == 0) rec2.type_array[i] = 3;
	if ((i % 2) == 0) rec2.type_array[i] = 2;
	rec2.value_array[i] = 2*i + 10;
    }

    if (!write_FFSfile(ffsfile, second_rec_ioformat, &rec2))
	printf("write failed\n");

    close_FFSfile(ffsfile);
    free_FFSfile(ffsfile);

    ffsfile = open_FFSfile(output_file, "r");
    FFSset_fixed_target(FFSContext_of_file(ffsfile), str_list);
    memset(&rec2, 0, sizeof(rec2));
    FFSread(ffsfile, &rec2);

/* change success conditions */
    if (rec2.len != 3) {
        printf("Expected len array truncated to 3 elements, found %d\n", rec2.len);
	return 1;
    }
    if (rec2.len2 != 3) {
        printf("Expected len2 array truncated to 3 elements, found %d\n", rec2.len);
	return 1;
    }
    for(i=0; i< rec2.len; i++) {
	if (rec2.type_array[i] != 2) 
	    printf("Expected type_array[%d] == 2, found %d\n", i, rec2.type_array[i]);
    }
    if (rec2.value_array[0] != 14) {
	printf("Expected element 0 == 14, got %g\n", rec2.value_array[0]);
	return 1;
    }
    if (rec2.value_array[1] != 18) {
	printf("Expected element 1 == 18, got %g\n", rec2.value_array[1]);
	return 1;
    }
    if (rec2.value_array[2] != 26) {
	printf("Expected element 2 == 26, got %g\n", rec2.value_array[2]);
	return 1;
    }

    FFSread(ffsfile, &rec2);
    if (rec2.len != 5) {
        printf("Expected len array truncated to 3 elements, found %d\n", rec2.len);
	return 1;
    }
    if (rec2.len2 != 5) {
        printf("Expected len2 array truncated to 3 elements, found %d\n", rec2.len);
	return 1;
    }
    for(i=0; i< rec2.len; i++) {
	if (rec2.type_array[i] != 2) 
	    printf("Expected type_array[%d] == 2, found %d\n", i, rec2.type_array[i]);
    }
    if (rec2.value_array[0] != 10) {
	printf("2Expected element 0 == 10, got %g\n", rec2.value_array[0]);
	return 1;
    }
    if (rec2.value_array[1] != 14) {
	printf("2Expected element 1 == 14, got %g\n", rec2.value_array[1]);
	return 1;
    }
    if (rec2.value_array[2] != 18) {
	printf("2Expected element 2 == 18, got %g\n", rec2.value_array[2]);
	return 1;
    }
    if (rec2.value_array[3] != 22) {
	printf("2Expected element 3 == 22, got %g\n", rec2.value_array[1]);
	return 1;
    }
    if (rec2.value_array[4] != 26) {
	printf("2Expected element 4 == 26, got %g\n", rec2.value_array[2]);
	return 1;
    }

    unlink(output_file);
    return 0;
}
