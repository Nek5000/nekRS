// This is a copied program from ffs_write.c

#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ffs.h"

typedef struct _second_rec {
    int integer_field;
    short short_field;
    long long_field;
    char *string;
    double double_field;
    char char_field;
} second_rec, *second_rec_ptr;

FMField field_list2[] = {
    {"integer field(-14)", "integer", 
       sizeof(int), FMOffset(second_rec_ptr, integer_field)},
    {"short field", "integer", 
       sizeof(short), FMOffset(second_rec_ptr, short_field)},
    {"long field(-9876129)", "integer", 
       sizeof(long), FMOffset(second_rec_ptr, long_field)},
    {"string field", "string",
       sizeof(char *), FMOffset(second_rec_ptr, string)},
    {"double field(2.717)", "float",
       sizeof(double), FMOffset(second_rec_ptr, double_field)},
    {"char field", "char",
       sizeof(char), FMOffset(second_rec_ptr, char_field)},
    { NULL, NULL, 0, 0}
};

char *string_xml = "\
<StringRecord integer_attribute=<FFS:data field_id=0>>\n\
    <Short value=<FFS:data field_id=1>>\n\
    <Long  value=<FFS:data field_id=2>>\n\
    <String>\n\
       <FFS:data field_id=3>\n\
    </String>\n\
    <Double value=<FFS:data field_id=4>> \n\
    <Character value=<FFS:data field_id=5>>\n\
</StringRecord>\n";

int
main(int argc, char **argv)
{
    FMContext src_context;
    FFSFile ffsfile;
    FMFormat second_rec_ioformat;

    second_rec rec2;

    FMStructDescRec str_list[5];
    FMOptInfo opt_info[2], opt_info2[2];
    int verbose = 0, indexed = 0, i;
    char *output_file = "test_output";

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
	    /* argument */
	    if (argv[i][1] == 'v') {
		verbose++;
	    } else if (argv[i][1] == 'i') {
		indexed++;
	    } else {
		printf("Unknown argument \"%s\"\n", argv[i]);
	    }
	} else {
	    output_file = argv[i];
	}
    }

    if (!indexed) {
	ffsfile = open_FFSfile(output_file, "a");
    } else {
	ffsfile = open_FFSfile(output_file, "ai");
    }

    src_context = create_local_FMcontext();

    str_list[0].format_name = "string format";
    str_list[0].field_list = field_list2;
    str_list[0].struct_size = sizeof(second_rec);
    str_list[0].opt_info = &opt_info[0];
    str_list[1].format_name = NULL;
    opt_info[0].info_len = strlen(string_xml) +1;
    opt_info[0].info_block = string_xml;
    second_rec_ioformat = register_data_format(src_context, str_list);

    memset((char *) &rec2, 0, sizeof(rec2));
    rec2.integer_field = 14;
    rec2.short_field = 27;
    rec2.long_field = 987234;
    rec2.string = "testing";
    rec2.double_field = 2.717;
    rec2.char_field = 'A';
    if (!write_FFSfile(ffsfile, second_rec_ioformat, &rec2))
	printf("write failed\n");

    write_comment_FFSfile(ffsfile, "this is another comment in the file");

    close_FFSfile(ffsfile);
    free_FFSfile(ffsfile);
    return 0;
}
