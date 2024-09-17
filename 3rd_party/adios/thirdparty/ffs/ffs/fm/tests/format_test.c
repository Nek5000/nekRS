
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include "fm.h"
#include "fm_internal.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define sleep(x) Sleep(1000*x)
#endif

#include "test_funcs.h"

int
main(int argc, char **argv)
{

    int continuous_test = 0;
    int restart_test = 0;
    int get_test = 0;
    int XML_test = 0;
    FMStructDescRec str_list[5];

    if (argc > 1) {
	if (strcmp(argv[1], "-c") == 0) {
	    continuous_test++;
	    printf("Doing continuous test\n");
	}
	if (strcmp(argv[1], "-r") == 0) {
	    restart_test++;
	    printf("Doing restart test\n");
	}
	if (strcmp(argv[1], "-g") == 0) {
	    get_test++;
	}
	if (strcmp(argv[1], "-x") == 0) {
	    XML_test++;
	    printf("Doing XML test\n");
	}
    }
    
    if (get_test) {
	FMFormat first_rec_ioformat;
	FMContext context = create_FMcontext();
	FMcontext_allow_self_formats(context);
	printf("Format server identifier is %x\n", FMcontext_get_format_server_identifier(context));
	unsigned char id[] = {02, 00, 00, 37, 103, 189, 231, 165, 33, 254, 42, 32};
	printf("Doing get test\n");
	first_rec_ioformat = FMformat_from_ID(context, (char *) &id[0]);
	printf("format is %p\n", first_rec_ioformat);
	return 0;
    }
    if (XML_test) {
	FMFormat first_rec_ioformat;
	FMContext context = create_FMcontext();
	FMcontext_allow_self_formats(context);
	FMOptInfo opt_info[2];
	printf("Format server identifier is %x\n", FMcontext_get_format_server_identifier(context));
	opt_info[0].info_type = 0x584D4C20;
	opt_info[0].info_len = 44;
	opt_info[0].info_block = "this is a bunch of text meant to be XML text";
	opt_info[1].info_type = 0;
	str_list[0].format_name = "first format";
	str_list[0].field_list = field_list;
	str_list[0].struct_size = sizeof(first_rec);
	str_list[0].opt_info = opt_info;
	str_list[1].format_name = NULL;
	first_rec_ioformat = register_data_format(context, str_list);
	dump_FMFormat(first_rec_ioformat);
    } else if (!continuous_test && !restart_test) {
	FMContext context = create_FMcontext();
	FMContext local_context = create_local_FMcontext();
	FMFormat first_rec_ioformat;
	FMcontext_allow_self_formats(context);
	str_list[0].format_name = "first format";
	str_list[0].field_list = field_list;
	str_list[0].struct_size = sizeof(first_rec);
	str_list[0].opt_info = NULL;
	str_list[1].format_name = NULL;
	first_rec_ioformat = register_data_format(context, str_list);
	printf("format is %p\n", first_rec_ioformat);
	if (first_rec_ioformat != NULL) {
	    int id_len, i;
	    char *id =get_server_ID_FMformat(first_rec_ioformat, &id_len);
	    printf("Remote context ID is : ");
	    print_server_ID((unsigned char*)id);
	    printf("Hex is : ");
	    for (i=0; i < id_len; i++) {
		printf("0x%02d, ", ((unsigned char*)id)[i]);
	    }
	    printf("\n");
	}
	first_rec_ioformat = register_data_format(local_context, str_list);
	if (first_rec_ioformat != NULL) {
	    int id_len;
	    char *id =get_server_ID_FMformat(first_rec_ioformat, &id_len);
	    printf("Local context ID is :  ");
	    print_server_ID((unsigned char *)id);
	}
    } else if (restart_test) {
	FMContext static_context = create_FMcontext();
	while (1) {
	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    sleep(120);
	}
    } else {
	FMContext static_context = create_FMcontext();
	while (1) {
	    FMContext context = create_FMcontext();
	    FMFormat first_rec_ioformat;

	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "first format";
	    str_list[0].field_list = field_list;
	    str_list[0].struct_size = sizeof(first_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("first format is %p\n", first_rec_ioformat);
	    sleep(10);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "string format";
	    str_list[0].field_list = field_list2;
	    str_list[0].struct_size = sizeof(second_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("second format is %p\n", first_rec_ioformat);
	    sleep(20);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "two string format";
	    str_list[0].field_list = field_list3;
	    str_list[0].struct_size = sizeof(third_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("third format is %p\n", first_rec_ioformat);
	    sleep(30);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "internal array format";
	    str_list[0].field_list = field_list4;
	    str_list[0].struct_size = sizeof(fourth_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("fourth format is %p\n", first_rec_ioformat);
	    sleep(40);

	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "embedded";
	    str_list[0].field_list = embedded_field_list;
	    str_list[0].struct_size = sizeof(embedded_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("embedded format is %p\n", first_rec_ioformat);
	    sleep(50);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "structured array format";
	    str_list[0].field_list = field_list5;
	    str_list[0].struct_size = sizeof(fifth_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = "embedded";
	    str_list[1].field_list = embedded_field_list;
	    str_list[1].struct_size = sizeof(embedded_rec);
	    str_list[1].opt_info = NULL;
	    str_list[2].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("fifth format is %p\n", first_rec_ioformat);
	    sleep(60);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "later format";
	    str_list[0].field_list = later_field_list;
	    str_list[0].struct_size = sizeof(later_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("later format is %p\n", first_rec_ioformat);
	    sleep(70);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "later format";
	    str_list[0].field_list = later_field_list2;
	    str_list[0].struct_size = sizeof(later_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("later2 format is %p\n", first_rec_ioformat);
	    sleep(80);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "nested format";
	    str_list[0].field_list = nested_field_list;
	    str_list[0].struct_size = sizeof(nested_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = "string format";
	    str_list[1].field_list = field_list2;
	    str_list[1].struct_size = sizeof(second_rec);
	    str_list[1].opt_info = NULL;
	    str_list[2].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("nested format is %p\n", first_rec_ioformat);
	    sleep(90);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "variant array format";
	    str_list[0].field_list = field_list6;
	    str_list[0].struct_size = sizeof(sixth_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = "string format";
	    str_list[1].field_list = field_list2;
	    str_list[1].struct_size = sizeof(second_rec);
	    str_list[1].opt_info = NULL;
	    str_list[2].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("sixth format is %p\n", first_rec_ioformat);
	    sleep(100);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("eventvec format is %p\n", first_rec_ioformat);
	    sleep(110);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "EventV";
	    str_list[0].field_list = field_list9;
	    str_list[0].struct_size = sizeof(ninth_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = "EventVecElem";
	    str_list[1].field_list = event_vec_elem_fields;
	    str_list[1].struct_size = sizeof(struct _io_encode_vec);
	    str_list[1].opt_info = NULL;
	    str_list[2].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("ninth format is %p\n", first_rec_ioformat);
	    sleep(120);


	    if (format_server_restarted(static_context)) {
		printf("Format server was restarted\n");
	    }
	    str_list[0].format_name = "string_array";
	    str_list[0].field_list = string_array_field_list;
	    str_list[0].struct_size = sizeof(string_array_rec);
	    str_list[0].opt_info = NULL;
	    str_list[1].format_name = NULL;
	    first_rec_ioformat = register_data_format(context, str_list);
	    printf("string array format is %p\n", first_rec_ioformat);
	    sleep(130);

	    free_FMcontext(context);
	}
	/* NOTREACHED */
	/* free_IOcontext(static_context);*/
    }
    return 0;
}
