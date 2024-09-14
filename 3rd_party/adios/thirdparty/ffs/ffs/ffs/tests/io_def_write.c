
#include "config.h"
#include <fcntl.h>
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include "ffs.h"
#include "unix_defs.h"

typedef struct _def_first_rec {
    double double_field;
} def_first_rec, *def_first_rec_ptr;

typedef struct _second_rec {
    short short_field;
    char *string;
    char char_field;
} second_rec, *second_rec_ptr;

typedef enum {Red_Stripe, Paulaner, Pilsner} enum_type;

typedef struct _third_rec {
    int integer_field;
    long long_field;
    char *string;
    char *string2;
    char char_field;
} third_rec, *third_rec_ptr;

#define ARRAY_SIZE 14

typedef struct _fourth_rec {
    long	ifield;
    int		int_array[ARRAY_SIZE];
    double double_array[2][2];
} fourth_rec, *fourth_rec_ptr;

typedef struct _embedded_rec {
    short	ifield;
    char	*string;
    double	dfield;
} embedded_rec, *embedded_rec_ptr;

typedef struct _fifth_rec {
    embedded_rec earray[4];
} fifth_rec, *fifth_rec_ptr;

typedef struct _later_rec {
#if SIZEOF_LONG_LONG > SIZEOF_LONG
    long long integer_field;
#else
    long integer_field;
#endif
    char *string;
    double double_field;
} later_rec, *later_rec_ptr;

typedef struct _later_rec2 {
#if SIZEOF_LONG_LONG > SIZEOF_LONG
    long long integer_field;
#else
    long integer_field;
#endif
    double double_field;
    char *string;
} later_rec2, *later_rec2_ptr;

typedef struct _nested_rec {
    int integer_field;
    second_rec  nested_rec;
    char *string;
} nested_rec, *nested_rec_ptr;

typedef struct _sixth_rec {
    char	*string;
    long	icount;
#if SIZEOF_LONG_LONG > SIZEOF_LONG
    long long	*var_int_array;
#else
    long	*var_int_array;
#endif
    second_rec	*var_string_array;
    double	dfield;
    double	*var_double_array;
} sixth_rec, *sixth_rec_ptr;

typedef struct _nested_variants {
    int vec_length;
    IOEncodeVector eventv;
} ninth_rec, *ninth_rec_ptr;

typedef struct _string_array {
    int array_len;
    char *base_string;
    char **array;
} string_array_rec, *string_array_rec_ptr;

typedef struct _ptr_rec {
    int integer_field;
    second_rec  *nested_rec;
    char *string;
} ptr_rec, *ptr_rec_ptr;

typedef struct _FormatListElement {
    char *format_name;
    int field_list_len;
    IOFieldList field_list;
} format_list_element;

typedef struct {
    int  len;
    char *channel;
} channel_ID_struct;

typedef struct _DeriveMsg
{
    char *chan_str;
    int cond;
    channel_ID_struct client_channel_id;
    char *client_contact_str;
    char *filter;
    int field_list_len;
    IOFieldList field_list;
    int format_list_len;
    format_list_element *format_list;

    char *init_data_block;
    int init_data_len;
} DeriveMsg, *DeriveMsgPtr;

typedef struct _multi_array_rec {
    long	ifield;
    double	double_array[2][2][2][2];
    int		(*int_array)[2];
    int		(*int_array2)[4];
    int		(*int_array3)[4][4];
} multi_array_rec, *multi_array_rec_ptr;


IOField def_field_list[] = {
    {"double field", "float",
       sizeof(double), IOOffset(def_first_rec_ptr, double_field)},
    { NULL, NULL, 0, 0}
};

IOField def_field_list2[] = {
    {"short field", "integer", 
       sizeof(short), IOOffset(second_rec_ptr, short_field)},
    {"string field", "string",
       sizeof(char *), IOOffset(second_rec_ptr, string)},
    {"char field", "char",
       sizeof(char), IOOffset(second_rec_ptr, char_field)},
    { NULL, NULL, 0, 0}
};

IOField def_field_list3[] = {
    {"integer field", "integer", 
       sizeof(int), IOOffset(third_rec_ptr, integer_field)},
    {"long field", "integer", 
       sizeof(long), IOOffset(third_rec_ptr, long_field)},
    {"string field", "string",
       sizeof(char *), IOOffset(third_rec_ptr, string)},
    {"string field2", "string",
       sizeof(char *), IOOffset(third_rec_ptr, string2)},
    {"char field", "char",
       sizeof(char), IOOffset(third_rec_ptr, char_field)},
    { NULL, NULL, 0, 0}
};

IOField def_field_list4[] = {
    {"ifield", "integer", 
       sizeof(long), IOOffset(fourth_rec_ptr, ifield)},
    {"int_array", IOArrayDecl(integer,ARRAY_SIZE), 
       sizeof(int), IOOffset(fourth_rec_ptr, int_array[0])},
    {"double field", "float[2][2]",
       sizeof(double), IOOffset(fourth_rec_ptr, double_array[0][0])},
    { NULL, NULL, 0, 0}
};

IOField def_embedded_field_list[] = {
    {"ifield", "integer",
       sizeof(short), IOOffset(embedded_rec_ptr, ifield)},
    {"string field", "string",
       sizeof(char *), IOOffset(embedded_rec_ptr, string)},
    {"dfield", "float",
       sizeof(double), IOOffset(embedded_rec_ptr, dfield)},
    {NULL, NULL, 0, 0}
};

IOField def_field_list5[] = {
    {"earray", "embedded[4]",
       sizeof(embedded_rec), IOOffset(fifth_rec_ptr, earray)},
    {NULL, NULL, 0, 0}
};

IOField def_later_field_list[] = {
    {"integer field", "integer", 
       sizeof(((later_rec_ptr)0)->integer_field), IOOffset(later_rec_ptr, integer_field)},
    {"string field", "string",
       sizeof(char *), IOOffset(later_rec_ptr, string)},
    {"double field", "float",
       sizeof(double), IOOffset(later_rec_ptr, double_field)},
    { NULL, NULL, 0, 0}
};

IOField def_later_field_list2[] = {
    {"integer field", "integer", 
       sizeof(((later_rec2_ptr)0)->integer_field), 
       IOOffset(later_rec2_ptr, integer_field)},
    {"string field", "string",
       sizeof(char *), IOOffset(later_rec2_ptr, string)},
    {"double field", "float",
       sizeof(double), IOOffset(later_rec2_ptr, double_field)},
    { NULL, NULL, 0, 0}
};

IOField def_nested_field_list[] = {
    {"integer field", "integer", 
       sizeof(((nested_rec_ptr)0)->integer_field), 
       IOOffset(nested_rec_ptr, integer_field)},
    {"nested record", "string format",
       sizeof(second_rec), IOOffset(nested_rec_ptr, nested_rec)},
    {"string field", "string",
       sizeof(char *), IOOffset(nested_rec_ptr, string)},
    { NULL, NULL, 0, 0}
};

IOField def_field_list6[] = {
    {"string field", "string",
       sizeof(char *), IOOffset(sixth_rec_ptr, string)},
    {"icount", "integer", 
       sizeof(long), IOOffset(sixth_rec_ptr, icount)},
    {"var_int_array", "integer[icount]",
       sizeof(((sixth_rec_ptr)0)->var_int_array[0]), IOOffset(sixth_rec_ptr, var_int_array)},
    {"var_string_array", "string format[icount]",
       sizeof(second_rec), IOOffset(sixth_rec_ptr, var_string_array)},
    {"double field", "float",
       sizeof(double), IOOffset(sixth_rec_ptr, dfield)},
    {"var_double_array", "float[icount]",
       sizeof(double), IOOffset(sixth_rec_ptr, var_double_array)},
    { NULL, NULL, 0, 0}
};

IOField def_event_vec_elem_fields[] =
{
    {"elem", "char[len]", sizeof(char), IOOffset(IOEncodeVector,iov_base)},
    {"len", "integer", sizeof(((IOEncodeVector)0)[0].iov_len), 
     IOOffset(IOEncodeVector, iov_len)},
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_field_list9[] =
{
    {"vec_length", "integer", sizeof(int),
     IOOffset(ninth_rec_ptr, vec_length)},
    {"eventv", "EventVecElem[vec_length]", sizeof(struct _io_encode_vec), 
     IOOffset(ninth_rec_ptr, eventv)},
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_string_array_field_list[] =
{
    {"array_len", "integer", sizeof(int),
     IOOffset(string_array_rec_ptr, array_len)},
    {"base_string", "string", sizeof(char*),
     IOOffset(string_array_rec_ptr, base_string)},
    {"array", "string[array_len]", sizeof(char*), 
     IOOffset(string_array_rec_ptr, array)},
    {(char *) 0, (char *) 0, 0, 0}

};

IOField def_channel_id_flds[] = {
    {"len", "integer", sizeof(int), IOOffset(channel_ID_struct *, len)},
    {"id", "char[len]", 1, IOOffset(channel_ID_struct *, channel)},
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_field_list_flds[] = {
    {"field_name", "string", sizeof(char *),
     IOOffset(IOFieldList, field_name)},
    {"field_type", "string", sizeof(char *),
     IOOffset(IOFieldList, field_type)},
    {"field_size", "integer", sizeof(int),
     IOOffset(IOFieldList, field_size)},
    {"field_offset", "integer", sizeof(int),
     IOOffset(IOFieldList, field_offset)},
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_format_list_field_list[] = {
    {"format_name", "string", sizeof(char *),
     IOOffset(format_list_element *, format_name)},
    {"field_list_len", "integer", sizeof(int),
     IOOffset(format_list_element *, field_list_len)},
    {"field_list", "IOfield_list[field_list_len]", sizeof(IOField),
     IOOffset(format_list_element *, field_list)}
    ,
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_derive_msg_field_list[] = {
    {"channel", "string", sizeof(char *),
     IOOffset(DeriveMsgPtr, chan_str)},
    {"condition", "integer", sizeof(int), IOOffset(DeriveMsgPtr, cond)},
    {"client_channel_id", "channel_ID", sizeof(channel_ID_struct),
     IOOffset(DeriveMsgPtr, client_channel_id)}
    ,
    {"client_contact_str", "string", sizeof(char *),
     IOOffset(DeriveMsgPtr, client_contact_str)},
    {"filter", "string", sizeof(char *), IOOffset(DeriveMsgPtr, filter)},
    {"field_list_len", "integer", sizeof(int),
     IOOffset(DeriveMsgPtr, field_list_len)},
    {"field_list", "IOfield_list[field_list_len]", sizeof(IOField),
     IOOffset(DeriveMsgPtr, field_list)}
    ,
    {"format_list_len", "integer", sizeof(int),
     IOOffset(DeriveMsgPtr, format_list_len)},
    {"format_list", "DEFormatList[format_list_len]",
     sizeof(format_list_element), IOOffset(DeriveMsgPtr, format_list)}
    ,
    {"init_data_block", "char[init_data_len]",
     1, IOOffset(DeriveMsgPtr, init_data_block)}
    ,
    {"init_data_len", "integer",
     sizeof(int), IOOffset(DeriveMsgPtr, init_data_len)},
    {(char *) 0, (char *) 0, 0, 0}
};

IOField def_multi_array_flds[] = {
    {"ifield", "integer", sizeof(long), IOOffset(multi_array_rec_ptr, ifield)},
    {"double_array", "float[2][2][2][2]", sizeof(double),
     IOOffset(multi_array_rec_ptr, double_array)},
    {"int_array", "integer[2][ifield]", sizeof(int),
     IOOffset(multi_array_rec_ptr, int_array)},
    {"int_array2", "integer[ifield][2]", sizeof(int),
     IOOffset(multi_array_rec_ptr, int_array2)},
    {"int_array3", "integer[ifield][ifield][ifield]", sizeof(int),
    IOOffset(multi_array_rec_ptr, int_array3)},
    {(char *) 0, (char *) 0, 0, 0}
};

int
main()
{
    IOFile iofile = open_IOfile("def_test_output", "w");
    IOFormat first_rec_ioformat, second_rec_ioformat, third_rec_ioformat;
    IOFormat fourth_rec_ioformat, later_ioformat, nested_ioformat;
    IOFormat embedded_rec_ioformat, fifth_rec_ioformat, sixth_rec_ioformat;
    IOFormat ninth_rec_ioformat;
    IOFormat string_array_ioformat;
    struct _pbiovec vector[5];
    def_first_rec rec1;
    second_rec rec2;
    third_rec rec3;
    fourth_rec rec4;
    later_rec rec5;
    later_rec2 rec6;
    nested_rec rec7;
    fifth_rec emb_array;
    sixth_rec var_array;
    ninth_rec var_var;
    string_array_rec str_array;
    int i, j;
    IOOptInfo opt_info[2];

    def_first_rec array1[10];

    opt_info[0].info_type = 0;
    opt_info[0].info_len = 0;
    opt_info[0].info_block = NULL;

    first_rec_ioformat = register_opt_format("first format",
					     def_field_list, opt_info,
					     (IOContext) iofile);
    second_rec_ioformat = register_opt_format("string format",
					      def_field_list2, opt_info, 
					      (IOContext) iofile);
    third_rec_ioformat = register_opt_format("two string format",
					     def_field_list3, opt_info, 
					     (IOContext) iofile);
    fourth_rec_ioformat = register_opt_format("internal array format",
					      def_field_list4, opt_info,
					      (IOContext) iofile);
    embedded_rec_ioformat = register_opt_format("embedded",
						def_embedded_field_list, opt_info,
						(IOContext) iofile);
    (void)embedded_rec_ioformat;
    fifth_rec_ioformat = register_opt_format("structured array format",
					     def_field_list5, opt_info, 
					     (IOContext) iofile);
    sixth_rec_ioformat = register_opt_format("variant array format",
					     def_field_list6, opt_info,
					     (IOContext) iofile);
    write_comment_IOfile(iofile, "this is a comment in the file");
    memset((char *) &rec1, 0, sizeof(rec1));
    rec1.double_field = 2.717;
    if (!write_IOfile(iofile, first_rec_ioformat, &rec1))
	IOperror(iofile, "write failed\n");
    memset((char *) &emb_array, 0, sizeof(emb_array));
    memset((char *) &var_array, 0, sizeof(var_array));
    emb_array.earray[0].dfield = 4.0;
    emb_array.earray[0].ifield = 4;
    emb_array.earray[0].string = (char *) malloc(10);
    sprintf(emb_array.earray[0].string, "string%d", emb_array.earray[0].ifield * 5);
    emb_array.earray[1].dfield = 3.0;
    emb_array.earray[1].ifield = 3;
    emb_array.earray[1].string = (char *) malloc(10);
    sprintf(emb_array.earray[1].string, "string%d", emb_array.earray[1].ifield * 5);
    emb_array.earray[2].dfield = 2.0;
    emb_array.earray[2].ifield = 2;
    emb_array.earray[2].string = (char *) malloc(10);
    sprintf(emb_array.earray[2].string, "string%d", emb_array.earray[2].ifield * 5);
    emb_array.earray[3].dfield = 1.0;
    emb_array.earray[3].ifield = 1;
    emb_array.earray[3].string = (char *) malloc(10);
    sprintf(emb_array.earray[3].string, "string%d", emb_array.earray[3].ifield * 5);
    if (!write_IOfile(iofile, fifth_rec_ioformat, &emb_array))
	IOperror(iofile, "write failed\n");
    IOfree_var_rec_elements(iofile, fifth_rec_ioformat, &emb_array);

    memset((char *) &rec2, 0, sizeof(rec2));
    rec2.short_field = 27;
    rec2.string = "testing";
    rec2.char_field = 'A';
    if (!write_IOfile(iofile, second_rec_ioformat, &rec2))
	IOperror(iofile, "write failed\n");
    rec2.short_field = 27;
    rec2.string = NULL;
    rec2.char_field = 'A';
    if (!write_IOfile(iofile, second_rec_ioformat, &rec2))
	IOperror(iofile, "write failed\n");
    rec1.double_field *= 3.0;
    write_IOfile(iofile, first_rec_ioformat, &rec1);
    rec2.short_field = 27;
    rec2.string = NULL;
    rec2.char_field = 'A';
    if (!write_IOfile(iofile, second_rec_ioformat, &rec2))
	IOperror(iofile, "write failed\n");
    rec1.double_field *= 2.717;
    write_comment_IOfile(iofile, "this is another comment in the file");
    write_IOfile(iofile, first_rec_ioformat, &rec1);
    memset((char *) &rec3, 0, sizeof(rec3));
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    if (!write_IOfile(iofile, third_rec_ioformat, &rec3))
	IOperror(iofile, "write failed\n");
    later_ioformat = register_opt_format("later format",
					 def_later_field_list, opt_info, 
					 (IOContext) iofile);
    nested_ioformat = register_opt_format("nested format",
					  def_nested_field_list, opt_info, 
					  (IOContext) iofile);
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    if (!write_IOfile(iofile, third_rec_ioformat, &rec3))
	IOperror(iofile, "write failed\n");
    memset((char *) &rec7, 0, sizeof(rec7));
    rec7.integer_field = 47;
    rec7.nested_rec.short_field = 27;
    rec7.nested_rec.string = "Another string";
    rec7.nested_rec.char_field = 'A';
    rec7.string = "Yet another string";
    if (!write_IOfile(iofile, nested_ioformat, &rec7))
	IOperror(iofile, "write failed\n");
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    if (!write_IOfile(iofile, third_rec_ioformat, &rec3))
	IOperror(iofile, "write failed\n");
    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = NULL;
    rec3.string2 = NULL;
    rec3.char_field = 'A';
    if (!write_IOfile(iofile, third_rec_ioformat, &rec3))
	IOperror(iofile, "write failed\n");
    memset((char *) &rec5, 0, sizeof(rec5));
    rec5.integer_field = 9872346;
    rec5.string = "ABCD";
    rec5.double_field = 3.14159265358797323;
    vector[0].data = &rec5;
    vector[0].format = later_ioformat;

    rec3.integer_field = 14;
    rec3.long_field = 987234;
    rec3.string = "testing";
    rec3.string2 = "jambalaya";
    rec3.char_field = 'A';
    vector[1].data = &rec3;
    vector[1].format = third_rec_ioformat;

    rec2.short_field = 27;
    rec2.string = "the end";
    rec2.char_field = 'A';
    vector[2].data = &rec2;
    vector[2].format = second_rec_ioformat;
    if (!writev_IOfile(iofile, &(vector[0]), 3))
	IOperror(iofile, "writev failed\n");
    for (i = 0; i < 10; i++) {
	memset((char *) &array1[i], 0, sizeof(array1[i]));
	array1[i].double_field = 2.717 * (i * i);
    }
    later_ioformat = register_opt_format("later format",
					 def_later_field_list2, opt_info,
					 (IOContext) iofile);
    if (!write_array_IOfile(iofile, first_rec_ioformat, &array1[0],
			    10, sizeof(array1[0])))
	IOperror(iofile, "write failed\n");

    memset((char *) &rec4, 0, sizeof(rec4));
    for (i = 0; i < ARRAY_SIZE; i++) {
	rec4.int_array[i] = 297 + i;
    }
    rec4.double_array[0][0] = 1.0;
    rec4.double_array[0][1] = 2.0;
    rec4.double_array[1][0] = 3.0;
    rec4.double_array[1][1] = 4.0;
    rec4.ifield = -rec4.int_array[ARRAY_SIZE - 1];
    vector[0].data = &rec4;
    vector[0].format = fourth_rec_ioformat;

    memset((char *) &rec6, 0, sizeof(rec6));
    rec6.integer_field = 23462346;
    rec6.string = "Efghij";
    rec6.double_field = 3.14159265358797323 * 2.0;
    vector[1].data = &rec6;
    vector[1].format = later_ioformat;
    if (!writev_IOfile(iofile, &(vector[0]), 2))
	IOperror(iofile, "writev failed\n");
    rec6.integer_field = 2346987;
    rec6.string = "Klmn";
    rec6.double_field = 3.14159265358797323 * 3.0;
    
    if (!write_IOfile(iofile, later_ioformat, &rec6))
	IOperror(iofile, "write failed\n");

    for (i = 1; i < 20; i += 5) {
        memset((char *) &var_array, 0, sizeof(var_array));
	var_array.string = malloc(10);
	sprintf(var_array.string, "variant%d", i);
	var_array.icount = 2 * i;
	var_array.var_int_array = malloc(sizeof(((sixth_rec_ptr) 0)->var_int_array[0]) * var_array.icount);
	var_array.var_double_array = malloc(sizeof(double) * var_array.icount);
	var_array.var_string_array = malloc(sizeof(second_rec) * var_array.icount);
        memset((char *) var_array.var_string_array, 0, 
	       sizeof(second_rec) * var_array.icount);
	for (j = 0; j < var_array.icount; j++) {
	    var_array.var_int_array[j] = 297 + j;
	    var_array.var_double_array[j] = 2.717 * j;
	    var_array.var_string_array[j].short_field = j;
	    var_array.var_string_array[j].string = malloc(15);
	    sprintf(var_array.var_string_array[j].string,
		    "substring%d", j);
	    var_array.var_string_array[j].char_field = 'a' + 2 * j;
	}
	if (!write_IOfile(iofile, sixth_rec_ioformat, &var_array))
	    IOperror(iofile, "write failed");
	IOfree_var_rec_elements(iofile, sixth_rec_ioformat, &var_array);
    }

    register_opt_format("EventVecElem", def_event_vec_elem_fields, opt_info,
			(IOContext) iofile);

    ninth_rec_ioformat = register_opt_format("EventV",
					     def_field_list9, opt_info,
					     (IOContext) iofile);
    string_array_ioformat = register_IOrecord_format("string_array",
						     def_string_array_field_list,
						iofile);
    for (i = 1; i < 10; i += 2) {
        memset((char *) &var_var, 0, sizeof(var_var));
	memset((char *) &str_array, 0, sizeof(str_array));
	var_var.vec_length = i;
	var_var.eventv = malloc(sizeof(((ninth_rec_ptr) 0)->eventv[0]) * var_var.vec_length);

	str_array.array_len = i;
	str_array.array = malloc(sizeof(char*) * i);

	for (j = 0; j < var_var.vec_length; j++) {
	    int k;
	    var_var.eventv[j].iov_len = j + i;
	    var_var.eventv[j].iov_base = malloc(j + i);
	    str_array.array[j] = malloc(i + j + 2);
	    for (k=0; k<j+i; k++) {
		((char*)var_var.eventv[j].iov_base)[k] = 'A' + k + i/5;
		str_array.array[j][k] = 'a' + k + i;
	    }
	    str_array.array[j][k] = 0;
	}
	if ((i % 4) == 1) {
	    str_array.base_string = NULL;
	} else {
	    str_array.base_string = strdup("Whoa there!");
	}
	if (!write_IOfile(iofile, ninth_rec_ioformat, &var_var))
	    IOperror(iofile, "write failed");
	if (!write_IOfile(iofile, string_array_ioformat, &str_array))
	    IOperror(iofile, "write failed");
	IOfree_var_rec_elements(iofile, string_array_ioformat, &str_array);
	IOfree_var_rec_elements(iofile, ninth_rec_ioformat, &var_var);
    }

    close_IOfile(iofile);
    free_IOfile(iofile);
    return 0;
}


    
