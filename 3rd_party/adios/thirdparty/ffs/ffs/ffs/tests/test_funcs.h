#ifdef _MSC_VER
#define strdup(s) _strdup(s)
#define open(x,y,z) _open(x,y,z)
#define read(x,y,z) _read(x,y,z)
#define write(x,y,z) _write(x,y,z)
#define close(x) _close(x)
#define unlink(x) _unlink(x)
#include <io.h>
#endif
typedef struct _newer_rec {
    int ganzzahl;
    double gleitkommazahl;
    char zeichen;
} newer_rec, *newer_rec_ptr;

typedef struct _first_rec {
    int integer_field;
    double double_field;
    char char_field;
} first_rec, *first_rec_ptr;

typedef struct _second_rec {
    int integer_field;
    short short_field;
    long long_field;
    char *string;
    double double_field;
    char char_field;
} second_rec, *second_rec_ptr;

typedef enum {Red_Stripe, Paulaner, Pilsner} enum_type;

typedef struct _third_rec {
    int integer_field;
    long long_field;
    unsigned int uint_field;
    unsigned long ulong_field;
    char *string;
    double double_field;
    char *string2;
    char char_field;
    enum_type enum_field;
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

typedef struct _io_encode_vec {
     void *iov_base;
     long  iov_len;
} *IOEncodeVector;

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
    FMFieldList field_list;
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
    FMFieldList field_list;
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

typedef struct _multi_array_rec2 {
    long	ifield;
    double	double_array[2][2][2][2];
    int		(*int_array)[2][4];
    int		(*int_array2)[4][2];
    int		(*int_array3)[4][4][4];
} multi_array_rec2, *multi_array_rec2_ptr;


typedef struct
{
  int id;
  int data_size;
    int old_data_size;
  int compressed_size;
  int * compressed_data;

  int corner1x;
  int corner1y;
  int corner1z;
  int corner2x;
  int corner2y;
  int corner2z;

  int atom_type;
} compressed_mesh_param, *compressed_mesh_param_ptr;

/* form pipelined server to pipelined client */
typedef struct
{
  int corner1x;
  int corner1y;
  int corner1z;
  int corner2x;
  int corner2y;
  int corner2z;

  int compression_type;

  int codebook_size;
  int* codebook_data;

  int timestamp;
  int nonce;
  int num_meshes;
  compressed_mesh_param* mesh_data;
} triangle_param, *triangle_param_ptr;


typedef struct _XMLFormatListElement {
    char *format_name;
    int field_list_len;
    FMFieldList field_list;
    char *xml_markup;
} msg_format_list_element;

typedef struct {
    int action;
    char *in_format_name;
    char *func_str;
    int format_count;
    msg_format_list_element *out_formats;
} add_rec, *add_rec_ptr;

typedef struct node {
    int node_num;
    struct node *link1;
    struct node *link2;
} *node_ptr;

typedef struct visit_table {
    int node_count;
    void *nodes[100];
} *visit_table;

typedef struct pointer_to_static_array {
    int (*int_array)[4][4];
} *pointer_to_static_array_ptr;

extern int calc_signature(node_ptr n, visit_table v);

typedef struct _cp_init_info {
    char *contact_info;
    int target_stone;
    void *reader_ID;
} *cp_init_info;

struct _reader_register_msg {
    void *writer_file;
    int writer_response_condition;
    int reader_cohort_size;
    cp_init_info *CP_reader_info;
};

extern FMField cp_reader_init_list[];
extern FMField cp_reader_register_list[];
extern FMStructDescRec reader_register_format_list[];

extern FMField field_list[];
extern FMField newer_field_list[];
extern FMField field_list2[];
extern FMField field_list3[];
extern FMField field_list4[];
extern FMField embedded_field_list[];
extern FMField field_list5[];
extern FMField later_field_list[];
extern FMField later_field_list2[];
extern FMField nested_field_list[];
extern FMField field_list6[];
extern FMField event_vec_elem_fields[];
extern FMField field_list9[];
extern FMField string_array_field_list[];
extern FMField derive_msg_field_list[];
extern FMField format_list_field_list[];
extern FMField field_list_flds[];
extern FMField channel_id_flds[];
extern FMField multi_array_flds[];
extern FMField compressed_mesh[];
extern FMField triangle_field[];
extern FMField xml_format_list_flds[];
extern FMField add_field_list[];
extern FMField node_field_list[];
extern FMField pointer_to_static_field_list[];

extern FMStructDescRec first_format_list[];
extern FMStructDescRec string_format_list[];
extern FMStructDescRec structured_format_list[];
extern FMStructDescRec two_string_format_list[];
extern FMStructDescRec fourth_format_list[];
extern FMStructDescRec later_format_list[];
extern FMStructDescRec nested_format_list[];
extern FMStructDescRec embedded_format_list[];
extern FMStructDescRec variant_format_list[];
extern FMStructDescRec ninth_format_list[];
extern FMStructDescRec string_array_format_list[];
extern FMStructDescRec derive_format_list[];
extern FMStructDescRec multi_array_format_list[];
extern FMStructDescRec triangle_format_list[];
extern FMStructDescRec add_action_format_list[];
extern FMStructDescRec node_format_list [];
extern FMStructDescRec pointer_to_static_format_list [];

extern void init_written_data();
extern void free_written_data();

extern char *comment_array[];

extern first_rec rec1_array[13];
extern int rec1_attr_val[13];

extern second_rec rec2_array[4];
extern int rec2_attr_val[4];

extern third_rec rec3_array[5];
extern int rec3_attr_val[5];

extern fourth_rec rec4;

extern fifth_rec rec5;
extern int rec5_attr_val;

extern sixth_rec rec6_array[4];

extern nested_rec rec7_array[1];

extern later_rec rec8_array[3];

extern ninth_rec rec9_array[6];

extern string_array_rec string_array_array[6];

extern DeriveMsg derive;

extern multi_array_rec multi_array;

extern multi_array_rec2 multi_array2;

extern multi_array_rec fortran_array;

extern triangle_param triangle;

extern add_rec add_action_record;

extern struct pointer_to_static_array psa;

extern struct _reader_register_msg reader_register;

extern int first_rec_eq(first_rec *r1, first_rec *r2);
extern int second_rec_eq(second_rec *r1, second_rec *r2);
extern int third_rec_eq(third_rec *r1, third_rec *r2);
extern int fourth_rec_eq(fourth_rec *r1, fourth_rec *r2);
extern int emb_rec_eq(embedded_rec *r1, embedded_rec *r2);
extern int fifth_rec_eq(fifth_rec *r1, fifth_rec *r2);
extern int sixth_rec_eq(sixth_rec *r1, sixth_rec *r2);
extern int nested_rec_eq(nested_rec *r1, nested_rec *r2);
extern int later_rec_eq(later_rec *r1, later_rec *r2);
extern int ninth_rec_eq(ninth_rec *r1, ninth_rec *r2);
extern int string_array_eq(string_array_rec *r1, string_array_rec *r2);
extern int derive_eq(DeriveMsgPtr d1, DeriveMsgPtr d2);
extern int multi_array_eq(multi_array_rec *r1, multi_array_rec *r2);
extern int triangle_param_eq(triangle_param *r1, triangle_param *r2);
extern int add_rec_eq(add_rec_ptr r1, add_rec_ptr r2);
extern int
pointer_to_static_rec_eq(pointer_to_static_array_ptr r1, pointer_to_static_array_ptr r2);
extern int
reader_register_rec_eq(struct _reader_register_msg *r1, struct _reader_register_msg *r2);

extern const char *first_xml;
extern char *string_xml;
extern char *third_xml;
extern char *fourth_xml;
extern char *embedded_xml;
extern char *struct_array_xml;
extern char *var_array_xml;
extern char *later_xml;
extern char *later2_xml;
extern char *nested_xml;
extern char *event_xml;
extern char *event_vec_xml;
extern char *multi_array_xml;
