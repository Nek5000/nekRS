#ifndef ATL_H
#define ATL_H
#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif
#if defined(FUNCPROTO) || defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
#ifndef	NULL
#define NULL	((void *) 0)
#endif
#else
#ifndef	NULL
#define NULL	0
#endif
#endif

#ifndef _STDINT_H
#include <stdint.h>
#endif
#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <sys/types.h>
#endif

typedef int atom_t;

/* opaque type for atom server handle */
typedef struct _atom_server *atom_server;

extern
atom_t
atom_from_string (atom_server as, char *str);

extern
char *
string_from_atom (atom_server as, atom_t atom);

typedef enum {
    no_atom_cache, atom_cache, prefill_atom_cache
} atom_cache_type;

extern
atom_server
init_atom_server (atom_cache_type cache_style);

extern void
free_atom_server(atom_server as);

extern
char *
get_server_id (atom_server as);

typedef enum _attr_value_type { Attr_Undefined, Attr_Int4, Attr_Int8, 
				Attr_String, Attr_Opaque, Attr_Atom, 
				Attr_List, Attr_Float16, Attr_Float8, 
				Attr_Float4 } attr_value_type;

/* opaque type for attr_lists */
typedef struct _attr_list_struct *attr_list;

typedef int64_t attr_value;

typedef struct Attr_tmp_buffer *AttrBuffer;

#define ATL_CHAR_CONS(a,b,c,d) ((((unsigned char)a) << 24) + (((unsigned char)b) << 16) + (((unsigned char)c) << 8) + ((unsigned char)d))

typedef struct attr_opaque {
    int length;
    void *buffer;
} attr_opaque, *attr_opaque_p;

/* operations on attr_lists */
extern attr_list create_attr_list(void);

extern void free_attr_list (attr_list list);

extern attr_list add_ref_attr_list (attr_list list);

extern int attr_list_ref_count (attr_list list);

extern attr_list attr_join_lists (attr_list list1, attr_list list2);

extern attr_list attr_add_list (attr_list list1, attr_list list2);

extern void attr_merge_lists (attr_list list1, attr_list list2);

extern attr_list attr_copy_list (attr_list list);

extern int add_attr (attr_list attrs, atom_t attr_id, 
			  attr_value_type val_type,
			  attr_value value );

extern int add_float_attr (attr_list attrs, atom_t attr_id, 
				 double value );

extern int add_double_attr (attr_list attrs, atom_t attr_id, 
				 double value );

extern int add_int_attr (attr_list attrs, atom_t attr_id, 
			      int value );

extern int add_long_attr (attr_list attrs, atom_t attr_id, 
			      long value );

extern int add_string_attr (attr_list attrs, atom_t attr_id, 
			      char *value );

extern int add_opaque_attr (attr_list attrs, atom_t attr_id, 
				 int length, char *buffer );

extern int set_attr (attr_list attrs, atom_t attr_id, 
			  attr_value_type val_type,
			  attr_value value );

extern int set_float_attr (attr_list attrs, atom_t attr_id, 
				double value );

extern int set_double_attr (attr_list attrs, atom_t attr_id, 
				 double value );

extern int set_int_attr (attr_list attrs, atom_t attr_id, 
			      int value );

extern int set_long_attr (attr_list attrs, atom_t attr_id, 
			      ssize_t value );

extern int set_string_attr (attr_list attrs, atom_t attr_id, 
				 char *string );

extern int set_opaque_attr (attr_list attrs, atom_t attr_id, 
				 int length, char *buffer );

extern int replace_attr (attr_list attrs, atom_t attr_id, 
			      attr_value_type val_type, attr_value value );

extern int replace_float_attr (attr_list attrs, atom_t attr_id, 
				    double value );

extern int replace_double_attr (attr_list attrs, atom_t attr_id, 
				 double value );

extern int replace_int_attr (attr_list attrs, atom_t attr_id, 
			      int value );

extern int replace_long_attr (attr_list attrs, atom_t attr_id, 
			      long value );

extern int replace_string_attr (attr_list attrs, atom_t attr_id, 
				 char *string );

extern int replace_opaque_attr (attr_list attrs, atom_t attr_id, 
				     int length, char *buffer );

extern int query_attr ( attr_list attrs, atom_t attr_id, 
			    attr_value_type *val_type_p, attr_value *value_p);

extern int get_int_attr ( attr_list attrs, atom_t attr_id, int *value_p);
extern int get_long_attr ( attr_list attrs, atom_t attr_id, ssize_t *value_p);
extern int get_double_attr ( attr_list attrs, atom_t attr_id, double *value_p);
extern int get_float_attr ( attr_list attrs, atom_t attr_id, float *value_p);
extern int get_string_attr ( attr_list attrs, atom_t attr_id, char **value_p);

extern int get_opaque_attr ( attr_list attrs, atom_t attr_id, int *length_p, char **value_p);

extern void dump_attr_list ( attr_list attrs );
extern void fdump_attr_list ( void *file, attr_list attrs );

extern int get_attr_id(attr_list list, int item_no, atom_t *item);

extern int get_attr (attr_list list,int index, atom_t *name,
			  attr_value_type *val_type, attr_value *value);

typedef void (*atl_lock_func)(void *client_data);
extern void
atl_install_mutex_funcs(atl_lock_func lock, atl_lock_func unlock, void *client_data);

extern
atom_t
attr_atom_from_string (const char *str);

extern
char *
attr_string_from_atom (atom_t atom);

extern int attr_count (attr_list attrs);

extern
char *
attr_list_to_string (attr_list attrs);

extern
attr_list
attr_list_from_string (const char * str);

extern
AttrBuffer create_AttrBuffer(void);

void
free_AttrBuffer(AttrBuffer buf);

extern
int
attr_list_subset(attr_list l1, attr_list l2);

extern void *
encode_attr_for_xmit(attr_list l, AttrBuffer b, int *length);

extern attr_list 
decode_attr_from_xmit(void * buf);

extern int 
atl_base64_decode(unsigned char *input, unsigned char *output);

extern char *
atl_base64_encode(char *buffer, unsigned int len);

#if defined(__cplusplus) || defined(c_plusplus)
	   }
#endif

#endif
