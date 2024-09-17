typedef struct {
    int line;
    int character;
} srcpos;

typedef enum {op_modulus, op_plus, op_minus, op_leq, op_lt, op_geq, op_gt, op_eq, op_neq, op_log_neg, op_log_or, op_log_and, op_arith_and, op_arith_or, op_arith_xor, op_left_shift, op_right_shift, op_mult, op_div, op_deref, op_inc, op_dec, op_address, op_sizeof, op_not} operator_t;

typedef struct {
    srcpos lx_srcpos;
    char *string;
    operator_t op;
    int type_stack_count;
} lx_info;

typedef struct {
    int static_size;
    sm_ref control_field;
} dimen_s;

typedef struct dimen_p{
    int dimen_count;
    dimen_s dimens[1];
} *dimen_p;

typedef struct enc_struct {
    int byte_order;
    int float_order;
    int pointer_size;
} *enc_info;

struct cod_exec_struct {
    void *ec;
    void *static_data;
    cod_code gen_code;
};

extern void cod_print_operator_t(operator_t o);

extern void cod_print_srcpos(srcpos pos);

extern void cod_print_dimen_p(dimen_p dimen);

extern void cod_print_enc_info(enc_info enc);

extern void free_enc_info(enc_info enc);

extern void *cod_cg_net(sm_ref net, int ret_type, unsigned int *offset_p,
			cod_code code_descriptor);
extern int cod_sm_get_type(sm_ref node);
extern int cod_expr_is_string(sm_ref expr);
extern int is_control_value(sm_ref expr, sm_ref strct);
extern sm_ref get_complex_type(cod_parse_context context, sm_ref node);
extern int are_compatible_ptrs(sm_ref left, sm_ref right);
extern void cod_add_standard_elements(cod_parse_context context);
extern void cod_process_include(char *name, cod_parse_context context);
extern void cod_add_defined_type(char *id, cod_parse_context context);
typedef struct list_struct *sm_list;
extern void cod_remove_defined_types(cod_parse_context context, int count);
extern void cod_add_decl_to_scope(char *id, sm_ref node, cod_parse_context context);
extern int cod_semanticize_added_decls(cod_parse_context context);
extern void cod_swap_decls_to_standard(cod_parse_context context);
extern void
cod_add_decl_to_parse_context(const char *name, sm_ref item, cod_parse_context context);
extern sm_ref
evaluate_constant_return_expr(cod_parse_context context, sm_ref expr, int *free_result);
