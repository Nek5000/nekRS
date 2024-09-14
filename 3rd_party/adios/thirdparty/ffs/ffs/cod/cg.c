#include "config.h"

#ifdef LINUX_KERNEL_MODULE
#ifndef MODULE
#define MODULE
#endif
#ifndef __KERNEL__
#define __KERNEL__
#endif
#include <linux/kernel.h>
#include <linux/module.h>
#endif

#undef NDEBUG
#include "assert.h"
#ifndef LINUX_KERNEL_MODULE
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdio.h>
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#endif
#ifdef HAVE_DILL_H
#include "dill.h"
#define dill_mark_label_type int
#else
/*stuff*/
#define dill_reg int
#define dill_mark_label_type int
#define dill_stream void *
enum {
    DILL_C,    /* char */
    DILL_UC,   /* unsigned char */
    DILL_S,    /* short */
    DILL_US,   /* unsigned short */
    DILL_I,    /* int */
    DILL_U,    /* unsigned */
    DILL_L,    /* long */
    DILL_UL,   /* unsigned long */
    DILL_P,    /* pointer */
    DILL_F,    /* floating */
    DILL_D,    /* double */
    DILL_V,    /* void */
    DILL_B,    /* block structure */
    DILL_ERR   /* no type */
};
#define dill_type_size(c, t)  0
#define dill_type_align(c, t) 1
#define dill_alloc_label(c, n) 0
#endif
#include "cod.h"
#include "cod_internal.h"
#include "structs.h"
#include "string.h"
#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 2215)
#endif
#include "cod.tab.h"
#ifndef FALSE
#define FALSE 0
#endif

#if defined(_MSC_VER)
#define strdup _strdup
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

typedef struct _opnd_enc_info {
    int is_encoded;
    int byte_swap_on_fetch;
    dill_reg string_base;
} operand_enc_info;

typedef struct _oprnd {
    dill_reg reg;
    int is_addr;
    int in_kernel;
    int offset;
    operand_enc_info enc;
} operand;

#if defined(HAVE_DILL_H)
static void cg_compound_statement(dill_stream s, sm_ref stmt, cod_code descr);
static void cg_decl(dill_stream s, sm_ref decl, cod_code descr);
static void cg_statement(dill_stream s, sm_ref stmt, cod_code descr);
static operand cg_expr(dill_stream s, sm_ref expr, int left, cod_code descr);
static void cg_selection_statement(dill_stream s, sm_ref stmt, cod_code descr);
static void cg_iteration_statement(dill_stream s, sm_ref stmt, cod_code descr);
static void cg_return_statement(dill_stream s, sm_ref stmt, cod_code descr);
static void cg_jump_statement(dill_stream s, sm_ref stmt, cod_code descr);
static void gen_mov(dill_stream s, operand left, dill_reg right, int type);
static void gen_load(dill_stream s, dill_reg left, operand right, int type);
static void gen_encoded_field_load(dill_stream s, dill_reg ret, operand base, int load_type, sm_ref field);
static dill_reg coerce_type(dill_stream s, dill_reg local, int target_type, int item_type);
static int is_comparison_operator(sm_ref expr);
static void cg_branch_if_false(dill_stream s, sm_ref pred, dill_mark_label_type label, cod_code descr, int reverse);
static int is_complex_type(sm_ref expr);
static int is_static_var(sm_ref expr);
extern int cod_sm_get_type(sm_ref node);
#endif
int cg_get_size(dill_stream s, sm_ref node);
extern int is_array(sm_ref expr);

static int inst_count_guess = 0;
static void
init_operand(operand *op) 
{
    op->reg = 0;
    op->is_addr = 0;
    op->in_kernel = 0;
    op->offset = 0;
    op->enc.is_encoded = 0;
}

static int
cg_required_align(dill_stream s, sm_ref node)
{
  int required_align;
  switch (node->node_type) {
  case cod_declaration:
    if (node->node.declaration.sm_complex_type == NULL) {
      required_align = dill_type_align(s, node->node.declaration.cg_type);
    } else {
      required_align = cg_required_align(s, node->node.declaration.sm_complex_type);
    }
    break;
  case cod_field:
    if (node->node.field.sm_complex_type == NULL) {
      required_align = dill_type_align(s, node->node.field.cg_type);
    } else {
      required_align = cg_required_align(s, node->node.field.sm_complex_type);
    }
    break;
  case cod_array_type_decl:
    if (node->node.array_type_decl.sm_complex_element_type == NULL) {
      required_align = dill_type_align(s, node->node.array_type_decl.cg_element_type);
    } else {
      required_align = cg_required_align(s, node->node.array_type_decl.sm_complex_element_type);
    }
    break;
  case cod_struct_type_decl:
    required_align = dill_type_align(s, DILL_D);
    break;
  case cod_reference_type_decl:
    required_align = dill_type_align(s, DILL_P);
    break;
  default:
   required_align = 0; //silence warning
    assert(0);
  }
  return required_align;
}

static void
cg_preprocess(sm_ref node, void *data) {
    cod_code code_descriptor = (cod_code) data;

    switch(node->node_type) {
    case cod_jump_statement:{
	inst_count_guess += 2;
	break;
    }
    case cod_operator: {
	inst_count_guess += 6;	/* 5 probable worst case */
	break;
    }
    case cod_cast: {
	inst_count_guess += 2;
	break;
    }
    case cod_reference_type_decl: {
	break;
    }
    case cod_type_specifier: {
	break;
    }
    case cod_subroutine_call: {
	sm_list args = node->node.subroutine_call.arguments;
	while (args != NULL) {
	    inst_count_guess++;
	    args = args->next;
	}
	inst_count_guess+= 3;		/* call plus a few extra */
	break;
    }
    case cod_selection_statement: {
	inst_count_guess += 8;
	break;
    }
    case cod_constant: {
	inst_count_guess += 3;	/* probably 2 (strings??) */
	break;
    }
    case cod_compound_statement: {
	break;
    }
    case cod_label_statement: {
	dill_stream s = code_descriptor->drisc_context;
	dill_mark_label_type stmt_label = dill_alloc_label(s, node->node.label_statement.name);
	node->node.label_statement.cg_label = stmt_label;
	break;
    }
    case cod_assignment_expression: {
	inst_count_guess += 3;	/* coerce plus move */
	break;
    }
    case cod_declaration: {
	inst_count_guess += 4;
/*	if (decl->node.declaration.static_var) {
	    code_description->static_block_address_register = 1;
	    }*/
	break;
    }
    case cod_element_ref: {
	inst_count_guess += 6;
	break;
    }
    case cod_return_statement: {
	inst_count_guess += 2;
	break;
    }
    case cod_expression_statement: {
	break;
    }
    case cod_field: {
	break;
    }
    case cod_iteration_statement: {
	inst_count_guess += 6;
	break;
    }
    case cod_identifier: {
	break;
    }
    case cod_field_ref: {
	inst_count_guess += 6;
	break;
    }
    case cod_conditional_operator: {
	inst_count_guess += 8;
	break;
    }
    case cod_enumerator: {
	inst_count_guess += 1;
	break;
    }
    case cod_comma_expression:
    case cod_array_type_decl:
    case cod_enum_type_decl:
	break;
    case cod_struct_type_decl:
	if (node->node.struct_type_decl.fields && 
	    (node->node.struct_type_decl.fields->node->node.field.cg_offset == -1)) {
	    /* we need to do structure layout */
	    dill_stream s = code_descriptor->drisc_context;
	    int last_offset = 0;
	    int off_by;
	    sm_list fields = node->node.struct_type_decl.fields;
	    while(fields != NULL) {
		sm_ref field = fields->node;
		int required_align = cg_required_align(s, field);
		int field_size;
		off_by = last_offset % required_align;
		if (off_by != 0) {
		    last_offset += required_align - off_by;
		}
		field->node.field.cg_offset = last_offset;
		if (field->node.field.sm_complex_type != NULL) {
		    field_size = cg_get_size(s, field->node.field.sm_complex_type);
		} else {
		    field_size = dill_type_size(s, field->node.field.cg_type);
		}

		last_offset += field_size;
		field->node.field.cg_size = field_size;
		fields = fields->next;
	    }
	    last_offset += (dill_type_align(s, DILL_D) - (last_offset % dill_type_align(s, DILL_D))) % dill_type_align(s, DILL_D);
	    node->node.struct_type_decl.cg_size = last_offset;
	}
	break;
    case cod_initializer_list:
    case cod_initializer:
    case cod_designator:
	break;
    default:
	cod_print(node);
	assert(FALSE);
    }
}
    
static char *arg_type_str[] = { "c", "uc", "s", "us", "i", "u", "l", "ul", "p", "f", "d", "v", "p", "EC"};

char *
generate_arg_str(sm_ref net)
{
    sm_list decls = NULL;
    char *arg_str = malloc(4);  /* 4 makes valgrind happier */
    signed char *arg_types = malloc(1);
    int arg_count = 0;
    int i;
    if (net->node_type == cod_compound_statement) {
	decls = net->node.compound_statement.decls;
    }
    arg_str[0] = 0;
    while (decls != NULL) {
	sm_ref decl = decls->node;
	
	if (decl->node_type == cod_declaration) {
	    if (decl->node.declaration.param_num != -1) {
		/* parameter */
		int param_num = decl->node.declaration.param_num;
		if (param_num >= arg_count) {
		    arg_types = realloc(arg_types, param_num + 1);
		    memset(arg_types + arg_count, -1, 
			   param_num - arg_count);
		    arg_count = param_num + 1;
		}
		arg_types[param_num] = decl->node.declaration.cg_type;
	    }
	} else if (decl->node_type == cod_array_type_decl) {
	    while(decl->node_type == cod_array_type_decl)
		decl = decl->node.array_type_decl.element_ref;
	    if (decl->node.declaration.param_num != -1) {
		int param_num = decl->node.declaration.param_num;
		if (param_num >= arg_count) {
		    arg_types = realloc(arg_types, param_num + 1);
		    memset(arg_types + arg_count, -1, 
			   param_num - arg_count);
		    arg_count = param_num + 1;
		}
		arg_types[param_num] = DILL_P;
	    }
	}
	    
	decls = decls->next;
    }
    if ((arg_count == 1) && (arg_types[0] == DILL_V)) {
	arg_count = 0;
    }
    for (i = 0; i < arg_count ; i++) {
	size_t len;
	if (arg_types[i] == -1) {
	    printf("Arg %d not declared\n", i);
	    return arg_str;
	}
	len = strlen(arg_str) + 8;
	arg_str = realloc(arg_str, len);
	strncat(arg_str, "%", len);
	strncat(arg_str, arg_type_str[arg_types[i]], len);
    }
    free(arg_types);
    return arg_str;
}

#ifdef HAVE_DILL_H
int
cg_get_size(dill_stream s, sm_ref node) {

    sm_ref ref = node;
    
    switch(ref->node_type) {
    case cod_identifier:
	if(ref->node.identifier.sm_declaration) {
	    return cg_get_size(s, ref->node.identifier.sm_declaration);
	}
	return dill_type_size(s, ref->node.identifier.cg_type);

    case cod_declaration:
	if(ref->node.declaration.sm_complex_type) {
	    return cg_get_size(s, ref->node.declaration.sm_complex_type);
	}
	return dill_type_size(s, ref->node.declaration.cg_type);
	
    case cod_reference_type_decl:
	return dill_type_size(s, DILL_P);

    case cod_field:
	return ref->node.field.cg_size;
		
    case cod_struct_type_decl:{
	if ((ref->node.struct_type_decl.cg_size % dill_type_align(s, DILL_D)) != 0) {
	    int struct_size = ref->node.struct_type_decl.cg_size;
	    struct_size += (dill_type_align(s, DILL_D) - (struct_size % dill_type_align(s, DILL_D))) % dill_type_align(s, DILL_D);
	    ref->node.struct_type_decl.cg_size = struct_size;
	}
	return ref->node.struct_type_decl.cg_size;
    }
    case cod_cast:
	if (ref->node.cast.sm_complex_type != NULL) {
	    return cg_get_size(s, ref->node.cast.sm_complex_type);
	} else {
	    return dill_type_size(s, ref->node.cast.cg_type);
	};
    case cod_array_type_decl: {
	/* handle dynamic here */
	return ref->node.array_type_decl.cg_static_size * ref->node.array_type_decl.cg_element_size;
    }
    case cod_assignment_expression:
    case cod_field_ref:
    case cod_element_ref:
    case cod_subroutine_call:
    case cod_operator: {
	sm_ref typ = get_complex_type(NULL, ref);
	if (typ != NULL) {
	    return cg_get_size(s, typ);
	}
	return dill_type_size(s, cod_sm_get_type(ref));
    }
    case cod_enum_type_decl:
	return dill_type_size(s, DILL_I);
    case cod_constant:
    case cod_iteration_statement:
    case cod_expression_statement:
    case cod_return_statement:
    case cod_compound_statement:
    case cod_type_specifier:
    case cod_selection_statement:
    default:
	/* shouldn't happen */
	assert(FALSE);
    }
    return 0;
}
#else
int cg_get_size(dill_stream s, sm_ref node)
{
    return 0;
}
#endif

static void cg_generate_static_block(dill_stream s, cod_code descr);

extern void *
cod_cg_net(sm_ref net, int ret_type, unsigned int *offset_p, cod_code code_descriptor)
{
    void *init_func;
    static int debug_cg = -1;
    dill_stream s = NULL;
    char *arg_str;
#if defined(HAVE_DILL_H)
    dill_exec_handle handle;
#endif
    if (debug_cg == -1) {
	debug_cg = (int)(intptr_t)(getenv("COD_DEBUG"));
    }
#if defined(HAVE_DILL_H)
    arg_str = generate_arg_str(net);
    s = dill_create_stream();
    code_descriptor->drisc_context = s;
    inst_count_guess = 0;
    code_descriptor->static_size_required = 0;
    code_descriptor->static_block_address_register = 0;
    dill_start_proc(s, "no name", ret_type, arg_str);
    cod_apply(net, cg_preprocess, NULL, NULL, code_descriptor);
    free(arg_str);
    code_descriptor->static_block_address_register = -1;
    cg_compound_statement(s, net, code_descriptor);
    if (code_descriptor->static_block_address_register != -1) {
	dill_begin_prefix_code(s);
	cg_generate_static_block(s, code_descriptor);
    }
    if (debug_cg) {
	printf("Virtual insn dump\n");
	dill_dump(s);
	printf("\n");
    }
    handle = dill_finalize(s);
    code_descriptor->execution_handle = handle;
    init_func = dill_get_fp(handle);
    if (debug_cg) {
	printf("Native insn dump\n");
	dill_dump(s);
	printf("\n");
    }
    code_descriptor->code_memory_block = NULL;
    return init_func;
#else
    printf("No dcg available\n");
    return NULL;
#endif
}

#if defined(HAVE_DILL_H)
static void
cg_decls_stmts(dill_stream s, sm_list list, cod_code descr)
{
    while(list != NULL) {
	sm_ref item = list->node;
	switch(item->node_type) {
	case cod_declaration: 
	case cod_struct_type_decl:
	case cod_array_type_decl:
	case cod_reference_type_decl:
	case cod_constant:
	case cod_enum_type_decl:
	    cg_decl(s, item, descr);
	    break;
	default:
	    cg_statement(s, item, descr);
	    break;
	}
	list = list->next;
    }
}

static void
cg_compound_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    cg_decls_stmts(s, stmt->node.compound_statement.decls, descr);
    cg_decls_stmts(s, stmt->node.compound_statement.statements, descr);
}

static dill_reg 
cg_get_static_block(dill_stream s, cod_code descr)
{
    dill_reg ret;
    if (descr->static_block_address_register != -1) {
	return descr->static_block_address_register;
    } else {
	ret = dill_getreg(s, DILL_P);
	descr->static_block_address_register = ret;
	return ret;
    }
}

static void
cg_generate_static_block(dill_stream s, cod_code descr)
{
    if (descr->has_exec_ctx) {
	dill_reg cod_ctx = dill_param_reg(s, 0);
	dill_reg ec = dill_getreg(s, DILL_P);
	dill_reg ret;
	dill_ldpi(s, ec, cod_ctx, FMOffset(struct cod_exec_struct *, ec));
	ret = dill_scallp(s, (void*)&dill_get_client_data, "dill_get_client_data", "%p%I", ec, 0x23234);
	dill_movp(s, descr->static_block_address_register, ret);
    } else {
	dill_setp(s, descr->static_block_address_register, descr->data);
    }
}

static int
evaluate_constant_expr(dill_stream s, sm_ref expr, long*value)
{
    switch(expr->node_type) {
    case cod_constant: {
	char *val = expr->node.constant.const_val;
	long l;
	if (val[0] == '0') {
	    /* hex or octal */
	    if (val[1] == 'x') {
		/* hex */
		if (sscanf(val+2, "%lx", &l) != 1) 
		    printf("sscanf failed\n");
	    } else {
		if (sscanf(val, "%lo", &l) != 1) 
		    printf("sscanf failed\n");
	    }
	} else {
	    if (sscanf(val, "%ld", &l) != 1)
		printf("sscanf failed\n");
	}
	*value = l;
	return 1;
    }
    case cod_identifier:
	return evaluate_constant_expr(s, expr->node.identifier.sm_declaration, value);
    case cod_declaration:
	if (!expr->node.declaration.const_var) return 0;
	return evaluate_constant_expr(s, expr->node.declaration.init_value, value);
    case cod_operator: {
	long left, right;
	if (expr->node.operator.left != NULL) {
	    if (!evaluate_constant_expr(s, expr->node.operator.left, &left)) return 0;
	}
	if (expr->node.operator.op == op_sizeof) {
	    *value = cg_get_size(s, expr->node.operator.right);
	    return 1;
	}
	if (expr->node.operator.right != NULL) {
	    if (!evaluate_constant_expr(s, expr->node.operator.right, &right)) return 0;
	}
	switch(expr->node.operator.op) {
	case  op_modulus:
	    *value = left % right;
	    break;
	case  op_plus:
	    *value = left + right;
	    break;
	case  op_minus:
	    *value = left - right;
	    break;
	case  op_leq:
	    *value = left <= right;
	    break;
	case  op_lt:
	    *value = left < right;
	    break;
	case  op_geq:
	    *value = left >= right;
	    break;
	case  op_gt:
	    *value = left > right;
	    break;
	case  op_eq:
	    *value = left = right;
	    break;
	case  op_neq:
	    *value = left != right;
	    break;
	case  op_log_or:
	    *value = left || right;
	    break;
	case  op_arith_or:
	    *value = left | right;
	    break;
	case  op_arith_xor:
	    *value = left ^ right;
	    break;
	case  op_log_and:
	    *value = left && right;
	    break;
	case  op_arith_and:
	    *value = left & right;
	    break;
	case  op_mult:
	    *value = left * right;
	    break;
	case  op_div:
	    *value = left / right;
	    break;
	case op_log_neg:
	    *value = !right;
	    break;
	case op_not:
	    *value = ~right;
	    break;
	case op_left_shift:
	    *value = left << right;
	    break;
	case op_right_shift:
	    *value = left >> right;
	    break;
	case  op_deref:
	case  op_address:
	case op_inc:
	case op_dec:
	case op_sizeof:
	    assert(FALSE);
	}
	return 1;
    }
    case cod_cast:
	return evaluate_constant_expr(s, expr->node.cast.expression, value);
    case cod_assignment_expression:
    case cod_field_ref:
    case cod_element_ref:
    case cod_subroutine_call:
	assert(FALSE);
    default:
	assert(FALSE);
    }
    return 0;
}

static void *
add_global(cod_code descr, char *name, int size)
{
    return malloc(size);
}
	
static void
assign_enum_value(dill_stream s, sm_list enums, int *enum_value_p)
{
    sm_ref en;
    if (enums == NULL) return;
    en = enums->node;
    assign_enum_value(s, enums->next, enum_value_p);
    if (en->node.enumerator.const_expression) {
	long tmp;
	assert(evaluate_constant_expr(s, en->node.enumerator.const_expression, &tmp));
	*enum_value_p = tmp;
    }
    en->node.enumerator.enum_value = (*enum_value_p)++;
    return;
}

static void
cg_enum_type_decl(dill_stream s, sm_ref decl, cod_code descr)
{
    int enum_value = 0;
    assign_enum_value(s, decl->node.enum_type_decl.enums, &enum_value);
}

static void
set_dimen_values(dill_stream s, sm_ref top, sm_ref this, int dimension)
{
    long size = -1;
    if (this->node_type != cod_array_type_decl) return;
    if (this->node.array_type_decl.size_expr) {
	evaluate_constant_expr(s, this->node.array_type_decl.size_expr, &size);
	top->node.array_type_decl.dimensions->dimens[dimension].static_size = size;
    }
    set_dimen_values(s, top, this->node.array_type_decl.element_ref, dimension+1);
}

static void
set_complex_type_sizes(dill_stream s, sm_ref typ)
{
    switch(typ->node_type) {
    case cod_array_type_decl: {
	int i, element_count = 1;
	if (typ->node.array_type_decl.dimensions != NULL) {
	    set_dimen_values(s, typ, typ, 0);
	}
	if (typ->node.array_type_decl.sm_complex_element_type) {
	    set_complex_type_sizes(s, typ->node.array_type_decl.sm_complex_element_type);
	}
	typ->node.array_type_decl.cg_element_size = 
	    dill_type_size(s, typ->node.array_type_decl.cg_element_type);
	if (dill_type_align(s, typ->node.array_type_decl.cg_element_type) >
	    typ->node.array_type_decl.cg_element_size) 
	    typ->node.array_type_decl.cg_element_size = 
		dill_type_align(s, typ->node.array_type_decl.cg_element_type);
	if (typ->node.array_type_decl.dimensions) {
	    for (i=0; i < typ->node.array_type_decl.dimensions->dimen_count; i++) {
		element_count *= typ->node.array_type_decl.dimensions->dimens[i].static_size;
	    }
	}
	typ->node.array_type_decl.cg_static_size = element_count;
	break;
    }
    case cod_struct_type_decl:
    case cod_enum_type_decl:
	break;
    case cod_declaration:
	if (typ->node.declaration.sm_complex_type) {
	    set_complex_type_sizes(s, typ->node.declaration.sm_complex_type);
	}
	break;
    case cod_reference_type_decl:
	if (typ->node.reference_type_decl.sm_complex_referenced_type) {
	    set_complex_type_sizes(s, typ->node.reference_type_decl.sm_complex_referenced_type);
	}
	break;
    default:
	printf("unhandled case in set_complex_type_sizes\n");
	cod_print(typ);
	assert(0);
    }
}

static void
evaluate_simple_init_and_assign(dill_stream s, sm_ref init, int cg_type, void *var_base);

static
char *generate_block_init_value(dill_stream s, sm_ref decl)
{
    sm_ref init = decl->node.declaration.init_value;
    sm_ref typ = get_complex_type(NULL, decl);
    if (init->node_type == cod_constant) {
	/* should be string! */
	return strdup(init->node.constant.const_val);
    }
    if (typ->node_type == cod_array_type_decl) {
	int size = cg_get_size(s, decl);
	char *ret = malloc(size), *tmp;
	sm_list init_list;
	memset(ret, 0, size);
	assert(init->node_type == cod_initializer_list);
	init_list = init->node.initializer_list.initializers;
	tmp = ret;
	while (init_list) {
	    sm_ref initializer = init_list->node;
	    evaluate_simple_init_and_assign(s, initializer->node.initializer.initializer, 
					    typ->node.array_type_decl.cg_element_type, tmp);
	    tmp += typ->node.array_type_decl.cg_element_size;
	    init_list = init_list->next;
	}
	return ret;
    } else if (typ->node_type == cod_struct_type_decl) {
	int size = cg_get_size(s, decl);
	char *ret = malloc(size);
	sm_list init_list, fields;
	memset(ret, 0, size);
	assert(init->node_type == cod_initializer_list);
	init_list = init->node.initializer_list.initializers;
	fields = typ->node.struct_type_decl.fields;
	while (init_list) {
	    sm_ref initializer = init_list->node;
	    sm_ref field = fields->node;
	    evaluate_simple_init_and_assign(s, initializer->node.initializer.initializer, 
					    field->node.field.cg_type, ret + field->node.field.cg_offset);
	    init_list = init_list->next;
	    fields = fields->next;
	}
	return ret;
    }
    cod_print(decl->node.declaration.init_value);
    return NULL;
}

static void
evaluate_simple_init_and_assign(dill_stream s, sm_ref init, int cg_type, void *var_base)
{
    sm_ref const_val;
    int free_expr = 0;
    const_val = evaluate_constant_return_expr(NULL, init, &free_expr);
    assert(const_val->node_type == cod_constant);
    if (const_val->node.constant.token == string_constant) {
	/* leaky */
	*(char**)(var_base) = strdup(const_val->node.constant.const_val);
    } else if (const_val->node.constant.token == floating_constant) {
	double f;
	
	sscanf(const_val->node.constant.const_val, "%lf", &f);
	switch (cg_type) {
	case DILL_C:
	    *(char *)(var_base) = (int)f;
	    break;
	case DILL_UC:
	    *(unsigned char *)(var_base) = (int)f;
	    break;
	case DILL_S:
	    *(short *)(var_base) = (int)f;
	    break;
	case DILL_US:
	    *(unsigned short *)(var_base) = (int)f;
	    break;
	case DILL_I:
	    *(int *)(var_base) = (int)f;
	    break;
	case DILL_U:
	    *(unsigned int *)(var_base) = (int)f;
	    break;
	case DILL_L:
	    *(long *)(var_base) = (long)f;
	    break;
	case DILL_UL:
	    *(unsigned long *)(var_base) = (long)f;
	    break;
	case DILL_F:
	    *(float*)(var_base) = (float)f;
	    break;
	case DILL_D:
	    *(double*)(var_base) = (double)f;
	    break;
	default:
	    assert(FALSE);
	}
    } else {
	ssize_t l;
	char *val = const_val->node.constant.const_val;
	if (const_val->node.constant.token == character_constant) {
	    if (*val == 'L') val++ ;
	    if (*val == 'u') val++ ;
	    if (*val == 'U') val++ ;
	    val++ ; /* skip ' */
	    if (*val != '\\') {
		l = *val;
	    } else {
		val++;
		switch(*val) {
		case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
		    if (sscanf(val, "%zo", &l) != 1) 
			printf("octal char sscanf failed %s\n", val);
		    break;
		case 'x':
		    if (sscanf(val+1, "%zx", &l) != 1) 
			printf("hex char sscanf failed, %s\n", val);
		    break;
		case '\\':
		    l = 0134;
		    break;
		case '\'':
		    l = 047;
		    break;
		case '"':
		    l = 042;
		    break;
		case '?':
		    l = 077;
		    break;
		case 'a':
		    l = 07;
		    break;
		case 'b':
		    l = 010;
		    break;
		case 'f':
		    l = 06;
		    break;
		case 'n':
		    l = 012;
		    break;
		case 'r':
		    l = 015;
		    break;
		case 't':
		    l = 011;
		    break;
		default:
		    printf("Bad character constant %s\n", val);
		}
	    }
	} else {
	    if (val[0] == '0') {
		/* hex or octal */
		if (val[1] == 'x') {
		    /* hex */
		    if (sscanf(val+2, "%zx", &l) != 1) 
			printf("sscanf failed\n");
		} else {
		    if (sscanf(val, "%zo", &l) != 1) 
			printf("sscanf failed\n");
		}
	    } else {
		if (sscanf(val, "%zd", &l) != 1) 
		    printf("sscanf failed\n");
	    }
	}
	switch (cg_type) {
	case DILL_C:
	    *(char *)(var_base) = (char)l;
	    break;
	case DILL_UC:
	    *(unsigned char *)(var_base) = (unsigned char)l;
	    break;
	case DILL_S:
	    *(short *)(var_base) = (short)l;
	    break;
	case DILL_US:
	    *(unsigned short *)(var_base) = (unsigned short)l;
	    break;
	case DILL_I:
	    *(int *)(var_base) = (int)l;
	    break;
	case DILL_U:
	    *(unsigned int *)(var_base) = (unsigned int)l;
	    break;
	case DILL_L:
	    *(ssize_t *)(var_base) = l;
	    break;
	case DILL_UL:
	    *(size_t *)(var_base) = l;
	    break;
	case DILL_P:
	    *(void **)(var_base) = (void*)l;
	    break;
	case DILL_F:
	    *(float*)(var_base) = (float)l;
	    break;
	case DILL_D:
	    *(double*)(var_base) = (double)l;
	    break;
	default:
	    assert(FALSE);
	}
	if (free_expr) cod_free(const_val);
    }
}

static void 
cg_decl(dill_stream s, sm_ref decl, cod_code descr)
{
    dill_reg lvar = -1;
    switch(decl->node_type) {
    case cod_declaration: {
	void *var_base = NULL;
	sm_ref ctype = decl->node.declaration.sm_complex_type;
	int is_block_type = 0;
	if ((ctype != NULL) && ((ctype->node_type == cod_array_type_decl) || (ctype->node_type == cod_struct_type_decl))) {
	    is_block_type = 1;
	}
	if (ctype && (decl->node.declaration.is_typedef)) {
	    cg_decl(s, decl->node.declaration.sm_complex_type, descr);
	}
	if (ctype) {
	    set_complex_type_sizes(s, ctype);
	}
	if (decl->node.declaration.const_var && !ctype) {
	    /* Nothing to do for SIMPLE const vars.  We'll CG them as consts on reference. */
	    return;
	}
	if (decl->node.declaration.is_extern) {
	    if (!decl->node.declaration.cg_address) {
		/* GLOBAL.  Need to allocate */
		int size = cg_get_size(s, decl);
		void *extern_value = add_global(descr, decl->node.declaration.id, size);
		decl->node.declaration.cg_address = extern_value;
		decl->node.declaration.is_extern = 1;
		decl->node.declaration.static_var = 0;
		var_base = extern_value;
	    } else {
		/* if extern and already has address, we've initialized it already */
		/* so, nothing remains to be done and we're done here... */
		return;
	    }
	}
	if (decl->node.declaration.static_var && !is_block_type) {
	    /* do SIMPLE statics */
	    decl->node.declaration.cg_address = 
		(void*)(ssize_t)descr->static_size_required;
	    descr->static_size_required += 8; 
	    if (descr->data == NULL) {
		descr->data = malloc(descr->static_size_required);
	    } else {
		descr->data = realloc(descr->data, descr->static_size_required);
	    }
	    var_base = (char*)descr->data + (intptr_t)decl->node.declaration.cg_address;
	}
	if (var_base) {
	    if (decl->node.declaration.init_value == NULL) {
		/* init to zero's */
		memset(var_base, 0, cg_get_size(s, decl));
	    } else {
		if (is_block_type) {
		    /* this is really an array or structure being initialized */
		    char * init_value = generate_block_init_value(s, decl);
		    memcpy(var_base, init_value, cg_get_size(s, decl));
		    free(init_value);
		} else {
		    evaluate_simple_init_and_assign(s, decl->node.declaration.init_value, decl->node.declaration.cg_type, var_base);
		}
	    }
	    return;
	}
	if (decl->node.declaration.param_num != -1) {
	    /* PARAMETER */
	    sm_ref typ = decl->node.declaration.sm_complex_type;
	    if (decl->node.declaration.cg_type == DILL_V) {
		return;
	    }
	    lvar = dill_param_reg(s, decl->node.declaration.param_num);
	    if (typ && (typ->node_type == cod_struct_type_decl) &&
	        (typ->node.struct_type_decl.encode_info != NULL)) {
		/* for encoded parameters, skip header */
		dill_addpi(s, lvar, lvar, 24);
	    }
	    if ((ctype != NULL) && (ctype->node_type == cod_array_type_decl)){
		/* array param */
		long i;
		char *val;
		sm_ref const_val = ctype->node.array_type_decl.size_expr;
		if (const_val != NULL) {
		    assert(const_val->node_type == cod_constant);
		    val = const_val->node.constant.const_val;
		    if (val[0] == '0') {
			/* hex or octal */
			if (val[1] == 'x') {
			    /* hex */
			    if (sscanf(val+2, "%lx", &i) != 1) 
				printf("sscanf failed\n");
			} else {
			    if (sscanf(val, "%lo", &i) != 1) 
				printf("sscanf failed\n");
			}
		    } else {
			if (sscanf(val, "%ld", &i) != 1) 
			    printf("sscanf failed\n");
		    }
		    ctype->node.array_type_decl.cg_static_size = i;
		}
	    }
	} else {
	    if (decl->node.declaration.is_subroutine) {
		lvar = 0;
	    } else {
		if ((ctype != NULL) && (ctype->node_type == cod_array_type_decl)){
		    /* array */
		    if (decl->node.declaration.static_var) {
			decl->node.declaration.cg_address = 
			    (void*)(ssize_t)descr->static_size_required;
			descr->static_size_required += ctype->node.array_type_decl.cg_static_size * ctype->node.array_type_decl.cg_element_size;
			if (descr->data == NULL) {
			    descr->data = malloc(descr->static_size_required);
			} else {
			    descr->data = realloc(descr->data, descr->static_size_required);
			}
			var_base = (char*)descr->data + (intptr_t)decl->node.declaration.cg_address;
			lvar = dill_getreg(s, DILL_P);
			dill_addpi(s, lvar, cg_get_static_block(s, descr), (intptr_t)decl->node.declaration.cg_address);  /* op_i_addpi */
		    } else {
			dill_reg block =
			    dill_getvblock(s, ctype->node.array_type_decl.cg_static_size * ctype->node.array_type_decl.cg_element_size);
			lvar = dill_getreg(s, DILL_P);
			dill_virtual_lea(s, lvar, block);	/* op_i_leai */
		    }
		} else if (!is_block_type) {
		    if (decl->node.declaration.static_var) {
			decl->node.declaration.cg_address = 
			    (void*)(ssize_t)descr->static_size_required;
			descr->static_size_required += cg_get_size(s, decl);
			if (descr->data == NULL) {
			    descr->data = malloc(descr->static_size_required);
			} else {
			    descr->data = realloc(descr->data, descr->static_size_required);
			}
			var_base = (char*)descr->data + (intptr_t)decl->node.declaration.cg_address;
		    } else if ((decl->node.declaration.addr_taken) || 
			       (ctype && (ctype->node_type == cod_struct_type_decl))) {
		        /* make sure it's in memory if its address is taken */
			dill_reg block = dill_getvblock(s, 8);
			lvar = dill_getreg(s, DILL_P);
			dill_virtual_lea(s, lvar, block);	/* op_i_leai */
		    } else {
		        lvar = dill_getreg(s, decl->node.declaration.cg_type);
		    }
		} else {
		    sm_ref struct_type = decl->node.declaration.sm_complex_type;
		    assert(struct_type->node_type == cod_struct_type_decl);
		    if (decl->node.declaration.static_var) {
			decl->node.declaration.cg_address = 
			    (void*)(ssize_t)descr->static_size_required;
			descr->static_size_required += cg_get_size(s, decl);
			if (descr->data == NULL) {
			    descr->data = malloc(descr->static_size_required);
			} else {
			    descr->data = realloc(descr->data, descr->static_size_required);
			}
			var_base = (char*)descr->data + (intptr_t)decl->node.declaration.cg_address;
		    } else {
			if ((struct_type->node.struct_type_decl.cg_size % dill_type_align(s, DILL_D)) != 0) {
			    int struct_size = struct_type->node.struct_type_decl.cg_size;
			    struct_size += (dill_type_align(s, DILL_D) - (struct_size % dill_type_align(s, DILL_D))) % dill_type_align(s, DILL_D);
			    struct_type->node.struct_type_decl.cg_size = struct_size;
			}
			dill_reg block =
			    dill_getvblock(s, struct_type->node.struct_type_decl.cg_size);
			lvar = dill_getreg(s, DILL_P);
			dill_virtual_lea(s, lvar, block);	/* op_i_leai */
		    }
		}
	    }
	}
	if (!decl->node.declaration.is_typedef) {
	    decl->node.declaration.cg_oprnd = lvar;
	}
	if (decl->node.declaration.init_value != NULL) {
	    operand left, right;
	    int assign_type = cod_sm_get_type(decl);
	    init_operand(&left);
	    if (is_block_type) {
		/* this is really an array or structure being initialized */
		char * init_value = generate_block_init_value(s, decl);
		/* leak init_value */
		(void) dill_scallv(s, (void*)memcpy, "memcpy", "%p%P%I", lvar, init_value, cg_get_size(s, decl));
	    } else {
		right = cg_expr(s, decl->node.declaration.init_value, 0, descr);
		if (right.is_addr != 0) {
		    fprintf(stderr, "Failure in processing declaration init value, expr is :\n");
		    cod_print(decl->node.declaration.init_value);
		    fprintf(stderr, "Please report this bug and include as much information as possible for reproduction\n");
		    exit(1);
		}
		right.reg = coerce_type(s, right.reg, assign_type, 
					cod_sm_get_type(decl->node.declaration.init_value));
		if (decl->node.declaration.addr_taken) {
		    init_operand(&left);
		    left.reg = lvar;
		    left.is_addr = 1;
		    left.in_kernel = 0; /* we have no access to kernel structures */
		    left.offset = 0;
		} else {
		    left.reg = lvar;
		    left.is_addr = 0;
		    left.in_kernel = 0; /* we have no access to kernel structures */
		    left.offset = 0;
		}
		gen_mov(s, left, right.reg, assign_type);
	    }
	} else if (is_block_type &&
		   (decl->node.declaration.param_num == -1)) {
	    /* init structure to zero's */
	    if (decl->node.declaration.static_var) {
		memset(var_base, 0, cg_get_size(s, decl));
	    } else {
		(void) dill_scallv(s, (void*)memset, "memset", "%p%I%I", lvar, 0, cg_get_size(s, decl));
	    }
	}
	break;
    }
    case cod_struct_type_decl: {
	    sm_list fields = decl->node.struct_type_decl.fields;
	    int do_cg = 0, last_offset = 0;
	    int off_by;
	    while (fields != NULL) {
		if (fields->node->node.field.cg_size <= 0) {
		    do_cg = 1;
		}
		fields = fields->next;
	    }
	    if (do_cg == 0) return;
	    fields = decl->node.struct_type_decl.fields;
	    while(fields != NULL) {
		sm_ref field = fields->node;
		int field_size = cg_get_size(s, field);
		int field_align = cg_required_align(s, field);
		off_by = last_offset % field_align;
		if (off_by != 0) {
		    last_offset += field_align - off_by;
		}
		field->node.field.cg_offset = last_offset;
		field->node.field.cg_size = field_size;
		last_offset += field_size;
		fields = fields->next;
	    }
	    last_offset += (dill_type_align(s, DILL_D) - (last_offset % dill_type_align(s, DILL_D))) % dill_type_align(s, DILL_D);
	    decl->node.struct_type_decl.cg_size = last_offset;
        }
	break;
    case cod_enum_type_decl:
	cg_enum_type_decl(s, decl, descr);
	break;
    case cod_reference_type_decl:
	break;
    case cod_array_type_decl:
	if (decl->node.array_type_decl.dimensions != NULL) {
	    set_dimen_values(s, decl, decl, 0);
	}
	cg_decl(s, decl->node.array_type_decl.element_ref, descr);
	break;
    case cod_constant:
        /* don't need to do anything */
        break;
    default:
	printf("Unexpected case in cg_decl\n");
	cod_print(decl);
    }
}

static void
cg_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    if (!stmt) return;
    switch(stmt->node_type) {
    case cod_selection_statement:
	cg_selection_statement(s, stmt, descr);
	break;
    case cod_iteration_statement:
	cg_iteration_statement(s, stmt, descr);
	break;
    case cod_expression_statement: 
	(void) cg_expr(s, stmt->node.expression_statement.expression, 0, descr);
	break;
    case cod_compound_statement:
	cg_compound_statement(s, stmt, descr);
	break;
    case cod_return_statement:
	cg_return_statement(s, stmt, descr);
	break;
    case cod_jump_statement:
	cg_jump_statement(s, stmt, descr);
	break;
    case cod_label_statement: {
	dill_mark_label(s, stmt->node.label_statement.cg_label);
	cg_statement(s, stmt->node.label_statement.statement, descr);
	break;
    }
    default:
	printf("unhandled case in cg_statement\n");
    }
}

static int
cod_streq(char *str1, char *str2)
{
    int ret;
    if ((str1 == NULL) && (str2 == NULL)) return 0;
    if ((str1 == NULL) || (str2 == NULL)) return 1;
    ret = (strcmp(str1, str2) != 0);
/*    printf("Comparing %s and %s, returning %d\n", str1, str2, ret);*/
    return ret;
}

static char*
cod_strdup(char *str)
{
    if (str == NULL) return NULL;
/*    printf("Duplicating str %s\n", str);*/
    return strdup(str);
}

static void
gen_bz(dill_stream s, int conditional, int target_label, int op_type);
static void
gen_bnz(dill_stream s, int conditional, int target_label, int op_type);

/*
  Operator_prep generally just code generates the left and right operands to an operator, 
  however it does more stuff for some operators.  In particular, because of the need to 
  generate short-circuit code for logical_or and logical_and (I.E don't evaluate the 
  right operand), we generate the complete operation here and return the result as *lp.  
  The generated code looks like this:
AND case

eval left
  set result 0
  if (left eq 0) goto return_result
eval right
  if right eq 0 goto return_result
  set result 1
label return_result


OR case
eval left
  set result 1
  if left ne 0 goto return_result
eval right
  if right ne 0 goto return_result
  set result 0
label return_result

*/
static void
generate_short_circuit(dill_stream s, dill_reg operand, dill_reg *result_p, operator_t op, 
		       int end_label, int op_type, int right);

static void
operator_prep(dill_stream s, sm_ref expr, dill_reg *rp, dill_reg *lp, cod_code descr)
{
    operand right_op, left_op;
    dill_reg tmp_right = 0, right = 0, left = 0, result;
    int op_type = expr->node.operator.operation_type;
    int string_op = FALSE;
    dill_mark_label_type short_circuit = 0;
    if (expr->node.operator.op == op_address) {
	right_op = cg_expr(s, expr->node.operator.right, 1, descr);
	assert((right_op.is_addr == 1) || (cod_sm_get_type(expr->node.operator.right) == DILL_B));
	if (right_op.offset != 0) {
	    dill_reg result = dill_getreg(s, DILL_P);
	    dill_addpi(s, result, right_op.reg, right_op.offset);
	    right_op.reg = result;
	}
	*rp = right_op.reg;
	return;
    }
	
    if (expr->node.operator.left != NULL) {
	int left_cg_type = cod_sm_get_type(expr->node.operator.left);
	left_op = cg_expr(s, expr->node.operator.left, 0, descr);
	if (left_op.is_addr && is_array(expr->node.operator.left)) {
	    if (left_op.offset != 0) {
		dill_reg result = dill_getreg(s, DILL_P);
		dill_addpi(s, result, left_op.reg, left_op.offset);
		left_op.reg = result;
	    }
	    left_op.is_addr = 0;
	}

	if (left_op.is_addr != 0) {
	    fprintf(stderr, "Failure in processing LHS of operator, problem expr is :\n");
	    cod_print(expr->node.operator.left);
	    fprintf(stderr, "Please report this bug and include as much information as possible for reproduction\n");
	    assert(left_op.is_addr == 0);
	}
	switch(left_cg_type) {
	case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	    /* do integer promotion */
	    left = coerce_type(s, left_op.reg, DILL_I, left_cg_type);
	    left_cg_type = DILL_I;
	    left_op.reg = left;
	}
	left = left_op.reg;
	if (!string_op && op_type != DILL_P) {
	    left = coerce_type(s, left_op.reg, op_type, left_cg_type);
	}
    }
    if ((expr->node.operator.op == op_log_and) || (expr->node.operator.op == op_log_or)) {
	result = dill_getreg(s, op_type);
	short_circuit = dill_alloc_label(s, "end_of_and/or");
	generate_short_circuit(s, left, &result, expr->node.operator.op, short_circuit, op_type, 0);
    }
    if (expr->node.operator.right != NULL) {
	int right_cg_type = cod_sm_get_type(expr->node.operator.right);
	string_op = cod_expr_is_string(expr->node.operator.right);

	right_op = cg_expr(s, expr->node.operator.right, 0, descr);
	if (right_op.is_addr != 0) {
	    fprintf(stderr, "Failure in processing RHS of operator, expr is :\n");
	    cod_print(expr->node.operator.right);
	    fprintf(stderr, "Please report this bug and include as much information as possible for reproduction\n");
	    assert(right_op.is_addr == 0);
	}
	switch(right_cg_type) {
	case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	    /* do integer promotion */
	    right = coerce_type(s, right_op.reg, DILL_I, right_cg_type);
	    right_cg_type = DILL_I;
	    right_op.reg = right;
	}

	tmp_right = right_op.reg;
	if (!string_op && op_type != DILL_P) 
	    tmp_right = coerce_type(s, right_op.reg, op_type, right_cg_type);
	right = tmp_right;
	if ((expr->node.operator.op == op_log_and) || (expr->node.operator.op == op_log_or)) {
	    generate_short_circuit(s, right, &result, expr->node.operator.op, short_circuit, op_type, 1);
	    dill_mark_label(s, short_circuit);
	    /* return result as left operand */
	    left = result;
	}
    }
    
    *rp = right;
    *lp = left;
}

static void
generate_short_circuit(dill_stream s, dill_reg operand, dill_reg *result_p, operator_t op, 
		       int target_label, int op_type, int right)
{
    if (!right) { /* left (first evaluated) operand, prepare result */
	long init_result_val;
	if (op == op_log_or) {
	    init_result_val = 1;
	} else {
	    init_result_val = 0;
	}
	dill_piset(s, op_type, *result_p, init_result_val);
    }
    if (op == op_log_or) {
	gen_bnz(s, operand, target_label, op_type);
    } else {
	gen_bz(s, operand, target_label, op_type);
    }
    if (right) { /* right (last evaluated) operand, finalize if we haven't skipped forward */
	long final_result_val;
	if (op == op_log_or) {
	    final_result_val = 0;
	} else {
	    final_result_val = 1;
	}
	dill_piset(s, op_type, *result_p, final_result_val);
    }
}

#ifndef HAVE_DILL_H
static void
dill_pmov(dill_stream s, int cg_type, dill_reg dest, dill_reg src)
{
    switch(cg_type) {
    case DILL_C:
    case DILL_UC:
    case DILL_S:
    case DILL_US:
    case DILL_I:
	dill_movi(s, dest, src);	/* op_i_movi */
	break;
    case DILL_U:
	dill_movu(s, dest, src);	/* op_i_movu */
	break;
    case DILL_L:
	dill_movl(s, dest, src);	/* op_i_movl */
	break;
    case DILL_UL:
	dill_movul(s, dest, src);	/* op_i_movul */
	break;
    case DILL_P:
	dill_movp(s, dest, src);	/* op_i_movp */
	break;
    case DILL_D:
	dill_movd(s, dest, src);	/* op_i_movd */
	break;
    case DILL_F:
	dill_movf(s, dest, src);	/* op_i_movf */
	break;
    default:
	fprintf(stderr, "unhandled case in dill_pmov\n");
	break;
    }
}
#endif

static void
poly_addi(dill_stream s, int op_type, dill_reg result, dill_reg left,
	  long imm)
{
    switch (op_type) {
    case DILL_I:
	dill_addii(s, result, left, imm);	/* op_i_addii */
	break;
    case DILL_U:
	dill_addui(s, result, left, imm);	/* op_i_addui */
	break;
    case DILL_L:
	dill_addli(s, result, left, imm);	/* op_i_addli */
	break;
    case DILL_UL:
	dill_adduli(s, result, left, imm);	/* op_i_adduli */
	break;
    case DILL_P:
	dill_addpi(s, result, left, imm);	/* op_i_addpi */
	break;
    case DILL_D: {
	dill_reg one = 	dill_getreg(s, op_type);
	dill_setd(s, one, 1.0);
	dill_Padd(s, op_type, result, left, one);
	break;
    }
    case DILL_F: {
	dill_reg one = 	dill_getreg(s, op_type);
	dill_setf(s, one, 1.0);
	dill_Padd(s, op_type, result, left, one);
	break;
    }
    default:
	assert(0);
    }
}

static void
poly_subi(dill_stream s, int op_type, dill_reg result, dill_reg left,
	  long imm)
{
    switch (op_type) {
    case DILL_I:
	dill_subii(s, result, left, imm);	/* op_i_subii */
	break;
    case DILL_U:
	dill_subui(s, result, left, imm);	/* op_i_subui */
	break;
    case DILL_L:
	dill_subli(s, result, left, imm);	/* op_i_subli */
	break;
    case DILL_UL:
	dill_subuli(s, result, left, imm);	/* op_i_subuli */
	break;
    case DILL_P:
	dill_subpi(s, result, left, imm);	/* op_i_subpi */
	break;
    case DILL_D: {
	dill_reg one = 	dill_getreg(s, op_type);
	dill_setd(s, one, 1.0);
	dill_Psub(s, op_type, result, left, one);
	break;
    }
    case DILL_F: {
	dill_reg one = 	dill_getreg(s, op_type);
	dill_setf(s, one, 1.0);
	dill_Psub(s, op_type, result, left, one);
	break;
    }
    default:
	assert(0);
    }
}

static operand
handle_incdec(dill_stream s, sm_ref expr, cod_code descr)
{
    dill_reg right = 0, opnd, left = 0;
    operand ret_op, store_op;
    operand right_op, left_op;
    int op_type = expr->node.operator.operation_type;
    dill_reg result = dill_getreg(s, expr->node.operator.result_type);
    operator_t op = expr->node.operator.op;

    sm_ref ref = expr->node.operator.right;
    sm_ref typ;
    int size = 1;
    if (ref == NULL) {
	ref = expr->node.operator.left;
    }
    
    if (expr->node.operator.right != NULL) {
	right_op = cg_expr(s, expr->node.operator.right, 1, descr);
	if (right_op.is_addr == 1) {
	    right = dill_getreg(s, op_type);
	    gen_load(s, right, right_op, op_type);
	} else {
	    right = right_op.reg;
	}
    }
    if (expr->node.operator.left != NULL) {
	left_op = cg_expr(s, expr->node.operator.left, 1, descr);
	if (left_op.is_addr == 1) {
	    left = dill_getreg(s, op_type);
	    gen_load(s, left, left_op, op_type);
	} else {
	    left = left_op.reg;
	}
    }

    typ = get_complex_type(NULL, ref);
    
    if(typ && typ->node_type == cod_reference_type_decl) {
	if(typ->node.reference_type_decl.sm_complex_referenced_type) {
	    size = cg_get_size(s, typ->node.reference_type_decl.sm_complex_referenced_type);
	} else {
	    size = dill_type_size(s, typ->node.reference_type_decl.cg_referenced_type);
	}
    }
    
    if (expr->node.operator.right == NULL) {
	/* post-increment/decrement, grab the result first */
	dill_pmov(s, op_type, result, left);
	opnd = left;
	store_op = left_op;
    } else {
	opnd = right;
	store_op = right_op;
    }
    if (op == op_inc) {
	poly_addi(s, op_type, opnd, opnd, size);
    } else {
	poly_subi(s, op_type, opnd, opnd, size);
    }
    if (expr->node.operator.left == NULL) {
	/* pre-increment/decrement, grab the result last */
	dill_pmov(s, op_type, result, right);	/* op_i_movp */
    }
    if (store_op.is_addr == 1) {
	/* store the result home either way if we didn't have an assignable */
	gen_mov(s, store_op, opnd, op_type);
    }
    init_operand(&ret_op);
    ret_op.is_addr = 0;
    ret_op.reg = result;
    return ret_op;
}	

static operand
execute_operator_cg(dill_stream s, operator_t op, int op_type, dill_reg result, dill_reg left, dill_reg right, sm_ref sm_left, sm_ref sm_right, int need_assignable, int result_type);

static operand
cg_operator(dill_stream s, sm_ref expr, int need_assignable, cod_code descr)
{
    dill_reg right, left = 0;
    operand ret_op;
    int op_type = expr->node.operator.operation_type;
    dill_reg result = dill_getreg(s, expr->node.operator.result_type);
    operator_t op = expr->node.operator.op;

    init_operand(&ret_op);
    ret_op.is_addr = 0;
    ret_op.offset = 0;

    if (is_comparison_operator(expr)) {
	dill_mark_label_type false_label = dill_alloc_label(s, "compare start");
	dill_seti(s, result, 0);	/* op_i_seti */
	cg_branch_if_false(s, expr, false_label, descr, 0);
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, false_label);
	ret_op.reg = result;
	return ret_op;
    }

    if ((op == op_inc) || (op == op_dec)) {
	assert(!need_assignable);  /* if we do, this is wrong */
	return handle_incdec(s, expr, descr);
    }
    if (op == op_sizeof) {
	int size = cg_get_size(s, expr->node.operator.right);
	dill_seti(s, result, size);
	ret_op.reg = result;
	return ret_op;
    }

    operator_prep(s, expr, &right, &left, descr);

    if (op == op_log_neg) {
	op = op_neq;
	left = dill_getreg(s, op_type);
	switch (op_type) {
	case DILL_I: case DILL_U: case DILL_L: case DILL_UL: case DILL_P:
	    dill_piset(s, op_type, left, 0); break;
	case DILL_F: dill_setf(s, left, 0.0); break;	/* op_i_setf */
	case DILL_D: dill_setd(s, left, 0.0); break;	/* op_i_setd */
	}
    }
    if ((expr->node.operator.op == op_log_and) || (expr->node.operator.op == op_log_or)) {
	/* everything done by operator_prep */
	ret_op.reg = left;
	return ret_op;
    }
    ret_op = 
	execute_operator_cg(s, op, op_type, result, left, right, expr->node.operator.left,
			    expr->node.operator.right, need_assignable, expr->node.operator.result_type);
    return ret_op;
}

static operand
execute_operator_cg(dill_stream s, operator_t op, int op_type, dill_reg result, dill_reg left, dill_reg right, sm_ref sm_left, sm_ref sm_right, int need_assignable, int result_type)
{
    operand ret_op;
    init_operand(&ret_op);
    ret_op.is_addr = 0;
    ret_op.offset = 0;
    ret_op.reg = result;

    int string_op = FALSE;
    if (sm_right) {
	string_op = cod_expr_is_string(sm_right);
    }
    switch(op) {
    case  op_modulus:
        dill_Pmod(s, op_type, result, left, right);
	break;
    case  op_plus:
	if (sm_left == NULL) {
	    ret_op.reg = right;
	    return ret_op;
	}
	switch (op_type) {
	case DILL_P:
	{
	    sm_ref ptr, arg;
	    sm_ref typ;
	    int size;
	    dill_reg opPtr, opArg;
	    dill_reg offset = dill_getreg(s, DILL_UL);	    

	    /* Figure out which arg is ptr and which is integral */
	    if((typ = get_complex_type(NULL, sm_left))) {
		if((typ->node_type == cod_reference_type_decl) || is_array(sm_left)) {
		    ptr = sm_left;  opPtr = left;
		    arg = sm_right; opArg = right;
		} else {
		    arg = sm_left;  opArg = left;
		    ptr = sm_right; opPtr = right;
		}
	    } else {
		if(cod_sm_get_type(sm_left) == DILL_P) {
		    ptr = sm_left;  opPtr = left; 
		    arg = sm_right; opArg = right;
		} else {
		    arg = sm_left;  opArg = left; 
		    ptr = sm_right; opPtr = right;
		}
	    }

	    /* Get the size of the referenced type */
	    typ = get_complex_type(NULL, ptr);
	    if (typ) {
		if (typ->node_type == cod_reference_type_decl) {
		    if(typ->node.reference_type_decl.sm_complex_referenced_type) {
			if (is_array(typ->node.reference_type_decl.sm_complex_referenced_type)) {
			    sm_ref arr = typ->node.reference_type_decl.sm_complex_referenced_type;
			    size = arr->node.array_type_decl.cg_element_size;
			} else {
			    size = cg_get_size(s, typ->node.reference_type_decl.sm_complex_referenced_type);
			}
		    } else {
			size = dill_type_size(s, typ->node.reference_type_decl.cg_referenced_type);
		    }
		} else if (typ->node_type == cod_array_type_decl) {
		    size = typ->node.array_type_decl.cg_element_size;
		} else {
		    size = 1; // avoid warning
		    assert(0);
		}
	    } else if ((ptr->node_type == cod_operator) && (ptr->node.operator.op == op_address)) {
		size = cg_get_size(s, ptr->node.operator.right);
	    } else {
		if (cod_expr_is_string(ptr)) {
		    size = 1;
		} else {
		    size = 1;  // silence warning
		    assert(FALSE);
		}
	    }
	    if (size != 1) {
		switch(cod_sm_get_type(arg)) {
		case DILL_I:
		    dill_mulii(s, offset, opArg, size);  /* op_i_mulii */
		    dill_cvi2l(s, offset, offset);       /* op_i_cvi2l */
		    break;
		case DILL_U:
		    dill_mului(s, offset, opArg, size);  /* op_i_mului */
		    dill_cvu2l(s, offset, offset);       /* op_i_cvi2l */
		    break;
		case DILL_L:
		    dill_mulli(s, offset, opArg, size);  /* op_i_mulli */
		    break;
		case DILL_UL:
		    dill_mululi(s, offset, opArg, size); /* op_i_mululi */
		    break;
		}
	    } else {
		offset = opArg;
	    }
	    dill_addp(s, result, opPtr, offset);           /* op_i_addp */ 
	}
	break;
	default:
	    dill_Padd(s, op_type, result, left, right);
	    break;
	}
	break;
    case  op_minus:
	if (sm_left == NULL) {
	    left = dill_getreg(s, op_type);
	    switch (op_type) {
	    case DILL_I: dill_seti(s, left, 0); break;	/* op_i_seti */
	    case DILL_U: dill_setu(s, left, 0); break;	/* op_i_setu */
	    case DILL_L: dill_setl(s, left, 0); break;	/* op_i_setl */
	    case DILL_UL: dill_setul(s, left, 0); break;	/* op_i_setul */
	    case DILL_F: dill_setf(s, left, 0.0); break;	/* op_i_setf */
	    case DILL_D: dill_setd(s, left, 0.0); break;	/* op_i_setd */
            case DILL_P:
                /* Unary minus not applicable to pointers */
                assert(FALSE);
                break;
	    }
	}
	switch (op_type) {
	case DILL_P:
	{
	    sm_ref ltyp, rtyp;
	    char   lIsPtr, rIsPtr;
	    int size;
	    
	    /* Figure out if left arg is a pointer or integral */
	    if((ltyp = get_complex_type(NULL, sm_left))) {
		if((ltyp->node_type == cod_reference_type_decl) || is_array(sm_left)) {
		    lIsPtr = 1;
		} else {
		    lIsPtr = 0;
		}
	    } else {
		if(cod_sm_get_type(sm_left) == DILL_P) {
		    lIsPtr = 1;
		} else {
		    lIsPtr = 0;
		}
	    }

	    /* Figure out if right arg is a pointer or integral */
	    if((rtyp = get_complex_type(NULL, sm_right))) {
		if((rtyp->node_type == cod_reference_type_decl) || is_array(sm_right)) {
		    rIsPtr = 1;
		} else {
		    rIsPtr = 0;
		}
	    } else {
		if(cod_sm_get_type(sm_right) == DILL_P) {
		    rIsPtr = 1;
		} else {
		    rIsPtr = 0;
		}
	    }
	    
	    /*
	     * The semanticizer would only allow op_minus w/ a pointer arg if:
	     * (1) it is between two pointers referencing compatible types
	     * (2) it is between a pointer and an integral type
	     */
	    if(lIsPtr && rIsPtr) {
		/* Case #1: binary minus between two compatible pointers */

		/* Figure out referenced type size */
		if(ltyp->node.reference_type_decl.sm_complex_referenced_type) {
		    size = cg_get_size(s, ltyp->node.reference_type_decl.sm_complex_referenced_type);
		} else if (ltyp->node_type == cod_array_type_decl) {
		    size = ltyp->node.array_type_decl.cg_element_size;
		} else {
		    size = dill_type_size(s, ltyp->node.reference_type_decl.cg_referenced_type);
		}
		
		/* Perform subtraction */
		dill_subp(s, result, left, right);           /* op_i_subp   */

		if(size != 1) {
		    dill_divli(s, result, result, size);     /* op_i_divli */
		}
	    } else {
		/* Case #2: binary minus between a pointer and an integral */
		sm_ref ptr, arg;
		sm_ref typ;
		dill_reg opPtr, opArg;
		dill_reg offset = dill_getreg(s, DILL_UL);	    

		if(lIsPtr) {
		    /* Left arg is pointer, right arg is integral */ 
		    ptr = sm_left;  opPtr = left;
		    arg = sm_right; opArg = right;
		} else {
		    /* Left arg is integral, right arg is pointer */ 
		    arg = sm_left;  opArg = left;
		    ptr = sm_right; opPtr = right;
		}

		/* Get the size of the referenced type */
		typ = get_complex_type(NULL, ptr);
		if (typ) {
		    if(typ->node.reference_type_decl.sm_complex_referenced_type) {
			size = cg_get_size(s, typ->node.reference_type_decl.sm_complex_referenced_type);
		    } else if (typ->node_type == cod_array_type_decl) {
			size = typ->node.array_type_decl.cg_element_size;
		    } else {
			size = dill_type_size(s, typ->node.reference_type_decl.cg_referenced_type);
		    }
		} else if ((ptr->node_type == cod_operator) && (ptr->node.operator.op == op_address)) {
		    size = cg_get_size(s, ptr->node.operator.right);
		} else {
		    if (cod_expr_is_string(ptr)) {
			size = 1;
		    } else {
			size = 1; // avoid warning
			assert(FALSE);
		    }
		}
		if (size != 1) {
		    switch(cod_sm_get_type(arg)) {
		    case DILL_I:
			dill_mulii(s, offset, opArg, size);  /* op_i_mulii */
			dill_cvi2l(s, offset, offset);       /* op_i_cvi2l */
			break;
		    case DILL_U:
			dill_mului(s, offset, opArg, size);  /* op_i_mului */
			dill_cvu2l(s, offset, offset);       /* op_i_cvi2l */
			break;
		    case DILL_L:
			dill_mulli(s, offset, opArg, size);  /* op_i_mulli */
			break;
		    case DILL_UL:
			dill_mululi(s, offset, opArg, size); /* op_i_mululi */
			break;
		    }
		} else {
		    offset = opArg;
		}
		dill_subp(s, result, opPtr, offset);           /* op_i_subp */ 
	    }
	}
	break;
	default:
	    dill_Psub(s, op_type, result, left, right);
	    break;
	}
	ret_op.reg = result;
	return ret_op;
    case  op_leq:{
	dill_mark_label_type true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_blei(s, left, right, true_label);	/* op_i_blei */
	    break;
	case DILL_U:
	    dill_bleu(s, left, right, true_label);	/* op_i_bleu */
	    break;
	case DILL_L:
	    dill_blel(s, left, right, true_label);	/* op_i_blel */
	    break;
	case DILL_UL:
	    dill_bleul(s, left, right, true_label);	/* op_i_bleul */
	    break;
	case DILL_F:
	    dill_blef(s, left, right, true_label);	/* op_i_blef */
	    break;
	case DILL_D:
	    dill_bled(s, left, right, true_label);	/* op_i_bled */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_lt:{
	dill_mark_label_type true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_blti(s, left, right, true_label);	/* op_i_blti */
	    break;
	case DILL_U:
	    dill_bltu(s, left, right, true_label);	/* op_i_bltu */
	    break;
	case DILL_L:
	    dill_bltl(s, left, right, true_label);	/* op_i_bltl */
	    break;
	case DILL_UL:
	    dill_bltul(s, left, right, true_label);	/* op_i_bltul */
	    break;
	case DILL_F:
	    dill_bltf(s, left, right, true_label);	/* op_i_bltf */
	    break;
	case DILL_D:
	    dill_bltd(s, left, right, true_label);	/* op_i_bltd */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_geq:{
	dill_mark_label_type true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_bgei(s, left, right, true_label);	/* op_i_bgei */
	    break;
	case DILL_U:
	    dill_bgeu(s, left, right, true_label);	/* op_i_bgeu */
	    break;
	case DILL_L:
	    dill_bgel(s, left, right, true_label);	/* op_i_bgel */
	    break;
	case DILL_UL:
	    dill_bgeul(s, left, right, true_label);	/* op_i_bgeul */
	    break;
	case DILL_F:
	    dill_bgef(s, left, right, true_label);	/* op_i_bgef */
	    break;
	case DILL_D:
	    dill_bged(s, left, right, true_label);	/* op_i_bged */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_gt:{
	dill_mark_label_type true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_bgti(s, left, right, true_label);	/* op_i_bgti */
	    break;
	case DILL_U:
	    dill_bgtu(s, left, right, true_label);	/* op_i_bgtu */
	    break;
	case DILL_L:
	    dill_bgtl(s, left, right, true_label);	/* op_i_bgtl */
	    break;
	case DILL_UL:
	    dill_bgtul(s, left, right, true_label);	/* op_i_bgtul */
	    break;
	case DILL_F:
	    dill_bgtf(s, left, right, true_label);	/* op_i_bgtf */
	    break;
	case DILL_D:
	    dill_bgtd(s, left, right, true_label);	/* op_i_bgtd */
							/* op_i_bgtp */
							/* op_i_bnep */
							/* op_i_beqp */
							/* op_i_blep */
							/* op_i_bgep */
							/* op_i_bltp */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_eq:{
	dill_mark_label_type true_label;
	if (string_op) {
	    /* this has to be strcmp */
	    dill_reg rv;
	    dill_push_init(s);
	    if (dill_do_reverse_vararg_push(s)) {
		dill_push_argp(s, right);
		dill_push_argp(s, left);
	    } else {
		dill_push_argp(s, left);
		dill_push_argp(s, right);
	    }
#ifdef i_callii
	    rv = dill_getreg(s, DILL_I);
	    i_callii(rv, (void *)cod_streq);
#else
	    rv = dill_calli(s, (void *)cod_streq, "cod_streq");
#endif
	    dill_movi(s, result, rv);
	    break;
	}
	true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_beqi(s, left, right, true_label);	/* op_i_beqi */
	    break;
	case DILL_U:
	    dill_bequ(s, left, right, true_label);	/* op_i_bequ */
	    break;
	case DILL_L:
	    dill_beql(s, left, right, true_label);	/* op_i_beql */
	    break;
	case DILL_UL:
	    dill_bequl(s, left, right, true_label);	/* op_i_bequl */
	    break;
	case DILL_F:
	    dill_beqf(s, left, right, true_label);	/* op_i_beqf */
	    break;
	case DILL_D:
	    dill_beqd(s, left, right, true_label);	/* op_i_beqd */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_neq: {
	dill_mark_label_type true_label = dill_alloc_label(s, "true");
	dill_seti(s, result, 0);	/* op_i_seti */
	switch (op_type) {
	case DILL_I:
	    dill_bnei(s, left, right, true_label);	/* op_i_bnei */
	    break;
	case DILL_U:
	    dill_bneu(s, left, right, true_label);	/* op_i_bneu */
	    break;
	case DILL_L:
	    dill_bnel(s, left, right, true_label);	/* op_i_bnel */
	    break;
	case DILL_UL:
	    dill_bneul(s, left, right, true_label);	/* op_i_bneul */
	    break;
	case DILL_P:
	    dill_bnep(s, left, right, true_label);	/* op_i_bnep */
	    break;
	case DILL_F:
	    dill_bnef(s, left, right, true_label);	/* op_i_bnef */
	    break;
	case DILL_D:
	    dill_bned(s, left, right, true_label);	/* op_i_bned */
	    break;
	default:
	    assert(FALSE);
	}
	dill_seti(s, result, 1);	/* op_i_seti */
	dill_mark_label(s, true_label);
	break;
    }
    case  op_log_or:
    case  op_arith_or:
	switch (op_type) {
	case DILL_P:
	    if(op == op_arith_or) {
		/* Arithmetic OR not allowed on pointers. */
		assert(FALSE);
	    } else {
		dill_cvul2p(s, left, left);	/* op_i_cvul2p */
		dill_cvul2p(s, right, right);	/* op_i_cvul2p */
		dill_orul(s, result, left, right);
	    }
	    break;
	default:
	    dill_Por(s, op_type, result, left, right);
	}
	break;
    case op_not:
	left = dill_getreg(s, op_type);
	switch (op_type) {
	case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	case DILL_I: case DILL_U: case DILL_L: case DILL_UL:
	    dill_piset(s, op_type, left, ~0); break;
	}
	/* falling through */
    case op_arith_xor:
        dill_Pxor(s, op_type, result, left, right);
	break;
    case  op_log_and:
    case  op_arith_and:
	switch (op_type) {
	case DILL_P:
	    if(op == op_arith_and) {
		/* Arithmetic AND not allowed on pointers. */
		assert(FALSE);
	    } else {
		dill_cvul2p(s, left, left);	/* op_i_cvul2p */
		dill_cvul2p(s, right, right);	/* op_i_cvul2p */
		dill_andul(s, result, left, right);
	    }
	    break;
	default:
	    dill_Pand(s, op_type, result, left, right);
	}
	break;
    case op_right_shift:
        dill_Prsh(s, op_type, result, left, right);
	break;
    case op_left_shift:
        dill_Plsh(s, op_type, result, left, right);
	break;
    case  op_mult:
        dill_Pmul(s, op_type, result, left, right);
	break;
    case  op_div:
        dill_Pdiv(s, op_type, result, left, right);
	break;
    case  op_deref:
    {
	sm_ref typ = get_complex_type(NULL, sm_right);
	operand right_op;
	init_operand(&right_op);
	right_op.reg = right;
	right_op.is_addr = 1;
	right_op.in_kernel = 0;
	right_op.offset = 0;

	if(typ && typ->node_type == cod_reference_type_decl) {
	    if(typ->node.reference_type_decl.kernel_ref)
		right_op.in_kernel = 1;
	    else
		right_op.in_kernel = 0;
	}

	if (!need_assignable) {
	    if(right_op.in_kernel) {
		dill_special(s, DILL_SEGMENTED_FOLLOWS, DILL_X86_GS_PREFIX);
	    }
	    gen_load(s, result, right_op, result_type);
	    right_op.reg = result;
	    right_op.is_addr = 0;
	    right_op.in_kernel = 0;
	}
	return right_op;
    }
    case  op_address:
	ret_op.reg = right;
	return ret_op;
    case op_log_neg:
	/* handled by if at top */
	assert(FALSE);
	break;
    case op_sizeof:
    case op_inc:
    case op_dec: 
	/* handled at top */
	assert(FALSE);
	break;
    }
    ret_op.reg = result;
    return ret_op;
}

static int
next_formal_is_drisc_exec_ctx(sm_list formals)
{
    sm_ref formal, ct;
    if (formals == NULL) return 0;
    formal = formals->node;
    ct = formal->node.declaration.sm_complex_type;
    if (ct == NULL) return 0;
    if ((ct->node_type == cod_reference_type_decl) &&
	(ct->node.reference_type_decl.name != NULL) &&
	(strcmp(ct->node.reference_type_decl.name, "cod_exec_context") == 0)) {
	return 1;
    }
    return 0;
}

static void
cod_do_noop_statement(dill_stream s, sm_ref func_ref, dill_reg duration)
{
    dill_mark_label_type top_label = dill_alloc_label(s, "noop_loop");
    dill_reg count = dill_getreg(s, DILL_I);
    dill_seti(s, count, 0);
    dill_mark_label(s, top_label);
    dill_addii(s, count, count, 1);
    dill_blei(s, count, duration, top_label);
}

#ifdef NOTDEF
static int
next_formal_is_cod_type_spec(sm_list formals)
{
    sm_ref formal, ct;
    if (formals == NULL) return 0;
    formal = formals->node;
    ct = formal->node.declaration.sm_complex_type;
    if (ct == NULL) return 0;
    if ((ct->node_type == cod_reference_type_decl) &&
	(strcmp(ct->node.reference_type_decl.name, "cod_type_spec")
	 == 0)) {
	return 1;
    }
    return 0;
}
#endif    

#define MAX_ARG 128

static operand
cg_subroutine_call(dill_stream s, sm_ref expr, cod_code descr)
{
    sm_ref func_ref = expr->node.subroutine_call.sm_func_ref;
    sm_list formals = func_ref->node.declaration.params;
    sm_list args = expr->node.subroutine_call.arguments;
    sm_ref arg, formal;
    dill_reg args_array[MAX_ARG];  /* should be large enough */
    int types_array[MAX_ARG];  /* should be large enough */
    int start, i, direction;
    int arg_count = 0;

    operand ret;

    init_operand(&ret);
    ret.is_addr = 0;
    ret.offset = 0;
    if (next_formal_is_drisc_exec_ctx(formals)) {
	if (arg_count < MAX_ARG) {
	    args_array[arg_count] = 0; /* context must be first param */
	    types_array[arg_count] = DILL_P;
	}
	arg_count++;
	formals = formals->next;
    }
    while (args != NULL) {
	operand param;
	int formal_type;
	arg = args->node;
	if (formals != NULL) {
	    formal = formals->node;
	    if (strcmp(formal->node.declaration.id, "...") == 0) {
		formal = NULL;
	    }
	} else {
	    formal = NULL;
	}
	if (formal != NULL) {
	    int actual_type = cod_sm_get_type(arg);
	    formal_type = formal->node.declaration.cg_type;
	    if ((actual_type == DILL_B) || is_array(arg)) {
		/* structure param */
		param = cg_expr(s, arg, 1, descr);
		if (is_complex_type(arg) && (param.is_addr == 0)) {
		    param.is_addr = 1;
		}
		if (param.offset != 0) {
		    dill_reg tmp = dill_getreg(s, DILL_P);
		    dill_addpi(s, tmp, param.reg, param.offset);
		    param.reg = tmp;
		    param.offset = 0;
		}
		formal_type = DILL_P;
	    } else {
		param = cg_expr(s, arg, 0, descr);
		param.reg = coerce_type(s, param.reg, formal_type,
					actual_type);
	    }
	} else {
	    int actual_type = cod_sm_get_type(arg);
	    if ((actual_type == DILL_B) || is_array(arg)) {
		/* structure param */
		param = cg_expr(s, arg, 1, descr);
		if (is_complex_type(arg) && (param.is_addr == 0)) {
		    param.is_addr = 1;
		}
		if (param.offset != 0) {
		    dill_reg tmp = dill_getreg(s, DILL_P);
		    dill_addpi(s, tmp, param.reg, param.offset);
		    param.reg = tmp;
		    param.offset = 0;
		}
		formal_type = DILL_P;
	    } else {
		param = cg_expr(s, arg, 0, descr);
		formal_type = cod_sm_get_type(arg);
	    }
	}
	if (arg_count < MAX_ARG) {
	    args_array[arg_count] = param.reg;
	    types_array[arg_count] = formal_type;
	}
	arg_count++;
	args = args->next;
	if (formals) formals = formals->next;
	if (next_formal_is_drisc_exec_ctx(formals)) {
	    if (arg_count < MAX_ARG) {
		args_array[arg_count] = 0; /* context must be first param */
		types_array[arg_count] = DILL_P;
	    }
	    arg_count++;
	    formals = formals->next;
	}
    }	

    if (dill_do_reverse_vararg_push(s)) {
	start = arg_count - 1;
	direction = -1;
    } else {
	start = 0;
	direction = 1;
    }
    if (strcmp(func_ref->node.declaration.id, "cod_NoOp") == 0) {
	cod_do_noop_statement(s, func_ref, args_array[0]);
	return ret;
    }
    if (func_ref->node.declaration.varidiac_subroutine_param_count == -1) {
	dill_push_init(s);
    } else {
	dill_push_varidiac_init(s, func_ref->node.declaration.varidiac_subroutine_param_count);
    }
    for (i=0; i< arg_count; i++) {
	int index = start + i * direction;
	dill_push_arg(s, types_array[index], args_array[index]);
    }
    ret.reg = dill_pcall(s, cod_sm_get_type(expr),
			 func_ref->node.declaration.cg_address, 
			 strdup(func_ref->node.declaration.id));
    return ret;
}

static void
cod_expand_dyn_array(void *base_addr, long new_size, long old_size, long struct_size)
{

    static int debug_cg = -1;

    if (debug_cg == -1) {
	debug_cg = (int)(intptr_t)(getenv("COD_DEBUG"));
    }
    if (debug_cg) {
	printf("cod_expand_dyn_array, base_addr %p, old_base %p, new_size %ld, old_size %ld, struct_size %ld\n",
	       base_addr, *(void**)base_addr, new_size, old_size, 
	       struct_size);
    }
    if (*(void**)base_addr == NULL) {
	*(void**)base_addr = malloc(new_size * struct_size);
	memset(*(void**)base_addr, 0, new_size *struct_size);
    } else {
	if (new_size > old_size) {
	    long malloc_size = new_size * struct_size;
	    long memset_size = (new_size - old_size) * struct_size;
	    void *cur_base = *(char **)base_addr;
	    void *new_base = realloc(cur_base, malloc_size);
	    memset((char*)new_base + (old_size *struct_size), 0,  memset_size);
	    *(void**)base_addr = new_base;
	}
    }	
    if (debug_cg) {
	printf("\tnew base %p, new size %ld, end %p\n", *(void**)base_addr, new_size*struct_size, *(char**)base_addr + new_size*struct_size);
    }
}

static void do_var_array_adjust(dill_stream s, sm_ref expr, operand old_val, 
				operand new_val, cod_code descr)
{
    sm_ref strct;
    sm_ref struct_decl;
    sm_list fields;
    assert(expr->node_type == cod_field_ref);

    strct = expr->node.field_ref.struct_ref;
    expr = expr->node.field_ref.sm_field_ref;
    struct_decl = get_complex_type(NULL, strct);
    if (struct_decl->node_type == cod_reference_type_decl) {
	struct_decl = struct_decl->node.reference_type_decl.sm_complex_referenced_type;
    }
    fields = struct_decl->node.struct_type_decl.fields;
    while(fields != NULL) {
	sm_ref tmp, darray_field;
	int size_reg, old_size;
	operand array_base;
	int dimen = -1, i = 0;
	sm_ref typ = fields->node->node.field.sm_complex_type;
	if ((typ == NULL) || (typ->node_type != cod_reference_type_decl)) {
	    fields = fields->next;
	    continue;
	}
	typ = typ->node.reference_type_decl.sm_complex_referenced_type;
	int static_size_multiplier = 1;
	while (typ != NULL) {
	    if (typ->node_type != cod_array_type_decl) {
		typ = NULL;
		continue;
	    }
	    if (typ->node.array_type_decl.sm_dynamic_size == expr) {
		dimen = i;
	    }
	    if (typ->node.array_type_decl.sm_dynamic_size == NULL) {
		static_size_multiplier *= typ->node.array_type_decl.cg_static_size;
	    }
	    i++;
	    typ = typ->node.array_type_decl.sm_complex_element_type;
	}
	if (dimen == -1) {
	    fields = fields->next;
	    continue;
	}
	darray_field = fields->node;
	tmp = malloc(sizeof(*tmp));
	tmp->node_type = cod_field_ref;
	tmp->node.field_ref.struct_ref = strct;
	tmp->node.field_ref.lx_field = NULL;
	tmp->node.field_ref.sm_field_ref = darray_field;
	array_base = cg_expr(s, tmp, 1, descr);
	free(tmp);
	fields = fields->next;
	if (array_base.offset != 0) {
	    int new_base = dill_getreg(s, DILL_P);
	    dill_addpi(s, new_base, array_base.reg, array_base.offset); /* op_i_addpi */
	    array_base.reg = new_base;
	    array_base.offset = 0;
	}
	size_reg = dill_getreg(s, DILL_I);
	old_size = dill_getreg(s, cod_sm_get_type(expr));
	gen_load(s, old_size, old_val, DILL_I);
	dill_seti(s, size_reg, darray_field->node.field.cg_size * static_size_multiplier);

	dill_push_init(s);
	if (dill_do_reverse_vararg_push(s)) {
	    dill_push_argl(s, size_reg);
	    dill_push_argl(s, old_size);
	    dill_push_argl(s, new_val.reg);
	    dill_push_argp(s, array_base.reg);
	} else {	    
	    dill_push_argp(s, array_base.reg);
	    dill_push_argl(s, new_val.reg);
	    dill_push_argl(s, old_size);
	    dill_push_argl(s, size_reg);
	}
	dill_callv(s, (void*)&cod_expand_dyn_array, "cod_expand_dyn_array");  /* op_i_callvi */
    }
}

static int 
is_var_array(sm_ref expr)
{
    sm_ref typ;
    if (expr->node_type == cod_field_ref) {
	return is_var_array(expr->node.field_ref.sm_field_ref);
    }
    if (expr->node_type == cod_field) {
	sm_ref ref = expr->node.field.sm_complex_type;
	if (ref == NULL) return 0;
	if (ref->node_type != cod_reference_type_decl) return 0;
	if (ref->node.reference_type_decl.sm_complex_referenced_type == NULL)
	    return 0;
	return (ref->node.reference_type_decl.sm_complex_referenced_type->node_type == cod_array_type_decl);
    }
    typ = get_complex_type(NULL, expr);
    if ((typ == NULL) || (typ->node_type != cod_array_type_decl)) return 0;
    if (typ->node.array_type_decl.cg_static_size == -1) {
	/* this array is varying */
	return 1;
    }
    if (typ->node.array_type_decl.cg_element_size == -1) {
	/* subarray is varying */
	return 1;
    }

    return 0;
}

static int 
is_static_var(sm_ref expr)
{
    if (expr->node_type == cod_identifier) {
	return is_static_var(expr->node.identifier.sm_declaration);
    } else if (expr->node_type == cod_declaration) {
	return expr->node.declaration.static_var;
    } else {
	assert(0);
    }
    return 0;
}

static dill_reg
load_dynamic_array_dimension(dill_stream s, sm_ref containing, sm_ref field, cod_code descr)
{
    sm_ref tmp;
    operand dsize;
    tmp = malloc(sizeof(*tmp));
    tmp->node_type = cod_field_ref;
    tmp->node.field_ref.struct_ref = containing;
    tmp->node.field_ref.lx_field = NULL;
    tmp->node.field_ref.sm_field_ref = field;
    dsize = cg_expr(s, tmp, 0, descr);
    free(tmp);
    return dsize.reg;
}

static int  /* return true if this is a structured type or a reference type */
is_complex_type(sm_ref expr)
{
    return (get_complex_type(NULL, expr) != NULL);
}

extern int type_of_int_const_string(char *val);

static operand
cg_expr(dill_stream s, sm_ref expr, int need_assignable, cod_code descr)
{
    operand oprnd;

    init_operand(&oprnd);
    switch(expr->node_type) {
    case cod_identifier:
	return cg_expr(s, expr->node.identifier.sm_declaration, need_assignable, descr);
    case cod_declaration:
    {
	sm_ref typ = expr->node.declaration.sm_complex_type;
	if (typ && (typ->node_type == cod_struct_type_decl) &&
	    (typ->node.struct_type_decl.encode_info != NULL)) {
	    /* this is a PBIO-encoded structure */
	    oprnd.enc.is_encoded = 1;
	    oprnd.enc.byte_swap_on_fetch = 0;
	    if (dill_target_byte_order(s) != typ->node.struct_type_decl.encode_info->byte_order)
		oprnd.enc.byte_swap_on_fetch = 1;
	}
	if (expr->node.declaration.const_var && expr->node.declaration.init_value && !typ) {
	    oprnd = cg_expr(s, expr->node.declaration.init_value, 0, descr);
	    return oprnd;
	}
	if (expr->node.declaration.static_var || 
	    expr->node.declaration.is_extern) {
	    dill_reg lvar = dill_getreg(s, DILL_P);
	    dill_reg static_block;
	    if (expr->node.declaration.static_var) {
		static_block = cg_get_static_block(s, descr);
		dill_addpi(s, lvar, static_block, (intptr_t)expr->node.declaration.cg_address);  /* op_i_addpi */
	    } else {
		dill_setp(s, lvar, expr->node.declaration.cg_address);  /* op_i_addpi */
	    }
	    assert(oprnd.enc.is_encoded == 0);
	    oprnd.reg = lvar;
	    oprnd.offset = 0;
	    oprnd.is_addr = 1;
	    if (need_assignable == 1) {
		return oprnd;
	    } else if (expr->node.declaration.cg_type == DILL_B) {
		oprnd.is_addr = 0;
		return oprnd;
	    } else {
		dill_reg ret = dill_getreg(s, expr->node.declaration.cg_type);
		gen_load(s, ret, oprnd, expr->node.declaration.cg_type);
		oprnd.reg = ret;
		oprnd.is_addr = 0;
		oprnd.offset = 0;
		return oprnd;
	    }
	}
	if (need_assignable) {
	    if (expr->node.declaration.param_num != -1) {
		/* it's a parameter */
		int is_reference_type = 0;
		if (expr->node.declaration.sm_complex_type &&
		    (expr->node.declaration.sm_complex_type->node_type ==
		    cod_reference_type_decl)) {
		    is_reference_type++;
		}
		if ((expr->node.declaration.cg_type != DILL_P) ||
		    (is_reference_type)) {
		    /* passed by value */
		    oprnd.reg = expr->node.declaration.cg_oprnd;
		    oprnd.is_addr = 0;
		    oprnd.offset = 0;
		} else {
		    /* passed by reference */
		    oprnd.reg = expr->node.declaration.cg_oprnd;
		    oprnd.is_addr = 1;
		    oprnd.offset = 0;
		    oprnd.enc.string_base = oprnd.reg;
		}
	    } else {
		/* it's a local */
		if ((expr->node.declaration.cg_type != DILL_B)  &&
		    ! expr->node.declaration.addr_taken) {
		    /* value is assignable */
		    oprnd.reg = expr->node.declaration.cg_oprnd;
		    oprnd.is_addr = 0;
		    oprnd.offset = 0;
		} else {
		    oprnd.reg = expr->node.declaration.cg_oprnd;
		    oprnd.is_addr = 1;
		    oprnd.offset = 0;
		}
	    }
	} else {
	    if (expr->node.declaration.addr_taken) {
		dill_reg ret = dill_getreg(s, expr->node.declaration.cg_type);
		oprnd.reg = expr->node.declaration.cg_oprnd;
		oprnd.is_addr = 1;
		oprnd.offset = 0;
		oprnd.enc.string_base = oprnd.reg;
		gen_load(s, ret, oprnd, expr->node.declaration.cg_type);
		oprnd.reg = ret;
		oprnd.is_addr = 0;
		oprnd.offset = 0;
		oprnd.enc.is_encoded = 0;
	    } else {
		oprnd.reg = expr->node.declaration.cg_oprnd;
		oprnd.is_addr = 0;
		oprnd.offset = 0;
		oprnd.enc.string_base = oprnd.reg;
	    }
	}
	return oprnd;
    }
    case cod_comma_expression:{
	(void) cg_expr(s, expr->node.comma_expression.left, 0, descr);
	return cg_expr(s, expr->node.comma_expression.right, need_assignable, descr);
    }
    case cod_operator: {
	oprnd = cg_operator(s, expr, need_assignable, descr);
	return oprnd;
    }
    case cod_cast: {
	operand right;
	int expr_type = cod_sm_get_type(expr->node.cast.expression);
	right = cg_expr(s, expr->node.cast.expression, 0, descr);
	if (expr->node.cast.cg_type == DILL_V) return right;  /* likely ignored */
	if (right.is_addr != 0) {
	    fprintf(stderr, "Failure in processing of cast, expr is :\n");
	    cod_print(expr->node.cast.expression);
	    fprintf(stderr, "Please report this bug and include as much information as possible for reproduction\n");
	    assert(right.is_addr == 0);
	}
	right.reg = coerce_type(s, right.reg, expr->node.cast.cg_type, 
				expr_type);
	return right;
    }
    case cod_constant: {
	dill_reg lvar = -1;
	if (expr->node.constant.token == string_constant) {
	    lvar = dill_getreg(s, DILL_P);
	    /* memory leak of constant value here.... */
	    dill_setp(s, lvar, strdup(expr->node.constant.const_val));  /* op_i_setp */
	} else if (expr->node.constant.token == floating_constant) {
	    double d;

	    lvar = dill_getreg(s, DILL_D);
	    sscanf(expr->node.constant.const_val, "%lf", &d);
	    dill_setd(s, lvar, d);	/* op_i_setd */
	} else if (expr->node.constant.token == character_constant) {
	    long l;
	    char *val = expr->node.constant.const_val;
	    lvar = dill_getreg(s, DILL_I);
	    if (*val == 'L') val++ ;
	    if (*val == 'u') val++ ;
	    if (*val == 'U') val++ ;
	    val++ ; /* skip ' */
	    if (*val != '\\') {
		l = *val;
	    } else {
		val++;
		switch(*val) {
		case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
		    if (sscanf(val, "%lo", &l) != 1) 
			printf("octal char sscanf failed %s\n", val);
		    break;
		case 'x':
		    if (sscanf(val+1, "%lx", &l) != 1) 
			printf("hex char sscanf failed, %s\n", val);
		    break;
		case '\\':
		    l = 0134;
		    break;
		case '\'':
		    l = 047;
		    break;
		case '"':
		    l = 042;
		    break;
		case '?':
		    l = 077;
		    break;
		case 'a':
		    l = 07;
		    break;
		case 'b':
		    l = 010;
		    break;
		case 'f':
		    l = 06;
		    break;
		case 'n':
		    l = 012;
		    break;
		case 'r':
		    l = 015;
		    break;
		case 't':
		    l = 011;
		    break;
		default:
		    l = 0;
		    printf("Bad character constant %s\n", val);
		}
	    }
	    dill_seti(s, lvar, l);	/* op_i_setc */
	} else {   /* integer_constant */
	    long i;
	    char *val = expr->node.constant.const_val;
	    int typ = type_of_int_const_string(val);
	    lvar = dill_getreg(s, typ);
	    if (val[0] == '0') {
		/* hex or octal */
		if (val[1] == 'x') {
		    /* hex */
		    if (sscanf(val+2, "%lx", &i) != 1) 
			printf("hex sscanf failed, %s\n", val);
		} else if (val[1] == 'b') {
		    /* binary */
		    int j = 2;
		    i = 0;
		    while (val[j]) {
			i <<= 1;
			if (val[j] == '1') {
			    i += 1;
			}
			j++;
		    }
		} else {
		    if (sscanf(val, "%lo", &i) != 1) 
			printf("octal sscanf failed %s\n", val);
		}
	    } else {
		if (sscanf(val, "%ld", &i) != 1) 
		    printf("decimal sscanf failed %s\n", val);
	    }
	    switch (typ) {
	    case DILL_I:
		dill_seti(s, lvar, i);	/* op_i_setl */
		break;
	    case DILL_U:
		dill_setu(s, lvar, i);	/* op_i_setl */
		break;
	    default:
	    case DILL_L:
		dill_setl(s, lvar, i);	/* op_i_setl */
		break;
	    case DILL_UL:
		dill_setu(s, lvar, i);	/* op_i_setl */
		break;
	    }
	}
	oprnd.reg = lvar;
	oprnd.is_addr = 0;
	oprnd.offset = 0;
	return oprnd;
    }
    case cod_enumerator:{
	dill_reg lvar = -1;
	lvar = dill_getreg(s, DILL_I);
	dill_seti(s, lvar, expr->node.enumerator.enum_value);	/* op_i_seti */
	oprnd.reg = lvar;
	oprnd.is_addr = 0;
	oprnd.offset = 0;
	return oprnd;
    }
    case cod_assignment_expression: {
	operand right;
	operand left;
	int assign_type = cod_sm_get_type(expr);
	right = cg_expr(s, expr->node.assignment_expression.right, 0, descr);
	left = cg_expr(s, expr->node.assignment_expression.left, 1, descr);
	if (right.is_addr != 0) {
	    fprintf(stderr, "Failure in processing RHS of assignment, expr is :\n");
	    cod_print(expr->node.assignment_expression.right);
	    fprintf(stderr, "Please report this bug and include as much information as possible for reproduction\n");
	    assert(right.is_addr == 0);
	}
	right.reg = coerce_type(s, right.reg, assign_type, 
				cod_sm_get_type(expr->node.assignment_expression.right));
	if (cod_expr_is_string(expr->node.assignment_expression.right)) {
	    dill_push_init(s);
	    dill_push_argp(s, right.reg);
#ifdef i_callii
	    right.reg = dill_getreg(s, DILL_P);
	    i_callli(rv, (void *)cod_strdup);
#else
	    right.reg = dill_calll(s, (void *)cod_strdup, "cod_strdup");
#endif
	}

	if(left.is_addr && left.in_kernel) {
	    dill_special(s, DILL_SEGMENTED_FOLLOWS, DILL_X86_GS_PREFIX);
	}
	
	if (expr->node.assignment_expression.op == op_eq) {
	    sm_ref left_val = expr->node.assignment_expression.left;
	    if (is_control_value(left_val, NULL)) {
		do_var_array_adjust(s, expr->node.assignment_expression.left, left,
				    right, descr);
	    }
	    if ((assign_type == DILL_B) || is_array(left_val)) {
		if (!right.is_addr) {
		    right.is_addr = 1;
		}
		(void) dill_scallv(s, (void*)memcpy, "memcpy", "%p%p%I", left.reg, right.reg, cg_get_size(s, left_val));
	    } else {
		gen_mov(s, left, right.reg, assign_type);
	    }
	} else {
	    dill_reg left_value;
	    operand result;
	    /* we have an operator-assignment */
	    if (left.is_addr == 1) {
		left_value = dill_getreg(s, assign_type);
		gen_load(s, left_value, left, assign_type);
	    } else {
		left_value = left.reg;
	    }
	    init_operand(&result);
	    result.reg = dill_getreg(s, assign_type);
	    result = execute_operator_cg(s, expr->node.assignment_expression.op, assign_type, result.reg, left_value, right.reg, 
					     expr->node.assignment_expression.left,
					     expr->node.assignment_expression.right, 0, assign_type);
	    if (is_control_value(expr->node.assignment_expression.left, NULL)) {
		do_var_array_adjust(s, expr->node.assignment_expression.left, left,
				    result, descr);
	    }
	    gen_mov(s, left, result.reg, assign_type);
	}
	if ((need_assignable == 1) || (!left.is_addr)) {
	    return left;
	} else {
	    dill_reg ret = dill_getreg(s, assign_type);
	    gen_load(s, ret, left, assign_type);
	    oprnd.reg = ret;
	    oprnd.is_addr = 0;
	    oprnd.offset = 0;
	    return oprnd;
	}
    }
    case cod_field_ref: {
	operand base = cg_expr(s, expr->node.field_ref.struct_ref, 1, descr);
	sm_ref field = expr->node.field_ref.sm_field_ref;
	if (base.is_addr && get_complex_type(NULL, expr->node.field_ref.struct_ref)->node_type == cod_reference_type_decl) {
	    dill_reg ret = dill_getreg(s, DILL_P);
	    gen_load(s, ret, base, DILL_P);
	    base.reg = ret;
	    base.is_addr = 1;
	}
	base.offset += field->node.field.cg_offset;
	if (is_var_array(field) && (need_assignable == 0)) {
	    /* variable array */
	    dill_reg ret = dill_getreg(s, DILL_P);
	    if (base.enc.is_encoded) {
		gen_encoded_field_load(s, ret, base, DILL_P, expr);
	    } else {
		gen_load(s, ret, base, DILL_P);
	    }
	    if (base.enc.is_encoded) {
		dill_addp(s, ret, ret, base.enc.string_base);
	    }
	    oprnd.reg = ret;
	    oprnd.is_addr = 1;
	    oprnd.offset = 0;
	    return oprnd;
	}
	{
	    sm_ref base_expr = expr->node.field_ref.struct_ref;
	    if (is_complex_type(base_expr) && (base.is_addr == 0)) {
		base.is_addr = 1;
	    }
	}

	if (need_assignable == 1) {
	    return base;
	} else {
	    dill_reg ret = dill_getreg(s, field->node.field.cg_type);
	    if (base.enc.is_encoded) {
		gen_encoded_field_load(s, ret, base, field->node.field.cg_type, expr);
	    } else {
		gen_load(s, ret, base, field->node.field.cg_type);
	    }
	    if (base.enc.is_encoded) {
		if (field->node.field.cg_type == DILL_P) {
		    dill_addp(s, ret, ret, base.enc.string_base);
		}
	    }
	    oprnd.reg = ret;
	    oprnd.is_addr = 0;
	    oprnd.offset = 0;
	    return oprnd;
	}
    }
    case cod_element_ref: {
	operand base;
	operand index[100];
	int array_dimensions = 0;
	int cg_type;
	sm_ref arr;
	int element_size;
	sm_ref tmp = expr;
	int i;
	sm_ref containing = expr->node.element_ref.sm_containing_structure_ref;

	/* rightmost index is encountered first and is index[0] */
	while(tmp->node_type == cod_element_ref) {
	    index[array_dimensions] = cg_expr(s, tmp->node.element_ref.expression, 0, descr);
	    tmp = tmp->node.element_ref.array_ref;
	    array_dimensions++;
	}
	/* reverse index operands so that the leftmost is index[0] */
	for (i=0; i< array_dimensions/2; i++) {
	    operand tmp = index[i];
	    index[i] = index[array_dimensions - i - 1];
	    index[array_dimensions - i - 1] = tmp;
	}
	arr = get_complex_type(NULL, tmp);
	
	base = cg_expr(s, tmp, 1, descr);

	if (arr && (arr->node_type == cod_reference_type_decl)
	    && (is_array(arr) || is_static_var(tmp))) {
	    /* we didn't load the address */
	    int load_type = DILL_P;
	    dill_reg ret = dill_getreg(s, DILL_P);
	    if (base.enc.is_encoded) {
		gen_encoded_field_load(s, ret, base, DILL_P, expr->node.element_ref.array_ref);
	    } else {
		gen_load(s, ret, base, load_type);
	    }
	    if (base.enc.is_encoded) {
		dill_addp(s, ret, ret, base.enc.string_base);
	    }
	    base.reg = ret;
	    base.is_addr = 1;
	    base.offset = 0;
	    arr = arr->node.reference_type_decl.sm_complex_referenced_type;
	}
	/* base is valid at this point */
	cg_type = expr->node.element_ref.cg_element_type;
	if (cg_type != DILL_B) {
	    element_size  = dill_type_size(s, cg_type);
	} else {
	    element_size = cg_get_size(s, expr->node.element_ref.sm_complex_element_type);
	}
/*
ROW MAJOR  (C-style):

2D = index1 * dim2 * element_size + index2 * element_size
3D = index1 * dim2 * dim3 * element_size + index2 * dim3 * element_size + index3*element_size
4D = index1 * dim2 * dim3 * dim4 * element_size + index2 * dim3 * dim4 * element_size + index3 * dim4 * element_size + index4*element_size
*/
	if (arr && (arr->node_type == cod_array_type_decl)) {
	    dill_reg size_accum = dill_getreg(s, DILL_L);
	    dimen_p d = arr->node.array_type_decl.dimensions;
	    dill_reg new_base = dill_getreg(s, DILL_P);
	    for (i = array_dimensions - 1; i >= 0; i--) {
		dill_reg tmp;
		if (i == array_dimensions -1) {
		    dill_seti(s, size_accum, element_size);
		    dill_movp(s, new_base, base.reg);
		    base.reg = new_base;
		} else {
		    if (d->dimens[i+1].static_size != -1) {
			dill_mulii(s, size_accum, size_accum, d->dimens[i+1].static_size);
		    } else {
			dill_reg dsize = load_dynamic_array_dimension(s, containing, d->dimens[i+1].control_field, descr);
			dill_muli(s, size_accum, size_accum, dsize);
		    }
		}
		tmp = dill_getreg(s, DILL_L);
		dill_muli(s, tmp, index[i].reg, size_accum);
		dill_cvi2l(s, tmp, tmp);       /* op_i_cvi2l */
		dill_addp(s, base.reg, base.reg, tmp);
	    }
	} else {
	    /* must be non-standard array */
	    dill_reg size_accum = dill_getreg(s, DILL_L);
	    dill_reg new_base = dill_getreg(s, DILL_P);
	    dill_seti(s, size_accum, element_size);
	    dill_movp(s, new_base, base.reg);
	    base.reg = new_base;
	    dill_muli(s, size_accum, index[i].reg, size_accum);
	    dill_addp(s, base.reg, base.reg, size_accum);
	    base.is_addr = 1;
	}
	if (need_assignable == 1) {
	    return base;
	} else {
	    dill_reg ret = dill_getreg(s, cg_type);
	    gen_load(s, ret, base, cg_type);
	    oprnd.reg = ret;
	    oprnd.is_addr = 0;
	    oprnd.offset = 0;
	    return oprnd;
	}
    }
    case cod_subroutine_call:
	return cg_subroutine_call(s, expr, descr);
    case cod_conditional_operator: {
	operand ret_op;
	dill_reg result = dill_getreg(s, expr->node.conditional_operator.result_type);
	int result_type = expr->node.conditional_operator.result_type;
	operand e1_op, e2_op;

	init_operand(&ret_op);
	ret_op.is_addr = 0;
	ret_op.offset = 0;

	dill_mark_label_type false_label = dill_alloc_label(s, "conditional e2");
	dill_mark_label_type end_label = dill_alloc_label(s, "end conditional");
	cg_branch_if_false(s, expr->node.conditional_operator.condition, false_label, descr, 0);
	e1_op = cg_expr(s, expr->node.conditional_operator.e1, need_assignable, descr);
	dill_pmov(s, result_type, result, e1_op.reg);	/* op_i_movp */
	dill_jv(s, end_label);
	dill_mark_label(s, false_label);
	e2_op = cg_expr(s, expr->node.conditional_operator.e2, need_assignable, descr);
	e2_op.reg = coerce_type(s, e2_op.reg, result_type, 
				cod_sm_get_type(expr->node.conditional_operator.e2));
	dill_pmov(s, result_type, result, e2_op.reg);	/* op_i_movp */
	dill_mark_label(s, end_label);
	ret_op.reg = result;
	return ret_op;
    }
    default:
	fprintf(stderr, "Unknown case in cg_expression\n");
	cod_print(expr);
    }
    oprnd.reg = -1;
    oprnd.is_addr = -1;
    oprnd.offset = -1;
    return oprnd;
}

static dill_reg
coerce_type(dill_stream s, dill_reg obj, int target_type, int obj_type)
{
    dill_reg ret;

    if (target_type == obj_type) return obj;

    ret = dill_getreg(s, target_type);
    switch(target_type) {
    case DILL_C:
    case DILL_UC:
	switch (obj_type) {
	case DILL_UC:
	    dill_orii(s, ret, obj, 0xff);
	    break;
	case DILL_C:
	    dill_cvuc2i(s, ret, obj);
	    break;
	case DILL_I:
	    dill_cvi2c(s, ret, obj);	/* op_i_cvi2c */
	    break;
	case DILL_U:
	    dill_cvu2c(s, ret, obj);	/* op_i_cvu2c */
	    break;
	case DILL_L:
	    dill_cvl2c(s, ret, obj);	/* op_i_cvl2c */
	    break;
	case DILL_UL:
	    dill_cvul2c(s, ret, obj);	/* op_i_cvul2c */
	    break;
	case DILL_F:
	    if (target_type == DILL_UC) {
	        dill_cvf2u(s, ret, obj);	/* op_i_cvf2u */
		dill_cvu2c(s, ret, ret);	/* op_i_cvu2c */
	    } else {
	        dill_cvf2i(s, ret, obj);	/* op_i_cvf2i */
		dill_cvi2c(s, ret, ret);	/* op_i_cvi2c */
	    }	    
	    break;
	case DILL_D:
	    if (target_type == DILL_UC) {
	        dill_cvd2u(s, ret, obj);	/* op_i_cvd2u */
		dill_cvu2c(s, ret, ret);	/* op_i_cvu2c */
	    } else {
	        dill_cvd2i(s, ret, obj);	/* op_i_cvd2i */
		dill_cvi2c(s, ret, ret);	/* op_i_cvi2c */
	    }
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    dill_cvul2c(s, ret, ret);	/* op_i_cvul2c */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_S:
    case DILL_US:
	switch (obj_type) {
	case DILL_C:
	    if (target_type == DILL_S) {
		dill_cvc2i(s, ret, obj);	/* op_i_cvi2s */
		dill_cvi2s(s, ret, ret);
	    } else {
		dill_cvc2i(s, ret, obj);	/* op_i_cvi2us */
		dill_cvi2us(s, ret, ret);
	    }
	    break;
	case DILL_S:
	    if (target_type == DILL_S) {
		dill_movs(s, ret, obj);	/* op_i_cvi2s */
	    } else {
		dill_cvs2i(s, ret, obj);	/* op_i_cvi2us */
		dill_cvi2us(s, ret, ret);
	    }
	    break;
	case DILL_I:
	    if (target_type == DILL_S) {
		dill_cvi2s(s, ret, obj);	/* op_i_cvi2s */
	    } else {
		dill_cvi2us(s, ret, obj);	/* op_i_cvi2us */
	    }
	    break;
	case DILL_U:
	    if (target_type == DILL_S) {
		dill_cvu2s(s, ret, obj);	/* op_i_cvu2s */
	    } else {
		dill_cvu2us(s, ret, obj);	/* op_i_cvu2us */
	    }
	    break;
	case DILL_L:
	    if (target_type == DILL_S) {
		dill_cvl2s(s, ret, obj);	/* op_i_cvl2s */
	    } else {
		dill_cvl2us(s, ret, obj);	/* op_i_cvl2us */
	    }
	    break;
	case DILL_UL:
	    if (target_type == DILL_S) {
		dill_cvul2s(s, ret, obj);	/* op_i_cvul2s */
	    } else {
		dill_cvul2us(s, ret, obj);	/* op_i_cvul2us */
	    }
	    break;
	case DILL_F:
	    if (target_type == DILL_US) {
	        dill_cvf2u(s, ret, obj);	/* op_i_cvf2u */
		dill_cvu2s(s, ret, ret);	/* op_i_cvu2s */
	    } else {
	        dill_cvf2i(s, ret, obj);	/* op_i_cvf2i */
		dill_cvi2s(s, ret, ret);	/* op_i_cvi2s */
	    }
	    break;
	case DILL_D:
	    if (target_type == DILL_US) {
	        dill_cvd2u(s, ret, obj);	/* op_i_cvd2u */
		dill_cvu2s(s, ret, ret);	/* op_i_cvu2s */
	    } else {
	        dill_cvd2i(s, ret, obj);	/* op_i_cvd2i */
		dill_cvi2s(s, ret, ret);	/* op_i_cvi2s */
	    }
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    dill_cvul2s(s, ret, ret);	/* op_i_cvul2s */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_I:
	switch (obj_type) {
	case DILL_C:
	    dill_cvc2i(s, ret, obj);	/* op_i_cvc2i */
	    break;
	case DILL_UC:
	    dill_andii(s, ret, obj, 0xff); /* op_i_andii */
	    break;
	case DILL_S:
	    dill_cvs2i(s, ret, obj);	/* op_i_cvs2i */
	    break;
	case DILL_US:
	    dill_cvus2i(s, ret, obj);	/* op_i_cvus2i */
	    break;
	case DILL_I:
	    assert(FALSE);
	    break;
	case DILL_U:
	    dill_cvu2i(s, ret, obj);	/* op_i_cvu2i */
	    break;
	case DILL_L:
	    dill_cvl2i(s, ret, obj);	/* op_i_cvl2i */
	    break;
	case DILL_UL:
	    dill_cvul2i(s, ret, obj);	/* op_i_cvul2i */
	    break;
	case DILL_F:
	    dill_cvf2i(s, ret, obj);	/* op_i_cvf2i */
	    break;
	case DILL_D:
	    dill_cvd2i(s, ret, obj);	/* op_i_cvd2i */
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    dill_cvul2i(s, ret, ret);	/* op_i_cvul2i */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_U:
	switch (obj_type) {
	case DILL_C:
	    dill_cvc2u(s, ret, obj);	/* op_i_cvc2u */
	    break;
	case DILL_UC:
	    dill_andii(s, ret, obj, 0xff); /* op_i_andii */
	    break;
	case DILL_S:
	    dill_cvs2u(s, ret, obj);	/* op_i_cvs2u */
	    break;
	case DILL_US:
	    dill_cvus2u(s, ret, obj);	/* op_i_cvus2u */
	    break;
	case DILL_I:
	    dill_cvi2u(s, ret, obj);	/* op_i_cvi2u */
	    break;
	case DILL_U:
	    assert(FALSE);
	    break;
	case DILL_L:
	    dill_cvl2u(s, ret, obj);	/* op_i_cvl2u */
	    break;
	case DILL_UL:
	    dill_cvul2u(s, ret, obj);	/* op_i_cvul2u */
	    break;
	case DILL_F:
	    dill_cvf2u(s, ret, obj);	/* op_i_cvf2u */
	    break;
	case DILL_D:
	    dill_cvd2u(s, ret, obj);	/* op_i_cvd2u */
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    dill_cvul2u(s, ret, ret);	/* op_i_cvul2u */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_L:
	switch (obj_type) {
	case DILL_C:
	    dill_cvc2l(s, ret, obj);	/* op_i_cvc2l */
	    break;
	case DILL_UC:
	    dill_andii(s, ret, obj, 0xff); /* op_i_andii */
	    break;
	case DILL_S:
	    dill_cvs2l(s, ret, obj);	/* op_i_cvs2l */
	    break;
	case DILL_US:
	    dill_cvus2l(s, ret, obj);	/* op_i_cvus2l */
	    break;
	case DILL_I:
	    dill_cvi2l(s, ret, obj);	/* op_i_cvi2l */
	    break;
	case DILL_U:
	    dill_cvu2l(s, ret, obj);	/* op_i_cvu2l */
	    break;
	case DILL_L:
	    assert(FALSE);
	    break;
	case DILL_UL:
	    dill_cvul2l(s, ret, obj);	/* op_i_cvul2l */
	    break;
	case DILL_F:
	    dill_cvf2l(s, ret, obj);	/* op_i_cvf2l */
	    break;
	case DILL_D:
	    dill_cvd2l(s, ret, obj);	/* op_i_cvd2l */
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    dill_cvul2l(s, ret, ret);	/* op_i_cvul2l */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_UL:
	switch (obj_type) {
	case DILL_C:
	    dill_cvc2ul(s, ret, obj);	/* op_i_cvc2ul */
	    break;
	case DILL_UC:
	    dill_andii(s, ret, obj, 0xff); /* op_i_andii */
	    break;
	case DILL_S:
	    dill_cvs2ul(s, ret, obj);	/* op_i_cvs2ul */
	    break;
	case DILL_US:
	    dill_cvus2ul(s, ret, obj);	/* op_i_cvus2ul */
	    break;
	case DILL_I:
	    dill_cvi2ul(s, ret, obj);	/* op_i_cvi2ul */
	    break;
	case DILL_U:
	    dill_cvu2ul(s, ret, obj);	/* op_i_cvu2ul */
	    break;
	case DILL_L:
	    dill_cvl2ul(s, ret, obj);	/* op_i_cvl2ul */
	    break;
	case DILL_UL:
	    assert(FALSE);
	    break;
	case DILL_F:
	    dill_cvf2ul(s, ret, obj);	/* op_i_cvf2ul */
	    break;
	case DILL_D:
	    dill_cvd2ul(s, ret, obj);	/* op_i_cvd2ul */
	    break;
	case DILL_P:
	    dill_cvp2ul(s, ret, obj);	/* op_i_cvp2ul */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_F:
	switch (obj_type) {
	case DILL_S:
	case DILL_US:
	case DILL_UC:
	case DILL_C:
	    dill_cvi2f(s, ret, coerce_type(s, obj, DILL_I, obj_type));	/* op_i_cvi2f */
	    break;
	case DILL_I:
	    dill_cvi2f(s, ret, obj);	/* op_i_cvi2f */
	    break;
	case DILL_U:
	    dill_cvu2f(s, ret, obj);	/* op_i_cvu2f */
	    break;
	case DILL_L:
	    dill_cvl2f(s, ret, obj);	/* op_i_cvl2f */
	    break;
	case DILL_UL:
	    dill_cvul2f(s, ret, obj);	/* op_i_cvul2f */
	    break;
	case DILL_F:
	    assert(FALSE);
	    break;
	case DILL_D:
	    dill_cvd2f(s, ret, obj);	/* op_i_cvd2f */
	    break;
	default:
	    assert(FALSE);
	}
	break;
    case DILL_D:
	switch (obj_type) {
	case DILL_S:
	case DILL_US:
	case DILL_UC:
	case DILL_C:
	    dill_cvi2d(s, ret, coerce_type(s, obj, DILL_I, obj_type));	/* op_i_cvi2d */
	    break;
	case DILL_I:
	    dill_cvi2d(s, ret, obj);	/* op_i_cvi2d */
	    break;
	case DILL_U:
	    dill_cvu2d(s, ret, obj);	/* op_i_cvu2d */
	    break;
	case DILL_L:
	    dill_cvl2d(s, ret, obj);	/* op_i_cvl2d */
	    break;
	case DILL_UL:
	    dill_cvul2d(s, ret, obj);	/* op_i_cvul2d */
	    break;
	case DILL_F:
	    dill_cvf2d(s, ret, obj);	/* op_i_cvf2d */
	    break;
	case DILL_D:
	    assert(FALSE);
	    break;
	default:
	    assert(FALSE);
	}
	break;

    case DILL_P:
	switch (obj_type) {
	case DILL_L:
	    dill_cvl2ul(s, ret, obj);	/* op_i_cvl2ul */
	    dill_cvul2p(s, ret, ret);	/* op_i_cvul2p */
	    break;
	case DILL_UL:
	    dill_cvul2p(s, ret, obj);	/* op_i_cvul2p */
	    break;
	case DILL_P:
	    assert(FALSE);
	    break;
	case DILL_I:
	    dill_cvi2ul(s, ret, obj);	/* op_i_cvul2p */
	    dill_cvul2p(s, ret, ret);	/* op_i_cvul2p */
	    break;
	default:
	    assert(FALSE);
	}
	break;

    default:
	assert(FALSE);
    }
    return ret;
}

		
static void
gen_encoded_field_load(dill_stream s, dill_reg ret, operand base, int load_type, sm_ref field_ref)
{
    sm_ref typ = NULL;
    sm_ref field = field_ref->node.field_ref.sm_field_ref;
    int field_size = field->node.field.cg_size;
    if (field_ref->node_type == cod_field_ref)
	typ = get_complex_type(NULL, field_ref->node.field_ref.struct_ref);
    if (is_var_array(field)) 
	field_size = typ->node.struct_type_decl.encode_info->pointer_size;
    if (base.enc.is_encoded && (field_size != dill_type_size(s, load_type))) {
	/* encoded and native type sizes don't match */
	/* do some PBIO-like things here */
	if (field_size > dill_type_size(s, load_type)) {
	    /* field type is larger than we have, get low bits */
	    switch(load_type) {
	    case DILL_L:
	    case DILL_P:
		if (typ->node.struct_type_decl.encode_info->byte_order
		    == 1 /* Format_Integer_bigendian */) {
		    /* use low-end bytes */
		    base.offset += field_size - dill_type_size(s, load_type);
		}
	    }
	} else {
	    /* field type is smaller than we have, probably pointer */
	    if (field_size == dill_type_size(s, DILL_L)) {
		load_type = DILL_L;
	    } else if (field_size == dill_type_size(s, DILL_I)) {
		load_type = DILL_I;
	    } else if (field_size == dill_type_size(s, DILL_S)) {
		load_type = DILL_S;
	    } else {
		printf("Bad type size\n");
	    }
	}
    }
    gen_load(s, ret, base, load_type);
}

static void
gen_load(dill_stream s, dill_reg left, operand right, int type)
{
    int bswap = 0;
    if (right.enc.is_encoded && right.enc.byte_swap_on_fetch) {
	bswap++;
    }
    if (!dill_has_ldbs(s) || !bswap) {
	/* no byte-swapping load or not byte-swapping, do a simple load */
	switch(type) {
	case DILL_C:
	    dill_ldci(s, left, right.reg, right.offset);	/* op_i_ldci */
	    break;
	case DILL_UC:
	    dill_lduci(s, left, right.reg, right.offset);	/* op_i_lduci */
	    break;
	case DILL_S:
	    dill_ldsi(s, left, right.reg, right.offset);	/* op_i_ldsi */
	    break;
	case DILL_US:
	    dill_ldusi(s, left, right.reg, right.offset);	/* op_i_ldusi */
	    break;
	case DILL_I:
	    dill_ldii(s, left, right.reg, right.offset);	/* op_i_ldii */
	    break;
	case DILL_U:
	    dill_ldui(s, left, right.reg, right.offset);	/* op_i_ldui */
	    break;
	case DILL_L:
	    dill_ldli(s, left, right.reg, right.offset);	/* op_i_ldli */
	    break;
	case DILL_UL:
	    dill_lduli(s, left, right.reg, right.offset);	/* op_i_lduli */
	    break;
	case DILL_F:
	    dill_ldfi(s, left, right.reg, right.offset);	/* op_i_ldfi */
	    break;
	case DILL_D:
	    dill_lddi(s, left, right.reg, right.offset);	/* op_i_lddi */
	    break;
	case DILL_P:
	    dill_ldpi(s, left, right.reg, right.offset);  /* op_i_ldpi */
	    break;
	default:
	    fprintf(stderr, "unhandled case in gen_load\n");
	}
	if (bswap) {
	    switch (type) {
	    case DILL_C:
		break;
	    case DILL_UC:
		break;
	    case DILL_S:
		dill_bswaps(s, left, left);	/* op_i_bswaps */
		break;
	    case DILL_US:
		dill_bswapus(s, left, left);	/* op_i_bswapus */
		break;
	    case DILL_I:
		dill_bswapi(s, left, left);	/* op_i_bswapi */
		break;
	    case DILL_U:
		dill_bswapu(s, left, left);	/* op_i_bswapu */
		break;
	    case DILL_L:
	    case DILL_P:
		dill_bswapl(s, left, left);	/* op_i_bswapl */
		break;
	    case DILL_UL:
		dill_bswapul(s, left, left);	/* op_i_bswapul */
		break;
	    case DILL_F:
		dill_bswapf(s, left, left);	/* op_i_bswapf */
		break;
	    case DILL_D:
		dill_bswapd(s, left, left);	/* op_i_bswapd */
		break;
	    default:
		fprintf(stderr, "unhandled case in gen_load\n");
	    }
	}
    } else {
	/* must be byte-swapping with special load available */
	switch(type) {
	case DILL_C:
	    dill_ldbsci(s, left, right.reg, right.offset);
	    break;
	case DILL_UC:
	    dill_ldbsuci(s, left, right.reg, right.offset);
	    break;
	case DILL_S:
	    dill_ldbssi(s, left, right.reg, right.offset);
	    break;
	case DILL_US:
	    dill_ldbsusi(s, left, right.reg, right.offset);
	    break;
	case DILL_I:
	    dill_ldbsii(s, left, right.reg, right.offset);
	    break;
	case DILL_U:
	    dill_ldbsui(s, left, right.reg, right.offset);
	    break;
	case DILL_L:
	    dill_ldbsli(s, left, right.reg, right.offset);
	    break;
	case DILL_UL:
	    dill_ldbsuli(s, left, right.reg, right.offset);
	    break;
	case DILL_F:
	    dill_ldbsfi(s, left, right.reg, right.offset);
	    break;
	case DILL_D:
	    dill_ldbsdi(s, left, right.reg, right.offset);
	    break;
	case DILL_P:
	    dill_ldbspi(s, left, right.reg, right.offset);
	    break;
	default:
	    fprintf(stderr, "unhandled case in gen_load\n");
	}
    }
}

static void
gen_mov(dill_stream s, operand left, dill_reg right, int type)
{
    if (left.is_addr == 0) {
	switch(type) {
	case DILL_C:
	    dill_movc(s, left.reg, right);	/* op_i_movc */
	    break;
	case DILL_UC:
	    dill_movuc(s, left.reg, right);	/* op_i_movuc */
	    break;
	case DILL_S:
	    dill_movs(s, left.reg, right);	/* op_i_movs */
	    break;
	case DILL_US:
	    dill_movus(s, left.reg, right);	/* op_i_movus */
	    break;
	case DILL_I:
	    dill_movi(s, left.reg, right);	/* op_i_movi */
	    break;
	case DILL_U:
	    dill_movu(s, left.reg, right);	/* op_i_movu */
	    break;
	case DILL_L:
	    dill_movl(s, left.reg, right);	/* op_i_movl */
	    break;
	case DILL_UL:
	    dill_movul(s, left.reg, right);	/* op_i_movul */
	    break;
	case DILL_F:
	    dill_movf(s, left.reg, right);	/* op_i_movf */
	    break;
	case DILL_D:
	    dill_movd(s, left.reg, right);	/* op_i_movd */
	    break;
	case DILL_P:
	    dill_movp(s, left.reg, right);	/* op_i_movp */
	    break;
	default:
	    fprintf(stderr, "unhandled case in gen_mov, mov side\n");
	}
    } else {
	/* left is addr */
	switch(type) {
	case DILL_C:
	    dill_stci(s, right, left.reg, left.offset);	/* op_i_stci */
	    break;
	case DILL_UC:
	    dill_stuci(s, right, left.reg, left.offset);	/* op_i_stuci */
	    break;
	case DILL_S:
	    dill_stsi(s, right, left.reg, left.offset);	/* op_i_stsi */
	    break;
	case DILL_US:
	    dill_stusi(s, right, left.reg, left.offset);	/* op_i_stusi */
	    break;
	case DILL_I:
	    dill_stii(s, right, left.reg, left.offset);	/* op_i_stii */
	    break;
	case DILL_U:
	    dill_stui(s, right, left.reg, left.offset);	/* op_i_stui */
	    break;
	case DILL_L:
	    dill_stli(s, right, left.reg, left.offset);	/* op_i_stli */
	    break;
	case DILL_UL:
	    dill_stuli(s, right, left.reg, left.offset);	/* op_i_stuli */
	    break;
	case DILL_F:
	    dill_stfi(s, right, left.reg, left.offset);	/* op_i_stfi */
	    break;
	case DILL_D:
	    dill_stdi(s, right, left.reg, left.offset);	/* op_i_stdi */
	    break;
	case DILL_P:
	    dill_stpi(s, right, left.reg, left.offset);	/* op_i_stpi */
	    break;
	default:
	    fprintf(stderr, "unhandled case in gen_mov st side\n");
	}
    }	
}

static void
gen_bnz(dill_stream s, int conditional, int target_label, int op_type)
{
    switch(op_type) {
    case DILL_I:
	dill_bneii(s, conditional, 0, target_label);	/* op_i_bneii */
	break;
    case DILL_U:
	dill_bneui(s, conditional, 0, target_label);	/* op_i_bneui */
	break;
    case DILL_L:
	dill_bneli(s, conditional, 0, target_label);	/* op_i_bneli */
	break;
    case DILL_P:
	dill_bnepi(s, conditional, 0, target_label);	/* op_i_bnepi */
	break;
    case DILL_UL:
	dill_bneuli(s, conditional, 0, target_label);	/* op_i_bneuli */
	break;
    case DILL_F: {
	dill_reg tmp = dill_getreg(s, DILL_F);
	dill_setf(s, tmp, 0.0);	/* op_i_setf */
	dill_bnef(s, conditional, tmp, target_label);
	break;
    }
    case DILL_D: {
	dill_reg tmp = dill_getreg(s, DILL_D);
	dill_setd(s, tmp, 0.0);
	dill_bned(s, conditional, tmp, target_label);	/* op_i_bned */
	break;
    }
    default:
	fprintf(stderr, "unhandled case in gen_bnz\n");
    }
}

static void
gen_bz(dill_stream s, int conditional, int target_label, int op_type)
{
    switch(op_type) {
    case DILL_C: case DILL_S:
    case DILL_I:
	dill_beqii(s, conditional, 0, target_label);	/* op_i_beqii */
	break;
    case DILL_UC: case DILL_US:
    case DILL_U:
	dill_bequi(s, conditional, 0, target_label);	/* op_i_bequi */
	break;
    case DILL_L:
    case DILL_P:
	dill_beqli(s, conditional, 0, target_label);	/* op_i_beqli */
	break;
    case DILL_UL:
	dill_bequli(s, conditional, 0, target_label);	/* op_i_bequli */
	break;
    case DILL_F: {
	dill_reg tmp = dill_getreg(s, DILL_F);
	dill_setf(s, tmp, 0.0);	/* op_i_setf */
	dill_beqf(s, conditional, tmp, target_label);	/* op_i_beqf */
	break;
    }
    case DILL_D: {
	dill_reg tmp = dill_getreg(s, DILL_D);
	dill_setd(s, tmp, 0.0);	/* op_i_setd */
	dill_beqd(s, conditional, tmp, target_label);	/* op_i_beqd */
	break;
    }
    default:
	fprintf(stderr, "unhandled case in gen_bz op_type %d\n", op_type);
    }
}

static void cg_jump_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    if (stmt->node.jump_statement.goto_target) {
	/* this is a goto */
	sm_ref label_stmt = stmt->node.jump_statement.sm_target_stmt;
	dill_jv(s, label_stmt->node.label_statement.cg_label);
	return;
    } else {
	/* continue or break */
	sm_ref iterator = stmt->node.jump_statement.sm_target_stmt;
	if (stmt->node.jump_statement.continue_flag == 1) {
	    /* continue - jump to loop start (iterator evaluation) */
	    dill_jv(s, iterator->node.iteration_statement.cg_iter_label);
	} else {
	    /* break  -  jump to loop end */
	    dill_jv(s, iterator->node.iteration_statement.cg_end_label);
	}
    }
}

static void cg_return_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    operand ret_val;
    int expr_cg_type;
    int func_cg_type = stmt->node.return_statement.cg_func_type;
    if (func_cg_type == DILL_V) {
	dill_retii(s, 0);
	return;
    }
    expr_cg_type = cod_sm_get_type(stmt->node.return_statement.expression);
    ret_val = cg_expr(s, stmt->node.return_statement.expression, 0, descr);
    if (ret_val.is_addr) {
	/* returning a reference */
	ret_val.is_addr = 0;
    }

    ret_val.reg = coerce_type(s, ret_val.reg, func_cg_type, expr_cg_type);
#ifdef HAVE_DILL_H
    dill_pret(s, func_cg_type, ret_val.reg);
#else
    switch(func_cg_type) {
    case DILL_C:
    case DILL_UC:
    case DILL_S:
    case DILL_US:
    case DILL_I:
	dill_reti(s, ret_val.reg);	/* op_i_reti */
	break;
    case DILL_U:
	dill_retu(s, ret_val.reg);	/* op_i_retu */
	break;
    case DILL_L:
	dill_retl(s, ret_val.reg);	/* op_i_retl */
	break;
    case DILL_UL:
	dill_retul(s, ret_val.reg);	/* op_i_retul */
	break;
    case DILL_D:
	dill_retd(s, ret_val.reg);	/* op_i_retd */
	break;
    case DILL_F:
	dill_retf(s, ret_val.reg);	/* op_i_retf */
	break;
    default:
	fprintf(stderr, "unhandled case in return\n");
	break;
    }
#endif
}

static int
is_comparison_operator(sm_ref expr)
{
    if (expr->node_type != cod_operator) return 0;
    switch (expr->node.operator.op) {
    case op_neq: case op_eq: case op_gt: case op_geq: case op_lt: case op_leq:
	return 1;
    default:
	return 0;
    }
}

extern cod_exec_context
cod_create_exec_context(cod_code code)
{
    cod_exec_context cod_ctx = malloc(sizeof(struct cod_exec_struct));
    cod_ctx->ec = dill_get_exec_context(code->drisc_context);
    cod_ctx->gen_code = code;
    if (code->has_exec_ctx && code->static_size_required > 0) {
	void *block = malloc(code->static_size_required);
	/* copy initialized data */
	memcpy(block, code->data, code->static_size_required);
	dill_assoc_client_data(cod_ctx->ec, 0x23234, (intptr_t)block);
	cod_ctx->static_data = block;
    } else {
	cod_ctx->static_data = NULL;
    }
    return cod_ctx;
}


extern void
cod_exec_context_free(cod_exec_context ec)
{
    dill_free_exec_context(ec->ec);
    if (ec->static_data) {
	free(ec->static_data);
    }
    free(ec);
}

extern void *cod_extract_state(cod_exec_context ec, int *length_p)
{
    cod_code code = ec->gen_code;
    void *ret = malloc(code->static_size_required);
    /*    
    int i;
    printf("Static block size is %d\n", code->static_size_required);
    for (i=0; i < code->static_size_required; i++) {
	printf(" %02x", ((unsigned char*)ec->static_data)[i]);
    }
    printf("\n");*/
    memcpy(ret, ec->static_data, code->static_size_required);
    *length_p = code->static_size_required;
    /*    *ioformats = code->static_formats;*/
    
    return ret;
}

extern int cod_install_state(cod_exec_context ec, void *state, int length)
{
    memcpy(ec->static_data, state, length);
    return 1;
}

extern void
cod_assoc_client_data(cod_exec_context ec, int key, intptr_t value)
{
    dill_assoc_client_data(ec->ec, key, value);
}

extern intptr_t
cod_get_client_data(cod_exec_context ec, int key)
{
    return dill_get_client_data(ec->ec, key);
}

static int
reverse_operator(int cod_op)
{
    switch(cod_op) {
    case op_neq:   return dill_eq_code;
    case op_eq:   return dill_ne_code;
    case op_gt:   return dill_le_code;
    case op_geq:   return dill_lt_code;
    case op_lt:   return dill_ge_code;
    case op_leq:   return dill_gt_code;
    default:
	assert(0);
    }
    return -1;
}

static int
operator_to_dill_op(operator_t cod_op)
{
    switch(cod_op) {
    case op_neq:   return dill_ne_code;
    case op_eq:   return dill_eq_code;
    case op_gt:   return dill_gt_code;
    case op_geq:   return dill_gt_code;
    case op_lt:   return dill_lt_code;
    case op_leq:   return dill_le_code;
    default:
	assert(0);
    }
    return -1;
}

static void
cg_branch_if_false(dill_stream s, sm_ref pred, dill_mark_label_type label, 
		   cod_code descr, int reverse)
{
    operator_t op;
    int string_op;
    int branch_operator;
    dill_reg right = 0, left = 0;

    if (!is_comparison_operator(pred)) {
	operand conditional;
	conditional = cg_expr(s, pred, 0, descr);
	assert(conditional.is_addr == 0);
	if (reverse) {
	    gen_bnz(s, conditional.reg, label, cod_sm_get_type(pred));
	} else {
	    gen_bz(s, conditional.reg, label, cod_sm_get_type(pred));
	}
	return;
    }

    op = pred->node.operator.op;
    string_op = cod_expr_is_string(pred->node.operator.right);

    operator_prep(s, pred, &right, &left, descr);

    if (string_op) {
	/* this has to be strcmp */
	dill_reg rv;
	dill_push_init(s);
	if (dill_do_reverse_vararg_push(s)) {
	    dill_push_argp(s, right);
	    dill_push_argp(s, left);
	} else {
	    dill_push_argp(s, left);
	    dill_push_argp(s, right);
	}
	rv = dill_calli(s, (void *)cod_streq, "cod_streq");
	if (reverse) {
	    gen_bz(s, rv, label, DILL_I);
	} else {
	    gen_bnz(s, rv, label, DILL_I);
	}
	return;
    }

    if (reverse) {
	branch_operator = operator_to_dill_op(op);
    } else {
	branch_operator = reverse_operator(op);
    }
    dill_pbr(s, branch_operator, pred->node.operator.operation_type, 
	   left, right, label);
}
    
static void cg_selection_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    dill_mark_label_type else_label = dill_alloc_label(s, "else");
    
    cg_branch_if_false(s, stmt->node.selection_statement.conditional, 
		       else_label, descr, 0);
    cg_statement(s, stmt->node.selection_statement.then_part, descr);
    if (stmt->node.selection_statement.else_part != NULL) {
	dill_mark_label_type end_label = dill_alloc_label(s, "if-end");
	dill_jv(s, end_label);	/* op_i_jpi */
	dill_mark_label(s, else_label);
	cg_statement(s, stmt->node.selection_statement.else_part, descr);
	dill_mark_label(s, end_label);
    } else {
	dill_mark_label(s, else_label);
    }
}

static void cg_iteration_statement(dill_stream s, sm_ref stmt, cod_code descr)
{
    dill_mark_label_type begin_label = dill_alloc_label(s, "loop begin"), 
	end_label = dill_alloc_label(s, "loop end"),
	iter_label = dill_alloc_label(s, "loop iteration");
    stmt->node.iteration_statement.cg_end_label = end_label;
    stmt->node.iteration_statement.cg_iter_label = iter_label;
    if (stmt->node.iteration_statement.init_expr != NULL) {
	(void) cg_expr(s, stmt->node.iteration_statement.init_expr, 0, descr);
    }
    dill_mark_label(s, begin_label);
    if (stmt->node.iteration_statement.test_expr != NULL) {
	cg_branch_if_false(s, stmt->node.iteration_statement.test_expr, 
			   end_label, descr, 0);
    }
    cg_statement(s, stmt->node.iteration_statement.statement, descr);
    dill_mark_label(s, iter_label);
    if (stmt->node.iteration_statement.iter_expr != NULL) {
	(void) cg_expr(s, stmt->node.iteration_statement.iter_expr, 0, descr);
    }
    if (stmt->node.iteration_statement.post_test_expr == NULL) {
	dill_jv(s, begin_label);	/* op_i_jpi */
    } else {
	/* specify reverse decision, so this is cg_branch_if_*true* */
	cg_branch_if_false(s, stmt->node.iteration_statement.post_test_expr, 
			   begin_label, descr, 1);
    }
    dill_mark_label(s, end_label);
    
}
#endif
