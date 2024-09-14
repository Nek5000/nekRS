%{
#include "config.h"
#ifdef __NVCOMPILER
#pragma diag_suppress 550, 111, 941
#endif
#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 2215)
#endif
#ifdef SEGMENTED_POINTERS
int cod_segmented_pointers = 1;
#else
int cod_segmented_pointers = 0;
#endif
#ifdef KPLUGINS_INTEGRATION
int cod_kplugins_integration = 1;
#else
int cod_kplugins_integration = 0;
#endif
#ifndef LINUX_KERNEL_MODULE
#include "stdio.h"
#endif
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
#ifndef LINUX_KERNEL_MODULE
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#endif
#undef NDEBUG
#include "assert.h"
#ifndef LINUX_KERNEL_MODULE
#include <ctype.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#else
#include <linux/ctype.h>
#include <linux/string.h>
#include <linux/mm.h>
#endif
#include "float.h"
#ifdef DBL_DECIMAL_DIG
  #define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else  
  #ifdef DECIMAL_DIG
    #define OP_DBL_Digs (DECIMAL_DIG)
  #else  
    #define OP_DBL_Digs (DBL_DIG + 3)
  #endif
#endif
#include "fm.h"
#include "fm_internal.h"
#include "cod.h"
#include "cod_internal.h"
#include "structs.h"
#ifdef HAVE_DILL_H
#include "dill.h"
#else
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
    DILL_EC,
    DILL_ERR   /* no type */
};
typedef void *dill_stream;
#define dill_create_stream() 0
#define dill_type_size(c, s) 0
#endif
#if defined(_MSC_VER)
#define strdup _strdup
#define isatty _isatty
#define fileno _fileno
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#endif
#ifndef LINUX_KERNEL_MODULE
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#else
#include "kecl.h"
#define malloc (void *)DAllocMM
#define free(a) DFreeMM((addrs_t)a)
#define realloc(a,b) DReallocMM((addrs_t)a,b)
#define fprintf(fmt, args...) printk(args);
#define printf printk
char *strdup(const char *s)
{
	char *p;

	p = (char *)kmalloc(strlen(s)+1, GFP_KERNEL);
	if (p != NULL)
		strcpy(p,s);
	return p;
}
#endif
#ifdef _MSC_VER
#undef strncpy
#endif
#define YY_NO_INPUT

static char*
gen_anon()
{
    static int anon_count = 0;
    char *ret = malloc(40);
    snprintf(ret, 40, "Anonymous-%d", anon_count++);
    return ret;
}

#define yyparse cod_yyparse
#define yylex cod_yylex
#define yyrestart cod_yyrestart
#define yywrap cod_yywrap
#define yyerror cod_yyerror
#define yylineno cod_yylineno
#define yy_flex_debug cod_yy_flex_debug
#define yy_create_buffer cod_yy_create_buffer
#define yy_delete_buffer cod_yy_delete_buffer
#define yy_flush_buffer cod_yy_flush_buffer
#define yy_init_buffer cod_yy_init_buffer
#define yy_load_buffer_state cod_yy_load_buffer_state
#define yy_scan_buffer cod_yy_scan_buffer
#define yy_scan_bytes cod_yy_scan_bytes
#define yy_scan_string cod_yy_scan_string
#define yy_switch_to_buffer cod_yy_switch_to_buffer
#define yychar cod_yychar
#define yyin cod_yyin
#define yyleng cod_yyleng
#define yylval cod_yylval
#define yynerrs cod_yynerrs
#define yyout cod_yyout
#define yytext cod_yytext
#define yyset_out cod_yyset_out
#define yyset_lineno cod_yyset_lineno
#define yyset_in cod_yyset_in
#define yyset_debug cod_yyset_debug
#define yyrealloc cod_yyrealloc
#define yyalloc cod_yyalloc
#define yyfree cod_yyfree
#define yypush_buffer_state cod_yypush_buffer_state
#define yypop_buffer_state cod_yypop_buffer_state
#define yylex_destroy cod_yylex_destroy
#define yyget_out cod_yyget_out
#define yyget_lineno cod_yyget_lineno
#define yyget_in cod_yyget_in
#define yyget_debug cod_yyget_debug
#define yyget_text cod_yyget_text
#define yyget_leng cod_yyget_leng

static char *create_string_from_yytext();
extern int yylex();
extern int yyparse();
static sm_ref yyparse_value;
static int yyerror_count = 1;
extern void yyerror(char *str);
static int parsing_type = 0;
static int parsing_param_spec = 0;
static cod_parse_context yycontext;
static sm_ref cod_build_parsed_type_node(cod_parse_context c, char *name, sm_list l);
static sm_list
cod_dup_list(sm_list list)
{
    sm_list ret_list, new_list = NULL;
    sm_list *last_p = &ret_list;
    while (list != NULL) {
	*last_p = new_list = malloc(sizeof(struct list_struct));
	last_p = &(new_list->next);
	new_list->node = cod_copy(list->node);
	list = list->next;
    }
    *last_p = NULL;
    return ret_list;
}
%}

%expect 1
%union {
    lx_info info;
    sm_ref reference;
    operator_t operator;
    sm_list list;
    char *string;
};

%token <info> ARROW
%token <info> LPAREN
%token <info> RPAREN
%token <info> LCURLY
%token <info> RCURLY
%token <info> COLON
%token <info> QUESTION
%token <info> LBRACKET
%token <info> RBRACKET
%token <info> DOT
%token <info> STAR
%token <info> AT
%token <info> SLASH
%token <info> MODULUS
%token <info> PLUS
%token <info> MINUS
%token <info> TILDE
%token <info> LEQ
%token <info> LT
%token <info> GEQ
%token <info> GT
%token <info> EQ
%token <info> NEQ
%token <info> LEFT_SHIFT
%token <info> RIGHT_SHIFT
%token <info> ASSIGN
%token <info> MUL_ASSIGN
%token <info> DIV_ASSIGN
%token <info> MOD_ASSIGN
%token <info> ADD_ASSIGN
%token <info> SUB_ASSIGN
%token <info> LEFT_ASSIGN
%token <info> RIGHT_ASSIGN
%token <info> AND_ASSIGN
%token <info> XOR_ASSIGN
%token <info> OR_ASSIGN
%token <info> LOG_OR
%token <info> LOG_AND
%token <info> ARITH_OR
%token <info> ARITH_AND
%token <info> ARITH_XOR
%token <info> INC_OP
%token <info> DEC_OP
%token <info> BANG
%token <info> SEMI
%token <info> IF
%token <info> ELSE
%token <info> FOR
%token <info> DO
%token <info> WHILE
%token <info> CHAR
%token <info> SHORT
%token <info> INT
%token <info> LONG
%token <info> UNSIGNED
%token <info> SIGNED
%token <info> FLOAT
%token <info> DOUBLE
%token <info> VOID
%token <info> STRING
%token <info> STATIC
%token <info> EXTERN_TOKEN
%token <info> STRUCT
%token <info> ENUM
%token <info> UNION
%token <info> CONST
%token <info> SIZEOF
%token <info> TYPEDEF
%token <info> RETURN_TOKEN
%token <info> CONTINUE
%token <info> BREAK
%token <info> GOTO
%token <info> PRINT
%token <info> COMMA
%token <info> DOTDOTDOT
%token <info> integer_constant
%token <info> character_constant
%token <info> string_constant
%token <info> floating_constant
%token <info> identifier_ref
%token <info> type_identifier
%token <info> enumeration_constant

%type <info> struct_or_union;
%type <info> assignment_operator
%type <info> unary_operator
%type <reference> compound_statement;
%type <reference> labeled_statement;
%type <list> declaration_list decls_stmts_list;
%type <list> declaration;
%type <reference> statement;
%type <reference> init_declarator;
%type <list> init_declarator_list;
%type <reference> declarator;
%type <reference> direct_declarator;
%type <list> enumerator_list;
%type <reference> enumerator;
%type <list> declaration_specifiers;
%type <list> type_qualifier_list;
%type <list> specifier_qualifier_list;
%type <list> struct_declaration_list;
%type <list> struct_declaration;
%type <list> struct_declarator_list;
%type <reference> struct_declarator;
%type <list> pointer;
%type <list> abstract_declarator;
%type <list> type_name;
%type <list> parameter_list;
%type <list> parameter_type_list;
%type <list> argument_expression_list;
%type <reference> struct_or_union_specifier
%type <reference> enum_specifier
%type <reference> parameter_declaration;
%type <reference> type_specifier;
%type <reference> type_qualifier;
%type <reference> storage_class_specifier;
%type <reference> initializer;
%type <list> initializer_list;
%type <list> designator_list;
%type <list> designation;
%type <reference> designator;
%type <reference> assignment_expression;
%type <reference> expression_statement;
%type <reference> selection_statement;
%type <reference> iteration_statement;
%type <reference> jump_statement;
%type <reference> conditional_expression;
%type <reference> constant_expression;
%type <reference> unary_expression;
%type <reference> cast_expression;
%type <reference> postfix_expression;
%type <reference> primary_expression;
%type <reference> expression;
%type <reference> expression_opt;
%type <reference> constant;
%type <reference> multiplicative_expression;
%type <reference> additive_expression;
%type <reference> shift_expression;
%type <reference> relational_expression;
%type <reference> equality_expression;
%type <reference> logical_or_expression;
%type <reference> logical_and_expression;
%type <reference> inclusive_or_expression;
%type <reference> exclusive_or_expression;
%type <reference> and_expression;
%type <reference> start
%%

start :
	declaration_list {
	    yyparse_value = (sm_ref)$1;
	}
	|
	compound_statement {
	    yyparse_value = $1;
	}
	;

primary_expression:
	identifier_ref {
	    $$ = cod_new_identifier();
	    $$->node.identifier.id = $1.string;
	    $$->node.identifier.lx_srcpos = $1.lx_srcpos;
	}
	|
	constant
	|
	LPAREN expression RPAREN
	{ $$ = $2; }
	;

/* missing ->  	| '(' type_name ')' '{' initializer_list '}'
	| '(' type_name ')' '{' initializer_list ',' '}'  */
postfix_expression:
	primary_expression
	|
	postfix_expression LBRACKET expression RBRACKET {
	    $$ = cod_new_element_ref();
	    $$->node.element_ref.lx_srcpos = $2.lx_srcpos;
	    $$->node.element_ref.expression = $3;
	    $$->node.element_ref.array_ref = $1;
	}
	|
	postfix_expression DOT identifier_ref {
	    $$ = cod_new_field_ref();
	    $$->node.field_ref.lx_srcpos = $2.lx_srcpos;
	    $$->node.field_ref.lx_field = $3.string;
	    $$->node.field_ref.struct_ref = $1;
	}
	|
	postfix_expression LPAREN argument_expression_list RPAREN {
	    $$ = cod_new_subroutine_call();
	    $$->node.subroutine_call.lx_srcpos = $2.lx_srcpos;
	    $$->node.subroutine_call.arguments = $3;
	    $$->node.subroutine_call.sm_func_ref = $1;
	}
	|
	postfix_expression LPAREN  RPAREN {
	    $$ = cod_new_subroutine_call();
	    $$->node.subroutine_call.lx_srcpos = $2.lx_srcpos;
	    $$->node.subroutine_call.arguments = NULL;
	    $$->node.subroutine_call.sm_func_ref = $1;
	}
	| postfix_expression ARROW identifier_ref {
	    $$ = cod_new_field_ref();
	    $$->node.field_ref.lx_srcpos = $2.lx_srcpos;
	    $$->node.field_ref.lx_field = $3.string;
	    $$->node.field_ref.struct_ref = $1;
	}
	| postfix_expression INC_OP {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_inc;
	    $$->node.operator.right = NULL;
	    $$->node.operator.left = $1;
	}
        | postfix_expression DEC_OP {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_dec;
	    $$->node.operator.right = NULL;
	    $$->node.operator.left = $1;
	}
	;

argument_expression_list:
	assignment_expression {
		$$ = malloc(sizeof(struct list_struct));
		$$->node = $1;
		$$->next = NULL;
	}
	|
	argument_expression_list COMMA assignment_expression {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = $3;
	    tmp->next->next = NULL;
	    $$ = $1;
	};


/* missing ALIGNOF */
unary_expression:
	postfix_expression
        | INC_OP unary_expression {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $1.lx_srcpos;
	    $$->node.operator.op = op_inc;
	    $$->node.operator.right = $2;
	    $$->node.operator.left = NULL;
	}
        | DEC_OP unary_expression {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $1.lx_srcpos;
	    $$->node.operator.op = op_dec;
	    $$->node.operator.right = $2;
	    $$->node.operator.left = NULL;
	}
	| unary_operator cast_expression {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $1.lx_srcpos;
	    $$->node.operator.op = $1.op;
	    $$->node.operator.right = $2;
	    $$->node.operator.left = NULL;
	}
	| SIZEOF unary_expression {
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $1.lx_srcpos;
	    $$->node.operator.op = op_sizeof;
	    $$->node.operator.right = $2;
	    $$->node.operator.left = NULL;
	}
	| SIZEOF LPAREN type_name RPAREN {
	    /* dummy up a cast to hold the sm_list of the type */
	    sm_ref cast = cod_new_cast();
	    cast->node.cast.lx_srcpos = $1.lx_srcpos;
	    cast->node.cast.type_spec = $3;
	    cast->node.cast.expression = NULL;

	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $1.lx_srcpos;
	    $$->node.operator.op = op_sizeof;
	    $$->node.operator.right = cast;
	    $$->node.operator.left = NULL;
	};

unary_operator
	: ARITH_AND {
	    $$.op = op_address;
	}
	| STAR {
	    $$.op = op_deref;
	}
	| PLUS {
	    $$.op = op_plus;
	}
	| MINUS {
	    $$.op = op_minus;
	}
	| TILDE {
	    $$.op = op_not;
	}
	| BANG {
	    $$.op = op_log_neg;
	  }
	;

cast_expression:
	unary_expression
	| LPAREN type_name RPAREN cast_expression {
	    $$ = cod_new_cast();
	    $$->node.cast.lx_srcpos = $1.lx_srcpos;
	    $$->node.cast.type_spec = $2;
	    $$->node.cast.expression = $4;
	}
	;

multiplicative_expression:
	cast_expression
	|
	multiplicative_expression STAR cast_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_mult;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	multiplicative_expression SLASH cast_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_div;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	multiplicative_expression MODULUS cast_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_modulus;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

additive_expression:
	multiplicative_expression
	|
	additive_expression PLUS multiplicative_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.op = op_plus;
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	additive_expression MINUS multiplicative_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_minus;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

shift_expression:
	additive_expression
	|
	shift_expression LEFT_SHIFT additive_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_left_shift;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	shift_expression RIGHT_SHIFT additive_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_right_shift;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

relational_expression:
	shift_expression
	|
	relational_expression LT shift_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_lt;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	relational_expression GT shift_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_gt;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	relational_expression LEQ shift_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_leq;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	relational_expression GEQ shift_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_geq;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

equality_expression:
	relational_expression
	|
	equality_expression EQ relational_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_eq;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	|
	equality_expression NEQ relational_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_neq;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

and_expression:
	equality_expression
	|
	and_expression ARITH_AND equality_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_arith_and;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

exclusive_or_expression:
	and_expression
	|
	exclusive_or_expression ARITH_XOR and_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_arith_xor;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

inclusive_or_expression:
	exclusive_or_expression
	|
	inclusive_or_expression ARITH_OR exclusive_or_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_arith_or;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

logical_and_expression:
	inclusive_or_expression
	|
	logical_and_expression LOG_AND inclusive_or_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_log_and;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

logical_or_expression:
	logical_and_expression
	|
	logical_or_expression LOG_OR logical_and_expression
	{
	    $$ = cod_new_operator();
	    $$->node.operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.operator.op = op_log_or;
	    $$->node.operator.right = $3;
	    $$->node.operator.left = $1;
	}
	;

conditional_expression:
	logical_or_expression
	| 
	logical_or_expression QUESTION expression COLON conditional_expression
	{
	    $$ = cod_new_conditional_operator();
	    $$->node.conditional_operator.lx_srcpos = $2.lx_srcpos;
	    $$->node.conditional_operator.condition = $1;
	    $$->node.conditional_operator.e1 = $3;
	    $$->node.conditional_operator.e2 = $5;
	}
	;

assignment_operator
	: ASSIGN
	{ $$ = $1; $$.op = op_eq;} 
	| MUL_ASSIGN
	{ $$ = $1; $$.op = op_mult;} 
	| DIV_ASSIGN
	{ $$ = $1; $$.op = op_div;} 
	| MOD_ASSIGN
	{ $$ = $1; $$.op = op_modulus;} 
	| ADD_ASSIGN
	{ $$ = $1; $$.op = op_plus;} 
	| SUB_ASSIGN
	{ $$ = $1; $$.op = op_minus;} 
	| LEFT_ASSIGN
	{ $$ = $1; $$.op = op_left_shift;} 
	| RIGHT_ASSIGN
	{ $$ = $1; $$.op = op_right_shift;} 
	| AND_ASSIGN
	{ $$ = $1; $$.op = op_arith_and;} 
	| XOR_ASSIGN
	{ $$ = $1; $$.op = op_arith_xor;} 
	| OR_ASSIGN
	{ $$ = $1; $$.op = op_arith_or;} 
	;

assignment_expression:
	conditional_expression
	{ $$ = $1;} 
	|
	unary_expression assignment_operator assignment_expression
	{
	    $$ = cod_new_assignment_expression();
	    $$->node.assignment_expression.lx_srcpos = $2.lx_srcpos;
	    $$->node.assignment_expression.left = $1;
	    $$->node.assignment_expression.right = $3;
	    $$->node.assignment_expression.op = $2.op;
	}
	;

expression
	: assignment_expression
	    {$$ = $1;}
	| expression COMMA assignment_expression
	{
	    $$ = cod_new_comma_expression();
	    $$->node.comma_expression.lx_srcpos = $2.lx_srcpos;
	    $$->node.comma_expression.left = $1;
	    $$->node.comma_expression.right = $3;
	}


constant_expression
        : conditional_expression
        ;

init_declarator_list
	: init_declarator {
		$$ = malloc(sizeof(struct list_struct));
		$$->node = $1;
		$$->next = NULL;
	}
	| init_declarator_list COMMA init_declarator {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp = tmp->next;
	    tmp->node = $3;
	    tmp->next = NULL;
	    $$ = $1;
	}
	;

/*
missing :
	| static_assert_declaration
*/
declaration
	: declaration_specifiers 
	     { 
		 if (parsing_type) {
		     yyparse_value = (sm_ref) $1;
		     YYACCEPT;
		 }
	     }
	init_declarator_list
	    {  /* stop here if we're just doing a proc decl */
		if (parsing_param_spec) {
		    $<reference>$ = $3->node;
		    if ($<reference>$->node_type == cod_declaration) {
			if  ($<reference>$->node.declaration.type_spec == NULL) {
			    $<reference>$->node.declaration.type_spec = $1;
			} else {
			    /* 
			     * the pointer type list (with the declarator)
			     * goes at the end 
			     */
			    sm_list tmp = $1;
			    while (tmp->next != NULL) {
				tmp = tmp->next;
			    }
			    tmp->next = $<reference>$->node.declaration.type_spec;
			    $<reference>$->node.declaration.type_spec = $1;
			}
		    } else {
		        printf("unexpected node in init_declarator\n");
			cod_print($<reference>$);
		    }
		    yyparse_value = $3->node;
		    free($3);
		    YYACCEPT;
		}
	    }
	SEMI
	    {
		$$ = $3;
		sm_list dtmp = $3;
		while (dtmp) {
		    sm_list type_spec;
		    if (dtmp->next != NULL) {
			type_spec = cod_dup_list($1);
		    } else {
			type_spec = $1;
		    }
		    sm_ref decl = dtmp->node;
		    if (decl->node_type == cod_declaration) {
			if  (decl->node.declaration.type_spec == NULL) {
			    decl->node.declaration.type_spec = type_spec;
			} else {
			    /* 
			     * the pointer type list (with the declarator)
			     * goes at the end 
			     */
			    sm_list tmp = type_spec;
			    while (tmp->next != NULL) {
				tmp = tmp->next;
			    }
			    tmp->next = decl->node.declaration.type_spec;
			    decl->node.declaration.type_spec = type_spec;
			}
		    } else if (decl->node_type == cod_array_type_decl) {
			if  (decl->node.array_type_decl.type_spec == NULL) {
			    decl->node.array_type_decl.type_spec = type_spec;
			} else {
			    /* 
			     * the pointer type list (with the declarator)
			     * goes at the end 
			     */
			    sm_list tmp = type_spec;
			    while (tmp->next != NULL) {
				tmp = tmp->next;
			    }
			    tmp->next = decl->node.array_type_decl.type_spec;
			    decl->node.array_type_decl.type_spec = type_spec;
			}
		    } else {
			printf("Unknown decl entry\n");
			cod_print(decl);
		    }
		    while (type_spec != NULL) {
			if (type_spec->node->node.type_specifier.token == TYPEDEF) {
			    cod_add_defined_type(decl->node.declaration.id, yycontext);
			}
			type_spec = type_spec->next;
		    }
		    dtmp = dtmp->next;
		}
		(void)$<reference>4;
	    }
	| declaration_specifiers SEMI {
	    $$ = $1;
	}
	;

declaration_specifiers
	: storage_class_specifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| storage_class_specifier declaration_specifiers {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = $2;
	    $$ = tmp;
	}
	| type_specifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| type_specifier declaration_specifiers {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = $2;
	    $$ = tmp;
	}
	| type_qualifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| type_qualifier declaration_specifiers {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = $2;
	    $$ = tmp;
	};

init_declarator:
	declarator
	|
	declarator ASSIGN initializer
	    {
		if ($1->node_type == cod_declaration) {
		    $1->node.declaration.init_value = $3;
		} else if ($1->node_type == cod_array_type_decl) {
		    sm_ref tmp = $1->node.array_type_decl.element_ref;
		    while (tmp->node_type == cod_array_type_decl) {
			tmp = tmp->node.array_type_decl.element_ref;
		    }
		    assert(tmp->node_type == cod_declaration);
		    tmp->node.declaration.init_value = $3;
		}
	    }
	;

storage_class_specifier
	: TYPEDEF {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = TYPEDEF;
	}
	| STATIC {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = STATIC;
	}
	| EXTERN_TOKEN {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = EXTERN_TOKEN;
	}
	;

/* missing BOOL  COMPLEX IMAGINARY */
type_specifier:
	CHAR {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = CHAR;
	}
	| SHORT {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = SHORT;
	}
	| INT {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = INT;
	}
	| LONG {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = LONG;
	}
	| FLOAT {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = FLOAT;
	}
	| DOUBLE {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = DOUBLE;
	}
	| VOID {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = VOID;
	}
	| SIGNED {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = SIGNED;
	}
	| UNSIGNED {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = UNSIGNED;
	}
	| STRING {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = STRING;
	}
	| type_identifier {
	    $$ = cod_new_identifier();
	    $$->node.identifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.identifier.id = $1.string;
	}
	| struct_or_union_specifier {
	    $$ = $1;
	}
	| enum_specifier {
	    $$ = $1;
	}
	;

struct_or_union_specifier
	: struct_or_union identifier_ref LCURLY struct_declaration_list RCURLY {
	    $$ = cod_build_parsed_type_node(yycontext, $2.string, $4);
	}
	| struct_or_union LCURLY struct_declaration_list RCURLY {
	    $$ = cod_build_parsed_type_node(yycontext, strdup("anon"), $3);
	}
	| struct_or_union identifier_ref {
	    $$ = cod_build_parsed_type_node(yycontext, $2.string, NULL);
	}
	;

struct_or_union
	: STRUCT
	| UNION {
            yyerror("UNIONs not supported!");
	}
	;

struct_declaration_list
	: struct_declaration
	| struct_declaration_list struct_declaration {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next =$2;
	    $$ = $1;
	}
	;

/* missing static_assert_declaration */
struct_declaration
	: specifier_qualifier_list SEMI { }
	| specifier_qualifier_list struct_declarator_list SEMI {
	    sm_list type_spec = $1;
	    sm_list decl_list = $2;
 	    $$ = $2;
/******** GSE This isn't right.  Reusing potentially modified type spec */
	    while (decl_list != NULL) {
		sm_ref decl = decl_list->node;
		if (decl->node_type == cod_declaration) {
		    if  (decl->node.declaration.type_spec == NULL) {
			decl->node.declaration.type_spec = type_spec;
		    } else {
			/* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
			sm_list tmp = $1;
			while (tmp->next != NULL) {
			    tmp = tmp->next;
			}
			tmp->next = decl->node.declaration.type_spec;
			decl->node.declaration.type_spec = type_spec;
		    }
		} else if (decl->node_type == cod_array_type_decl) {
		    if  (decl->node.array_type_decl.type_spec == NULL) {
			decl->node.array_type_decl.type_spec = type_spec;
		    } else {
			/* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
			sm_list tmp = type_spec;
			while (tmp->next != NULL) {
			    tmp = tmp->next;
			}
			tmp->next = decl->node.array_type_decl.type_spec;
			decl->node.array_type_decl.type_spec = type_spec;
		    }
		} else {
		    printf("Unknown decl entry\n");
		    cod_print(decl);
		}
		decl_list = decl_list->next;
		if (decl_list != NULL) {
		    type_spec = cod_dup_list(type_spec);
		}
	    }
	}
	;


struct_declarator_list
 	: struct_declarator {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| struct_declarator_list COMMA struct_declarator {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = $3;
	    tmp->next->next = NULL;
	    $$ = $1;
	}
	;

struct_declarator : declarator ;

specifier_qualifier_list
	: type_specifier specifier_qualifier_list {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = $2;
	    $$ = tmp;
	}
	| type_specifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| type_qualifier specifier_qualifier_list {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = $2;
	    $$ = tmp;
	}
	| type_qualifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	;

enum_specifier
	: ENUM LCURLY enumerator_list RCURLY {
	    $$ = cod_new_enum_type_decl();
	    $$->node.enum_type_decl.id = gen_anon();
	    $$->node.enum_type_decl.enums = $3;
	    $$->node.enum_type_decl.lx_srcpos = $1.lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
	| ENUM LCURLY enumerator_list COMMA RCURLY {
	    $$ = cod_new_enum_type_decl();
	    $$->node.enum_type_decl.id = gen_anon();
	    $$->node.enum_type_decl.enums = $3;
	    $$->node.enum_type_decl.lx_srcpos = $1.lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
	| ENUM identifier_ref LCURLY enumerator_list RCURLY {
	    $$ = cod_new_enum_type_decl();
	    $$->node.enum_type_decl.id = $2.string;
	    $$->node.enum_type_decl.enums = $4;
	    $$->node.enum_type_decl.lx_srcpos = $1.lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
	| ENUM identifier_ref LCURLY enumerator_list COMMA RCURLY {
	    $$ = cod_new_enum_type_decl();
	    $$->node.enum_type_decl.id = $2.string;
	    $$->node.enum_type_decl.enums = $4;
	    $$->node.enum_type_decl.lx_srcpos = $1.lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
	| ENUM identifier_ref {
	    $$ = cod_new_enum_type_decl();
	    $$->node.enum_type_decl.id = $2.string;
	    $$->node.enum_type_decl.enums = NULL;
	    $$->node.enum_type_decl.lx_srcpos = $1.lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
	;

enumerator_list
	: enumerator {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = NULL;
	    $$ = tmp;
	}
	| enumerator_list COMMA enumerator {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $3;
	    tmp->next = $1;
	    $$ = tmp;
	}

	;

enumerator	/* identifiers must be flagged as ENUMERATION_CONSTANT */
	: identifier_ref ASSIGN constant_expression {
	    $$ = cod_new_enumerator();
	    $$->node.enumerator.id = $1.string;
	    $$->node.enumerator.const_expression = $3;
	}

	| identifier_ref {
	    $$ = cod_new_enumerator();
	    $$->node.enumerator.id = $1.string;
	    $$->node.enumerator.const_expression = NULL;
	}
	;

type_qualifier
	: CONST {
	    $$ = cod_new_type_specifier();
	    $$->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    $$->node.type_specifier.token = CONST;
	}
	;

declarator:
	direct_declarator
	|
	pointer direct_declarator {
	    $$ = $2;
	    if ($$->node_type == cod_declaration) {
		$$->node.declaration.type_spec = $1;
	    } else if ($$->node_type == cod_array_type_decl) {
		$$->node.array_type_decl.type_spec = $1;
	    } else {
		printf("Unknown direct_declarator entry\n");
		cod_print($$);
	    }
	}
	;

direct_declarator:
	identifier_ref	    {
		$$ = cod_new_declaration();
		$$->node.declaration.param_num = -1;
		$$->node.declaration.id = $1.string;
		$$->node.declaration.init_value = NULL;
		$$->node.declaration.lx_srcpos = $1.lx_srcpos;
		$$->node.declaration.is_subroutine = 0;
		$$->node.declaration.params = NULL;
	    }
	| LPAREN declarator RPAREN {
	    $$ = $2;
	}	
	| identifier_ref LPAREN parameter_type_list RPAREN {
		$$ = cod_new_declaration();
		$$->node.declaration.param_num = -1;
		$$->node.declaration.id = $1.string;
		$$->node.declaration.init_value = NULL;
		$$->node.declaration.lx_srcpos = $1.lx_srcpos;
		$$->node.declaration.is_subroutine = 1;
		$$->node.declaration.params = $3;
	}
	| identifier_ref LPAREN RPAREN {
		$$ = cod_new_declaration();
		$$->node.declaration.param_num = -1;
		$$->node.declaration.id = $1.string;
		$$->node.declaration.init_value = NULL;
		$$->node.declaration.lx_srcpos = $1.lx_srcpos;
		$$->node.declaration.is_subroutine = 1;
		$$->node.declaration.params = NULL;
	}
	| direct_declarator LBRACKET constant_expression RBRACKET {
		$$ = cod_new_array_type_decl();
		$$->node.array_type_decl.lx_srcpos = $2.lx_srcpos;
		$$->node.array_type_decl.size_expr = $3;
		$$->node.array_type_decl.element_ref = $1;
		$$->node.array_type_decl.sm_dynamic_size = NULL;
	}
	| direct_declarator LBRACKET RBRACKET {
		$$ = cod_new_array_type_decl();
		$$->node.array_type_decl.lx_srcpos = $2.lx_srcpos;
		$$->node.array_type_decl.size_expr = NULL;
		$$->node.array_type_decl.element_ref = $1;
		$$->node.array_type_decl.sm_dynamic_size = NULL;
	}
	;

pointer :
	STAR {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = NULL;
	}
	| STAR type_qualifier_list {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	| STAR pointer {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	| STAR type_qualifier_list pointer {
	    sm_list tmp = $2;
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = STAR;

	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = $3;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	| AT {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) { 
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = AT;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = NULL;
	}
	| AT type_qualifier_list {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = AT;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	| AT pointer {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = AT;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	| AT type_qualifier_list pointer {
	    sm_list tmp = $2;
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = $1.lx_srcpos;
	    star->node.type_specifier.token = AT;

	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = $3;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = star;
	    $$->next = $2;
	}
	;

type_qualifier_list
	: type_qualifier {
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = $1;
	    $$->next = NULL;
	}
	| type_qualifier_list type_qualifier {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = $2;
	    tmp->next->next = NULL;
	    $$ = $1;
	}
	;

parameter_type_list:
	parameter_list |
	parameter_list COMMA DOTDOTDOT {
	    sm_list tmp = $1;
	    sm_ref id = cod_new_declaration();
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = id;
	    tmp->next->next = NULL;
	    id->node.declaration.id = strdup("...");
	    $$ = $1;
	}
	;

parameter_list:
	parameter_declaration {
		$$ = malloc(sizeof(struct list_struct));
		$$->node = $1;
		$$->next = NULL;
	}
	|
	parameter_list COMMA parameter_declaration {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = $3;
	    tmp->next->next = NULL;
	    $$ = $1;
	};
	
/* missing
 	| declaration_specifiers abstract_declarator
*/
parameter_declaration
	: declaration_specifiers {
	    $$ = cod_new_declaration();
	    $$->node.declaration.param_num = -1;
	    $$->node.declaration.id = gen_anon();
	    $$->node.declaration.init_value = NULL;
	    $$->node.declaration.is_subroutine = 0;
	    $$->node.declaration.params = NULL;
	    $$->node.declaration.type_spec = $1;
	}
	| declaration_specifiers declarator {
		$$ = $2;
		if ($$->node_type == cod_declaration) {
		    $$->node.declaration.static_var = 0;
		    if  ($$->node.declaration.type_spec == NULL) {
		        $$->node.declaration.type_spec = $1;
		    } else {
		        /* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
		      sm_list tmp = $1;
		      while (tmp->next != NULL) {
			  tmp = tmp->next;
		      }
		      tmp->next = $$->node.declaration.type_spec;
		      $$->node.declaration.type_spec = $1;
		    }
		} else if ($$->node_type == cod_array_type_decl) {
		    if  ($$->node.array_type_decl.type_spec == NULL) {
		        $$->node.array_type_decl.type_spec = $1;
		    } else {
		        /* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
		      sm_list tmp = $1;
		      while (tmp->next != NULL) {
			  tmp = tmp->next;
		      }
		      tmp->next = $$->node.array_type_decl.type_spec;
		      $$->node.array_type_decl.type_spec = $1;
		    }
		} else {
		    printf("unexpected node in parameter_declaration");
		}
	};

type_name: 
	specifier_qualifier_list
	| specifier_qualifier_list abstract_declarator {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = $2;
	    $$ = $1;
	}
	;

/* missing 
	: pointer direct_abstract_declarator
	| direct_abstract_declarator

identifier_list
	: IDENTIFIER
	| identifier_list COMMA IDENTIFIER
	;

direct_abstract_declarator
	: '(' abstract_declarator ')'
	| '[' ']'
	| '[' '*' ']'
	| '[' STATIC type_qualifier_list assignment_expression ']'
	| '[' STATIC assignment_expression ']'
	| '[' type_qualifier_list STATIC assignment_expression ']'
	| '[' type_qualifier_list assignment_expression ']'
	| '[' type_qualifier_list ']'
	| '[' assignment_expression ']'
	| direct_abstract_declarator '[' ']'
	| direct_abstract_declarator '[' '*' ']'
	| direct_abstract_declarator '[' STATIC type_qualifier_list assignment_expression ']'
	| direct_abstract_declarator '[' STATIC assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list STATIC assignment_expression ']'
	| direct_abstract_declarator '[' type_qualifier_list ']'
	| direct_abstract_declarator '[' assignment_expression ']'
	| '(' ')'
	| '(' parameter_type_list ')'
	| direct_abstract_declarator '(' ')'
	| direct_abstract_declarator '(' parameter_type_list ')'
	;

*/

abstract_declarator:
	pointer
	;

initializer:
	LCURLY initializer_list RCURLY 
	{ 
	    $$ = cod_new_initializer_list();
	    $$->node.initializer_list.initializers = $2;
	}
	| LCURLY initializer_list COMMA RCURLY
	{ 
	    $$ = cod_new_initializer_list();
	    $$->node.initializer_list.initializers = $2;
	}
	| assignment_expression { $$ = $1;}
	;


initializer_list :
	designation initializer {
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = $1;
	    initializer->node.initializer.initializer = $2;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = initializer;
	    $$->next = NULL;
	}
	| initializer {
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = NULL;
	    initializer->node.initializer.initializer = $1;
	    $$ = malloc(sizeof(struct list_struct));
	    $$->node = initializer;
	    $$->next = NULL;
	}
	| initializer_list COMMA designation initializer {
	    sm_list tmp = $1;
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = $3;
	    initializer->node.initializer.initializer = $4;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = initializer;
	    tmp->next->next = NULL;
	    $$ = $1;
	}
	| initializer_list COMMA initializer {
	    sm_list tmp = $1;
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = NULL;
	    initializer->node.initializer.initializer = $3;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = initializer;
	    tmp->next->next = NULL;
	    $$ = $1;
	}
	;

designation
	: designator_list ASSIGN
	{ $$ = $1;}
	;

designator_list
	: designator {
		$$ = malloc(sizeof(struct list_struct));
		$$->node = $1;
		$$->next = NULL;
	}
	| designator_list designator {
	    sm_list tmp = $1;
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = $2;
	    tmp->next->next = NULL;
	    $$ = $1;
	}
	;

designator: 
	LBRACKET constant_expression RBRACKET
	{ 
	    $$ = cod_new_designator();
	    $$->node.designator.expression = $2;
	    $$->node.designator.id = NULL;
	}
	| DOT identifier_ref
	{ 
	    $$ = cod_new_designator();
	    $$->node.designator.expression = NULL;
	    $$->node.designator.id = $2.string;
	}
	;

decls_stmts_list:
	statement {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $1;
	    tmp->next = NULL;
	    $$ = tmp;
	} 
	|  declaration {
	    $$ = $1;
	   }
	| error SEMI {
	      $$ = NULL;
	  }
	| decls_stmts_list statement {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = $2;
	    tmp->next = NULL;
	    $$ = cod_append_list($1, tmp);
	}
	| decls_stmts_list declaration {
	    $$ = cod_append_list($1, $2);
	};

statement:
	labeled_statement
	| compound_statement
	| expression_statement
	| selection_statement
	| iteration_statement
	| jump_statement
	;

/* missing
	| CASE constant_expression COLON statement
	| DEFAULT COLON statement*/
labeled_statement:
	identifier_ref COLON statement {
	    $$ = cod_new_label_statement();
	    $$->node.label_statement.name =  $1.string;
	    $$->node.label_statement.statement = $3;
	};

compound_statement:
	LCURLY RCURLY {
	    $$ = cod_new_compound_statement();
	}
	| LCURLY decls_stmts_list RCURLY {
	    int count = $1.type_stack_count;
	    $$ = cod_new_compound_statement();
	    $$->node.compound_statement.decls = $2;
	    cod_remove_defined_types(yycontext, count);
	};

declaration_list:
	declaration { $$ = $1; }
	|
	declaration_list declaration {
	    if ($1 == NULL) {
		$$ = $2;
	    } else {
		sm_list tmp = $1;
		while (tmp->next != NULL) {
		    tmp = tmp->next;
		}
		tmp->next = $2;
		$$ = $1;
	    }
	};

jump_statement:
	RETURN_TOKEN expression SEMI {
	    $$ = cod_new_return_statement();
	    $$->node.return_statement.expression = $2;
	    $$->node.return_statement.lx_srcpos = $1.lx_srcpos;
	}
	| RETURN_TOKEN SEMI {
	    $$ = cod_new_return_statement();
	    $$->node.return_statement.expression = NULL;
	    $$->node.return_statement.lx_srcpos = $1.lx_srcpos;
	}
	| CONTINUE SEMI {
	    $$ = cod_new_jump_statement();
	    $$->node.jump_statement.continue_flag = 1;
	    $$->node.jump_statement.goto_target = NULL;
	    $$->node.jump_statement.lx_srcpos = $1.lx_srcpos;
	}
	| BREAK SEMI {
	    $$ = cod_new_jump_statement();
	    $$->node.jump_statement.continue_flag = 0;
	    $$->node.jump_statement.goto_target = NULL;
	    $$->node.jump_statement.lx_srcpos = $1.lx_srcpos;
	}
	| GOTO identifier_ref SEMI{
	    $$ = cod_new_jump_statement();
	    $$->node.jump_statement.continue_flag = 0;
	    $$->node.jump_statement.goto_target = $2.string;
	    $$->node.jump_statement.lx_srcpos = $1.lx_srcpos;
	}
	;

expression_statement:
	SEMI {
	    $$ = NULL;
	}
	|
	expression SEMI	{ 
	    $$ = cod_new_expression_statement();
	    $$->node.expression_statement.expression = $1;
	}
	;

/* missing switch 
	| SWITCH LPAREN expression RPAREN statement
*/
selection_statement:
	IF LPAREN expression RPAREN statement
	{ 
	    $$ = cod_new_selection_statement();
	    $$->node.selection_statement.lx_srcpos = $1.lx_srcpos;
	    $$->node.selection_statement.conditional = $3;
	    $$->node.selection_statement.then_part = $5;
	    $$->node.selection_statement.else_part = NULL;
	}
	| 
	IF LPAREN expression RPAREN statement ELSE statement
	{ 
	    $$ = cod_new_selection_statement();
	    $$->node.selection_statement.lx_srcpos = $1.lx_srcpos;
	    $$->node.selection_statement.conditional = $3;
	    $$->node.selection_statement.then_part = $5;
	    $$->node.selection_statement.else_part = $7;
	}
	;

/* Could do this instead 
	| FOR '(' expression_statement expression_statement ')' statement
	| FOR '(' expression_statement expression_statement expression ')' statement
	| FOR '(' declaration expression_statement ')' statement
	| FOR '(' declaration expression_statement expression ')' statement
*/
iteration_statement:
	FOR LPAREN expression_opt SEMI expression_opt SEMI expression_opt RPAREN statement
	{ 
	    $$ = cod_new_iteration_statement();
	    $$->node.iteration_statement.lx_srcpos = $1.lx_srcpos;
	    $$->node.iteration_statement.init_expr = $3;
	    $$->node.iteration_statement.test_expr = $5;
	    $$->node.iteration_statement.iter_expr = $7;
	    $$->node.iteration_statement.statement = $9;
	} 
	|
	WHILE LPAREN expression RPAREN statement
	{ 
	    $$ = cod_new_iteration_statement();
	    $$->node.iteration_statement.lx_srcpos = $1.lx_srcpos;
	    $$->node.iteration_statement.init_expr = NULL;
	    $$->node.iteration_statement.test_expr = $3;
	    $$->node.iteration_statement.iter_expr = NULL;
	    $$->node.iteration_statement.statement = $5;
	} 
	|
	DO statement WHILE LPAREN expression RPAREN SEMI
	{ 
	    $$ = cod_new_iteration_statement();
	    $$->node.iteration_statement.lx_srcpos = $1.lx_srcpos;
	    $$->node.iteration_statement.init_expr = NULL;
	    $$->node.iteration_statement.test_expr = NULL;
	    $$->node.iteration_statement.post_test_expr = $5;
	    $$->node.iteration_statement.iter_expr = NULL;
	    $$->node.iteration_statement.statement = $2;
	} 

	;

expression_opt:
	{ $$ = NULL; }
	/* null */
	| expression;

constant :
	integer_constant {
	    $$ = cod_new_constant();
	    $$->node.constant.token = integer_constant;
	    $$->node.constant.const_val = $1.string;
	    $$->node.constant.lx_srcpos = $1.lx_srcpos;
	}
	|
	floating_constant {
	    $$ = cod_new_constant();
	    $$->node.constant.token = floating_constant;
	    $$->node.constant.const_val = $1.string;
	    $$->node.constant.lx_srcpos = $1.lx_srcpos;
	}
	|
	string_constant {
	    $$ = cod_new_constant();
	    $$->node.constant.token = string_constant;
	    $$->node.constant.const_val = $1.string;
	    $$->node.constant.lx_srcpos = $1.lx_srcpos;
	}
	|
	character_constant {
	    $$ = cod_new_constant();
	    $$->node.constant.token = character_constant;
	    $$->node.constant.const_val = $1.string;
	    $$->node.constant.lx_srcpos = $1.lx_srcpos;
	}
	|
	enumeration_constant {
	    $$ = cod_new_constant();
	    $$->node.constant.token = character_constant;
	    $$->node.constant.const_val = $1.string;
	    $$->node.constant.lx_srcpos = $1.lx_srcpos;
	}
	;

%%
#ifdef _MSC_VER
#define YY_NO_UNISTD_H
#endif
#include "lex.yy.c"

typedef struct scope *scope_ptr;

struct parse_struct {
    sm_list decls;
    sm_list standard_decls;
    scope_ptr scope;
    char **defined_types;
    char **enumerated_constants;
    err_out_func_t error_func;
    void *client_data;
    sm_list return_type_list;
    int return_cg_type;
    sm_ref freeable_declaration;
    int has_exec_context;
    int dont_coerce_return;
    int alloc_globals;
};

static int
semanticize_compound_statement(cod_parse_context context, sm_ref compound, 
			       scope_ptr containing_scope, int require_last_return);
static int semanticize_decls_list(cod_parse_context context, sm_list decls, 
				  scope_ptr scope);
static int semanticize_array_type_node(cod_parse_context context,
				       sm_ref array, scope_ptr scope);
static int semanticize_reference_type_node(cod_parse_context context,
					   sm_ref decl, scope_ptr scope);
static void add_decl(char *id, sm_ref node, scope_ptr scope);
static sm_ref find_complex_type(sm_ref node, scope_ptr scope);
static const char *cod_code_string;
static int is_string(sm_ref expr);


int
cod_semanticize_added_decls(cod_parse_context context)
{
    return semanticize_decls_list(context, context->decls, context->scope);
}

extern void
cod_swap_decls_to_standard(cod_parse_context context)
{
    context->standard_decls = context->decls;
    context->decls = NULL;
}

void cod_set_error_func(cod_parse_context context, err_out_func_t err_func)
{
    context->error_func = err_func;
}

void cod_set_dont_coerce_return(cod_parse_context context, int value)
{
    context->dont_coerce_return = value;
}

static 
char *
cod_preprocessor(char *input, cod_parse_context context, int*white)
{
    char *out;
    char *ptr;
    if (strchr(input, '#') == NULL) return NULL;
    out = strdup(input);
    ptr = out;
    *white = 0;
    while (ptr && (*ptr)) {
	if (isspace(*ptr)) ptr++;
	if (*ptr == '#') {
	    char *start = ptr;
	    char *line_end;
	    if ((strncmp(ptr, "#include", 8) == 0) && isspace(*(ptr+8))) {
		/* got a #include */
		char *include_end;
		ptr += 8;
		while(isspace(*ptr)) ptr++;
		line_end = strchr(ptr, '\n');
		if (line_end) *line_end = 0;
		if ((*ptr == '<') || (*ptr == '"')) {
		    include_end = (*ptr == '<') ? strchr(ptr, '>') : strchr((ptr+1), '"');
		    if (!include_end) {
			printf("improper #include, \"%s\"\n", ptr);
			goto skip;
		    }
		    *include_end = 0;
		    cod_process_include((ptr+1), context);
		} else {
		    printf("improper #include, \"%s\"\n", ptr);
		    goto skip;
		}
		/* good #include, replace with spaces */
		if (line_end) *line_end = '\n';
		*include_end = ' ';
		while ((start != include_end) && (*start)) {
		    *(start++) = ' ';
		}
	    }
	}
    skip:
	/* skip to next line */
	ptr = strchr(ptr, '\n');
	while (ptr && (*(ptr - 1) == '\'')) {
	    /* continued line */
	    ptr = strchr(ptr, '\n');
	}
    }
    { 
	char *tmp = out;
	while(isspace(*tmp)) tmp++;
	if(*tmp == 0) {
	    free(out);
	    *white = 1;
	}
    }
    return out;
}
	
int
cod_parse_for_globals(char *code, cod_parse_context context)
{
    int ret;
    context->alloc_globals = 1;
    ret = cod_parse_for_context(code, context);
    context->alloc_globals = 0;
    return ret;
}
int
cod_parse_for_context(char *code, cod_parse_context context)
{
    sm_list decls;
    int ret;
    int all_whitespace = 0;
    char *freeable_code = NULL;
#if defined(YYDEBUG) && (YYDEBUG != 0)
    extern int yydebug;
    yydebug = 1;
#endif
    freeable_code = cod_preprocessor(code, context, &all_whitespace);
    if (all_whitespace) return 1;
    if (freeable_code) {
	code = freeable_code;
    }
    if (code != NULL) {
	setup_for_string_parse(code, context->defined_types, context->enumerated_constants);
	cod_code_string = code;
    }
    yyerror_count = 0;
    yycontext = context;
    yyparse();
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	if (freeable_code) free(freeable_code);
	return 0;
    }

    decls = (sm_list) yyparse_value;
    if (context->decls) {
	sm_list last = context->decls;
	while (last->next != NULL)
	    last = last->next;
	last->next = decls;
    } else {
	context->decls = decls;
    }
    ret = semanticize_decls_list(context, decls, context->scope);
    if (ret == 0) {
	cod_rfree_list(decls, NULL);
	context->decls = NULL;
    }
    if (freeable_code) free(freeable_code);
    return ret;
}

static int
semanticize_gotos(cod_parse_context context, sm_ref stmt, sm_list function_context);
static int semanticize_decl(cod_parse_context context, sm_ref decl, 
			    scope_ptr scope);

static int include_prefix(char *code)
{
    char *tmp = code;
    int not_done = 1;
    while (not_done) {
	while(isspace(*tmp)) tmp++;
	if (*tmp == '#') {
	    /* skip this line */
	    while(*tmp != '\n') tmp++;
	} else if (*tmp == '{') {
	    break;
	}
    }
    return (int)(intptr_t)(tmp - code);
}
cod_code
cod_code_gen(char *code, cod_parse_context context)
{
    sm_ref tmp, tmp2;
    cod_code ret_code;
    unsigned int offset;
    void *func;
    int bracket = 0;

    if (code != NULL) {
	if ((bracket = include_prefix(code))) {
	    char *prefix = malloc(bracket+1), *tmp;
	    strncpy(prefix, code, bracket + 1);
	    prefix[bracket] = 0;
	    tmp = prefix;
	    while(isspace(*tmp)) tmp++;
	    if (strlen(tmp) > 0) {
	        cod_parse_for_globals(tmp, context);
	    }
	    free(prefix);
	    code += bracket;
	}
	setup_for_string_parse(code, context->defined_types, context->enumerated_constants);
	cod_code_string = code;
    }

    yyerror_count = 0;
    yycontext = context;
    yyparse();
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	return 0;
    }
    tmp = cod_new_compound_statement();
    tmp->node.compound_statement.decls = context->decls;
    tmp->node.compound_statement.statements = NULL;
    tmp->node.compound_statement.statements =
	malloc(sizeof(struct list_struct));
    tmp->node.compound_statement.statements->next = NULL;
    tmp->node.compound_statement.statements->node = yyparse_value;
    tmp2 = cod_new_compound_statement();
    tmp2->node.compound_statement.decls = context->standard_decls;
    tmp2->node.compound_statement.statements =
	malloc(sizeof(struct list_struct));
    tmp2->node.compound_statement.statements->next = NULL;
    tmp2->node.compound_statement.statements->node = tmp;
    if (!semanticize_gotos(context, tmp, tmp2->node.compound_statement.statements) ||
	!semanticize_compound_statement(context, tmp, context->scope, (context->return_cg_type != DILL_V))) {
	tmp->node.compound_statement.decls = NULL;
	tmp2->node.compound_statement.decls = NULL;
	cod_rfree(tmp2);
	return NULL;
    }
    ret_code = malloc(sizeof(struct _cod_code_struct));
    memset(ret_code, 0, sizeof(struct _cod_code_struct));
    ret_code->code_memory_block = NULL;
    ret_code->data = NULL;
    ret_code->has_exec_ctx = context->has_exec_context;
    ret_code->static_block_address_register = -1;
    func = cod_cg_net(tmp, context->return_cg_type, &offset, ret_code);
    tmp->node.compound_statement.decls = NULL;
    tmp2->node.compound_statement.decls = NULL;
    cod_rfree(tmp2);
    ret_code->func = (void(*)(void))(intptr_t)func;
    return ret_code;
}

void 
cod_dump(cod_code code)
{
    printf("ECL CODE structure %p - \n", code);
    printf("  function pointer %p, code memory block %p, data %p, static size %d\n",
	   code->func, code->code_memory_block,
	   code->data, code->static_size_required);
#ifdef HAVE_DILL_H
    dill_dump((dill_stream) code->drisc_context);
#endif
}
    

int
cod_code_verify(char *code, cod_parse_context context)
{
    sm_ref tmp;

    if (code != NULL) {
	setup_for_string_parse(code, context->defined_types, context->enumerated_constants);
	cod_code_string = code;
    }

    yyerror_count = 0;
    yycontext = context;
    yyparse();
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	if (yyparse_value) {
	    cod_rfree(yyparse_value);
	}
	return 0;
    }

    tmp = cod_new_compound_statement();
    tmp->node.compound_statement.decls = context->decls;
    tmp->node.compound_statement.statements =
	malloc(sizeof(struct list_struct));
    tmp->node.compound_statement.statements->next = NULL;
    tmp->node.compound_statement.statements->node = yyparse_value;
    if (semanticize_compound_statement(context, tmp, context->scope, (context->return_cg_type != DILL_V)) == 0) {
	tmp->node.compound_statement.decls = NULL;
	cod_rfree(tmp);
	return 0;
    }
    tmp->node.compound_statement.decls = NULL;
    cod_rfree(tmp);
    return 1;
}

extern void 
cod_code_free(cod_code code)
{
    if (code->code_memory_block) free(code->code_memory_block);
    if (code->data) free(code->data);
#if defined(HAVE_DILL_H)
    if (code->drisc_context) {
	dill_free_stream((dill_stream) code->drisc_context);
    }
    if (code->execution_handle) {
	dill_free_handle((dill_exec_handle) code->execution_handle);
    }
#endif
    free(code);
}

static char *
copy_line(const char *line_begin)
{
    const char *line_end;
    if ((line_end = strchr(line_begin, 10)) == NULL) {
	/* no CR */
	return strdup(line_begin);
    } else {
	char *tmp = malloc(line_end - line_begin + 1);
	strncpy(tmp, line_begin, line_end - line_begin);
	tmp[line_end - line_begin] = 0;
	return tmp;
    }
}

static void
default_error_out(void *client_data, char *string)
{
    fprintf(stderr, "%s", string);
}

static void
print_context(cod_parse_context context, int line, int character)
{
    const char *tmp = cod_code_string;
    const char *line_begin = cod_code_string;
    char *line_copy = NULL;
    int i, line_len, offset = 0;

    while (line > 1) {
	switch(*tmp) {
	case 10:
	    line_begin = tmp + 1;
	    line--;
	    break;
	case 0:
	    line = 1;   /* end of src */
	    break;
	}
	tmp++;
    }
    if (character > 40) {
	offset = character - 40;
    }
    line_copy = copy_line(line_begin + offset);
    line_len = (int)strlen(line_copy);
    if (line_len > 60) {
	line_copy[60] = 0;
    }
    context->error_func(context->client_data, line_copy);
    context->error_func(context->client_data, "\n");
    free(line_copy);
    for(i=offset + 1; i< character; i++) {
	if (line_begin[i-1] == '\t') {
	    context->error_func(context->client_data, "\t");
	} else {
	    context->error_func(context->client_data, " ");
	}
    }
    context->error_func(context->client_data, "^\n");
}

void yyerror(char *str)
{
    char tmp_str[100];
    sprintf(tmp_str, "## Error %s\n", str);
    yycontext->error_func(yycontext->client_data, tmp_str);
    yycontext->error_func(yycontext->client_data, "## While parsing near ");
    yycontext->error_func(yycontext->client_data, yytext);
    sprintf(tmp_str, ", offset = %d, line = %d ####\n",lex_offset,line_count);
    yycontext->error_func(yycontext->client_data, tmp_str);
    print_context(yycontext, line_count, lex_offset);
    yyerror_count++;
}

#ifdef STDC_HEADERS
static void
cod_src_error(cod_parse_context context, sm_ref expr, char *format, ...)
#else
static void
cod_src_error(context, expr, format, va_alist)
cod_parse_context context;
sm_ref expr;
char *format;
va_dcl
#endif
{

    va_list ap;
    char *tmp = malloc(10240); /* arbitrarily large */
    srcpos lx_srcpos = {0,0};
#ifdef STDC_HEADERS
    va_start(ap, format);
#else
    va_start(ap);
#endif
    if (expr) lx_srcpos = cod_get_srcpos(expr);
    context->error_func(context->client_data, "## Ecode Error:  ");
    vsprintf(tmp, format, ap);
    context->error_func(context->client_data, tmp);
    sprintf(tmp, " at line %d, char %d\n", lx_srcpos.line, lx_srcpos.character);
    context->error_func(context->client_data, tmp);
    free(tmp);
    print_context(context, lx_srcpos.line, lx_srcpos.character);
}

extern void
cod_print_dimen_p(dimen_p d)
{
    int i;
    if (!d) {
	printf("DIMENS NOT SET YET\n");
	return;
    }
    for (i=0; i < d->dimen_count; i++) {
	if (d->dimens[i].static_size != -1) {
	    printf("[%d]", d->dimens[i].static_size);
	} else {
	    sm_ref field = d->dimens[i].control_field;
	    printf("[%s]", field->node.field.name);
	}
    }
    printf("\n");
}

extern void
cod_print_operator_t(operator_t o)
{
    switch (o) {
    case  op_modulus:
	printf("MODULUS");
	break;
    case  op_plus:
	printf("PLUS");
	break;
    case  op_minus:
	printf("MINUS");
	break;
    case  op_leq:
	printf("LEQ");
	break;
    case  op_lt:
	printf("LESS THAN");
	break;
    case  op_geq:
	printf("GEQ");
	break;
    case  op_gt:
	printf("GREATER THAN");
	break;
    case  op_eq:
	printf("EQUAL");
	break;
    case  op_neq:
	printf("NOT EQUAL");
	break;
    case  op_log_or:
	printf("LOGICAL OR");
	break;
    case  op_log_and:
	printf("LOGICAL AND");
	break;
    case op_log_neg:
	printf("LOGICAL NEGATION");
	break;
    case op_arith_and:
	printf("ARITH AND");
	break;
    case op_arith_or:
	printf("ARITH OR");
	break;
    case op_arith_xor:
	printf("ARITH XOR");
	break;
    case op_left_shift:
	printf("LEFT SHIFT");
	break;
    case op_right_shift:
	printf("RIGHT SHIFT");
	break;
    case  op_mult:
	printf("MULTIPLY");
	break;
    case  op_div:
	printf("DIVISION");
	break;
    case  op_deref:
	printf("DEREFERENCE");
	break;
    case  op_inc:
	printf("INCREMENT");
	break;
    case  op_not:
	printf("BITWISE NOT");
	break;
    case  op_dec:
	printf("DECREMENT");
	break;
    case  op_address:
	printf("ADDRESS");
	break;
    case op_sizeof:
	printf("SIZEOF");
	break;
    }
}

extern void
cod_print_srcpos(srcpos pos)
{
    printf("line %d, char %d", pos.line, pos.character);
}

extern void
cod_print_enc_info(enc_info enc)
{
    if (enc == NULL) {
	printf("Not encoded");
    } else {
	switch(enc->byte_order) {
	case 1:
	    printf("Bigendian");
	    break;
	case 2:
	    printf("Littleendian");
	    break;
	}
    }
}

extern void
free_enc_info(enc_info enc)
{
    free(enc);
}

enum namespace { NS_DEFAULT, NS_STRUCT, NS_ENUM };

char *namespace_str[] = {"DEFAULT", "STRUCT"};

typedef struct st_entry {
    char *id;
    sm_ref node;
    enum namespace ns;
    struct st_entry *next;
} *st_entry;

struct scope {
    cod_extern_list externs;
    struct st_entry *entry_list;
    sm_ref code_container;
    struct scope *containing_scope;
};


extern cod_parse_context
cod_copy_context(cod_parse_context context)
{
    int i, count;
    int type_count = 0;
    cod_parse_context new_context = new_cod_parse_context();
    new_context->has_exec_context = context->has_exec_context;
    new_context->decls = cod_copy_list(context->decls);
    count = 0;
    while (context->scope->externs && context->scope->externs[count].extern_value) count++;
    i=0;
    while(new_context->scope->externs[i].extern_name) free(new_context->scope->externs[i++].extern_name);
    free(new_context->scope->externs);
    new_context->scope->externs = malloc(sizeof(context->scope->externs[0]) *
					 (count+1));
    for (i=0; i < count; i++) {
      new_context->scope->externs[i].extern_name = strdup(context->scope->externs[i].extern_name);
      new_context->scope->externs[i].extern_value = context->scope->externs[i].extern_value;
    }
    new_context->scope->externs[count].extern_name = NULL;
    new_context->scope->externs[count].extern_value = NULL;

    new_context->error_func = context->error_func;
    new_context->client_data = context->client_data;
    semanticize_decls_list(new_context, new_context->decls, 
			   new_context->scope);
    free(new_context->defined_types);
    while(context->defined_types && context->defined_types[type_count]) type_count++;
    new_context->defined_types = malloc(sizeof(char*) * (type_count + 2));
    for (i=0; i<= type_count; i++) {
	new_context->defined_types[i] = context->defined_types[i];
    }
    return new_context;
}

extern void dump_scope(scope_ptr scope);

extern cod_parse_context
cod_copy_globals(cod_parse_context context)
{
    int i, count;
    int type_count = 0;
    cod_parse_context new_context = new_cod_parse_context();
    new_context->has_exec_context = context->has_exec_context;
    new_context->decls = cod_copy_list(context->decls);
    sm_list new_decls = new_context->decls;
    sm_list old_decls = context->decls;
    sm_list *last_new_p;
    last_new_p = &new_context->decls;
    while(new_decls != NULL) {
	sm_ref new_decl = new_decls->node;
	sm_ref old_decl = old_decls->node;
	switch(new_decl->node_type) {
	case cod_declaration:
	    if ((old_decl->node.declaration.param_num != -1) || 
		(old_decl->node.declaration.cg_address == (void*)-1)){
		/* this is an old parameter or subroutine, we have to kill it */
		*last_new_p = new_decls->next;
		new_decls = new_decls->next;
		old_decls = old_decls->next;
		continue;
	    }
	    new_decl->node.declaration.cg_address = old_decl->node.declaration.cg_address;
	    new_decl->node.declaration.sm_complex_type = NULL;
	    break;
	case cod_array_type_decl:
	    new_decl->node.array_type_decl.element_ref->node.declaration.cg_address = 
		old_decl->node.array_type_decl.element_ref->node.declaration.cg_address;
	    break;
	default:
	    break;
	}
	last_new_p = &new_decls->next;
	new_decls = new_decls->next;
	old_decls = old_decls->next;
    }
    count = 0;
    while (context->scope->externs && context->scope->externs[count].extern_value) count++;
    i=0;
    while(new_context->scope->externs[i].extern_name) free(new_context->scope->externs[i++].extern_name);
    free(new_context->scope->externs);
    new_context->scope->externs = malloc(sizeof(context->scope->externs[0]) *
					 (count+1));
    for (i=0; i < count; i++) {
      new_context->scope->externs[i].extern_name = strdup(context->scope->externs[i].extern_name);
      new_context->scope->externs[i].extern_value = context->scope->externs[i].extern_value;
    }
    new_context->scope->externs[count].extern_name = NULL;
    new_context->scope->externs[count].extern_value = NULL;

    new_context->error_func = context->error_func;
    new_context->client_data = context->client_data;
    semanticize_decls_list(new_context, new_context->decls, 
			   new_context->scope);
    free(new_context->defined_types);
    while(context->defined_types && context->defined_types[type_count]) type_count++;
    new_context->defined_types = malloc(sizeof(char*) * (type_count + 2));
    for (i=0; i<= type_count; i++) {
	new_context->defined_types[i] = context->defined_types[i];
    }
    return new_context;
}

static sm_ref
find_containing_iterator(scope_ptr scope)
{
    if (scope == NULL) return NULL;
    if ((scope->code_container != NULL) &&
	(scope->code_container->node_type == cod_iteration_statement)) {
	return scope->code_container;
    }
    return find_containing_iterator(scope->containing_scope);
}

static void *
resolve_extern(char *id, scope_ptr scope)
{
    if (scope == NULL) return NULL;
    if (scope->externs != NULL) {
	cod_extern_list externs = scope->externs;
	while(externs->extern_name != NULL) {
	    if (strcmp(id, externs->extern_name) == 0) {
		return externs->extern_value;
	    }
	    externs++;
	}
    }
    return resolve_extern(id, scope->containing_scope);
}

static scope_ptr
push_scope_container(scope_ptr containing_scope, sm_ref container)
{
    scope_ptr new_scope = malloc(sizeof(*new_scope));
    new_scope->externs = NULL;
    new_scope->entry_list = NULL;
    new_scope->code_container = container;
    new_scope->containing_scope = containing_scope;
    return new_scope;
}

static scope_ptr
push_scope(scope_ptr containing_scope)
{
    scope_ptr new_scope = malloc(sizeof(*new_scope));
    new_scope->externs = NULL;
    new_scope->entry_list = NULL;
    new_scope->code_container = NULL;
    new_scope->containing_scope = containing_scope;
    return new_scope;
}

static void
pop_scope(scope_ptr scope)
{
    st_entry list = scope->entry_list;
    while (list != NULL) {
	st_entry tmp = list->next;
	free(list);
	list = tmp;
    }
    free(scope);
}

extern void
dump_scope(scope_ptr scope)
{
    printf("Containing_scope is %p\n", scope->containing_scope);
    printf("Extern list:");
    if (scope->externs != NULL) {
	int i = 0;
	while (scope->externs[i].extern_name != NULL) {
	    printf("\t\"%s\" -> 0x%p\n", scope->externs[i].extern_name,
		   scope->externs[i].extern_value);
	    i++;
	}
    }
    printf("Symbol list:");
    if (scope->entry_list != NULL) {
	st_entry e = scope->entry_list;
	while (e != NULL) {
	    printf("\t\"%s\" -> 0x%p   [%s]\n", e->id, e->node, namespace_str[e->ns]);
	    cod_print(e->node);
	    e = e->next;
	}
    }
}	
static void
add_decl(char *id, sm_ref node, scope_ptr scope)
{
    st_entry entry = malloc(sizeof(*entry));
    entry->node = node;
    entry->id = id;
    entry->ns = NS_DEFAULT;
    entry->next = scope->entry_list;
    scope->entry_list = entry;
}

static void
add_decl_ns(char *id, sm_ref node, scope_ptr scope, enum namespace ns)
{
    st_entry entry = malloc(sizeof(*entry));
    entry->node = node;
    entry->id = id;
    entry->ns = ns;
    entry->next = scope->entry_list;
    scope->entry_list = entry;
}

extern void
cod_add_decl_to_scope(char *id, sm_ref node, cod_parse_context context)
{
    add_decl(id, node, context->scope);
}

static sm_ref
resolve_local(char *id, scope_ptr scope)
{
    st_entry list = scope->entry_list;
    while(list != NULL) {
	if (strcmp(list->id, id) == 0) {
	    return list->node;
	}
	list = list->next;
    }
    return NULL;
}

static sm_ref
resolve(char *id, scope_ptr scope)
{
    sm_ref tmp;
    if (scope == NULL) return NULL;
    tmp = resolve_local(id, scope);
    if (tmp != NULL) {
	return tmp;
    }
    return resolve(id, scope->containing_scope);
}

static int
determine_unary_type(cod_parse_context context, sm_ref expr, sm_ref right)
{
    int right_type = cod_sm_get_type(right);
    operator_t op = expr->node.operator.op;
    switch(right_type) {
    case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	right_type = DILL_I;   /* integer promotion */
    }
    if (op == op_minus) {
	switch(right_type) {
	case DILL_U:
	    return DILL_I;
	case DILL_UL:
	    return DILL_L;
	}
    }
    return right_type;
}

static int 
determine_op_type(cod_parse_context context, sm_ref expr, 
		  sm_ref left, sm_ref right)
{
    int unsigned_used = 0;
    int left_type  = cod_sm_get_type(left);
    int right_type = cod_sm_get_type(right);
    operator_t op = expr->node.operator.op;

    if (left_type == DILL_P) {
	sm_ref ctype;
	if (is_string(expr->node.operator.left) &&
	    is_string(expr->node.operator.right) &&
	    (op == op_eq)) {
	    return DILL_I;
	}

	ctype = get_complex_type(context, left);
	if(ctype && (ctype->node_type == cod_struct_type_decl)) {
	    cod_src_error(context, expr, 
			  "Illegal arithmetic. Left side is a structure.");
	    return DILL_ERR;
	}

	switch(right_type) {
	case DILL_P:
	case DILL_C:
	case DILL_UC:
	case DILL_S:
	case DILL_US:
	case DILL_I:
	case DILL_U:
	case DILL_L:
	case DILL_UL:
	case DILL_B:
	    return DILL_P;
	    
	default:
	    cod_src_error(context, expr, 
			  "Illegal pointer arithmetic. Right side is of incompatible type.");
	    return DILL_ERR;
	}
    }
    if (right_type == DILL_P) {
	sm_ref right_complex = get_complex_type(context, right);
	if(right_complex && (right->node_type == cod_struct_type_decl)) {
	    cod_src_error(context, expr, 
			  "Illegal arithmetic. Right side is a structure.");
	    return DILL_ERR;
	}

	switch(left_type) {
	case DILL_P:
	case DILL_C:
	case DILL_UC:
	case DILL_S:
	case DILL_US:
	case DILL_I:
	case DILL_U:
	case DILL_L:
	case DILL_UL:
	    return DILL_P;
	    
	default:
	    cod_src_error(context, expr, 
			  "Illegal pointer arithmetic. Left side is of incompatible type.");
	    return DILL_ERR;
	}	
	
    }
    if (left_type == DILL_B) {
	cod_src_error(context, expr, 
		      "Illegal arithmetic.  Left side is a structured type");
	return DILL_ERR;
    }
    if (right_type == DILL_B) {
	cod_src_error(context, expr, 
		      "Illegal arithmetic.  Right side is a structured type");
	return DILL_ERR;
    }
    if ((left_type == DILL_D) || (right_type == DILL_D)) {
	if ((op == op_modulus) || (op == op_log_or) || (op == op_log_and)) {
	    cod_src_error(context, expr, "Operands must be integral.");
	    return DILL_ERR;
	} else {
	    return DILL_D;
	}
    }
    if ((left_type == DILL_F) || (right_type == DILL_F)) {
	if ((op == op_modulus) || (op == op_log_or) || (op == op_log_and)) {
	    cod_src_error(context, expr, "Operands must be integral.");
	    return DILL_ERR;

	} else {
	    return DILL_F;
	}
    }
    switch(left_type) {
    case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	left_type = DILL_I;   /* integer promotion */
    }
    switch(right_type) {
    case DILL_C: case DILL_UC: case DILL_S: case DILL_US:
	right_type = DILL_I;   /* integer promotion */
    }
    switch(left_type) {
    case DILL_UC: case DILL_US: case DILL_U: case DILL_UL:
	unsigned_used++;
    }
    switch(right_type) {
    case DILL_UC: case DILL_US: case DILL_U: case DILL_UL:
	unsigned_used++;
    }
    if ((op == op_left_shift) || (op == op_right_shift)) return left_type;
    if ((left_type == DILL_UL) || (right_type == DILL_UL)) return DILL_UL;
    if ((left_type == DILL_L) || (right_type == DILL_L)) {
	/* GSE -bug  This test should be for *generated* target, not host */
	if (sizeof(intptr_t) > sizeof(unsigned int)) {
	    /* Long can represent all values of unsigned int */
	    return DILL_L;
	} else {
	    return unsigned_used? DILL_UL : DILL_L;
	}
    }
    if ((left_type == DILL_U) || (right_type == DILL_U)) return DILL_U;
    return unsigned_used? DILL_U: DILL_I;
}

static sm_ref reduce_type_list(cod_parse_context context, sm_list type_list, 
			       int *cg_type, scope_ptr scope, int *is_typedef,
			       sm_ref *freeable_type);
static int 
assignment_types_match(cod_parse_context context, sm_ref left, sm_ref right, int strict);

#ifdef NOTDEF
static int
is_n_dimen_array(int dimen, sm_ref expr)
{
    if ((dimen == 0) && (expr == NULL)) return 1;
    if ((dimen > 0) && (expr == NULL)) return 0;
    if (dimen == 0) {
	if (expr->node_type == cod_array_type_decl) {
	    return 0;
	} else {
	    return 1;
	}
    }
    if (expr->node_type == cod_field_ref) {
	return is_n_dimen_array(dimen, expr->node.field_ref.sm_field_ref);
    }
    if (expr->node_type == cod_element_ref) {
	return is_n_dimen_array(dimen + 1, expr);
    }
    if (expr->node_type != cod_field) return 0;
    if (expr->node_type == cod_field) {
	return is_n_dimen_array(dimen, expr->node.field.sm_complex_type);
    }
    /* ought to recurse or handle above */
    assert(0);
    return 0;
}
#endif

static int 
is_string(sm_ref expr)
{
    if (expr->node_type == cod_field) {
	return (expr->node.field.string_type && (strcmp(expr->node.field.string_type, "string") == 0));
    } else if (expr->node_type == cod_field_ref) {
	return is_string(expr->node.field_ref.sm_field_ref);
    } else if (expr->node_type == cod_identifier) {
	return is_string(expr->node.identifier.sm_declaration);
    } else if (expr->node_type == cod_conditional_operator) {
	return is_string(expr->node.conditional_operator.e1);  /* e2 must be similar */
    } else if (expr->node_type == cod_declaration) {
	if (expr->node.declaration.cg_type != DILL_P) return 0;
	if (expr->node.declaration.sm_complex_type != NULL) return 0;
	/* only strings have pointers without complex types */
	return 1;
    } else if (expr->node_type == cod_constant) {
	return (expr->node.constant.token == string_constant);
    }
    return 0;
}

static int 
is_const(sm_ref expr)
{
    switch(expr->node_type) {
    case cod_field_ref:
	return is_const(expr->node.field_ref.struct_ref);
    case cod_element_ref:
	return is_const(expr->node.element_ref.array_ref);
    case cod_cast:
	return is_const(expr->node.cast.expression);
    case cod_identifier:
	return is_const(expr->node.identifier.sm_declaration);
    case cod_declaration:
	return expr->node.declaration.const_var;
    case cod_constant:
	return 1;
    case cod_operator:
	/* most operator errors handled elsewhere.  Here we're only looking for dereference */
	if (expr->node.operator.op == op_deref) {
	    return is_const(expr->node.operator.right);
	}
	return 0;
    default:
	printf("Unhandled case in is_const()\n");
	cod_print(expr);
	assert(0);
    }
    return 0;
}

extern int 
cod_expr_is_string(sm_ref expr)
{
    return is_string(expr);
}

extern int 
is_control_value(sm_ref expr, sm_ref strct)
{
    sm_list fields;
    if (expr->node_type == cod_field_ref) {
	return is_control_value(expr->node.field_ref.sm_field_ref,
				expr->node.field_ref.struct_ref);
    }
    if (expr->node_type != cod_field) return 0;
    assert(strct != NULL);
    strct = get_complex_type(0, strct);	
    if (strct->node_type == cod_reference_type_decl) {
	strct = strct->node.reference_type_decl.sm_complex_referenced_type;
    }
    if (strct->node_type == cod_declaration) {
	strct = strct->node.declaration.sm_complex_type;
    }
    assert(strct->node_type == cod_struct_type_decl);
    fields =  strct->node.struct_type_decl.fields;
    while(fields != NULL) {
	sm_ref ctype = fields->node->node.field.sm_complex_type;
	if ((ctype != NULL) && (ctype->node_type == cod_reference_type_decl))
	    ctype = ctype->node.reference_type_decl.sm_complex_referenced_type;
	while (ctype != NULL) {
	    if (ctype->node_type == cod_array_type_decl) {
		if (ctype->node.array_type_decl.sm_dynamic_size == expr) {
		    return 1;
		}
		ctype = ctype->node.array_type_decl.sm_complex_element_type;
	    } else {
		ctype = NULL;
	    }
	}
	fields = fields->next;
    }
    return 0;
}

#ifndef FALSE
#define FALSE 0
#endif

static char*
type_list_to_string(cod_parse_context context, sm_list type_list, int *size)
{
    sm_list orig_list = type_list;
    int short_appeared = 0;
    int long_appeared = 0;
    int long_long_appeared = 0;
    int int_appeared = 0;
    int double_appeared = 0;
    int float_appeared = 0;
    int char_appeared = 0;
    int signed_appeared = 0;
    int unsigned_appeared = 0;
    int void_appeared = 0;
    int string_appeared = 0;
    int spec_count = 0;
    int prefix_end = 0;
    int type_found = 0;
    int cg_type;

    cg_type = DILL_ERR;
    while ((type_list != NULL) && (prefix_end == 0)) {
	sm_ref node = type_list->node;
	int typ = -1;
	if (node->node_type == cod_type_specifier) {
	    typ = type_list->node->node.type_specifier.token;
	    if ((typ == STAR) || (typ == AT)) {
		prefix_end = 1;
		type_list = type_list->next;
		continue;
	    }
	}
	if (node->node_type != cod_type_specifier) {
	    if (node->node_type == cod_identifier) {
		return NULL;
	    } else if (node->node_type == cod_struct_type_decl) {
		return NULL;
	    } else {
		printf("Unknown node type in type_list_to_string\n");
		break;
	    }
	} else {
	    spec_count++;
	    switch (typ) {
	    case INT:
		int_appeared++;
		break;
	    case LONG:
		long_appeared++;
		break;
	    case SHORT:
		short_appeared++;
		break;
	    case DOUBLE:
		double_appeared++;
		break;
	    case STRING:
		string_appeared++;
		break;
	    case VOID:
		void_appeared++;
		break;
	    case FLOAT:
		float_appeared++;
		break;
	    case CHAR:
		char_appeared++;
		break;
	    case SIGNED:
		signed_appeared++;
		break;
	    case UNSIGNED:
		unsigned_appeared++;
		break;
	    case TYPEDEF:
		spec_count--;
		break;
	    case STATIC:
		spec_count--;
		break;
	    case EXTERN_TOKEN:
		spec_count--;
		break;
	    case CONST:
		spec_count--;
		break;
	    default:
		printf("Unknown type\n");
	    }
	    type_list = type_list->next;
	}
    }
    if (spec_count == 0) {
	if (type_list == NULL) cg_type = DILL_I;   /* default to int */
	goto finalize;
    }
    if (void_appeared && (spec_count > 1)) {
	cod_src_error(context, orig_list->node, 
		      "Void type may not appear with other specifiers");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (string_appeared && (spec_count > 1)) {
	cod_src_error(context, orig_list->node, 
		      "String type may not appear with other specifiers");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (void_appeared) {
	cg_type = DILL_V;
	goto finalize;
    }
    if (string_appeared) {
	cg_type = DILL_P;
	goto finalize;
    }
    if (short_appeared && long_appeared) {
	cod_src_error(context, orig_list->node, 
		      "Only one of long or short permitted");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (short_appeared && (double_appeared + float_appeared)) {
	cod_src_error(context, orig_list->node, 
		      "Short may not be specified with double or float");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (double_appeared + float_appeared) {
	if (double_appeared + float_appeared + short_appeared + signed_appeared + unsigned_appeared + char_appeared + int_appeared > 1) {
	    cod_src_error(context, orig_list->node, "Bad type spec");
	    cg_type = DILL_ERR;
	    return NULL;
	} else {
	    /* not handling LONG plus one of these */
	    if (double_appeared) {
		cg_type = DILL_D;
		goto finalize;
	    } else {
		cg_type = DILL_F;
		goto finalize;
	    }
	}
    }

    /* neither float or double appeared */
    if (long_appeared == 2) {
	long_long_appeared++;
	long_appeared = 0;
    }
    if (short_appeared + char_appeared + long_appeared + long_long_appeared >= 2) {
	cod_src_error(context, orig_list->node, 
		      "Only one integer size spec may be specified");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (unsigned_appeared + signed_appeared > 1) {
	cod_src_error(context, orig_list->node, "Bad type spec");
	cg_type = DILL_ERR;
	return NULL;
    }
    if (unsigned_appeared) {
	if (char_appeared) {
	    cg_type = DILL_UC;
	    goto finalize;
	} else if (short_appeared) {
	    cg_type = DILL_US;
	    goto finalize;
	} else if (long_appeared || long_long_appeared) {
	    cg_type = DILL_UL;
	    goto finalize;
	} else {
	    cg_type = DILL_U;
	    goto finalize;
	}
    } else {
	if (char_appeared) {
	    cg_type = DILL_C;
	    goto finalize;
	} else if (short_appeared) {
	    cg_type = DILL_S;
	    goto finalize;
	} else if (long_appeared || long_long_appeared) {
	    cg_type = DILL_L;
	    goto finalize;
	} else {
	    cg_type = DILL_I;
	    goto finalize;
	}
    }
 finalize:
    if (cg_type != DILL_ERR) {
	type_found++;
    }
    switch(cg_type) {
    case DILL_C: 
	*size = sizeof(char);
	return strdup("integer");
    case DILL_UC: 
	*size = sizeof(char);
	return strdup("unsigned integer");
    case DILL_I: 
	*size = sizeof(int);
	return strdup("integer");
    case DILL_L:
	*size = sizeof(intptr_t);
	return strdup("integer");
    case DILL_S:
	*size = sizeof(short);
	return strdup("integer");
    case DILL_U: 
	*size = sizeof(int);
	return strdup("unsigned integer");
    case DILL_UL:
	*size = sizeof(intptr_t);
	return strdup("unsigned integer");
    case DILL_US:
	*size = sizeof(short);
	return strdup("unsigned integer");
    case DILL_F:
	*size = sizeof(float);
	return strdup("float");
    case DILL_D:
	*size = sizeof(double);
	return strdup("float");
    }
    return NULL;
}

static sm_ref
cod_build_parsed_type_node(cod_parse_context c, char *name, sm_list l)
{
    sm_ref decl = cod_new_struct_type_decl();
    sm_list *end_ptr = &decl->node.struct_type_decl.fields;

    sm_list tmp = l;
    sm_list last_type = NULL;
    int field_count = 0;
    decl->node.struct_type_decl.id = name;
    
     while(tmp != NULL) {
	sm_ref node = tmp->node;
	sm_list typ = NULL;
	sm_list new_elem;
	new_elem = malloc(sizeof(*new_elem));
	new_elem->next = NULL;
	new_elem->node = cod_new_field();
	if (node->node_type == cod_declaration) {
	    typ = cod_dup_list(node->node.declaration.type_spec);
	    new_elem->node->node.field.name = strdup(node->node.declaration.id);
	    new_elem->node->node.field.string_type  = 
		type_list_to_string(c, typ, &new_elem->node->node.field.cg_size);
	} else if (node->node_type == cod_array_type_decl) {
	    sm_ref base_decl = node->node.array_type_decl.element_ref;
	    sm_ref size = node->node.array_type_decl.size_expr;
	    char *base_string_type = NULL;
	    char *size_str = NULL, *final_type;
	    typ = cod_dup_list(node->node.array_type_decl.type_spec);
	    if (base_decl->node_type != cod_declaration) {
		printf("Array base type must be a simple type\n");
		return NULL;
	    }
	    new_elem->node->node.field.name = strdup(base_decl->node.declaration.id);
	    base_string_type  = 
		type_list_to_string(c, typ, &new_elem->node->node.field.cg_size);
	    if (size->node_type == cod_identifier) {
		size_str = size->node.identifier.id;
	    } else {
		int free_val = 0;
		sm_ref constant = evaluate_constant_return_expr(c, size, &free_val);
		if (constant->node_type == cod_constant) {
		    if (constant->node.constant.token != integer_constant) {
			printf("Array size constant is non-integer\n");
			return NULL;
		    } else {
			size_str = constant->node.constant.const_val;
		    }
		    if (free_val) free(constant);
		} else {
		    printf("Unexpected value for array size\n");
		    return NULL;
		}
	    }
	    if (base_string_type) {
		final_type = malloc(strlen(base_string_type) + 
				    strlen(size_str) + 3);
		sprintf(final_type, "%s[%s]", base_string_type, size_str);
		new_elem->node->node.field.string_type = final_type;
		free(base_string_type);
	    } else {
		new_elem->node->node.field.string_type = NULL;
	    }
	}
	new_elem->node->node.field.cg_offset = -1;
	new_elem->node->node.field.cg_type = DILL_ERR;
	new_elem->node->node.field.type_spec = typ;
	cod_rfree(node);
	field_count++;
	last_type = tmp;
	tmp = tmp->next;
	free(last_type);
	*end_ptr = new_elem;
	end_ptr = &new_elem->next;
    }
    return decl;
}

int 
is_array(sm_ref expr)
{
    sm_ref typ;
    if (expr->node_type == cod_field_ref) {
	return is_array(expr->node.field_ref.sm_field_ref);
    }
    if (expr->node_type == cod_identifier) {
	return is_array(expr->node.identifier.sm_declaration);
    }
    if (expr->node_type == cod_declaration) {
	sm_ref ctype = expr->node.declaration.sm_complex_type;
	if ((ctype != NULL) && (ctype->node_type == cod_array_type_decl)) {
	    return 1;
	}
    }
    typ = get_complex_type(NULL, expr);
    if (typ == NULL) return 0;

    if (typ->node_type == cod_array_type_decl) {
	return 1;
    }

    if (typ->node_type == cod_reference_type_decl) {
	sm_ref ctype = 
	    typ->node.reference_type_decl.sm_complex_referenced_type;
	if (ctype == NULL) return 0;
	if (ctype->node_type == cod_array_type_decl) {
	    return 1;
	}
    }
    return 0;
}
    
static sm_ref
get_containing_structure(sm_ref expr)
{
    switch(expr->node_type) {
    case cod_element_ref:
	return get_containing_structure(expr->node.element_ref.array_ref);
    case cod_field_ref:
	return expr->node.field_ref.struct_ref;
    default:
	return NULL;
    }
}

    
static void
add_field_list(int *format_count_p, FMStructDescList *format_list_p, sm_ref typ)
{
    sm_list fields =  typ->node.struct_type_decl.fields;
    FMFieldList field_list = malloc(sizeof(field_list[0]) * 2);
    int field_count = 0;
    int my_format_num = (*format_count_p)++;
    *format_list_p = realloc(*format_list_p, sizeof(*format_list_p[0]) * (*format_count_p + 1));
    while(fields != NULL) {
	sm_ref typ = fields->node->node.field.sm_complex_type;
	field_list = realloc(field_list, (sizeof(field_list[0]) * (field_count +2)));
	field_list[field_count].field_name = strdup(fields->node->node.field.name);
	field_list[field_count].field_type = strdup(fields->node->node.field.string_type);
	field_list[field_count].field_size = fields->node->node.field.cg_size;
	field_list[field_count].field_offset = fields->node->node.field.cg_offset;
	while((typ != NULL) && ((typ->node_type == cod_reference_type_decl) || (typ->node_type == cod_declaration)
				|| (typ->node_type == cod_array_type_decl))) {
	    if (typ->node_type == cod_reference_type_decl) {
		typ = typ->node.reference_type_decl.sm_complex_referenced_type;
	    } else if (typ->node_type == cod_array_type_decl) {
		typ = typ->node.array_type_decl.sm_complex_element_type;
	    } else if (typ->node_type == cod_declaration) {
		typ = typ->node.declaration.sm_complex_type;
	    }
	}
	if ((typ != NULL) && (typ->node_type == cod_struct_type_decl)) {
	    add_field_list(format_count_p, format_list_p, typ);
	}
	field_count++;
	fields = fields->next;
    }
    field_list[field_count].field_name = field_list[field_count].field_type = NULL;
    field_list[field_count].field_size = field_list[field_count].field_offset = 0;
    (*format_list_p)[my_format_num].format_name = strdup(typ->node.struct_type_decl.id);
    (*format_list_p)[my_format_num].field_list = field_list;
    (*format_list_p)[my_format_num].struct_size = typ->node.struct_type_decl.cg_size;
    (*format_list_p)[my_format_num].opt_info = NULL;
}

static FMStructDescList
build_format_list(cod_parse_context context, sm_ref expr)
{
    sm_ref typ = get_complex_type(context, expr);
    FMStructDescList formats = malloc(sizeof(formats[0]) * 2);
    int format_count = 0;
    if (typ == NULL) {
	cod_src_error(context, expr->node.field_ref.struct_ref, 
		      "Reference must be structured type", 
		      expr->node.field_ref.lx_field);
	return 0;
    }
    if (typ->node_type == cod_reference_type_decl) {
	typ = typ->node.reference_type_decl.sm_complex_referenced_type;
    }
    if (typ->node_type == cod_declaration) {
	typ = typ->node.declaration.sm_complex_type;
    }
    add_field_list(&format_count, &formats, typ);
    formats[format_count].format_name = NULL;
    formats[format_count].field_list = NULL;
    return formats;
}

static int is_left_hand_side(sm_ref expr);

static int semanticize_expr(cod_parse_context context, sm_ref expr, 
			    scope_ptr scope) 
{
    switch(expr->node_type) {
    case cod_identifier: {
	sm_ref tmp = resolve(expr->node.identifier.id, scope);
	if (tmp != NULL) {
            if (tmp->node_type == cod_constant) {
                srcpos old_srcpos = expr->node.identifier.lx_srcpos;
		free(expr->node.identifier.id);
                /* morph identifier into constant */
                expr->node_type = cod_constant;
                expr->node.constant.token = tmp->node.constant.token;
                expr->node.constant.const_val = strdup(tmp->node.constant.const_val);
		expr->node.constant.freeable_name = NULL;
                expr->node.constant.lx_srcpos = old_srcpos;
                return semanticize_expr(context, expr, scope);
            } else {
                expr->node.identifier.sm_declaration = tmp;
                return 1;
            }
	} else {
	    cod_src_error(context, expr,
			  "Undefined Symbol \"%s\"", 
			  expr->node.identifier.id);
	    return 0;
	}
    }
    case cod_comma_expression:
	if (!semanticize_expr(context, expr->node.comma_expression.left, scope)) 
	    return 0;
	if (!semanticize_expr(context, expr->node.comma_expression.right, scope)) 
	    return 0;
	return 1;
    case cod_cast: {
	int cg_type;
	sm_ref typ;
	if (expr->node.cast.expression &&
	    !semanticize_expr(context, expr->node.cast.expression, scope)) {
	    return 0;
	}

	typ = reduce_type_list(context, expr->node.cast.type_spec, &cg_type, 
			       scope, NULL, NULL);
	if ((cg_type == DILL_ERR) && (typ == NULL)) {
	    cod_src_error(context, expr, "Illegal cast");
	    return 0;
	}
	expr->node.cast.cg_type = cg_type;
	expr->node.cast.sm_complex_type = typ;
	return 1;
    }
    case cod_operator: {
	int ret = 1;
	if (expr->node.operator.left != NULL) {
	    if (!semanticize_expr(context, expr->node.operator.left, scope)) {
		ret = 0;
	    }
	}
	if (expr->node.operator.right != NULL) {
	    if (!semanticize_expr(context, expr->node.operator.right, scope)) {
		ret = 0;
	    }
	}
	if (ret == 0) return 0;
	if ((expr->node.operator.left != NULL) && 
	    (expr->node.operator.right != NULL)) {
	    expr->node.operator.operation_type = 
		determine_op_type(context, expr, 
				  expr->node.operator.left,
				  expr->node.operator.right);
	    if (expr->node.operator.operation_type == DILL_ERR) {
		return 0;
	    }
	} else if (expr->node.operator.right != NULL) {
	    expr->node.operator.operation_type = 
		determine_unary_type(context, expr, expr->node.operator.right);
	} else if (expr->node.operator.left != NULL) {
	    expr->node.operator.operation_type = 
		determine_unary_type(context, expr, expr->node.operator.left);
	}
	switch (expr->node.operator.op) {
	case op_leq: case op_lt: case op_geq: case op_gt: case op_neq: 
	case op_eq:  case op_log_neg: case op_log_or: case op_log_and:
	case op_sizeof:
	    expr->node.operator.result_type = DILL_I;
	    break;
	case op_address:
	    expr->node.operator.result_type = DILL_P;
	    if (expr->node.operator.right->node_type == cod_identifier) {
		sm_ref decl = expr->node.operator.right->node.identifier.sm_declaration;
		if (decl->node_type == cod_declaration) {
		    if (decl->node.declaration.param_num != -1) {
			if (decl->node.declaration.sm_complex_type == NULL) {
			    cod_src_error(context, expr, "Cannot take address of a pass-by-value parameter");
			    return 0;
			}
		    }
		    decl->node.declaration.addr_taken = 1;
		}
	    } else {
		if (!is_left_hand_side(expr->node.operator.right)) {
		    cod_src_error(context, expr, "Invalid operand to address operator");
		    return 0;
		}
	    }
	    break;
	case op_deref: {
	    sm_ref typ = get_complex_type(context, expr->node.operator.right);
	    if (!typ || ((typ->node_type != cod_reference_type_decl) && 
			 (typ->node_type != cod_array_type_decl))) {
		cod_src_error(context, expr, "Cannot dereference a non-reference type");
		return 0;
	    } else if (typ->node_type == cod_reference_type_decl) {
		expr->node.operator.result_type =
		    typ->node.reference_type_decl.cg_referenced_type;
	    } else if (typ->node_type == cod_array_type_decl) {
		expr->node.operator.result_type =
		    typ->node.array_type_decl.cg_element_type;
	    } else {
		assert(0);
	    }
	    break;
	}
	default:
	    /* Operator applied to pointer types? Check compatibility... */
	    if(expr->node.operator.operation_type == DILL_P) {
		
		switch(expr->node.operator.op) {
		case op_inc:
		case op_dec:
		    break;
		    
		case op_plus:
		{
		    sm_ref left  = expr->node.operator.left;
		    sm_ref right = expr->node.operator.right;
		    
		    sm_ref lcplx = NULL;
		    sm_ref rcplx = NULL;

		    if(!left) {
			cod_src_error(context, expr,
				      "Invalid operand to unary plus\n");
			return 0;
		    }
		    
		    /* Extract complex types, if any */
		    lcplx = get_complex_type(context, left);
		    rcplx = get_complex_type(context, right);
		 
		    /* Pointers do not add with complex types */
		    if(lcplx && rcplx) {
			cod_src_error(context, expr,
				      "Invalid operands to binary plus");
			return 0;
		    }

		    /*
		     * We're ok if we reach this, since that implies we have subtraction 
		     * between a pointer and an integral type. The suitability of the
		     * integral type has been checked in determine_op_type() already.
		     */
		}
		break;

		case op_minus:
		{
		    sm_ref left  = expr->node.operator.left;
		    sm_ref right = expr->node.operator.right;
		    
		    sm_ref lcplx = NULL;
		    sm_ref rcplx = NULL;

		    if(!left) {
			cod_src_error(context, expr,
				      "Invalid operand to unary minus\n");
			return 0;
		    }

		    /* Extract complex types, if any */
		    lcplx = get_complex_type(context, left);
		    rcplx = get_complex_type(context, right);

		    
		    /* If both are complex types... */
		    if(lcplx && rcplx) {

			/* If both are pointers... */
			if(((lcplx->node_type == cod_reference_type_decl) || (lcplx->node_type == cod_array_type_decl)) &&
			   ((rcplx->node_type == cod_reference_type_decl) || (rcplx->node_type == cod_array_type_decl))) {
			    /* Check if the argument pointers are compatible */
			    if(!are_compatible_ptrs(lcplx, rcplx)) {
				cod_src_error(context, expr,
					      "Incompatible pointer arguments to binary minus");
				return 0;
			    } else {
				/*
				 * Binary minus between two compatible pointers is allowed,
				 * but it produces an integral type, so we fix that here.
				 */
				expr->node.operator.result_type=DILL_L;
				
				/*
				 * NOTE how we return success directly from here and do not
				 * break from the switch through the line below setting the
				 * result_type from the operation_type... In this case this
				 * would cause problems. We want operation_type to stay DILL_P
				 * but the final result_type to be a DILL_L.
				 */
				return 1;
			    }
			} else {
			    /*
			     * Pointers and other complex types do not subtract.
			     * Arithmetic canno be done on non-pointer complex types.
			     */
			    cod_src_error(context, expr,
					  "Incompatible arguments to binary minus");
			    return 0;
			}
		    }

		    /*
		     * We're ok if we reach this, since that implies we have subtraction 
		     * between a pointer and an integral type. The suitability of the
		     * integral type has been checked in determine_op_type() already.
		     */
		}
		break;
		
		default:
		    cod_src_error(context, expr,
				  "Operator cannot be applied to pointer types!\n");
		    return 0;
		}
	    }
	    /*
	     * NOTE: If anything here changes, one (potentially) has to
	     * update the code above which deals with binary minus
	     * between two compatible pointers, changes the result to
	     * an integral type, and returns directly without going
	     * through this code (as all other cases do).
	     */
	    expr->node.operator.result_type=expr->node.operator.operation_type;
	}
	return ret;
    }
    case cod_constant:
	return 1;
    case cod_assignment_expression: {
	int ret = 1;
	if (!semanticize_expr(context, expr->node.assignment_expression.left, scope)) {
	    ret = 0;
	} else {
	    expr->node.assignment_expression.cg_type = 
		cod_sm_get_type(expr->node.assignment_expression.left);
	}
	if (expr->node.assignment_expression.left && is_const(expr->node.assignment_expression.left)) {
	    cod_src_error(context, expr->node.assignment_expression.left, "Invalid assignment, left side is const");
	    ret = 0;
	}
	if (!semanticize_expr(context, expr->node.assignment_expression.right, scope)){
	    ret = 0;
	} else {
	    int right_type = 
		cod_sm_get_type(expr->node.assignment_expression.right);
	    if ((right_type == DILL_P) && 
		(is_string(expr->node.assignment_expression.right))) {
		if (expr->node.assignment_expression.cg_type != DILL_P) {
		    cod_src_error(context, expr, "assignment mixes string and non-string types");
		    ret = 0;
		}
	    } else if ((right_type == DILL_B) || (right_type == DILL_ERR)) {
		cod_src_error(context, expr->node.assignment_expression.right, "Invalid assignment, right side must be simple type");
		ret = 0;
	    }
	}
	if ((expr->node.assignment_expression.cg_type == DILL_P) ||
	    (expr->node.assignment_expression.cg_type == DILL_ERR)) {
	    sm_ref ltyp = 
		get_complex_type(context, 
				 expr->node.assignment_expression.left);
//	    sm_ref rtyp = 
//		get_complex_type(context, 
//				 expr->node.assignment_expression.right);
	    if (ltyp == NULL) {
		if (!is_string(expr->node.assignment_expression.left)) {
		    cod_src_error(context, expr->node.assignment_expression.left, "Invalid assignment, left side must be simple, non-pointer type");
		    ret = 0;
		}
	    } else {
		if ((ltyp->node_type == cod_struct_type_decl) || (ltyp->node_type == cod_array_type_decl) || (ltyp->node_type == cod_enum_type_decl)) {
		    /* maybe OK */
		} else if (ltyp->node_type != cod_reference_type_decl) {
		    cod_src_error(context, expr->node.assignment_expression.left, "Invalid assignment, left side must be simple, non-pointer type");
		    ret = 0;
		}
	    }
	}
	if (ret == 1) {
	    ret = assignment_types_match(context,
					 expr->node.assignment_expression.left,
					 expr->node.assignment_expression.right, 
					 /* strict */ (expr->node.assignment_expression.op == op_eq));
	}
	return ret;
    }
    case cod_field_ref: {
	sm_ref typ;
	sm_list fields;
	if (!semanticize_expr(context, expr->node.field_ref.struct_ref, scope)) {
	    return 0;
	}
	typ = get_complex_type(context, expr->node.field_ref.struct_ref);
	if (typ == NULL) {
	    cod_src_error(context, expr->node.field_ref.struct_ref, 
			  "Reference must be structured type", 
			  expr->node.field_ref.lx_field);
	    return 0;
	}
	if (typ->node_type == cod_reference_type_decl) {
	    typ = typ->node.reference_type_decl.sm_complex_referenced_type;
	}
	if (typ->node_type == cod_declaration) {
	    typ = typ->node.declaration.sm_complex_type;
	}
	fields =  typ->node.struct_type_decl.fields;
	while(fields != NULL) {
	    if (strcmp(expr->node.field_ref.lx_field,
		       fields->node->node.field.name) == 0) {
		break;
	    }
	    fields = fields->next;
	}
	if (fields == NULL) {
	    cod_src_error(context, expr, 
			  "Unknown field reference, \"%s\".",
			  expr->node.field_ref.lx_field);
	    return 0;
	}
	expr->node.field_ref.sm_field_ref = fields->node;
	return 1;
    }
    case cod_element_ref: {
	if (semanticize_expr(context, expr->node.element_ref.array_ref, scope)) {
	    int cg_type;
	    sm_ref arr = get_complex_type(NULL, expr->node.element_ref.array_ref);
	    if (is_string(expr->node.element_ref.array_ref)) {
		expr->node.element_ref.this_index_dimension = 0;
		expr->node.element_ref.sm_complex_element_type = NULL;
		expr->node.element_ref.cg_element_type = DILL_C;
		expr->node.element_ref.sm_containing_structure_ref =
		    get_containing_structure(expr->node.element_ref.array_ref);
	    } else if (is_array(expr->node.element_ref.array_ref)) {
		if (arr->node_type == cod_reference_type_decl) {
		    arr = arr->node.reference_type_decl.sm_complex_referenced_type;
		}
		if (expr->node.element_ref.array_ref->node_type != cod_element_ref) {
		    /* bottom level of recursion, we're the left-most array index */
		    expr->node.element_ref.this_index_dimension = 0;
		} else {
		    sm_ref subindex = expr->node.element_ref.array_ref;
		    expr->node.element_ref.this_index_dimension = subindex->node.element_ref.this_index_dimension + 1;
		}
		expr->node.element_ref.sm_complex_element_type = 
		    arr->node.array_type_decl.sm_complex_element_type;

		expr->node.element_ref.cg_element_type = 
		    arr->node.array_type_decl.cg_element_type;
		expr->node.element_ref.sm_containing_structure_ref =
		    get_containing_structure(expr->node.element_ref.array_ref);
	    } else if (arr && (arr->node_type == cod_reference_type_decl)) {
		expr->node.element_ref.sm_complex_element_type =
		    arr->node.reference_type_decl.sm_complex_referenced_type;
		expr->node.element_ref.cg_element_type = 
		    arr->node.reference_type_decl.cg_referenced_type;
		expr->node.element_ref.sm_containing_structure_ref =
		    get_containing_structure(expr->node.element_ref.array_ref);
	    } else {
		cod_src_error(context, expr, "Indexed element must be array, string or reference type.");
		return 0;
	    }

	    if (!semanticize_expr(context, expr->node.element_ref.expression, scope)) {
		return 0;
	    }
	    
	    cg_type = cod_sm_get_type(expr->node.element_ref.expression);
	    switch(cg_type) {
	    case DILL_C:
	    case DILL_UC:
	    case DILL_S:
	    case DILL_US:
	    case DILL_I: 
	    case DILL_U: 
	    case DILL_L:
	    case DILL_UL:
		return 1;
		break;
	    }
	    cod_src_error(context, expr, 
			  "Index for element reference must be integer type");
	    return 0;
	}
	return 0;
    }
    case cod_subroutine_call: {
	sm_ref func_ref = expr->node.subroutine_call.sm_func_ref;
	char *id;
	sm_ref tmp;
	sm_list args;
	sm_list formals, tmp_formals, tmp_args;
	int ret = 1;
	if (func_ref->node_type == cod_identifier) {
	    id = func_ref->node.identifier.id;
	} else {
	    id = func_ref->node.declaration.id;
	}
	tmp = resolve(id, scope);
	args = expr->node.subroutine_call.arguments;
	int done;
	if (tmp != NULL) {
	    if ((tmp->node_type != cod_declaration) ||
		!tmp->node.declaration.is_subroutine) {
		cod_src_error(context, expr, 
			      "Identifier is not subroutine \"%s\".", 
			func_ref->node.identifier.id);
		return 0;
	    }
	    free(func_ref->node.identifier.id);
	    free(func_ref);
	    expr->node.subroutine_call.sm_func_ref = func_ref = tmp;
	    formals = func_ref->node.declaration.params;
	} else {
	    cod_src_error(context, func_ref, "Undefined Subroutine \"%s\".", 
			  func_ref->node.identifier.id);

	    return 0;
	}
	tmp_formals = formals;
	tmp_args = args;
	sm_list *last_arg_p = &args;
	/* add closure args if required */
	while (tmp_formals != NULL) {
	    sm_ref formal = tmp_formals->node;
	    if (formal && (formal->node.declaration.sm_complex_type != NULL)) {
		sm_ref ct = formal->node.declaration.sm_complex_type;
		if ((ct->node_type == cod_reference_type_decl) &&
		    (strcmp(ct->node.reference_type_decl.name, "cod_closure_context") == 0)) {
		    sm_list new_arg = malloc(sizeof(struct list_struct));
		    char tmp[30];
		    new_arg->next = tmp_args;
		    if (func_ref->node.declaration.closure_id == NULL) {
			strcpy(tmp, "0");
		    } else {
			sprintf(tmp, "%p", func_ref->node.declaration.closure_id);
			if (strncmp(tmp, "0x", 2) != 0) {
			    sprintf(tmp, "0x%p", func_ref->node.declaration.closure_id);
			}
		    }

		    new_arg->node = cod_new_constant();
		    new_arg->node->node.constant.token = integer_constant;
		    new_arg->node->node.constant.const_val = strdup(tmp);
		    *last_arg_p = new_arg;
		    tmp_args = new_arg;
		} else if ((ct->node_type == cod_reference_type_decl) &&
			   (strcmp(ct->node.reference_type_decl.name, "cod_exec_context") == 0)) {
		    tmp_formals = tmp_formals->next;
		    continue;
		}

	    }
	    tmp_formals = tmp_formals->next;
	    if (tmp_args) {
		last_arg_p = &tmp_args->next;
		tmp_args = tmp_args->next;
	    }
	}
	/* must do this assigment, in case things changed from the loop above */
	expr->node.subroutine_call.arguments = args;

	done = 0;
	while (!done) {
	    sm_ref arg = NULL;
	    sm_ref formal = NULL;
	    if (formals != NULL) {
		formal = formals->node;
	    }
	    if (args != NULL) {
		arg = args->node;
	    }
	    if (formal && (formal->node.declaration.sm_complex_type != NULL)) {
		sm_ref ct = formal->node.declaration.sm_complex_type;
		if ((ct->node_type == cod_reference_type_decl) &&
		    (ct->node.reference_type_decl.name != NULL)) {
		    if (strcmp(ct->node.reference_type_decl.name, "cod_exec_context") == 0) {
                        if (context->has_exec_context == 0) {
                            cod_src_error(context, arg, "Calling subroutine has no cod_exec_context");
                            return 0;
                        }
                        /* swallow next formal, we'll fill that in ourselves */
                        formals = formals->next;
                        continue;
                    }
                }
	    }
	    if ((args == NULL) && (formals != NULL)) {
		if (strcmp(formal->node.declaration.id, "...") != 0) {
		    cod_src_error(context, arg, "Too few arguments to function");
		    ret = 0;
		}
	    }
	    if (args == NULL) {
		done++;
		continue;
	    }
	    if (!semanticize_expr(context, arg, scope) ) {
		args = args->next;
		continue;
	    }
	    if (formal == NULL) {
		cod_src_error(context, arg, "Too many arguments to subroutine");
		ret = 0;
		return ret;
	    }
	    if (strcmp(formal->node.declaration.id, "...") != 0) {
		/* we've got a real formal to check against */
		/* do some checking... */
		int mismatch = 0;
		switch (cod_sm_get_type(arg)) {
		case DILL_D: case DILL_F:
		    if (formal->node.declaration.cg_type >= DILL_V) {
			mismatch++;
		    }
		    break;
		case DILL_I: case DILL_U:
		case DILL_L: case DILL_UL:
		    if (formal->node.declaration.cg_type == DILL_P) {
			sm_ref ct = formal->node.declaration.sm_complex_type;
			if (!ct || 
			    (ct->node_type != cod_reference_type_decl) ||
			    ((strcmp(ct->node.reference_type_decl.name, "cod_type_spec") != 0) &&
			     (strcmp(ct->node.reference_type_decl.name, "cod_closure_context") != 0))) {
			    if ((arg->node_type != cod_constant) || 
				(arg->node.constant.token != integer_constant)) {
				mismatch++;
			    } else {
				int tmp = -1;
				sscanf(arg->node.constant.const_val, "%d", &tmp);
				/* zero is an acceptable pointer */
				if (tmp != 0) {
				    mismatch++;
				}
			    }
			}
		    }
		    break;
		case DILL_P:
		    if (formal->node.declaration.cg_type != DILL_P) {
			if (!(formal->node.declaration.sm_complex_type &&
			      (formal->node.declaration.sm_complex_type->node_type ==
			       cod_reference_type_decl))) {
			    mismatch++;
			}
		    }
		    break;
		}

		if (mismatch) {
		    cod_src_error(context, arg, 
				  "Type mismatch, parameter \"%s\".",
			    formal->node.declaration.id);
		    ret = 0;
		}
	    }
	    if ((formals != NULL) &&
		(strcmp(formal->node.declaration.id, "...") != 0)) {
		formals = formals->next;
		formal = NULL;
		if (formals != NULL) formal = formals->node;
	    }
	    /* look ahead to next formal and insert an arg if it's cod_type_spec */
	    if (formal &&
		formal->node.declaration.sm_complex_type != NULL) {
		sm_ref ct = formal->node.declaration.sm_complex_type;
		if ((ct->node_type == cod_reference_type_decl) &&
		    (strcmp(ct->node.reference_type_decl.name, "cod_type_spec")
		     == 0)) {
		    /* swallow next formal, we'll fill that in ourselves */
		    sm_list tmp_args = malloc(sizeof(struct list_struct));
		    FMStructDescList list = build_format_list(context,  arg);
		    char tmp[30];
		    sprintf(&tmp[0], "0x%p", list);
		    tmp_args->node = cod_new_constant();
		    tmp_args->node->node.constant.token = integer_constant;
		    tmp_args->node->node.constant.const_val = strdup(tmp);
		    tmp_args->next = args->next;
		    args->next = tmp_args;
		}
	    }
	    args = args->next;
	    if ((args == NULL) && (formals != NULL)) {
		if (strcmp(formal->node.declaration.id, "...") != 0) {
		    cod_src_error(context, arg, "Too few arguments to function");
		    ret = 0;
		}
	    }
	}
	return ret;
    }
    case cod_conditional_operator: {
	int ret = 1;
	if (expr->node.conditional_operator.condition != NULL) {
	    if (!semanticize_expr(context, expr->node.conditional_operator.condition, scope)) {
		ret = 0;
	    }
	}
	if (expr->node.conditional_operator.e1 != NULL) {
	    if (!semanticize_expr(context, expr->node.conditional_operator.e1, scope)) {
		ret = 0;
	    }
	}
	if (expr->node.conditional_operator.e2 != NULL) {
	    if (!semanticize_expr(context, expr->node.conditional_operator.e2, scope)) {
		ret = 0;
	    }
	}
	expr->node.conditional_operator.result_type = 
	    determine_unary_type(context, expr, expr->node.conditional_operator.e1);
	return ret;
    }
    case cod_initializer_list: {
	sm_list items = expr->node.initializer_list.initializers;
	int ret = 1;
	while (items) {
	    if (!semanticize_expr(context, items->node, scope)) ret = 0;
	    items = items->next;
	}
	return ret;
    }
    case cod_initializer: {
	if (!semanticize_expr(context, expr->node.initializer.initializer, scope))
	    return 0;
	if (expr->node.initializer.designation) {
//	    if (!semanticize_expr(context, expr->node.initializer.designation, scope)) return 0;
	}
	return 1;
    }
    default:
	fprintf(stderr, "Unknown case in semanticize_expression\n");
	cod_print(expr);
    }
    return 0;
}

static int
is_left_hand_side(sm_ref expr)
{
    switch(expr->node_type) {
    case cod_identifier:
	return 1;
    case cod_operator:
	return 0;
    case cod_cast:
	return is_left_hand_side(expr->node.cast.expression);
    case cod_assignment_expression:
	return expr->node.assignment_expression.cg_type;
    case cod_declaration:
	return 1;
    case cod_constant:
	return 0;
    case cod_field_ref:
	return 1;
    case cod_element_ref:
	return 1;
    case cod_subroutine_call:
	return 0;
    default:
	fprintf(stderr, "Unknown case in is_left_hand_side()\n");
	cod_print(expr);
	assert(0);
    }
    return 0;
}

int
type_of_int_const_string(char *val)
{
/*
For decimal, it is the first type the value can fit in: int, long, long long

For hexadecimal, it is the first type the value can fit in: int, unsigned int, long, 
unsigned long, long long, unsigned long long
*/

    long i;
    int len = (int)strlen(val);
    int hex = 0;
    int specified_unsgned = 0, specified_lng = 0;
    if (val[0] == '0') {
	/* hex or octal */
	hex++;
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
    switch(val[len-1]) {
    case 'U':
    case 'u':
	specified_unsgned++;
	break;
    case 'l':
    case 'L':
	specified_lng++;
	break;
    }
    if (len > 2) 
	switch(val[len-2]) {
	case 'U':
	case 'u':
	    specified_unsgned++;
	    break;
	case 'l':
	case 'L':
	    specified_lng++;
	    break;
	}
    if (len > 3) 
	switch(val[len-3]) {
	case 'U':
	case 'u':
	    specified_unsgned++;
	    break;
	case 'l':
	case 'L':
	    specified_lng++;
	    break;
	}
    if (specified_lng == 0) {
	/* unspecified */
	if (hex) {
	    if (i == (int)i) return DILL_I;
	    if (i == (unsigned)i) return DILL_U;
	    if (i == (long)i) return DILL_L;
	    if (i == (unsigned)i) return DILL_UL;
	    return DILL_UL;  /* don't do long long now */
	} else {
	    if (i == (int)i) return DILL_I;
	    if (i == (long)i) return DILL_L;
	    return DILL_L; /* don't do long long now */
	}
    }
    /* must have specified long */
    if (specified_unsgned) {
	return DILL_UL;
    } else {
	return DILL_L;
    }
}

extern int
cod_sm_get_type(sm_ref node)
{
    switch(node->node_type) {
    case cod_identifier:
	if (node->node.identifier.sm_declaration != NULL) {
	    return cod_sm_get_type(node->node.identifier.sm_declaration);
	}
	return node->node.identifier.cg_type;
    case cod_enumerator:
	return DILL_I;
    case cod_operator:
	return node->node.operator.result_type;
    case cod_conditional_operator:
	return node->node.conditional_operator.result_type;
    case cod_cast:
	return node->node.cast.cg_type;
    case cod_assignment_expression:
	return node->node.assignment_expression.cg_type;
    case cod_declaration:
	if (is_array(node)) {
	    return DILL_P;
	} else {
	    return node->node.declaration.cg_type;
	}
    case cod_constant:
	/* need to handle bigger constants */
	if (node->node.constant.token == string_constant) {
	    return DILL_P;
	} else if (node->node.constant.token == floating_constant) {
	    return DILL_D;
	} else if (node->node.constant.token == character_constant) {
	    return DILL_C;
	} else {
	    return type_of_int_const_string(node->node.constant.const_val);
	}
    case cod_field_ref:
	return cod_sm_get_type(node->node.field_ref.sm_field_ref);
    case cod_element_ref:
	return node->node.element_ref.cg_element_type;
    case cod_field:
	if (is_array(node)) {
	    return DILL_P;
	} else {
	    return node->node.field.cg_type;
	}
    case cod_initializer_list:
	return DILL_ERR;
    case cod_subroutine_call:
	return cod_sm_get_type(node->node.subroutine_call.sm_func_ref);
    case cod_comma_expression:
	return cod_sm_get_type(node->node.comma_expression.right);
    default:
	fprintf(stderr, "Unknown case in cod_sm_get_type()\n");
	cod_print(node);
    }
    return DILL_ERR;
}

extern int
are_compatible_ptrs(sm_ref left, sm_ref right) {
    sm_ref lTyp = NULL, rTyp = NULL;
    int lcgTyp = -1, rcgTyp = -1;

    /* Sanity check */
    if (left->node_type == cod_reference_type_decl) {
	lTyp = left->node.reference_type_decl.sm_complex_referenced_type;
	lcgTyp = left->node.reference_type_decl.cg_referenced_type;
    } else if (left->node_type == cod_array_type_decl) {
	lTyp = left->node.array_type_decl.sm_complex_element_type;
	lcgTyp = left->node.array_type_decl.cg_element_type;
    } else {
	return 0;
    }
    if (right->node_type == cod_reference_type_decl) {
	rTyp = right->node.reference_type_decl.sm_complex_referenced_type;
	rcgTyp = right->node.reference_type_decl.cg_referenced_type;
    } else if (right->node_type == cod_array_type_decl) {
	rTyp = right->node.array_type_decl.sm_complex_element_type;
	rcgTyp = right->node.array_type_decl.cg_element_type;
    } else {
	return 0;
    }

    if(lTyp && rTyp) {
	/* Two complex referenced types */
	if(((lTyp->node_type == cod_reference_type_decl) || (lTyp->node_type == cod_array_type_decl)) &&
	   ((rTyp->node_type == cod_reference_type_decl) || (rTyp->node_type == cod_array_type_decl))) {
	    /* Recurse if both are pointers */
	    return are_compatible_ptrs(lTyp, rTyp);
	}
	return (lTyp == rTyp)?1:0;
    }

    if(!lTyp && !rTyp) {
	/* Two integral referenced types */
	return (rcgTyp == lcgTyp)?1:0;
    }

    /* Mix of a pointer to a complex type and a pointer to an integral type */
    return 0;
}

extern sm_ref
get_complex_type(cod_parse_context context, sm_ref node)
{
    if (!node) return NULL;
    switch(node->node_type) {
    case cod_array_type_decl:
    case cod_reference_type_decl:
    case cod_struct_type_decl:
    case cod_enum_type_decl:
	return node;
    case cod_subroutine_call:
	return get_complex_type(context,
				node->node.subroutine_call.sm_func_ref);
    case cod_identifier:
	return get_complex_type(context, 
				node->node.identifier.sm_declaration);
    case cod_element_ref:
	return node->node.element_ref.sm_complex_element_type;
    case cod_field:
	return node->node.field.sm_complex_type;
    case cod_declaration:
	return get_complex_type(context, node->node.declaration.sm_complex_type);
    case cod_field_ref:{
	sm_ref typ;
	sm_list fields;
	typ = get_complex_type(context, node->node.field_ref.struct_ref);
	if (typ->node_type == cod_reference_type_decl) {
	    typ = typ->node.reference_type_decl.sm_complex_referenced_type;
	}
	if (typ->node_type == cod_declaration) {
	    typ = typ->node.declaration.sm_complex_type;
	}
	fields =  typ->node.struct_type_decl.fields;
	while ((fields != NULL) && 
	       (strcmp(node->node.field_ref.lx_field,
		       fields->node->node.field.name) != 0)) {
	    fields = fields->next;
	}
	if (fields == NULL) {
	    cod_src_error(context, node, "Unknown field reference \"%s\".",
		    node->node.field_ref.lx_field);
	    return NULL;
	}
	return get_complex_type(context, fields->node->node.field.sm_complex_type);
    }
    case cod_conditional_operator:
	return NULL;
    case cod_constant:
	return NULL;
    case cod_operator:
	switch (node->node.operator.op) {
	case op_deref: {
	    sm_ref right = get_complex_type(NULL, node->node.operator.right);
	    if ((right != NULL) && 
		(right->node_type == cod_reference_type_decl)) {
		sm_ref typ = right->node.reference_type_decl.sm_complex_referenced_type;
		if (typ && (typ->node_type == cod_declaration)) {
		    return get_complex_type(context, typ);
		} else {
		    return typ;
		}
	    }
	    return NULL;
	}
	case op_plus: case op_minus: case op_inc: case op_dec: {
	    sm_ref right = NULL;
	    sm_ref left  = NULL;
	    if (node->node.operator.right)
		right = get_complex_type(NULL, node->node.operator.right);
	    if (node->node.operator.left)
		left = get_complex_type(NULL, node->node.operator.left);
	    if (right && (left == NULL)) return right;
	    if (left && (right == NULL)) return left;
	    if ((left == NULL) && (right == NULL)) return NULL;
	    /*
	     * GANEV: op_minus can be applied to two pointers,
	     * i.e. two complex types => both left _and_ right can be
	     * non-NULL, hence this code. This shouldn't happen in
	     * other cases... (I think).
	     */
	    if(node->node.operator.op == op_minus && right && left) {
		if(left->node_type == cod_reference_type_decl &&
		   right->node_type == cod_reference_type_decl) {
		    /* Ok, so it's op_minus between two pointers, then check compatibility */
		    if(are_compatible_ptrs(left, right)) {
			return left;
		    } else {
			cod_src_error(context, node, 
				      "Incompatible pointer args to binary minus");
			return NULL;
		    }
		}
	    }
	    cod_src_error(context, node, "Incompatible pointer arguments to operator");
	    return NULL;
	}
	default:
	    return NULL;
	}

    case cod_cast:
	return node->node.cast.sm_complex_type;
	break;
    case cod_initializer_list:
    case cod_enumerator:
	return NULL;
    case cod_assignment_expression:
	return get_complex_type(context, node->node.assignment_expression.left);
    default:
	fprintf(stderr, "Unknown case in get_complex_type()\n");
	cod_print(node);
    }
    return NULL;
}

static sm_ref
reduce_type_list(cod_parse_context context, sm_list type_list, int *cg_type,
		 scope_ptr scope, int*is_typedef, sm_ref *freeable_type) 
{
    sm_list orig_list = type_list;
    int short_appeared = 0;
    int long_appeared = 0;
    int long_long_appeared = 0;
    int int_appeared = 0;
    int double_appeared = 0;
    int float_appeared = 0;
    int char_appeared = 0;
    int signed_appeared = 0;
    int unsigned_appeared = 0;
    int void_appeared = 0;
    int string_appeared = 0;
    int spec_count = 0;
    int prefix_end = 0;
    int type_found = 0;
    sm_ref complex_return_type = NULL;;

    *cg_type = DILL_ERR;
    while ((type_list != NULL) && (prefix_end == 0)) {
	int typ = type_list->node->node.type_specifier.token;
	if ((type_list->node->node_type != cod_type_specifier) ||
	    (typ == STAR) || (typ == AT)) {
	    prefix_end = 1;
	} else {
	    spec_count++;
	    switch (typ) {
	    case INT:
		int_appeared++;
		break;
	    case LONG:
		long_appeared++;
		break;
	    case SHORT:
		short_appeared++;
		break;
	    case DOUBLE:
		double_appeared++;
		break;
	    case STRING:
		string_appeared++;
		break;
	    case VOID:
		void_appeared++;
		break;
	    case FLOAT:
		float_appeared++;
		break;
	    case CHAR:
		char_appeared++;
		break;
	    case SIGNED:
		signed_appeared++;
		break;
	    case UNSIGNED:
		unsigned_appeared++;
		break;
	    case TYPEDEF:
		if (is_typedef) (*is_typedef)++;
		spec_count--;
		break;
	    case STATIC:
		spec_count--;
		break;
	    case EXTERN_TOKEN:
		spec_count--;
		break;
	    case CONST:
		spec_count--;
		break;
	    default:
		printf("Unknown type\n");
	    }
	    type_list = type_list->next;
	}
    }
    if (spec_count == 0) {
	if (type_list == NULL) *cg_type = DILL_I;   /* default to int */
	goto finalize;
    }
    if (void_appeared && (spec_count > 1)) {
	cod_src_error(context, orig_list->node, 
		      "Void type may not appear with other specifiers");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (string_appeared && (spec_count > 1)) {
	cod_src_error(context, orig_list->node, 
		      "String type may not appear with other specifiers");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (void_appeared) {
	*cg_type = DILL_V;
	goto finalize;
    }
    if (string_appeared) {
	*cg_type = DILL_P;
	goto finalize;
    }
    if (short_appeared && long_appeared ) {
	cod_src_error(context, orig_list->node, 
		      "Only one of long or short permitted");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (short_appeared && (double_appeared + float_appeared)) {
	cod_src_error(context, orig_list->node, 
		      "Short may not be specified with double or float");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (double_appeared + float_appeared) {
	if (double_appeared + float_appeared + short_appeared + signed_appeared + unsigned_appeared + char_appeared + int_appeared > 1) {
	    cod_src_error(context, orig_list->node, "Bad type spec");
	    *cg_type = DILL_ERR;
	    return NULL;
	} else {
	    /* not handling LONG plus one of these */
	    if (double_appeared) {
		*cg_type = DILL_D;
		goto finalize;
	    } else {
		*cg_type = DILL_F;
		goto finalize;
	    }
	}
    }

    /* neither float or double appeared */
    if (long_appeared == 2) {
	long_long_appeared++;
	long_appeared = 0;
    }
    if (short_appeared + char_appeared + long_appeared + long_long_appeared >= 2) {
	cod_src_error(context, orig_list->node, 
		      "Only one integer size spec may be specified");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (unsigned_appeared + signed_appeared > 1) {
	cod_src_error(context, orig_list->node, "Bad type spec");
	*cg_type = DILL_ERR;
	return NULL;
    }
    if (unsigned_appeared) {
	if (char_appeared) {
	    *cg_type = DILL_UC;
	    goto finalize;
	} else if (short_appeared) {
	    *cg_type = DILL_US;
	    goto finalize;
	} else if (long_appeared || long_long_appeared) {
	    *cg_type = DILL_UL;
	    goto finalize;
	} else {
	    *cg_type = DILL_U;
	    goto finalize;
	}
    } else {
	if (char_appeared) {
	    *cg_type = DILL_C;
	    goto finalize;
	} else if (short_appeared) {
	    *cg_type = DILL_S;
	    goto finalize;
	} else if (long_appeared || long_long_appeared) {
	    *cg_type = DILL_L;
	    goto finalize;
	} else if (int_appeared) {
	    *cg_type = DILL_I;
	    goto finalize;
	}
    }
 finalize:
    if (type_list == NULL) {
	/* no error and no more to process */
	return NULL;  
    }
    if (*cg_type != DILL_ERR) {
	type_found++;
    }
    while (type_list != NULL) {
	sm_ref node = type_list->node;
	switch (node->node_type) {
	case cod_identifier: 
	{
	    if (type_found != 0) {
		cod_src_error(context, node, 
			      "Type identifier cannot follow prior identifiers");
		*cg_type = DILL_ERR;
		return NULL;
	    }
	    complex_return_type = find_complex_type(node, scope);
	    if ((complex_return_type != NULL)&&
		(complex_return_type->node_type == cod_declaration)) {
		if (complex_return_type->node.declaration.sm_complex_type) {
		    complex_return_type = complex_return_type->node.declaration.sm_complex_type;
		} else {
		    *cg_type = complex_return_type->node.declaration.cg_type;
		}
	    }
	    if ((complex_return_type != NULL)&&
		(complex_return_type->node_type == cod_reference_type_decl)) {
		*cg_type = DILL_P;
	    }
	    if ((complex_return_type != NULL)&&
		(complex_return_type->node_type == cod_struct_type_decl)) {
		*cg_type = DILL_B;
	    }
	    if ((complex_return_type != NULL)&&
		(complex_return_type->node_type == cod_enum_type_decl)) {
		*cg_type = DILL_I;
	    }
	    if ((complex_return_type == NULL) && node->node.identifier.id &&
		((strcmp(node->node.identifier.id, "cod_type_spec") == 0) ||
                 (strcmp(node->node.identifier.id, "cod_closure_context") == 0) ||
		 (strcmp(node->node.identifier.id, "cod_exec_context") == 0))) {
		/* special ECL type information for prior arg */
		sm_ref typ = cod_new_reference_type_decl();
		typ->node.reference_type_decl.name = strdup(node->node.identifier.id);
		if (strcmp(node->node.identifier.id, "cod_type_spec") == 0) {
		    *cg_type = DILL_P;
                } else if (strcmp(node->node.identifier.id, "cod_closure_context") == 0) {
                    *cg_type = DILL_P;
		} else {
		    context->has_exec_context = 1;
		    *cg_type = DILL_P;
		}
		typ->node.reference_type_decl.cg_referenced_type = *cg_type;
		typ->node.reference_type_decl.sm_complex_referenced_type =
		    complex_return_type;
		typ->node.reference_type_decl.kernel_ref = 0;
		complex_return_type = typ;
		if (*freeable_type) {
		    cod_rfree(*freeable_type);
		    *freeable_type = NULL;
		}
		*freeable_type = typ;
	    }
	    assert((complex_return_type != NULL) || (*cg_type != DILL_ERR));
	    type_found++;
	}
	break;
	case cod_type_specifier:
	    switch (node->node.type_specifier.token) {
	    case STAR:
	      {
		  if (node->node.type_specifier.created_type_decl == NULL) {
		      /* GSE create anon-type */
		      sm_ref typ = cod_new_reference_type_decl();
		      typ->node.reference_type_decl.name = gen_anon();
		      typ->node.reference_type_decl.cg_referenced_type = *cg_type;
		      *cg_type = DILL_P;
		      typ->node.reference_type_decl.sm_complex_referenced_type =
			  complex_return_type;
		      typ->node.reference_type_decl.kernel_ref = 0;
		      complex_return_type = typ;
		      node->node.type_specifier.created_type_decl = typ;
		  } else {
		      complex_return_type = node->node.type_specifier.created_type_decl;
		      *cg_type = DILL_P;
		  }
	      }
	      break;
	    case AT:
	      {
		  /* GSE create anon-type */
		  sm_ref typ = cod_new_reference_type_decl();
		  typ->node.reference_type_decl.name = gen_anon();
		  typ->node.reference_type_decl.cg_referenced_type = *cg_type;
		  *cg_type = DILL_P;
		  typ->node.reference_type_decl.sm_complex_referenced_type =
		      complex_return_type;
		  typ->node.reference_type_decl.kernel_ref = 1;
		  complex_return_type = typ;
	      }
	      break;
	    default:
		if (type_found != 0) {
		    cod_src_error(context, node, 
				  "Only '*', '@', and CONST can follow valid type");
		    *cg_type = DILL_ERR;
		    return NULL;
		}
	    }
	    break;
	case cod_struct_type_decl: {
	    if (node->node.struct_type_decl.fields != NULL) {
		semanticize_decl(context, node, scope);
		complex_return_type = node;
	    } else {
		complex_return_type = resolve(node->node.struct_type_decl.id, scope);
		if (complex_return_type == NULL) {
		    cod_src_error(context, node,
				  "Struct declaration not found");
		    return NULL;
		}
	    }
	    *cg_type = DILL_B;
	    break;
	}
	case cod_enum_type_decl: {
	    if (node->node.enum_type_decl.enums != NULL) {
		semanticize_decl(context, node, scope);
		complex_return_type = node;
	    } else {
		complex_return_type = resolve(node->node.enum_type_decl.id, scope);
		if (complex_return_type == NULL) {
		    cod_src_error(context, node,
				  "Enum declaration not found");
		    return NULL;
		}
	    }
	    *cg_type = DILL_I;
	    break;
	}
	default:
	    printf("Unexpected node in reduce_type_list\n");
	    return NULL;
	}
	type_list = type_list->next;
    }
    return complex_return_type;
}

static int 
assignment_types_match(cod_parse_context context, sm_ref left, sm_ref right, int strict)
{
    sm_ref left_smt, right_smt;
    int left_cgt, right_cgt;
    left_smt = get_complex_type(context, left);
    right_smt = get_complex_type(context, right);
    left_cgt = cod_sm_get_type(left);
    right_cgt = cod_sm_get_type(right);
    if ((left_smt == NULL) && (right_smt == NULL)) {
	/* just check cgts */
	/* don't mix DILL_P, DILL_B and anything else */
	switch (left_cgt) {
	case DILL_P: 
	    switch (right_cgt) {
	    case DILL_P:
	    case DILL_L:
	    case DILL_UL:
		return 1;
		break;
	    default:
		cod_src_error(context, left, "Trying to assign a pointer variable with a non-pointer value.");
		return 0;
	    }
	default:
	    switch (right_cgt) {
	    case DILL_P:
		cod_src_error(context, left, "Trying to assign pointer to an incompatible variable.");
		return 0;
	    default:
		return 1;
		break;
	    }
	}
    }
    if ((left_smt != NULL) && 
	((left_smt->node_type != cod_reference_type_decl) &&
	 (left_smt->node_type != cod_array_type_decl) &&
	 (left_smt->node_type != cod_struct_type_decl) &&
	 (left_smt->node_type != cod_enum_type_decl))) {
	if ((left_cgt == DILL_P) || (left_cgt == DILL_B)) {
	    cod_src_error(context, left, "Only pointer, array, struct or enum complex types allowed as LHS in assignment");
	    return 0;
	}
    }
    if ((right_smt != NULL) && 
	((right_smt->node_type != cod_reference_type_decl) &&
	 (right_smt->node_type != cod_array_type_decl) &&
	 (right_smt->node_type != cod_struct_type_decl) &&
	 (right_smt->node_type != cod_enum_type_decl))) {
	if ((right_cgt == DILL_P) || (right_cgt == DILL_B)) {
	    cod_src_error(context, right, "Only pointer, array, struct or enum complex types allowed as RHS in assignment");
	    return 0;
	}
    }
    if (left_smt && (left_smt->node_type == cod_reference_type_decl) &&
	(right_smt == NULL)) {

	switch(right_cgt) {
	case DILL_P:
	case DILL_L:
	case DILL_UL:
	    return 1;
	case DILL_I:
	case DILL_U:
	    if (!strict) return 1;
	    if ((right->node_type == cod_constant) &&
		(right->node.constant.token == integer_constant)) {
		int i = -1;
		sscanf(right->node.constant.const_val, "%d", &i);
		if (i== 0) return 1;
	    }
	    /* falling through */
	default:
	    cod_src_error(context, right, "Right hand side must be pointer type");
	    return 0;
	}
    }
    if (right_smt && (left_smt == NULL)) {
	switch(left_cgt) {
	case DILL_C:
	case DILL_UC:
	case DILL_S:
	case DILL_US:
	case DILL_I:
	case DILL_U:
	case DILL_L:
	case DILL_UL:
	case DILL_P:
	    /* GANEV: should we have a warning here? */
	    return 1;

	default:
	    cod_src_error(context, right, "Pointer converted without explicit cast");
	    return 0;
	}
    }
    return 1;	
}

static int semanticize_struct_type_node(cod_parse_context context, sm_ref decl, 
					scope_ptr scope);

static int semanticize_enum_type_node(cod_parse_context context, sm_ref decl, 
				      scope_ptr scope);

static int
is_constant_expr(sm_ref expr)
{
    switch(expr->node_type) {
    case cod_constant: {
	return 1;
	break;
    }
    case cod_identifier:
	if (!expr->node.identifier.sm_declaration) return 0;
	return is_constant_expr(expr->node.identifier.sm_declaration);
    case cod_declaration:
	if (!expr->node.declaration.const_var) return 0;
	return is_constant_expr(expr->node.declaration.init_value);
    case cod_operator: {
	if (expr->node.operator.left != NULL) {
	    if (!is_constant_expr(expr->node.operator.left)) return 0;
	}
	if (expr->node.operator.op == op_sizeof) {
	    return 1;
	}
	if (expr->node.operator.right != NULL) {
	    if (!is_constant_expr(expr->node.operator.right)) return 0;
	}
	switch(expr->node.operator.op) {
	case  op_modulus:
	case  op_not:
	case  op_plus:
	case  op_minus:
	case  op_leq:
	case  op_lt:
	case  op_geq:
	case  op_gt:
	case  op_eq:
	case  op_neq:
	case  op_log_or:
	case  op_arith_or:
	case  op_arith_xor:
	case  op_log_and:
	case  op_arith_and:
	case  op_mult:
	case  op_div:
	case op_log_neg:
	case op_left_shift:
	case op_right_shift:
	    return 1;
	    break;
	case op_deref:
	case op_address:
	case op_inc:
	case op_dec:
	case op_sizeof:
	    return 0;
	}
	return 1;
    }
    case cod_cast:
	return is_constant_expr(expr->node.cast.expression);
    case cod_assignment_expression:
    case cod_field_ref:
    case cod_element_ref:
    case cod_subroutine_call:
	return 0;
    default:
	assert(0);
    }
    return 0;
}

static int
possibly_set_sizes_to_match(cod_parse_context context, sm_ref decl, sm_ref init_value)
{
    sm_ref array_type = get_complex_type(context, decl);
    if (array_type->node.array_type_decl.size_expr) return 1;
    if ((init_value->node_type == cod_constant) && 
	(init_value->node.constant.token == string_constant)) {
	/* init value is a string, set the array size to strlen + 1 */
	sm_ref size_expr = cod_new_constant();
	char *str = malloc(40); /* plenty */
	size_expr->node.constant.token = integer_constant;
	sprintf(str, "%ld\n", (long) strlen(init_value->node.constant.const_val) + 1);
	size_expr->node.constant.const_val = str;
	array_type->node.array_type_decl.size_expr = size_expr;
	return 1;
    }

    if (is_array(decl)) {
	sm_list items;
	long size = 0;
	assert(init_value->node_type == cod_initializer_list);
	items = init_value->node.initializer_list.initializers;
	/* init value is a list of initializers, count */
	while (items) {
	    size++;
	    items = items->next;
	}
	sm_ref size_expr = cod_new_constant();
	char *str = malloc(40); /* plenty */
	size_expr->node.constant.token = integer_constant;
	sprintf(str, "%ld\n", size);
	size_expr->node.constant.const_val = str;
	array_type->node.array_type_decl.size_expr = size_expr;
	return 1;
    }
    printf("Decl is : \n"); cod_print(decl);
    printf("init_value is : \n"); cod_print(init_value);
    return 1;
}
static int semanticize_decl(cod_parse_context context, sm_ref decl, 
			    scope_ptr scope)
{
    switch(decl->node_type) {
    case cod_declaration: {
	sm_ref ctype;
	int is_block_type = 0;

	if (resolve_local(decl->node.declaration.id, scope) != NULL) {
	    if (resolve_local(decl->node.declaration.id, scope) != decl) {
		cod_src_error(context, decl, "Duplicate Symbol \"%s\"", 
			      decl->node.declaration.id);
		return 0;
	    } else {
		/* been here, done that */
		return 1;
	    }
	} else {
	    add_decl(decl->node.declaration.id, decl, scope);
	}
	if (scope->containing_scope == NULL) {
	    /* must be external variable */
	    void *extern_value = 
		resolve_extern(decl->node.declaration.id, scope);
	    if ((extern_value == NULL) && (context->alloc_globals)) {
		;
	    } else if ((extern_value == NULL) && (decl->node.declaration.cg_address == NULL) &&
		       (decl->node.declaration.const_var == 0)) {
		cod_src_error(context, decl, 
			      "External symbol lacking address \"%s\"", 
			decl->node.declaration.id);
		return  0;
	    }
	    if (extern_value) {
		decl->node.declaration.cg_address = extern_value;
	    }
	    decl->node.declaration.is_extern = 1;
	}
	if (decl->node.declaration.type_spec != NULL) {
	    sm_list l = decl->node.declaration.type_spec;
	    if ((l->node->node_type == cod_type_specifier) && 
		((l->node->node.type_specifier.token == STATIC) ||
		 (l->node->node.type_specifier.token == CONST))) {
		if ((l->node->node.type_specifier.token == STATIC) &&
		    !decl->node.declaration.is_subroutine) {
		    decl->node.declaration.static_var = 1;
		} else {
		    decl->node.declaration.const_var = 1;
		}
		decl->node.declaration.type_spec = l->next;
		free(l->node);
		free(l);
	    }
	}
	if (decl->node.declaration.static_var) {
	    if (decl->node.declaration.init_value != NULL) {
		sm_ref const_val = decl->node.declaration.init_value;
		if (const_val->node_type == cod_initializer_list) {
		    sm_list items = const_val->node.initializer_list.initializers;
		    /* init value is a list of initializers, count */
		    while (items) {
			sm_ref sub_init = items->node->node.initializer.initializer;
			if (!is_constant_expr(sub_init)) {
			    cod_src_error(context, sub_init, 
					  "Static initializer not constant. Variable \"%s\"",
					  decl->node.declaration.id);
			    return 0;
			}
			items = items->next;
		    }
		    
		} else if (!is_constant_expr(const_val)) {
		    cod_src_error(context, const_val, 
				  "Static initializer not constant. Variable \"%s\"",
			decl->node.declaration.id);
		    return 0;
		}
	    }
	}
	if ((decl->node.declaration.sm_complex_type != NULL) &&
	    (decl->node.declaration.param_num != -1)) {
	    /* complex type + param, must be pass by reference */
	    sm_ref type = decl->node.declaration.sm_complex_type;
	    decl->node.declaration.cg_type = DILL_P;
	    if (type->node_type == cod_array_type_decl) {
		int ret = semanticize_array_type_node(context, type, 
						      scope);
		if (ret == 0) return ret;
	    }
	}
	/* some array decls have sm_complex_type set already */
	if (decl->node.declaration.sm_complex_type == NULL) {
	    sm_ref typ = NULL;
	    int cg_type = DILL_I;
	    if (decl->node.declaration.type_spec != NULL) {
		int type_def = 0;
		typ = reduce_type_list(context, decl->node.declaration.type_spec,
				       &cg_type, scope, &type_def, &decl->node.declaration.freeable_complex_type);
		if (type_def) {
		    decl->node.declaration.is_typedef = 1;
		}
	    } else {
		sm_ref arr = decl->node.declaration.sm_complex_type;
		if ((arr != NULL) && 
		    (arr->node_type == cod_array_type_decl)) {
		    typ = reduce_type_list(context, 
					    arr->node.array_type_decl.type_spec, 
					   &cg_type, scope, NULL, &decl->node.declaration.freeable_complex_type);
		} else if ((arr != NULL) && (arr->node_type == cod_enum_type_decl)) {
		    cg_type = DILL_I;
		}
	    }
	    if ((typ == NULL) && (cg_type == DILL_ERR)) return 0;
	    decl->node.declaration.cg_type = cg_type;
	    decl->node.declaration.sm_complex_type = typ;
	}
	ctype = decl->node.declaration.sm_complex_type;
	if ((ctype != NULL) && ((ctype->node_type == cod_array_type_decl) || (ctype->node_type == cod_struct_type_decl))) {
	    is_block_type = 1;
	}
	if (decl->node.declaration.init_value != NULL) {
	    int ret;
	    ret = semanticize_expr(context, decl->node.declaration.init_value, 
				   scope);
	    if (ret == 0) return ret;
	    if (is_array(decl)) {
		ret = possibly_set_sizes_to_match(context, decl, decl->node.declaration.init_value);
	    }
	    ret = assignment_types_match(context, decl, 
					 decl->node.declaration.init_value, 1);
	    return ret;
	}
	if (decl->node.declaration.is_subroutine) {
	    int ret;
	    int param_count = 0;
	    sm_list params = decl->node.declaration.params;
	    scope_ptr sub_scope = push_scope(scope);
	    ret = semanticize_decls_list(context,
					 decl->node.declaration.params, 
					 sub_scope);
	    decl->node.declaration.varidiac_subroutine_param_count = -1;
	    while(params) {
		sm_ref formal = params->node;
		while(formal->node_type == cod_array_type_decl) {
		    formal = formal->node.array_type_decl.element_ref;
		}
		if (strcmp(formal->node.declaration.id, "...") == 0) {
		    decl->node.declaration.varidiac_subroutine_param_count = param_count;		    
		}
		params = params->next;
		param_count++;
	    }
	    pop_scope(sub_scope);
	    return ret;
	}
	return 1;
	break;
    }
    case cod_struct_type_decl:
	return semanticize_struct_type_node(context, decl, scope);
	break;
    case cod_array_type_decl:
	if (decl->node.array_type_decl.type_spec != NULL) {
	    sm_list l = decl->node.array_type_decl.type_spec;
	    if ((l->node->node_type == cod_type_specifier) && 
		(l->node->node.type_specifier.token == STATIC)) {
		decl->node.array_type_decl.type_spec = l->next;
		decl->node.array_type_decl.element_ref->node.declaration.static_var = 1;
		free(l->node);
		free(l);
	    }
	}
	return semanticize_array_type_node(context, decl, scope);
	break;
    case cod_reference_type_decl:
	return semanticize_reference_type_node(context, decl, scope);
	break;
    case cod_constant:
        return 1;
        break;
    case cod_enum_type_decl:
	return semanticize_enum_type_node(context, decl, scope);
	break;
    default:
	printf("Unhandled case in semanticize decls_list\n");
	cod_print(decl);
    }
    return 0;
}

static int semanticize_statement(cod_parse_context context, sm_ref stmt, 
				 scope_ptr scope);
static int 
check_last_statement_return_list(cod_parse_context context, sm_list stmts);

static int 
check_last_statement_return(cod_parse_context context, sm_ref stmt)
{
    switch (stmt->node_type) {
    case cod_selection_statement:
	if (!check_last_statement_return(context, stmt->node.selection_statement.then_part)) return 0;
	if (stmt->node.selection_statement.else_part && 
	    !check_last_statement_return(context, stmt->node.selection_statement.else_part)) return 0;
	return 1;
    case cod_compound_statement: {
	sm_list list = stmt->node.compound_statement.statements;
	if (!list) list = stmt->node.compound_statement.decls;
	if (list) return check_last_statement_return_list(context, list);
	return 1;
    }
    case cod_return_statement:
	return 1;
    case cod_expression_statement:
	return check_last_statement_return(context, stmt->node.expression_statement.expression);
    case cod_label_statement:
	return check_last_statement_return(context, stmt->node.label_statement.statement);
    case cod_subroutine_call: {
	sm_ref func_ref = stmt->node.subroutine_call.sm_func_ref;
	char *id;
	if (func_ref->node_type == cod_identifier) {
	    id = func_ref->node.identifier.id;
	} else {
	    id = func_ref->node.declaration.id;
	}
	if (strcmp(id, "exit") == 0) return 1;
	if (strcmp(id, "abort") == 0) return 1;
	return 0;
    }
    default:
	return 0;
    }
}

static int 
check_last_statement_return_list(cod_parse_context context, sm_list stmts)
{
    sm_ref stmt;
    while (stmts != NULL) {
	stmt = stmts->node;
	stmts = stmts->next;
    }
    if (!stmt) return 0;
    return check_last_statement_return(context, stmt);
}

static int
semanticize_decls_stmts_list(cod_parse_context context, sm_list decls_stmts, scope_ptr scope)
{
    int ret = 1;
    while (decls_stmts != NULL) {
	sm_ref item = decls_stmts->node;
	switch(item->node_type) {
	case cod_declaration: 
	case cod_struct_type_decl:
	case cod_array_type_decl:
	case cod_reference_type_decl:
	case cod_enum_type_decl:
	case cod_constant:
	    if (!semanticize_decl(context, item, scope)) {
		ret = 0;
	    }
	    break;
	default: {
	    int t = semanticize_statement(context, item, scope);
	    if (!t) {
		ret = 0;
	    }
	}
	}
	decls_stmts = decls_stmts->next;
    }
    return ret;
}

static int
semanticize_decls_list(cod_parse_context context, sm_list decls, 
		       scope_ptr scope)
{
    int ret = 1;
    while (decls != NULL) {
	if (!semanticize_decl(context, decls->node, scope)) {
	    ret = 0;
	}
	decls = decls->next;
    }
    return ret;
}

static int
semanticize_selection_statement(cod_parse_context context, sm_ref selection,
				scope_ptr scope)
{
    int ret = 1;
    if (!semanticize_expr(context, 
			  selection->node.selection_statement.conditional,
			  scope)) {
	ret = 0;
    }
    if (!semanticize_statement(context,
			       selection->node.selection_statement.then_part,
			       scope)) {
	ret = 0;
    }
    if (selection->node.selection_statement.else_part) {
	if (!semanticize_statement(context,
				   selection->node.selection_statement.else_part,
				   scope)) {
	    ret = 0;
	}
    }
    return ret;
}

static int
semanticize_iteration_statement(cod_parse_context context, sm_ref iteration,
				scope_ptr scope)
{
    int ret = 1;
    if (iteration->node.iteration_statement.init_expr != NULL) {
	if (!semanticize_expr(context, 
			      iteration->node.iteration_statement.init_expr,
			      scope)) {
	    ret = 0;
	}
    }

    if (iteration->node.iteration_statement.test_expr != NULL) {
	if (!semanticize_expr(context, 
			      iteration->node.iteration_statement.test_expr,
			      scope)) {
	    ret = 0;
	}
    }

    if (iteration->node.iteration_statement.iter_expr != NULL) {
	if (!semanticize_expr(context, 
			      iteration->node.iteration_statement.iter_expr,
			      scope)) {
	    ret = 0;
	}
    }

    if (iteration->node.iteration_statement.statement != NULL) {
	scope_ptr sub_scope = push_scope_container(scope, iteration);
	if (!semanticize_statement(context,
				   iteration->node.iteration_statement.statement,
				   sub_scope)) {
	    ret = 0;
	}
	pop_scope(sub_scope);
    }
    if (iteration->node.iteration_statement.post_test_expr != NULL) {
	if (!semanticize_expr(context, 
			      iteration->node.iteration_statement.post_test_expr,
			      scope)) {
	    ret = 0;
	}
    }

    return ret;
}

static int 
semanticize_statement(cod_parse_context context, sm_ref stmt, 
		      scope_ptr scope)
{
    if (!stmt) return 1;
    switch (stmt->node_type) {
    case cod_selection_statement:
	return semanticize_selection_statement(context, stmt, scope);
    case cod_iteration_statement:
	return semanticize_iteration_statement(context, stmt, scope);
    case cod_expression_statement: {
	return semanticize_expr(context, 
				stmt->node.expression_statement.expression,
				scope);
    }	
    case cod_compound_statement:
	return semanticize_compound_statement(context, stmt, scope, 0);
    case cod_return_statement:{
	int expr_type;
	stmt->node.return_statement.cg_func_type = context->return_cg_type;
	if (stmt->node.return_statement.cg_func_type == DILL_V) {
	    if (stmt->node.return_statement.expression != NULL) {
		cod_src_error(context, stmt, 
			      "Return value supplied in subroutine declared to return VOID");
		return 0;
	    }
	} else {
	    if (stmt->node.return_statement.expression == NULL) {
		cod_src_error(context, stmt, 
			      "Return value missing in non-VOID subroutine");
		return 0;
	    }
	}	    
	if (stmt->node.return_statement.expression == NULL) return 1;
	if (!semanticize_expr(context, stmt->node.return_statement.expression,
			      scope)) return 0;
	expr_type = cod_sm_get_type(stmt->node.return_statement.expression);
	if (context->dont_coerce_return) {
	    int type_failure = 0;
	    switch (stmt->node.return_statement.cg_func_type) {
	    case DILL_C: case DILL_UC:  case DILL_S: case DILL_US: case DILL_I: case DILL_U: case DILL_L: case DILL_UL:
		if (expr_type > DILL_UL) type_failure++;
		break;
	    case DILL_F: case DILL_D:
		if ((expr_type != DILL_F) && (expr_type != DILL_D)) type_failure++;
		break;
	    }
	    if (type_failure) {
		cod_src_error(context, stmt, 
			      "Return value doesn't match procedure type declaration and now allowed to use coercion");
		return 0;
	    }
	}
	return 1;
    }
    case cod_label_statement:{
//	add_decl(stmt->node.label_statement.name, stmt, scope);
	return semanticize_statement(context, stmt->node.label_statement.statement, scope);
    }
    case cod_jump_statement:{
	if (stmt->node.jump_statement.goto_target != NULL) {
	    if (!stmt->node.jump_statement.sm_target_stmt) {
		cod_src_error(context, stmt, 
			      "Label \"%s\" not found.  Goto has no target.", stmt->node.jump_statement.goto_target);
		return 0;
	    }
	} else {
	    /* this is a continue or a break */
	    sm_ref tmp = find_containing_iterator(scope);
	    if (!tmp) {
		cod_src_error(context, stmt, 
			      "Continue or Break statement not contained inside an iterator.");
		return 0;
	    }
	    stmt->node.jump_statement.sm_target_stmt = tmp;
	}
	return 1;
	break;
    }
    default:
	printf("unhandled case in semanticize statement\n");
	return 1;
    }
    return 1;
}

typedef struct goto_semantic_state {
    int backward_jump;
    int passed_init_decl;
    int already_found;
} *goto_state;

static int semanticize_goto_l(cod_parse_context context, sm_ref this_goto, sm_list stmts, goto_state gs);

static int semanticize_goto(cod_parse_context context, sm_ref this_goto, sm_ref stmt, goto_state gs)
{
    int ret = 1;
    if (!stmt) return ret;
    switch (stmt->node_type) {
    case cod_declaration: 
	if (!gs->backward_jump && stmt->node.declaration.init_value) {
	    gs->passed_init_decl = 1;
	}
	break;
    case cod_struct_type_decl:
    case cod_array_type_decl:
    case cod_reference_type_decl:
    case cod_enum_type_decl:
    case cod_constant:
	/* no action for decls */
	break;
    case cod_selection_statement:
	ret &= semanticize_goto(context, this_goto, stmt->node.selection_statement.then_part, gs);
	if (stmt->node.selection_statement.else_part)
	    ret &= semanticize_goto(context, this_goto, stmt->node.selection_statement.else_part, gs);
	break;
    case cod_iteration_statement:
	ret &= semanticize_goto(context, this_goto, stmt->node.iteration_statement.statement, gs);
	break;
    case cod_compound_statement:
	ret &= semanticize_goto_l(context, this_goto, stmt->node.compound_statement.decls, gs);
	ret &= semanticize_goto_l(context, this_goto, stmt->node.compound_statement.statements, gs);
	break;
    case cod_return_statement:
	break;
    case cod_label_statement:
	if (strcmp(this_goto->node.jump_statement.goto_target, stmt->node.label_statement.name) == 0) {
	    /* found target */
	    if (!gs->backward_jump && gs->passed_init_decl) {
		cod_src_error(context, stmt, "Goto jumps over initialized declaration, illegal forward jump.");
		ret = 0;
	    } else if (gs->already_found) {
		cod_src_error(context, stmt, "Duplicate label \"%s\".", stmt->node.label_statement.name);
		ret = 0;
	    } else {
		this_goto->node.jump_statement.sm_target_stmt = stmt;
		gs->already_found = 1;
	    }
	}
	ret &= semanticize_goto(context, this_goto, stmt->node.label_statement.statement, gs);
	break;
    case cod_jump_statement:
	if (stmt == this_goto) {
	    gs->backward_jump = 0;
	}
	break;
    case cod_expression_statement:
	break;
    default:
	printf("unhandled case in semanticize goto\n");
	return 0;
    }
    return ret;
}

static int semanticize_goto_l(cod_parse_context context, sm_ref this_goto, sm_list stmts, goto_state gs)
{
    int saved_passed_init_decl = gs->passed_init_decl;
    int ret = 1;
    while(stmts) {
	ret &= semanticize_goto(context, this_goto, stmts->node, gs);
	stmts = stmts->next;
    }
    gs->passed_init_decl = saved_passed_init_decl;
    return ret;
}
    
static int
semanticize_gotos_list(cod_parse_context context, sm_list stmts, sm_list function_context)
{
    int ret = 1;
    while(stmts) {
	ret &= semanticize_gotos(context, stmts->node, function_context);
	stmts = stmts->next;
    }
    return ret;
}

static int
semanticize_gotos(cod_parse_context context, sm_ref stmt, sm_list function_context)
{
    /* 
     * recursive descent looking for goto's, followed by a recursive descent in 
     * the entire scope looking for their target 
     */
    int ret = 1;
    if (!stmt) return 1;
    switch (stmt->node_type) {
    case cod_declaration: 
    case cod_struct_type_decl:
    case cod_array_type_decl:
    case cod_reference_type_decl:
    case cod_enum_type_decl:
    case cod_constant:
	/* no action for most decls */
	break;
    case cod_selection_statement:
	ret &= semanticize_gotos(context, stmt->node.selection_statement.then_part, function_context);
	if (stmt->node.selection_statement.else_part)
	    ret &= semanticize_gotos(context, stmt->node.selection_statement.else_part, function_context);
	break;
    case cod_iteration_statement:
	ret &= semanticize_gotos(context, stmt->node.iteration_statement.statement, function_context);
	break;
    case cod_compound_statement:
	ret &= semanticize_gotos_list(context, stmt->node.compound_statement.decls, function_context);
	ret &= semanticize_gotos_list(context, stmt->node.compound_statement.statements, function_context);
	break;
    case cod_return_statement:
	break;
    case cod_label_statement:
	ret &= semanticize_gotos(context, stmt->node.label_statement.statement, function_context);
	break;
    case cod_jump_statement:
	if (stmt->node.jump_statement.goto_target != NULL) {
	    /* this is a goto */
	    struct goto_semantic_state gs;
	    gs.backward_jump = 1;
	    gs.passed_init_decl = 0;
	    gs.already_found = 0;
	    ret &= semanticize_goto_l(context, stmt, function_context, &gs);
	}
	break;
    case cod_expression_statement:
	break;
    default:
	printf("unhandled case in semanticize gotos\n");
	return 0;
    }
    return ret;
}

extern int
semanticize_compound_statement(cod_parse_context context, sm_ref compound, 
			       scope_ptr containing_scope, int require_last_return)
{
    int ret = 1;
    scope_ptr current_scope = push_scope(containing_scope);

    ret &= semanticize_decls_stmts_list(context,
					compound->node.compound_statement.decls,
					current_scope);
    ret &= semanticize_decls_stmts_list(context,
					compound->node.compound_statement.statements,
					current_scope);
    if (ret && require_last_return) {
	int tmp;
	sm_list list = compound->node.compound_statement.statements;
	if (!list) list = compound->node.compound_statement.decls;
	tmp = check_last_statement_return_list(context, compound->node.compound_statement.statements);
	if (!tmp) {
	    cod_src_error(context, NULL, 
			  "Control reaches end of non-void function.");
	}
	ret &= tmp;
    }
    pop_scope(current_scope);
    return ret;
}

extern sm_ref
cod_build_type_node(const char *name, FMFieldList field_list)
{
    sm_ref decl = cod_new_struct_type_decl();
    sm_list *end_ptr = &decl->node.struct_type_decl.fields;

    decl->node.struct_type_decl.id = strdup(name);
    while ((field_list != NULL) && (field_list->field_name != NULL)) {
	sm_list new_elem;
	new_elem = malloc(sizeof(*new_elem));
	new_elem->next = NULL;
	new_elem->node = cod_new_field();
	new_elem->node->node.field.name = strdup(field_list->field_name);
	new_elem->node->node.field.string_type = strdup(field_list->field_type);
	new_elem->node->node.field.cg_size = field_list->field_size;
	new_elem->node->node.field.cg_offset = field_list->field_offset;
	new_elem->node->node.field.cg_type = DILL_ERR;
	*end_ptr = new_elem;
	end_ptr = &new_elem->next;
	field_list++;
    }
    return decl;
}


extern void
cod_remove_defined_types(cod_parse_context context, int count)
{
    char **types = context->defined_types;
    while(types && types[count]) types[count++] = NULL;
}

void
cod_add_defined_type(char *id, cod_parse_context context)
{
    int count = 0;
    while(context->defined_types && context->defined_types[count]) count++;
    if (count == 0) {
	context->defined_types = malloc(sizeof(char*) * 2);
    } else {
	context->defined_types = realloc(context->defined_types,
					 (count+2)*sizeof(char*));
    }
    context->defined_types[count] = id;
    context->defined_types[count+1] = NULL;
    reset_types_table(context->defined_types, context->enumerated_constants);
		      
}

void
cod_add_enum_const(char *id, cod_parse_context context)
{
    int count = 0;
    while(context->enumerated_constants && context->enumerated_constants[count]) count++;
    if (count == 0) {
	context->enumerated_constants = malloc(sizeof(char*) * 2);
    } else {
	context->enumerated_constants = realloc(context->enumerated_constants,
					 (count+2)*sizeof(char*));
    }
    context->enumerated_constants[count] = id;
    context->enumerated_constants[count+1] = NULL;
    reset_types_table(context->defined_types, context->enumerated_constants);
}

extern void
cod_add_simple_struct_type(const char *name, FMFieldList field_list, 
		    cod_parse_context context)
{
    sm_ref node = cod_build_type_node(name, field_list);
    cod_add_decl_to_parse_context(name, node, context);
    cod_add_decl_to_scope((char*)name, node, context);
}

extern void
cod_add_struct_type(FMStructDescList format_list, 
		    cod_parse_context context)
{
    int count=0;
    while(format_list && format_list[count].format_name) {
	count++;
    }
    count = count-1;
    for ( ; count >= 0; count--) {
	cod_add_simple_struct_type(format_list[count].format_name,
				   format_list[count].field_list,
				   context);
    }
}

static int
str_to_data_type(char *str, int size)
{
    char *tmp = malloc(strlen(str) + 1);
    char *free_str = tmp;
    strcpy(tmp, str);
    str = tmp;			/* make a copy of str parameter */

    while (isspace((int)*str) || (*str == '*') || (*str == '(')) {	/* skip preceeding space */
	str++;
    }
    tmp = str + strlen(str) - 1;
    while (isspace((int)*tmp) || (*tmp == ')')) {  /* test trailing space */
	*tmp = 0;
	tmp--;
    }
    tmp = str;
    while (*tmp) {		/* map to lower case */
	*tmp = tolower(*tmp);
	tmp++;
    }
    if ((strcmp(str, "integer") == 0) || (strcmp(str, "enumeration") == 0)) {
	free(free_str);
	if (size == sizeof(intptr_t)) {
	    return DILL_L;
	} else if (size == sizeof(int)) {
	    return DILL_I;
	} else if (size == sizeof(short)) {
	    return DILL_S;
	} else if (size == sizeof(char)) {
	    return DILL_C;
	} else {
	    return DILL_L;
	}
    } else if (strcmp(str, "unsigned integer") == 0) {
	free(free_str);
	if (size == sizeof(intptr_t)) {
	    return DILL_UL;
	} else if (size == sizeof(int)) {
	    return DILL_U;
	} else if (size == sizeof(short)) {
	    return DILL_US;
	} else if (size == sizeof(char)) {
	    return DILL_UC;
	} else {
	    return DILL_UL;
	}
    } else if ((strcmp(str, "float") == 0) || (strcmp(str, "double") == 0)) {
	free(free_str);
	if (size == sizeof(double)) {
	    return DILL_D;
	} else if (size == sizeof(float)) {
	    return DILL_F;
	} else {
	    fprintf(stderr, "unsupported float size %d\n", size);
	    return DILL_D;
	}
    } else if (strcmp(str, "char") == 0) {
	free(free_str);
	assert(size == 1);
	return DILL_C;
    } else if (strcmp(str, "string") == 0) {
	free(free_str);
	return DILL_P;
    } else {
	free(free_str);
	return DILL_ERR;
    }
}

static int
array_str_to_data_type(char *str, int size)
{
    int ret_type;
    char field_type[1024];
    char *left_paren;
    if ((left_paren = strchr(str, '[')) == NULL) {
	ret_type = str_to_data_type(str, size);
    } else {
	char *tmp = str;
	int i = 0;
	for( ; tmp < left_paren; tmp++) {
	    field_type[i++] = *tmp;
	}
	field_type[i] = 0;
	ret_type = str_to_data_type(field_type, size);
    }
    return ret_type;
}

static sm_ref
build_subtype_nodes(cod_parse_context context, sm_ref decl, field* f, FMTypeDesc *desc,
		    int *err, scope_ptr scope, int *must_free_p)
{
    sm_ref ret = NULL;
    sm_ref subtype = NULL;
    int must_free_flag = 0;
    if (desc->next != NULL) {
	subtype = build_subtype_nodes(context, decl, f, desc->next, err, scope, &must_free_flag);
	if (*err != 0) {
	    printf("Subtype node failure\n");
	    return NULL;
	}
    }
    switch (desc->type) {
    case FMType_array: {
	sm_list fields = decl->node.struct_type_decl.fields;
	sm_ref cf;
	int i;
	ret = cod_new_array_type_decl();
	*must_free_p = 1;
	ret->node.array_type_decl.cg_static_size = desc->static_size;
	if (desc->static_size == 0) {
	    ret->node.array_type_decl.cg_static_size = -1;
	}
	ret->node.array_type_decl.cg_element_type = DILL_B;
	ret->node.array_type_decl.sm_complex_element_type = subtype;
	if (must_free_flag) {
	    if (ret->node.array_type_decl.freeable_complex_element_type) {
	        cod_rfree(ret->node.array_type_decl.freeable_complex_element_type);
	    }
	    ret->node.array_type_decl.freeable_complex_element_type = subtype;
	}
	if (subtype == NULL) {
	    ret->node.array_type_decl.cg_element_type = 
		array_str_to_data_type(f->string_type, f->cg_size);
	    ret->node.array_type_decl.cg_element_size = f->cg_size;
	    ret->node.array_type_decl.dimensions = malloc(sizeof(struct dimen_p));
	    ret->node.array_type_decl.dimensions->dimen_count = 1;
	} else {
	    if (subtype->node_type == cod_array_type_decl) {
		int sub_size = subtype->node.array_type_decl.cg_static_size;
		int sub_dimensions = subtype->node.array_type_decl.dimensions->dimen_count;
		if (sub_size == -1) {
		    /* element of *this* array has varying elements */
		    ret->node.array_type_decl.cg_element_size = -1;
		} else {
		    ret->node.array_type_decl.cg_element_size = 
			sub_size * subtype->node.array_type_decl.cg_element_size;;
		}
		
		ret->node.array_type_decl.dimensions = malloc(sizeof(struct dimen_p) + sub_dimensions * sizeof(dimen_s));
		ret->node.array_type_decl.dimensions->dimen_count = sub_dimensions+1;;
		memcpy(&ret->node.array_type_decl.dimensions->dimens[1], &subtype->node.array_type_decl.dimensions->dimens[0], sub_dimensions * sizeof(dimen_s));
	    } else {
		ret->node.array_type_decl.cg_element_size = f->cg_size;
		ret->node.array_type_decl.dimensions = malloc(sizeof(struct dimen_p));
		ret->node.array_type_decl.dimensions->dimen_count = 1;
		if (subtype->node_type == cod_reference_type_decl) {
		    ret->node.array_type_decl.cg_element_type = DILL_P;
		}
	    }
	}
	    
	if (ret->node.array_type_decl.cg_static_size != -1) {
	    ret->node.array_type_decl.sm_dynamic_size = NULL;
	    ret->node.array_type_decl.dimensions->dimens[0].static_size = ret->node.array_type_decl.cg_static_size;
	    ret->node.array_type_decl.dimensions->dimens[0].control_field = NULL;
	} else {
	    for (i=0; i < desc->control_field_index; i++) {
		fields = fields->next;
	    }
	    cf = fields->node;
	    switch (str_to_data_type(cf->node.field.string_type, 
				     (int)sizeof(int))) {
	    case DILL_C: case DILL_UC: case DILL_S: case DILL_US: 
	    case DILL_I: case DILL_U: case DILL_L: case DILL_UL:
		break;
	    default:
		cod_src_error(context, NULL, 
			      "Variable length control field \"%s\"not of integer type.", cf->node.field.string_type);
		*err = 1;
		return NULL;
		break;
	    }
	    ret->node.array_type_decl.sm_dynamic_size = cf;
	    ret->node.array_type_decl.dimensions->dimens[0].static_size = -1;
	    ret->node.array_type_decl.dimensions->dimens[0].control_field = cf;
	}
	break;
    }
    case FMType_pointer:
	ret = cod_new_reference_type_decl();
	*must_free_p = 1;
	ret->node.reference_type_decl.name = gen_anon();
	ret->node.reference_type_decl.cg_referenced_type = DILL_ERR;
	ret->node.reference_type_decl.sm_complex_referenced_type = subtype;
	if (must_free_flag) {
	    if (ret->node.reference_type_decl.freeable_complex_referenced_type) {
	        cod_rfree(ret->node.reference_type_decl.freeable_complex_referenced_type);
	    }
	    ret->node.reference_type_decl.freeable_complex_referenced_type = subtype;
	}
	ret->node.reference_type_decl.cg_referenced_size = -1;
	break;
    case FMType_subformat: {
	char *tmp_str = FMbase_type(f->string_type);
	ret = resolve(tmp_str, scope);
	free(tmp_str);
	if (ret == NULL) {
	    printf("Didn't find base type %s\n", tmp_str);
	    *err = 1;
	}
	break;
    }
    case FMType_simple:
    case FMType_string:
	ret = NULL;
	break;
    }
    return ret;

}

static void
build_type_nodes(cod_parse_context context, sm_ref decl, field* f, sm_list fields,
		 int cg_size, int cg_type, FMTypeDesc* desc, int *err, scope_ptr scope)
{
    int must_free_flag = 0;
    sm_ref complex_type = build_subtype_nodes(context, decl, f, desc, err, scope, &must_free_flag);
    f->sm_complex_type = complex_type;
    if (must_free_flag) {
        if (f->freeable_complex_type) {
	    cod_rfree(f->freeable_complex_type);
	}
	f->freeable_complex_type = complex_type;
    }
}

static int
semanticize_array_element_node(cod_parse_context context, sm_ref array, sm_ref super_type, sm_list base_type_spec, scope_ptr scope)
{
    if (array->node.array_type_decl.size_expr != NULL) {
	if (!is_constant_expr(array->node.array_type_decl.size_expr)) {
	    cod_src_error(context, array, 
			  "Array size expression must be constant.");
		return 0;
	}
	if (semanticize_expr(context,
			     array->node.array_type_decl.size_expr, scope) == 0) {
	    return 0;
	}
    } else {
	sm_ref element_ref = array->node.array_type_decl.element_ref;
	/* allow NULL array type sizes */
	if (element_ref->node_type != cod_declaration) {
	    cod_src_error(context, element_ref, 
			  "Null sizes only allowed in parameter contexts");
	    return 0;
	}
	    
    }
    dimen_p d = super_type->node.array_type_decl.dimensions;
    d->dimen_count++;
    d = realloc(d, sizeof(*d) + (sizeof(d->dimens[0]) * d->dimen_count));
    d->dimens[d->dimen_count].control_field = NULL;
    super_type->node.array_type_decl.dimensions = d;

    if (array->node.array_type_decl.element_ref->node_type 
	== cod_declaration) {
	/* we're the last in line */
	sm_ref typ = NULL;
	int cg_type = DILL_ERR;
	int ret;
	sm_ref decl = array->node.array_type_decl.element_ref;
	decl->node.declaration.sm_complex_type = super_type;
	decl->node.declaration.cg_type = DILL_B;
	ret = semanticize_decl(context, decl, scope);
	if (ret == 0) return 0;

	if (decl->node.declaration.type_spec != NULL) {
	    typ = reduce_type_list(context, decl->node.declaration.type_spec,
				   &cg_type, scope, NULL, &decl->node.declaration.freeable_complex_type);
	} else {
	    sm_ref arr = decl->node.declaration.sm_complex_type;
	    if ((arr != NULL) && 
		(arr->node_type == cod_array_type_decl)) {
		typ = reduce_type_list(context, 
					arr->node.array_type_decl.type_spec, 
					&cg_type, scope, NULL, &decl->node.declaration.freeable_complex_type);
	    }
	}
	if ((typ == NULL) && (cg_type == DILL_ERR)) return 0;
	array->node.array_type_decl.cg_element_type = cg_type;
	array->node.array_type_decl.sm_complex_element_type = typ;
	super_type->node.array_type_decl.cg_element_type = cg_type;
    } else {
	assert(array->node.array_type_decl.element_ref->node_type == cod_array_type_decl);
	array->node.array_type_decl.sm_complex_element_type = array->node.array_type_decl.element_ref;
	return semanticize_array_element_node(context, 
					      array->node.array_type_decl.element_ref,
					      array,  
					      base_type_spec, scope);
    }
    return 1;
}	

static int
semanticize_array_type_node(cod_parse_context context, sm_ref array, scope_ptr scope)
{
    if (!array->node.array_type_decl.dimensions) {
        array->node.array_type_decl.dimensions = malloc(sizeof(dimen_s));
	memset(array->node.array_type_decl.dimensions, 0, sizeof(dimen_s));
    }
    array->node.array_type_decl.dimensions->dimen_count = 0;
    return semanticize_array_element_node(context, array, array,  
					  array->node.array_type_decl.type_spec,
					  scope);
}

#define Max(i,j) ((i<j) ? j : i)

extern FMTypeDesc*
gen_FMTypeDesc(FMFieldList fl, int field, const char *typ);

static void
free_FMTypeDesc(FMTypeDesc *desc)
{
    while (desc) {
	FMTypeDesc *tmp = desc->next;
	free(desc);
	desc = tmp;
    }
    return;
}

static int
semanticize_struct_type_node(cod_parse_context context, sm_ref decl, 
		      scope_ptr scope)
{
    FMFieldList fl = malloc(sizeof(fl[0]));
    int field_num = 0;
    int ret = 1;
    int struct_size = 0;
    sm_list fields = decl->node.struct_type_decl.fields;
    add_decl_ns(decl->node.struct_type_decl.id, decl, scope, NS_STRUCT);
    while(fields != NULL) {
	field *f = &fields->node->node.field;
	fl[field_num].field_name = f->name;
	fl[field_num].field_type = f->string_type;
	fl = realloc(fl, sizeof(fl[0]) * (field_num + 2));
	field_num++;
	fields = fields->next;
    }
    fl[field_num].field_name = NULL;
    fl[field_num].field_type = NULL;
    field_num = 0;
    fields = decl->node.struct_type_decl.fields;
    while(fields != NULL) {
	field *f = &fields->node->node.field;
	int err = 0;
	int field_size = f->cg_size;

	if (f->string_type != NULL) {
	    /* FFS-compatible field type */
	    FMTypeDesc* desc = gen_FMTypeDesc(fl, field_num, f->string_type);
	    if (desc == NULL) {
		cod_src_error(context, decl, 
			      "Field \"%s\" has unknown type \"%s\".",
			      f->name, f->string_type);
		ret = 0;
	    }
	    build_type_nodes(context, decl, f, fields, f->cg_size, f->cg_type,
			     desc, &err, scope);
	    
	    free_FMTypeDesc(desc);
	    f->cg_type = str_to_data_type(f->string_type, f->cg_size);
	    field_size = f->cg_size;
	    if (f->sm_complex_type) {
		if (f->sm_complex_type->node_type == cod_reference_type_decl) {
		    field_size = sizeof(char*);
		} else if (f->sm_complex_type->node_type == cod_array_type_decl) {
		    sm_ref arr = f->sm_complex_type;
		    while (arr && (arr->node_type == cod_array_type_decl)) {
			if (arr->node.array_type_decl.cg_static_size != -1) {
			    field_size *= arr->node.array_type_decl.cg_static_size;
			}
			arr = arr->node.array_type_decl.sm_complex_element_type;
		    }
		}
	    }
	} else {
	    /* not able to get a FFS-compatible form */
	    int type_def = 0;
	    int cg_type;
	    sm_ref typ;
	    typ = reduce_type_list(context, f->type_spec, &cg_type, scope, 
				   &type_def, &f->freeable_complex_type);
	    f->sm_complex_type = typ;
	    f->cg_type = cg_type;
	    field_size = -1;
	}
	if (err == 1) ret = 0;
	struct_size = Max(struct_size,
			  (f->cg_offset + field_size));
	fields = fields->next;
	field_num++;
    }
    free(fl);
    decl->node.struct_type_decl.cg_size = struct_size;
    return ret;
}

static int
semanticize_enum_type_node(cod_parse_context context, sm_ref decl, 
		      scope_ptr scope)
{
    sm_list enums = decl->node.enum_type_decl.enums;
    while(enums != NULL) {
	if (enums->node->node.enumerator.const_expression) {
	    if (!is_constant_expr(enums->node->node.enumerator.const_expression)) {
		cod_src_error(context, enums->node, 
			      "Enumerator value expression must be constant.");
		return 0;
	    }
	}
	add_decl(enums->node->node.enumerator.id, enums->node, context->scope);
	enums = enums->next;
    }
    add_decl_ns(decl->node.enum_type_decl.id, decl, scope, NS_ENUM);
    return 1;
}

static int
semanticize_reference_type_node(cod_parse_context context, sm_ref decl, 
				scope_ptr scope)
{
    int ret = 1;
    if (decl->node.reference_type_decl.name != NULL) {
	add_decl(decl->node.reference_type_decl.name, decl, scope);
    }
    return ret;
}

extern sm_ref
cod_build_param_node(const char *id, sm_ref typ, int param_num)
{
    sm_ref node = cod_new_declaration();
    sm_ref ident;
    node->node.declaration.param_num = param_num;
    node->node.declaration.id = strdup(id);
    node->node.declaration.sm_complex_type = typ;
    if (typ != NULL) {
	ident = cod_new_identifier();
	node->node.declaration.type_spec = malloc(sizeof(struct list_struct));
	node->node.declaration.type_spec->next = NULL;
	node->node.declaration.type_spec->node = ident;
	ident->node.identifier.id = strdup(typ->node.struct_type_decl.id);
    }
    return node;
}

extern
void
get_FMformat_characteristics(FMFormat format, FMfloat_format *ff, FMinteger_format *intf, int *column_major, int *pointer_size);

static
sm_ref cod_build_type_node_FMformat(FMFormat format, cod_parse_context context)
{
    sm_ref decl = cod_new_struct_type_decl();
    sm_list *end_ptr = &decl->node.struct_type_decl.fields;
    FMfloat_format data_float;
    FMinteger_format data_int;
    int column_major;
    int pointer_size;
    FMFieldList field_list = format->field_list;
    get_FMformat_characteristics(format, &data_float, &data_int, &column_major, &pointer_size);

    decl->node.struct_type_decl.id = strdup(name_of_FMformat(format));
    decl->node.struct_type_decl.encode_info = malloc(sizeof(struct enc_struct));
    decl->node.struct_type_decl.encode_info->byte_order = data_int;
    decl->node.struct_type_decl.encode_info->float_order = data_float;
    decl->node.struct_type_decl.encode_info->pointer_size = pointer_size;
    while ((field_list != NULL) && (field_list->field_name != NULL)) {
	sm_list new_elem;
	char *colon = strchr(field_list->field_type, ':');
	char *bracket = strchr(field_list->field_type, '[');

	if (colon != NULL) {
	    *colon = 0;
	    if (bracket != NULL) strcpy(colon, bracket);
	}
	new_elem = malloc(sizeof(*new_elem));
	new_elem->next = NULL;
	new_elem->node = cod_new_field();
	new_elem->node->node.field.name = strdup(field_list->field_name);
	new_elem->node->node.field.string_type = strdup(field_list->field_type);
	new_elem->node->node.field.cg_size = field_list->field_size;
	new_elem->node->node.field.cg_offset = field_list->field_offset;
	new_elem->node->node.field.cg_type = DILL_ERR;
	*end_ptr = new_elem;
	end_ptr = &new_elem->next;
	field_list++;
    }
    return decl;
}

extern void
cod_add_encoded_param(const char *id, char *data, int param_num, 
		      FMContext c, cod_parse_context context)
{
    int i = 0;
    FMFormat format = FMformat_from_ID(c, data);
    FMFormat *formats;
    sm_ref top_type = NULL, param_node;
    sm_ref node;
    if (format == NULL) {
	printf("No FMFormat ID found in buffer supplied to cod_add_encoded_param()\n");
	printf("No parameter added\n");
	return;
    }
    formats = format->subformats;
    while (formats[i] != NULL) {
	node = cod_build_type_node_FMformat(formats[i], context);
	cod_add_decl_to_parse_context(name_of_FMformat(formats[i]), node, context); 
	cod_add_decl_to_scope(name_of_FMformat(formats[i]), node, context); 
	top_type = node;
	i++;
    }
    
    node = cod_build_type_node_FMformat(format, context);
    cod_add_decl_to_parse_context(name_of_FMformat(format), node, context); 
    cod_add_decl_to_scope(name_of_FMformat(format), node, context); 
    top_type = node;
    param_node = cod_build_param_node(id, NULL, param_num);
    param_node->node.declaration.sm_complex_type = top_type;
    cod_add_decl_to_parse_context(id, param_node, context);
}

extern void
cod_add_param(const char *id, const char *typ, int param_num, 
	      cod_parse_context context)
{
    sm_list type_list;
    sm_ref node;
    setup_for_string_parse(typ, context->defined_types, context->enumerated_constants);
    cod_code_string = typ;
    parsing_type = 1;
    yyerror_count = 0;
    yycontext = context;
    yyparse();
    parsing_type = 0;
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	return;
    }
    type_list = (sm_list) yyparse_value;

    node = cod_build_param_node(id, NULL, param_num);
    node->node.declaration.type_spec = type_list;
    cod_add_decl_to_parse_context(id, node, context);
}

extern void
cod_subroutine_declaration(const char *decl, cod_parse_context context)
{
    sm_list type_list, params;
    sm_ref complex_type, declaration, freeable_complex_type = NULL;
    int cg_type, param_num;

    setup_for_string_parse(decl, context->defined_types, context->enumerated_constants);
    cod_code_string = decl;
    parsing_param_spec = 1;
    yyerror_count = 0;
    yycontext = context;
    yyparse();
    parsing_param_spec = 0;
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	return;
    }
    declaration = yyparse_value;
    context->freeable_declaration = declaration;
    type_list = declaration->node.declaration.type_spec;

    /* handle return type */
    complex_type = reduce_type_list(context, type_list, &cg_type, context->scope, NULL, &freeable_complex_type);
    if (freeable_complex_type) cod_rfree(freeable_complex_type);
 /* context->return_type_list = type_list; - Free'd as part of parse, not here*/
    if (complex_type != NULL) {
	cg_type = DILL_P;
    }
    context->return_cg_type = cg_type;

    /* handle params */
    params = declaration->node.declaration.params;
    param_num = 0;
    while (params != NULL) {
	sm_ref param = NULL;
	switch (params->node->node_type) {
	case cod_declaration:
	    param = params->node;
	    break;
	case cod_array_type_decl:
	    param = params->node->node.array_type_decl.element_ref;
	    param->node.declaration.sm_complex_type = params->node;
	    break;
	default:
	    printf("unhandled case in cod_subroutine_declaration\n");
	}
	param->node.declaration.param_num = param_num;
	cod_add_decl_to_parse_context(param->node.declaration.id,
				      cod_copy(params->node), context);
	param_num++;
	params = params->next;
    }
}

extern void
cod_set_return_type(char *typ, cod_parse_context context)
{
    sm_list type_list;
    int cg_type;
    sm_ref complex_type, freeable_complex_type = NULL;
    setup_for_string_parse(typ, context->defined_types, context->enumerated_constants);
    cod_code_string = typ;
    parsing_type = 1;
    yyerror_count = 0;
    yycontext = context;
    yyparse();
    parsing_type = 0;
    terminate_string_parse();

    if ((yyparse_value == NULL) || (yyerror_count != 0)) {
	return;
    }
    type_list = (sm_list) yyparse_value;

    complex_type = reduce_type_list(context, type_list, &cg_type, context->scope, NULL, &freeable_complex_type);
    context->return_type_list = type_list;
    if (complex_type != NULL) {
	cg_type = DILL_P;
	if (freeable_complex_type) cod_rfree(freeable_complex_type);
    }
    context->return_cg_type = cg_type;
}

static sm_ref
find_complex_type(sm_ref node, scope_ptr scope)
{
    assert(node->node_type == cod_identifier);
    return resolve(node->node.identifier.id, scope);
}

extern cod_parse_context
new_cod_parse_context()
{
    cod_parse_context context = malloc(sizeof(struct parse_struct));
    context->decls = NULL;
    context->standard_decls = NULL;
    context->scope = push_scope(NULL);
    context->defined_types = NULL;
    context->enumerated_constants = NULL;
    context->error_func = default_error_out;
    context->client_data = NULL;
    context->return_type_list = NULL;
    context->return_cg_type = DILL_I;
    context->freeable_declaration = NULL;
    context->has_exec_context = 0;
    context->dont_coerce_return = 0;
    context->alloc_globals = 0;
    cod_add_standard_elements(context);
    return context;
}

extern void
cod_free_parse_context(cod_parse_context parse_context)
{
    if (parse_context->scope->externs) {
        int i = 0;
	while(parse_context->scope->externs[i].extern_name) {
	  free(parse_context->scope->externs[i].extern_name);
	  i++;
	}
	free(parse_context->scope->externs);
    }
    pop_scope(parse_context->scope);
    if (parse_context->defined_types) {
	free(parse_context->defined_types);
    }
    if (parse_context->decls) {
	cod_rfree_list(parse_context->decls, NULL);
    }
    if (parse_context->return_type_list) {
	cod_rfree_list(parse_context->return_type_list, NULL);
    }
    if (parse_context->standard_decls) {
	cod_rfree_list(parse_context->standard_decls, NULL);
    }
    if (parse_context->freeable_declaration) {
        cod_rfree(parse_context->freeable_declaration);
    }
    free(parse_context);
}

extern void
cod_assoc_externs(cod_parse_context context, cod_extern_list externs)
{
    int new_count = 0;

    while(externs[new_count].extern_value) new_count++;

    if (context->scope->externs == NULL) {
        int i;
	context->scope->externs = malloc((new_count+1) * sizeof(externs[0]));
	for (i=0; i < new_count; i++) {
	  context->scope->externs[i].extern_name = strdup(externs[i].extern_name);
	  context->scope->externs[i].extern_value = externs[i].extern_value;
	}
	context->scope->externs[new_count].extern_name = NULL;
	context->scope->externs[new_count].extern_value = NULL;
    } else {
	int old_count = 0;
        int i;
	while(context->scope->externs[old_count++].extern_value);
	context->scope->externs = realloc(context->scope->externs, (new_count + old_count) * sizeof(externs[0]));
	
	for (i=0; i < new_count; i++) {
	    int j;
	    for (j=0; j < old_count-1; j++) {
		if (strcmp(externs[i].extern_name, context->scope->externs[j].extern_name) == 0) {
		    context->scope->externs[j].extern_value = externs[i].extern_value;
		}
	    }
	    context->scope->externs[old_count + i - 1].extern_name = strdup(externs[i].extern_name);
	    context->scope->externs[old_count + i - 1].extern_value = externs[i].extern_value;
	}
	context->scope->externs[new_count + old_count -1].extern_name = NULL;
	context->scope->externs[new_count + old_count -1].extern_value = NULL;
    }
}

extern void
cod_add_decl_to_parse_context(const char *name, sm_ref item, cod_parse_context context)
{
    sm_list *last_ptr = &context->decls;
    sm_list list = context->decls;
    while (list != NULL) {
	last_ptr = &list->next;
	list = list->next;
    }
    *last_ptr = malloc(sizeof(*list));
    (*last_ptr)->next = NULL;
    (*last_ptr)->node = item;
    if (item->node_type == cod_struct_type_decl) {
	cod_add_defined_type((char *)name, context);
    }
}

extern void
cod_add_int_constant_to_parse_context(const char *const_name, int value, cod_parse_context context)
{
    sm_ref constant;
    char str_value[64];
    char *name = strdup(const_name);
    sprintf(str_value, "%d", value);
    constant = cod_new_constant();
    constant->node.constant.token = integer_constant;
    constant->node.constant.const_val = strdup(str_value);
    constant->node.constant.freeable_name = name;
    cod_add_decl_to_scope((char*) name, constant, context);
    cod_add_decl_to_parse_context(name, constant, context);
}

extern void
cod_set_closure(char *name, void* closure_context, cod_parse_context context)
{
    sm_ref decl;
    decl = resolve(name, context->scope);
    assert(decl->node_type == cod_declaration);
    assert(decl->node.declaration.is_subroutine);
    decl->node.declaration.closure_id = closure_context;
}

static void
space_to_underscore(char *str){
    while(*str != '\0'){
	if(isspace(*str))
	    *str = '_';
	    str++;
    }
}

static void
purify_name(FMStructDescList list){
    int i,j;
    for(i=0; list[i].format_name; i++){
	FMFieldList fl = list[i].field_list;
	space_to_underscore((char*)list[i].format_name);
	for(j=0; fl[j].field_name; j++){
	    space_to_underscore((char*)fl[j].field_name);
	    space_to_underscore((char*)fl[j].field_type);
	}
    }
}

static void
uniqueify_names(FMStructDescList list, char *prefix)
{
    int i = 0;
    int prefix_len = (int)strlen(prefix);
    while (list[i].format_name != NULL) {
	int j = 0;
	FMFieldList fl = list[i].field_list;
	char *new_name =
	    malloc(strlen(list[i].format_name) + prefix_len + 1);
	strcpy(new_name, prefix);
	strcpy(new_name + prefix_len, list[i].format_name);
	free((char*)list[i].format_name);
	list[i].format_name = new_name;
	while (fl[j].field_name != 0) {
	    int field_type_len = (int)strlen(fl[j].field_type);
	    char *bracket = strchr(fl[j].field_type, '[');
	    int k;
	    if (bracket != NULL) {
		field_type_len = (int)((intptr_t) bracket - (intptr_t) fl[j].field_type);
	    }
	    for (k = 0; k < i; k++) {
		char *new_type;
		if (strncmp
		    (fl[j].field_type, list[k].format_name + prefix_len,
		     field_type_len) != 0) {
		    /* don't match in N chars */
		    continue;
		}
		if ((list[k].format_name + prefix_len)[field_type_len] !=
		    0) {
		    /* list name is longer */
		    continue;
		}
		new_type =
		    malloc(strlen(fl[j].field_type) + prefix_len + 1);
		strcpy(new_type, prefix);
		strcpy(new_type + prefix_len, fl[j].field_type);
		free((void *) fl[j].field_type);
		fl[j].field_type = new_type;
		break;
	    }
	    j++;
	}
	i++;
    }
    purify_name(list);
}

/* Returns the ecode function which will do message format conversion */
extern cod_code
gen_rollback_code(FMStructDescList format1, FMStructDescList format2, char *xform_code)
{
    /* setup ECL generation */
    /* 
     *  NOTE:  We have to make the type names (structure names)
     *  in format1 and format2 to be unique. 
     *  Because of the nature of the problem, they are likely to be 
     *  identical and this may cause namespace collisions in ECL.
     */
    cod_code code;
    cod_parse_context parse_context = new_cod_parse_context();

    int i = 0;
    uniqueify_names(format1, "f0_");
    while (format1[i].format_name != NULL) {
	cod_add_simple_struct_type(format1[i].format_name,
			    format1[i].field_list, parse_context);
	i++;
    }
    cod_add_param("new", format1[i - 1].format_name, 0, parse_context);

    i = 0;
    uniqueify_names(format2, "f1_");
    while (format2[i].format_name != NULL) {
	cod_add_simple_struct_type(format2[i].format_name,
			    format2[i].field_list, parse_context);
	i++;
    }
    cod_add_param("old", format2[i - 1].format_name, 1, parse_context);

    code = cod_code_gen(xform_code, parse_context);
    cod_free_parse_context(parse_context);

    /* 
     * the "code" structure is the only output from this block,
     * all else is free'd.
     */
    return code;
}

static double
get_constant_float_value(cod_parse_context context, sm_ref expr)
{
    double result;
    switch(expr->node.constant.token) {
    case integer_constant:
    case floating_constant:
	sscanf(expr->node.constant.const_val, "%lg", &result);
	return result;
    case string_constant:
	return 0.0;
    case character_constant:
	return (double)(unsigned char)expr->node.constant.const_val[0];
    default:
	assert(FALSE);
    }
	return 0.0;
}

static intptr_t
get_constant_long_value(cod_parse_context context, sm_ref expr)
{
    double dresult;
    long result;
    switch(expr->node.constant.token) {
    case integer_constant:
	sscanf(expr->node.constant.const_val, "%ld", &result);
	return result;
    case floating_constant:
	sscanf(expr->node.constant.const_val, "%lg", &dresult);
	return (long)dresult;
    case string_constant:
	return -1;
    case character_constant:
	return (long)(unsigned char)expr->node.constant.const_val[0];
    default:
	assert(FALSE);
    }
	return -1;
}

extern sm_ref
evaluate_constant_return_expr(cod_parse_context context, sm_ref expr, int *free_result)
{
    switch(expr->node_type) {
    case cod_constant:
	*free_result = 0;
	return expr;
    case cod_identifier:
	return evaluate_constant_return_expr(context, expr->node.identifier.sm_declaration, free_result);
    case cod_declaration:
	if (!expr->node.declaration.const_var) return NULL;
	return evaluate_constant_return_expr(context, expr->node.identifier.sm_declaration, free_result);
    case cod_operator: {
	sm_ref left = NULL, right = NULL, ret;
	int free_left = 0, free_right = 0;
	int left_token, right_token;
	if (expr->node.operator.left != NULL) {
	    if (!(left = evaluate_constant_return_expr(context, expr->node.operator.left, &free_left))) return NULL;
	    left_token = left->node.constant.token;
	}
	if (expr->node.operator.op == op_sizeof) {
	    int cg_type;
	    sm_ref cast = expr->node.operator.right;
	    sm_ref typ;
	    long size;
	    assert(cast->node_type == cod_cast);
	    typ = reduce_type_list(context, cast->node.cast.type_spec, &cg_type, context? context->scope:NULL, NULL, NULL);
	    static dill_stream s = NULL;
	    char str_val[40];
	    extern int cg_get_size(dill_stream s, sm_ref node);
    
	    if (s == NULL) {
		s = dill_create_stream();
	    }
	    if (typ == NULL) {
		size = dill_type_size(s, cg_type);
	    } else {
		size = cg_get_size(s, cast);
	    }
	    ret = cod_new_constant();
	    ret->node.constant.token = integer_constant;
	    sprintf(str_val, "%ld", size);
	    ret->node.constant.const_val = strdup(str_val);
	    *free_result = 1;
	    return ret;
	}	    
	if (expr->node.operator.right != NULL) {
	    if (!(right = evaluate_constant_return_expr(context, expr->node.operator.right, &free_right))) return NULL;
	    right_token = right->node.constant.token;
	    if (!expr->node.operator.left) {
		left_token = right_token;
	    }
	}
	if ((left_token == string_constant) || 
	    (right_token == string_constant)) {
	    /* blech.  No operators on strings. */
	    return NULL;
	}
	if ((left_token == floating_constant) || 
	    (right_token == floating_constant)) {
	    double left_val, right_val, fvalue;
	    int ivalue, is_ivalue = 0;
	    char str_val[40];
	    if (left)
		left_val = get_constant_float_value(context, left);
	    if (right)
		right_val = get_constant_float_value(context, right);
	    switch(expr->node.operator.op) {
	    case  op_plus:
		fvalue = left_val + right_val;
		break;
	    case  op_minus:
		fvalue = left_val - right_val;
	    break;
	    case  op_leq:
		ivalue = left_val <= right_val;
		is_ivalue=1;
		break;
	    case  op_lt:
		ivalue = left_val < right_val;
		is_ivalue=1;
		break;
	    case  op_geq:
		ivalue = left_val >= right_val;
		is_ivalue=1;
		break;
	    case  op_gt:
		ivalue = left_val > right_val;
		is_ivalue=1;
		break;
	    case  op_eq:
		ivalue = left_val == right_val;
		is_ivalue=1;
		break;
	    case  op_neq:
		ivalue = left_val != right_val;
		is_ivalue=1;
		break;
	    case  op_mult:
		fvalue = left_val * right_val;
		break;
	    case  op_div:
		if (right_val == 0) {
		    return NULL;
		}
		fvalue = left_val / right_val;
		break;
	    case op_log_neg:
	    case op_not:
	    case op_left_shift:
	    case op_right_shift:
	    case  op_modulus:
	    case  op_log_or:
	    case  op_arith_or:
	    case  op_arith_xor:
	    case  op_log_and:
	    case  op_arith_and:
	    case  op_deref:
	    case  op_address:
	    case op_inc:
	    case op_dec:
	    case op_sizeof:
		assert(FALSE);
	    }
	    ret = cod_new_constant();
	    if (is_ivalue) {
		ret->node.constant.token = integer_constant;
		sprintf(str_val, "%d", ivalue);
	    } else {
		ret->node.constant.token = floating_constant;
		sprintf(str_val, "%.*e\n", OP_DBL_Digs - 1, fvalue);
	    }
	    ret->node.constant.const_val = strdup(str_val);
	    *free_result = 1;
	} else {
	    /* we get an integer result */
	    intptr_t left_val = 0, right_val = 0, value;
	    char str_val[40];
	    if (expr->node.operator.left)
		left_val = get_constant_long_value(context, left);
	    if (expr->node.operator.right)
		right_val = get_constant_long_value(context, right);
	    switch(expr->node.operator.op) {
	    case  op_modulus:
		if (right_val == 0) {
		    return NULL;
		}
		value = left_val % right_val;
		break;
	    case  op_plus:
		value = left_val + right_val;
		break;
	    case  op_minus:
		value = left_val - right_val;
		break;
	    case  op_leq:
		value = left_val <= right_val;
		break;
	    case  op_lt:
		value = left_val < right_val;
		break;
	    case  op_geq:
		value = left_val >= right_val;
		break;
	    case  op_gt:
		value = left_val > right_val;
		break;
	    case  op_eq:
		value = left_val = right_val;
		break;
	    case  op_neq:
		value = left_val != right_val;
		break;
	    case  op_log_or:
		value = left_val || right_val;
		break;
	    case  op_arith_or:
		value = left_val | right_val;
		break;
	    case  op_arith_xor:
		value = left_val ^ right_val;
		break;
	    case  op_log_and:
		value = left_val && right_val;
		break;
	    case  op_arith_and:
		value = left_val & right_val;
		break;
	    case  op_mult:
		value = left_val * right_val;
		break;
	    case  op_div:
		value = left_val / right_val;
		break;
	    case op_log_neg:
		value = !right_val;
		break;
	    case op_not:
		value = ~right_val;
		break;
	    case op_left_shift:
		value = left_val << right_val;
		break;
	    case op_right_shift:
		value = left_val >> right_val;
		break;
	    case op_deref:
	    case op_address:
	    case op_inc:
	    case op_dec:
	    case op_sizeof:
		assert(FALSE);
	    }
	    ret = cod_new_constant();
	    ret->node.constant.token = integer_constant;
	    sprintf(str_val, "%zd", value);
	    ret->node.constant.const_val = strdup(str_val);
	    *free_result = 1;
	}
	if (free_left) {
	    /* do stuff*/ 
	} 
	if (free_right) {
	    /* do stuff*/ 
	} 
	return ret;
    }
    case cod_cast:
	return evaluate_constant_return_expr(context, expr->node.cast.expression, free_result);
    case cod_assignment_expression:
    case cod_field_ref:
    case cod_element_ref:
    case cod_subroutine_call:
	assert(FALSE);
    default:
	assert(FALSE);
    }
	return NULL;
}
	
