/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "cod.y"

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

#line 270 "cod.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "cod.tab.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_ARROW = 3,                      /* ARROW  */
  YYSYMBOL_LPAREN = 4,                     /* LPAREN  */
  YYSYMBOL_RPAREN = 5,                     /* RPAREN  */
  YYSYMBOL_LCURLY = 6,                     /* LCURLY  */
  YYSYMBOL_RCURLY = 7,                     /* RCURLY  */
  YYSYMBOL_COLON = 8,                      /* COLON  */
  YYSYMBOL_QUESTION = 9,                   /* QUESTION  */
  YYSYMBOL_LBRACKET = 10,                  /* LBRACKET  */
  YYSYMBOL_RBRACKET = 11,                  /* RBRACKET  */
  YYSYMBOL_DOT = 12,                       /* DOT  */
  YYSYMBOL_STAR = 13,                      /* STAR  */
  YYSYMBOL_AT = 14,                        /* AT  */
  YYSYMBOL_SLASH = 15,                     /* SLASH  */
  YYSYMBOL_MODULUS = 16,                   /* MODULUS  */
  YYSYMBOL_PLUS = 17,                      /* PLUS  */
  YYSYMBOL_MINUS = 18,                     /* MINUS  */
  YYSYMBOL_TILDE = 19,                     /* TILDE  */
  YYSYMBOL_LEQ = 20,                       /* LEQ  */
  YYSYMBOL_LT = 21,                        /* LT  */
  YYSYMBOL_GEQ = 22,                       /* GEQ  */
  YYSYMBOL_GT = 23,                        /* GT  */
  YYSYMBOL_EQ = 24,                        /* EQ  */
  YYSYMBOL_NEQ = 25,                       /* NEQ  */
  YYSYMBOL_LEFT_SHIFT = 26,                /* LEFT_SHIFT  */
  YYSYMBOL_RIGHT_SHIFT = 27,               /* RIGHT_SHIFT  */
  YYSYMBOL_ASSIGN = 28,                    /* ASSIGN  */
  YYSYMBOL_MUL_ASSIGN = 29,                /* MUL_ASSIGN  */
  YYSYMBOL_DIV_ASSIGN = 30,                /* DIV_ASSIGN  */
  YYSYMBOL_MOD_ASSIGN = 31,                /* MOD_ASSIGN  */
  YYSYMBOL_ADD_ASSIGN = 32,                /* ADD_ASSIGN  */
  YYSYMBOL_SUB_ASSIGN = 33,                /* SUB_ASSIGN  */
  YYSYMBOL_LEFT_ASSIGN = 34,               /* LEFT_ASSIGN  */
  YYSYMBOL_RIGHT_ASSIGN = 35,              /* RIGHT_ASSIGN  */
  YYSYMBOL_AND_ASSIGN = 36,                /* AND_ASSIGN  */
  YYSYMBOL_XOR_ASSIGN = 37,                /* XOR_ASSIGN  */
  YYSYMBOL_OR_ASSIGN = 38,                 /* OR_ASSIGN  */
  YYSYMBOL_LOG_OR = 39,                    /* LOG_OR  */
  YYSYMBOL_LOG_AND = 40,                   /* LOG_AND  */
  YYSYMBOL_ARITH_OR = 41,                  /* ARITH_OR  */
  YYSYMBOL_ARITH_AND = 42,                 /* ARITH_AND  */
  YYSYMBOL_ARITH_XOR = 43,                 /* ARITH_XOR  */
  YYSYMBOL_INC_OP = 44,                    /* INC_OP  */
  YYSYMBOL_DEC_OP = 45,                    /* DEC_OP  */
  YYSYMBOL_BANG = 46,                      /* BANG  */
  YYSYMBOL_SEMI = 47,                      /* SEMI  */
  YYSYMBOL_IF = 48,                        /* IF  */
  YYSYMBOL_ELSE = 49,                      /* ELSE  */
  YYSYMBOL_FOR = 50,                       /* FOR  */
  YYSYMBOL_DO = 51,                        /* DO  */
  YYSYMBOL_WHILE = 52,                     /* WHILE  */
  YYSYMBOL_CHAR = 53,                      /* CHAR  */
  YYSYMBOL_SHORT = 54,                     /* SHORT  */
  YYSYMBOL_INT = 55,                       /* INT  */
  YYSYMBOL_LONG = 56,                      /* LONG  */
  YYSYMBOL_UNSIGNED = 57,                  /* UNSIGNED  */
  YYSYMBOL_SIGNED = 58,                    /* SIGNED  */
  YYSYMBOL_FLOAT = 59,                     /* FLOAT  */
  YYSYMBOL_DOUBLE = 60,                    /* DOUBLE  */
  YYSYMBOL_VOID = 61,                      /* VOID  */
  YYSYMBOL_STRING = 62,                    /* STRING  */
  YYSYMBOL_STATIC = 63,                    /* STATIC  */
  YYSYMBOL_EXTERN_TOKEN = 64,              /* EXTERN_TOKEN  */
  YYSYMBOL_STRUCT = 65,                    /* STRUCT  */
  YYSYMBOL_ENUM = 66,                      /* ENUM  */
  YYSYMBOL_UNION = 67,                     /* UNION  */
  YYSYMBOL_CONST = 68,                     /* CONST  */
  YYSYMBOL_SIZEOF = 69,                    /* SIZEOF  */
  YYSYMBOL_TYPEDEF = 70,                   /* TYPEDEF  */
  YYSYMBOL_RETURN_TOKEN = 71,              /* RETURN_TOKEN  */
  YYSYMBOL_CONTINUE = 72,                  /* CONTINUE  */
  YYSYMBOL_BREAK = 73,                     /* BREAK  */
  YYSYMBOL_GOTO = 74,                      /* GOTO  */
  YYSYMBOL_PRINT = 75,                     /* PRINT  */
  YYSYMBOL_COMMA = 76,                     /* COMMA  */
  YYSYMBOL_DOTDOTDOT = 77,                 /* DOTDOTDOT  */
  YYSYMBOL_integer_constant = 78,          /* integer_constant  */
  YYSYMBOL_character_constant = 79,        /* character_constant  */
  YYSYMBOL_string_constant = 80,           /* string_constant  */
  YYSYMBOL_floating_constant = 81,         /* floating_constant  */
  YYSYMBOL_identifier_ref = 82,            /* identifier_ref  */
  YYSYMBOL_type_identifier = 83,           /* type_identifier  */
  YYSYMBOL_enumeration_constant = 84,      /* enumeration_constant  */
  YYSYMBOL_YYACCEPT = 85,                  /* $accept  */
  YYSYMBOL_start = 86,                     /* start  */
  YYSYMBOL_primary_expression = 87,        /* primary_expression  */
  YYSYMBOL_postfix_expression = 88,        /* postfix_expression  */
  YYSYMBOL_argument_expression_list = 89,  /* argument_expression_list  */
  YYSYMBOL_unary_expression = 90,          /* unary_expression  */
  YYSYMBOL_unary_operator = 91,            /* unary_operator  */
  YYSYMBOL_cast_expression = 92,           /* cast_expression  */
  YYSYMBOL_multiplicative_expression = 93, /* multiplicative_expression  */
  YYSYMBOL_additive_expression = 94,       /* additive_expression  */
  YYSYMBOL_shift_expression = 95,          /* shift_expression  */
  YYSYMBOL_relational_expression = 96,     /* relational_expression  */
  YYSYMBOL_equality_expression = 97,       /* equality_expression  */
  YYSYMBOL_and_expression = 98,            /* and_expression  */
  YYSYMBOL_exclusive_or_expression = 99,   /* exclusive_or_expression  */
  YYSYMBOL_inclusive_or_expression = 100,  /* inclusive_or_expression  */
  YYSYMBOL_logical_and_expression = 101,   /* logical_and_expression  */
  YYSYMBOL_logical_or_expression = 102,    /* logical_or_expression  */
  YYSYMBOL_conditional_expression = 103,   /* conditional_expression  */
  YYSYMBOL_assignment_operator = 104,      /* assignment_operator  */
  YYSYMBOL_assignment_expression = 105,    /* assignment_expression  */
  YYSYMBOL_expression = 106,               /* expression  */
  YYSYMBOL_constant_expression = 107,      /* constant_expression  */
  YYSYMBOL_init_declarator_list = 108,     /* init_declarator_list  */
  YYSYMBOL_declaration = 109,              /* declaration  */
  YYSYMBOL_110_1 = 110,                    /* $@1  */
  YYSYMBOL_111_2 = 111,                    /* @2  */
  YYSYMBOL_declaration_specifiers = 112,   /* declaration_specifiers  */
  YYSYMBOL_init_declarator = 113,          /* init_declarator  */
  YYSYMBOL_storage_class_specifier = 114,  /* storage_class_specifier  */
  YYSYMBOL_type_specifier = 115,           /* type_specifier  */
  YYSYMBOL_struct_or_union_specifier = 116, /* struct_or_union_specifier  */
  YYSYMBOL_struct_or_union = 117,          /* struct_or_union  */
  YYSYMBOL_struct_declaration_list = 118,  /* struct_declaration_list  */
  YYSYMBOL_struct_declaration = 119,       /* struct_declaration  */
  YYSYMBOL_struct_declarator_list = 120,   /* struct_declarator_list  */
  YYSYMBOL_struct_declarator = 121,        /* struct_declarator  */
  YYSYMBOL_specifier_qualifier_list = 122, /* specifier_qualifier_list  */
  YYSYMBOL_enum_specifier = 123,           /* enum_specifier  */
  YYSYMBOL_enumerator_list = 124,          /* enumerator_list  */
  YYSYMBOL_enumerator = 125,               /* enumerator  */
  YYSYMBOL_type_qualifier = 126,           /* type_qualifier  */
  YYSYMBOL_declarator = 127,               /* declarator  */
  YYSYMBOL_direct_declarator = 128,        /* direct_declarator  */
  YYSYMBOL_pointer = 129,                  /* pointer  */
  YYSYMBOL_type_qualifier_list = 130,      /* type_qualifier_list  */
  YYSYMBOL_parameter_type_list = 131,      /* parameter_type_list  */
  YYSYMBOL_parameter_list = 132,           /* parameter_list  */
  YYSYMBOL_parameter_declaration = 133,    /* parameter_declaration  */
  YYSYMBOL_type_name = 134,                /* type_name  */
  YYSYMBOL_abstract_declarator = 135,      /* abstract_declarator  */
  YYSYMBOL_initializer = 136,              /* initializer  */
  YYSYMBOL_initializer_list = 137,         /* initializer_list  */
  YYSYMBOL_designation = 138,              /* designation  */
  YYSYMBOL_designator_list = 139,          /* designator_list  */
  YYSYMBOL_designator = 140,               /* designator  */
  YYSYMBOL_decls_stmts_list = 141,         /* decls_stmts_list  */
  YYSYMBOL_statement = 142,                /* statement  */
  YYSYMBOL_labeled_statement = 143,        /* labeled_statement  */
  YYSYMBOL_compound_statement = 144,       /* compound_statement  */
  YYSYMBOL_declaration_list = 145,         /* declaration_list  */
  YYSYMBOL_jump_statement = 146,           /* jump_statement  */
  YYSYMBOL_expression_statement = 147,     /* expression_statement  */
  YYSYMBOL_selection_statement = 148,      /* selection_statement  */
  YYSYMBOL_iteration_statement = 149,      /* iteration_statement  */
  YYSYMBOL_expression_opt = 150,           /* expression_opt  */
  YYSYMBOL_constant = 151                  /* constant  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  88
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   985

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  85
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  67
/* YYNRULES -- Number of rules.  */
#define YYNRULES  206
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  321

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   339


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   358,   358,   362,   368,   374,   376,   383,   385,   392,
     399,   406,   412,   418,   425,   435,   441,   455,   456,   463,
     470,   477,   484,   499,   502,   505,   508,   511,   514,   520,
     521,   530,   532,   541,   550,   561,   563,   572,   583,   585,
     594,   605,   607,   616,   625,   634,   645,   647,   656,   667,
     669,   680,   682,   693,   695,   706,   708,   719,   721,   732,
     734,   745,   747,   749,   751,   753,   755,   757,   759,   761,
     763,   765,   770,   773,   784,   786,   796,   800,   805,   824,
     831,   823,   914,   920,   925,   931,   936,   942,   947,   955,
     957,   973,   978,   983,   992,   997,  1002,  1007,  1012,  1017,
    1022,  1027,  1032,  1037,  1042,  1047,  1050,  1056,  1059,  1062,
    1068,  1069,  1075,  1076,  1088,  1089,  1140,  1145,  1157,  1160,
    1166,  1171,  1177,  1185,  1192,  1199,  1206,  1213,  1223,  1229,
    1239,  1245,  1253,  1261,  1263,  1277,  1286,  1289,  1298,  1307,
    1314,  1324,  1332,  1340,  1348,  1362,  1373,  1384,  1395,  1415,
    1420,  1433,  1434,  1449,  1455,  1470,  1479,  1518,  1519,  1565,
    1569,  1574,  1579,  1584,  1592,  1600,  1613,  1629,  1634,  1639,
    1652,  1658,  1667,  1673,  1676,  1679,  1685,  1690,  1691,  1692,
    1693,  1694,  1695,  1702,  1709,  1712,  1720,  1722,  1736,  1741,
    1746,  1752,  1758,  1767,  1771,  1781,  1790,  1807,  1817,  1827,
    1841,  1843,  1846,  1853,  1860,  1867,  1874
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "ARROW", "LPAREN",
  "RPAREN", "LCURLY", "RCURLY", "COLON", "QUESTION", "LBRACKET",
  "RBRACKET", "DOT", "STAR", "AT", "SLASH", "MODULUS", "PLUS", "MINUS",
  "TILDE", "LEQ", "LT", "GEQ", "GT", "EQ", "NEQ", "LEFT_SHIFT",
  "RIGHT_SHIFT", "ASSIGN", "MUL_ASSIGN", "DIV_ASSIGN", "MOD_ASSIGN",
  "ADD_ASSIGN", "SUB_ASSIGN", "LEFT_ASSIGN", "RIGHT_ASSIGN", "AND_ASSIGN",
  "XOR_ASSIGN", "OR_ASSIGN", "LOG_OR", "LOG_AND", "ARITH_OR", "ARITH_AND",
  "ARITH_XOR", "INC_OP", "DEC_OP", "BANG", "SEMI", "IF", "ELSE", "FOR",
  "DO", "WHILE", "CHAR", "SHORT", "INT", "LONG", "UNSIGNED", "SIGNED",
  "FLOAT", "DOUBLE", "VOID", "STRING", "STATIC", "EXTERN_TOKEN", "STRUCT",
  "ENUM", "UNION", "CONST", "SIZEOF", "TYPEDEF", "RETURN_TOKEN",
  "CONTINUE", "BREAK", "GOTO", "PRINT", "COMMA", "DOTDOTDOT",
  "integer_constant", "character_constant", "string_constant",
  "floating_constant", "identifier_ref", "type_identifier",
  "enumeration_constant", "$accept", "start", "primary_expression",
  "postfix_expression", "argument_expression_list", "unary_expression",
  "unary_operator", "cast_expression", "multiplicative_expression",
  "additive_expression", "shift_expression", "relational_expression",
  "equality_expression", "and_expression", "exclusive_or_expression",
  "inclusive_or_expression", "logical_and_expression",
  "logical_or_expression", "conditional_expression", "assignment_operator",
  "assignment_expression", "expression", "constant_expression",
  "init_declarator_list", "declaration", "$@1", "@2",
  "declaration_specifiers", "init_declarator", "storage_class_specifier",
  "type_specifier", "struct_or_union_specifier", "struct_or_union",
  "struct_declaration_list", "struct_declaration",
  "struct_declarator_list", "struct_declarator",
  "specifier_qualifier_list", "enum_specifier", "enumerator_list",
  "enumerator", "type_qualifier", "declarator", "direct_declarator",
  "pointer", "type_qualifier_list", "parameter_type_list",
  "parameter_list", "parameter_declaration", "type_name",
  "abstract_declarator", "initializer", "initializer_list", "designation",
  "designator_list", "designator", "decls_stmts_list", "statement",
  "labeled_statement", "compound_statement", "declaration_list",
  "jump_statement", "expression_statement", "selection_statement",
  "iteration_statement", "expression_opt", "constant", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-239)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-1)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     808,   232,  -239,  -239,  -239,  -239,  -239,  -239,  -239,  -239,
    -239,  -239,  -239,  -239,  -239,    46,  -239,  -239,  -239,  -239,
      64,  -239,    13,   352,   352,  -239,    47,  -239,   352,  -239,
     352,    39,   385,  -239,  -239,  -239,  -239,  -239,  -239,   653,
     653,  -239,  -239,    71,    84,   466,   146,   723,   539,   108,
     118,    85,  -239,  -239,  -239,  -239,   162,  -239,  -239,   170,
     194,   744,  -239,    32,    25,    11,   140,   128,   129,   132,
     144,   147,     7,  -239,  -239,   -18,  -239,   313,  -239,  -239,
    -239,  -239,  -239,  -239,  -239,  -239,   107,   184,  -239,  -239,
      17,  -239,  -239,   902,   186,  -239,  -239,  -239,  -239,    21,
     902,    42,   902,   188,   744,  -239,  -239,   744,   744,   142,
     744,   385,  -239,  -239,   -11,  -239,  -239,   148,   466,   114,
     583,   744,   116,  -239,  -239,  -239,  -239,  -239,  -239,  -239,
    -239,  -239,  -239,  -239,  -239,  -239,   744,  -239,  -239,   744,
     744,   744,   744,   744,   744,   744,   744,   744,   744,   744,
     744,   744,   744,   744,   744,   744,   744,   744,  -239,   744,
    -239,  -239,  -239,   172,    59,  -239,   107,    17,   135,   135,
     201,   130,  -239,   183,   203,    20,   717,  -239,    -3,   902,
    -239,  -239,  -239,  -239,  -239,   744,    28,   159,   169,   233,
      49,   235,  -239,  -239,  -239,  -239,  -239,    54,  -239,    65,
    -239,  -239,  -239,  -239,  -239,    32,    32,    25,    25,    11,
      11,    11,    11,   140,   140,   128,   129,   132,   144,    27,
     147,  -239,   744,  -239,     0,    60,   236,  -239,  -239,   135,
    -239,   135,   776,    17,   195,   631,   674,   203,  -239,  -239,
    -239,    -7,  -239,  -239,   840,  -239,   466,   744,   744,   466,
    -239,  -239,   744,  -239,   744,  -239,  -239,  -239,  -239,  -239,
      18,  -239,  -239,  -239,  -239,  -239,    17,   238,   168,  -239,
    -239,  -239,   561,  -239,  -239,  -239,   237,  -239,    17,  -239,
     197,   200,    56,  -239,  -239,  -239,  -239,  -239,  -239,   871,
     744,   171,  -239,    61,   631,    29,  -239,  -239,  -239,   466,
     744,   205,  -239,  -239,   244,  -239,  -239,   482,  -239,  -239,
    -239,  -239,   252,  -239,  -239,  -239,  -239,   631,   466,  -239,
    -239
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,    94,    95,    96,    97,   102,   101,    98,    99,
     100,   103,    92,    93,   110,     0,   111,   132,    91,   104,
       0,   186,    79,    83,    85,   105,     0,   106,    87,     3,
       2,     0,     0,   184,    24,    25,    26,    27,    23,     0,
       0,    28,   193,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   202,   205,   204,   203,     4,   206,     7,    17,
      29,     0,    31,    35,    38,    41,    46,    49,    51,    53,
      55,    57,    59,    72,    74,     0,   173,     0,   172,   177,
     178,   182,   179,   180,   181,     5,     0,   127,     1,    82,
       0,    84,    86,     0,   109,    88,   187,   174,     4,     0,
     120,   157,   122,     0,     0,    18,    19,     0,   200,     0,
       0,     0,    21,   189,     0,   190,   191,     0,     0,     0,
       0,     0,     0,    13,    14,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,     0,    29,    20,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   194,     0,
     185,   176,   175,   131,     0,   128,     0,     0,   141,   145,
     135,    80,    77,    89,   133,     0,     0,   112,     0,     0,
       6,   119,   159,   158,   121,     0,     0,   201,     0,     0,
       0,     0,   188,   192,   183,    12,    11,     0,    15,     0,
       9,    73,    32,    33,    34,    36,    37,    39,    40,    44,
      42,    45,    43,    47,    48,    50,    52,    54,    56,     0,
      58,    75,     0,   123,     0,     0,     0,   149,   143,   142,
     147,   146,     0,     0,     0,     0,     0,   134,   108,   113,
     114,     0,   116,   118,     0,    30,     0,   200,     0,     0,
      22,    10,     0,     8,     0,    76,   130,   124,   129,   125,
       0,   136,   150,   144,   148,   138,   155,     0,   151,   153,
      78,    81,     0,   162,    90,   140,     0,   115,     0,   107,
     195,     0,     0,   198,    16,    60,   126,   156,   137,     0,
       0,     0,   164,     0,     0,     0,   168,   139,   117,     0,
     200,     0,   152,   154,     0,   171,   160,     0,   163,   167,
     169,   196,     0,   199,   170,   161,   166,     0,     0,   165,
     197
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -239,  -239,  -239,  -239,  -239,   -34,  -239,   -46,   -71,    12,
      30,    33,   106,   109,   105,   110,   103,  -239,   -82,  -239,
    -108,   -30,   -93,  -239,    19,  -239,  -239,   -20,    34,  -239,
      31,  -239,  -239,    82,  -106,  -239,   -15,   -13,  -239,    98,
    -210,   -10,  -144,    91,   -84,    99,  -239,  -239,   -17,   160,
    -239,  -126,  -239,   -26,  -239,    14,  -239,   -45,  -239,   307,
    -239,  -239,  -239,  -239,  -239,  -238,  -239
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,    20,    58,    59,   197,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,   136,
      74,    75,   256,   171,    21,    90,   234,    22,   172,    23,
      24,    25,    26,   176,   177,   241,   242,   178,    27,   164,
     165,    28,   173,   174,   175,   229,   267,   268,   269,   103,
     183,   274,   293,   294,   295,   296,    77,    78,    79,    80,
      30,    81,    82,    83,    84,   188,    85
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     109,   167,    99,    91,    92,   105,   106,   257,    95,   281,
     168,   169,   198,   112,   258,   138,   156,   182,   114,   101,
      76,   167,   102,   226,   167,   286,   180,   137,   201,   158,
     168,   169,   162,   246,   243,   254,   192,   144,   145,   290,
     277,   291,   142,   143,   240,   139,   157,   140,   141,    96,
     258,   221,    86,    93,   249,   168,   169,   309,   159,   251,
      89,   301,   312,   100,    88,   159,   223,   259,   306,   278,
     239,   205,   206,   194,    99,   107,   253,   186,   187,   170,
     190,    99,   163,   102,   228,   230,    97,   181,   108,   184,
     102,   199,   102,   202,   203,   204,   161,   159,   101,   170,
     163,   102,   170,   159,   159,   137,   137,   137,   137,   137,
     137,   137,   137,   137,   137,   137,   137,   137,   137,   137,
     137,   137,   287,   137,   100,   159,   219,   273,    87,    94,
     252,   100,   159,   100,   243,   224,   260,   307,   239,   245,
     255,   159,   100,   276,   284,   263,   292,   264,   168,   169,
     110,   137,   150,   151,   255,   115,   207,   208,   227,   227,
     146,   147,   148,   149,   273,   116,   102,   117,   308,   102,
     118,   152,   285,   119,   120,   153,   209,   210,   211,   212,
     121,   316,   122,   213,   214,   154,   273,   155,   137,   163,
     166,   319,   179,   185,   189,   193,   195,   304,   200,   273,
     222,   280,   137,    17,   283,   232,   233,   100,   255,   273,
     100,   235,   266,   236,   123,   124,   247,   187,   282,   262,
     137,   262,   125,   126,   127,   128,   129,   130,   131,   132,
     133,   134,   135,    31,   102,   159,    32,   248,     1,    33,
     250,   261,   271,   288,   289,    34,   299,   300,   297,    35,
      36,    37,   313,   305,   311,   314,   137,   318,   215,   217,
     220,   244,   216,   298,   225,   218,   237,   270,   231,   266,
     187,   191,   303,   320,    38,   100,    39,    40,    41,    42,
      43,   317,    44,    45,    46,     2,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    47,    18,    48,    49,    50,    51,    29,     0,   310,
      52,    53,    54,    55,    56,    19,    57,    32,     0,     1,
     160,     0,     0,     0,     0,     0,    34,     0,     0,     0,
      35,    36,    37,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    38,     0,    39,    40,    41,
      42,    43,     0,    44,    45,    46,     2,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    47,    18,    48,    49,    50,    51,     0,    32,
       0,    52,    53,    54,    55,    56,    19,    57,    34,     0,
       0,     0,    35,    36,    37,     2,     3,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,     0,    18,     0,     0,     0,     0,    38,     0,    39,
      40,    41,     0,     0,     0,    19,     0,     0,     2,     3,
       4,     5,     6,     7,     8,     9,    10,    11,     0,     0,
      14,    15,    16,    17,    47,     0,     0,     0,     0,     0,
       0,     0,     0,    52,    53,    54,    55,    98,    19,    57,
      32,     0,     1,     0,     0,     0,     0,     0,     0,    34,
       0,     0,     0,    35,    36,    37,    32,     0,   272,   315,
       0,     0,   290,     0,   291,    34,     0,     0,     0,    35,
      36,    37,     0,     0,     0,     0,     0,     0,    38,     0,
      39,    40,    41,    42,    43,     0,    44,    45,    46,     0,
       0,     0,     0,     0,    38,     0,    39,    40,    41,     0,
       0,     0,     0,     0,     0,    47,     0,    48,    49,    50,
      51,     0,     0,    32,    52,    53,    54,    55,    56,     0,
      57,    47,    34,     0,     0,     0,    35,    36,    37,     0,
      52,    53,    54,    55,    98,    32,    57,   272,     0,     0,
       0,   290,     0,   291,    34,     0,     0,     0,    35,    36,
      37,    38,     0,    39,    40,    41,   113,    32,   196,     0,
       0,     0,     0,     0,     0,     0,    34,     0,     0,     0,
      35,    36,    37,    38,     0,    39,    40,    41,    47,     0,
       0,     0,     0,     0,     0,     0,     0,    52,    53,    54,
      55,    98,     0,    57,     0,    38,     0,    39,    40,    41,
      47,     0,     0,     0,     0,    32,     0,   272,     0,    52,
      53,    54,    55,    98,    34,    57,     0,     0,    35,    36,
      37,     0,    47,     0,     0,     0,     0,   104,     0,     0,
       0,    52,    53,    54,    55,    98,    34,    57,     0,     0,
      35,    36,    37,    38,     0,    39,    40,    41,    32,     0,
       0,     0,     0,     0,     0,   275,     0,    34,     0,     0,
       0,    35,    36,    37,     0,    38,     0,    39,    40,    41,
      47,     0,     0,     0,     0,     0,     0,     0,     0,    52,
      53,    54,    55,    98,     0,    57,    38,     0,    39,    40,
      41,     0,    47,     0,   238,     0,     0,   111,     0,     0,
       0,    52,    53,    54,    55,    98,    34,    57,     0,     0,
      35,    36,    37,    47,     0,     0,     0,     0,    32,     0,
       0,     0,    52,    53,    54,    55,    98,    34,    57,     0,
       0,    35,    36,    37,     0,    38,     0,    39,    40,    41,
       2,     3,     4,     5,     6,     7,     8,     9,    10,    11,
       0,   265,    14,    15,    16,    17,    38,     0,    39,    40,
      41,     0,    47,     0,     0,     0,     0,     0,     0,     0,
      19,    52,    53,    54,    55,    98,     0,    57,     0,     0,
       0,     0,     0,    47,     1,     0,     0,     0,     0,     0,
       0,     0,    52,    53,    54,    55,    98,     0,    57,     2,
       3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,     0,    18,   279,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    19,
       0,     2,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,     0,    18,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    19,     0,     2,     3,     4,     5,     6,     7,     8,
       9,    10,    11,     0,     0,    14,    15,    16,    17,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    19,     2,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
       0,    18,     0,     0,     0,     0,     0,     0,   302,     0,
       0,     0,     0,     0,    19,     2,     3,     4,     5,     6,
       7,     8,     9,    10,    11,     0,     0,    14,    15,    16,
      17,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    19
};

static const yytype_int16 yycheck[] =
{
      45,     4,    32,    23,    24,    39,    40,     7,    28,   247,
      13,    14,   120,    47,   224,    61,     9,   101,    48,    32,
       1,     4,    32,   167,     4,     7,     5,    61,   136,    47,
      13,    14,    77,     5,   178,     8,    47,    26,    27,    10,
      47,    12,    17,    18,    47,    13,    39,    15,    16,    30,
     260,   159,     6,     6,     5,    13,    14,    28,    76,     5,
      47,     5,   300,    32,     0,    76,     7,     7,     7,    76,
     176,   142,   143,   118,   104,     4,    11,   107,   108,    82,
     110,   111,    82,    93,   168,   169,    47,   100,     4,   102,
     100,   121,   102,   139,   140,   141,    77,    76,   111,    82,
      82,   111,    82,    76,    76,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   152,   153,
     154,   155,   266,   157,    93,    76,   156,   235,    82,    82,
      76,   100,    76,   102,   278,    76,    76,    76,   244,   185,
     222,    76,   111,   236,   252,   229,   272,   231,    13,    14,
       4,   185,    24,    25,   236,    47,   144,   145,   168,   169,
      20,    21,    22,    23,   272,    47,   176,    82,   294,   179,
       8,    42,   254,     3,     4,    43,   146,   147,   148,   149,
      10,   307,    12,   150,   151,    41,   294,    40,   222,    82,
       6,   317,     6,     5,    52,    47,    82,   290,    82,   307,
      28,   246,   236,    68,   249,     4,    76,   176,   290,   317,
     179,    28,   232,    10,    44,    45,    47,   247,   248,   229,
     254,   231,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,     1,   244,    76,     4,     4,     6,     7,
       5,     5,    47,     5,    76,    13,    49,    47,    11,    17,
      18,    19,    47,    82,   299,    11,   290,     5,   152,   154,
     157,   179,   153,   278,   166,   155,   175,   233,   169,   289,
     300,   111,   289,   318,    42,   244,    44,    45,    46,    47,
      48,   307,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,     0,    -1,   295,
      78,    79,    80,    81,    82,    83,    84,     4,    -1,     6,
       7,    -1,    -1,    -1,    -1,    -1,    13,    -1,    -1,    -1,
      17,    18,    19,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    42,    -1,    44,    45,    46,
      47,    48,    -1,    50,    51,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    -1,     4,
      -1,    78,    79,    80,    81,    82,    83,    84,    13,    -1,
      -1,    -1,    17,    18,    19,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    -1,    70,    -1,    -1,    -1,    -1,    42,    -1,    44,
      45,    46,    -1,    -1,    -1,    83,    -1,    -1,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    -1,    -1,
      65,    66,    67,    68,    69,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    78,    79,    80,    81,    82,    83,    84,
       4,    -1,     6,    -1,    -1,    -1,    -1,    -1,    -1,    13,
      -1,    -1,    -1,    17,    18,    19,     4,    -1,     6,     7,
      -1,    -1,    10,    -1,    12,    13,    -1,    -1,    -1,    17,
      18,    19,    -1,    -1,    -1,    -1,    -1,    -1,    42,    -1,
      44,    45,    46,    47,    48,    -1,    50,    51,    52,    -1,
      -1,    -1,    -1,    -1,    42,    -1,    44,    45,    46,    -1,
      -1,    -1,    -1,    -1,    -1,    69,    -1,    71,    72,    73,
      74,    -1,    -1,     4,    78,    79,    80,    81,    82,    -1,
      84,    69,    13,    -1,    -1,    -1,    17,    18,    19,    -1,
      78,    79,    80,    81,    82,     4,    84,     6,    -1,    -1,
      -1,    10,    -1,    12,    13,    -1,    -1,    -1,    17,    18,
      19,    42,    -1,    44,    45,    46,    47,     4,     5,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    13,    -1,    -1,    -1,
      17,    18,    19,    42,    -1,    44,    45,    46,    69,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    78,    79,    80,
      81,    82,    -1,    84,    -1,    42,    -1,    44,    45,    46,
      69,    -1,    -1,    -1,    -1,     4,    -1,     6,    -1,    78,
      79,    80,    81,    82,    13,    84,    -1,    -1,    17,    18,
      19,    -1,    69,    -1,    -1,    -1,    -1,     4,    -1,    -1,
      -1,    78,    79,    80,    81,    82,    13,    84,    -1,    -1,
      17,    18,    19,    42,    -1,    44,    45,    46,     4,    -1,
      -1,    -1,    -1,    -1,    -1,    11,    -1,    13,    -1,    -1,
      -1,    17,    18,    19,    -1,    42,    -1,    44,    45,    46,
      69,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    78,
      79,    80,    81,    82,    -1,    84,    42,    -1,    44,    45,
      46,    -1,    69,    -1,     7,    -1,    -1,     4,    -1,    -1,
      -1,    78,    79,    80,    81,    82,    13,    84,    -1,    -1,
      17,    18,    19,    69,    -1,    -1,    -1,    -1,     4,    -1,
      -1,    -1,    78,    79,    80,    81,    82,    13,    84,    -1,
      -1,    17,    18,    19,    -1,    42,    -1,    44,    45,    46,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      -1,     5,    65,    66,    67,    68,    42,    -1,    44,    45,
      46,    -1,    69,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      83,    78,    79,    80,    81,    82,    -1,    84,    -1,    -1,
      -1,    -1,    -1,    69,     6,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    78,    79,    80,    81,    82,    -1,    84,    53,
      54,    55,    56,    57,    58,    59,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    -1,    70,     7,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    83,
      -1,    53,    54,    55,    56,    57,    58,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    -1,    70,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    83,    -1,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    -1,    -1,    65,    66,    67,    68,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    83,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    63,    64,    65,    66,    67,    68,
      -1,    70,    -1,    -1,    -1,    -1,    -1,    -1,    77,    -1,
      -1,    -1,    -1,    -1,    83,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    -1,    -1,    65,    66,    67,
      68,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    83
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     6,    53,    54,    55,    56,    57,    58,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    70,    83,
      86,   109,   112,   114,   115,   116,   117,   123,   126,   144,
     145,     1,     4,     7,    13,    17,    18,    19,    42,    44,
      45,    46,    47,    48,    50,    51,    52,    69,    71,    72,
      73,    74,    78,    79,    80,    81,    82,    84,    87,    88,
      90,    91,    92,    93,    94,    95,    96,    97,    98,    99,
     100,   101,   102,   103,   105,   106,   109,   141,   142,   143,
     144,   146,   147,   148,   149,   151,     6,    82,     0,    47,
     110,   112,   112,     6,    82,   112,   109,    47,    82,   106,
     115,   122,   126,   134,     4,    90,    90,     4,     4,   142,
       4,     4,    90,    47,   106,    47,    47,    82,     8,     3,
       4,    10,    12,    44,    45,    28,    29,    30,    31,    32,
      33,    34,    35,    36,    37,    38,   104,    90,    92,    13,
      15,    16,    17,    18,    26,    27,    20,    21,    22,    23,
      24,    25,    42,    43,    41,    40,     9,    39,    47,    76,
       7,   109,   142,    82,   124,   125,     6,     4,    13,    14,
      82,   108,   113,   127,   128,   129,   118,   119,   122,     6,
       5,   122,   129,   135,   122,     5,   106,   106,   150,    52,
     106,   134,    47,    47,   142,    82,     5,    89,   105,   106,
      82,   105,    92,    92,    92,    93,    93,    94,    94,    95,
      95,    95,    95,    96,    96,    97,    98,    99,   100,   106,
     101,   105,    28,     7,    76,   124,   127,   126,   129,   130,
     129,   130,     4,    76,   111,    28,    10,   128,     7,   119,
      47,   120,   121,   127,   118,    92,     5,    47,     4,     5,
       5,     5,    76,    11,     8,   103,   107,     7,   125,     7,
      76,     5,   126,   129,   129,     5,   112,   131,   132,   133,
     113,    47,     6,   105,   136,    11,   107,    47,    76,     7,
     142,   150,   106,   142,   105,   103,     7,   127,     5,    76,
      10,    12,   136,   137,   138,   139,   140,    11,   121,    49,
      47,     5,    77,   133,   107,    82,     7,    76,   136,    28,
     140,   142,   150,    47,    11,     7,   136,   138,     5,   136,
     142
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_uint8 yyr1[] =
{
       0,    85,    86,    86,    87,    87,    87,    88,    88,    88,
      88,    88,    88,    88,    88,    89,    89,    90,    90,    90,
      90,    90,    90,    91,    91,    91,    91,    91,    91,    92,
      92,    93,    93,    93,    93,    94,    94,    94,    95,    95,
      95,    96,    96,    96,    96,    96,    97,    97,    97,    98,
      98,    99,    99,   100,   100,   101,   101,   102,   102,   103,
     103,   104,   104,   104,   104,   104,   104,   104,   104,   104,
     104,   104,   105,   105,   106,   106,   107,   108,   108,   110,
     111,   109,   109,   112,   112,   112,   112,   112,   112,   113,
     113,   114,   114,   114,   115,   115,   115,   115,   115,   115,
     115,   115,   115,   115,   115,   115,   115,   116,   116,   116,
     117,   117,   118,   118,   119,   119,   120,   120,   121,   122,
     122,   122,   122,   123,   123,   123,   123,   123,   124,   124,
     125,   125,   126,   127,   127,   128,   128,   128,   128,   128,
     128,   129,   129,   129,   129,   129,   129,   129,   129,   130,
     130,   131,   131,   132,   132,   133,   133,   134,   134,   135,
     136,   136,   136,   137,   137,   137,   137,   138,   139,   139,
     140,   140,   141,   141,   141,   141,   141,   142,   142,   142,
     142,   142,   142,   143,   144,   144,   145,   145,   146,   146,
     146,   146,   146,   147,   147,   148,   148,   149,   149,   149,
     150,   150,   151,   151,   151,   151,   151
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     1,     1,     1,     3,     1,     4,     3,
       4,     3,     3,     2,     2,     1,     3,     1,     2,     2,
       2,     2,     4,     1,     1,     1,     1,     1,     1,     1,
       4,     1,     3,     3,     3,     1,     3,     3,     1,     3,
       3,     1,     3,     3,     3,     3,     1,     3,     3,     1,
       3,     1,     3,     1,     3,     1,     3,     1,     3,     1,
       5,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     3,     1,     3,     1,     1,     3,     0,
       0,     5,     2,     1,     2,     1,     2,     1,     2,     1,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     5,     4,     2,
       1,     1,     1,     2,     2,     3,     1,     3,     1,     2,
       1,     2,     1,     4,     5,     5,     6,     2,     1,     3,
       3,     1,     1,     1,     2,     1,     3,     4,     3,     4,
       3,     1,     2,     2,     3,     1,     2,     2,     3,     1,
       2,     1,     3,     1,     3,     1,     2,     1,     2,     1,
       3,     4,     1,     2,     1,     4,     3,     2,     1,     2,
       3,     2,     1,     1,     2,     2,     2,     1,     1,     1,
       1,     1,     1,     3,     2,     3,     1,     2,     3,     2,
       2,     2,     3,     1,     2,     5,     7,     9,     5,     7,
       0,     1,     1,     1,     1,     1,     1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* start: declaration_list  */
#line 358 "cod.y"
                         {
	    yyparse_value = (sm_ref)(yyvsp[0].list);
	}
#line 1802 "cod.tab.c"
    break;

  case 3: /* start: compound_statement  */
#line 362 "cod.y"
                           {
	    yyparse_value = (yyvsp[0].reference);
	}
#line 1810 "cod.tab.c"
    break;

  case 4: /* primary_expression: identifier_ref  */
#line 368 "cod.y"
                       {
	    (yyval.reference) = cod_new_identifier();
	    (yyval.reference)->node.identifier.id = (yyvsp[0].info).string;
	    (yyval.reference)->node.identifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 1820 "cod.tab.c"
    break;

  case 6: /* primary_expression: LPAREN expression RPAREN  */
#line 377 "cod.y"
        { (yyval.reference) = (yyvsp[-1].reference); }
#line 1826 "cod.tab.c"
    break;

  case 8: /* postfix_expression: postfix_expression LBRACKET expression RBRACKET  */
#line 385 "cod.y"
                                                        {
	    (yyval.reference) = cod_new_element_ref();
	    (yyval.reference)->node.element_ref.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	    (yyval.reference)->node.element_ref.expression = (yyvsp[-1].reference);
	    (yyval.reference)->node.element_ref.array_ref = (yyvsp[-3].reference);
	}
#line 1837 "cod.tab.c"
    break;

  case 9: /* postfix_expression: postfix_expression DOT identifier_ref  */
#line 392 "cod.y"
                                              {
	    (yyval.reference) = cod_new_field_ref();
	    (yyval.reference)->node.field_ref.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.field_ref.lx_field = (yyvsp[0].info).string;
	    (yyval.reference)->node.field_ref.struct_ref = (yyvsp[-2].reference);
	}
#line 1848 "cod.tab.c"
    break;

  case 10: /* postfix_expression: postfix_expression LPAREN argument_expression_list RPAREN  */
#line 399 "cod.y"
                                                                  {
	    (yyval.reference) = cod_new_subroutine_call();
	    (yyval.reference)->node.subroutine_call.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	    (yyval.reference)->node.subroutine_call.arguments = (yyvsp[-1].list);
	    (yyval.reference)->node.subroutine_call.sm_func_ref = (yyvsp[-3].reference);
	}
#line 1859 "cod.tab.c"
    break;

  case 11: /* postfix_expression: postfix_expression LPAREN RPAREN  */
#line 406 "cod.y"
                                          {
	    (yyval.reference) = cod_new_subroutine_call();
	    (yyval.reference)->node.subroutine_call.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.subroutine_call.arguments = NULL;
	    (yyval.reference)->node.subroutine_call.sm_func_ref = (yyvsp[-2].reference);
	}
#line 1870 "cod.tab.c"
    break;

  case 12: /* postfix_expression: postfix_expression ARROW identifier_ref  */
#line 412 "cod.y"
                                                  {
	    (yyval.reference) = cod_new_field_ref();
	    (yyval.reference)->node.field_ref.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.field_ref.lx_field = (yyvsp[0].info).string;
	    (yyval.reference)->node.field_ref.struct_ref = (yyvsp[-2].reference);
	}
#line 1881 "cod.tab.c"
    break;

  case 13: /* postfix_expression: postfix_expression INC_OP  */
#line 418 "cod.y"
                                    {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_inc;
	    (yyval.reference)->node.operator.right = NULL;
	    (yyval.reference)->node.operator.left = (yyvsp[-1].reference);
	}
#line 1893 "cod.tab.c"
    break;

  case 14: /* postfix_expression: postfix_expression DEC_OP  */
#line 425 "cod.y"
                                    {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_dec;
	    (yyval.reference)->node.operator.right = NULL;
	    (yyval.reference)->node.operator.left = (yyvsp[-1].reference);
	}
#line 1905 "cod.tab.c"
    break;

  case 15: /* argument_expression_list: assignment_expression  */
#line 435 "cod.y"
                              {
		(yyval.list) = malloc(sizeof(struct list_struct));
		(yyval.list)->node = (yyvsp[0].reference);
		(yyval.list)->next = NULL;
	}
#line 1915 "cod.tab.c"
    break;

  case 16: /* argument_expression_list: argument_expression_list COMMA assignment_expression  */
#line 441 "cod.y"
                                                             {
	    sm_list tmp = (yyvsp[-2].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = (yyvsp[0].reference);
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 1930 "cod.tab.c"
    break;

  case 18: /* unary_expression: INC_OP unary_expression  */
#line 456 "cod.y"
                                  {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_inc;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = NULL;
	}
#line 1942 "cod.tab.c"
    break;

  case 19: /* unary_expression: DEC_OP unary_expression  */
#line 463 "cod.y"
                                  {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_dec;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = NULL;
	}
#line 1954 "cod.tab.c"
    break;

  case 20: /* unary_expression: unary_operator cast_expression  */
#line 470 "cod.y"
                                         {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = (yyvsp[-1].info).op;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = NULL;
	}
#line 1966 "cod.tab.c"
    break;

  case 21: /* unary_expression: SIZEOF unary_expression  */
#line 477 "cod.y"
                                  {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_sizeof;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = NULL;
	}
#line 1978 "cod.tab.c"
    break;

  case 22: /* unary_expression: SIZEOF LPAREN type_name RPAREN  */
#line 484 "cod.y"
                                         {
	    /* dummy up a cast to hold the sm_list of the type */
	    sm_ref cast = cod_new_cast();
	    cast->node.cast.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
	    cast->node.cast.type_spec = (yyvsp[-1].list);
	    cast->node.cast.expression = NULL;

	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_sizeof;
	    (yyval.reference)->node.operator.right = cast;
	    (yyval.reference)->node.operator.left = NULL;
	}
#line 1996 "cod.tab.c"
    break;

  case 23: /* unary_operator: ARITH_AND  */
#line 499 "cod.y"
                    {
	    (yyval.info).op = op_address;
	}
#line 2004 "cod.tab.c"
    break;

  case 24: /* unary_operator: STAR  */
#line 502 "cod.y"
               {
	    (yyval.info).op = op_deref;
	}
#line 2012 "cod.tab.c"
    break;

  case 25: /* unary_operator: PLUS  */
#line 505 "cod.y"
               {
	    (yyval.info).op = op_plus;
	}
#line 2020 "cod.tab.c"
    break;

  case 26: /* unary_operator: MINUS  */
#line 508 "cod.y"
                {
	    (yyval.info).op = op_minus;
	}
#line 2028 "cod.tab.c"
    break;

  case 27: /* unary_operator: TILDE  */
#line 511 "cod.y"
                {
	    (yyval.info).op = op_not;
	}
#line 2036 "cod.tab.c"
    break;

  case 28: /* unary_operator: BANG  */
#line 514 "cod.y"
               {
	    (yyval.info).op = op_log_neg;
	  }
#line 2044 "cod.tab.c"
    break;

  case 30: /* cast_expression: LPAREN type_name RPAREN cast_expression  */
#line 521 "cod.y"
                                                  {
	    (yyval.reference) = cod_new_cast();
	    (yyval.reference)->node.cast.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
	    (yyval.reference)->node.cast.type_spec = (yyvsp[-2].list);
	    (yyval.reference)->node.cast.expression = (yyvsp[0].reference);
	}
#line 2055 "cod.tab.c"
    break;

  case 32: /* multiplicative_expression: multiplicative_expression STAR cast_expression  */
#line 533 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_mult;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2067 "cod.tab.c"
    break;

  case 33: /* multiplicative_expression: multiplicative_expression SLASH cast_expression  */
#line 542 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_div;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2079 "cod.tab.c"
    break;

  case 34: /* multiplicative_expression: multiplicative_expression MODULUS cast_expression  */
#line 551 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_modulus;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2091 "cod.tab.c"
    break;

  case 36: /* additive_expression: additive_expression PLUS multiplicative_expression  */
#line 564 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.op = op_plus;
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2103 "cod.tab.c"
    break;

  case 37: /* additive_expression: additive_expression MINUS multiplicative_expression  */
#line 573 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_minus;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2115 "cod.tab.c"
    break;

  case 39: /* shift_expression: shift_expression LEFT_SHIFT additive_expression  */
#line 586 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_left_shift;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2127 "cod.tab.c"
    break;

  case 40: /* shift_expression: shift_expression RIGHT_SHIFT additive_expression  */
#line 595 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_right_shift;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2139 "cod.tab.c"
    break;

  case 42: /* relational_expression: relational_expression LT shift_expression  */
#line 608 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_lt;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2151 "cod.tab.c"
    break;

  case 43: /* relational_expression: relational_expression GT shift_expression  */
#line 617 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_gt;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2163 "cod.tab.c"
    break;

  case 44: /* relational_expression: relational_expression LEQ shift_expression  */
#line 626 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_leq;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2175 "cod.tab.c"
    break;

  case 45: /* relational_expression: relational_expression GEQ shift_expression  */
#line 635 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_geq;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2187 "cod.tab.c"
    break;

  case 47: /* equality_expression: equality_expression EQ relational_expression  */
#line 648 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_eq;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2199 "cod.tab.c"
    break;

  case 48: /* equality_expression: equality_expression NEQ relational_expression  */
#line 657 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_neq;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2211 "cod.tab.c"
    break;

  case 50: /* and_expression: and_expression ARITH_AND equality_expression  */
#line 670 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_arith_and;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2223 "cod.tab.c"
    break;

  case 52: /* exclusive_or_expression: exclusive_or_expression ARITH_XOR and_expression  */
#line 683 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_arith_xor;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2235 "cod.tab.c"
    break;

  case 54: /* inclusive_or_expression: inclusive_or_expression ARITH_OR exclusive_or_expression  */
#line 696 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_arith_or;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2247 "cod.tab.c"
    break;

  case 56: /* logical_and_expression: logical_and_expression LOG_AND inclusive_or_expression  */
#line 709 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_log_and;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2259 "cod.tab.c"
    break;

  case 58: /* logical_or_expression: logical_or_expression LOG_OR logical_and_expression  */
#line 722 "cod.y"
        {
	    (yyval.reference) = cod_new_operator();
	    (yyval.reference)->node.operator.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.operator.op = op_log_or;
	    (yyval.reference)->node.operator.right = (yyvsp[0].reference);
	    (yyval.reference)->node.operator.left = (yyvsp[-2].reference);
	}
#line 2271 "cod.tab.c"
    break;

  case 60: /* conditional_expression: logical_or_expression QUESTION expression COLON conditional_expression  */
#line 735 "cod.y"
        {
	    (yyval.reference) = cod_new_conditional_operator();
	    (yyval.reference)->node.conditional_operator.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
	    (yyval.reference)->node.conditional_operator.condition = (yyvsp[-4].reference);
	    (yyval.reference)->node.conditional_operator.e1 = (yyvsp[-2].reference);
	    (yyval.reference)->node.conditional_operator.e2 = (yyvsp[0].reference);
	}
#line 2283 "cod.tab.c"
    break;

  case 61: /* assignment_operator: ASSIGN  */
#line 746 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_eq;}
#line 2289 "cod.tab.c"
    break;

  case 62: /* assignment_operator: MUL_ASSIGN  */
#line 748 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_mult;}
#line 2295 "cod.tab.c"
    break;

  case 63: /* assignment_operator: DIV_ASSIGN  */
#line 750 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_div;}
#line 2301 "cod.tab.c"
    break;

  case 64: /* assignment_operator: MOD_ASSIGN  */
#line 752 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_modulus;}
#line 2307 "cod.tab.c"
    break;

  case 65: /* assignment_operator: ADD_ASSIGN  */
#line 754 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_plus;}
#line 2313 "cod.tab.c"
    break;

  case 66: /* assignment_operator: SUB_ASSIGN  */
#line 756 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_minus;}
#line 2319 "cod.tab.c"
    break;

  case 67: /* assignment_operator: LEFT_ASSIGN  */
#line 758 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_left_shift;}
#line 2325 "cod.tab.c"
    break;

  case 68: /* assignment_operator: RIGHT_ASSIGN  */
#line 760 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_right_shift;}
#line 2331 "cod.tab.c"
    break;

  case 69: /* assignment_operator: AND_ASSIGN  */
#line 762 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_arith_and;}
#line 2337 "cod.tab.c"
    break;

  case 70: /* assignment_operator: XOR_ASSIGN  */
#line 764 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_arith_xor;}
#line 2343 "cod.tab.c"
    break;

  case 71: /* assignment_operator: OR_ASSIGN  */
#line 766 "cod.y"
        { (yyval.info) = (yyvsp[0].info); (yyval.info).op = op_arith_or;}
#line 2349 "cod.tab.c"
    break;

  case 72: /* assignment_expression: conditional_expression  */
#line 771 "cod.y"
        { (yyval.reference) = (yyvsp[0].reference);}
#line 2355 "cod.tab.c"
    break;

  case 73: /* assignment_expression: unary_expression assignment_operator assignment_expression  */
#line 774 "cod.y"
        {
	    (yyval.reference) = cod_new_assignment_expression();
	    (yyval.reference)->node.assignment_expression.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.assignment_expression.left = (yyvsp[-2].reference);
	    (yyval.reference)->node.assignment_expression.right = (yyvsp[0].reference);
	    (yyval.reference)->node.assignment_expression.op = (yyvsp[-1].info).op;
	}
#line 2367 "cod.tab.c"
    break;

  case 74: /* expression: assignment_expression  */
#line 785 "cod.y"
            {(yyval.reference) = (yyvsp[0].reference);}
#line 2373 "cod.tab.c"
    break;

  case 75: /* expression: expression COMMA assignment_expression  */
#line 787 "cod.y"
        {
	    (yyval.reference) = cod_new_comma_expression();
	    (yyval.reference)->node.comma_expression.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    (yyval.reference)->node.comma_expression.left = (yyvsp[-2].reference);
	    (yyval.reference)->node.comma_expression.right = (yyvsp[0].reference);
	}
#line 2384 "cod.tab.c"
    break;

  case 77: /* init_declarator_list: init_declarator  */
#line 800 "cod.y"
                          {
		(yyval.list) = malloc(sizeof(struct list_struct));
		(yyval.list)->node = (yyvsp[0].reference);
		(yyval.list)->next = NULL;
	}
#line 2394 "cod.tab.c"
    break;

  case 78: /* init_declarator_list: init_declarator_list COMMA init_declarator  */
#line 805 "cod.y"
                                                     {
	    sm_list tmp = (yyvsp[-2].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp = tmp->next;
	    tmp->node = (yyvsp[0].reference);
	    tmp->next = NULL;
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 2410 "cod.tab.c"
    break;

  case 79: /* $@1: %empty  */
#line 824 "cod.y"
             { 
		 if (parsing_type) {
		     yyparse_value = (sm_ref) (yyvsp[0].list);
		     YYACCEPT;
		 }
	     }
#line 2421 "cod.tab.c"
    break;

  case 80: /* @2: %empty  */
#line 831 "cod.y"
            {  /* stop here if we're just doing a proc decl */
		if (parsing_param_spec) {
		    (yyval.reference) = (yyvsp[0].list)->node;
		    if ((yyval.reference)->node_type == cod_declaration) {
			if  ((yyval.reference)->node.declaration.type_spec == NULL) {
			    (yyval.reference)->node.declaration.type_spec = (yyvsp[-2].list);
			} else {
			    /* 
			     * the pointer type list (with the declarator)
			     * goes at the end 
			     */
			    sm_list tmp = (yyvsp[-2].list);
			    while (tmp->next != NULL) {
				tmp = tmp->next;
			    }
			    tmp->next = (yyval.reference)->node.declaration.type_spec;
			    (yyval.reference)->node.declaration.type_spec = (yyvsp[-2].list);
			}
		    } else {
		        printf("unexpected node in init_declarator\n");
			cod_print((yyval.reference));
		    }
		    yyparse_value = (yyvsp[0].list)->node;
		    free((yyvsp[0].list));
		    YYACCEPT;
		}
	    }
#line 2453 "cod.tab.c"
    break;

  case 81: /* declaration: declaration_specifiers $@1 init_declarator_list @2 SEMI  */
#line 859 "cod.y"
            {
		(yyval.list) = (yyvsp[-2].list);
		sm_list dtmp = (yyvsp[-2].list);
		while (dtmp) {
		    sm_list type_spec;
		    if (dtmp->next != NULL) {
			type_spec = cod_dup_list((yyvsp[-4].list));
		    } else {
			type_spec = (yyvsp[-4].list);
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
		(void)(yyvsp[-1].reference);
	    }
#line 2513 "cod.tab.c"
    break;

  case 82: /* declaration: declaration_specifiers SEMI  */
#line 914 "cod.y"
                                      {
	    (yyval.list) = (yyvsp[-1].list);
	}
#line 2521 "cod.tab.c"
    break;

  case 83: /* declaration_specifiers: storage_class_specifier  */
#line 920 "cod.y"
                                  {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2531 "cod.tab.c"
    break;

  case 84: /* declaration_specifiers: storage_class_specifier declaration_specifiers  */
#line 925 "cod.y"
                                                         {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[-1].reference);
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = tmp;
	}
#line 2542 "cod.tab.c"
    break;

  case 85: /* declaration_specifiers: type_specifier  */
#line 931 "cod.y"
                         {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2552 "cod.tab.c"
    break;

  case 86: /* declaration_specifiers: type_specifier declaration_specifiers  */
#line 936 "cod.y"
                                                {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[-1].reference);
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = tmp;
	}
#line 2563 "cod.tab.c"
    break;

  case 87: /* declaration_specifiers: type_qualifier  */
#line 942 "cod.y"
                         {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2573 "cod.tab.c"
    break;

  case 88: /* declaration_specifiers: type_qualifier declaration_specifiers  */
#line 947 "cod.y"
                                                {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[-1].reference);
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = tmp;
	}
#line 2584 "cod.tab.c"
    break;

  case 90: /* init_declarator: declarator ASSIGN initializer  */
#line 958 "cod.y"
            {
		if ((yyvsp[-2].reference)->node_type == cod_declaration) {
		    (yyvsp[-2].reference)->node.declaration.init_value = (yyvsp[0].reference);
		} else if ((yyvsp[-2].reference)->node_type == cod_array_type_decl) {
		    sm_ref tmp = (yyvsp[-2].reference)->node.array_type_decl.element_ref;
		    while (tmp->node_type == cod_array_type_decl) {
			tmp = tmp->node.array_type_decl.element_ref;
		    }
		    assert(tmp->node_type == cod_declaration);
		    tmp->node.declaration.init_value = (yyvsp[0].reference);
		}
	    }
#line 2601 "cod.tab.c"
    break;

  case 91: /* storage_class_specifier: TYPEDEF  */
#line 973 "cod.y"
                  {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = TYPEDEF;
	}
#line 2611 "cod.tab.c"
    break;

  case 92: /* storage_class_specifier: STATIC  */
#line 978 "cod.y"
                 {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = STATIC;
	}
#line 2621 "cod.tab.c"
    break;

  case 93: /* storage_class_specifier: EXTERN_TOKEN  */
#line 983 "cod.y"
                       {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = EXTERN_TOKEN;
	}
#line 2631 "cod.tab.c"
    break;

  case 94: /* type_specifier: CHAR  */
#line 992 "cod.y"
             {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = CHAR;
	}
#line 2641 "cod.tab.c"
    break;

  case 95: /* type_specifier: SHORT  */
#line 997 "cod.y"
                {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = SHORT;
	}
#line 2651 "cod.tab.c"
    break;

  case 96: /* type_specifier: INT  */
#line 1002 "cod.y"
              {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = INT;
	}
#line 2661 "cod.tab.c"
    break;

  case 97: /* type_specifier: LONG  */
#line 1007 "cod.y"
               {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = LONG;
	}
#line 2671 "cod.tab.c"
    break;

  case 98: /* type_specifier: FLOAT  */
#line 1012 "cod.y"
                {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = FLOAT;
	}
#line 2681 "cod.tab.c"
    break;

  case 99: /* type_specifier: DOUBLE  */
#line 1017 "cod.y"
                 {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = DOUBLE;
	}
#line 2691 "cod.tab.c"
    break;

  case 100: /* type_specifier: VOID  */
#line 1022 "cod.y"
               {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = VOID;
	}
#line 2701 "cod.tab.c"
    break;

  case 101: /* type_specifier: SIGNED  */
#line 1027 "cod.y"
                 {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = SIGNED;
	}
#line 2711 "cod.tab.c"
    break;

  case 102: /* type_specifier: UNSIGNED  */
#line 1032 "cod.y"
                   {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = UNSIGNED;
	}
#line 2721 "cod.tab.c"
    break;

  case 103: /* type_specifier: STRING  */
#line 1037 "cod.y"
                 {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = STRING;
	}
#line 2731 "cod.tab.c"
    break;

  case 104: /* type_specifier: type_identifier  */
#line 1042 "cod.y"
                          {
	    (yyval.reference) = cod_new_identifier();
	    (yyval.reference)->node.identifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.identifier.id = (yyvsp[0].info).string;
	}
#line 2741 "cod.tab.c"
    break;

  case 105: /* type_specifier: struct_or_union_specifier  */
#line 1047 "cod.y"
                                    {
	    (yyval.reference) = (yyvsp[0].reference);
	}
#line 2749 "cod.tab.c"
    break;

  case 106: /* type_specifier: enum_specifier  */
#line 1050 "cod.y"
                         {
	    (yyval.reference) = (yyvsp[0].reference);
	}
#line 2757 "cod.tab.c"
    break;

  case 107: /* struct_or_union_specifier: struct_or_union identifier_ref LCURLY struct_declaration_list RCURLY  */
#line 1056 "cod.y"
                                                                               {
	    (yyval.reference) = cod_build_parsed_type_node(yycontext, (yyvsp[-3].info).string, (yyvsp[-1].list));
	}
#line 2765 "cod.tab.c"
    break;

  case 108: /* struct_or_union_specifier: struct_or_union LCURLY struct_declaration_list RCURLY  */
#line 1059 "cod.y"
                                                                {
	    (yyval.reference) = cod_build_parsed_type_node(yycontext, strdup("anon"), (yyvsp[-1].list));
	}
#line 2773 "cod.tab.c"
    break;

  case 109: /* struct_or_union_specifier: struct_or_union identifier_ref  */
#line 1062 "cod.y"
                                         {
	    (yyval.reference) = cod_build_parsed_type_node(yycontext, (yyvsp[0].info).string, NULL);
	}
#line 2781 "cod.tab.c"
    break;

  case 111: /* struct_or_union: UNION  */
#line 1069 "cod.y"
                {
            yyerror("UNIONs not supported!");
	}
#line 2789 "cod.tab.c"
    break;

  case 113: /* struct_declaration_list: struct_declaration_list struct_declaration  */
#line 1076 "cod.y"
                                                     {
	    sm_list tmp = (yyvsp[-1].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next =(yyvsp[0].list);
	    (yyval.list) = (yyvsp[-1].list);
	}
#line 2802 "cod.tab.c"
    break;

  case 114: /* struct_declaration: specifier_qualifier_list SEMI  */
#line 1088 "cod.y"
                                        { }
#line 2808 "cod.tab.c"
    break;

  case 115: /* struct_declaration: specifier_qualifier_list struct_declarator_list SEMI  */
#line 1089 "cod.y"
                                                               {
	    sm_list type_spec = (yyvsp[-2].list);
	    sm_list decl_list = (yyvsp[-1].list);
 	    (yyval.list) = (yyvsp[-1].list);
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
			sm_list tmp = (yyvsp[-2].list);
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
#line 2860 "cod.tab.c"
    break;

  case 116: /* struct_declarator_list: struct_declarator  */
#line 1140 "cod.y"
                            {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2870 "cod.tab.c"
    break;

  case 117: /* struct_declarator_list: struct_declarator_list COMMA struct_declarator  */
#line 1145 "cod.y"
                                                         {
	    sm_list tmp = (yyvsp[-2].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = (yyvsp[0].reference);
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 2885 "cod.tab.c"
    break;

  case 119: /* specifier_qualifier_list: type_specifier specifier_qualifier_list  */
#line 1160 "cod.y"
                                                  {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[-1].reference);
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = tmp;
	}
#line 2896 "cod.tab.c"
    break;

  case 120: /* specifier_qualifier_list: type_specifier  */
#line 1166 "cod.y"
                         {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2906 "cod.tab.c"
    break;

  case 121: /* specifier_qualifier_list: type_qualifier specifier_qualifier_list  */
#line 1171 "cod.y"
                                                  {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[-1].reference);
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = tmp;
	}
#line 2917 "cod.tab.c"
    break;

  case 122: /* specifier_qualifier_list: type_qualifier  */
#line 1177 "cod.y"
                         {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 2927 "cod.tab.c"
    break;

  case 123: /* enum_specifier: ENUM LCURLY enumerator_list RCURLY  */
#line 1185 "cod.y"
                                             {
	    (yyval.reference) = cod_new_enum_type_decl();
	    (yyval.reference)->node.enum_type_decl.id = gen_anon();
	    (yyval.reference)->node.enum_type_decl.enums = (yyvsp[-1].list);
	    (yyval.reference)->node.enum_type_decl.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
#line 2939 "cod.tab.c"
    break;

  case 124: /* enum_specifier: ENUM LCURLY enumerator_list COMMA RCURLY  */
#line 1192 "cod.y"
                                                   {
	    (yyval.reference) = cod_new_enum_type_decl();
	    (yyval.reference)->node.enum_type_decl.id = gen_anon();
	    (yyval.reference)->node.enum_type_decl.enums = (yyvsp[-2].list);
	    (yyval.reference)->node.enum_type_decl.lx_srcpos = (yyvsp[-4].info).lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
#line 2951 "cod.tab.c"
    break;

  case 125: /* enum_specifier: ENUM identifier_ref LCURLY enumerator_list RCURLY  */
#line 1199 "cod.y"
                                                            {
	    (yyval.reference) = cod_new_enum_type_decl();
	    (yyval.reference)->node.enum_type_decl.id = (yyvsp[-3].info).string;
	    (yyval.reference)->node.enum_type_decl.enums = (yyvsp[-1].list);
	    (yyval.reference)->node.enum_type_decl.lx_srcpos = (yyvsp[-4].info).lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
#line 2963 "cod.tab.c"
    break;

  case 126: /* enum_specifier: ENUM identifier_ref LCURLY enumerator_list COMMA RCURLY  */
#line 1206 "cod.y"
                                                                  {
	    (yyval.reference) = cod_new_enum_type_decl();
	    (yyval.reference)->node.enum_type_decl.id = (yyvsp[-4].info).string;
	    (yyval.reference)->node.enum_type_decl.enums = (yyvsp[-2].list);
	    (yyval.reference)->node.enum_type_decl.lx_srcpos = (yyvsp[-5].info).lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
#line 2975 "cod.tab.c"
    break;

  case 127: /* enum_specifier: ENUM identifier_ref  */
#line 1213 "cod.y"
                              {
	    (yyval.reference) = cod_new_enum_type_decl();
	    (yyval.reference)->node.enum_type_decl.id = (yyvsp[0].info).string;
	    (yyval.reference)->node.enum_type_decl.enums = NULL;
	    (yyval.reference)->node.enum_type_decl.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    // cod_add_defined_type(decl->node.declaration.id, yycontext);
	}
#line 2987 "cod.tab.c"
    break;

  case 128: /* enumerator_list: enumerator  */
#line 1223 "cod.y"
                     {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[0].reference);
	    tmp->next = NULL;
	    (yyval.list) = tmp;
	}
#line 2998 "cod.tab.c"
    break;

  case 129: /* enumerator_list: enumerator_list COMMA enumerator  */
#line 1229 "cod.y"
                                           {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[0].reference);
	    tmp->next = (yyvsp[-2].list);
	    (yyval.list) = tmp;
	}
#line 3009 "cod.tab.c"
    break;

  case 130: /* enumerator: identifier_ref ASSIGN constant_expression  */
#line 1239 "cod.y"
                                                    {
	    (yyval.reference) = cod_new_enumerator();
	    (yyval.reference)->node.enumerator.id = (yyvsp[-2].info).string;
	    (yyval.reference)->node.enumerator.const_expression = (yyvsp[0].reference);
	}
#line 3019 "cod.tab.c"
    break;

  case 131: /* enumerator: identifier_ref  */
#line 1245 "cod.y"
                         {
	    (yyval.reference) = cod_new_enumerator();
	    (yyval.reference)->node.enumerator.id = (yyvsp[0].info).string;
	    (yyval.reference)->node.enumerator.const_expression = NULL;
	}
#line 3029 "cod.tab.c"
    break;

  case 132: /* type_qualifier: CONST  */
#line 1253 "cod.y"
                {
	    (yyval.reference) = cod_new_type_specifier();
	    (yyval.reference)->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    (yyval.reference)->node.type_specifier.token = CONST;
	}
#line 3039 "cod.tab.c"
    break;

  case 134: /* declarator: pointer direct_declarator  */
#line 1263 "cod.y"
                                  {
	    (yyval.reference) = (yyvsp[0].reference);
	    if ((yyval.reference)->node_type == cod_declaration) {
		(yyval.reference)->node.declaration.type_spec = (yyvsp[-1].list);
	    } else if ((yyval.reference)->node_type == cod_array_type_decl) {
		(yyval.reference)->node.array_type_decl.type_spec = (yyvsp[-1].list);
	    } else {
		printf("Unknown direct_declarator entry\n");
		cod_print((yyval.reference));
	    }
	}
#line 3055 "cod.tab.c"
    break;

  case 135: /* direct_declarator: identifier_ref  */
#line 1277 "cod.y"
                            {
		(yyval.reference) = cod_new_declaration();
		(yyval.reference)->node.declaration.param_num = -1;
		(yyval.reference)->node.declaration.id = (yyvsp[0].info).string;
		(yyval.reference)->node.declaration.init_value = NULL;
		(yyval.reference)->node.declaration.lx_srcpos = (yyvsp[0].info).lx_srcpos;
		(yyval.reference)->node.declaration.is_subroutine = 0;
		(yyval.reference)->node.declaration.params = NULL;
	    }
#line 3069 "cod.tab.c"
    break;

  case 136: /* direct_declarator: LPAREN declarator RPAREN  */
#line 1286 "cod.y"
                                   {
	    (yyval.reference) = (yyvsp[-1].reference);
	}
#line 3077 "cod.tab.c"
    break;

  case 137: /* direct_declarator: identifier_ref LPAREN parameter_type_list RPAREN  */
#line 1289 "cod.y"
                                                           {
		(yyval.reference) = cod_new_declaration();
		(yyval.reference)->node.declaration.param_num = -1;
		(yyval.reference)->node.declaration.id = (yyvsp[-3].info).string;
		(yyval.reference)->node.declaration.init_value = NULL;
		(yyval.reference)->node.declaration.lx_srcpos = (yyvsp[-3].info).lx_srcpos;
		(yyval.reference)->node.declaration.is_subroutine = 1;
		(yyval.reference)->node.declaration.params = (yyvsp[-1].list);
	}
#line 3091 "cod.tab.c"
    break;

  case 138: /* direct_declarator: identifier_ref LPAREN RPAREN  */
#line 1298 "cod.y"
                                       {
		(yyval.reference) = cod_new_declaration();
		(yyval.reference)->node.declaration.param_num = -1;
		(yyval.reference)->node.declaration.id = (yyvsp[-2].info).string;
		(yyval.reference)->node.declaration.init_value = NULL;
		(yyval.reference)->node.declaration.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
		(yyval.reference)->node.declaration.is_subroutine = 1;
		(yyval.reference)->node.declaration.params = NULL;
	}
#line 3105 "cod.tab.c"
    break;

  case 139: /* direct_declarator: direct_declarator LBRACKET constant_expression RBRACKET  */
#line 1307 "cod.y"
                                                                  {
		(yyval.reference) = cod_new_array_type_decl();
		(yyval.reference)->node.array_type_decl.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
		(yyval.reference)->node.array_type_decl.size_expr = (yyvsp[-1].reference);
		(yyval.reference)->node.array_type_decl.element_ref = (yyvsp[-3].reference);
		(yyval.reference)->node.array_type_decl.sm_dynamic_size = NULL;
	}
#line 3117 "cod.tab.c"
    break;

  case 140: /* direct_declarator: direct_declarator LBRACKET RBRACKET  */
#line 1314 "cod.y"
                                              {
		(yyval.reference) = cod_new_array_type_decl();
		(yyval.reference)->node.array_type_decl.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
		(yyval.reference)->node.array_type_decl.size_expr = NULL;
		(yyval.reference)->node.array_type_decl.element_ref = (yyvsp[-2].reference);
		(yyval.reference)->node.array_type_decl.sm_dynamic_size = NULL;
	}
#line 3129 "cod.tab.c"
    break;

  case 141: /* pointer: STAR  */
#line 1324 "cod.y"
             {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = NULL;
	}
#line 3142 "cod.tab.c"
    break;

  case 142: /* pointer: STAR type_qualifier_list  */
#line 1332 "cod.y"
                                   {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[0].list);
	}
#line 3155 "cod.tab.c"
    break;

  case 143: /* pointer: STAR pointer  */
#line 1340 "cod.y"
                       {
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    star->node.type_specifier.token = STAR;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[0].list);
	}
#line 3168 "cod.tab.c"
    break;

  case 144: /* pointer: STAR type_qualifier_list pointer  */
#line 1348 "cod.y"
                                           {
	    sm_list tmp = (yyvsp[-1].list);
	    sm_ref star = cod_new_type_specifier();
	    star->node.type_specifier.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	    star->node.type_specifier.token = STAR;

	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[-1].list);
	}
#line 3187 "cod.tab.c"
    break;

  case 145: /* pointer: AT  */
#line 1362 "cod.y"
             {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) { 
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	    star->node.type_specifier.token = AT;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = NULL;
	}
#line 3203 "cod.tab.c"
    break;

  case 146: /* pointer: AT type_qualifier_list  */
#line 1373 "cod.y"
                                 {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    star->node.type_specifier.token = AT;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[0].list);
	}
#line 3219 "cod.tab.c"
    break;

  case 147: /* pointer: AT pointer  */
#line 1384 "cod.y"
                     {
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	    star->node.type_specifier.token = AT;
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[0].list);
	}
#line 3235 "cod.tab.c"
    break;

  case 148: /* pointer: AT type_qualifier_list pointer  */
#line 1395 "cod.y"
                                         {
	    sm_list tmp = (yyvsp[-1].list);
	    sm_ref star = cod_new_type_specifier();
	    if(!cod_segmented_pointers) {
                yyerror("Segmented pointers disabled!");
	    }
	    star->node.type_specifier.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	    star->node.type_specifier.token = AT;

	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = star;
	    (yyval.list)->next = (yyvsp[-1].list);
	}
#line 3257 "cod.tab.c"
    break;

  case 149: /* type_qualifier_list: type_qualifier  */
#line 1415 "cod.y"
                         {
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = (yyvsp[0].reference);
	    (yyval.list)->next = NULL;
	}
#line 3267 "cod.tab.c"
    break;

  case 150: /* type_qualifier_list: type_qualifier_list type_qualifier  */
#line 1420 "cod.y"
                                             {
	    sm_list tmp = (yyvsp[-1].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = (yyvsp[0].reference);
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-1].list);
	}
#line 3282 "cod.tab.c"
    break;

  case 152: /* parameter_type_list: parameter_list COMMA DOTDOTDOT  */
#line 1434 "cod.y"
                                       {
	    sm_list tmp = (yyvsp[-2].list);
	    sm_ref id = cod_new_declaration();
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = id;
	    tmp->next->next = NULL;
	    id->node.declaration.id = strdup("...");
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 3299 "cod.tab.c"
    break;

  case 153: /* parameter_list: parameter_declaration  */
#line 1449 "cod.y"
                              {
		(yyval.list) = malloc(sizeof(struct list_struct));
		(yyval.list)->node = (yyvsp[0].reference);
		(yyval.list)->next = NULL;
	}
#line 3309 "cod.tab.c"
    break;

  case 154: /* parameter_list: parameter_list COMMA parameter_declaration  */
#line 1455 "cod.y"
                                                   {
	    sm_list tmp = (yyvsp[-2].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = (yyvsp[0].reference);
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 3324 "cod.tab.c"
    break;

  case 155: /* parameter_declaration: declaration_specifiers  */
#line 1470 "cod.y"
                                 {
	    (yyval.reference) = cod_new_declaration();
	    (yyval.reference)->node.declaration.param_num = -1;
	    (yyval.reference)->node.declaration.id = gen_anon();
	    (yyval.reference)->node.declaration.init_value = NULL;
	    (yyval.reference)->node.declaration.is_subroutine = 0;
	    (yyval.reference)->node.declaration.params = NULL;
	    (yyval.reference)->node.declaration.type_spec = (yyvsp[0].list);
	}
#line 3338 "cod.tab.c"
    break;

  case 156: /* parameter_declaration: declaration_specifiers declarator  */
#line 1479 "cod.y"
                                            {
		(yyval.reference) = (yyvsp[0].reference);
		if ((yyval.reference)->node_type == cod_declaration) {
		    (yyval.reference)->node.declaration.static_var = 0;
		    if  ((yyval.reference)->node.declaration.type_spec == NULL) {
		        (yyval.reference)->node.declaration.type_spec = (yyvsp[-1].list);
		    } else {
		        /* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
		      sm_list tmp = (yyvsp[-1].list);
		      while (tmp->next != NULL) {
			  tmp = tmp->next;
		      }
		      tmp->next = (yyval.reference)->node.declaration.type_spec;
		      (yyval.reference)->node.declaration.type_spec = (yyvsp[-1].list);
		    }
		} else if ((yyval.reference)->node_type == cod_array_type_decl) {
		    if  ((yyval.reference)->node.array_type_decl.type_spec == NULL) {
		        (yyval.reference)->node.array_type_decl.type_spec = (yyvsp[-1].list);
		    } else {
		        /* 
			 * the pointer type list (with the declarator)
			 * goes at the end 
			 */
		      sm_list tmp = (yyvsp[-1].list);
		      while (tmp->next != NULL) {
			  tmp = tmp->next;
		      }
		      tmp->next = (yyval.reference)->node.array_type_decl.type_spec;
		      (yyval.reference)->node.array_type_decl.type_spec = (yyvsp[-1].list);
		    }
		} else {
		    printf("unexpected node in parameter_declaration");
		}
	}
#line 3380 "cod.tab.c"
    break;

  case 158: /* type_name: specifier_qualifier_list abstract_declarator  */
#line 1519 "cod.y"
                                                       {
	    sm_list tmp = (yyvsp[-1].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = (yyvsp[0].list);
	    (yyval.list) = (yyvsp[-1].list);
	}
#line 3393 "cod.tab.c"
    break;

  case 160: /* initializer: LCURLY initializer_list RCURLY  */
#line 1570 "cod.y"
        { 
	    (yyval.reference) = cod_new_initializer_list();
	    (yyval.reference)->node.initializer_list.initializers = (yyvsp[-1].list);
	}
#line 3402 "cod.tab.c"
    break;

  case 161: /* initializer: LCURLY initializer_list COMMA RCURLY  */
#line 1575 "cod.y"
        { 
	    (yyval.reference) = cod_new_initializer_list();
	    (yyval.reference)->node.initializer_list.initializers = (yyvsp[-2].list);
	}
#line 3411 "cod.tab.c"
    break;

  case 162: /* initializer: assignment_expression  */
#line 1579 "cod.y"
                                { (yyval.reference) = (yyvsp[0].reference);}
#line 3417 "cod.tab.c"
    break;

  case 163: /* initializer_list: designation initializer  */
#line 1584 "cod.y"
                                {
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = (yyvsp[-1].list);
	    initializer->node.initializer.initializer = (yyvsp[0].reference);
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = initializer;
	    (yyval.list)->next = NULL;
	}
#line 3430 "cod.tab.c"
    break;

  case 164: /* initializer_list: initializer  */
#line 1592 "cod.y"
                      {
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = NULL;
	    initializer->node.initializer.initializer = (yyvsp[0].reference);
	    (yyval.list) = malloc(sizeof(struct list_struct));
	    (yyval.list)->node = initializer;
	    (yyval.list)->next = NULL;
	}
#line 3443 "cod.tab.c"
    break;

  case 165: /* initializer_list: initializer_list COMMA designation initializer  */
#line 1600 "cod.y"
                                                         {
	    sm_list tmp = (yyvsp[-3].list);
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = (yyvsp[-1].list);
	    initializer->node.initializer.initializer = (yyvsp[0].reference);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = initializer;
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-3].list);
	}
#line 3461 "cod.tab.c"
    break;

  case 166: /* initializer_list: initializer_list COMMA initializer  */
#line 1613 "cod.y"
                                             {
	    sm_list tmp = (yyvsp[-2].list);
	    sm_ref initializer = cod_new_initializer();
	    initializer->node.initializer.designation = NULL;
	    initializer->node.initializer.initializer = (yyvsp[0].reference);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = initializer;
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-2].list);
	}
#line 3479 "cod.tab.c"
    break;

  case 167: /* designation: designator_list ASSIGN  */
#line 1630 "cod.y"
        { (yyval.list) = (yyvsp[-1].list);}
#line 3485 "cod.tab.c"
    break;

  case 168: /* designator_list: designator  */
#line 1634 "cod.y"
                     {
		(yyval.list) = malloc(sizeof(struct list_struct));
		(yyval.list)->node = (yyvsp[0].reference);
		(yyval.list)->next = NULL;
	}
#line 3495 "cod.tab.c"
    break;

  case 169: /* designator_list: designator_list designator  */
#line 1639 "cod.y"
                                     {
	    sm_list tmp = (yyvsp[-1].list);
	    while (tmp->next != NULL) {
		tmp = tmp->next;
	    }
	    tmp->next = malloc(sizeof(struct list_struct));
	    tmp->next->node = (yyvsp[0].reference);
	    tmp->next->next = NULL;
	    (yyval.list) = (yyvsp[-1].list);
	}
#line 3510 "cod.tab.c"
    break;

  case 170: /* designator: LBRACKET constant_expression RBRACKET  */
#line 1653 "cod.y"
        { 
	    (yyval.reference) = cod_new_designator();
	    (yyval.reference)->node.designator.expression = (yyvsp[-1].reference);
	    (yyval.reference)->node.designator.id = NULL;
	}
#line 3520 "cod.tab.c"
    break;

  case 171: /* designator: DOT identifier_ref  */
#line 1659 "cod.y"
        { 
	    (yyval.reference) = cod_new_designator();
	    (yyval.reference)->node.designator.expression = NULL;
	    (yyval.reference)->node.designator.id = (yyvsp[0].info).string;
	}
#line 3530 "cod.tab.c"
    break;

  case 172: /* decls_stmts_list: statement  */
#line 1667 "cod.y"
                  {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[0].reference);
	    tmp->next = NULL;
	    (yyval.list) = tmp;
	}
#line 3541 "cod.tab.c"
    break;

  case 173: /* decls_stmts_list: declaration  */
#line 1673 "cod.y"
                       {
	    (yyval.list) = (yyvsp[0].list);
	   }
#line 3549 "cod.tab.c"
    break;

  case 174: /* decls_stmts_list: error SEMI  */
#line 1676 "cod.y"
                     {
	      (yyval.list) = NULL;
	  }
#line 3557 "cod.tab.c"
    break;

  case 175: /* decls_stmts_list: decls_stmts_list statement  */
#line 1679 "cod.y"
                                     {
	    sm_list tmp = malloc(sizeof(struct list_struct));
	    tmp->node = (yyvsp[0].reference);
	    tmp->next = NULL;
	    (yyval.list) = cod_append_list((yyvsp[-1].list), tmp);
	}
#line 3568 "cod.tab.c"
    break;

  case 176: /* decls_stmts_list: decls_stmts_list declaration  */
#line 1685 "cod.y"
                                       {
	    (yyval.list) = cod_append_list((yyvsp[-1].list), (yyvsp[0].list));
	}
#line 3576 "cod.tab.c"
    break;

  case 183: /* labeled_statement: identifier_ref COLON statement  */
#line 1702 "cod.y"
                                       {
	    (yyval.reference) = cod_new_label_statement();
	    (yyval.reference)->node.label_statement.name =  (yyvsp[-2].info).string;
	    (yyval.reference)->node.label_statement.statement = (yyvsp[0].reference);
	}
#line 3586 "cod.tab.c"
    break;

  case 184: /* compound_statement: LCURLY RCURLY  */
#line 1709 "cod.y"
                      {
	    (yyval.reference) = cod_new_compound_statement();
	}
#line 3594 "cod.tab.c"
    break;

  case 185: /* compound_statement: LCURLY decls_stmts_list RCURLY  */
#line 1712 "cod.y"
                                         {
	    int count = (yyvsp[-2].info).type_stack_count;
	    (yyval.reference) = cod_new_compound_statement();
	    (yyval.reference)->node.compound_statement.decls = (yyvsp[-1].list);
	    cod_remove_defined_types(yycontext, count);
	}
#line 3605 "cod.tab.c"
    break;

  case 186: /* declaration_list: declaration  */
#line 1720 "cod.y"
                    { (yyval.list) = (yyvsp[0].list); }
#line 3611 "cod.tab.c"
    break;

  case 187: /* declaration_list: declaration_list declaration  */
#line 1722 "cod.y"
                                     {
	    if ((yyvsp[-1].list) == NULL) {
		(yyval.list) = (yyvsp[0].list);
	    } else {
		sm_list tmp = (yyvsp[-1].list);
		while (tmp->next != NULL) {
		    tmp = tmp->next;
		}
		tmp->next = (yyvsp[0].list);
		(yyval.list) = (yyvsp[-1].list);
	    }
	}
#line 3628 "cod.tab.c"
    break;

  case 188: /* jump_statement: RETURN_TOKEN expression SEMI  */
#line 1736 "cod.y"
                                     {
	    (yyval.reference) = cod_new_return_statement();
	    (yyval.reference)->node.return_statement.expression = (yyvsp[-1].reference);
	    (yyval.reference)->node.return_statement.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	}
#line 3638 "cod.tab.c"
    break;

  case 189: /* jump_statement: RETURN_TOKEN SEMI  */
#line 1741 "cod.y"
                            {
	    (yyval.reference) = cod_new_return_statement();
	    (yyval.reference)->node.return_statement.expression = NULL;
	    (yyval.reference)->node.return_statement.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	}
#line 3648 "cod.tab.c"
    break;

  case 190: /* jump_statement: CONTINUE SEMI  */
#line 1746 "cod.y"
                        {
	    (yyval.reference) = cod_new_jump_statement();
	    (yyval.reference)->node.jump_statement.continue_flag = 1;
	    (yyval.reference)->node.jump_statement.goto_target = NULL;
	    (yyval.reference)->node.jump_statement.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	}
#line 3659 "cod.tab.c"
    break;

  case 191: /* jump_statement: BREAK SEMI  */
#line 1752 "cod.y"
                     {
	    (yyval.reference) = cod_new_jump_statement();
	    (yyval.reference)->node.jump_statement.continue_flag = 0;
	    (yyval.reference)->node.jump_statement.goto_target = NULL;
	    (yyval.reference)->node.jump_statement.lx_srcpos = (yyvsp[-1].info).lx_srcpos;
	}
#line 3670 "cod.tab.c"
    break;

  case 192: /* jump_statement: GOTO identifier_ref SEMI  */
#line 1758 "cod.y"
                                  {
	    (yyval.reference) = cod_new_jump_statement();
	    (yyval.reference)->node.jump_statement.continue_flag = 0;
	    (yyval.reference)->node.jump_statement.goto_target = (yyvsp[-1].info).string;
	    (yyval.reference)->node.jump_statement.lx_srcpos = (yyvsp[-2].info).lx_srcpos;
	}
#line 3681 "cod.tab.c"
    break;

  case 193: /* expression_statement: SEMI  */
#line 1767 "cod.y"
             {
	    (yyval.reference) = NULL;
	}
#line 3689 "cod.tab.c"
    break;

  case 194: /* expression_statement: expression SEMI  */
#line 1771 "cod.y"
                        { 
	    (yyval.reference) = cod_new_expression_statement();
	    (yyval.reference)->node.expression_statement.expression = (yyvsp[-1].reference);
	}
#line 3698 "cod.tab.c"
    break;

  case 195: /* selection_statement: IF LPAREN expression RPAREN statement  */
#line 1782 "cod.y"
        { 
	    (yyval.reference) = cod_new_selection_statement();
	    (yyval.reference)->node.selection_statement.lx_srcpos = (yyvsp[-4].info).lx_srcpos;
	    (yyval.reference)->node.selection_statement.conditional = (yyvsp[-2].reference);
	    (yyval.reference)->node.selection_statement.then_part = (yyvsp[0].reference);
	    (yyval.reference)->node.selection_statement.else_part = NULL;
	}
#line 3710 "cod.tab.c"
    break;

  case 196: /* selection_statement: IF LPAREN expression RPAREN statement ELSE statement  */
#line 1791 "cod.y"
        { 
	    (yyval.reference) = cod_new_selection_statement();
	    (yyval.reference)->node.selection_statement.lx_srcpos = (yyvsp[-6].info).lx_srcpos;
	    (yyval.reference)->node.selection_statement.conditional = (yyvsp[-4].reference);
	    (yyval.reference)->node.selection_statement.then_part = (yyvsp[-2].reference);
	    (yyval.reference)->node.selection_statement.else_part = (yyvsp[0].reference);
	}
#line 3722 "cod.tab.c"
    break;

  case 197: /* iteration_statement: FOR LPAREN expression_opt SEMI expression_opt SEMI expression_opt RPAREN statement  */
#line 1808 "cod.y"
        { 
	    (yyval.reference) = cod_new_iteration_statement();
	    (yyval.reference)->node.iteration_statement.lx_srcpos = (yyvsp[-8].info).lx_srcpos;
	    (yyval.reference)->node.iteration_statement.init_expr = (yyvsp[-6].reference);
	    (yyval.reference)->node.iteration_statement.test_expr = (yyvsp[-4].reference);
	    (yyval.reference)->node.iteration_statement.iter_expr = (yyvsp[-2].reference);
	    (yyval.reference)->node.iteration_statement.statement = (yyvsp[0].reference);
	}
#line 3735 "cod.tab.c"
    break;

  case 198: /* iteration_statement: WHILE LPAREN expression RPAREN statement  */
#line 1818 "cod.y"
        { 
	    (yyval.reference) = cod_new_iteration_statement();
	    (yyval.reference)->node.iteration_statement.lx_srcpos = (yyvsp[-4].info).lx_srcpos;
	    (yyval.reference)->node.iteration_statement.init_expr = NULL;
	    (yyval.reference)->node.iteration_statement.test_expr = (yyvsp[-2].reference);
	    (yyval.reference)->node.iteration_statement.iter_expr = NULL;
	    (yyval.reference)->node.iteration_statement.statement = (yyvsp[0].reference);
	}
#line 3748 "cod.tab.c"
    break;

  case 199: /* iteration_statement: DO statement WHILE LPAREN expression RPAREN SEMI  */
#line 1828 "cod.y"
        { 
	    (yyval.reference) = cod_new_iteration_statement();
	    (yyval.reference)->node.iteration_statement.lx_srcpos = (yyvsp[-6].info).lx_srcpos;
	    (yyval.reference)->node.iteration_statement.init_expr = NULL;
	    (yyval.reference)->node.iteration_statement.test_expr = NULL;
	    (yyval.reference)->node.iteration_statement.post_test_expr = (yyvsp[-2].reference);
	    (yyval.reference)->node.iteration_statement.iter_expr = NULL;
	    (yyval.reference)->node.iteration_statement.statement = (yyvsp[-5].reference);
	}
#line 3762 "cod.tab.c"
    break;

  case 200: /* expression_opt: %empty  */
#line 1841 "cod.y"
        { (yyval.reference) = NULL; }
#line 3768 "cod.tab.c"
    break;

  case 202: /* constant: integer_constant  */
#line 1846 "cod.y"
                         {
	    (yyval.reference) = cod_new_constant();
	    (yyval.reference)->node.constant.token = integer_constant;
	    (yyval.reference)->node.constant.const_val = (yyvsp[0].info).string;
	    (yyval.reference)->node.constant.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 3779 "cod.tab.c"
    break;

  case 203: /* constant: floating_constant  */
#line 1853 "cod.y"
                          {
	    (yyval.reference) = cod_new_constant();
	    (yyval.reference)->node.constant.token = floating_constant;
	    (yyval.reference)->node.constant.const_val = (yyvsp[0].info).string;
	    (yyval.reference)->node.constant.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 3790 "cod.tab.c"
    break;

  case 204: /* constant: string_constant  */
#line 1860 "cod.y"
                        {
	    (yyval.reference) = cod_new_constant();
	    (yyval.reference)->node.constant.token = string_constant;
	    (yyval.reference)->node.constant.const_val = (yyvsp[0].info).string;
	    (yyval.reference)->node.constant.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 3801 "cod.tab.c"
    break;

  case 205: /* constant: character_constant  */
#line 1867 "cod.y"
                           {
	    (yyval.reference) = cod_new_constant();
	    (yyval.reference)->node.constant.token = character_constant;
	    (yyval.reference)->node.constant.const_val = (yyvsp[0].info).string;
	    (yyval.reference)->node.constant.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 3812 "cod.tab.c"
    break;

  case 206: /* constant: enumeration_constant  */
#line 1874 "cod.y"
                             {
	    (yyval.reference) = cod_new_constant();
	    (yyval.reference)->node.constant.token = character_constant;
	    (yyval.reference)->node.constant.const_val = (yyvsp[0].info).string;
	    (yyval.reference)->node.constant.lx_srcpos = (yyvsp[0].info).lx_srcpos;
	}
#line 3823 "cod.tab.c"
    break;


#line 3827 "cod.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 1882 "cod.y"

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
	    intptr_t left_val = 0, right_val = 0, value = 0;
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
	
