/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

#ifndef YY_YY_COD_TAB_H_INCLUDED
# define YY_YY_COD_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    ARROW = 258,                   /* ARROW  */
    LPAREN = 259,                  /* LPAREN  */
    RPAREN = 260,                  /* RPAREN  */
    LCURLY = 261,                  /* LCURLY  */
    RCURLY = 262,                  /* RCURLY  */
    COLON = 263,                   /* COLON  */
    QUESTION = 264,                /* QUESTION  */
    LBRACKET = 265,                /* LBRACKET  */
    RBRACKET = 266,                /* RBRACKET  */
    DOT = 267,                     /* DOT  */
    STAR = 268,                    /* STAR  */
    AT = 269,                      /* AT  */
    SLASH = 270,                   /* SLASH  */
    MODULUS = 271,                 /* MODULUS  */
    PLUS = 272,                    /* PLUS  */
    MINUS = 273,                   /* MINUS  */
    TILDE = 274,                   /* TILDE  */
    LEQ = 275,                     /* LEQ  */
    LT = 276,                      /* LT  */
    GEQ = 277,                     /* GEQ  */
    GT = 278,                      /* GT  */
    EQ = 279,                      /* EQ  */
    NEQ = 280,                     /* NEQ  */
    LEFT_SHIFT = 281,              /* LEFT_SHIFT  */
    RIGHT_SHIFT = 282,             /* RIGHT_SHIFT  */
    ASSIGN = 283,                  /* ASSIGN  */
    MUL_ASSIGN = 284,              /* MUL_ASSIGN  */
    DIV_ASSIGN = 285,              /* DIV_ASSIGN  */
    MOD_ASSIGN = 286,              /* MOD_ASSIGN  */
    ADD_ASSIGN = 287,              /* ADD_ASSIGN  */
    SUB_ASSIGN = 288,              /* SUB_ASSIGN  */
    LEFT_ASSIGN = 289,             /* LEFT_ASSIGN  */
    RIGHT_ASSIGN = 290,            /* RIGHT_ASSIGN  */
    AND_ASSIGN = 291,              /* AND_ASSIGN  */
    XOR_ASSIGN = 292,              /* XOR_ASSIGN  */
    OR_ASSIGN = 293,               /* OR_ASSIGN  */
    LOG_OR = 294,                  /* LOG_OR  */
    LOG_AND = 295,                 /* LOG_AND  */
    ARITH_OR = 296,                /* ARITH_OR  */
    ARITH_AND = 297,               /* ARITH_AND  */
    ARITH_XOR = 298,               /* ARITH_XOR  */
    INC_OP = 299,                  /* INC_OP  */
    DEC_OP = 300,                  /* DEC_OP  */
    BANG = 301,                    /* BANG  */
    SEMI = 302,                    /* SEMI  */
    IF = 303,                      /* IF  */
    ELSE = 304,                    /* ELSE  */
    FOR = 305,                     /* FOR  */
    DO = 306,                      /* DO  */
    WHILE = 307,                   /* WHILE  */
    CHAR = 308,                    /* CHAR  */
    SHORT = 309,                   /* SHORT  */
    INT = 310,                     /* INT  */
    LONG = 311,                    /* LONG  */
    UNSIGNED = 312,                /* UNSIGNED  */
    SIGNED = 313,                  /* SIGNED  */
    FLOAT = 314,                   /* FLOAT  */
    DOUBLE = 315,                  /* DOUBLE  */
    VOID = 316,                    /* VOID  */
    STRING = 317,                  /* STRING  */
    STATIC = 318,                  /* STATIC  */
    EXTERN_TOKEN = 319,            /* EXTERN_TOKEN  */
    STRUCT = 320,                  /* STRUCT  */
    ENUM = 321,                    /* ENUM  */
    UNION = 322,                   /* UNION  */
    CONST = 323,                   /* CONST  */
    SIZEOF = 324,                  /* SIZEOF  */
    TYPEDEF = 325,                 /* TYPEDEF  */
    RETURN_TOKEN = 326,            /* RETURN_TOKEN  */
    CONTINUE = 327,                /* CONTINUE  */
    BREAK = 328,                   /* BREAK  */
    GOTO = 329,                    /* GOTO  */
    PRINT = 330,                   /* PRINT  */
    COMMA = 331,                   /* COMMA  */
    DOTDOTDOT = 332,               /* DOTDOTDOT  */
    integer_constant = 333,        /* integer_constant  */
    character_constant = 334,      /* character_constant  */
    string_constant = 335,         /* string_constant  */
    floating_constant = 336,       /* floating_constant  */
    identifier_ref = 337,          /* identifier_ref  */
    type_identifier = 338,         /* type_identifier  */
    enumeration_constant = 339     /* enumeration_constant  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 201 "cod.y"

    lx_info info;
    sm_ref reference;
    operator_t operator;
    sm_list list;
    char *string;

#line 156 "cod.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_COD_TAB_H_INCLUDED  */
