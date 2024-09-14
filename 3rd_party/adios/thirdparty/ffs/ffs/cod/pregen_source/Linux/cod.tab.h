/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ARROW = 258,
     LPAREN = 259,
     RPAREN = 260,
     LCURLY = 261,
     RCURLY = 262,
     COLON = 263,
     QUESTION = 264,
     LBRACKET = 265,
     RBRACKET = 266,
     DOT = 267,
     STAR = 268,
     AT = 269,
     SLASH = 270,
     MODULUS = 271,
     PLUS = 272,
     MINUS = 273,
     TILDE = 274,
     LEQ = 275,
     LT = 276,
     GEQ = 277,
     GT = 278,
     EQ = 279,
     NEQ = 280,
     LEFT_SHIFT = 281,
     RIGHT_SHIFT = 282,
     ASSIGN = 283,
     MUL_ASSIGN = 284,
     DIV_ASSIGN = 285,
     MOD_ASSIGN = 286,
     ADD_ASSIGN = 287,
     SUB_ASSIGN = 288,
     LEFT_ASSIGN = 289,
     RIGHT_ASSIGN = 290,
     AND_ASSIGN = 291,
     XOR_ASSIGN = 292,
     OR_ASSIGN = 293,
     LOG_OR = 294,
     LOG_AND = 295,
     ARITH_OR = 296,
     ARITH_AND = 297,
     ARITH_XOR = 298,
     INC_OP = 299,
     DEC_OP = 300,
     BANG = 301,
     SEMI = 302,
     IF = 303,
     ELSE = 304,
     FOR = 305,
     DO = 306,
     WHILE = 307,
     CHAR = 308,
     SHORT = 309,
     INT = 310,
     LONG = 311,
     UNSIGNED = 312,
     SIGNED = 313,
     FLOAT = 314,
     DOUBLE = 315,
     VOID = 316,
     STRING = 317,
     STATIC = 318,
     EXTERN_TOKEN = 319,
     STRUCT = 320,
     ENUM = 321,
     UNION = 322,
     CONST = 323,
     SIZEOF = 324,
     TYPEDEF = 325,
     RETURN_TOKEN = 326,
     CONTINUE = 327,
     BREAK = 328,
     GOTO = 329,
     PRINT = 330,
     COMMA = 331,
     DOTDOTDOT = 332,
     integer_constant = 333,
     character_constant = 334,
     string_constant = 335,
     floating_constant = 336,
     identifier_ref = 337,
     type_identifier = 338,
     enumeration_constant = 339
   };
#endif
/* Tokens.  */
#define ARROW 258
#define LPAREN 259
#define RPAREN 260
#define LCURLY 261
#define RCURLY 262
#define COLON 263
#define QUESTION 264
#define LBRACKET 265
#define RBRACKET 266
#define DOT 267
#define STAR 268
#define AT 269
#define SLASH 270
#define MODULUS 271
#define PLUS 272
#define MINUS 273
#define TILDE 274
#define LEQ 275
#define LT 276
#define GEQ 277
#define GT 278
#define EQ 279
#define NEQ 280
#define LEFT_SHIFT 281
#define RIGHT_SHIFT 282
#define ASSIGN 283
#define MUL_ASSIGN 284
#define DIV_ASSIGN 285
#define MOD_ASSIGN 286
#define ADD_ASSIGN 287
#define SUB_ASSIGN 288
#define LEFT_ASSIGN 289
#define RIGHT_ASSIGN 290
#define AND_ASSIGN 291
#define XOR_ASSIGN 292
#define OR_ASSIGN 293
#define LOG_OR 294
#define LOG_AND 295
#define ARITH_OR 296
#define ARITH_AND 297
#define ARITH_XOR 298
#define INC_OP 299
#define DEC_OP 300
#define BANG 301
#define SEMI 302
#define IF 303
#define ELSE 304
#define FOR 305
#define DO 306
#define WHILE 307
#define CHAR 308
#define SHORT 309
#define INT 310
#define LONG 311
#define UNSIGNED 312
#define SIGNED 313
#define FLOAT 314
#define DOUBLE 315
#define VOID 316
#define STRING 317
#define STATIC 318
#define EXTERN_TOKEN 319
#define STRUCT 320
#define ENUM 321
#define UNION 322
#define CONST 323
#define SIZEOF 324
#define TYPEDEF 325
#define RETURN_TOKEN 326
#define CONTINUE 327
#define BREAK 328
#define GOTO 329
#define PRINT 330
#define COMMA 331
#define DOTDOTDOT 332
#define integer_constant 333
#define character_constant 334
#define string_constant 335
#define floating_constant 336
#define identifier_ref 337
#define type_identifier 338
#define enumeration_constant 339




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 201 "cod/cod.y"
{
    lx_info info;
    sm_ref reference;
    operator_t operator;
    sm_list list;
    char *string;
}
/* Line 1529 of yacc.c.  */
#line 225 "/Users/eisen/prog/ffs/build/cod.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

