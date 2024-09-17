%{
#include <cerrno>
#include <cstdlib>
#include <limits.h>
#include <cstring> // strerror
#include <string>
#include "ASTDriver.h"
#include "parser.h"
#if defined(_MSC_VER)
#include <io.h>
#define YY_NO_UNISTD_H
#define strdup _strdup
#define isatty _isatty
#define fileno _fileno
#include <sys/types.h>
#include <sys/stat.h>
#endif
%}

%{
#if defined __clang__
# define CLANG_VERSION (__clang_major__ * 100 + __clang_minor__)
#endif

// Clang and ICC like to pretend they are GCC.
#if defined __GNUC__ && !defined __clang__ && !defined __ICC
# define GCC_VERSION (__GNUC__ * 100 + __GNUC_MINOR__)
#endif

// Pacify warnings in yy_init_buffer (observed with Flex 2.6.4)
// and GCC 6.4.0, 7.3.0 with -O3.
#if defined GCC_VERSION && 600 <= GCC_VERSION
# pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

// This example uses Flex's C back end, yet compiles it as C++.
// So expect warnings about C style casts and NULL.
#if defined CLANG_VERSION && 500 <= CLANG_VERSION
# pragma clang diagnostic ignored "-Wold-style-cast"
# pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined GCC_VERSION && 407 <= GCC_VERSION
# pragma GCC diagnostic ignored "-Wold-style-cast"
# pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#define FLEX_VERSION (YY_FLEX_MAJOR_VERSION * 100 + YY_FLEX_MINOR_VERSION)

// Old versions of Flex (2.5.35) generate an incomplete documentation comment.
//
//  In file included from src/scan-code-c.c:3:
//  src/scan-code.c:2198:21: error: empty paragraph passed to '@param' command
//        [-Werror,-Wdocumentation]
//   * @param line_number
//     ~~~~~~~~~~~~~~~~~^
//  1 error generated.
#if FLEX_VERSION < 206 && defined CLANG_VERSION
# pragma clang diagnostic ignored "-Wdocumentation"
#endif

// Old versions of Flex (2.5.35) use 'register'.  Warnings introduced in
// GCC 7 and Clang 6.
#if FLEX_VERSION < 206
# if defined CLANG_VERSION && 600 <= CLANG_VERSION
#  pragma clang diagnostic ignored "-Wdeprecated-register"
# elif defined GCC_VERSION && 700 <= GCC_VERSION
#  pragma GCC diagnostic ignored "-Wregister"
# endif
#endif

#if FLEX_VERSION < 206
# if defined CLANG_VERSION
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wdocumentation"
#  pragma clang diagnostic ignored "-Wshorten-64-to-32"
#  pragma clang diagnostic ignored "-Wsign-conversion"
# elif defined GCC_VERSION
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
# endif
#endif

// Flex 2.6.4, GCC 9
// warning: useless cast to type 'int' [-Wuseless-cast]
// 1361 |   YY_CURRENT_BUFFER_LVALUE->yy_buf_size = (int) (new_size - 2);
//      |                                                 ^
#if defined GCC_VERSION && 900 <= GCC_VERSION
# pragma GCC diagnostic ignored "-Wuseless-cast"
#endif
%}

%option noyywrap nounput noinput batch

%{
  // A number symbol corresponding to the value in S.
  adios2::detail::parser::symbol_type
  make_INT (const std::string &s, const adios2::detail::parser::location_type& loc);
%}

op    [-+*/^]
id    [a-zA-Z][a-zA-Z_0-9_]*
path  [a-zA-Z][a-zA-Z0-9_\\./-]*
int   [0-9]+
blank [ \t\r]

%{
  // Code run each time a pattern is matched.
  # define YY_USER_ACTION  loc.columns (yyleng);
%}
%%
%{
  // A handy shortcut to the location held by the adios2::detail::ASTDriver.
  adios2::detail::location& loc = drv.location;
  // Code run each time yylex is called.
  loc.step ();
%}
{blank}+   loc.step ();
\n+        loc.lines (yyleng); loc.step ();

"="        return adios2::detail::parser::make_ASSIGN  (loc);
","        return adios2::detail::parser::make_COMMA   (loc);
":"        return adios2::detail::parser::make_COLON   (loc);
"("        return adios2::detail::parser::make_L_PAREN (loc);
")"        return adios2::detail::parser::make_R_PAREN (loc);
"["        return adios2::detail::parser::make_L_BRACE (loc);
"]"        return adios2::detail::parser::make_R_BRACE (loc);

{int}      return make_INT (yytext, loc);
{op}       return adios2::detail::parser::make_OPERATOR (yytext, loc);
{id}       return adios2::detail::parser::make_IDENTIFIER (yytext, loc);
{path}     return adios2::detail::parser::make_VARNAME (yytext, loc);
.          {
             throw adios2::detail::parser::syntax_error
               (loc, "invalid character: " + std::string(yytext));
}
<<EOF>>    return adios2::detail::parser::make_YYEOF (loc);
%%

adios2::detail::parser::symbol_type
make_INT (const std::string &s, const adios2::detail::parser::location_type& loc)
{
  errno = 0;
  long n = strtol (s.c_str(), NULL, 10);
  if (! (INT_MIN <= n && n <= INT_MAX && errno != ERANGE))
    throw adios2::detail::parser::syntax_error (loc, "integer is out of range: " + s);
  return adios2::detail::parser::make_INT ((int) n, loc);
}

void
adios2::detail::ASTDriver::parse (const std::string input)
{
  adios2::detail::parser parse (*this);
  yy_scan_string(input.c_str());
  parse.set_debug_level (trace_parsing);
  parse ();
}
