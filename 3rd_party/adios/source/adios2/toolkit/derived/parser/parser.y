%skeleton "lalr1.cc"
%require "3.8.2"
%header

%define api.token.raw
%define api.namespace {adios2::detail}
%define api.token.constructor
%define api.value.type variant
%define parse.assert

%code requires {
  #include <tuple>
  #include <vector>
  #include <string>
  namespace adios2
  {
    namespace detail
    {
      class ASTDriver;
    }
  }
}

// The parsing context.
%param { ASTDriver& drv }

%locations

%define parse.trace
%define parse.error detailed
%define parse.lac full

%code {
#include "ASTDriver.h"
#include "ASTNode.h"
}

%define api.token.prefix {TOK_}
%token
  ASSIGN  "="
  COMMA   ","
  COLON   ":"
  L_PAREN "("
  R_PAREN ")"
  L_BRACE "["
  R_BRACE "]"
;

%token <std::string> OPERATOR
%token <std::string> IDENTIFIER "identifier"
%token <std::string> VARNAME
%token <int> INT "number"
%nterm <int> list
%nterm <std::vector<std::tuple<int, int, int>>> indices_list
%nterm <std::tuple<int, int, int>> index
%left OPERATOR

%%
%start lines;
lines:
  assignment lines {}
| exp {}
;

assignment:
  IDENTIFIER ASSIGN VARNAME { drv.add_lookup_entry($1,  $3); }
| IDENTIFIER ASSIGN IDENTIFIER { drv.add_lookup_entry($1,  $3); }
| IDENTIFIER ASSIGN VARNAME L_BRACE indices_list R_BRACE { drv.add_lookup_entry($1, $3, $5); }
| IDENTIFIER ASSIGN IDENTIFIER L_BRACE indices_list R_BRACE { drv.add_lookup_entry($1, $3, $5); }
;

exp:
  "number" {  }
| exp OPERATOR exp { drv.createNode($2, 2); }
| "(" exp ")" {  }
| IDENTIFIER "(" list ")" { drv.createNode($1, $3); }
| IDENTIFIER "[" indices_list "]" { drv.createNode($1, $3); }
| IDENTIFIER  { drv.createNode($1); }
;


indices_list:
  %empty { $$ = {}; }
| indices_list COMMA index { $1.push_back($3); $$ = $1; }
| index { $$ = {$1}; }
;

index:
  %empty                  { $$ = {-1, -1,  1}; }
| INT COLON INT COLON INT { $$ = {$1, $3, $5}; }
| COLON INT COLON INT     { $$ = {-1, $2, $4}; }
| INT COLON COLON INT     { $$ = {$1, -1, $4}; }
| INT COLON INT COLON     { $$ = {$1, $3,  1}; }
| INT COLON INT           { $$ = {$1, $3,  1}; }
| COLON COLON INT         { $$ = {-1, -1, $3}; }
| COLON INT COLON         { $$ = {-1, $2,  1}; }
| COLON INT               { $$ = {-1, $2,  1}; }
| INT COLON COLON         { $$ = {$1, -1,  1}; }
| INT COLON               { $$ = {$1, -1,  1}; }
| INT                     { $$ = {$1, $1,  1}; }
;

list:
  %empty { $$ = 0; }
| exp COMMA list { $$ = $3 + 1; }
| exp { $$ = 1; }
%%

void
adios2::detail::parser::error (const location_type& l, const std::string& m)
{
  std::cerr << l << ": " << m << '\n';
}
