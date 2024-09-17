// A Bison parser, made by GNU Bison 3.8.2.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2021 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
// especially those whose name start with YY_ or yy_.  They are
// private implementation details that can be changed or removed.





#include "parser.h"


// Unqualified %code blocks.
#line 33 "../parser.y"

#include "ASTDriver.h"
#include "ASTNode.h"

#line 51 "parser.cpp"


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif


// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yy_stack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YY_USE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 6 "../parser.y"
namespace adios2 { namespace detail {
#line 144 "parser.cpp"

  /// Build a parser object.
  parser::parser (ASTDriver& drv_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      yy_lac_established_ (false),
      drv (drv_yyarg)
  {}

  parser::~parser ()
  {}

  parser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------.
  | symbol.  |
  `---------*/



  // by_state.
  parser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  parser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  parser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  parser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  parser::symbol_kind_type
  parser::by_state::kind () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return symbol_kind::S_YYEMPTY;
    else
      return YY_CAST (symbol_kind_type, yystos_[+state]);
  }

  parser::stack_symbol_type::stack_symbol_type ()
  {}

  parser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_INT: // "number"
      case symbol_kind::S_list: // list
        value.YY_MOVE_OR_COPY< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_OPERATOR: // OPERATOR
      case symbol_kind::S_IDENTIFIER: // "identifier"
      case symbol_kind::S_VARNAME: // VARNAME
        value.YY_MOVE_OR_COPY< std::string > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_index: // index
        value.YY_MOVE_OR_COPY< std::tuple<int, int, int> > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_indices_list: // indices_list
        value.YY_MOVE_OR_COPY< std::vector<std::tuple<int, int, int>> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  parser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.location))
  {
    switch (that.kind ())
    {
      case symbol_kind::S_INT: // "number"
      case symbol_kind::S_list: // list
        value.move< int > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_OPERATOR: // OPERATOR
      case symbol_kind::S_IDENTIFIER: // "identifier"
      case symbol_kind::S_VARNAME: // VARNAME
        value.move< std::string > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_index: // index
        value.move< std::tuple<int, int, int> > (YY_MOVE (that.value));
        break;

      case symbol_kind::S_indices_list: // indices_list
        value.move< std::vector<std::tuple<int, int, int>> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

    // that is emptied.
    that.kind_ = symbol_kind::S_YYEMPTY;
  }

#if YY_CPLUSPLUS < 201103L
  parser::stack_symbol_type&
  parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_INT: // "number"
      case symbol_kind::S_list: // list
        value.copy< int > (that.value);
        break;

      case symbol_kind::S_OPERATOR: // OPERATOR
      case symbol_kind::S_IDENTIFIER: // "identifier"
      case symbol_kind::S_VARNAME: // VARNAME
        value.copy< std::string > (that.value);
        break;

      case symbol_kind::S_index: // index
        value.copy< std::tuple<int, int, int> > (that.value);
        break;

      case symbol_kind::S_indices_list: // indices_list
        value.copy< std::vector<std::tuple<int, int, int>> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }

  parser::stack_symbol_type&
  parser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    switch (that.kind ())
    {
      case symbol_kind::S_INT: // "number"
      case symbol_kind::S_list: // list
        value.move< int > (that.value);
        break;

      case symbol_kind::S_OPERATOR: // OPERATOR
      case symbol_kind::S_IDENTIFIER: // "identifier"
      case symbol_kind::S_VARNAME: // VARNAME
        value.move< std::string > (that.value);
        break;

      case symbol_kind::S_index: // index
        value.move< std::tuple<int, int, int> > (that.value);
        break;

      case symbol_kind::S_indices_list: // indices_list
        value.move< std::vector<std::tuple<int, int, int>> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  parser::yy_print_ (std::ostream& yyo, const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YY_USE (yyoutput);
    if (yysym.empty ())
      yyo << "empty symbol";
    else
      {
        symbol_kind_type yykind = yysym.kind ();
        yyo << (yykind < YYNTOKENS ? "token" : "nterm")
            << ' ' << yysym.name () << " ("
            << yysym.location << ": ";
        YY_USE (yykind);
        yyo << ')';
      }
  }
#endif

  void
  parser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  parser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  parser::yypop_ (int n) YY_NOEXCEPT
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  parser::debug_level_type
  parser::debug_level () const
  {
    return yydebug_;
  }

  void
  parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  parser::state_type
  parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - YYNTOKENS] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - YYNTOKENS];
  }

  bool
  parser::yy_pact_value_is_default_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  parser::yy_table_value_is_error_ (int yyvalue) YY_NOEXCEPT
  {
    return yyvalue == yytable_ninf_;
  }

  int
  parser::operator() ()
  {
    return parse ();
  }

  int
  parser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // Discard the LAC context in case there still is one left from a
    // previous invocation.
    yy_lac_discard_ ("init");

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';
    YY_STACK_PRINT ();

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token\n";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            symbol_type yylookahead (yylex (drv));
            yyla.move (yylookahead);
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    if (yyla.kind () == symbol_kind::S_YYerror)
    {
      // The scanner already issued an error message, process directly
      // to error recovery.  But do not keep the error token as
      // lookahead, it is too special and may lead us to an endless
      // loop in error recovery. */
      yyla.kind_ = symbol_kind::S_YYUNDEF;
      goto yyerrlab1;
    }

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.kind ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.kind ())
      {
        if (!yy_lac_establish_ (yyla.kind ()))
          goto yyerrlab;
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        if (!yy_lac_establish_ (yyla.kind ()))
          goto yyerrlab;

        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    yy_lac_discard_ ("shift");
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
      switch (yyr1_[yyn])
    {
      case symbol_kind::S_INT: // "number"
      case symbol_kind::S_list: // list
        yylhs.value.emplace< int > ();
        break;

      case symbol_kind::S_OPERATOR: // OPERATOR
      case symbol_kind::S_IDENTIFIER: // "identifier"
      case symbol_kind::S_VARNAME: // VARNAME
        yylhs.value.emplace< std::string > ();
        break;

      case symbol_kind::S_index: // index
        yylhs.value.emplace< std::tuple<int, int, int> > ();
        break;

      case symbol_kind::S_indices_list: // indices_list
        yylhs.value.emplace< std::vector<std::tuple<int, int, int>> > ();
        break;

      default:
        break;
    }


      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 2: // lines: assignment lines
#line 61 "../parser.y"
                   {}
#line 644 "parser.cpp"
    break;

  case 3: // lines: exp
#line 62 "../parser.y"
      {}
#line 650 "parser.cpp"
    break;

  case 4: // assignment: "identifier" "=" VARNAME
#line 66 "../parser.y"
                            { drv.add_lookup_entry(yystack_[2].value.as < std::string > (),  yystack_[0].value.as < std::string > ()); }
#line 656 "parser.cpp"
    break;

  case 5: // assignment: "identifier" "=" "identifier"
#line 67 "../parser.y"
                               { drv.add_lookup_entry(yystack_[2].value.as < std::string > (),  yystack_[0].value.as < std::string > ()); }
#line 662 "parser.cpp"
    break;

  case 6: // assignment: "identifier" "=" VARNAME "[" indices_list "]"
#line 68 "../parser.y"
                                                         { drv.add_lookup_entry(yystack_[5].value.as < std::string > (), yystack_[3].value.as < std::string > (), yystack_[1].value.as < std::vector<std::tuple<int, int, int>> > ()); }
#line 668 "parser.cpp"
    break;

  case 7: // assignment: "identifier" "=" "identifier" "[" indices_list "]"
#line 69 "../parser.y"
                                                            { drv.add_lookup_entry(yystack_[5].value.as < std::string > (), yystack_[3].value.as < std::string > (), yystack_[1].value.as < std::vector<std::tuple<int, int, int>> > ()); }
#line 674 "parser.cpp"
    break;

  case 8: // exp: "number"
#line 73 "../parser.y"
           {  }
#line 680 "parser.cpp"
    break;

  case 9: // exp: exp OPERATOR exp
#line 74 "../parser.y"
                   { drv.createNode(yystack_[1].value.as < std::string > (), 2); }
#line 686 "parser.cpp"
    break;

  case 10: // exp: "(" exp ")"
#line 75 "../parser.y"
              {  }
#line 692 "parser.cpp"
    break;

  case 11: // exp: "identifier" "(" list ")"
#line 76 "../parser.y"
                          { drv.createNode(yystack_[3].value.as < std::string > (), yystack_[1].value.as < int > ()); }
#line 698 "parser.cpp"
    break;

  case 12: // exp: "identifier" "[" indices_list "]"
#line 77 "../parser.y"
                                  { drv.createNode(yystack_[3].value.as < std::string > (), yystack_[1].value.as < std::vector<std::tuple<int, int, int>> > ()); }
#line 704 "parser.cpp"
    break;

  case 13: // exp: "identifier"
#line 78 "../parser.y"
              { drv.createNode(yystack_[0].value.as < std::string > ()); }
#line 710 "parser.cpp"
    break;

  case 14: // indices_list: %empty
#line 83 "../parser.y"
         { yylhs.value.as < std::vector<std::tuple<int, int, int>> > () = {}; }
#line 716 "parser.cpp"
    break;

  case 15: // indices_list: indices_list "," index
#line 84 "../parser.y"
                           { yystack_[2].value.as < std::vector<std::tuple<int, int, int>> > ().push_back(yystack_[0].value.as < std::tuple<int, int, int> > ()); yylhs.value.as < std::vector<std::tuple<int, int, int>> > () = yystack_[2].value.as < std::vector<std::tuple<int, int, int>> > (); }
#line 722 "parser.cpp"
    break;

  case 16: // indices_list: index
#line 85 "../parser.y"
        { yylhs.value.as < std::vector<std::tuple<int, int, int>> > () = {yystack_[0].value.as < std::tuple<int, int, int> > ()}; }
#line 728 "parser.cpp"
    break;

  case 17: // index: %empty
#line 89 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {-1, -1,  1}; }
#line 734 "parser.cpp"
    break;

  case 18: // index: "number" ":" "number" ":" "number"
#line 90 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[4].value.as < int > (), yystack_[2].value.as < int > (), yystack_[0].value.as < int > ()}; }
#line 740 "parser.cpp"
    break;

  case 19: // index: ":" "number" ":" "number"
#line 91 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {-1, yystack_[2].value.as < int > (), yystack_[0].value.as < int > ()}; }
#line 746 "parser.cpp"
    break;

  case 20: // index: "number" ":" ":" "number"
#line 92 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[3].value.as < int > (), -1, yystack_[0].value.as < int > ()}; }
#line 752 "parser.cpp"
    break;

  case 21: // index: "number" ":" "number" ":"
#line 93 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[3].value.as < int > (), yystack_[1].value.as < int > (),  1}; }
#line 758 "parser.cpp"
    break;

  case 22: // index: "number" ":" "number"
#line 94 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[2].value.as < int > (), yystack_[0].value.as < int > (),  1}; }
#line 764 "parser.cpp"
    break;

  case 23: // index: ":" ":" "number"
#line 95 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {-1, -1, yystack_[0].value.as < int > ()}; }
#line 770 "parser.cpp"
    break;

  case 24: // index: ":" "number" ":"
#line 96 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {-1, yystack_[1].value.as < int > (),  1}; }
#line 776 "parser.cpp"
    break;

  case 25: // index: ":" "number"
#line 97 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {-1, yystack_[0].value.as < int > (),  1}; }
#line 782 "parser.cpp"
    break;

  case 26: // index: "number" ":" ":"
#line 98 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[2].value.as < int > (), -1,  1}; }
#line 788 "parser.cpp"
    break;

  case 27: // index: "number" ":"
#line 99 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[1].value.as < int > (), -1,  1}; }
#line 794 "parser.cpp"
    break;

  case 28: // index: "number"
#line 100 "../parser.y"
                          { yylhs.value.as < std::tuple<int, int, int> > () = {yystack_[0].value.as < int > (), yystack_[0].value.as < int > (),  1}; }
#line 800 "parser.cpp"
    break;

  case 29: // list: %empty
#line 104 "../parser.y"
         { yylhs.value.as < int > () = 0; }
#line 806 "parser.cpp"
    break;

  case 30: // list: exp "," list
#line 105 "../parser.y"
                 { yylhs.value.as < int > () = yystack_[0].value.as < int > () + 1; }
#line 812 "parser.cpp"
    break;

  case 31: // list: exp
#line 106 "../parser.y"
      { yylhs.value.as < int > () = 1; }
#line 818 "parser.cpp"
    break;


#line 822 "parser.cpp"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        context yyctx (*this, yyla);
        std::string msg = yysyntax_error_ (yyctx);
        error (yyla.location, YY_MOVE (msg));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.kind () == symbol_kind::S_YYEOF)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    YY_STACK_PRINT ();
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    // Pop stack until we find a state that shifts the error token.
    for (;;)
      {
        yyn = yypact_[+yystack_[0].state];
        if (!yy_pact_value_is_default_ (yyn))
          {
            yyn += symbol_kind::S_YYerror;
            if (0 <= yyn && yyn <= yylast_
                && yycheck_[yyn] == symbol_kind::S_YYerror)
              {
                yyn = yytable_[yyn];
                if (0 < yyn)
                  break;
              }
          }

        // Pop the current state because it cannot handle the error token.
        if (yystack_.size () == 1)
          YYABORT;

        yyerror_range[1].location = yystack_[0].location;
        yy_destroy_ ("Error: popping", yystack_[0]);
        yypop_ ();
        YY_STACK_PRINT ();
      }
    {
      stack_symbol_type error_token;

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      yy_lac_discard_ ("error recovery");
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
    goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
  `-------------------------------------*/
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;


  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
  `-----------------------------------*/
  yyabortlab:
    yyresult = 1;
    goto yyreturn;


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    YY_STACK_PRINT ();
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

  const char *
  parser::symbol_name (symbol_kind_type yysymbol)
  {
    static const char *const yy_sname[] =
    {
    "end of file", "error", "invalid token", "=", ",", ":", "(", ")", "[",
  "]", "OPERATOR", "identifier", "VARNAME", "number", "$accept", "lines",
  "assignment", "exp", "indices_list", "index", "list", YY_NULLPTR
    };
    return yy_sname[yysymbol];
  }



  // parser::context.
  parser::context::context (const parser& yyparser, const symbol_type& yyla)
    : yyparser_ (yyparser)
    , yyla_ (yyla)
  {}

  int
  parser::context::expected_tokens (symbol_kind_type yyarg[], int yyargn) const
  {
    // Actual number of expected tokens
    int yycount = 0;

#if YYDEBUG
    // Execute LAC once. We don't care if it is successful, we
    // only do it for the sake of debugging output.
    if (!yyparser_.yy_lac_established_)
      yyparser_.yy_lac_check_ (yyla_.kind ());
#endif

    for (int yyx = 0; yyx < YYNTOKENS; ++yyx)
      {
        symbol_kind_type yysym = YY_CAST (symbol_kind_type, yyx);
        if (yysym != symbol_kind::S_YYerror
            && yysym != symbol_kind::S_YYUNDEF
            && yyparser_.yy_lac_check_ (yysym))
          {
            if (!yyarg)
              ++yycount;
            else if (yycount == yyargn)
              return 0;
            else
              yyarg[yycount++] = yysym;
          }
      }
    if (yyarg && yycount == 0 && 0 < yyargn)
      yyarg[0] = symbol_kind::S_YYEMPTY;
    return yycount;
  }




  bool
  parser::yy_lac_check_ (symbol_kind_type yytoken) const
  {
    // Logically, the yylac_stack's lifetime is confined to this function.
    // Clear it, to get rid of potential left-overs from previous call.
    yylac_stack_.clear ();
    // Reduce until we encounter a shift and thereby accept the token.
#if YYDEBUG
    YYCDEBUG << "LAC: checking lookahead " << symbol_name (yytoken) << ':';
#endif
    std::ptrdiff_t lac_top = 0;
    while (true)
      {
        state_type top_state = (yylac_stack_.empty ()
                                ? yystack_[lac_top].state
                                : yylac_stack_.back ());
        int yyrule = yypact_[+top_state];
        if (yy_pact_value_is_default_ (yyrule)
            || (yyrule += yytoken) < 0 || yylast_ < yyrule
            || yycheck_[yyrule] != yytoken)
          {
            // Use the default action.
            yyrule = yydefact_[+top_state];
            if (yyrule == 0)
              {
                YYCDEBUG << " Err\n";
                return false;
              }
          }
        else
          {
            // Use the action from yytable.
            yyrule = yytable_[yyrule];
            if (yy_table_value_is_error_ (yyrule))
              {
                YYCDEBUG << " Err\n";
                return false;
              }
            if (0 < yyrule)
              {
                YYCDEBUG << " S" << yyrule << '\n';
                return true;
              }
            yyrule = -yyrule;
          }
        // By now we know we have to simulate a reduce.
        YYCDEBUG << " R" << yyrule - 1;
        // Pop the corresponding number of values from the stack.
        {
          std::ptrdiff_t yylen = yyr2_[yyrule];
          // First pop from the LAC stack as many tokens as possible.
          std::ptrdiff_t lac_size = std::ptrdiff_t (yylac_stack_.size ());
          if (yylen < lac_size)
            {
              yylac_stack_.resize (std::size_t (lac_size - yylen));
              yylen = 0;
            }
          else if (lac_size)
            {
              yylac_stack_.clear ();
              yylen -= lac_size;
            }
          // Only afterwards look at the main stack.
          // We simulate popping elements by incrementing lac_top.
          lac_top += yylen;
        }
        // Keep top_state in sync with the updated stack.
        top_state = (yylac_stack_.empty ()
                     ? yystack_[lac_top].state
                     : yylac_stack_.back ());
        // Push the resulting state of the reduction.
        state_type state = yy_lr_goto_state_ (top_state, yyr1_[yyrule]);
        YYCDEBUG << " G" << int (state);
        yylac_stack_.push_back (state);
      }
  }

  // Establish the initial context if no initial context currently exists.
  bool
  parser::yy_lac_establish_ (symbol_kind_type yytoken)
  {
    /* Establish the initial context for the current lookahead if no initial
       context is currently established.

       We define a context as a snapshot of the parser stacks.  We define
       the initial context for a lookahead as the context in which the
       parser initially examines that lookahead in order to select a
       syntactic action.  Thus, if the lookahead eventually proves
       syntactically unacceptable (possibly in a later context reached via a
       series of reductions), the initial context can be used to determine
       the exact set of tokens that would be syntactically acceptable in the
       lookahead's place.  Moreover, it is the context after which any
       further semantic actions would be erroneous because they would be
       determined by a syntactically unacceptable token.

       yy_lac_establish_ should be invoked when a reduction is about to be
       performed in an inconsistent state (which, for the purposes of LAC,
       includes consistent states that don't know they're consistent because
       their default reductions have been disabled).

       For parse.lac=full, the implementation of yy_lac_establish_ is as
       follows.  If no initial context is currently established for the
       current lookahead, then check if that lookahead can eventually be
       shifted if syntactic actions continue from the current context.  */
    if (yy_lac_established_)
      return true;
    else
      {
#if YYDEBUG
        YYCDEBUG << "LAC: initial context established for "
                 << symbol_name (yytoken) << '\n';
#endif
        yy_lac_established_ = true;
        return yy_lac_check_ (yytoken);
      }
  }

  // Discard any previous initial lookahead context.
  void
  parser::yy_lac_discard_ (const char* event)
  {
   /* Discard any previous initial lookahead context because of Event,
      which may be a lookahead change or an invalidation of the currently
      established initial context for the current lookahead.

      The most common example of a lookahead change is a shift.  An example
      of both cases is syntax error recovery.  That is, a syntax error
      occurs when the lookahead is syntactically erroneous for the
      currently established initial context, so error recovery manipulates
      the parser stacks to try to find a new initial context in which the
      current lookahead is syntactically acceptable.  If it fails to find
      such a context, it discards the lookahead.  */
    if (yy_lac_established_)
      {
        YYCDEBUG << "LAC: initial context discarded due to "
                 << event << '\n';
        yy_lac_established_ = false;
      }
  }


  int
  parser::yy_syntax_error_arguments_ (const context& yyctx,
                                                 symbol_kind_type yyarg[], int yyargn) const
  {
    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
         In the first two cases, it might appear that the current syntax
         error should have been detected in the previous state when
         yy_lac_check was invoked.  However, at that time, there might
         have been a different syntax error that discarded a different
         initial context during error recovery, leaving behind the
         current lookahead.
    */

    if (!yyctx.lookahead ().empty ())
      {
        if (yyarg)
          yyarg[0] = yyctx.token ();
        int yyn = yyctx.expected_tokens (yyarg ? yyarg + 1 : yyarg, yyargn - 1);
        return yyn + 1;
      }
    return 0;
  }

  // Generate an error message.
  std::string
  parser::yysyntax_error_ (const context& yyctx) const
  {
    // Its maximum.
    enum { YYARGS_MAX = 5 };
    // Arguments of yyformat.
    symbol_kind_type yyarg[YYARGS_MAX];
    int yycount = yy_syntax_error_arguments_ (yyctx, yyarg, YYARGS_MAX);

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
      default: // Avoid compiler warnings.
        YYCASE_ (0, YY_("syntax error"));
        YYCASE_ (1, YY_("syntax error, unexpected %s"));
        YYCASE_ (2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_ (3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_ (4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_ (5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    std::ptrdiff_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += symbol_name (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char parser::yypact_ninf_ = -7;

  const signed char parser::yytable_ninf_ = -1;

  const signed char
  parser::yypact_[] =
  {
       6,     9,    22,    -7,     5,     6,    13,    27,    -6,    -4,
       9,    -3,    -7,    -7,     9,    -7,    30,    31,    14,    24,
      -2,    35,    12,    -7,    -7,    -3,    -3,     9,    -7,    28,
      37,     1,    -3,    -7,    23,    25,    -7,    -7,    32,    33,
      38,    -7,    -7,    -7,    -7,    -7,    34,    -7
  };

  const signed char
  parser::yydefact_[] =
  {
       0,     0,    13,     8,     0,     0,     3,    13,     0,     0,
      29,    14,     1,     2,     0,    10,     5,     4,    31,     0,
       0,    28,     0,    16,     9,    14,    14,    29,    11,     0,
      25,    27,    17,    12,     0,     0,    30,    23,    24,    26,
      22,    15,     7,     6,    19,    20,    21,    18
  };

  const signed char
  parser::yypgoto_[] =
  {
      -7,    39,    -7,    -1,    11,    16,    26
  };

  const signed char
  parser::yydefgoto_[] =
  {
       0,     4,     5,     6,    22,    23,    19
  };

  const signed char
  parser::yytable_[] =
  {
       8,    15,    20,    29,    14,    12,    39,    16,    17,    18,
      21,    30,     1,    24,    40,     1,    32,     2,    27,     3,
       7,    33,     3,    14,    14,     9,    18,    32,    10,    32,
      11,    28,    42,    10,    43,    11,    34,    35,    25,    26,
      31,    37,    38,    46,    13,    44,    45,    47,    41,     0,
       0,     0,     0,    36
  };

  const signed char
  parser::yycheck_[] =
  {
       1,     7,     5,     5,    10,     0,     5,    11,    12,    10,
      13,    13,     6,    14,    13,     6,     4,    11,     4,    13,
      11,     9,    13,    10,    10,     3,    27,     4,     6,     4,
       8,     7,     9,     6,     9,     8,    25,    26,     8,     8,
       5,    13,     5,     5,     5,    13,    13,    13,    32,    -1,
      -1,    -1,    -1,    27
  };

  const signed char
  parser::yystos_[] =
  {
       0,     6,    11,    13,    15,    16,    17,    11,    17,     3,
       6,     8,     0,    15,    10,     7,    11,    12,    17,    20,
       5,    13,    18,    19,    17,     8,     8,     4,     7,     5,
      13,     5,     4,     9,    18,    18,    20,    13,     5,     5,
      13,    19,     9,     9,    13,    13,     5,    13
  };

  const signed char
  parser::yyr1_[] =
  {
       0,    14,    15,    15,    16,    16,    16,    16,    17,    17,
      17,    17,    17,    17,    18,    18,    18,    19,    19,    19,
      19,    19,    19,    19,    19,    19,    19,    19,    19,    20,
      20,    20
  };

  const signed char
  parser::yyr2_[] =
  {
       0,     2,     2,     1,     3,     3,     6,     6,     1,     3,
       3,     4,     4,     1,     0,     3,     1,     0,     5,     4,
       4,     4,     3,     3,     3,     2,     3,     2,     1,     0,
       3,     1
  };




#if YYDEBUG
  const signed char
  parser::yyrline_[] =
  {
       0,    61,    61,    62,    66,    67,    68,    69,    73,    74,
      75,    76,    77,    78,    83,    84,    85,    89,    90,    91,
      92,    93,    94,    95,    96,    97,    98,    99,   100,   104,
     105,   106
  };

  void
  parser::yy_stack_print_ () const
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  void
  parser::yy_reduce_print_ (int yyrule) const
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG


#line 6 "../parser.y"
} } // adios2::detail
#line 1406 "parser.cpp"

#line 107 "../parser.y"


void
adios2::detail::parser::error (const location_type& l, const std::string& m)
{
  std::cerr << l << ": " << m << '\n';
}
