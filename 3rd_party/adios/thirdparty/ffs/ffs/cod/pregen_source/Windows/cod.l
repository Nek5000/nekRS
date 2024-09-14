O   [0-7]
D   [0-9]
NZ  [1-9]
L   [a-zA-Z_]
A   [a-zA-Z_0-9]
H   [a-fA-F0-9]
HP  (0[xX])
E   ([Ee][+-]?{D}+)
P   ([Pp][+-]?{D}+)
FS  (f|F|l|L)
IS  (((u|U)(l|L|ll|LL)?)|((l|L|ll|LL)(u|U)?))
CP  (u|U|L)
SP  (u8|u|U|L)
ES  (\\(['"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
WS  [ \t\v\n\f]

%{
static	int lex_offset = 1;
static  int line_count = 1;
static char *create_string_from_yytext();
#define RETURN(val) {yylval.info.lx_srcpos.line = line_count; yylval.info.lx_srcpos.character = lex_offset; lex_offset +=yyleng; return(val);}
static int type_count;
static char**types;
static char**enums;
#ifdef input
/* lex, not flex */
static char *saved_input_str = NULL;
static void 
terminate_string_parse()
{
    if (saved_input_str) {
	free(saved_input_str);
        saved_input_str = NULL;
    }
}

static char *input_str = NULL;
static void
setup_for_string_parse(string, defined_types, enum_constants)
const char *string;
char **defined_types;
char **enum_constants;
{
    int len = strlen(string);

    types = defined_types;
    enums = enum_constants;
    if (saved_input_str != NULL) free(input_str);
    input_str = malloc(len + 2);
    saved_input_str = input_str;
    strcpy(input_str, string);
    input_str[len] = 0;
    input_str[len + 1] = 0;
    lex_offset = 1;
    line_count = 1;
}
#undef input
#define input() (((yytchar=*input_str)==0)?0:(input_str++,yytchar==10?(yylineno++,yytchar):yytchar))
#undef unput
#define unput(c) (input_str--,c==10?yylineno--:yylineno)
#else
/* flex */

#ifdef YY_INPUT
#undef YY_INPUT
extern int my_yy_input();
#define YY_INPUT(buf,x,y) x=my_yy_input(buf,x,y)
#endif

/* end of lex ifdef */
#endif
static
int
is_defined_type(char *id)
{
    int i = 0;
    while(types && types[i]) {
	if (strcmp(id, types[i]) == 0) return 1;
	i++;
    }
    return 0;
}
static
int
is_enumeration_constant(char *id)
{
    int i = 0;
    while(enums && enums[i]) {
	if (strcmp(id, enums[i]) == 0) return 1;
	i++;
    }
    return 0;
}
static void check_strbuf();
static int buffer_len;
static char *string_buffer;
static char *string_buf_ptr;
%}
%x string_cond comment

%%


"->"			{RETURN(ARROW);}
"("			{RETURN(LPAREN);}
")"			{RETURN(RPAREN);}
"["			{RETURN(LBRACKET);}
"]"			{RETURN(RBRACKET);}
"..."			{RETURN(DOTDOTDOT);}
"."			{RETURN(DOT);}
","			{RETURN(COMMA);}
"\*"			{RETURN(STAR);}
"\@"			{RETURN(AT);}
"/"			{RETURN(SLASH);}
"%"			{RETURN(MODULUS);}
"+"			{RETURN(PLUS);}
"-"			{RETURN(MINUS);}
"~"			{RETURN(TILDE);}
"<="			{RETURN(LEQ);}
"<"			{RETURN(LT);}
">="			{RETURN(GEQ);}
">"			{RETURN(GT);}
"<<"			{RETURN(LEFT_SHIFT);}
">>"			{RETURN(RIGHT_SHIFT);}
"=="			{RETURN(EQ);}
"!="			{RETURN(NEQ);}
"="			{RETURN(ASSIGN);}
"*="			{RETURN(MUL_ASSIGN);}
"/="			{RETURN(DIV_ASSIGN);}
"%="			{RETURN(MOD_ASSIGN);}
"+="			{RETURN(ADD_ASSIGN);}
"-="			{RETURN(SUB_ASSIGN);}
"<<="			{RETURN(LEFT_ASSIGN);}
">>="			{RETURN(RIGHT_ASSIGN);}
"&="			{RETURN(AND_ASSIGN);}
"^="			{RETURN(XOR_ASSIGN);}
"|="			{RETURN(OR_ASSIGN);}
"||"			{RETURN(LOG_OR);}
"&&"			{RETURN(LOG_AND);}
"|"			{RETURN(ARITH_OR);}
"&"			{RETURN(ARITH_AND);}
"^"			{RETURN(ARITH_XOR);}
"++"			{RETURN(INC_OP);}
"--"			{RETURN(DEC_OP);}
";"			{RETURN(SEMI);}
"if"			{RETURN(IF);}
"else"			{RETURN(ELSE);}
"for"			{RETURN(FOR);}
"while"			{RETURN(WHILE);}
"do"			{RETURN(DO);}
"unsigned"		{RETURN(UNSIGNED);}
"signed"		{RETURN(SIGNED);}
"short"			{RETURN(SHORT);}
"int"			{RETURN(INT);}
"long"			{RETURN(LONG);}
"char"			{RETURN(CHAR);}
"string"		{RETURN(STRING);}
"float"			{RETURN(FLOAT);}
"double"		{RETURN(DOUBLE);}
"void"			{RETURN(VOID);}
static			{RETURN(STATIC);}
extern			{RETURN(EXTERN_TOKEN);}
typedef			{RETURN(TYPEDEF);}
continue		{RETURN(CONTINUE);}
break			{RETURN(BREAK);}
goto			{RETURN(GOTO);}
const			{RETURN(CONST);}
sizeof			{RETURN(SIZEOF);}
struct			{RETURN(STRUCT);}
enum			{RETURN(ENUM);}
union			{RETURN(UNION);}
return			{RETURN(RETURN_TOKEN);}
"{"			{
			    int count = 0;
			    while(types && types[count]) count++;
			    yylval.info.type_stack_count = count;
			    RETURN(LCURLY);
			}
"}"			{RETURN(RCURLY);}
"!"			{RETURN(BANG);}
":"			{RETURN(COLON);}
"?"			{RETURN(QUESTION);}
[A-Za-z_][A-Za-z0-9_]*	{
				yylval.info.string = create_string_from_yytext();
				if (is_defined_type(yylval.info.string)) {
				    RETURN(type_identifier);
				} else if (is_enumeration_constant(yylval.info.string)) {
				    RETURN(enumeration_constant);
				} else {
				    RETURN(identifier_ref);
				}
			}

\"      {
	  buffer_len = 20;
	  string_buffer = malloc(20 + 1);
	  string_buf_ptr = string_buffer; BEGIN(string_cond);
}

<string_cond>\"        { /* saw closing quote - all done */
        BEGIN(INITIAL);
        *string_buf_ptr = '\0';
        /* return string constant token type and
         * value to parser
         */
	yylval.info.string = string_buffer;
	RETURN(string_constant);
}

<string_cond>\n        {
	yyerror("Unterminated string constant");
}

<string_cond>\\x[0-9A-Fa-f]+ {
        /* hex escape sequence */
        int result;

        (void) sscanf( yytext + 2, "%x", &result );

        if ( result > 0xff ) {
	   yyerror("bad hex escape character");
	}

        check_strbuf();
	*string_buf_ptr++ = result;
}

<string_cond>\\[0-7]{1,3} {
        /* octal escape sequence */
        int result;

        (void) sscanf( yytext + 1, "%o", &result );

        if ( result > 0xff ) {
	   yyerror("bad octal escape character");
	}

        check_strbuf();
	*string_buf_ptr++ = result;
}

<string_cond>\\[0-9]+ {
	 yyerror("bad character escape");
}

<string_cond>\\n  		{check_strbuf();*string_buf_ptr++ = '\n';}
<string_cond>\\t  		{check_strbuf();*string_buf_ptr++ = '\t';}
<string_cond>\\r  		{check_strbuf();*string_buf_ptr++ = '\r';}
<string_cond>\\b  		{check_strbuf();*string_buf_ptr++ = '\b';}
<string_cond>\\f  		{check_strbuf();*string_buf_ptr++ = '\f';}

<string_cond>\\(.|\n)		{check_strbuf();*string_buf_ptr++ = yytext[1];}

<string_cond>[^\\\n\"]+        {
        char *yptr = yytext;

        while ( *yptr ) {
		check_strbuf();
                *string_buf_ptr++ = *yptr++;
        }
}



"/*"         BEGIN(comment);

<comment>[^*\n]*        {lex_offset += yyleng;} /* eat anything that's not a '*' */
<comment>"*"+[^*/\n]*   {lex_offset += yyleng;} /* eat up '*'s not followed by '/'s */
<comment>\n             {++line_count;lex_offset = 1;}
<comment>"*"+"/"        {lex_offset += yyleng;BEGIN(INITIAL);}
"//".*                                    { /* consume //-comment */ }

{CP}?"'"([^'\\\n]|{ES})+"'"		{
				yylval.info.string = create_string_from_yytext();
				RETURN(character_constant);
			}
[0-9][0-9]*{IS}?	{
				yylval.info.string = create_string_from_yytext();
				RETURN(integer_constant);
			}
0x[0-9a-fA-F]+{IS}?	{
				yylval.info.string = create_string_from_yytext();
				RETURN(integer_constant);
			}

0b[01]+{IS}?			{
				yylval.info.string = create_string_from_yytext();
				RETURN(integer_constant);
			}

{HP}{H}+{P}{FS}?  	{
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}
{HP}{H}*"."{H}+{P}{FS}? {
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}
{HP}{H}+"."{P}{FS}? 	{
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}
{D}+{E}{FS}?  		{
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}

{D}*"."{D}+({E})?{FS}?  {
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}

{D}+"."{D}*({E})?{FS}?  {
				yylval.info.string = create_string_from_yytext();
				RETURN(floating_constant);
			}

[ \t]			{lex_offset += yyleng;}
[\n]			{lex_offset = 1; line_count++;}
%%
static char *create_string_from_yytext()
{
  char *st = (char *) malloc((yyleng+1)*sizeof(char));
  strcpy (st, yytext);
  return(st);
}

static void check_strbuf()
{
    intptr_t cur_len = string_buf_ptr - string_buffer;
    if ((cur_len + 1) == buffer_len) {
	buffer_len += 20;
	string_buffer = realloc(string_buffer, buffer_len + 1);
	string_buf_ptr = string_buffer + cur_len;
    }
}

#define yy_size_t int
#ifndef yyconst
#define yyconst
#endif 

#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
#define YY_USE_PROTOS
#endif

#ifndef YY_PROTO
#ifdef YY_USE_PROTOS
#define YY_PROTO(proto) proto
#else
#define YY_PROTO(proto) ()
#endif
#endif

extern int
yywrap YY_PROTO(( void ))
{
    return 1;
}


/* flex, not lex */
void yy_delete_buffer YY_PROTO((YY_BUFFER_STATE b));


static YY_BUFFER_STATE bb = NULL;

static void
reset_types_table(char **defined_types, char **enumerated_constants)
{
    types = defined_types;
    enums = enumerated_constants;
}

static void
setup_for_string_parse(const char *string, char **defined_types, char **enum_constants)
{
    types = defined_types;
    enums = enum_constants;

    if ((bb = yy_scan_string(string)) == NULL) {
       fprintf(stderr, "yyscan_buffer_failed\n");
     }
    lex_offset = 1;
    line_count = 1;
}

static void 
terminate_string_parse()
{
    if (bb) {
		yy_delete_buffer(bb);
		bb = NULL;
	}
}

