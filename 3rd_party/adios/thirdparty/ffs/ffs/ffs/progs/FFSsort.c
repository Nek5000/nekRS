#include "config.h"
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "ffs.h"
#include "cod.h"

char *usage = "\
Usage: FFSsort [<option>] <in filename> <out filename>\n\
       FFSsort sortsFFS files. \n\
       Options are:\n\
		-e <expression>   --  Expression to sort on\n\
		-v		  --  Turn on verbose output\n\
		-q		  --  Turn off warnings\n\
\n\
\n\
	The <expression> value can have a \"<format_name>:\" prefix, in which case it \n\
	is applied only to record formats with a matching name.  <Expression> will\n\
	be used to create a CoD subroutine which extracts a sort value from the \n\
	incoming record.  In particular, the incoming record will be provided as a \n\
	parameter named \"input\" to the extraction routine, and the <expression>\n\
	value must create a legal CoD subroutine when inserted into one of the following \n\
	templates:\n\
		{ return input.<expression>; }     -or-\n\
		{ return (<expression>); }\n\
\n\
	The resulting expression may evaluate to either an integer, float or string data type,\n\
	and FFSsort will attempt to sort based on that value.  Note that expressions applied \n\
	to different record types may result in different return types.  If all sort expressions\n\
	result in integer types, FFSsort will sort based upon integer values.  If sort \n\
	expressions are a mix of integer and float types, or are solely float, FFSsort\n\
	will coerce integers to floats and sort based upon floats.  String types may not be\n\
	mixed with numeric types, and their sort is based upon strcmp() results.\n\
\n";

static int verbose = 0;
static int quiet = 0;

extern char *
dump_raw_FMrecord_to_string(FMContext fmc, FMFormat format, void *data);
extern int
write_encoded_FFSfile(FFSFile f, void *data, DATA_LEN_TYPE byte_size, FFSContext c,
		      attr_list attrs);
extern int
FFSread_raw_header(FFSFile file, void *dest, int buffer_size, FFSTypeHandle *fp);

typedef struct {
    void *data;
    long size;
    double double_sort_entry;
    long int_sort_entry;
    char *string_sort_entry;
} data_entry;

void
write_data(data_entry *data_values, FFSContext c, FFSFile out_file)
{
    int i = 0;
    while (data_values[i].data != NULL) {
	write_encoded_FFSfile(out_file, data_values[i].data, 
			      data_values[i].size, c, NULL);
	free(data_values[i].data);
	i++;
    }
    free(data_values[i].data);
}

typedef double (*double_extraction_func)(void*data);
typedef int (*int_extraction_func)(void*data);
typedef char* (*string_extraction_func)(void*data);
typedef struct {
    FFSTypeHandle format;
    double_extraction_func dbl_func;
    int_extraction_func int_func;
    string_extraction_func str_func;
} func_table_entry;

static
void parse_expression(char *expr, char **format_name, char **expression)
{
    int not_simple = 0;
    int index;
    char *colon = strchr(expr, ':');
    if (colon) {
	*colon = 0;
	*format_name = strdup(expr);
	*colon = ':';
	expr = colon+1;
    }
    while (isspace(expr[0])) expr++;
    while(isspace(expr[strlen(expr)-1])) expr[strlen(expr)-1] = 0;
    if (!isalpha(expr[0])) not_simple++;
    index = 0;
    while ((!not_simple) && (expr[index] != 0)) {
	if (!isalnum(expr[index]) && (expr[index] != '_')) {
	    not_simple++;
	}
	index++;
    }
    *expression = expr;
}

static void null_err_func(void *client_data, char *string)
{
    /* do nothing */
}

static char* templates[] = {
    "{\nreturn input.%s;\n}",
    "{\nreturn (%s);\n}",
    NULL};

static int
generate_single_handler(FFSTypeHandle format, char *expression, func_table_entry *entry)
{
    char *prog = NULL;
    int i, j;
    FMFormat fmf = FMFormat_of_original(format);
    j = 0;
    while (templates[j] != NULL) {
	for (i=0 ; i < 3; i++) {
	    cod_parse_context context = new_cod_parse_context();
	    FMContext c = FMContext_from_FMformat(fmf);
	    cod_code gen_code;
	    int id_len;
	    char *decl;
	    switch(i) {
	    case 0:
		decl = "int proc()";
		break;
	    case 1:
		decl = "double proc()";
		break;
	    case 2:
		decl = "string proc()";
		break;
	    }
	    cod_subroutine_declaration(decl, context);
	    cod_add_encoded_param("input", get_server_ID_FMformat(fmf, &id_len), 0, c, context);
	    cod_set_error_func(context, null_err_func);
	    cod_set_dont_coerce_return(context, 1);
	    prog = malloc(strlen(templates[j]) + strlen(expression));
	    sprintf(prog, templates[j], expression);
	    gen_code = cod_code_gen(prog, context);
	    free(prog);
	    switch(i) {
	    case 0:
		if (gen_code) entry->int_func = (int_extraction_func) gen_code->func;
		break;
	    case 1:
		if (gen_code) entry->dbl_func = (double_extraction_func) gen_code->func;
		break;
	    case 2:
		if (gen_code) entry->str_func = (string_extraction_func) gen_code->func;
		break;
	    }
	    /* GSE fix leaks stuff! */
	    if (gen_code) {
		/* if anything worked, we're good, return */
		return 1;
	    }
	}
	j++;
    }
    return 0;
}

static void
generate_sort_handler(func_table_entry *entry, char **expressions, char** names)
{
    int i = 0;
    FFSTypeHandle format = entry->format;
    char *format_name = name_of_FMformat(FMFormat_of_original(format));
    
    /* first try all format-specific expressions */
    while(expressions[i] != NULL) {
	if (names[i] != NULL) {
	    if (strcmp(names[i], format_name) == 0) {
		if (!generate_single_handler(format, expressions[i], entry)) {
		    if (!quiet)
			printf("Warning, format specific handler \"%s\", failed to compile for format \"%s\"\n",
			       expressions[i], names[i]);
		}
		return;
	    }
	}
	i++;
    }
    /* first try all generic expressions (in order) */
    i = 0;
    while(expressions[i] != NULL) {
	if (names[i] == NULL) {
	    if (generate_single_handler(format, expressions[i], entry)) {
		return;
	    }
	}
	i++;
    }
    if (!quiet)	printf("No valid expression for format %s\n", format_name);
    return;
}

static func_table_entry *func_table = NULL;
static int func_table_count = 0;

static int
get_sort_handler_index(FFSTypeHandle format, char **expressions, char **names)
{
    int i;
    if (func_table == NULL) {
	func_table = malloc(sizeof(func_table[0]));
    }
    for (i=0; i < func_table_count; i++) {
	if (format == func_table[i].format) return i;
    }
    func_table = realloc(func_table, sizeof(func_table[0])*(func_table_count+2));
    func_table[func_table_count].format = format;
    generate_sort_handler(&func_table[func_table_count], expressions, names);
    return func_table_count++;
}

/* qsort struct comparision function (double field) */ 
int double_compar(const void *a, const void *b) 
{ 
    data_entry *ia = (data_entry *)a;
    data_entry *ib = (data_entry *)b;
    if (ia->double_sort_entry < ib->double_sort_entry) return -1;
    if (ia->double_sort_entry > ib->double_sort_entry) return 1;
    return 0;
} 
 
/* qsort struct comparision function (double field) */ 
int int_compar(const void *a, const void *b) 
{ 
    data_entry *ia = (data_entry *)a;
    data_entry *ib = (data_entry *)b;
    if (ia->int_sort_entry < ib->int_sort_entry) return -1;
    if (ia->int_sort_entry > ib->int_sort_entry) return 1;
    return 0;
} 
 
/* qsort struct comparision function (double field) */ 
int string_compar(const void *a, const void *b) 
{ 
    data_entry *ia = (data_entry *)a;
    data_entry *ib = (data_entry *)b;
    return strcmp(ia->string_sort_entry, ib->string_sort_entry);
} 
 
typedef enum {
    NO_SORT, INT_SORT, FLOAT_SORT, STRING_SORT
} sort_t;

#define combine(a,b) ((int)(a) * 10 + (int)(b))

int
main(int argc, char **argv)
{
    FFSFile in_file = NULL, out_file = NULL;
    int buffer_size = 1024;
    char *buffer = NULL;
    int i;
    int bitmap;
    char *in_filename = NULL, *out_filename = NULL;
    char **expressions = malloc(sizeof(expressions[0]));
    char **names = malloc(sizeof(expressions[0]));
    int expr_count = 0;
    int count = 0;
    data_entry *data_values = malloc(sizeof(data_entry));
    sort_t sort_type = NO_SORT;

    for (i = 1; i < argc; i++) {
	if (argv[i][0] == '-') {
    	    if (strcmp(argv[i], "-v") == 0) {
		verbose++;
	    } else if (strcmp(argv[i], "-q") == 0) {
		quiet++;
	    } else if (strcmp(argv[i], "-e") == 0) {
		if (argc <= i+1) {
		    fprintf(stderr, "Argument must follow -e\n");
		    fprintf(stderr, "%s", usage);
		    exit(0);
		}
		expressions = realloc(expressions, sizeof(expressions[0])*(expr_count+2));
		names = realloc(names, sizeof(names[0])*(expr_count+2));
		parse_expression(argv[i+1], &names[expr_count], &expressions[expr_count]);
		expr_count++;
		i++;
	    } else {
		fprintf(stderr, "Unknown argument \"%s\"\n%s\n", argv[i], usage);
	    }
	} else {
	    if (in_filename == NULL) {
		in_filename = argv[i];
	    } else if (out_filename == NULL) {
		out_filename = argv[i];
	    } else {
		fprintf(stderr, "Extra argument specified \"%s\"\n%s\n",
			argv[i], usage);
		exit(1);
	    }
	}
    }
    if (expr_count == 0) {
	fprintf(stderr, "At least one sort expression must be specified\n");
	fprintf(stderr, "%s", usage);
	exit(1);
    }
    if (!in_filename || !out_filename) {
	fprintf(stderr, "%s", usage);
	exit(1);
    }

    in_file = open_FFSfile(in_filename, "R");
    out_file = open_FFSfile(out_filename, "w");

    if (in_file == NULL) {
	fprintf(stderr, "File Open Failure \"%s\"", in_filename);
	perror("Opening input file");
	exit(1);
    }

    if (out_file == NULL) {
	fprintf(stderr, "File Open Failure \"%s\"", in_filename);
	perror("Opening output file");
	exit(1);
    }

    bitmap = FFSdata | FFSformat;
    FFSset_visible(in_file, bitmap);

    while (1) {
        switch (FFSnext_record_type(in_file)) {
	case FFSformat:{
	    FFSTypeHandle format;
	    int index;
	    sort_t this_item = NO_SORT;
	    format = FFSread_format(in_file);
	    index = get_sort_handler_index(format, expressions, names);
	    if (func_table[index].int_func) this_item = INT_SORT;
	    if (func_table[index].dbl_func) this_item = FLOAT_SORT;
	    if (func_table[index].str_func) this_item = STRING_SORT;
	    switch(combine(this_item, sort_type)) {
	    case combine(INT_SORT, NO_SORT):
	    case combine(FLOAT_SORT, NO_SORT):
	    case combine(STRING_SORT, NO_SORT):
		/* first format, go with what this one wants */
		sort_type = this_item;
		break;
	    case combine(INT_SORT, INT_SORT):
	    case combine(FLOAT_SORT, FLOAT_SORT):
	    case combine(STRING_SORT, STRING_SORT):
		/* not first, but no change */
		break;
	    case combine(STRING_SORT, INT_SORT):
	    case combine(STRING_SORT, FLOAT_SORT):
	    case combine(INT_SORT, STRING_SORT):
	    case combine(FLOAT_SORT, STRING_SORT):
		fprintf(stderr, "Incompatible SORT types.  Can't mix string expressions with non-string expressions.\n  SORT FAILED\n");
		exit(1);
		break;
	    case combine(FLOAT_SORT, INT_SORT):
		/* we were sorting based on INT values, but we're mixing doubles in so use those */
		sort_type = this_item;
		if (verbose) printf("Switching from INT to FLOAT sort\n");
		break;
	    case combine(INT_SORT, FLOAT_SORT):
		/* we were sorting based on double values, now adding ints, but we'll convert all ints to doubles */
		break;
	    case combine(NO_SORT, NO_SORT):
	    case combine(NO_SORT, INT_SORT):
	    case combine(NO_SORT, FLOAT_SORT):
	    case combine(NO_SORT, STRING_SORT):
		/* new item has no expression, no change */
		break;
	    default:
		printf("Didn't think of this case %d\n", combine(this_item, sort_type));
	    }
	}
        case FFSdata:{
	    FFSTypeHandle format;
	    int size = FFSnext_data_length(in_file);
	    int index;
	    buffer = malloc(size);
	    FFSread_raw_header(in_file, buffer, buffer_size, &format);
	    index = get_sort_handler_index(format, expressions, names);
	    if (func_table[index].int_func) {
		data_values[count].int_sort_entry = func_table[index].int_func(buffer);
		data_values[count].double_sort_entry = data_values[count].int_sort_entry;
		data_values[count].string_sort_entry = NULL;
		if (verbose) printf("Int sort value is %ld\n", data_values[count].int_sort_entry);
	    } else if (func_table[index].dbl_func) {
		data_values[count].double_sort_entry = func_table[index].dbl_func(buffer);
		data_values[count].int_sort_entry = data_values[count].double_sort_entry;
		data_values[count].string_sort_entry = NULL;
		if (verbose) printf("Double sort value is %g\n", data_values[count].double_sort_entry);
	    } else if (func_table[index].str_func) {
		data_values[count].int_sort_entry = 0;
		data_values[count].double_sort_entry = 0;
		data_values[count].string_sort_entry = func_table[index].str_func(buffer);
		if (verbose) printf("String sort value is \"%s\"\n", data_values[count].string_sort_entry);
	    } else {
		if (verbose) printf("No sort value found\n");
	    }
	    data_values[count].data = buffer;
	    data_values[count].size = size;
	    data_values = realloc(data_values, sizeof(data_values[0])*(count+2));
	    count++;
	    data_values[count].data = NULL;
	    break;
        }
        case FFSerror:
        case FFSend:
	    switch(sort_type) {
	    case INT_SORT:
		qsort(data_values, count, sizeof(data_values[0]), int_compar);
		break;
	    case FLOAT_SORT:
		qsort(data_values, count, sizeof(data_values[0]), double_compar);
		break;
	    case STRING_SORT:
		qsort(data_values, count, sizeof(data_values[0]), string_compar);
		break;
	    default:
	      assert(0);
	    }
	    write_data(data_values, FFSContext_of_file(in_file), out_file);
            close_FFSfile(in_file);
            close_FFSfile(out_file);
            exit(0);
	default:
	    break;
        }
    }
}
