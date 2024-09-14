#include "config.h"
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
#include "cod.h"
#include "cod_internal.h"
#include "structs.h"
#undef NDEBUG
#include "assert.h"
#include <ctype.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>

#ifndef _MSC_VER
#include <sys/time.h>
#else
#define strdup _strdup
#endif
#ifndef LINUX_KERNEL_MODULE
#ifdef HAVE_ATL_H
#include "atl.h"
#endif
#endif
#include "ffs.h"

#ifndef LINUX_KERNEL_MODULE
#ifdef HAVE_ATL_H
static int
attr_set(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    attr_value_type junk;
    attr_value junk2;
    if (atom == 0 ) return 0;
    
    return query_attr(l, atom, &junk, &junk2);
}

static attr_list
attr_create_list()
{
    return create_attr_list();
}

static void
attr_free_list(attr_list l)
{
    free_attr_list(l);
}

static void
std_set_int_attr(attr_list l, char *name, int value)
{
    atom_t atom = attr_atom_from_string(name);
    if (atom == 0 ) return;

    set_int_attr(l, atom, value);
}

static void
std_set_long_attr(attr_list l, char *name, long value)
{
    atom_t atom = attr_atom_from_string(name);
    if (atom == 0 ) return;

    set_long_attr(l, atom, value);
}

static void
std_set_double_attr(attr_list l, char *name, double value)
{
    atom_t atom = attr_atom_from_string(name);
    if (atom == 0 ) return;

    set_double_attr(l, atom, value);
}

static void
std_set_float_attr(attr_list l, char *name, float value)
{
    atom_t atom = attr_atom_from_string(name);
    if (atom == 0 ) return;

    set_float_attr(l, atom, value);
}

static void
std_set_string_attr(attr_list l, char *name, char *value)
{
    atom_t atom = attr_atom_from_string(name);
    if (atom == 0 ) return;

    set_string_attr(l, atom, value);
}

static int
attr_ivalue(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    int i = 0;
    if (atom == 0 ) return 0;
    
    get_int_attr(l, atom, &i);
    return i;
}

static ssize_t
attr_lvalue(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    ssize_t lo = 0;
    if (atom == 0 ) return 0;
    
    get_long_attr(l, atom, &lo);
    return lo;
}

static double
attr_dvalue(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    double d;
    if (atom == 0 ) return 0;
    
    get_double_attr(l, atom, &d);
    return d;
}

static float
attr_fvalue(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    float f;
    if (atom == 0 ) return 0;
    
    get_float_attr(l, atom, &f);
    return f;
}

static char *
attr_svalue(attr_list l, char *name)
{
    atom_t atom = attr_atom_from_string(name);
    char *s;
    if (atom == 0 ) return 0;
    
    get_string_attr(l, atom, &s);
    return strdup(s);
}
#endif

static FFSFile 
open_ffs_file(char * fname, char * mode)
{
    FFSFile temp;
    temp = open_FFSfile(fname, mode);
    if(!temp) {
	fprintf(stderr, "Could not open FFSfile from CoD\n");
    }
    return temp;
}

static void close_ffs_file(FFSFile fname)
{
    close_FFSfile(fname);
}

typedef struct chr_time {
    double d1;
    double d2;
    double d3;
} chr_time;

#ifdef _MSC_VER
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 

#ifndef _WINSOCKAPI_
// MSVC defines this in winsock2.h!?
typedef struct timeval {
    long tv_sec;
    long tv_usec;
} timeval;
#endif
typedef struct timezone {
    int tz_minuteswest;     /* minutes west of Greenwich */
    int tz_dsttime;         /* type of DST correction */
} timezone;

int gettimeofday(struct timeval* tp, struct timezone* tzp)
{
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    // This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
    // until 00:00:00 January 1, 1970 
    static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

    SYSTEMTIME  system_time;
    FILETIME    file_time;
    uint64_t    time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    time = ((uint64_t)file_time.dwLowDateTime);
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec = (long)((time - EPOCH) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
}
#else
#include <sys/time.h>
#endif

int
gettimeofday_wrapper(struct timeval* tp)
{
    int ret = gettimeofday(tp, NULL);
    return ret;
}

static void
chr_get_time( chr_time *time)
{
    gettimeofday((struct timeval *) time, NULL);
}

static void
chr_timer_start( chr_time *time)
{
    chr_get_time(time);
}

static void chr_timer_diff( chr_time *diff_time, chr_time *src1, chr_time *src2);
static void
chr_timer_stop( chr_time *time)
{
    struct timeval now;
    struct timeval duration;

    gettimeofday(&now, NULL);
    chr_timer_diff((chr_time*)&duration, (chr_time*)&now, time);
    *((struct timeval *) time) = duration;
}

static int
chr_timer_eq_zero (chr_time *time)
{
    struct timeval *t = (struct timeval *) time; 
    return ((t->tv_sec == 0) && (t->tv_usec == 0));
}

static void
chr_timer_diff( chr_time *diff, chr_time *src1, chr_time *src2)
{
    struct timeval d;
    struct timeval *s1 = (struct timeval *)src1;
    struct timeval *s2 = (struct timeval *)src2;
    d.tv_sec = s1->tv_sec - s2->tv_sec;
    d.tv_usec = s1->tv_usec - s2->tv_usec;
    if (d.tv_usec < 0) {
	d.tv_usec += 1000000;
	d.tv_sec--;
    }
    *((struct timeval*)diff) = d;
}

static void
chr_timer_sum( chr_time *sum, chr_time *src1, chr_time *src2)
{
    struct timeval s;
    struct timeval *s1 = (struct timeval *)src1;
    struct timeval *s2 = (struct timeval *)src2;
    s.tv_sec = s1->tv_sec + s2->tv_sec;
    s.tv_usec = s1->tv_usec + s2->tv_usec;
    if (s.tv_usec > 1000000) {
	s.tv_usec -= 1000000;
	s.tv_sec++;
    }
    *((struct timeval*)sum) = s;
}


static double
chr_time_to_secs(chr_time *time)
{
    return (double)((struct timeval*)time)->tv_sec + 
	((double)((struct timeval*)time)->tv_usec)/1000000.0;
}

static double
chr_time_to_millisecs(chr_time *time)
{
    return ((double)((struct timeval*)time)->tv_sec)*1000.0 + 
	((double)((struct timeval*)time)->tv_usec)/1000.0;
}

static double
chr_time_to_microsecs(chr_time *time)
{
    return ((double)((struct timeval*)time)->tv_sec)*1000000.0 + 
	((double)((struct timeval*)time)->tv_usec);
}

static double
chr_time_to_nanosecs(chr_time *time)
{
    return ((double)((struct timeval*)time)->tv_sec)*1000000000.0 + 
	((double)((struct timeval*)time)->tv_usec*1000.0);
}

static double
chr_approx_resolution()
{
    struct timeval start, stop, diff;
    gettimeofday(&start, NULL);
    gettimeofday(&stop, NULL);
    while(start.tv_usec == stop.tv_usec) {
	gettimeofday(&stop, NULL);
    }
    chr_timer_diff((chr_time*)&diff, (chr_time*)&stop, (chr_time*)&start);
    return chr_time_to_secs((chr_time*)&diff);
}

static char atl_extern_string[] = "\n\
	int attr_set(attr_list l, string name);\n\
	attr_list create_attr_list();\n\
	attr_list copy_attr_list(attr_list l);\n\
	void free_attr_list(attr_list l);\n					\
	void set_long_attr(attr_list l, string name, long value);\n\
	void set_float_attr(attr_list l, string name, double value);\n\
	void set_double_attr(attr_list l, string name, double value);\n\
	void set_int_attr(attr_list l, string name, int value);\n\
	void set_string_attr(attr_list l, string name, string value);\n\
	int attr_ivalue(attr_list l, string name);\n\
	long attr_lvalue(attr_list l, string name);\n\
	double attr_dvalue(attr_list l, string name);\n\
	double attr_fvalue(attr_list l, string name);\n\
	char* attr_svalue(attr_list l, string name);\n";
static char chr_extern_string[] = "\n\
        void chr_get_time( chr_time *time);\n\
        void chr_timer_diff( chr_time *diff_time, chr_time *src1, chr_time *src2);\n\
	int chr_timer_eq_zero( chr_time *time);\n\
	void chr_timer_sum( chr_time *sum_time, chr_time *src1, chr_time *src2);\n\
	void chr_timer_start (chr_time *timer);\n\
	void chr_timer_stop (chr_time *timer);\n\
	double chr_time_to_nanosecs (chr_time *time);\n\
	double chr_time_to_microsecs (chr_time *time);\n\
	double chr_time_to_millisecs (chr_time *time);\n\
	double chr_time_to_secs (chr_time *time);\n";

static char basic_extern_string[] = "\n\
	double chr_approx_resolution();\n\
	int gettimeofday(timeval *tp);\n\
	ffs_file open_ffs(char * fname, char * mode);\n\
	void close_ffs(ffs_file fname);\n";

static char internals[] = "\n\
	void cod_NoOp(int duration);\n";

static cod_extern_entry internal_externs[] = 
{
    {"cod_NoOp", (void*)(long)0xdeadbeef},    /* value is unimportant, but can't be NULL */
    {NULL, NULL}
};

static cod_extern_entry externs[] = 
{
#ifdef HAVE_ATL_H
    {"attr_set", (void*)(intptr_t)attr_set},
    {"create_attr_list", (void*)(intptr_t)attr_create_list},
    {"copy_attr_list", (void*)(intptr_t)attr_copy_list},
    {"free_attr_list", (void*)(intptr_t)attr_free_list},
    {"set_int_attr", (void*)(intptr_t)std_set_int_attr},
    {"set_long_attr", (void*)(intptr_t)std_set_long_attr},
    {"set_double_attr", (void*)(intptr_t)std_set_double_attr},
    {"set_float_attr", (void*)(intptr_t)std_set_float_attr},
    {"set_string_attr", (void*)(intptr_t)std_set_string_attr},
    {"attr_ivalue", (void*)(intptr_t)attr_ivalue},
    {"attr_lvalue", (void*)(intptr_t)attr_lvalue},
    {"attr_dvalue", (void*)(intptr_t)attr_dvalue},
    {"attr_fvalue", (void*)(intptr_t)attr_fvalue},
    {"attr_svalue", (void*)(intptr_t)attr_svalue},
#endif
    {"chr_get_time", (void*)(intptr_t)chr_get_time},
    {"chr_timer_diff", (void*)(intptr_t)chr_timer_diff},
    {"chr_timer_eq_zero", (void*)(intptr_t)chr_timer_eq_zero},
    {"chr_timer_sum", (void*)(intptr_t)chr_timer_sum},
    {"chr_timer_start", (void*)(intptr_t)chr_timer_start},
    {"chr_timer_stop", (void*)(intptr_t)chr_timer_stop},
    {"chr_time_to_nanosecs", (void*)(intptr_t)chr_time_to_nanosecs},
    {"chr_time_to_microsecs", (void*)(intptr_t)chr_time_to_microsecs},
    {"chr_time_to_millisecs", (void*)(intptr_t)chr_time_to_millisecs},
    {"chr_time_to_secs", (void*)(intptr_t)chr_time_to_secs},
    {"chr_approx_resolution", (void*)(intptr_t)chr_approx_resolution},
    {"gettimeofday", (void*)(intptr_t)gettimeofday_wrapper},
    {"open_ffs", (void*)(intptr_t)open_ffs_file},
    {"close_ffs", (void*)(intptr_t)close_ffs_file},
    {(void*)0, (void*)0}
};

FMField chr_time_list[] = {
    {"d1", "double", sizeof(double), FMOffset(chr_time*, d1)}, 
    {"d2", "double", sizeof(double), FMOffset(chr_time*, d2)}, 
    {"d3", "double", sizeof(double), FMOffset(chr_time*, d3)}, 
    {NULL, NULL, 0, 0}};

FMField timeval_list[] = {
    {"tv_sec", "integer", sizeof(((struct timeval*)0)->tv_sec), FMOffset(struct timeval *, tv_sec)}, 
    {"tv_usec", "integer", sizeof(((struct timeval*)0)->tv_usec), FMOffset(struct timeval *, tv_usec)}, 
    {NULL, NULL, 0, 0}};

extern void
cod_add_standard_elements(cod_parse_context context)
{
    cod_assoc_externs(context, externs);
#ifdef HAVE_ATL_H
    sm_ref attr_node = cod_new_reference_type_decl();
    attr_node->node.reference_type_decl.name = strdup("attr_list");
    cod_add_decl_to_parse_context("attr_list", attr_node, context);
    cod_add_decl_to_scope("attr_list", attr_node, context);
    cod_add_defined_type("attr_list", context);
    cod_parse_for_context(atl_extern_string, context);
#endif
    sm_ref ffs_node = cod_new_reference_type_decl();
    ffs_node->node.reference_type_decl.name = strdup("ffs_file");
    cod_add_decl_to_parse_context("ffs_file", ffs_node, context);
    cod_add_decl_to_scope("ffs_file", ffs_node, context);
    cod_add_defined_type("ffs_file", context);

    cod_add_int_constant_to_parse_context("NULL", 0, context);
    cod_add_simple_struct_type("chr_time", chr_time_list, context);
    cod_parse_for_context(chr_extern_string, context);
    cod_add_simple_struct_type("timeval", timeval_list, context);
    cod_add_defined_type("cod_type_spec", context);
    cod_add_defined_type("cod_exec_context", context);
    cod_add_defined_type("cod_closure_context", context);
    cod_semanticize_added_decls(context);
    
    cod_parse_for_context(basic_extern_string, context);

    cod_assoc_externs(context, internal_externs);
    cod_parse_for_context(internals, context);
    cod_swap_decls_to_standard(context);
}

#else /* LINUX_KERNEL_MODULE */

extern void
cod_add_standard_elements(cod_parse_context context)
{
}
#endif /* LINUX_KERNEL_MODULE */

#if defined(NO_DYNAMIC_LINKING) && !defined(_MSC_VER)
#define sym(x) (void*)(intptr_t)x
#else
#define sym(x) (void*)0
#endif

static cod_extern_entry string_externs[] = 
{
    {"memchr", (void*)(intptr_t)memchr},
    {"memcmp", (void*)(intptr_t)memcmp},
    {"memcpy", (void*)(intptr_t)memcpy},
    {"memmove", (void*)(intptr_t)memmove},
    {"memset", (void*)(intptr_t)memset},
    {"strcat", (void*)(intptr_t)strcat},
    {"strchr", (void*)(intptr_t)strchr},
    {"strcmp", (void*)(intptr_t)strcmp},
    {"strcoll", (void*)(intptr_t)strcoll},
    {"strcpy", (void*)(intptr_t)strcpy},
    {"strcspn", (void*)(intptr_t)strcspn},
    {"strerror", (void*)(intptr_t)strerror},
    {"strlen", (void*)(intptr_t)strlen},
    {"strncat", (void*)(intptr_t)strncat},
    {"strncmp", (void*)(intptr_t)strncmp},
    {"strncpy", (void*)(intptr_t)strncpy},
    {"strpbrk", (void*)(intptr_t)strpbrk},
    {"strrchr", (void*)(intptr_t)strrchr},
    {"strspn", (void*)(intptr_t)strspn},
    {"strstr", (void*)(intptr_t)strstr},
    {"strtok", (void*)(intptr_t)strtok},
    {"strxfrm", (void*)(intptr_t)strxfrm},
    {NULL, NULL}
};

static char string_extern_string[] = "\n\
void	*memchr(const void *s, int c, int size);\n\
int	 memcmp(const void *m, const void *s, int size);\n\
void	*memcpy(void *m, const void *s, int size);\n\
void	*memmove(void *m, const void *s, int size);\n\
void	*memset(void *m, int c, int size);\n\
char	*strcat(char *s1, const char *s2);\n\
char	*strchr(const char *s1, int c);\n\
int	 strcmp(const char *s1, const char *s2);\n\
int	 strcoll(const char *s1, const char *s2);\n\
char	*strcpy(char *s1, const char *s2);\n\
int	 strcspn(const char *s1, const char *s2);\n\
int	 strlen(const char *s);\n\
char	*strncat(char *s1, const char *s2, int s);\n\
int	 strncmp(const char *s1, const char *s2, int s);\n\
char	*strncpy(char *s1, const char *s2, int s);\n\
char	*strpbrk(const char *s1, const char *s2);\n\
char	*strrchr(const char *s1, int c);\n\
int	 strspn(const char *s1, const char *s2);\n\
char	*strstr(const char *s1, const char *s2);\n\
char	*strtok(char *s1, const char *s2);\n\
int	 strxfrm(char *s1, const char *s2, int size);\n\
";

static cod_extern_entry strings_externs[] = 
{
#ifndef _MSC_VER
    {"bcmp", (void*)(intptr_t)bcmp},
    {"bcopy", (void*)(intptr_t)bcopy},
    {"bzero", (void*)(intptr_t)bzero},
    {"index", (void*)(intptr_t)index},
    {"rindex", (void*)(intptr_t)rindex},
    {"ffs", (void*)(intptr_t)ffs},
    {"strcasecmp", (void*)(intptr_t)strcasecmp},
    {"strncasecmp", (void*)(intptr_t)strncasecmp},
#endif
    {NULL, NULL}
};

static char strings_extern_string[] = "\n\
int	 bcmp(const void *m1, const void *m2, int size);\n\
void	 bcopy(const void *m1, void *m2, int size);\n\
void	 bzero(void *m, int size);\n\
char	*index(const char *s1, int c);\n\
char	*rindex(const char *s1, int c);\n\
int	 ffs(int);\n\
int	 strcasecmp(const char *s1, const char *s2);\n\
int	 strncasecmp(const char *s1, const char *s2, int size);\n\
";


#include <math.h>

static cod_extern_entry math_externs[] = 
{
    {"acos", sym(acos)},
    {"asin", sym(asin)},
    {"atan", sym(atan)},
    {"atan2", sym(atan2)},
    {"cos", sym(cos)},
    {"sin", sym(sin)},
    {"tan", sym(tan)},
    {"acosh", sym(acosh)},
    {"asinh", sym(asinh)},
    {"atanh", sym(atanh)},
    {"cosh", sym(cosh)},
    {"sinh", sym(sinh)},
    {"tanh", sym(tanh)},
    {"exp", sym(exp)},
    {"exp2", sym(exp2)},
    {"expm1", sym(expm1)},
    {"log", sym(log)},
    {"log10", sym(log10)},
    {"log2", sym(log2)},
    {"log1p", sym(log1p)},
    {"logb", sym(logb)},
    {"modf", sym(modf)},
    {"ldexp", sym(ldexp)},
    {"frexp", sym(frexp)},
    {"ilogb", sym(ilogb)},
    {"scalbn", sym(scalbn)},
    {"scalbln", sym(scalbln)},
    {"fabs", sym(fabs)},
    {"cbrt", sym(cbrt)},
    {"hypot", sym(hypot)},
    {"pow", sym(pow)},
    {"sqrt", sym(sqrt)},
    {"erf", sym(erf)},
    {"erfc", sym(erfc)},
    {"lgamma", sym(lgamma)},
    {"tgamma", sym(tgamma)},
    {"ceil", sym(ceil)},
    {"floor", sym(floor)},
    {"nearbyint", sym(nearbyint)},
    {"rint", sym(rint)},
    {"lrint", sym(lrint)},
    {"round", sym(round)},
    {"lround", sym(lround)},
    {"trunc", sym(trunc)},
    {"fmod", sym(fmod)},
    {"remainder", sym(remainder)},
    {"remquo", sym(remquo)},
    {"copysign", sym(copysign)},
    {"nan", sym(nan)},
    {NULL, NULL}
};

static char math_extern_string[] = "\n\
double acos(double a);\n\
double asin(double a);\n\
double atan(double a);\n\
double atan2(double b, double a);\n\
double cos(double a);\n\
double sin(double a);\n\
double tan(double a);\n\
double acosh(double a);\n\
double asinh(double a);\n\
double atanh(double a);\n\
double cosh(double a);\n\
double sinh(double a);\n\
double tanh(double a);\n\
double exp(double a);\n\
double exp2(double a); \n\
double expm1(double a); \n\
double log(double a);\n\
double log10(double a);\n\
double log2(double a);\n\
double log1p(double a);\n\
double logb(double a);\n\
double modf(double b, double * a);\n\
double ldexp(double b, int a);\n\
double frexp(double b, int * a);\n\
int ilogb(double a);\n\
double scalbn(double b, int a);\n\
double scalbln(double b, long int a);\n\
double fabs(double a);\n\
double cbrt(double a);\n\
double hypot(double b, double a);\n\
double pow(double b, double a);\n\
double sqrt(double a);\n\
double erf(double a);\n\
double erfc(double a);\n\
double lgamma(double a);\n\
double tgamma(double a);\n\
double ceil(double a);\n\
double floor(double a);\n\
double nearbyint(double a);\n\
double rint(double a);\n\
long   lrint(double a);\n\
double round(double a);\n\
long   lround(double a);\n\
double trunc(double a);\n\
double fmod(double a, double b);\n\
double remainder(double a, double b);\n\
double remquo(double a, double b, int *c);\n\
double copysign(double a, double b);\n\
double nan(const char * a);\n\
";


#include <limits.h>

static char limits_extern_string[] = "\n\
const char SCHAR_MAX = 127;\n\
const char SCHAR_MIN = -128;\n\
\n\
const unsigned char UCHAR_MAX = 255;\n\
const char CHAR_MAX = 127;\n\
const char CHAR_MIN = (-128);\n\
\n\
const unsigned short USHRT_MAX = 65535;\n\
const short SHRT_MAX = 32767;\n\
const short SHRT_MIN = (-32768);\n\
\n\
const unsigned int	UINT_MAX = 0xffffffff;\n\
const int INT_MAX = 2147483647;\n\
const int INT_MIN = (-2147483647-1);\n\
const long LONG_MAX = 0x7fffffffffffffffL;\n\
const long LONG_MIN = (-0x7fffffffffffffffL-1);\n\
const unsigned long ULONG_MAX = 0xffffffffffffffffUL;\n\
";


static void dlload_externs(char *libname, cod_extern_entry *externs);

extern void
cod_process_include(char *name, cod_parse_context context)
{
    intptr_t char_count = strchr(name, '.') - name;
    if (char_count < 0) char_count = strlen(name);
    if (strncmp(name, "string", char_count) == 0) {
	cod_assoc_externs(context, string_externs);
	cod_parse_for_context(string_extern_string, context);
    } else if (strncmp(name, "strings", char_count) == 0) {
	cod_assoc_externs(context, strings_externs);
	cod_parse_for_context(strings_extern_string, context);
    } else if (strncmp(name, "math", char_count) == 0) {
	dlload_externs("libm", math_externs);
	cod_assoc_externs(context, math_externs);
	cod_parse_for_context(math_extern_string, context);
    } else if (strncmp(name, "limits", char_count) == 0) {
	cod_parse_for_context(limits_extern_string, context);
    }

}
#if !NO_DYNAMIC_LINKING
#include <dlfcn.h>
#endif

static void 
dlload_externs(char *libname, cod_extern_entry *externs)
{
#if NO_DYNAMIC_LINKING
    return;
#else
    int i = 0;
    char *name = malloc(strlen(libname) + strlen(LIBRARY_EXT) + 1);
    strcpy(name, libname);
    strcat(name, LIBRARY_EXT);
    void *handle = dlopen(name, RTLD_LAZY);
    free(name);
    while(externs[i].extern_name) {
	externs[i].extern_value = dlsym(handle, externs[i].extern_name);
	i++;
    }
#endif
}
