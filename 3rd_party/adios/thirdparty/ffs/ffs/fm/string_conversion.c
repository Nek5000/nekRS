#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "assert.h"
#include "fm.h"

#ifndef HAVE_STRTOF
static float
strtof(char *str, char **endptr)
{
    float float_part = 0.0;
    int integer_part = 0, i = 10, decimal = 0;
    float ret;
    char *s = str;
    int negative = 0;
    /* Skip white space.  */
    while (isspace(*s))
	++s;
    if (*s == '-') {
	negative = 1;
	++s;
    } else if (*s == '+')
	++s;

    while (*s != '\0') {
	int t = *s - '0';
	if (*s == '.') {
	    decimal = 1;
	    s++;
	    continue;
	}
	if (t < 0 || t > 9)
	    break;
	if (decimal) {
	    float_part = float_part + ((float) t / i);
	    i *= 10;
	} else
	    integer_part = (integer_part * 10) + t;
	s++;
    }
    *endptr = s;
    ret = (float) integer_part + float_part;
    return negative ? -ret : ret;
}
#endif

#ifndef HAVE_STRTOD
static double
strtod(char *str, char **endptr)
{
    double float_part = 0.0;
    int integer_part = 0, i = 10, decimal = 0;
    double ret;
    char *s = str;
    int negative = 0;
    /* Skip white space.  */
    while (isspace(*s))
	++s;
    if (*s == '-') {
	negative = 1;
	++s;
    } else if (*s == '+')
	++s;
    while (*s != '\0') {
	int t = *s - '0';
	if (*s == '.') {
	    decimal = 1;
	    s++;
	    continue;
	}
	if (t < 0 || t > 9)
	    break;
	if (decimal) {
	    float_part = float_part + ((double) t / i);
	    i *= 10;
	} else
	    integer_part = (integer_part * 10) + t;
	s++;
    }
    *endptr = s;
    ret = (double) integer_part + float_part;
    return negative ? -ret : ret;
}
#endif

/* TODO: Fix strtold -sandip */
#ifndef HAVE_STRTOLD
static long double
strtold(char *str, char **endptr)
{
    long double float_part = 0.0;
    long integer_part = 0;
    int i = 10, decimal = 0;
    long double ret;
    char *s = str;
    int negative = 0;
    /* Skip white space.  */
    while (isspace(*s))
	++s;
    if (*s == '-') {
	negative = 1;
	++s;
    } else if (*s == '+')
	++s;
    while (*s != '\0') {
	int t = *s - '0';
	if (*s == '.') {
	    decimal = 1;
	    s++;
	    continue;
	}
	if (t < 0 || t > 9)
	    break;
	if (decimal) {
	    float_part = float_part + ((long double) t / i);
	    i *= 10;
	} else
	    integer_part = (integer_part * 10) + t;
	s++;
    }
    *endptr = s;
    ret = (long double) integer_part + float_part;
    return negative ? -ret : ret;
}
#endif

extern void
str_to_val(char *val_str, FMdata_type typ, int dest_size, void **default_value)
{
    char *error_check;
    void *default_val = NULL;
    switch (typ) {
    case integer_type:
    case enumeration_type:
	{
	    long tmp = strtol(val_str, &error_check, 10);
	    if (!error_check || *error_check != '\0')
		break;
	    default_val = malloc(dest_size);
	    if (dest_size == sizeof(char)) {
		char *dest_field = (char *) default_val;
		*dest_field = (char) tmp;
	    } else if (dest_size == sizeof(short)) {
		short *dest_field = (short *) default_val;
		*dest_field = (short) tmp;
	    } else if (dest_size == sizeof(int)) {
		int *dest_field = (int *) default_val;
		*dest_field = tmp;
	    } else if (dest_size == sizeof(long)) {
		long *dest_field = (long *) default_val;
		*dest_field = tmp;
#if SIZEOF_LONG_LONG != 0
	    } else if (dest_size == sizeof(long long)) {
		long long lltmp = strtoll(val_str, &error_check, 10);
		if (!error_check || *error_check != '\0') {
		    free(default_val);
		    default_val = NULL;
		    break;
		}
		memcpy(default_val, &lltmp, sizeof(long long));
#endif
	    } else {
		printf("Size problem in Default value conversion\n");
	    }
	}
	break;
    case unsigned_type:
    case boolean_type:
	{
	    long tmp = strtoul(val_str, &error_check, 10);
	    if (!error_check || *error_check != '\0')
		break;
	    default_val = malloc(dest_size);
	    if (dest_size == sizeof(unsigned char)) {
		unsigned char *dest_field = (unsigned char *) default_val;
		*dest_field = (unsigned char) tmp;
	    } else if (dest_size == sizeof(unsigned short)) {
		unsigned short *dest_field = (unsigned short *) default_val;
		*dest_field = (unsigned short) tmp;
	    } else if (dest_size == sizeof(unsigned int)) {
		unsigned int *dest_field = (unsigned int *) default_val;
		*dest_field = tmp;
	    } else if (dest_size == sizeof(unsigned long)) {
		unsigned long *dest_field = (unsigned long *) default_val;
		*dest_field = tmp;
#if SIZEOF_LONG_LONG != 0
	    } else if (dest_size == sizeof(unsigned long long)) {
		unsigned long long *dest_field = (unsigned long long *) default_val;
		unsigned long long tmp = (unsigned long long) strtoll(val_str, &error_check, 10);
		if (!error_check || *error_check != '\0') {
		    free(default_val);
		    default_val = NULL;
		    break;
		}
		*dest_field = tmp;
#endif
	    } else {
		printf("Size problem in Default value conversion\n");
	    }
	}
	break;
    case float_type:
	{
	    if (dest_size == sizeof(float)) {
		float *dest_field;
		float tmp = strtof(val_str, &error_check);
		if (error_check && *error_check == '\0') {
		    default_val = malloc(dest_size);
		    dest_field = (float *) default_val;
		    *dest_field = tmp;
		}
	    } else if (dest_size == sizeof(double)) {
		double *dest_field;
		double tmp = strtod(val_str, &error_check);
		if (error_check && *error_check == '\0') {
		    default_val = malloc(dest_size);
		    dest_field = (double *) default_val;
		    *dest_field = tmp;
		}
#if SIZEOF_LONG_DOUBLE != 0
	    } else if (dest_size == sizeof(long double)) {
		long double *dest_field;
		long double tmp = strtof(val_str, &error_check);
		if (error_check && *error_check == '\0') {
		    default_val = malloc(dest_size);
		    dest_field = (long double *) default_val;
		    *dest_field = tmp;
		}
#endif
	    } else {
		printf("Size problem in Default value conversion\n");
	    }
	}
	break;
    case char_type:
	{
	    char *dest_field;
	    default_val = malloc(dest_size);
	    dest_field = (char *) default_val;
	    *dest_field = *val_str;
	    break;
	}
    case string_type:
	printf("Default values for string types not yet supported\n");
	break;
    case unknown_type:
	printf("Default values for complex types not yet supported\n");
	break;
    default:
	assert(0);
    }
    *default_value = default_val;
}
