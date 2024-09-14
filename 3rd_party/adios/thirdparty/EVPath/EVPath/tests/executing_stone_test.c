/*
 * executing_stone_test exists to test the EVexecuting_stone()
 * functionality.  I.E. it attempts to determine on which stone's behalf a
 * handler is executing.  We do use client_data to cross-check here, but
 * EVexecuting_stone() is designed for use where client_data wouldn't work.
 */

#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#ifdef HAVE_ARPA_INET_H
#include <arpa/inet.h>
#endif
#include "evpath.h"
#ifdef HAVE_SYS_WAIT_H
#include <sys/wait.h>
#endif
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#endif


typedef struct _complex_rec {
    double r;
    double i;
} complex, *complex_ptr;

typedef struct _nested_rec {
    complex item;
} nested, *nested_ptr;

static FMField nested_field_list[] =
{
    {"item", "complex", sizeof(complex), FMOffset(nested_ptr, item)},
    {NULL, NULL, 0, 0}
};

static FMField complex_field_list[] =
{
    {"r", "double", sizeof(double), FMOffset(complex_ptr, r)},
    {"i", "double", sizeof(double), FMOffset(complex_ptr, i)},
    {NULL, NULL, 0, 0}
};

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;

static FMField simple_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {"short_field", "integer",
     sizeof(short), FMOffset(simple_rec_ptr, short_field)},
    {"long_field", "integer",
     sizeof(long), FMOffset(simple_rec_ptr, long_field)},
    {"nested_field", "nested",
     sizeof(nested), FMOffset(simple_rec_ptr, nested_field)},
    {"double_field", "float",
     sizeof(double), FMOffset(simple_rec_ptr, double_field)},
    {"char_field", "char",
     sizeof(char), FMOffset(simple_rec_ptr, char_field)},
    {"scan_sum", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, scan_sum)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {"complex", complex_field_list, sizeof(complex), NULL},
    {"nested", nested_field_list, sizeof(nested), NULL},
    {NULL, NULL}
};

int quiet = 1;

EVsource source_handle[3];

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    (void)cm;
    if (EVexecuting_stone(cm) != *(int*)client_data) {
	printf("Inconsistency in EVexecuting_stone, expected %d, got %d!\n",
	       *(int*)client_data, EVexecuting_stone(cm));
	exit(1);
    }
    if (quiet <= 0) printf("In handler for stone %d\n", EVexecuting_stone(cm));
    if (event->integer_field != 0) {
	event->integer_field--;
	EVsubmit(source_handle[lrand48()%3], event, NULL);
    }
    return 0;
}

char *transport = NULL;
char *control = NULL;

#include "support.c"
int regression_master, regression;  /* not used in this test */

int
main(int argc, char **argv)
{
    CManager cm;
    EVstone sink1, sink2, sink3;
    simple_rec data;

    PARSE_ARGS();

    srand48(getpid());

    cm = CManager_create();
    (void) CMfork_comm_thread(cm);


    sink1 = EValloc_stone(cm);
    sink2 = EValloc_stone(cm);
    sink3 = EValloc_stone(cm);
    EVassoc_terminal_action(cm, sink1, simple_format_list, simple_handler, (void*)&sink1);
    EVassoc_terminal_action(cm, sink2, simple_format_list, simple_handler, (void*)&sink2);
    EVassoc_terminal_action(cm, sink3, simple_format_list, simple_handler, (void*)&sink3);

    source_handle[0] = EVcreate_submit_handle(cm, sink1, simple_format_list);
    source_handle[1] = EVcreate_submit_handle(cm, sink2, simple_format_list);
    source_handle[2] = EVcreate_submit_handle(cm, sink3, simple_format_list);
    
    memset(&data, 0, sizeof(data));
    data.integer_field = 10;   /* keep submitting to random stone until this is zero */
    if (quiet <= 0) printf("submitting %d\n", data.integer_field);
    EVsubmit(source_handle[lrand48()%3], &data, NULL);
    if (EVexecuting_stone(cm) != -1) {
	printf("Bad return from invalid call to EVexecuting_stone()\n");
	return 1;
    }

    EVfree_source(source_handle[0]);
    EVfree_source(source_handle[1]);
    EVfree_source(source_handle[2]);

    CManager_close(cm);
    return 0;
}

