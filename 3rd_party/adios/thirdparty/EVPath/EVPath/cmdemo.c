#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "cm.h"
#include "libltdl/ltdl.h"
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

static IOField nested_field_list[] =
{
    {"item", "complex", sizeof(complex), IOOffset(nested_ptr, item)},
    {NULL, NULL, 0, 0}
};

static IOField complex_field_list[] =
{
    {"r", "double", sizeof(double), IOOffset(complex_ptr, r)},
    {"i", "double", sizeof(double), IOOffset(complex_ptr, i)},
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

static IOField simple_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), IOOffset(simple_rec_ptr, integer_field)},
    {"short_field", "integer",
     sizeof(short), IOOffset(simple_rec_ptr, short_field)},
    {"long_field", "integer",
     sizeof(long), IOOffset(simple_rec_ptr, long_field)},
    {"nested_field", "nested",
     sizeof(nested), IOOffset(simple_rec_ptr, nested_field)},
    {"double_field", "float",
     sizeof(double), IOOffset(simple_rec_ptr, double_field)},
    {"char_field", "char",
     sizeof(char), IOOffset(simple_rec_ptr, char_field)},
    {"scan_sum", "integer",
     sizeof(int), IOOffset(simple_rec_ptr, scan_sum)},
    {NULL, NULL, 0, 0}
};

static CMFormatRec simple_format_list[] =
{
    {"complex", complex_field_list},
    {"nested", nested_field_list},
    {NULL, NULL}
};

static
void 
generate_record(event)
simple_rec_ptr event;
{
    long sum = 0;
    event->integer_field = (int) lrand48() % 100;
    sum += event->integer_field % 100;
    event->short_field = ((short) lrand48());
    sum += event->short_field % 100;
    event->long_field = ((long) lrand48());
    sum += event->long_field % 100;

    event->nested_field.item.r = drand48();
    sum += ((int) (event->nested_field.item.r * 100.0)) % 100;
    event->nested_field.item.i = drand48();
    sum += ((int) (event->nested_field.item.i * 100.0)) % 100;

    event->double_field = drand48();
    sum += ((int) (event->double_field * 100.0)) % 100;
    event->char_field = lrand48() % 128;
    sum += event->char_field;
    sum = sum % 100;
    event->scan_sum = (int) sum;
}

int quiet = -1;

static
void
simple_handler(cm, conn, vevent, client_data, attrs)
CManager cm;
CMConnection conn;
void *vevent;
void *client_data;
attr_list attrs;
{
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
    sum += event->integer_field % 100;
    sum += event->short_field % 100;
    sum += event->long_field % 100;
    sum += ((int) (event->nested_field.item.r * 100.0)) % 100;
    sum += ((int) (event->nested_field.item.i * 100.0)) % 100;
    sum += ((int) (event->double_field * 100.0)) % 100;
    sum += event->char_field;
    sum = sum % 100;
    scan_sum = event->scan_sum;
    if (sum != scan_sum) {
	printf("Received record checksum does not match. expected %d, got %d\n",
	       (int) sum, (int) scan_sum);
    }
    if ((quiet <= 0) || (sum != scan_sum)) {
	printf("In the handler, event data is :\n");
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
    }
    printf("Data was received with attributes : \n");
    dump_attr_list(attrs);
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	if (tmp > 0) {
	    *((int *) client_data) = tmp - 1;
	}
    }
}

int
main(argc, argv)
int argc;
char **argv;
{
    CManager cm;
    CMConnection conn = NULL;
    CMFormat format;

    LTDL_SET_PRELOADED_SYMBOLS();
    srand48(getpid());
    cm = CManager_create();
    (void) CMfork_comm_thread(cm);

    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	char *transport = NULL;
	if ((transport = getenv("CMTransport")) != NULL) {
	    listen_list = create_attr_list();
	    add_string_attr(listen_list, CM_TRANSPORT, strdup(transport));
	}
	CMlisten_specific(cm, listen_list);
	contact_list = CMget_contact_list(cm);
	printf("Contact list \"%s\"\n", attr_list_to_string(contact_list));
	format = CMregister_format(cm, "simple", simple_field_list,
				   simple_format_list);
	CMregister_handler(format, simple_handler, NULL);
	CMsleep(cm, 120);
    } else {
	simple_rec data;
	attr_list attrs;
	if (argc == 2) {
	    attr_list contact_list;
	    contact_list = attr_list_from_string(argv[1]);
	    conn = CMinitiate_conn(cm, contact_list);
	    if (conn == NULL) {
		printf("No connection, attr list was :");
		dump_attr_list(contact_list);
		printf("\n");
		exit(1);
	    }
	}
	format = CMregister_format(cm, "simple", simple_field_list,
				   simple_format_list);
	generate_record(&data);
	attrs = create_attr_list();
#define CMDEMO_TEST_ATOM ATL_CHAR_CONS('C','\115','\104','t')
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_int_attr(attrs, CMDEMO_TEST_ATOM, 45678);
	CMwrite_attr(conn, format, &data, attrs);
	free_attr_list(attrs);
    }
    CManager_close(cm);
    return 0;
}
