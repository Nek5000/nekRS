#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "evpath.h"
#define MAXPROC 8    /* Max number of procsses */
#define NAMELEN 80   /* Max length of machine name */
#define CONTACTLEN 1024
#define LENGTH 24    /* Lengt of send buffer is divisible by 2, 4, 6 and 8 */

int quiet = 1;
static atom_t CM_TRANSPORT;


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

static
void 
generate_record(simple_rec_ptr event)
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

void
simple_handler(CManager cm, CMConnection conn, void *vevent, void *client_data,
	       attr_list attrs)
{
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
    static int message_count = 0;
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);    /* Get nr of processes */
    (void)cm;
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
	printf("In the handler, connection is %lx, event data is :\n", (long)conn);
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
	printf("Data was received with attributes : \n");
	dump_attr_list(attrs);
    }
    message_count++;
    if (message_count == (np - 1)) {
	printf("All %d messages received\n", np - 1);
	CMCondition_signal(cm, (int)(long)client_data);
    }
}

int
main(int argc, char* argv[]) 
{
    int np, me;

    char master_contact[CONTACTLEN];             /* Local host name string */
    char *transport = NULL;
    attr_list listen_list = NULL, contact_list;
    char *string_list;
    CManager cm;

    MPI_Init(&argc, &argv);                /* Initialize MPI */
    MPI_Comm_size(MPI_COMM_WORLD, &np);    /* Get nr of processes */
    MPI_Comm_rank(MPI_COMM_WORLD, &me);    /* Get own identifier */
  
    while (argc > 1) {
	if (strcmp(argv[1], "-t") == 0) {
	    transport = argv[2];
	    argc--;
	    argv++;
	} else if (strcmp(argv[1], "-v") == 0) {
	    quiet--;
	} else {
	    printf("Unknown argument %s\n", argv[1]);
	}
	argv++;
	argc--;
    }
    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    if (transport || (transport = getenv("CMTransport")) != NULL) {
	if (listen_list == NULL) listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String,
		 (attr_value) strdup(transport));
    }
    cm = CManager_create();
    CMlisten_specific(cm, listen_list);
    contact_list = CMget_specific_contact_list(cm, listen_list);
    
    if (transport != NULL) {
      char *actual_transport = NULL;
      get_string_attr(contact_list, CM_TRANSPORT, &actual_transport);
      if (!actual_transport || (strcmp(actual_transport, transport) != 0)) {
	printf("Failed to load transport \"%s\"\n", transport);
	MPI_Finalize();
	exit(1);
      }
    }
    string_list = attr_list_to_string(contact_list);
    free_attr_list(contact_list);
    strcpy(master_contact, string_list);

    if (me == 0) {    /* Process 0 does this */
	CMFormat format;
	int message_wait_condition;
	time_t start, end;
	if (quiet <= 0) {
	    printf("Master contact is %s\n", master_contact);
	}
	format = CMregister_format(cm, simple_format_list);
	message_wait_condition = CMCondition_get(cm, NULL);
	CMregister_handler(format, simple_handler, (void*)(long)message_wait_condition);
	MPI_Bcast(master_contact,CONTACTLEN,MPI_CHAR,0,MPI_COMM_WORLD);
 	start = time(NULL);
	CMCondition_wait(cm, message_wait_condition);
	end = time(NULL);
	printf("Elapsed time was %d\n", (int)(end - start));
    } else { /* all other processes do this */
	attr_list contact_list;
	CMConnection conn = NULL;
	CMFormat format;
	simple_rec data;
	attr_list attrs;
	atom_t CMDEMO_TEST_ATOM;
	MPI_Bcast(master_contact,CONTACTLEN,MPI_CHAR,0,MPI_COMM_WORLD);
	if (quiet <= 0) {
	    printf("Node %d thinks master contact is %s\n", me, master_contact);
	}
	contact_list = attr_list_from_string(master_contact);
	conn = CMinitiate_conn(cm, contact_list);
	if (conn == NULL) {
	    printf("No connection, attr list was :");
	    dump_attr_list(contact_list);
	    printf("\n");
	    exit(1);
	}
	free_attr_list(contact_list);
	format = CMregister_format(cm, simple_format_list);
	generate_record(&data);
	attrs = create_attr_list();
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_attr(attrs, CMDEMO_TEST_ATOM, Attr_Int4, (attr_value)45678);
	CMwrite_attr(conn, format, &data, attrs);
	CMsleep(cm, 1);
	free_attr_list(attrs);
    }

    MPI_Finalize();
    exit(0);
}
