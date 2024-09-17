#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include "evpath.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0, 0, (DWORD)x),y)
#else
#include <sys/socket.h>
#include <sys/wait.h>
#include <arpa/inet.h>
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

int quiet = 1;

static
void
simple_handler(CManager cm, CMConnection conn, void *vevent, void *client_data,
	       attr_list attrs)
{
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
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
	printf("In the handler, connection is %p, event data is :\n", conn);
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
	printf("Data was received with attributes : \n");
	dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
}

static int do_regression_master_test();
static int regression = 1;

static atom_t CM_TRANSPORT;
static atom_t CM_NETWORK_POSTFIX;
static atom_t CM_MCAST_ADDR;
static atom_t CM_MCAST_PORT;

char *transport = NULL;
#include "support.c"

int
main(int argc, char **argv)
{
    CManager cm;
    CMConnection conn = NULL;
    CMFormat format;
    int regression_master = 1;

    PARSE_ARGS();

    srand48(getpid());
    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
    CM_MCAST_PORT = attr_atom_from_string("MCAST_PORT");
    CM_MCAST_ADDR = attr_atom_from_string("MCAST_ADDR");

    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
    (void) CMfork_comm_thread(cm);

    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	char *postfix = NULL;
	char *string_list;
	if ((transport = getenv("CMTransport")) != NULL) {
	    if (listen_list == NULL) listen_list = create_attr_list();
	    add_attr(listen_list, CM_TRANSPORT, Attr_String,
		     (attr_value) strdup(transport));
	}
	if ((postfix = getenv("CMNetworkPostfix")) != NULL) {
	    if (listen_list == NULL) listen_list = create_attr_list();
	    add_attr(listen_list, CM_NETWORK_POSTFIX, Attr_String,
		     (attr_value) strdup(postfix));
	}
	CMlisten(cm);
	CMlisten_specific(cm, listen_list);
	contact_list = CMget_specific_contact_list(cm, listen_list);
	if (transport != NULL) {
	    char *actual_transport = NULL;
	    get_string_attr(contact_list, CM_TRANSPORT, &actual_transport);
	    if (!actual_transport || (strncmp(actual_transport, transport, strlen(actual_transport)) != 0)) {
		printf("Failed to load transport \"%s\"\n", transport);
		exit(1);
	    }
	}
	if (contact_list) {
	    string_list = attr_list_to_string(contact_list);
	} else {
	    /* must be multicast, hardcode a contact list */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
	    int addr;
	    (void) inet_pton(AF_INET, HELLO_GROUP, (struct in_addr *)&addr);
	    contact_list = create_attr_list();
	    add_attr(contact_list, CM_MCAST_ADDR, Attr_Int4,
		     (attr_value) (long)addr);
	    add_attr(contact_list, CM_MCAST_PORT, Attr_Int4,
		     (attr_value) HELLO_PORT);
	    add_attr(contact_list, CM_TRANSPORT, Attr_String,
		     (attr_value) "multicast");
	    conn = CMinitiate_conn(cm, contact_list);
	    string_list = attr_list_to_string(contact_list);
	    free_attr_list(contact_list);
	}	
	printf("Contact list \"%s\"\n", string_list);
	format = CMregister_format(cm, simple_format_list);
	CMregister_handler(format, simple_handler, NULL);
	CMsleep(cm, 120);
    } else {
	simple_rec data;
	attr_list attrs;
	atom_t CMDEMO_TEST_ATOM;
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
	    free_attr_list(contact_list);
	}
	format = CMregister_format(cm, simple_format_list);
	memset(&data, 0, sizeof(data));
	generate_record(&data);
	attrs = create_attr_list();
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_attr(attrs, CMDEMO_TEST_ATOM, Attr_Int4, (attr_value)45678);
	CMwrite_attr(conn, format, &data, attrs);
	CMsleep(cm, 1);
	free_attr_list(attrs);
    }
    if (conn) {
	CMConnection_close(conn);
    }
    CManager_close(cm);
    return 0;
}

static pid_t subproc_proc = 0;

static void
fail_and_die(int signal)
{
    (void)signal;
    fprintf(stderr, "CMtest failed to complete in reasonable time\n");
    if (subproc_proc != 0) {
	kill(subproc_proc, 9);
    }
    exit(1);
}

static int
do_regression_master_test()
{
    CManager cm;
    char *args[] = {argv0, "-c", NULL, NULL};
    int exit_state;
    int forked = 0;
    attr_list contact_list, listen_list = NULL;
    char *string_list, *postfix;
    CMFormat format;
    int message_count = 0, i;
#ifdef HAVE_WINDOWS_H
    SetTimer(NULL, 5, 1000, (TIMERPROC) fail_and_die);
#else
    struct sigaction sigact;
    sigact.sa_flags = 0;
    sigact.sa_handler = fail_and_die;
    sigemptyset(&sigact.sa_mask);
    sigaddset(&sigact.sa_mask, SIGALRM);
    sigaction(SIGALRM, &sigact, NULL);
    alarm(300);
#endif
    cm = CManager_create();
    forked = CMfork_comm_thread(cm);
    if ((transport = getenv("CMTransport")) != NULL) {
	listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String,
		 (attr_value) strdup(transport));
    }
    if ((postfix = getenv("CMNetworkPostfix")) != NULL) {
	if (listen_list == NULL) listen_list = create_attr_list();
	add_attr(listen_list, CM_NETWORK_POSTFIX, Attr_String,
		 (attr_value) strdup(postfix));
    }
    CMlisten_specific(cm, listen_list);
    contact_list = CMget_contact_list(cm);
    if (transport != NULL) {
      char *actual_transport = NULL;
      get_string_attr(contact_list, CM_TRANSPORT, &actual_transport);
      if (!actual_transport || (strncmp(actual_transport, transport, strlen(actual_transport)) != 0)) {
	printf("Failed to load transport \"%s\"\n", transport);
	exit(1);
      }
    }
    if (contact_list) {
	string_list = attr_list_to_string(contact_list);
	free_attr_list(contact_list);
    } else {
	/* must be multicast, hardcode a contact list */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
	int addr;
	(void) inet_pton(AF_INET, HELLO_GROUP, (struct in_addr *)&addr);
	contact_list = create_attr_list();
	add_attr(contact_list, CM_MCAST_ADDR, Attr_Int4,
		 (attr_value) (long)addr);
	add_attr(contact_list, CM_MCAST_PORT, Attr_Int4,
		 (attr_value) HELLO_PORT);
	add_attr(contact_list, CM_TRANSPORT, Attr_String,
		 (attr_value) "multicast");
	(void) CMinitiate_conn(cm, contact_list);
	string_list = attr_list_to_string(contact_list);
	free_attr_list(contact_list);
    }	
    args[2] = string_list;

    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	} else {
	    printf("Doing non-threaded communication handling\n");
	}
    }
    srand48(1);

    format = CMregister_format(cm, simple_format_list);
    CMregister_handler(format, simple_handler, &message_count);
    subproc_proc = run_subprocess(args);

    /* give him time to start */
    for (i=0; i< 10; i++) {
	if (message_count == 1) break;
	CMsleep(cm, 1);
    }
/* stuff */
    if (quiet <= 0) {
	printf("Waiting for remote....\n");
    }
#ifdef HAVE_WINDOWS_H
    if (_cwait(&exit_state, subproc_proc, 0) == -1) {
	perror("cwait");
    }
    if (exit_state == 0) {
	if (quiet <= 0) 
	    printf("Passed single remote subproc test\n");
    } else {
	printf("Single remote subproc exit with status %d\n",
	       exit_state);
    }
#else
    if (waitpid(subproc_proc, &exit_state, 0) == -1) {
	perror("waitpid");
    }
    if (WIFEXITED(exit_state)) {
	if (WEXITSTATUS(exit_state) == 0) {
	    if (quiet <- 1) 
		printf("Passed single remote subproc test\n");
	} else {
	    printf("Single remote subproc exit with status %d\n",
		   WEXITSTATUS(exit_state));
	}
    } else if (WIFSIGNALED(exit_state)) {
	printf("Single remote subproc died with signal %d\n",
	       WTERMSIG(exit_state));
    }
#endif
    free(string_list);
    CManager_close(cm);
    if (message_count != 1) printf("Message count == %d\n", message_count);
    return !(message_count == 1);
}
