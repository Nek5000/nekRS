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
#include "revpath.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#include < ws2tcpip.h>
#include <winsock2.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0, 0, (DWORD)x),y)
#else
#include <sys/socket.h>
#include <arpa/inet.h>
#include <sys/wait.h>
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

typedef struct _client_rec {
    int output_index;
    int *message_count;
} *client_data_ptr;

static int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    client_data_ptr cdata = (client_data_ptr)client_data;
    int output_index = cdata->output_index;
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
	printf("Received data from output index %d\n", output_index);
	printf("In the handler, event data is :\n");
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (cdata->message_count != NULL) {
	int tmp = cdata->message_count[output_index];
	cdata->message_count[output_index] = tmp + 1;
    }
    return 0;
}

static int do_regression_master_test();

char *ECL_generate = "{\n\
    static int count = 0;\n\
    long sum = 0;\n\
    output.integer_field = (int) lrand48() % 100;\n\
    sum = sum + output.integer_field % 100;\n\
    output.short_field = ((short) lrand48());\n\
    sum = sum + output.short_field % 100;\n\
    output.long_field = ((long) lrand48());\n\
    sum = sum + output.long_field % 100;\n\
\n\
    output.nested_field.item.r = drand48();\n\
    sum = sum + ((int) (output.nested_field.item.r * 100.0)) % 100;\n\
    output.nested_field.item.i = drand48();\n\
    sum = sum + ((int) (output.nested_field.item.i * 100.0)) % 100;\n\
\n\
    output.double_field = drand48();\n\
    sum = sum + ((int) (output.double_field * 100.0)) % 100;\n\
    output.char_field = lrand48() % 128;\n\
    sum = sum + output.char_field;\n\
    sum = sum % 100;\n\
    output.scan_sum = (int) sum;\n\
    count++;\n\
    return count <= 10;\n\
}";

static int dont_fork = 0;
typedef struct {
    char *contact;
} alive_msg_t;
FMField alive_fields[] = {{"contact", "string", sizeof(char*), FMOffset(alive_msg_t*, contact)},
			  {NULL, NULL, 0, 0}};
FMStructDescRec alive_formats[]= {{"alive",alive_fields, sizeof(alive_msg_t), NULL}, {NULL,NULL,0,NULL}};

static atom_t CM_TRANSPORT;
static atom_t CM_NETWORK_POSTFIX;
static atom_t CM_MCAST_ADDR;
static atom_t CM_MCAST_PORT;

static CMConnection
handshake_with_parent(CManager cm, attr_list parent_contact_list)
{
    CMConnection conn = CMinitiate_conn(cm, parent_contact_list);
    alive_msg_t alive;
    attr_list tmp_list;
    CMFormat alive_format;
    alive.contact = attr_list_to_string(tmp_list = CMget_contact_list(cm));
    free_attr_list(tmp_list);
    free_attr_list(parent_contact_list);
    alive_format = CMregister_format(cm, alive_formats);
    CMwrite(conn, alive_format, &alive);
    return conn;
}

static int regression = 1;
static char* transport = NULL;
#include "../mtests/support.c"

int
main(int argc, char **argv)
{
    int regression_master = 1;
    PARSE_ARGS();
    (void)regression;

    srand48(getpid());
    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
    CM_MCAST_PORT = attr_atom_from_string("MCAST_PORT");
    CM_MCAST_ADDR = attr_atom_from_string("MCAST_ADDR");

    if (!regression_master) {
	CManager cm = CManager_create();
	attr_list parent_contact_list = attr_list_from_string(argv[1]);
	CMConnection conn = handshake_with_parent(cm, parent_contact_list);
/*    (void) CMfork_comm_thread(cm);*/
	CMsleep(cm, 20);
	CMConnection_close(conn);
	CManager_close(cm);
	return 0;
    } else {
	return do_regression_master_test();
    }
}

static pid_t subproc_proc = 0;

static void
fail_and_die(int signal)
{
    (void)signal;
    fprintf(stderr, "EVtest failed to complete in reasonable time\n");
    if (subproc_proc != 0) {
	kill(subproc_proc, 9);
    }
    exit(1);
}

static char *router_func = "{\n\
static int count = 0;\n\
count++;\n\
if (count < 6) { return 0;}\n\
if (count < 9) { return 1;}\n\
return 2;\n\
}\0\0";

static FMField filter_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec filter_format_list[] =
{
    {"simple", filter_field_list, sizeof(int), NULL},
    {NULL, NULL}
};

int message_counts[3];
struct _client_rec rec0, rec1, rec2;

static void
alive_handler(CManager cm, CMConnection conn, void *alive_v, 
	      void *client_data, attr_list attrs)
{
    EVaction action, faction;
    EVstone term0, term1, term2, bridge0, bridge1, bridge2, router_stone, source_stone;
    char *filter, *gen_spec;
    int count;
    attr_list self_list;

    message_counts[0] = message_counts[1] = message_counts[2] = 0;
    filter = create_router_action_spec(filter_format_list, router_func);
    gen_spec = create_transform_action_spec(NULL, simple_format_list, ECL_generate);
    

    (void)alive_v; (void)client_data; (void) attrs;

    
/* create local stones */
	term0 = EValloc_stone(cm);
	rec0.output_index = 0;
	rec0.message_count = &message_counts[0];
	EVassoc_terminal_action(cm, term0, simple_format_list, 
				simple_handler, (void*)&rec0);
	term1 = EValloc_stone(cm);
	rec1.output_index = 1;
	rec1.message_count = &message_counts[0];
	EVassoc_terminal_action(cm, term1, simple_format_list, 
				simple_handler, (void*)&rec1);
	term2 = EValloc_stone(cm);
	rec2.output_index = 2;
	rec2.message_count = &message_counts[0];
	EVassoc_terminal_action(cm, term2, simple_format_list, 
				simple_handler, (void*)&rec2);
/* done with local stones */
	
	bridge0 = REValloc_stone(conn);
	self_list = CMget_contact_list(cm);
	REVassoc_bridge_action(conn, bridge0, self_list, term0);
	bridge1 = REValloc_stone(conn);
	REVassoc_bridge_action(conn, bridge1, self_list, term1);
	bridge2 = REValloc_stone(conn);
	REVassoc_bridge_action(conn, bridge2, self_list, term2);
	free_attr_list(self_list);
	REVfreeze_stone(conn, bridge1);
	router_stone = REValloc_stone (conn);
	faction = REVassoc_immediate_action (conn, router_stone, filter);
	REVaction_set_output(conn, router_stone, faction, 0, bridge0);
	REVaction_set_output(conn, router_stone, faction, 1, bridge1);
	REVaction_set_output(conn, router_stone, faction, 2, bridge2);

	source_stone = REValloc_stone (conn);
	action = REVassoc_immediate_action (conn, source_stone, gen_spec);
	free(gen_spec);
	free(filter);
	REVaction_set_output(conn, source_stone, action, 0, router_stone);

	REVenable_auto_stone(conn, source_stone, 0, 100000);
	CMsleep(cm, 4);
	REVfreeze_stone(conn, bridge1);
	count = REVtransfer_events(conn, bridge1, bridge2);
	if (quiet <= 0) {
	    printf("Events transferred : %d\n", count);
	}
	REVdrain_stone(conn, bridge2);
}

static int
do_regression_master_test()
{
    CManager cm;
    char *args[] = {argv0, "-s", NULL, NULL};
    int exit_state;
    int forked = 0;
    int ret, i;
    attr_list contact_list, listen_list = NULL;
    char *string_list, *postfix;
    CMFormat alive_format;
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
    if (!transport) {
	transport = getenv("CMTransport");
    }
    if (transport != NULL) {
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

    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	} else {
	    printf("Doing non-threaded communication handling\n");
	}
    }
    srand48(1);

    args[2] = string_list;
    args[2] = malloc(10 + strlen(string_list));
    sprintf(args[2], "%s", string_list);
    alive_format = CMregister_format(cm, alive_formats);
    CMregister_handler(alive_format, alive_handler, NULL);
    subproc_proc = run_subprocess(args);

    /* give him time to start, run and send us data */
    for (i=0; i< 20; i++) {
	if ((message_counts[0] == 5) &&
	    (message_counts[1] == 0) &&
	    (message_counts[2] == 5)) break;
	CMsleep(cm, 1);
    }
    if (dont_fork) {
	CMsleep(cm, 300);
    } else {
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
    }
    free(string_list);
    free(args[2]);
    CManager_close(cm);
    ret = 0;
    if (message_counts[0] != 5) {
	printf("Message count[0] == %d, should be 5\n", message_counts[0]);
	ret = 1;
    }
    if (message_counts[1] != 0) {
	printf("Message count[1] == %d, should be 0\n", message_counts[1]);
	ret = 1;
    }
    if (message_counts[2] != 5) {
	printf("Message count[2] == %d, should be 5\n", message_counts[2]);
	ret = 1;
    }
    return ret;
}
