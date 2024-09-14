#if defined (__INTEL_COMPILER)
#pragma warning (disable:981)
#endif
#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include "evpath.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)
#else
#include <sys/wait.h>
#include <arpa/inet.h>
#endif

typedef struct _rec_a {
    int a_field;
    int sequence;
} rec_a, *rec_a_ptr;

typedef struct _rec_b {
    int b_field;
    int sequence;
} rec_b, *rec_b_ptr;

typedef struct _rec_c {
    int c_field;
    int sequence_a;
    int sequence_b;
} rec_c, *rec_c_ptr;

static FMField a_field_list[] =
{
    {"a_field", "integer",
     sizeof(int), FMOffset(rec_a_ptr, a_field)},
    {"sequence", "integer",
     sizeof(int), FMOffset(rec_a_ptr, sequence)},
    {NULL, NULL, 0, 0}
};

static FMField b_field_list[] =
{
    {"b_field", "integer",
     sizeof(int), FMOffset(rec_b_ptr, b_field)},
    {"sequence", "integer",
     sizeof(int), FMOffset(rec_b_ptr, sequence)},
    {NULL, NULL, 0, 0}
};

static FMField c_field_list[] =
{
    {"c_field", "integer",
     sizeof(int), FMOffset(rec_c_ptr, c_field)},
    {"sequence_a", "integer",
     sizeof(int), FMOffset(rec_c_ptr, sequence_a)},
    {"sequence_b", "integer",
     sizeof(int), FMOffset(rec_c_ptr, sequence_b)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec a_format_list[] =
{
    {"a_rec", a_field_list, sizeof(rec_a), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec b_format_list[] =
{
    {"b_rec", b_field_list, sizeof(rec_b), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec c_format_list[] =
{
    {"c_rec", c_field_list, sizeof(rec_c), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescList queue_list[] = {a_format_list, b_format_list, c_format_list, NULL};


static
void
generate_a_record(rec_a_ptr event)
{
    static int sequence = 0;
    /* always even */
    event->a_field = ((int) lrand48() % 50) * 2;
    event->sequence = sequence++;
}

static
void
generate_b_record(rec_b_ptr event)
{
    static int sequence = 0;
    /* always odd */
    event->b_field = ((int) lrand48() % 50) * 2 + 1;
    event->sequence = sequence++;
}

int quiet = 1;
static int repeat_count = 100;
char *a_map = NULL;
char *b_map = NULL;
static
int
output_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    rec_c_ptr event = vevent;
    (void)cm;
    if (event->c_field % 2 != 1) {
	printf("Received record should be odd, got %d\n", event->c_field);
    }

    if (a_map == NULL) {
	a_map = malloc(repeat_count);
	memset(a_map, 0, repeat_count);
	b_map = malloc(repeat_count);
	memset(b_map, 0, repeat_count);
    }
    a_map[event->sequence_a]++;
    b_map[event->sequence_b]++;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	c_field = %d\n", event->c_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
    return 0;
}

static int do_regression_master_test();
static int regression = 1;
static atom_t CM_TRANSPORT;
static atom_t CM_NETWORK_POSTFIX;
static atom_t CM_MCAST_ADDR;
static atom_t CM_MCAST_PORT;

static char *trans = "{\n\
    int found = 0;\n\
    a_rec *a;\n\
    b_rec *b;\n\
    c_rec c;\n\
    if (EVpresent(a_rec_ID, 0)) {\n\
	int seq = EVdata_a_rec(0).sequence;\n\
        a = EVdata_a_rec(0); ++found;\n\
    }\n\
    if (EVpresent(b_rec_ID, 0)) {\n\
        b = EVdata_b_rec(0); ++found;\n\
    }\n\
    if (found == 2) {\n\
	attr_list l = create_attr_list();\n\
	set_int_attr(l, \"CMdemo_test_atom\", 56789);\n\
        c.c_field = a.a_field + b.b_field;\n\
	c.sequence_a = a.sequence;\n\
	c.sequence_b = b.sequence;\n\
        if (!EVpresent_b_rec(0))\n\
            printf(\"??? <1> not present (1)\\n\");\n\
        EVdiscard_a_rec(0);\n\
        if (!EVpresent_b_rec(0))\n\
            printf(\"??? <2> not present (1)\\n\");\n\
        EVdiscard_b_rec(0);\n\
        EVsubmit_attr(0, c, l);\n\
        free_attr_list(l);\n\
    }\n\
}\0\0";


static void
data_free(void *event_data, void *client_data)
{
    (void) client_data;
    free(event_data);
}

char *transport = NULL;
char *control = NULL;
#include "support.c"

int
main(int argc, char **argv)
{
    CManager cm;
    int regression_master = 1;

    /* XXX for testing */ setbuf(stdout, NULL);

    PARSE_ARGS();

    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
    CM_MCAST_PORT = attr_atom_from_string("MCAST_PORT");
    CM_MCAST_ADDR = attr_atom_from_string("MCAST_ADDR");

    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
/*    (void) CMfork_comm_thread(cm);*/

    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	char *transport = NULL;
	char *postfix = NULL;
	char *string_list;
	char *filter;
	EVstone term, fstone;
	EVaction faction;
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
	CMlisten_specific(cm, listen_list);
	contact_list = CMget_contact_list(cm);
	if (contact_list) {
	    string_list = attr_list_to_string(contact_list);
	} else {
	    /* must be multicast, hardcode a contact list */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
	    int addr;
	    (void) inet_aton(HELLO_GROUP, (struct in_addr *)&addr);
	    contact_list = create_attr_list();
	    add_attr(contact_list, CM_MCAST_ADDR, Attr_Int4,
		     (attr_value) (long)addr);
	    add_attr(contact_list, CM_MCAST_PORT, Attr_Int4,
		     (attr_value) HELLO_PORT);
	    add_attr(contact_list, CM_TRANSPORT, Attr_String,
		     (attr_value) "multicast");
/*	    conn = CMinitiate_conn(cm, contact_list);*/
	    string_list = attr_list_to_string(contact_list);
	    free_attr_list(contact_list);
	}
	term = EValloc_stone(cm);
	EVassoc_terminal_action(cm, term, c_format_list, output_handler, NULL);
	EVassoc_terminal_action(cm, term, a_format_list, output_handler, NULL);
	filter = create_multityped_action_spec(queue_list, trans);

	fstone = EValloc_stone(cm);
	faction = EVassoc_multi_action(cm, fstone, filter, NULL);
	EVaction_set_output(cm, fstone, faction, 0, term);

	printf("Contact list \"%d:%s\"\n", fstone, string_list);
	CMsleep(cm, 120);
	free(filter);
    } else {
	attr_list attrs;
	int remote_stone, stone = 0;
	int count, i;
	char *map;
	int a_count = 0;
	int b_count = 0;
	EVsource a_handle, b_handle;
	atom_t CMDEMO_TEST_ATOM;
	if (argc == 2) {
	    attr_list contact_list;
	    char *list_str;
	    sscanf(argv[1], "%d:", &remote_stone);
	    list_str = strchr(argv[1], ':') + 1;
	    contact_list = attr_list_from_string(list_str);
	    stone = EValloc_stone(cm);
	    EVassoc_bridge_action(cm, stone, contact_list, remote_stone);
	    free_attr_list(contact_list);
	}
	attrs = create_attr_list();
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_attr(attrs, CMDEMO_TEST_ATOM, Attr_Int4, (attr_value)45678);
	a_handle = EVcreate_submit_handle_free(cm, stone, a_format_list,
					       data_free, NULL);
	b_handle = EVcreate_submit_handle_free(cm, stone, b_format_list,
					       data_free, NULL);
	count = repeat_count;
	map = malloc(count);
	memset(map, 0, count);
	/* setup map so that it is half ones and half zeroes */
	srand48(time(NULL));
	for (i=0; i < count / 2 ; i++) {
	    int try = lrand48() % count;
	    if (map[try] == 0) {
		map[try] = 1;
	    } else {
		i--;  /* try again */
	    }
	}
#ifdef NOT_DEF
	printf("Map : ");
	for (i=0; i < count ; i++) {
	    printf("%1d", map[i]);
	    if (map[i] == 1) {
		a_count++;
	    } else {
		b_count++;
	    }
	}
	printf("\n");
#endif
	if (a_count != b_count) printf("MAP FUNCTION FAILED TO GENERATE EQUAL MAP\n");
	for (i=0; i < count ; i++) {
	    if (map[i] == 1) {
		rec_a_ptr a = malloc(sizeof(*a));
		memset(a, 0, sizeof(*a));
		generate_a_record(a);
		if (quiet <=0) {printf("submitting a -> %d\n", a->a_field);}
		EVsubmit(a_handle, a, attrs);
	    } else {
		rec_b_ptr b = malloc(sizeof(*b));
		memset(b, 0, sizeof(*b));
		generate_b_record(b);
		if (quiet <=0) {printf("submitting b -> %d\n", b->b_field);}
		EVsubmit(b_handle, b, attrs);
	    }
	}
	EVfree_source(a_handle);
	EVfree_source(b_handle);
	CMsleep(cm, 1);
	free_attr_list(attrs);
	free(map);
    }
    CManager_close(cm);
    return 0;
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

static int
do_regression_master_test()
{
    CManager cm;
    char *args[] = {"multiq_test", "-c", NULL, NULL};
    char *filter;
    int exit_state;
    int forked = 0;
    attr_list contact_list, listen_list = NULL;
    char *string_list, *transport, *postfix;
    int message_count = 0;
    EVstone term, fstone;
    EVaction faction;
    int i;
#ifdef HAVE_WINDOWS_H
    SetTimer(NULL, 5, 300*1000, (TIMERPROC) fail_and_die);
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
    if (contact_list) {
	string_list = attr_list_to_string(contact_list);
	free_attr_list(contact_list);
    } else {
	/* must be multicast, hardcode a contact list */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
	int addr;
	(void) inet_aton(HELLO_GROUP, (struct in_addr *)&addr);
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

    term = EValloc_stone(cm);
    EVassoc_terminal_action(cm, term, c_format_list, output_handler, &message_count);

    filter = create_multityped_action_spec(queue_list, trans);

    fstone = EValloc_stone(cm);
    faction = EVassoc_multi_action(cm, fstone, filter, NULL);
    EVaction_set_output(cm, fstone, faction, 0, term);

    args[2] = string_list;
    args[2] = malloc(10 + strlen(string_list) + strlen(filter));
    sprintf(args[2], "%d:%s", fstone, string_list);
    subproc_proc = run_subprocess(args);
    free(filter);
    free(args[2]);

    /* give him time to start */
    for (i=0; i< 10; i++) {
	if (message_count == repeat_count / 2) break;
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
    if (message_count != repeat_count / 2) printf("Message count == %d\n", message_count);
    if (a_map == NULL) return !(message_count == repeat_count / 2);
    for (i=0; i < repeat_count / 2 ; i++) {
	if (a_map[i] != 1) printf("A_map[%d] = %d\n", i, a_map[i]);
	if (b_map[i] != 1) printf("B_map[%d] = %d\n", i, b_map[i]);
    }
    free(a_map);
    free(b_map);
    return !(message_count == repeat_count / 2);
}
