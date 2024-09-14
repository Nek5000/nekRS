#include "../config.h"

#include "chr_time.h"

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
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)
#else
#include <sys/wait.h>
#include <arpa/inet.h>
#endif

#define MSG_COUNT 30
static int msg_limit = MSG_COUNT;
static int message_count = 0;
static int expected_count = MSG_COUNT;

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

typedef struct response_rec {
    int condition;
} *response_rec_ptr;

typedef struct _simple_rec {
    int condition;
    char *contact_list;
    int target_stone;
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
    int vec_count;
    FFSEncodeVector vecs;
} simple_rec, *simple_rec_ptr;

FMField event_vec_elem_fields[] =
{
    {"len", "integer", sizeof(((FFSEncodeVector)0)[0].iov_len), 
     FMOffset(FFSEncodeVector, iov_len)},
    {"elem", "char[len]", sizeof(char), FMOffset(FFSEncodeVector,iov_base)},
    {(char *) 0, (char *) 0, 0, 0}
};

static FMField response_field_list[] =
{
    {"condition", "integer",
     sizeof(int), FMOffset(response_rec_ptr, condition)},
    {(char *) 0, (char *) 0, 0, 0}
};

static FMField simple_field_list[] =
{
    {"condition", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, condition)},
    {"contact_list", "string",
     sizeof(char*), FMOffset(simple_rec_ptr, contact_list)},
    {"target_stone", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, target_stone)},
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
    {"vec_count", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, vec_count)},
    {"vecs", "EventVecElem[vec_count]", sizeof(struct FFSEncodeVec), 
     FMOffset(simple_rec_ptr, vecs)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {"complex", complex_field_list, sizeof(complex), NULL},
    {"nested", nested_field_list, sizeof(nested), NULL},
    {"EventVecElem", event_vec_elem_fields, sizeof(struct FFSEncodeVec), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescRec response_format_list[] =
{
    {"response", response_field_list, sizeof(*((response_rec_ptr) NULL)), NULL},
    {NULL, NULL, 0, NULL}
};

static int size = 400;
static int vecs = 20;
static int quiet = 1;
static int request_response = 0;
static int print_bandwidth = 0;

static
void 
generate_record(simple_rec_ptr event)
{
    int i;
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
    event->vec_count = vecs;
    event->vecs = malloc(sizeof(event->vecs[0]) * vecs);
    if (quiet <= 0) printf("Sending %d vecs of size %d\n", vecs, size/vecs);
    for (i=0; i < vecs; i++) {
	event->vecs[i].iov_len = size/vecs;
	event->vecs[i].iov_base = malloc(event->vecs[i].iov_len);
    }
}

static int msg_count = 0;

static
int
response_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    CMCondition_signal(cm, ((response_rec_ptr)vevent)->condition);
    return 0;
}

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    static chr_time bandwidth_start_time;
    static EVsource response_handle;
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
    (void) cm;
    sum += event->integer_field % 100;
    sum += event->short_field % 100;
    sum += event->long_field % 100;
    sum += ((int) (event->nested_field.item.r * 100.0)) % 100;
    sum += ((int) (event->nested_field.item.i * 100.0)) % 100;
    sum += ((int) (event->double_field * 100.0)) % 100;
    sum += event->char_field;
    sum = sum % 100;
    scan_sum = event->scan_sum;
    if (msg_count == 0) {
	if (event->condition != -1) {
	    /* request_response! */
	    EVstone response_stone = EValloc_stone(cm);
	    attr_list contact_list = attr_list_from_string(event->contact_list);
	    EVassoc_bridge_action(cm, response_stone, contact_list, event->target_stone);
	    free_attr_list(contact_list);
	    response_handle = EVcreate_submit_handle(cm, response_stone, response_format_list);
	}
	chr_timer_start(&bandwidth_start_time);
    }
    if (sum != scan_sum) {
	printf("Received record checksum does not match. expected %d, got %d\n",
	       (int) sum, (int) scan_sum);
    }
    msg_count++;
//    usleep(10000);
    if ( quiet <= 0) {
      printf(".\n");
    }
    if ((quiet <= -1) || (sum != scan_sum)) {
	printf("In the handler, event data is :\n");
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
    if ((msg_count == msg_limit) && print_bandwidth) {
	chr_timer_stop(&bandwidth_start_time);
	double secs = chr_time_to_secs(&bandwidth_start_time);
	long data_size = (msg_limit-1) * size;
	double megabits = (double)data_size*8 / ((double)1000*1000);
	double megabits_sec = megabits / secs;
	printf("Megabits/sec is %g\n", megabits_sec);
	printf("transport = %s size = %d, count = %d, secs = %g, Mbps = %g\n",
		 "default", size, msg_count, secs, megabits_sec);
    }
    if (event->condition != -1) {
	struct response_rec response;
	response.condition = event->condition;
	EVsubmit(response_handle, &response, NULL);
    }
    return 0;
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
    int regression_master = 1;
    int forked = 0;

    argv0 = argv[0];\
    while (argv[1] && (argv[1][0] == '-')) {
	if (strcmp(&argv[1][1], "size") == 0) {
	    if (sscanf(argv[2], "%d", &size) != 1) {
		printf("Unparseable argument to -size, %s\n", argv[2]);
	    }
	    if (vecs == 0) { vecs = 1; printf("vecs not 1\n");}
	    argv++;
	    argc--;
	} else 	if (strcmp(&argv[1][1], "vecs") == 0) {
	    if (sscanf(argv[2], "%d", &vecs) != 1) {
		printf("Unparseable argument to -vecs, %s\n", argv[2]);
	    }
	    argv++;
	    argc--;
	} else 	if (strcmp(&argv[1][1], "msgs") == 0) {
	    if (sscanf(argv[2], "%d", &msg_limit) != 1) {
		printf("Unparseable argument to -msgs, %s\n", argv[2]);
	    }
	    expected_count = msg_limit;
	    argv++;
	    argc--;
	} else if (strcmp(&argv[1][1], "ssh") == 0) {
	    char *destination_host;
	    char *first_colon, *second_colon;
	    char *ssh_port = NULL;
	    if (!argv[2]) {
	        printf("Missing --ssh destination\n");
		usage();
	    }
	    first_colon = strchr(argv[2], ':');
	    if (first_colon) {
	        *first_colon = 0;
		second_colon = strchr(first_colon+1, ':');
	    } else {
	        second_colon = NULL;
	    }
	    destination_host = strdup(argv[2]);
	    if (first_colon) {
	        int ssh_port_int;
		if (second_colon) *second_colon = 0;
		if (sscanf(first_colon+1, "%d", &ssh_port_int) != 1) {
		    second_colon = first_colon;
		}  else {
		    ssh_port = first_colon + 1;
		}
	    }
	    if (second_colon) {
	        strcpy(remote_directory, second_colon+1);
	    }
	    if (strlen(SSH_PATH) == 0) {
		printf("SSH_PATH in config.h is empty!  Can't run ssh\n");
		exit(1);
	    }
	    ssh_args[0] = strdup(SSH_PATH);
	    ssh_args[1] = destination_host;
	    ssh_args[2] = NULL;
	    if (ssh_port != NULL) {
	        ssh_args[2] = "-p";
	        ssh_args[3] = ssh_port;
		ssh_args[4] = NULL;
	    }
	    argv++; argc--;
	} else if (argv[1][1] == 'c') {
	    regression_master = 0;
	} else if (argv[1][1] == 's') {
	    regression_master = 0;
	} else if (argv[1][1] == 'q') {
	    quiet++;
	} else if (argv[1][1] == 'v') {
	    quiet--;
	} else if (argv[1][1] == 'n') {
	    regression = 0;
	    quiet = -1;
	} else if (argv[1][1] == 'r') {
	    request_response = 1;
	} else if (argv[1][1] == 'b') {
	    print_bandwidth = 1;
	} else if (argv[1][1] == 't') {
	    transport = argv[2];
	    argv++;
	    argc--;
	}
	argv++;
	argc--;
    }
    srand48(getpid());
    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");
    CM_NETWORK_POSTFIX = attr_atom_from_string("CM_NETWORK_POSTFIX");
    CM_MCAST_PORT = attr_atom_from_string("MCAST_PORT");
    CM_MCAST_ADDR = attr_atom_from_string("MCAST_ADDR");

    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
    forked = CMfork_comm_thread(cm);
    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	} else {
	    printf("Doing non-threaded communication handling\n");
	}
    }

    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	char *postfix = NULL;
	char *string_list;
	EVstone stone;
	if (transport || ((transport = getenv("CMTransport")) != NULL)) {
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
	if (listen_list) free_attr_list(listen_list);
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
	stone = EValloc_stone(cm);
	EVassoc_terminal_action(cm, stone, simple_format_list, simple_handler, NULL);
	printf("Contact list \"%d:%s\"\n", stone, string_list);
	while(msg_count != msg_limit) {
	    CMsleep(cm, 20);
	    printf("Received %d messages\n", msg_count);
	}
    } else {
	simple_rec_ptr data;
	attr_list attrs;
	int i;
	int remote_stone, stone = 0;
	EVsource source_handle;
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
	data = malloc(sizeof(simple_rec));
	generate_record(data);
	attrs = create_attr_list();
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_attr(attrs, CMDEMO_TEST_ATOM, Attr_Int4, (attr_value)45678);
	source_handle = EVcreate_submit_handle(cm, stone, simple_format_list);
	for (i=0; i < msg_limit; i++) {
	    data->integer_field++;
	    data->long_field--;
	    data->contact_list = NULL;
	    data->target_stone = -1;
	    if (quiet <= 0) printf("Submitting %d of %d\n", i, msg_limit);
	    if (request_response) {
		data->condition = CMCondition_get(cm, NULL);
		if (i == 0) {
		    CMlisten(cm);
		    attr_list contact_list = CMget_contact_list(cm);
		    EVstone stone = EValloc_stone(cm);
		    data->target_stone = stone;
		    data->contact_list = attr_list_to_string(contact_list);
		    EVassoc_terminal_action(cm, stone, response_format_list, response_handler, NULL);
		}
	    } else {
		data->condition = -1;
	    }
	    EVsubmit(source_handle, data, attrs);
	    if (request_response) CMCondition_wait(cm, data->condition);
	}
	if (quiet <= 0) printf("Wrote %d messages\n", msg_limit);
	CMsleep(cm, 30);
	free_attr_list(attrs);
	for (i=0; i < vecs; i++) {
	    free(data->vecs[i].iov_base);
	}
	free(data->vecs);
	free(data);
	EVfree_source(source_handle);
    }
    CManager_close(cm);
    return 0;
}

static pid_t subproc_proc = 0;

static void
fail_and_die(int signal)
{
    (void)signal;
    fprintf(stderr, "bulktest failed to complete in reasonable time\n");
    if (message_count != expected_count) {
	printf ("failure, received %d messages instead of %d\n",
		message_count, expected_count);
    }
    if (subproc_proc != 0) {
	kill(subproc_proc, 9);
    }
    exit(1);
}

static int
do_regression_master_test()
{
    CManager cm;
    char* args[20] = { "bulktest", "-c", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
    int last_arg = 2;
    int exit_state;
    int forked = 0;
    attr_list contact_list, listen_list = NULL;
    char* string_list, * postfix;
    char size_str[40];
    char vec_str[40];
    char msg_str[40];
    EVstone handle;
    int done = 0;
#ifdef HAVE_WINDOWS_H
    SetTimer(NULL, 5, 300*1000, (TIMERPROC)fail_and_die);
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
    if (transport || ((transport = getenv("CMTransport")) != NULL)) {
	listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String,
	    (attr_value)strdup(transport));
    }
    if ((postfix = getenv("CMNetworkPostfix")) != NULL) {
	if (listen_list == NULL) listen_list = create_attr_list();
	add_attr(listen_list, CM_NETWORK_POSTFIX, Attr_String,
	    (attr_value)strdup(postfix));
    }
    CMlisten_specific(cm, listen_list);
    if (listen_list) free_attr_list(listen_list);
    contact_list = CMget_contact_list(cm);
    if (transport != NULL) {
	char* actual_transport = NULL;
	get_string_attr(contact_list, CM_TRANSPORT, &actual_transport);
	if (!actual_transport || (strncmp(actual_transport, transport, strlen(actual_transport)) != 0)) {
	    printf("Failed to load transport \"%s\"\n", transport);
	    exit(1);
	}
    }

    if (contact_list) {
	string_list = attr_list_to_string(contact_list);
	free_attr_list(contact_list);
    }
    else {
	/* must be multicast, hardcode a contact list */
#define HELLO_PORT 12345
#define HELLO_GROUP "225.0.0.37"
	int addr;
	(void)inet_aton(HELLO_GROUP, (struct in_addr*)&addr);
	contact_list = create_attr_list();
	add_attr(contact_list, CM_MCAST_ADDR, Attr_Int4,
	    (attr_value)(long)addr);
	add_attr(contact_list, CM_MCAST_PORT, Attr_Int4,
	    (attr_value)HELLO_PORT);
	add_attr(contact_list, CM_TRANSPORT, Attr_String,
	    (attr_value)"multicast");
	(void)CMinitiate_conn(cm, contact_list);
	string_list = attr_list_to_string(contact_list);
	free_attr_list(contact_list);
    }
    args[last_arg++] = "-size";
    sprintf(&size_str[0], "%d", size);
    args[last_arg++] = size_str;
    args[last_arg++] = "-vecs";
    sprintf(&vec_str[0], "%d", vecs);
    args[last_arg++] = vec_str;
    args[last_arg++] = "-msgs";
    sprintf(&msg_str[0], "%d", msg_limit);
    args[last_arg++] = msg_str;
    if (request_response) args[last_arg++] = "-r";
    args[last_arg] = malloc(strlen(string_list) + 10);

    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	}
	else {
	    printf("Doing non-threaded communication handling\n");
	}
    }
    srand48(1);

    handle = EValloc_stone(cm);
    EVassoc_terminal_action(cm, handle, simple_format_list, simple_handler, &message_count);
    sprintf(args[last_arg], "%d:%s", handle, string_list);
    subproc_proc = run_subprocess(args);

    if (quiet <= 0) {
	printf("Waiting for remote....\n");
    }
    while (!done) {
#ifdef HAVE_WINDOWS_H
	if (_cwait(&exit_state, subproc_proc, 0) == -1) {
	    perror("cwait");
	}
	if (exit_state == 0) {
	    if (quiet <= 0)
		printf("Subproc exitted\n");
	} else {
	    printf("Single remote subproc exit with status %d\n",
		exit_state);
	}
#else
	int result;
	if (quiet <= 0) {
	    printf(",");
	    fflush(stdout);
	}
	CMsleep(cm, 50);	done++;

	result = waitpid(subproc_proc, &exit_state, WNOHANG);
	if (result == -1) {
	    perror("waitpid");
	    done++;
	}
	if (result == subproc_proc) {
	    if (WIFEXITED(exit_state)) {
		if (WEXITSTATUS(exit_state) == 0) {
		    if (quiet <= 0)
			printf("Subproc exited\n");
		}
		else {
		    printf("Single remote subproc exit with status %d\n",
			WEXITSTATUS(exit_state));
		}
	    }
	    else if (WIFSIGNALED(exit_state)) {
		printf("Single remote subproc died with signal %d\n",
		    WTERMSIG(exit_state));
	    }
	    done++;
	}
#endif
    }
    if (msg_count != msg_limit) {
	int i = 10;
	while ((i >= 0) && (msg_count != msg_limit)) {
	    CMsleep(cm, 1);
	}
    }
    free(args[last_arg]);
    free(string_list);
    CManager_close(cm);
    if (message_count != expected_count) {
	printf("failure, received %d messages instead of %d\n",
	    message_count, expected_count);
    }
    return !(message_count == expected_count);
}
