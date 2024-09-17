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
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)

#else
#include <sys/wait.h>
#endif

#define MSG_COUNT 300

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
    size_t vec_count;
    FFSEncodeVector vecs;
} simple_rec, *simple_rec_ptr;

FMField event_vec_elem_fields[] =
{
    {"len", "integer", sizeof(((FFSEncodeVector)0)[0].iov_len), 
     FMOffset(FFSEncodeVector, iov_len)},
    {"elem", "char[len]", sizeof(char), FMOffset(FFSEncodeVector,iov_base)},
    {(char *) 0, (char *) 0, 0, 0}
};

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
    {"vec_count", "integer",
     sizeof(size_t), FMOffset(simple_rec_ptr, vec_count)},
    {"vecs", "EventVecElem[vec_count]", sizeof(struct FFSEncodeVec), 
     FMOffset(simple_rec_ptr, vecs)},
    {NULL, NULL, 0, 0},
    {NULL, NULL, 0, 0}
};

static FMField done_field_list[] =
{
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {NULL, NULL, 0, 0},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {"complex", complex_field_list, sizeof(complex), NULL},
    {"nested", nested_field_list, sizeof(nested), NULL},
    {"EventVecElem", event_vec_elem_fields, sizeof(struct FFSEncodeVec), NULL},
    {NULL, NULL}
};

static size_t size = 400;
static int vecs = 20;
int quiet = 1;

int message_count = 0;

static
void 
generate_record(simple_rec_ptr event)
{
    int i;
    long sum = 0;
    memset(event, 0, sizeof(*event));
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
    memset(event->vecs, 0, sizeof(event->vecs[0]) * vecs);
    if (quiet <= 0) printf("Sending %d vecs of size %zu\n", vecs, (size_t)size/vecs);
    for (i=0; i < vecs; i++) {
	event->vecs[i].iov_len = size/vecs;
	event->vecs[i].iov_base = malloc(event->vecs[i].iov_len);
	memset(event->vecs[i].iov_base, 0, event->vecs[i].iov_len);
    }
}

static
void 
free_record(simple_rec_ptr event)
{
    int i;
    for (i=0; i < event->vec_count; i++) {
      free(event->vecs[i].iov_base);
    }
    free(event->vecs);
}

static int done=0;
static
void
done_handler(CManager cm, CMConnection conn, void *vevent, void *client_data,
	     attr_list attrs)
{
    done++;
}
static int msg_count = 0;

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
    msg_count++;
#ifdef _MSC_VER
    Sleep(10);
#else
    usleep(10000);
#endif
    if ((quiet <= -1) || (sum != scan_sum)) {
	printf("In the handler, event data is :\n");
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
        printf(" 	Vec_count = %zu\n", event->vec_count);
        if (event->vec_count >= 1) {
            printf("	Vec size = %zu\n", event->vecs[0].iov_len);
        }
	printf("Data was received with attributes : \n");
	dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
	if (quiet <= 0) {
	    printf("message count = %d\n", *((int *) client_data));
	}
	if (msg_count == MSG_COUNT) {
	    int done_msg = 0;
	    CMFormat done_format = CMregister_simple_format(cm, "done", done_field_list, sizeof(int));
	    CMwrite(conn, done_format, &done_msg);
	}
    }
}

static int do_regression_master_test();
static int regression = 1;
static atom_t CM_TRANSPORT;
static char *transport;

#define PARSE_EXTRA_ARGS } else if (strcmp(&argv[1][1], "size") == 0) {\
	    if (sscanf(argv[2], "%zu", &size) != 1) {\
		printf("Unparseable argument to -size, %s\n", argv[2]);\
	    }\
	    if (vecs == 0) { vecs = 1; printf("vecs not 1\n");}\
	    argv++;\
	    argc--;\
	} else 	if (strcmp(&argv[1][1], "vecs") == 0) {\
	    if (sscanf(argv[2], "%d", &vecs) != 1) {\
		printf("Unparseable argument to -vecs, %s\n", argv[2]);\
	    }\
	    argv++;\
	    argc--;

#include "support.c"

int
main(int argc, char **argv)
{
    CManager cm;
    CMConnection conn = NULL;
    CMFormat format, done_format;
    int regression_master = 1;

    PARSE_ARGS();

    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");

    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
    (void) CMfork_comm_thread(cm);

    if (argc == 1) {
	attr_list contact_list, listen_list = NULL;
	if (transport == NULL) {
	    transport = getenv("CMTransport");
	}
	if (transport != NULL) {
	    listen_list = create_attr_list();
	    add_attr(listen_list, CM_TRANSPORT, Attr_String,
		     (attr_value) strdup(transport));
	}
	CMlisten_specific(cm, listen_list);
	contact_list = CMget_contact_list(cm);
	printf("Contact list \"%s\"\n", attr_list_to_string(contact_list));
	format = CMregister_format(cm, simple_format_list);
	CMregister_handler(format, simple_handler, NULL);
	while(msg_count != MSG_COUNT) {
	    CMsleep(cm, 20);
	    printf("Received %d messages\n", msg_count);
	}
    } else {
	simple_rec_ptr data;
	attr_list attrs;
	int i;
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
	data = malloc(sizeof(simple_rec));
	format = CMregister_format(cm, simple_format_list);
	done_format = CMregister_simple_format(cm, "done", done_field_list, sizeof(int));
	CMregister_handler(done_format, done_handler, NULL);
	generate_record(data);
	attrs = create_attr_list();
	CMDEMO_TEST_ATOM = attr_atom_from_string("CMdemo_test_atom");
	add_attr(attrs, CMDEMO_TEST_ATOM, Attr_Int4, (attr_value)45678);
	for (i=0; i < MSG_COUNT; i++) {
	    int block = 0;
	    if (CMConnection_write_would_block(conn)) {
		if (quiet <= 0) printf("Going to block\n");
		block = 1;
	    } else {
		if (quiet <= 0) {
		    printf(".");
		    fflush(stdout);
		}
	    }
	    if (quiet <= 0) {
		printf("Write message count = %d\n", i);
	    }
	    CMwrite_attr(conn, format, data, attrs);
	    if (block) {
		if (quiet <= 0) printf("going again\n");
		block = 0;
	    }
	}
	if (quiet <= 0) printf("Write %d messages\n", MSG_COUNT);
	free_attr_list(attrs);
	free_record(data);
	free(data);
    }
    while(!done) {
	CMsleep(cm, 1);
    }
    if (conn) CMConnection_close(conn);
    CManager_close(cm);
    return 0;
}

static pid_t subproc_proc = 0;
int expected_count = MSG_COUNT;

static void
fail_and_die(int signal)
{
    (void)signal;
    fprintf(stderr, "CMtest failed to complete in reasonable time\n");
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
    char *args[] = {argv0, "-c", NULL, NULL, NULL, NULL, NULL, NULL};
    int exit_state;
    int forked = 0;
    attr_list contact_list, listen_list = NULL;
    char *string_list;
    char size_str[40];
    char vec_str[40];
    CMFormat format;
    int done = 0;
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
    if (transport == NULL) {
	transport = getenv("CMTransport");
    }
    if (transport != NULL) {
	listen_list = create_attr_list();
	add_attr(listen_list, CM_TRANSPORT, Attr_String,
		 (attr_value) strdup(transport));
	printf("TRANSPORT IS %s\n", transport);
    }
    CMlisten_specific(cm, listen_list);
    contact_list = CMget_contact_list(cm);
    string_list = attr_list_to_string(contact_list);
    free_attr_list(contact_list);
    args[2] = "-size";
    sprintf(&size_str[0], "%zu", size);
    args[3] = size_str;
    args[4] = "-vecs";
    sprintf(&vec_str[0], "%d", vecs);
    args[5] = vec_str;
    args[6] = string_list;

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
	CMsleep(cm, 1);
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
		} else {
		    printf("Single remote subproc exit with status %d\n",
			   WEXITSTATUS(exit_state));
		}
	    } else if (WIFSIGNALED(exit_state)) {
		printf("Single remote subproc died with signal %d\n",
		       WTERMSIG(exit_state));
	    }
	    done++;
	}
#endif
    }
    if (msg_count != MSG_COUNT) {
	int i = 10;
	while ((i >= 0) && (msg_count != MSG_COUNT)) {
	    CMsleep(cm, 1);
	}
    }
    free(string_list);
    CManager_close(cm);
    if (message_count != expected_count) {
	printf ("failure, received %d messages instead of %d\n",
		message_count, expected_count);
    }
    return !(message_count == expected_count);
}
