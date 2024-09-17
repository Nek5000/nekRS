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
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)
#else
#include <sys/wait.h>
#endif

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
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
    {"thin_message", simple_field_list, sizeof(simple_rec), NULL},
    {NULL, NULL}
};

int quiet = 1;

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
    (void)cm;
    sum += event->integer_field % 100;
    sum += event->short_field % 100;
    sum += event->long_field % 100;
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
extern 	void EVthin_socket_listen(CManager cm,  char **hostname, int *thin_port);

char *transport = NULL;
char *control = NULL;
#include "support.c"

int
main(int argc, char **argv)
{
    CManager cm;
    int regression_master = 1;

    PARSE_ARGS();

    argv0 = "./thin_client";  /* override */
    srand48(getpid());
    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
/*    (void) CMfork_comm_thread(cm);*/

    if (argc == 1) {
	EVstone stone;
	int thin_port = 0;
	char *hostname;
	EVthin_socket_listen(cm,  &hostname, &thin_port);
	stone = EValloc_stone(cm);
	EVassoc_terminal_action(cm, stone, simple_format_list, simple_handler, NULL);
	printf("Thin client contact: %s %d %d\n", hostname, thin_port, stone);
	CMsleep(cm, 120);
    } else {
	printf("Don't run this as a client\n");
	return 1;
    }
    CManager_close(cm);
    return 0;
}

static pid_t subproc_proc = 0;

static void
fail_and_die(int signal)
{
    (void)signal;
    fprintf(stderr, "Thin_test failed to complete in reasonable time\n");
    if (subproc_proc != 0) {
	kill(subproc_proc, 9);
    }
    exit(1);
}

static int
do_regression_master_test()
{
    CManager cm;
    char *args[] = {"thin_client", NULL, NULL, NULL, NULL};
    int exit_state;
    int forked = 0;
    int message_count = 0;
    EVstone handle;
    char *hostname;
    int thin_port;
    int i;
    char port_str[10];
    char stone_str[10];
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
/*    forked = CMfork_comm_thread(cm);*/

    EVthin_socket_listen(cm,  &hostname, &thin_port);

    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	} else {
	    printf("Doing non-threaded communication handling\n");
	}
    }
    srand48(1);

    handle = EValloc_stone(cm);
    EVassoc_terminal_action(cm, handle, simple_format_list, simple_handler, &message_count);
    
    args[1] = hostname;
    args[2] = &port_str[0];
    args[3] = &stone_str[0];
    sprintf(args[2], "%d", thin_port);
    sprintf(args[3], "%d", handle);
    subproc_proc = run_subprocess(args);
    free(hostname);

    /* give him time to start */
    for (i=0; i< 100; i++) {
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
    EVfree_stone(cm, handle);
    CManager_close(cm);
    if (message_count != 1) printf("Message count == %d\n", message_count);
    return !(message_count == 1);
}
