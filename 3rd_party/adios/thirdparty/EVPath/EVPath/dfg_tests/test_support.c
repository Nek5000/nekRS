#include "../config.h"
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
#include "ev_dfg.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)
#define waitpid(x, tmp, z) {   WaitForSingleObject(OpenProcess(0,0,(DWORD)x), INFINITE);    GetExitCodeProcess(OpenProcess(0,0,(DWORD)x), tmp);}
#else
#include <sys/wait.h>
#endif

#include "test_support.h"

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

FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {"complex", complex_field_list, sizeof(complex), NULL},
    {"nested", nested_field_list, sizeof(nested), NULL},
    {NULL, NULL}
};

void 
generate_simple_record(simple_rec_ptr event)
{
    long sum = 0;
    memset(event, 0, sizeof(simple_rec));
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

extern int
checksum_simple_record(simple_rec_ptr event, attr_list attrs, int quiet)
{
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
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    return (sum == scan_sum);
}

static pid_t subproc_proc = 0;
int quiet = 1;
static int no_fork = 0;
void(*on_exit_handler)() = NULL;
static int doing_ssh = 0;

static void
fail_and_die(int signal)
{
    (void) signal;
    fprintf(stderr, "EVPath test failed to complete in reasonable time\n");
    if (on_exit_handler) on_exit_handler();
    if (subproc_proc != 0) {
	kill(subproc_proc, 9);
    }
    exit(1);
}

static
pid_t
run_subprocess(char **args)
{
    static int count = 0;
#ifdef HAVE_WINDOWS_H
    intptr_t child;
    child = _spawnv(_P_NOWAIT, "./evtest.exe", args);
    if (child == -1) {
	printf("failed for evtest\n");
	perror("spawnv");
    }
    return child;
#else
    pid_t child = -1;
    if (quiet <=0) {printf("Forking subprocess\n");}
    if (no_fork) {
	int i = 0;
	printf("Would have run :");
	while(args[i] != NULL) printf(" %s", args[i++]);
	printf("\n");
    } else {
	child = fork();
    }

    if (child == 0) {
	/* I'm the child */
	execv(args[0], args);
    }
    if((count++ > 0) && doing_ssh) usleep(500000);  /* don't go too fast */
    return child;
#endif
}

static char *argv0;
static pid_t *pid_list = NULL;
static int start_subproc_arg_count;
static char **subproc_args;
static int cur_subproc_arg;
static char remote_directory[1024];

extern void
test_fork_children(char **list, char *master_contact)
{
    char *args[20];
    int i = 0;
    int list_index = 1;
    static int pid_index = 1;
    /* assume that we are list[0] */
    while(subproc_args[start_subproc_arg_count + i] != NULL) {
        args[i] = subproc_args[start_subproc_arg_count+i];
	i++;
    }
    if (quiet < 1) {
	args[i++] = "-v";
    }	
    args[i++] = "-c";
    args[i+1] = master_contact;
    args[i+2] = NULL;
    if (!pid_list) pid_list = malloc(sizeof(pid_list[0]));
    while(list[list_index] != NULL) {
	args[i] = list[list_index];
	pid_list[pid_index -1 ] = run_subprocess(args);
	list_index++;
	pid_index++;
	pid_list = realloc(pid_list, sizeof(pid_list[0]) * pid_index);
    }
    pid_list[pid_index - 1] = 0;
}

static void
delay_fork_wrapper(CManager cm, void *client_data)
{
    delay_struct *str = (delay_struct*)client_data;
    (void) cm;
    test_fork_children(str->list, str->master_contact);
    free(str);
}

extern void
delayed_fork_children(CManager cm, char **list, char *master_contact, int delay_seconds)
{
    delay_struct *str = malloc(sizeof(delay_struct));
    str->list = list;
    str->master_contact = master_contact;
    CMTaskHandle handle = CMadd_delayed_task(cm, delay_seconds, 0, delay_fork_wrapper, (void*) str);
    free(handle);
}

int wait_for_children(char **list)
{
    int i=0, stat;
    (void)list;
    while(pid_list && pid_list[i] != 0) {
	waitpid(pid_list[i], &stat, 0);
	i++;
    }
    free(pid_list);
    pid_list = NULL;
    return 0;
}

static int regression = 1;

static void fail_and_die(int signal);

static void
usage()
{
    printf("Usage:  %s <options> \n", argv0);
    printf("  Options:\n");
    printf("\t-q  quiet\n");
    printf("\t-v  verbose\n");
    printf("\t-n  No regression test.  I.E. just run the master and print \n\t\twhat command would have been run for client.\n");
    printf("\t-ssh <hostname>:<ssh_port>:<remote directory>  parameters to use for remote client via ssh.\n");
    printf("\t-ssh <hostname>:<remote directory>  parameters to use for remote client via ssh.\n");

    exit(1);
}

int
main(int argc, char **argv)
{
    int regression_master = 1;
    int ret = 0;

    argv0 = argv[0];
    start_subproc_arg_count = 8; /* leave a few open at the beginning */
    subproc_args = calloc((argc + start_subproc_arg_count + 8), sizeof(argv[0]));
    cur_subproc_arg = start_subproc_arg_count;
    remote_directory[0] = 0;

    while (argv[1] && (argv[1][0] == '-')) {
	if (argv[1][1] == 'c') {
	    regression_master = 0;
	} else if (argv[1][1] == 'q') {
	    quiet++;
	} else if (argv[1][1] == 'v') {
	    quiet--;
	} else if (argv[1][1] == 'n') {
	    regression = 0;
	    quiet = -1;
	    no_fork = 1;
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
	    start_subproc_arg_count-=2;
	    if (first_colon) {
	        int ssh_port_int;
		if (second_colon) *second_colon = 0;
		if (sscanf(first_colon+1, "%d", &ssh_port_int) != 1) {
		    second_colon = first_colon;
		}  else {
		    start_subproc_arg_count-=2;
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
	    subproc_args[start_subproc_arg_count] = strdup(SSH_PATH);
	    subproc_args[start_subproc_arg_count+1] = destination_host;
	    if (ssh_port != NULL) {
	        subproc_args[start_subproc_arg_count+2] = "-p";
	        subproc_args[start_subproc_arg_count+3] = ssh_port;
	    }
	    doing_ssh++;
	    subproc_args[cur_subproc_arg] = NULL;
	    argv++; argc--;
	} else if (argv[1][1] == '-') {
	    argv++;
	    argc--;
	    break;
	}
	argv++;
	argc--;
    }
    if (remote_directory[0] != 0) {
        if (strrchr(argv0, '/')) argv0 = strrchr(argv0, '/') + 1;
        subproc_args[cur_subproc_arg] = malloc(strlen(remote_directory) + 
					       strlen(argv0) + 4);
	strcpy(subproc_args[cur_subproc_arg], remote_directory);
	if (remote_directory[strlen(remote_directory)-1] != '/')
	    strcat(subproc_args[cur_subproc_arg], "/");
	strcat(subproc_args[cur_subproc_arg], argv0);
    } else {
        subproc_args[cur_subproc_arg] = argv0;
    }
    subproc_args[cur_subproc_arg+1] = NULL;
    cur_subproc_arg++;
    srand48(getpid());

#ifdef HAVE_WINDOWS_H
    SetTimer(NULL, 5, 1000, (TIMERPROC) fail_and_die);
#else
    struct sigaction sigact;
    sigact.sa_flags = 0;
    sigact.sa_handler = fail_and_die;
    sigemptyset(&sigact.sa_mask);
    sigaddset(&sigact.sa_mask, SIGALRM);
    sigaction(SIGALRM, &sigact, NULL);
    if (regression == 0) {
	alarm(600);
    } else {
	alarm(60);
    }
#endif

    if (!regression_master) {
	ret = be_test_child(argc, argv);
    } else {
	ret = be_test_master(argc-1, &argv[1]);
    }
    free(subproc_args);
    return ret;
}
