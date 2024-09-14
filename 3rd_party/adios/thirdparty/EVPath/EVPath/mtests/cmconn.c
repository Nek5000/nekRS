/*
 *   CM connection test plan of action: 
 *
 *     We want to test to ensure that a passively-opened connection is not
 *     closed until the far side explicitly drops all references to it.
 *
 *  1. Child will initiate a connection with parent, sending it's context
 *   list in a message.
 *  2. Child then waits for message receipt.
 *  3. Master receives message, initiates a connection using the contact list (should be same connection), sends a message, and closes the connection.
 *  4. Child should receive the message, wait 5 seconds to give time for master to close, then send a 2nd message and close the connection (after brief pause).
 *  5. the master must receive the 2nd message for success.
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
#include <sys/wait.h>

#ifdef _MSC_VER
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#define kill(x,y) TerminateProcess(OpenProcess(0,0,(DWORD)x),y)
#endif
typedef struct _msg_rec {
    char *contact_list;
    int message_id;
} msgrec, *msg_rec_ptr;

static FMField msg_field_list[] =
{
    {"contact_list", "string",
     sizeof(char*), FMOffset(msg_rec_ptr, contact_list)},
    {"message_id", "integer",
     sizeof(int), FMOffset(msg_rec_ptr, message_id)},
    {NULL, NULL, 0, 0}
};

int quiet = 1;
int im_the_master = 0;
int master_success = 0;
CMFormat msg_format = NULL;
CMConnection conn_to_master = NULL;

static
void
msg_handler(CManager cm, CMConnection conn, void *vmsg, void *client_data,
	       attr_list attrs)
{
    msg_rec_ptr msg = vmsg;
    msgrec reply = {NULL, 0};
    (void)cm;
    switch(msg->message_id) {
    case 0:  { /* first message from client to master */
        CMConnection conn_to_client;
	attr_list client_contact_list;
	if (!quiet)
	    printf("Master received incoming message from client:\n  Initiating a new conn to him, writing a message and closing it.\n");
	reply.message_id = 1;
	client_contact_list = attr_list_from_string(msg->contact_list);
	conn_to_client = CMget_conn(cm, client_contact_list);
	if (conn != conn_to_client) printf("CONN_EQ MAY BE BROKEN\n");
	free_attr_list(client_contact_list);
	CMwrite(conn_to_client, msg_format, &reply);
	CMusleep(cm, 5000);
	CMConnection_close(conn_to_client);
	break;
    }
    case 1:   /* message from master to client */
	if (!quiet)
	    printf("Client received incoming message from the master  - waiting\n");
	CMsleep(cm, 5);   /* waiting for master to close his connection */
	reply.message_id = 2;
	if (!quiet)
	    printf("Client sending message to the master, he has to get this\n");
	CMwrite(conn_to_master, msg_format, &reply);
	CMusleep(cm, 5000);
	break;
    case 2:   /* last message from client to master */
	if (!quiet)
	    printf("Master got success message... \n");
	master_success = 1;
	break;
    }
}

static int do_regression_master_test();
static int regression = 1;

static atom_t CM_TRANSPORT;

char *transport = NULL;
#include "support.c"

void
do_master_stuff(CManager cm)
{
	attr_list contact_list, listen_list = NULL;
	char *string_list;
	if ((transport = getenv("CMTransport")) != NULL) {
	    if (listen_list == NULL) listen_list = create_attr_list();
	    add_attr(listen_list, CM_TRANSPORT, Attr_String,
		     (attr_value) strdup(transport));
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
	string_list = attr_list_to_string(contact_list);
	printf("Contact list \"%s\"\n", string_list);
	free(string_list);
	msg_format = CMregister_simple_format(cm, "conn_msg", msg_field_list, sizeof(msgrec));
	CMregister_handler(msg_format, msg_handler, NULL);
	CMsleep(cm, 120);
}

int
main(int argc, char **argv)
{
    CManager cm;
    int regression_master = 1;

    PARSE_ARGS();

    srand(getpid());
    CM_TRANSPORT = attr_atom_from_string("CM_TRANSPORT");

    if (regression && regression_master) {
	return do_regression_master_test();
    }
    cm = CManager_create();
    (void) CMfork_comm_thread(cm);

    if (argc == 1) {
	do_master_stuff(cm);
    } else {
	msgrec msg;
	attr_list tmp, listen_list = NULL;
	if (argc == 2) {
	    attr_list contact_list;
	    contact_list = attr_list_from_string(argv[1]);
	    conn_to_master = CMget_conn(cm, contact_list);
	    if (conn_to_master == NULL) {
		printf("No connection, attr list was :");
		dump_attr_list(contact_list);
		printf("\n");
		exit(1);
	    }
	    free_attr_list(contact_list);
	}
	msg_format = CMregister_simple_format(cm, "conn_msg", msg_field_list, sizeof(msgrec));
	CMregister_handler(msg_format, msg_handler, NULL);
	memset(&msg, 0, sizeof(msg));
	msg.message_id = 0;
	if ((transport = getenv("CMTransport")) != NULL) {
	    if (listen_list == NULL) listen_list = create_attr_list();
	    add_attr(listen_list, CM_TRANSPORT, Attr_String,
		     (attr_value) strdup(transport));
	}
	CMlisten(cm);
	CMlisten_specific(cm, listen_list);
	msg.contact_list = attr_list_to_string(tmp = CMget_specific_contact_list(cm, listen_list));
	free_attr_list(tmp);
	if (!quiet)
	    printf("Contact list sent to client is \"%s\"\n", msg.contact_list);
	CMwrite(conn_to_master, msg_format, &msg);
	free(msg.contact_list);
	CMsleep(cm, 20);
    }
    if (conn_to_master) {
	CMConnection_close(conn_to_master);
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
    char *string_list;
    int i;
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
    string_list = attr_list_to_string(contact_list);
    free_attr_list(contact_list);

    args[2] = string_list;

    if (quiet <= 0) {
	if (forked) {
	    printf("Forked a communication thread\n");
	} else {
	    printf("Doing non-threaded communication handling\n");
	}
    }
    srand48(1);

    msg_format = CMregister_simple_format(cm, "conn_msg", msg_field_list, sizeof(msgrec));
    CMregister_handler(msg_format, msg_handler, NULL);
    subproc_proc = run_subprocess(args);

    /* give him time to start */
    for (i=0; i< 10; i++) {
	if (master_success == 1) break;
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
    if (master_success != 1) printf("Master_success == %d\n", master_success);
    return !(master_success == 1);
}
