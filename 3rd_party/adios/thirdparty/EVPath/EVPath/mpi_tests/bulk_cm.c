#include "config.h"

#include <stdio.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include "evpath.h"
#define MAXPROC 8    /* Max number of procsses */
#define NAMELEN 80   /* Max length of machine name */
#define CONTACTLEN 1024
#define LENGTH 24    /* Length of send buffer is divisible by 2, 4, 6 and 8 */

int quiet = 1;
static atom_t CM_TRANSPORT;


#define MSG_COUNT 1000
#define MSG_SIZE  102400 

typedef struct _bulk_rec {
    long len;
    double *buf;
} bulk_rec, *bulk_rec_ptr;

FMField bulk_fields[] =
{
    {"len", "integer", sizeof(((bulk_rec_ptr)0)[0].len), 
     FMOffset(bulk_rec_ptr, len)},
    {"elem", "float[len]", sizeof(((bulk_rec_ptr)0)[0].buf), FMOffset(bulk_rec_ptr, buf)},
    {(char *) 0, (char *) 0, 0, 0}
};

static int msg_count = MSG_COUNT;
static int msg_size = MSG_SIZE;

static
void 
generate_record(bulk_rec_ptr event)
{
    memset(event, 0, sizeof(*event));
    event->len = msg_size/sizeof(event->buf[0]);
    event->buf = malloc(event->len * sizeof(event->buf[0]));
}

static
void
bulk_handler(CManager cm, CMConnection conn, void *vevent, void *client_data,
	       attr_list attrs)
{
    long sum = 0, scan_sum = 0;
    static int message_count = 0;
    int np;
    MPI_Comm_size(MPI_COMM_WORLD, &np);    /* Get nr of processes */
    (void)cm;
    if ((quiet <= 0) || (sum != scan_sum)) {
	printf("In the handler, connection is %lx\n", (long)conn);
    }
    message_count++;
    if (message_count == (np - 1) * msg_count) {
	printf("All messages received\n");
	CMCondition_signal(cm, (int)(long)client_data);
    }
}

int
main(int argc, char* argv[]) 
{
    int i, np, me;

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
	} else if (strcmp(argv[1], "-k") == 0) {
	    int size = 0;
	    sscanf(argv[2], "%d", &size);
	    msg_size = size * 1024;
	    argc--;
	    argv++;
	} else if (strcmp(argv[1], "-m") == 0) {
	    int size = 0;
	    sscanf(argv[2], "%d", &size);
	    msg_size = size * 1024 * 1024;
	    argc--;
	    argv++;
	} else if (strcmp(argv[1], "-c") == 0) {
	    int count = 0;
	    sscanf(argv[2], "%d", &count);
	    msg_count = count;
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
    CMfork_comm_thread(cm);
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
	MPI_Bcast(master_contact,CONTACTLEN,MPI_CHAR,0,MPI_COMM_WORLD);
	format = CMregister_simple_format(cm, "bulk", bulk_fields, sizeof(bulk_rec));
	message_wait_condition = CMCondition_get(cm, NULL);
	CMregister_handler(format, bulk_handler, (void*)(long)message_wait_condition);
	printf("Master waiting on Barrier\n");
	MPI_Barrier(MPI_COMM_WORLD);
 	start = time(NULL);
	CMCondition_wait(cm, message_wait_condition);
	end = time(NULL);
	if (start == end) end++;
	printf("Duration = %ld\n", (long)(end - start));
	printf("Aggregate gather bandwidth roughly %g Mbytes/sec\n", (double)np * msg_count * msg_size / ((double)(end - start) * 1000000));
    } else { /* all other processes do this */
	attr_list contact_list;
	CMConnection conn = NULL;
	CMFormat format;
	bulk_rec data;
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
	format = CMregister_simple_format(cm, "bulk", bulk_fields, sizeof(bulk_rec));
	generate_record(&data);
	printf("Child %d waiting on Barrier\n", me);
	MPI_Barrier(MPI_COMM_WORLD);
	for (i=0; i < msg_count; i++) {
	    CMwrite(conn, format, &data);
	}
	CMsleep(cm, 1);
    }

    MPI_Finalize();
    exit(0);
}
