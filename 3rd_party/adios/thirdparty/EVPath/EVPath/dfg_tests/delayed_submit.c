/*
 *  Test the ability of a CoD handler to reference other stones by name, see notes below for details
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef _MSC_VER
#include <sys/types.h>
#include <sys/timeb.h>
#endif
#ifndef timersub
#define timersub(a, b, result) \
        do { \
                (result)->tv_sec = (a)->tv_sec - (b)->tv_sec; \
                (result)->tv_usec = (a)->tv_usec - (b)->tv_usec; \
                if ((result)->tv_usec < 0) { \
                        --(result)->tv_sec; \
                        (result)->tv_usec += 1000000; \
                } \
        } while (0)
#endif // timersub

#include "evpath.h"
#include "ev_dfg.h"
#include "test_support.h"

typedef struct _outbound {
    struct timeval submit_time;
} outbound, *outbound_ptr;

static FMField time_list[] =
{
    {"submit_sec", "integer", sizeof(((struct timeval*)0)->tv_sec), FMOffset(struct timeval *, tv_sec)}, 
    {"submit_usec", "integer", sizeof(((struct timeval*)0)->tv_usec), FMOffset(struct timeval *, tv_usec)}, 
    {NULL, NULL, 0, 0}
};

static FMStructDescRec time_format_list[] =
{
    {"outbound", time_list, sizeof(struct timeval), NULL},
    {NULL, NULL, 0, NULL}
};

static FMStructDescList queue_list[] = {time_format_list, NULL};

static int status;
static EVclient test_client;

static
int
event_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    outbound_ptr event = vevent;
    static int message_count = 0;
    static int failure = 0;
    struct timeval now, delay;
    (void)cm;
    (void)client_data;
    if (quiet <= 0) printf("In handler for stone %d\n", EVexecuting_stone(cm));
#ifdef HAVE_GETTIMEOFDAY
    gettimeofday((struct timeval*)&now, NULL);
#else
    /* GSE...  No gettimeofday on windows.
     * Must use _ftime, get millisec time, convert to usec.  Bleh.
     */
    struct _timeb nowb;
    _ftime(&nowb);
    ((struct timeval*)&now)->tv_sec = (long)nowb.time;
    ((struct timeval*)&now)->tv_usec = nowb.millitm * 1000;
#endif
    timersub(&now, &event->submit_time, &delay);
    if (quiet <= 0) {
        printf("Now is %ld.%06d, sent %ld.%06d\n", (long)now.tv_sec, (int)now.tv_usec, (long)event->submit_time.tv_sec, (int)event->submit_time.tv_usec);
	printf("Delay is %ld.%06d\n", (long)delay.tv_sec, (int)delay.tv_usec);
    }
    if ((message_count % 2) == 0) {
	if (labs(delay.tv_sec *1000000 + delay.tv_usec - 1000000) > 100000) {
	    printf("Message delayed too much, %ld\n", labs(delay.tv_sec *1000000 + delay.tv_usec - 1000000));
	    failure = 1;
	}
    }
    if ((message_count % 2) == 1) {
	if (labs(delay.tv_sec *1000000 + delay.tv_usec - 1500000) > 100000) {
	    printf("Message delayed too much, %ld\n", labs(delay.tv_sec *1000000 + delay.tv_usec - 1500000));
	    failure = 1;
	}
    }
    message_count++;
    if (message_count >= 2) {
	EVclient_shutdown(test_client, failure);
    }
    return 0;
}

static
int
result_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    outbound *event = vevent;
    (void)cm;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	submit_time = %ld sec, %d usec\n", (long)event->submit_time.tv_sec, 
	       (int)event->submit_time.tv_usec);
    }
    return 0;
}

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", "b", NULL};
    CManager cm;
    char *str_contact;
    EVdfg_stone S1, T1, M1;
    EVsource source_handle;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;
    char *M1_action_spec;
    (void)argc; (void)argv;
    cm = CManager_create();
    CMlisten(cm);

/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/

    source_handle = EVcreate_submit_handle(cm, DFG_SOURCE, time_format_list);
    source_capabilities = EVclient_register_source("master_source", source_handle);
    (void)EVclient_register_sink_handler(cm, "event_handler", time_format_list,
				(EVSimpleHandlerFunc) event_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "result_handler", time_format_list,
				(EVSimpleHandlerFunc) result_handler, NULL);

/*
**  MASTER and DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);

    /* create:
       - one source (S1) to generate an event to trigger the next stone.
       - one multityped-stone (M1) whose action will be to snag the current time 
         via gettimeofday() and then submit a delayed event containing that time.
	 Actually, we'll do this twice, once submitting the incoming event, and 
	 once submitting a created event.
       - a terminal stone (T1) that will catch the events, compare the time 
         in the event with the current time and determine if an appropriate delay 
	 has been performed.
    */

char *COD_multi = "{\n\
    outbound m2, *m1;\n\
    timeval t;\n\
    timeval delay;\n\
    gettimeofday(&t);\n\
    m1 = EVdata_outbound(0);\n\
    m2.submit_sec = m1.submit_sec = t.tv_sec;\n\
    m2.submit_usec = m1.submit_usec = t.tv_usec;\n	       \
    delay.tv_sec = 1;\n\
    delay.tv_usec = 0;\n\
    EVsubmit_delayed(0, m1, (void*)0, &delay);\n\
    delay.tv_usec = 500000;\n\
    EVsubmit_delayed(0, m2, (void*)0, &delay);\n\
    EVdiscard_outbound(0);\n\
}";


    test_dfg = EVdfg_create(test_master);
    S1 = EVdfg_create_source_stone(test_dfg, "master_source");
    T1 = EVdfg_create_sink_stone(test_dfg, "event_handler");

    M1_action_spec = create_multityped_action_spec(queue_list, COD_multi);
    M1 = EVdfg_create_stone(test_dfg, M1_action_spec);
    free(M1_action_spec);
    EVdfg_assign_node(S1, "b");
    EVdfg_assign_node(T1, "a");
    EVdfg_assign_node(M1, "a");
    EVdfg_link_dest(S1, M1);
    EVdfg_link_dest(M1, T1);

    EVdfg_realize(test_dfg);

/* We're node 0 in the process group */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, source_capabilities, sink_capabilities);

/* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    free(str_contact);
    if (EVclient_ready_wait(test_client) != 1) {
	/* dfg initialization failed! */
	exit(1);
    }

    EVclient_ready_for_shutdown(test_client);
    if (EVclient_source_active(source_handle)) {
	while (!EVclient_test_for_shutdown(test_client)) {
	    outbound rec;
	    /* submit would be quietly ignored if source is not active */
	    EVsubmit(source_handle, &rec, NULL);
	    CMsleep(cm, 3);
	}
	/* exit if shutdown is imminent */
    }

    status = EVclient_wait_for_shutdown(test_client);

    wait_for_children(nodes);

    EVfree_source(source_handle);
    CManager_close(cm);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVsource src;
    EVclient_sinks sink_capabilities;
    EVclient_sources source_capabilities;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  delayed_submit  <nodename> <mastercontact>\n");
	exit(1);
    }

    src = EVcreate_submit_handle(cm, DFG_SOURCE, time_format_list);
    source_capabilities = EVclient_register_source("master_source", src);
    (void)EVclient_register_sink_handler(cm, "event_handler", time_format_list,
				(EVSimpleHandlerFunc) event_handler, NULL);
    sink_capabilities = EVclient_register_sink_handler(cm, "result_handler", time_format_list,
				(EVSimpleHandlerFunc) result_handler, NULL);
    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);
    EVclient_ready_wait(test_client);

    EVclient_ready_for_shutdown(test_client);
    if (EVclient_source_active(src)) {
	while (!EVclient_test_for_shutdown(test_client)) {
	    outbound rec;
	    /* submit would be quietly ignored if source is not active */
	    EVsubmit(src, &rec, NULL);
	    CMsleep(cm, 3);
	}
	/* exit if shutdown is imminent */
    }

    EVfree_source(src);
    return EVclient_wait_for_shutdown(test_client);
}
