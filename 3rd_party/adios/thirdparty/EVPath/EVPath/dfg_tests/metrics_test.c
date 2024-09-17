#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>

#include "cod.h"
#include "ev_dfg.h"
#include "test_support.h"

static int status;
static EVclient test_client;

typedef struct _metrics_rec {
    double    dtimeofday;
    int           hw_cpus;
    long           hw_cpu_min_freq;
    long           hw_cpu_max_freq;
    long           hw_cpu_curr_freq;
    char*       os_type;
    char*       os_release;
    char*       hostname;
    double    stat_uptime;
    double    stat_loadavg_one;
    double    stat_loadavg_five;
    double    stat_loadavg_fifteen;
    unsigned long    vm_mem_total;
    unsigned long    vm_mem_free;
    unsigned long    vm_swap_total;
    unsigned long    vm_swap_free;
} metrics_rec, *metrics_rec_ptr;

static FMField metrics_field_list[] =
{
    {"dtimeofday", "double", sizeof(double),
     FMOffset(metrics_rec_ptr, dtimeofday)},
    {"hw_cpus", "integer", sizeof(int),
     FMOffset(metrics_rec_ptr, hw_cpus)},
    {"hw_cpu_min_freq", "integer", sizeof(long),
     FMOffset(metrics_rec_ptr, hw_cpu_min_freq)},
    {"hw_cpu_max_freq", "integer", sizeof(long),
     FMOffset(metrics_rec_ptr, hw_cpu_max_freq)},
    {"hw_cpu_curr_freq", "integer", sizeof(long),
     FMOffset(metrics_rec_ptr, hw_cpu_curr_freq)},
    {"os_type", "string", sizeof(char*),
     FMOffset(metrics_rec_ptr, os_type)},
    {"os_release", "string", sizeof(char*),
     FMOffset(metrics_rec_ptr, os_release)},
    {"hostname", "string", sizeof(char*),
     FMOffset(metrics_rec_ptr, hostname)},
    {"stat_uptime", "double", sizeof(double),
     FMOffset(metrics_rec_ptr, stat_uptime)},
    {"stat_loadavg_one", "double", sizeof(double),
     FMOffset(metrics_rec_ptr, stat_loadavg_one)},
    {"stat_loadavg_five", "double", sizeof(double),
     FMOffset(metrics_rec_ptr, stat_loadavg_five)},
    {"stat_loadavg_fifteen", "double", sizeof(double),
     FMOffset(metrics_rec_ptr, stat_loadavg_fifteen)},
    {"vm_mem_total", "unsigned integer", sizeof(unsigned long),
     FMOffset(metrics_rec_ptr, vm_mem_total)},
    {"vm_mem_free", "unsigned integer", sizeof(unsigned long),
     FMOffset(metrics_rec_ptr, vm_mem_free)},
    {"vm_swap_total", "unsigned integer", sizeof(unsigned long),
     FMOffset(metrics_rec_ptr, vm_swap_total)},
    {"vm_swap_free", "unsigned integer", sizeof(unsigned long),
     FMOffset(metrics_rec_ptr, vm_swap_free)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec metrics_format_list[] =
{
    {"metrics", metrics_field_list, sizeof(metrics_rec), NULL},
    {NULL, NULL}
};

static
int
simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    static int first = 1;
    metrics_rec_ptr event = vevent;
    static metrics_rec first_event;
    if (quiet <= 0) {
	printf("In the handler, event data is :\n");
	printf("	dtimeofday = %g\n", event->dtimeofday);
	printf("	hw_cpus = %d\n", event->hw_cpus);
	printf("	hw_cpu_min_freq = %ld\n", event->hw_cpu_min_freq);
	printf("	hw_cpu_max_freq = %ld\n", event->hw_cpu_max_freq);
	printf("	hw_cpu_curr_freq = %ld\n", event->hw_cpu_curr_freq);
	printf("	os_type = %s\n", event->os_type);
	printf("	os_release = %s\n", event->os_release);
	printf("	hostname = %s\n", event->hostname);
	printf("	stat_uptime = %g\n", event->stat_uptime);
	printf("	stat_loadavg_one = %g\n", event->stat_loadavg_one);
	printf("	stat_loadavg_five = %g\n", event->stat_loadavg_five);
	printf("	stat_loadavg_fifteen = %g\n", event->stat_loadavg_fifteen);
	printf("	vm_mem_total = %ld\n", event->vm_mem_total);
	printf("	vm_mem_free = %ld\n", event->vm_mem_free);
	printf("	vm_swap_total = %ld\n", event->vm_swap_total);
	printf("	vm_swap_free = %ld\n", event->vm_swap_free);
    }
    if (first) {
	first = 0;
	first_event = *event;
	return 0;
    } else {
        if (first_event.dtimeofday >= event->dtimeofday) {
	    EVclient_shutdown(test_client, DFG_STATUS_FAILURE);
        } else {
	    EVclient_shutdown(test_client, DFG_STATUS_SUCCESS);
        }
    }
    (void)cm;
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
    return 0;
}

/* static */
/* int */
/* simple_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) */
/* { */
/*     (void)cm; */
/*     (void)client_data; */
/*     static int count = 0; */
/*     count++; */
/*     if (count == node_count) { */
/* 	EVclient_shutdown(test_client, 0); */
/*     } */
/*     return 0; */
/* } */

char *ECL_generate = "{\n\
    static int count = 0;\n\
    output.dtimeofday = dgettimeofday();\n\
    output.hw_cpus = hw_cpus();\n\
    output.hw_cpu_min_freq = hw_cpu_min_freq();\n\
    output.hw_cpu_max_freq = hw_cpu_max_freq();\n\
    output.hw_cpu_curr_freq = hw_cpu_curr_freq();\n\
    output.os_type = os_type();\n\
    output.os_release = os_release();          \n\
    output.hostname = hostname();          \n\
    output.stat_uptime = stat_uptime();           \n\
    output.stat_loadavg_one = stat_loadavg_one();       \n	\
    output.stat_loadavg_five = stat_loadavg_five();      \n	\
    output.stat_loadavg_fifteen = stat_loadavg_fifteen();   \n	\
    output.vm_mem_total = vm_mem_total();   \n	\
    output.vm_mem_free = vm_mem_free();    \n	\
    output.vm_swap_total = vm_swap_total();   \n	\
    output.vm_swap_free = vm_swap_free();     \n\
    count++;\n\
    return count <= 2;\n\
}";

extern int
be_test_master(int argc, char **argv)
{
    char *nodes[] = {"a", NULL};
    CManager cm;
    char *str_contact;
    char *chandle = NULL;
    EVmaster test_master;
    EVdfg test_dfg;
    EVclient_sinks sink_capabilities;
    char *auto_action_spec;
    EVdfg_stone src, sink;

    cm = CManager_create();
    CMlisten(cm);
    
/*
**  LOCAL DFG SUPPORT   Sources and sinks that might or might not be utilized.
*/
    
    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", metrics_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

/*
**  MASTER AND DFG CREATION
*/
    test_master = EVmaster_create(cm);
    str_contact = EVmaster_get_contact_list(test_master);
    EVmaster_register_node_list(test_master, &nodes[0]);
    test_dfg = EVdfg_create(test_master);
    
    auto_action_spec = create_transform_action_spec(NULL, metrics_format_list, ECL_generate);
    src = EVdfg_create_stone(test_dfg, auto_action_spec);
    free(auto_action_spec);

    EVdfg_assign_node(src, "a");
    EVdfg_enable_auto_stone(src, 1, 0);
    sink = EVdfg_create_sink_stone(test_dfg, "simple_handler");
    EVdfg_assign_node(sink, "a");
    EVdfg_link_dest(src, sink);

    EVdfg_realize(test_dfg);
    
    /* We're node 0 in the DFG */
    test_client = EVclient_assoc_local(cm, nodes[0], test_master, NULL, sink_capabilities);
    
    /* Fork the others */
    test_fork_children(&nodes[0], str_contact);

    if (EVclient_ready_wait(test_client) != 1) {
      /* dfg initialization failed! */
      exit(1);
    }

    
    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    status = EVclient_wait_for_shutdown(test_client);
    free(str_contact);
    wait_for_children(nodes);

    CManager_close(cm);
    if (chandle) free(chandle);
    return status;
}


extern int
be_test_child(int argc, char **argv)
{
    CManager cm;
    EVclient_sinks sink_capabilities = NULL;
    EVclient_sources source_capabilities = NULL;

    cm = CManager_create();
    if (argc != 3) {
	printf("Child usage:  evtest  <nodename> <mastercontact>\n");
	exit(1);
    }

    sink_capabilities = EVclient_register_sink_handler(cm, "simple_handler", metrics_format_list,
				(EVSimpleHandlerFunc) simple_handler, NULL);

    test_client = EVclient_assoc(cm, argv[1], argv[2], source_capabilities, sink_capabilities);

    EVclient_ready_wait(test_client);
    if (EVclient_active_sink_count(test_client) == 0) {
	EVclient_ready_for_shutdown(test_client);
    }

    return EVclient_wait_for_shutdown(test_client);
}
