#include "config.h"
#include <stdio.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>
#include <math.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#define __ANSI_CPP__
#else
#include <netinet/in.h>
#include <arpa/inet.h>
#endif
#include <ffs.h>
#include <atl.h>
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"

extern atom_t CM_REBWM_RLEN;
extern atom_t CM_REBWM_REPT;
extern atom_t CM_BW_MEASURE_INTERVAL;
extern atom_t CM_BW_MEASURE_TASK;
extern atom_t CM_BW_MEASURED_VALUE;
extern atom_t CM_BW_MEASURED_COF;
extern atom_t CM_BW_MEASURE_SIZE;
extern atom_t CM_BW_MEASURE_SIZEINC;
static void init_atoms();

static atom_t CM_TRANS_TEST_SIZE = -1;
static atom_t CM_TRANS_TEST_NODE = -1;
static atom_t CM_TRANS_TEST_VECS = -1;
static atom_t CM_TRANS_TEST_VERBOSE = -1;
static atom_t CM_TRANS_TEST_REPEAT = -1;
static atom_t CM_TRANS_TEST_REUSE_WRITE_BUFFER = -1;
static atom_t CM_TRANS_TEST_DURATION = -1;
static atom_t CM_TRANS_MEGABITS_SEC = -1;

#define CMPerfProbe (unsigned int) 0xf0
#define CMPerfProbeResponse (unsigned int) 0xf1
#define CMPerfBandwidthInit (unsigned int) 0xf2
#define CMPerfBandwidthBody (unsigned int) 0xf3
#define CMPerfBandwidthEnd  (unsigned int) 0xf4
#define CMPerfBandwidthResult  (unsigned int) 0xf5
#define CMPerfTestInit (unsigned int) 0xfa
#define CMPerfTestBody (unsigned int) 0xfb
#define CMPerfTestEnd  (unsigned int) 0xfc
#define CMPerfTestResult  (unsigned int) 0xfd

#define CMRegressivePerfBandwidthInit (unsigned int) 0xf6
#define CMRegressivePerfBandwidthBody (unsigned int) 0xf7
#define CMRegressivePerfBandwidthEnd  (unsigned int) 0xf8
#define CMRegressivePerfBandwidthResult  (unsigned int) 0xf9

void
CMdo_performance_response(CMConnection conn, size_t length, int func,
			  int byte_swap, char *buffer)
{
    /* part of length was read already */
    length += 8;
    CMtrace_out(conn->cm, CMControlVerbose, "CMDo_performance_response func %d \n", func);
    init_atoms();
    switch(func) {
    case CMPerfProbe:
	/* first half of latency probe arriving */
	{
	    struct FFSEncodeVec tmp_vec[2];
	    int32_t header[3];
	    long actual;
	    tmp_vec[0].iov_base = &header;
	    tmp_vec[0].iov_len = sizeof(header);
	    header[0] = 0x434d5000;  /* CMP\0 */
#if SIZEOF_LONG == 4
	    header[1] = length & 0xffffff;
	    header[2] = 0;
#else
	    header[1] = (length>>32) & 0xffffff;
	    header[2] = length & 0xffffffff;
#endif
	    header[1] = header[1] | (CMPerfProbeResponse << 24);

	    tmp_vec[1].iov_len = length - sizeof(header);
	    tmp_vec[1].iov_base = buffer;

	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - responding to latency probe of %zd bytes\n", length);
	    actual = INT_CMwrite_raw(conn, tmp_vec, tmp_vec + 1, 2, length, NULL, 0);
	    if (actual != 2) {
		printf("perf write failed\n");
	    }
	}
	break;
    case CMPerfProbeResponse:
	/* last half of latency probe arriving, probe completion*/
	{
	    int cond = *(int*)buffer;  /* first entry should be condition */
	    chr_time *timer = INT_CMCondition_get_client_data(conn->cm, cond);
	    
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - latency probe response, condition %d\n", cond);
	    chr_timer_stop(timer);
	    INT_CMCondition_signal(conn->cm, cond);
	}
	break;
    case CMPerfBandwidthInit:
	/* initiate bandwidth measure */
	chr_timer_start(&conn->bandwidth_start_time);
	CMtrace_out(conn->cm, CMTransportVerbose, "CM - Starting bandwidth probe\n");
	break;
    case CMPerfBandwidthBody:
	/* no activity for inner packets */
	CMtrace_out(conn->cm, CMTransportVerbose, "CM - bandwidth probe - body packet\n");
	break;
    case CMPerfBandwidthEnd:
	/* end bandwidth measure, send result back */
	{
	    int header[6];
	    int actual;
	    struct FFSEncodeVec tmp_vec[1];
	    chr_timer_stop(&conn->bandwidth_start_time);
	    union {
		double d;
		float f;
		int i[2];
	    } t;

	    header[0] = 0x434d5000;  /* CMP\0 */
	    header[1] = 0 | (CMPerfBandwidthResult << 24);
	    header[2] = sizeof(header);
	    header[3] = *(int*)buffer;  /* first entry should be condition */
	    t.d = chr_time_to_secs(&conn->bandwidth_start_time);
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Completing bandwidth probe - %g seconds to receive\n", t.d);
	    if (htonl(0xdeadbeef) == 0xdeadbeef) {
		header[4] = htonl(t.i[0]);
		header[5] = htonl(t.i[1]);
	    } else {
		header[4] = htonl(t.i[1]);
		header[5] = htonl(t.i[0]);
	    }
	    tmp_vec[0].iov_base = &header;
	    tmp_vec[0].iov_len = sizeof(header);
	    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, sizeof(header), NULL, 0);

	    if (actual != 1) {
		printf("perf write failed\n");
	    }
	}
	break;
    case CMPerfBandwidthResult:
	/* result of bandwidth measure arriving, wake waiting guy */
	{
	    int cond = *(int*)buffer;  /* first entry should be condition */
	    double *result_p = INT_CMCondition_get_client_data(conn->cm, cond);
	    union {
		double d;
		float f;
		int i[2];
	    } t;

	    if (htonl(0xdeadbeef) == 0xdeadbeef) {
		t.i[0] = ntohl(((int*)buffer)[1]);
		t.i[1] = ntohl(((int*)buffer)[2]);
	    } else {
		t.i[0] = ntohl(((int*)buffer)[2]);
		t.i[1] = ntohl(((int*)buffer)[1]);
	    }
	    if (result_p) *result_p = t.d;
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - bandwidth probe response, condition %d\n", cond);
	    INT_CMCondition_signal(conn->cm, cond);
	}
	break;
    case CMPerfTestInit: 
        {
	    /* initiate bandwidth measure */
	    int header_size = ((int*)buffer)[1];  /* first entry should be condition */
	    char *attr_string = buffer + header_size - 12;
	    attr_list list = attr_list_from_string(attr_string);
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Starting transport test\n");
	    if (conn->cm->perf_upcall)  {
		CManager_unlock(conn->cm);
		(void)conn->cm->perf_upcall(conn->cm, buffer, length - 8, 0, list);
		CManager_lock(conn->cm);
	    }
	    free_attr_list(list);
	    chr_timer_start(&conn->bandwidth_start_time);
	    break;
	}
    case CMPerfTestBody:
	/* no activity for inner packets */
	CMtrace_out(conn->cm, CMTransportVerbose, "CM - transport test - body packet\n");
	if (conn->cm->perf_upcall)  {
	    CManager_unlock(conn->cm);
	    (void)conn->cm->perf_upcall(conn->cm, buffer, length - 8, 1, NULL);
	    CManager_lock(conn->cm);
	}
	break;
    case CMPerfTestEnd:
	/* end bandwidth measure, send result back */
	{
	    int header[6];
	    int actual;
	    attr_list upcall_result, upcall_feed;
	    char *str_list = NULL;
	    struct FFSEncodeVec tmp_vec[2];
	    chr_timer_stop(&conn->bandwidth_start_time);

	    header[0] = 0x434d5000;  /* CMP\0 */
	    header[1] = 0 | (CMPerfTestResult << 24);
	    header[2] = sizeof(header);
	    header[3] = *(int*)buffer;  /* first entry should be condition */
	    header[4] = 0;
	    header[5] = 0;
	    tmp_vec[0].iov_base = &header;
	    tmp_vec[0].iov_len = sizeof(header);
	    tmp_vec[1].iov_base = NULL;
	    tmp_vec[1].iov_len = 0;
	    upcall_feed = create_attr_list();

	    set_double_attr(upcall_feed, CM_TRANS_TEST_DURATION,
			    chr_time_to_secs(&conn->bandwidth_start_time));
	    if (conn->cm->perf_upcall)  {
		CManager_unlock(conn->cm);
		upcall_result = conn->cm->perf_upcall(conn->cm, buffer, length - 8, 2, upcall_feed);
		CManager_lock(conn->cm);
		if (upcall_result) {
		    str_list = attr_list_to_string(upcall_result);
		    free_attr_list(upcall_result);
		    tmp_vec[1].iov_len = header[4] = (int)strlen(str_list) + 1;
		    tmp_vec[1].iov_base = str_list;
		    header[2] += header[4];
		}
	    }
	    free_attr_list(upcall_feed);
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - transport test response sent:");
	    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 2, sizeof(header) + tmp_vec[1].iov_len, NULL, 0);
	    if (str_list) free(str_list);
	    if (actual != 1) {
		printf("perf write failed\n");
	    }
	}
	break;
    case CMPerfTestResult:
	/* result of bandwidth measure arriving, wake waiting guy */
	{
	    int cond = *(int*)buffer;  /* first entry should be condition */
	    attr_list *result_p = INT_CMCondition_get_client_data(conn->cm, cond);

	    if (ntohl(((int*)buffer)[1]) != 0) {
		char *attr_string = (char*)&((int*)buffer)[3];
		attr_list upcall_result = attr_list_from_string(attr_string);
		if (result_p) *result_p = upcall_result;
	    }
	    CMtrace_out(conn->cm, CMConnectionVerbose, "CM - transport test response, condition %d\n", cond);
	    INT_CMCondition_signal(conn->cm, cond);
	}
	break;
    case CMRegressivePerfBandwidthInit:
	/* initiate bandwidth measure */
        CMtrace_out(conn->cm, CMConnectionVerbose, "CM - received CM bw measure initiate\n");
	chr_timer_start(&conn->regressive_bandwidth_start_time);
	break;
    case CMRegressivePerfBandwidthBody:
	/* no activity for inner packets */
	break;
    case CMRegressivePerfBandwidthEnd:
	/* first half of latency probe arriving */
	{
	    int header[5];
	    int actual;
	    struct FFSEncodeVec tmp_vec[1];
	    chr_timer_stop(&conn->regressive_bandwidth_start_time);

	    header[0] = 0x434d5000;  /* CMP\0 */
	    header[1] = 0 | (CMRegressivePerfBandwidthResult << 24);
	    header[2] = sizeof(header);
	    header[3] = *(int*)buffer;  /* first entry should be condition */
	    header[4] = (int) 
		chr_time_to_microsecs(&conn->regressive_bandwidth_start_time);
            CMtrace_out(conn->cm, CMConnectionVerbose, "CM - received CM bw measure end, condition %d\n", *(int*)buffer);
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Completing bandwidth probe - %d microseconds to receive\n", header[2]);

	    tmp_vec[0].iov_base = &header;
	    tmp_vec[0].iov_len = sizeof(header);
	    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, sizeof(header), NULL, 0);
	    if (actual != 1) {
		printf("perf write failed\n");
	    }
	}
	break;
    case CMRegressivePerfBandwidthResult:
	/* last half of latency probe arriving, probe completion*/
	{
	    int cond = *(int*)buffer;  /* first entry should be condition */
	    int time;
	    char *chr_time, tmp;
	    int *result_p = INT_CMCondition_get_client_data(conn->cm, cond);
	    
	    time = ((int*)buffer)[1];/* second entry should be condition */
	    if (byte_swap) {
		chr_time = (char*)&time;
		tmp = chr_time[0];
		chr_time[0] = chr_time[3];
		chr_time[3] = tmp;
		tmp = chr_time[1];
		chr_time[1] = chr_time[2];
		chr_time[2] = tmp;
	    }
	    *result_p = time;
	    CMtrace_out(conn->cm, CMTransportVerbose, "CM - bandwidth probe response, condition %d\n", cond);
	    INT_CMCondition_signal(conn->cm, cond);
	}
	break;
    default:
	printf("BAD!  unknown perf function %d\n", func);
    }
}

static long
do_single_probe(CMConnection conn, long size, attr_list attrs)
{
    int cond;
    static long max_block_size = 0;
    static char *block = NULL;
    chr_time round_trip_time;
    long actual;
    struct FFSEncodeVec tmp_vec[1];

    (void)attrs;
    cond = INT_CMCondition_get(conn->cm, conn);

    if (size < 12) size = 12;
    if (max_block_size == 0) {
	char *new_block = malloc(size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = size;
	memset(block, 0xef, size);
    } else if (size > max_block_size) {
	char *new_block = realloc(block, size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = size;
	memset(block, 0xef, size);
    }
    
    /* CMP\0 in first entry for CMPerformance message */
    ((int*)block)[0] = 0x434d5000;
    /* size in second entry, high byte gives CMPerf operation */
#if SIZEOF_LONG == 4
    ((int*)block)[1] = (size & 0xffffff) | (CMPerfProbe<<24);;
    ((int*)block)[2] = 0;
#else
    ((int*)block)[1] = (size >>32) & 0xffffff;
    ((int*)block)[2] = size;
#endif
    ((int*)block)[3] = cond;   /* condition value in fourth entry */
    
    INT_CMCondition_set_client_data( conn->cm, cond, &round_trip_time);

    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Initiating latency probe of %ld bytes\n", size);
    chr_timer_start(&round_trip_time);

    tmp_vec[0].iov_base = &block[0];
    tmp_vec[0].iov_len = size;
    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);

    if (actual != 1) return -1;

    INT_CMCondition_wait(conn->cm, cond);
    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Completed latency probe - result %g microseconds\n", chr_time_to_microsecs(&round_trip_time));
    return (long) chr_time_to_microsecs(&round_trip_time);
}

/* return units are microseconds */
extern long
INT_CMprobe_latency(CMConnection conn, int size, attr_list attrs)
{
    int i;
    long result = 0;
    int repeat_count = 5;
    for (i=0; i < 2; i++) {
	(void) do_single_probe(conn, size, attrs);
    }
    for (i=0; i < repeat_count; i++) {
	result += do_single_probe(conn, size, attrs);
    }
    result /= repeat_count;
    return result;
}

/* return units are Kbytes/sec */
extern double
INT_CMprobe_bandwidth(CMConnection conn, long size, attr_list attrs)
{
    int i;
    int cond;
    int repeat_count = 100000/size;  /* send about 100K */
    static long max_block_size = 0;
    static char *block = NULL;
    double secs_to_receive;
    long actual;
    double bandwidth;
    struct FFSEncodeVec tmp_vec[1];


    cond = INT_CMCondition_get(conn->cm, conn);

    if (size < 24) size = 24;
    if (repeat_count < 10) repeat_count = 10;
    if (max_block_size == 0) {
	char *new_block = malloc(size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = size;
	memset(block, 0xef, size);
    } else if (size > max_block_size) {
	char *new_block = realloc(block, size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = size;
	memset(block, 0xef, size);
    }
    
    /* CMP\0 in first entry for CMPerformance message */
    ((int*)block)[0] = 0x434d5000;
    /* size in second entry, high byte gives CMPerf operation */
#if SIZEOF_LONG == 4
    ((int*)block)[1] = 0 | (CMPerfBandwidthInit<<24);
    ((int*)block)[2] = size;
#else
    ((int*)block)[1] = ((size >> 32) &0xffffff) | (CMPerfBandwidthInit<<24);
    ((int*)block)[2] = (size & 0xffffffff);
#endif
    ((int*)block)[3] = cond;   /* condition value in third entry */
    
    INT_CMCondition_set_client_data( conn->cm, cond, &secs_to_receive);

    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Initiating bandwidth probe of %ld bytes, %d messages\n", size, repeat_count);
    tmp_vec[0].iov_base = &block[0];
    tmp_vec[0].iov_len = size;
    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);
    if (actual != 1) { 
	return -1;
    }

    ((int*)block)[1] = (((int*)block)[1]&0xffffff) | (CMPerfBandwidthBody<<24);
    for (i=0; i <(repeat_count-1); i++) {
	actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);
	if (actual != 1) {
	    return -1;
	}
    }

    ((int*)block)[1] = (((int*)block)[1]&0xffffff) | (CMPerfBandwidthEnd<<24);
    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);
    if (actual != 1) {
	return -1;
    }

    INT_CMCondition_wait(conn->cm, cond);
    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Completed bandwidth probe - result %g seconds\n", secs_to_receive);
    bandwidth = ((double) size * (double)repeat_count) / secs_to_receive;
    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Estimated bandwidth - %g Mbites/sec\n", bandwidth / 1000.0 * 1000.0 * 8);
    return  bandwidth;
}


static void
init_atoms()
{
    if (CM_TRANS_TEST_SIZE==-1) {
	CM_TRANS_TEST_SIZE = attr_atom_from_string("CM_TRANS_TEST_SIZE");
	CM_TRANS_TEST_NODE = attr_atom_from_string("CM_TRANS_TEST_NODE");
	CM_TRANS_TEST_VECS = attr_atom_from_string("CM_TRANS_TEST_VECS");
	CM_TRANS_TEST_VERBOSE = attr_atom_from_string("CM_TRANS_TEST_VERBOSE");
	CM_TRANS_TEST_REPEAT = attr_atom_from_string("CM_TRANS_TEST_REPEAT");
	CM_TRANS_TEST_REUSE_WRITE_BUFFER = attr_atom_from_string("CM_TRANS_TEST_REUSE_WRITE_BUFFER");
	CM_TRANS_TEST_DURATION = attr_atom_from_string("CM_TRANS_TEST_DURATION_SECS");
	CM_TRANS_MEGABITS_SEC = attr_atom_from_string("CM_TRANS_MEGABITS_SEC");
	(void)CM_TRANS_MEGABITS_SEC;
    }
}

struct _free_struct {
    int free_all;
    struct FFSEncodeVec *write_vec;
    int vecs;
    struct FFSEncodeVec *tmp_vec;
};

static void
write_is_done(void *vdata)
{
    struct _free_struct *data = (struct _free_struct *)vdata;
    int count;

    free(data->write_vec[0].iov_base);
    if (data->tmp_vec) {
	for (count = 0; count < data->vecs; count++) {
	    if (data->tmp_vec[count].iov_base == data->write_vec[0].iov_base) {
		continue;
	    }
	    free(data->tmp_vec[count].iov_base);
	    data->tmp_vec[count].iov_base = NULL;
	}
	free(data->tmp_vec);
    }
    free(data->write_vec);
    free(data);
}

extern attr_list
INT_CMtest_transport(CMConnection conn, attr_list how)
{
    int i;
    int cond;
    attr_list result = NULL;
    long actual;
    struct FFSEncodeVec *write_vec;
    struct FFSEncodeVec *tmp_vec, *header_vec;
    int header[6];
    ssize_t size;
    int vecs = 1;
    int verbose = 0;
    int repeat_count = 1;
    int reuse_write_buffer = 1;
    long start_size, count;
    init_atoms();
    cond = INT_CMCondition_get(conn->cm, conn);
    CManager cm = conn->cm;
    struct _free_struct *write_data;
    int node_id;

    if (!get_long_attr(how, CM_TRANS_TEST_SIZE, &size)) {
	printf("CM_TRANS_TEST_SIZE attr not found by CMtest_transport, required\n");
	return 0;
    }
    get_int_attr(how, CM_TRANS_TEST_VECS, &vecs);
    if (vecs < 1) {
	printf("Stupid vecs value in CMtest_transport, %d\n", vecs);
	return 0;
    }
    if (((float)size / (float) vecs) < 24.0) {
      vecs = 1;
      if (size < 24) {
	size = 24;
      }
    }
    get_int_attr(how, CM_TRANS_TEST_VERBOSE, &verbose);
    get_int_attr(how, CM_TRANS_TEST_REPEAT, &repeat_count);
    get_int_attr(how, CM_TRANS_TEST_REUSE_WRITE_BUFFER, &reuse_write_buffer);
    get_int_attr(how, CM_TRANS_TEST_NODE, &node_id);
    char *attr_str = attr_list_to_string(how);

    
    start_size = (int)strlen(attr_str) + (int)sizeof(header) + 1;
    /* CMP\0 in first entry for CMPerformance message */
    ((int*)header)[0] = 0x434d5000;
    /* size in second entry, high byte gives CMPerf operation */
#if SIZEOF_LONG == 4
    ((int*)header)[1] = (CMPerfTestInit<<24);
#else
    ((int*)header)[1] = ((start_size >> 32) &0xffffff) | (CMPerfTestInit<<24);
#endif
    ((int*)header)[2] = (start_size & 0xffffffff);
    ((int*)header)[3] = cond;   /* condition value in third entry */
    ((int*)header)[4] = sizeof(header);   /* header size 4th entry */
    ((int*)header)[5] = 0;
    INT_CMCondition_set_client_data( conn->cm, cond, &result);

    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Initiating transport test of %zd bytes, %d messages\n", size, repeat_count);

    CMtrace_out(conn->cm, CMTransportVerbose, "CM - transport test, sending first message\n");

    header_vec = malloc(sizeof(header_vec[0]) * (vecs + 1));  /* at least 2 */
    header_vec[0].iov_base = &header[0];
    header_vec[0].iov_len = sizeof(header);
    header_vec[1].iov_base = attr_str;
    header_vec[1].iov_len = strlen(attr_str) + 1; /* send NULL */
    actual = INT_CMwrite_raw(conn, header_vec, NULL, 2, header_vec[0].iov_len + header_vec[1].iov_len, NULL, 1);
    free(attr_str);
    if (actual != 1) { 
	free(header_vec);
	return NULL;
    }

    tmp_vec = NULL;
    size_t each = (size + vecs - 1) / vecs;
    for (i=0; i <repeat_count; i++) {
	if (tmp_vec == NULL) {
	    tmp_vec = malloc(sizeof(tmp_vec[0]) * (vecs + 2));  /* at least 2 */
	    tmp_vec[0].iov_len = 5 * sizeof(int); /* body header */
	    for (count = 0; count < vecs; count++) {
		tmp_vec[count+1].iov_base = calloc(each + repeat_count, 1);
		tmp_vec[count+1].iov_len = each;
	    }
	    for (count = 0; count < vecs; count++) {
		/* for each vector, give it unique data */
		int j;
		for (j=0; j < ((each + repeat_count) /sizeof(int)); j++) {
		    ((int*)tmp_vec[count+1].iov_base)[j] = rand();
		}
	    }
	    if (tmp_vec[1].iov_len > tmp_vec[0].iov_len) {
	      tmp_vec[1].iov_len -= tmp_vec[0].iov_len;
	    } else {
	      tmp_vec[1].iov_len = 1;  /* just so there's something */
	    }
	}
	tmp_vec[0].iov_base = malloc(5 * sizeof(int)); /* body header */
	((int*)tmp_vec[0].iov_base)[0] = 0x434d5000;
    /* size in second entry, high byte gives CMPerf operation */
#if SIZEOF_LONG == 4
	((int*)tmp_vec[0].iov_base)[1] = (CMPerfTestBody<<24);
#else
	((int*)tmp_vec[0].iov_base)[1] = ((size >> 32) &0xffffff) | (CMPerfTestBody<<24);
#endif
	((int*)tmp_vec[0].iov_base)[2] = (size & 0xffffffff);
	((int*)tmp_vec[0].iov_base)[3] = i;   /* sequence number */
	((int*)tmp_vec[0].iov_base)[4] = node_id;   /* node_id */
	if (vecs > 1) {
	    /* if more than one vec, handle rounding */
	    tmp_vec[vecs].iov_len = size - (each * (vecs-1));
	}
	write_vec = malloc(sizeof(write_vec[0]) * (vecs + 2));  /* at least 3 */
	memcpy(write_vec, tmp_vec, sizeof(write_vec[0]) * (vecs + 2));
	for (count = 0; count < vecs; count++) {
	    /* On each iteration, increment the write base for each buffer by 1 */
//	    write_vec[count+1].iov_base += i;
	}
	write_data = malloc(sizeof(struct _free_struct));
	write_data->write_vec = write_vec;
	if ((i == (repeat_count-1)) || (!reuse_write_buffer)){
	    /* free this when done */
	    write_data->tmp_vec = tmp_vec;
	    write_data->vecs = vecs;
	} else {
	    write_data->tmp_vec = NULL;
	}
	actual = INT_CMwrite_raw_notify(conn, write_vec, NULL, vecs+1, size, 
					NULL, 0, write_is_done, (void*)write_data);
	if ((i == (repeat_count-1)) || (!reuse_write_buffer)){
	    /* free this when done */
	    tmp_vec = NULL;
	}

	if (actual != 1) {
	    free(tmp_vec);
	    return NULL;
	}
	if (conn->write_pending) {
	    wait_for_pending_write(conn);
	}
    }
    ((int*)header)[1] = 0 | (CMPerfTestEnd<<24);
    ((int*)header)[2] = sizeof(header);
    if (!tmp_vec) tmp_vec = malloc(sizeof(tmp_vec[0]));
    tmp_vec[0].iov_base = &header[0];
    tmp_vec[0].iov_len = sizeof(header);
    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, sizeof(header), NULL, 0);
    free(tmp_vec);
    free(header_vec);
    if (actual != 1) {
	return NULL;
    }

    int ret = INT_CMCondition_wait(conn->cm, cond);
    if (ret) {
      CMtrace_out(cm, CMTransportVerbose, "CM - Completed transport test - result %p \n", result);
    } else {
      CMtrace_out(cm, CMTransportVerbose, "CM - Completed transport test CONNECTION FAILED- result %p \n", result);
    }
    return result;
}



/* matrix manipulation function for regression */

/***********************************************************
Name:     AtA
Passed:   **A, a matrix of dim m by n 
Returns:  **R, a matrix of dim m by p
***********************************************************/
static void 
AtA(double **A, int m, int n, double **R)
{
    int i, j, k;

    for (i = 0; i < n; i++)
	for (j = 0; j < n; j++) {
	    R[i][j] = 0.0;
	    for (k = 0; k < m; k++)
		R[i][j] += A[k][i] * A[k][j];
	}
}

/***********************************************************
Name:     Atf
Passed:   **A, a matrix of dim n by p
*f, a vector of dim n
Returns:  **R, a matrix of dim p by 1
***********************************************************/
static void 
Atf(double **A, double *f, int n, int p, double **R)
{
    int i, j, k;

    for (i = 0; i < p; i++)
	for (j = 0; j < 1; j++) {
	    R[i][j] = 0.0;
	    for (k = 0; k < n; k++)
		R[i][j] += A[k][i] * f[k];
	}
}
/*********************************************************************
 * dp_inv: inverse matrix A
 *              matrix A is replaced by its inverse
 *		returns ks, an integer
 *    A - square matrix of size n by n
 *    ks equals n when A is invertable
 * example usage:
	  if(dp_inv(N,u) != u){
            printf("Error inverting N\n");
            exit(-1);}
***********************************************************/
static int 
dp_inv(double **A, int n)
{
    int i, j, k, ks;

    ks = 0;
    for (k = 0; k < n; k++) {
	if (A[k][k] != 0.0) {
	    ks = ks + 1;
	    for (j = 0; j < n; j++) {
		if (j != k)
		    A[k][j] = A[k][j] / A[k][k];
	    }			/* end for j */
	    A[k][k] = 1.0 / A[k][k];
	}			/* end if */
	for (i = 0; i < n; i++) {
	    if (i != k) {
		for (j = 0; j < n; j++) {
		    if (j != k)
			A[i][j] = A[i][j] - A[i][k] * A[k][j];
		}		/* end for j */
		A[i][k] = -A[i][k] * A[k][k];
	    }			/* end if */
	}			/* end for i */
    }				/* end for k */
    return ks;
}				/* end MtxInverse() */


/***********************************************************
Name:     AtB
Passed:   **A, a matrix of dim m by n 
**B, a matrix of dim m by p
Returns:  **R, a matrix of dim n by p
***********************************************************/
static void 
AtB(double **A, double **B, int m, int n, int p, double **R)
{
    int i, j, k;

    for (i = 0; i < n; i++)
	for (j = 0; j < p; j++) {
	    R[i][j] = 0.0;
	    for (k = 0; k < m; k++)
		R[i][j] += A[k][i] * B[k][j];
	}
}

/************************************************************************/
/* Frees a double pointer array mtx[row][] */
/************************************************************************/
static void 
dub_dp_free(double **mtx, int row)
{
    int tmp_row;
    for (tmp_row = 0; tmp_row < row; tmp_row++)
	free(mtx[tmp_row]);
    free(mtx);
}

/************************************************************************/
/* Allocate memory for a double pointer array mtx[row][col], */
/* return double pointer  */
/************************************************************************/
static double **
dub_dp_mtxall(int row, int col)
{
    double **mtx;
    int tmp_row;

    /* Set Up Row of Pointers */
    mtx = (double **) malloc((unsigned) (row) * sizeof(double *));
    if (mtx == NULL)
	return NULL;

    /* Set Up Columns in Matrix */
    for (tmp_row = 0; tmp_row < row; tmp_row++) {
	mtx[tmp_row] = (double *) malloc((unsigned) (col) * sizeof(double));
	/* If could not Allocate All Free Memory */
	if (mtx[tmp_row] == NULL) {
	    dub_dp_free(mtx, row);
	    return NULL;
	}			/* Return Null Pointer */
    }
    return mtx;			/* Return Pointer to Matrix */
}


static double 
Regression(CMConnection conn, DelaySizeMtx *inputmtx)
{

    /* CMatrixUtils mtxutl; */
    double *Y, **X, **a;
    double **XtX, **XtY;
    int j;			/* , i, k; */
    /* char plinkid[9]; */
    int numofsizes;
    (void)conn;
    numofsizes = inputmtx->MsgNum;

    Y = (double *) malloc(sizeof(double) * inputmtx->MsgNum);
    X = dub_dp_mtxall(inputmtx->MsgNum, 2);
    a = dub_dp_mtxall(2, 1);
    XtX = dub_dp_mtxall(2, 2);
    XtY = dub_dp_mtxall(2, 1);

    /* regression estimate using least square method */
    for (j = 0; j < numofsizes; j++) {
	/* convert to one way delay in milliseconds */
	Y[j] = inputmtx->AveRTDelay[j] / 2 * 1000;
	/* convert to Kbytes */
	X[j][0] = inputmtx->MsgSize[j] / 1024;
	X[j][1] = 1;
    }
    AtA(X, numofsizes, 2, XtX);
    Atf(X, Y, numofsizes, 2, XtY);
    if (!dp_inv(XtX, 2)) {
	CMtrace_out(conn->cm, CMTransportVerbose, "CM - Regression()- Matrix XtX is not invertible\n");
	return -1;
    }
    AtB(XtX, XtY, 2, 2, 1, a);
    CMtrace_out(conn->cm, CMTransportVerbose,"CM - Regression():\nslope = %f (Bandwidth = %f Mbps), intercept = %f\n", a[0][0], 8 / a[0][0], a[1][0]);

    dub_dp_free(X, numofsizes);
    return 8 / a[0][0];
}

/* return units are Mbps */
extern double
INT_CMregressive_probe_bandwidth(CMConnection conn, long size, attr_list attrs)
{
    int i, j;
    int cond;
    int N = 9; /* send out N tcp streams with varied length*/
    int repeat_count = 100000/size;  /* send about 100K */
    repeat_count = 100;
    static long max_block_size = 0;
    static char *block = NULL;
    int microsecs_to_receive;
    int actual;
    double bandwidth;
    long biggest_size;
    DelaySizeMtx dsm;
    double ave_delay=0.0, var_delay=0.0;
    double ave_size=0.0, var_size=0.0;
    double EXY=0.0;
    double covXY=0.0, cofXY=0.0;
    struct FFSEncodeVec tmp_vec[1];

    if (size < 24) size = 24;

    if (attrs != NULL) {
	get_int_attr(attrs, CM_REBWM_RLEN, &N);
	
	get_int_attr(attrs, CM_REBWM_REPT, &repeat_count);
	CMtrace_out(conn->cm, CMTransportVerbose, "INT_CMregressive_probe_bandwidth: get from attr, N: %d, repeat_count: %d\n", N, repeat_count);
	if(N<6) N=6;
	if(repeat_count<3) repeat_count=3;
    } else {
	N=9;
	repeat_count=3;
    }
    CMtrace_out(conn->cm, CMConnectionVerbose, "CM - INITIATE BW MEASURE on CONN %p\n", conn);
    biggest_size=size*(N+1);

    if (max_block_size == 0) {
	char *new_block = malloc(biggest_size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = biggest_size;
	memset(block, 0xef, biggest_size);
    } else if (biggest_size > max_block_size) {
	char *new_block = realloc(block, biggest_size);
	if (new_block == NULL) return -1;
	block = new_block;
	max_block_size = biggest_size;
	memset(block, 0xef, biggest_size);
    }

    /* CMP\0 in first entry for CMPerformance message */
    ((int*)block)[0] = 0x434d5000;
    /* size in second entry, high byte gives CMPerf operation */

    dsm.MsgNum=N;
    dsm.AveRTDelay=malloc(sizeof(double)*N);
    dsm.MsgSize=malloc(sizeof(int)*N);
    
    for(i =0; i<N; i++){
	cond = INT_CMCondition_get(conn->cm, conn);
	((int*)block)[3] = cond;   /* condition value in third entry */
	INT_CMCondition_set_client_data( conn->cm, cond, &microsecs_to_receive);

	/* size in second entry, high byte gives CMPerf operation */
#if SIZEOF_LONG == 8	
	((int*)block)[1] = ((size>>32)&0xffffff) | (CMRegressivePerfBandwidthInit<<24);
	((int*)block)[2] = (size & 0xffffffff);
	((int*)block)[5] = size >> 32;
	((int*)block)[6] = size & 0xffffff;
#else
	((int*)block)[1] = 0 | (CMRegressivePerfBandwidthInit<<24);
	((int*)block)[2] = size;
	((int*)block)[5] = size;
	((int*)block)[6] = 0;
#endif

	CMtrace_out(conn->cm, CMTransportVerbose, "CM - Initiating bandwidth probe of %ld bytes, %d messages\n", size, repeat_count);
	tmp_vec[0].iov_base = &block[0];
	tmp_vec[0].iov_len = size;
	actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);

	if (actual != 1) {
	    return -1;
	}

	((int*)block)[1] = size | (CMRegressivePerfBandwidthBody<<24);
	for (j=0; j <(repeat_count-1); j++) {
	    actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);
	    if (actual != 1) {
		return -1;
	    }
	}
	((int*)block)[1] = size | (CMRegressivePerfBandwidthEnd <<24);
	actual = INT_CMwrite_raw(conn, tmp_vec, NULL, 1, size, NULL, 0);
	if (actual != 1) {
	    return -1;
	}

	if (INT_CMCondition_wait(conn->cm, cond) == 0) {
	    return 0.0;
	}
	bandwidth = ((double) size * (double)repeat_count * 1000.0) / 
	    (double)microsecs_to_receive;

	dsm.AveRTDelay[i]=(double)microsecs_to_receive*2.0/(double)repeat_count/1000000.0;
	dsm.MsgSize[i]=size;
	/*change size for the next round of bw measurment. */
	size+=biggest_size/(N+1);
	ave_delay+=dsm.AveRTDelay[i]*1000.0;
	ave_size+=dsm.MsgSize[i];
	EXY+=dsm.AveRTDelay[i]*1000.0*dsm.MsgSize[i];

	CMtrace_out(conn->cm, CMTransportVerbose, "CM - Partial Estimated bandwidth- %f Mbps, size: %ld, delay: %d, ave_delay+=%f\n", bandwidth*8.0 / 1000.0, size-biggest_size/(N+1), microsecs_to_receive, dsm.AveRTDelay[i]*1000.0);


    }
    bandwidth=Regression( conn, &dsm);

    ave_delay /= (double)N;
    ave_size /= (double)N;
    EXY /= (double)N;
    for(i=0;i<N; i++){
      var_delay += (dsm.AveRTDelay[i]*1000.0-ave_delay)*(dsm.AveRTDelay[i]*1000.0-ave_delay);
      var_size += (dsm.MsgSize[i]-ave_size)*(dsm.MsgSize[i]-ave_size);
    }
    var_delay /= (double)N;
    var_size /= (double)N;

    covXY=EXY-ave_delay*ave_size;
    cofXY=covXY/(sqrt(var_delay)*sqrt(var_size));
    
     CMtrace_out(conn->cm, CMTransportVerbose,"INT_CMregressive_probe_bandwidth: ave_delay: %f, ave_size: %f, var_delay: %f, var_size: %f, EXY: %f, covXY: %f, cofXY: %f\n", ave_delay, ave_size, var_delay, var_size, EXY, covXY, cofXY);
    
    CMtrace_out(conn->cm, CMTransportVerbose, "CM - Regressive Estimated bandwidth- %f Mbps, size: %ld\n", bandwidth, size);

    free(dsm.AveRTDelay);
    free(dsm.MsgSize);
  
    if(cofXY<0.97 && cofXY>-0.97)
	if(bandwidth>0) bandwidth*=-1; /*if the result is not reliable, return negative bandwidth*/
  
    if (conn->attrs == NULL) conn->attrs = CMcreate_attr_list(conn->cm);

    {
	int ibandwidth, icof;
	ibandwidth = (int) (bandwidth * 1000.0);
	icof = (int) (cofXY * 1000.0);

	CMtrace_out(conn->cm, CMTransportVerbose, "CM - Regressive Setting measures to BW %d kbps, COF %d\n", ibandwidth, icof);
	set_int_attr(conn->attrs, CM_BW_MEASURED_VALUE, ibandwidth);
	set_int_attr(conn->attrs, CM_BW_MEASURED_COF, icof);
     
    }
    return  bandwidth;

}

