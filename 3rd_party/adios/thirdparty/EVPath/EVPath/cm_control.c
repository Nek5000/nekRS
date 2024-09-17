#include "config.h"
#include <stdio.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <ffs.h>
#include <atl.h>
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"

typedef struct _CMCondition {
    CMCondition next;
    int condition_num;
    int waiting;
    int signaled;
    int failed;
    thr_condition_t cond_condition;
    CMConnection conn;
    void *client_data;
} CMCondition_s;

static int cm_control_debug_flag = -1;

static void set_debug_flag(CManager cm)
{
    (void)cm;
    if (cm_control_debug_flag == -1) {
	if (CMtrace_on(cm, CMLowLevelVerbose)) {
	    cm_control_debug_flag = 1;
	} else {
	    cm_control_debug_flag = 0;
	}
    }
}

static CMCondition
CMCondition_find(CMControlList cl, int condition)
{
    CMCondition next = cl->condition_list;
    while (next != NULL) {
	if (next->condition_num == condition) {
	    return next;
	}
	next = next->next;
    }
    fprintf(stderr, "Serious internal error.  Use of condition %d, no longer in control list\n", condition);
    return NULL;
}


extern int
INT_CMCondition_get(CManager cm, CMConnection conn)
{
    CMControlList cl = cm->control_list;
    CMCondition cond = INT_CMmalloc(sizeof(CMCondition_s));
    set_debug_flag(cm);
    cond->next = cl->condition_list;
    cl->condition_list = cond;
    cond->condition_num = cl->next_condition_num++;
    cond->conn = conn;
    if (cl->next_condition_num >= 0xffffff) {
	/* recycle at  (16 M - 1) [ Caution on number reuse ] */
	cl->next_condition_num = 0;
    }
    cond->waiting = 0;
    cond->signaled = 0;
    cond->failed = 0;
    if (conn && conn->closed) {
	cond->failed = 1;
    }
    thr_condition_init(cond->cond_condition);
    return cond->condition_num;
}

static void
CMCondition_trigger(CManager cm, CMCondition cond, CMControlList cl)
{
    (void)cl;
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel Triggering CMcondition %d\n", cond->condition_num);
    }
    if (cond->waiting) {
	if (cm_control_debug_flag) {
	    fprintf(cm->CMTrace_file, "CMLowLevel Triggering CMcondition %d\n", 
		   cond->condition_num);
	}
	thr_condition_signal(cond->cond_condition);
    }
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel After trigger for CMcondition %d\n", 
	       cond->condition_num);
    }
}

extern void
CMconn_fail_conditions(CMConnection conn)
{
    CMControlList cl = conn->cm->control_list;
    CMCondition cond_list;
    set_debug_flag(conn->cm);
    cond_list = cl->condition_list;
    while(cond_list != NULL) {
	if (cond_list->conn == conn) {
	    cond_list->failed = 1;
	    CMCondition_trigger(conn->cm, cond_list, cl);
	}
	cond_list = cond_list->next;
    }
    if (cl->cond_polling) {
	/* wake the server thread in case we're not him */
	CMwake_server_thread(conn->cm);
    }
}

extern void
INT_CMCondition_fail(CManager cm, int condition)
{
    CMControlList cl = cm->control_list;
    CMCondition cond = CMCondition_find(cl, condition);
    if (!cond) return;
    cond->failed = 1;
    CMCondition_trigger(cm, cond, cm->control_list);

    if (cl->cond_polling) {
	/* wake the server thread in case we're not him */
	CMwake_server_thread(cm);
    }
}

void
CMCondition_destroy(CMControlList cl, int condition)
{
    CMCondition cond = NULL, prev = NULL;
    CMCondition next = NULL;

    if (cl->condition_list) {
	if (cl->condition_list->condition_num == condition) {
	    cond = cl->condition_list;
	    cl->condition_list = cl->condition_list->next;
	} else {
	    prev = cl->condition_list;
	    next = cl->condition_list->next;
	    while (next != NULL) {
		if (next->condition_num == condition) {
		    cond = next;
		    prev->next = next->next;
		    break;
		}
		prev = next;
		next = next->next;
	    }
	}
    }
    if (cond == NULL) {
	fprintf(stderr, "Serious internal error.  Use of condition %d, no longer in control list\n", condition);
    } else {
	/* free internal elements */
	thr_condition_free(cond->cond_condition);
	INT_CMfree(cond);
    }
}

extern void
internal_condition_free(CMControlList cl)
{
    while(cl->condition_list) {
        CMCondition_destroy(cl, cl->condition_list->condition_num);
    }
}

extern int
INT_CMCondition_has_signaled(CManager cm, int condition)
{
    int retval;
    CMCondition cond;
    CMControlList cl = cm->control_list;
    set_debug_flag(cm);

    cond = CMCondition_find(cl, condition);
    if (!cond) return -1;
    retval = cond->signaled;
    
    return retval;
}

extern int
INT_CMCondition_has_failed(CManager cm, int condition)
{
    int retval;
    CMCondition cond;
    CMControlList cl = cm->control_list;
    set_debug_flag(cm);

    cond = CMCondition_find(cl, condition);
    if (!cond) return -1;
    retval = cond->failed;
    
    return retval;
}

extern int
INT_CMCondition_wait(CManager cm, int condition)
{
    CMCondition cond;
    CMControlList cl = cm->control_list;
    int result;

    assert(CManager_locked(cm));
    set_debug_flag(cm);
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel Waiting for CMcondition %d\n", condition);
    }
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel locked cl\n");
    }
    cond = CMCondition_find(cl, condition);

    if (cond == NULL) return -1;
    if (cond->signaled) {
	if (cm_control_debug_flag) {
	    fprintf(cm->CMTrace_file, "CMcondition %d already signalled\n", condition);
	}
	return 1;
    }
    if (cond->failed) {
	if (cm_control_debug_flag) {
	    fprintf(cm->CMTrace_file, "CMcondition %d already failed\n", condition);
	}
	return 0;
    }
    cond->waiting++;
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel In condition wait, server thread = %p\n", 
	       (void*)(intptr_t)cl->server_thread);
    }
    if (!cl->has_thread) {
	if ((cl->server_thread ==  (thr_thread_t) (intptr_t) NULL) || (cl->server_thread == thr_thread_self())) {
	    cl->cond_polling = 1;
	    while (!(cond->signaled || cond->failed)) {
		if (cm_control_debug_flag) {
		    fprintf(cm->CMTrace_file, "CMLowLevel  Polling for CMcondition %d\n", condition);
		}
		CMcontrol_list_wait(cl);
	    }
	    cl->cond_polling = 0;
	    if (cm_control_debug_flag) {
		fprintf(cm->CMTrace_file, "CMLowLevel  after Polling for CMcondition %d\n", condition);
	    }
	    /* the poll and handle will set cl->server_thread, restore it */
	    cl->server_thread =  (thr_thread_t) (intptr_t)NULL;
	    if (cm_control_debug_flag) {
		fprintf(cm->CMTrace_file, "CMLowLevel  In condition wait, reset server thread = %lx\n", 
		       (long)cl->server_thread);
	    }
	} else {
	    /* some other thread is servicing the network here 
	       hopefully they'll keep doing it */
	    /* some other thread is the server thread */
	    if (cm_control_debug_flag) {
		fprintf(cm->CMTrace_file, "CMLowLevel Waiting for CMcondition %d\n", 
		       condition);
	    }
	    assert(CManager_locked(cm));
	    cm->locked--;
	    thr_condition_wait(cond->cond_condition, cm->exchange_lock);
	    cm->locked++;
	    if (cm_control_debug_flag) {
		fprintf(cm->CMTrace_file, "CMLowLevel After wait for CMcondition %d\n", 
		       condition);
	    }
	}
    } else if (thr_thread_self() == cl->server_thread) {
	/* we're the server thread */
	cl->cond_polling = 1;
	while (!(cond->signaled || cond->failed)) {
	    if (cm_control_debug_flag) {
		fprintf(cm->CMTrace_file, "CMLowLevel polling for CMcondition %d\n", condition);
	    }
	    CMcontrol_list_wait(cl);
	    if (cl->closed) cond->failed = 1;
	}
	cl->cond_polling = 0;
    } else {
	/* some other thread is the server thread */
	if (cm_control_debug_flag) {
	    fprintf(cm->CMTrace_file, "CMLowLevel Waiting for CMcondition %d\n", 
		   condition);
	}
	assert(CManager_locked(cm));
	cm->locked--;
	thr_condition_wait(cond->cond_condition, cm->exchange_lock);
	cm->locked++;
	if (cm_control_debug_flag) {
	    fprintf(cm->CMTrace_file, "CMLowLevel After wait for CMcondition %d\n", 
		   condition);
	}
    }
    result = cond->signaled;
    CMCondition_destroy(cl, condition);
    if (cm_control_debug_flag) {
	fprintf(cm->CMTrace_file, "CMLowLevel Return from wait CMcondition %d\n", condition);
    }
    return result;
}

extern void
INT_CMCondition_signal(CManager cm, int condition)
{
    CMCondition cond;
    CMControlList cl = cm->control_list;
    if(!CManager_locked(cm)) {
	printf("Not LOCKED!\n");
    }
    set_debug_flag(cm);
    cond = CMCondition_find(cl, condition);
    if (!cond) return;
    cond->signaled = 1;
    CMCondition_trigger(cm, cond, cl);
    if (cl->has_thread == 0) cm->abort_read_ahead = 1;
    if (cl->cond_polling) {
	/* wake the server thread in case we're not him */
	CMwake_server_thread(cm);
    }
}

extern void
INT_CMCondition_set_client_data(CManager cm, int condition, void *client_data)
{
    CMCondition cond;
    CMControlList cl = cm->control_list;
    set_debug_flag(cm);
    cond = CMCondition_find(cl, condition);
    if (!cond) return;
    cond->client_data = client_data;
}

extern void *
INT_CMCondition_get_client_data(CManager cm, int condition)
{
    CMCondition cond;
    void *client_data;
    CMControlList cl = cm->control_list;
    set_debug_flag(cm);
    cond = CMCondition_find(cl, condition);
    if (cond) {
	client_data = cond->client_data;
	return client_data;
    } else {
	return NULL;
    }
}
