/***** Includes *****/
#include "config.h"
#include <sys/types.h>

#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#include <windows.h>
#include <sys/timeb.h>
#define getpid()	_getpid()
#else
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#endif
#include <sys/socket.h>
#ifdef HAVE_SYS_UN_H
#include <sys/un.h>
#endif
#ifdef HAVE_SYS_UIO_H
#include <sys/uio.h>
#endif
#ifdef HAVE_HOSTLIB_H
#include "hostLib.h"
#endif
#ifdef HAVE_STREAMS_UN_H
#include <streams/un.h>
#endif
#include <netinet/in.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#endif
#include <stdio.h>
#include <fcntl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#define assert(e)  \
    ((void) ((e) ? 0 : __assert (#e, __FILE__, __LINE__)))
#define __assert(e, file, line) \
    ((void)printf ("%s:%u: failed assertion `%s'\n", file, line, e), abort())
#ifdef HAVE_MEMORY_H
#include <memory.h>
#endif

#include <atl.h>
#include "evpath.h"
#include "cm_transport.h"
#include "ev_select.h"
#include <pthread.h>
#include <sched.h>
#define thr_thread_t pthread_t
#define thr_thread_self() pthread_self()
#define thr_thread_yield() sched_yield()

#include <sys/epoll.h>
#define MAX_EVENTS 32

#ifndef SOCKET_ERROR
#define SOCKET_ERROR -1
#endif
#if defined (__INTEL_COMPILER)
#  pragma warning (disable: 869)
#  pragma warning (disable: 310)
#  pragma warning (disable: 1418)
#  pragma warning (disable: 180)
#  pragma warning (disable: 177)
#endif

typedef struct func_list_item {
    select_list_func func;
    void *arg1;
    void *arg2;
} FunctionListElement;

typedef struct select_data {
    thr_thread_t server_thread;
    int epfd;

    int 	sel_item_max;
    FunctionListElement *select_items;
    FunctionListElement *write_items;

    periodic_task_handle periodic_task_list;

    int		closed;
    CManager	cm;
    int 	select_consistency_number;
    int 	wake_read_fd;
    int 	wake_write_fd;
} *select_data_ptr;

static void wake_server_thread(select_data_ptr socket_data);
static void setup_wake_mechanism(CMtrans_services svcs,
				       select_data_ptr *sdp);
static int remove_periodic_task(select_data_ptr sd,
				      periodic_task_handle handle);

#ifdef HAVE_WINDOWS_H
/* Winsock init stuff  */
/* ask for ver 1.1 */
static WORD wVersionRequested = MAKEWORD(1, 1);
static WSADATA wsaData;
int nErrorStatus;
static char*WSAerror_str(int err);
#define FD_SETSIZE 1024
#endif

static void
init_select_data(CMtrans_services svc, select_data_ptr *sdp, CManager cm)
{
    select_data_ptr sd = malloc(sizeof(struct select_data));
    *sdp = sd;
    sd->epfd = epoll_create(1);
    sd->server_thread =  (thr_thread_t) NULL;
    sd->closed = 0;
    sd->sel_item_max = 0;
    sd->select_items = (FunctionListElement *) svc->malloc_func(sizeof(FunctionListElement));
    sd->select_items[0].func = NULL;
    sd->select_items[0].arg1 = NULL;
    sd->select_items[0].arg2 = NULL;
    sd->write_items = (FunctionListElement *) svc->malloc_func(sizeof(FunctionListElement));
    sd->write_items[0].func = NULL;
    sd->write_items[0].arg1 = NULL;
    sd->write_items[0].arg2 = NULL;
    
    sd->periodic_task_list = NULL;
    sd->select_consistency_number = 0;
    sd->wake_read_fd = -1;
    sd->wake_write_fd = -1;
    if (cm != NULL) {
	sd->cm = cm;
    }
    setup_wake_mechanism(svc, sdp);
}

	
typedef struct _periodic_task {
    int period_sec;
    int period_usec;
    thr_thread_t executing;
    struct timeval next_time;
    select_list_func func;
    void *arg1;
    void *arg2;
    periodic_task_handle next;
} task_handle_s;

static void
free_epoll_data(CMtrans_services svc, select_data_ptr *sdp)
{
    periodic_task_handle tasks;
    select_data_ptr sd = *sdp;
    *sdp = NULL;
    tasks = sd->periodic_task_list;
    svc->free_func(sd->select_items);
    svc->free_func(sd->write_items);
    while (tasks != NULL) {
	periodic_task_handle next = tasks->next;
	svc->free_func(tasks);
	tasks = next;
    }
    svc->free_func(sd);
}

#undef timercmp
#define	timercmp(tvp, uvp, cmp) \
	/* CSTYLED */ \
	(((tvp)->tv_sec cmp (uvp)->tv_sec) || \
	((((tvp)->tv_sec == (uvp)->tv_sec) && \
	/* CSTYLED */ \
	((tvp)->tv_usec cmp (uvp)->tv_usec))))

static void
set_soonest_timeout(struct timeval *timeout, periodic_task_handle task_list, struct timeval now)
{
    struct timeval this_delay;
    if (task_list == NULL) return;
    this_delay.tv_sec = task_list->next_time.tv_sec - now.tv_sec;
    this_delay.tv_usec = task_list->next_time.tv_usec - now.tv_usec;
    if (task_list->executing == (thr_thread_t)-1) {
	/* this task not executing already, see when it needs to run  */
	if (this_delay.tv_usec < 0) {
	    this_delay.tv_sec--;
	    this_delay.tv_usec += 1000000;
	}
	if (this_delay.tv_sec < 0) {
	    this_delay.tv_sec = this_delay.tv_usec = 0;
	}
	if ((timeout->tv_sec == -1) || (timercmp(&this_delay, timeout, <))) {
	    *timeout = this_delay;
	}
    }
    set_soonest_timeout(timeout, task_list->next, now);
}

static void
increment_time(struct timeval *time, int increment_sec, int increment_usec)
{
    time->tv_usec += increment_usec;
    time->tv_sec += increment_sec;
    if (time->tv_usec >= 1000000) {
	time->tv_sec += (time->tv_usec / 1000000);
	time->tv_usec = (time->tv_usec % 1000000);
    }
}

static void
shutdown_wake_mechanism(select_data_ptr sd);

static void
socket_select(CMtrans_services svc, select_data_ptr sd, int timeout_sec, int timeout_usec)
{
    int i, res;
    int fd;
    struct epoll_event events[MAX_EVENTS];
    int ep_timeout;

    struct timeval timeout;
    int tmp_select_consistency_number = sd->select_consistency_number;

    if (sd->closed) {
	sd->server_thread =  (thr_thread_t) NULL; 
	return;
    }

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    if (sd->server_thread ==  (thr_thread_t) NULL) {
	/* no server thread set, must be this one */
	sd->server_thread = thr_thread_self();
    }
    if (sd->server_thread != thr_thread_self()) {
	/* What?  We're polling, but we're not the server thread? */
	fprintf(stderr, "Warning:  Multiple threads calling CManager_socket_select.\n");
	fprintf(stderr, "          This situation may result in unexpected I/O blocking.\n");
	fprintf(stderr, "          Server thread set to %lx.\n", (long) thr_thread_self());
	sd->server_thread = thr_thread_self();
    }
    if ((timeout_sec >= 0) || (sd->periodic_task_list != NULL)) {
	struct timeval now;
#ifndef HAVE_WINDOWS_H
	gettimeofday(&now, NULL);
#else
	/* GSE...  No gettimeofday on windows.  
	 * Must use _ftime, get millisec time, convert to usec.  Bleh.
	 */
	struct _timeb nowb;
	_ftime(&nowb);
	now.tv_sec = nowb.time;
	now.tv_usec = nowb.millitm * 1000;
#endif
	if (timeout_usec >= 1000000) {
	    timeout_sec += (timeout_usec / 1000000);
	    timeout_usec = timeout_usec % 1000000;
	}
	timeout.tv_sec = timeout_sec;
	timeout.tv_usec = timeout_usec;
	
	set_soonest_timeout(&timeout, sd->periodic_task_list, now);
        svc->verbose(sd->cm, CMSelectVerbose, "CMSelect with timeout %d sec, %d usec", 
		       timeout.tv_sec, timeout.tv_usec);

	if (timeout.tv_sec == -1) {
	    timeout.tv_usec = 1;
	    timeout.tv_sec = 0;
	}
	DROP_CM_LOCK(svc, sd->cm);
	ep_timeout = (1000 * timeout.tv_sec) + (timeout.tv_usec / 1000);
        res = epoll_wait(sd->epfd, events, MAX_EVENTS, ep_timeout);

	ACQUIRE_CM_LOCK(svc, sd->cm);
    } else {
	svc->verbose(sd->cm, CMSelectVerbose, "CMSelect blocking select");
	DROP_CM_LOCK(svc, sd->cm);
	res = epoll_wait(sd->epfd, events, MAX_EVENTS, -1);

	ACQUIRE_CM_LOCK(svc, sd->cm);
    }
    if (sd->closed) {
	sd->server_thread =  (thr_thread_t) NULL; 
	return;
    }
#ifndef HAVE_WINDOWS_H
    if (res == -1) {
	if (errno == EINTR) {
	    return;
	}
	/* 
	 * if upon returning from select the consistency number
	 * has changed, the only safe thing to do is return so 
	 * that a new select can happen with consistent info.
	 */
	if (sd->select_consistency_number != 
	    tmp_select_consistency_number) return;
	if (errno == 0) {
	    /* 
	     * odd, but x86 Solaris seems to return -1 from select but not set
	     * the errno in some circumstances.  Just return when this happens.
	     */
	    return;
	}
	if (errno == EBADF) {
	    fprintf(stderr, "The epoll fd is invalid. This is catastrophic.\n");
	} else if (errno != EAGAIN) {
#ifdef HAVE_FDS_BITS
	    fprintf(stderr, "select failed, errno %d, rd_set was %lx, %lx,%lx, %lx\n\n", errno,
		    (long) ((fd_set *) sd->fdset)->fds_bits[0],
		    (long) ((fd_set *) sd->fdset)->fds_bits[1],
		    (long) ((fd_set *) sd->fdset)->fds_bits[2],
		    (long) ((fd_set *) sd->fdset)->fds_bits[3]);
	    fprintf(stderr, "timeout was %d and %d\n", timeout_sec, timeout_usec);
#else
	    fprintf(stderr, "select failed, errno %d\n", errno);
#endif
	    exit(1);
	    return;
	}
    }
#else
    /* 
     * Have to ignore the bad invalue for select on NT because it can't
     * do selects on files.  Otherwise you will get a bunch of
     * irritating select errors that aren't really errors 
     */
    if (res == SOCKET_ERROR) {
	int errno_val;
	errno_val = WSAGetLastError();
	if (errno_val == WSAEINTR || errno_val == WSAEINVAL) {
	    return;
	} else {
	    fprintf(stderr, "select failed, errno %d\n", 
		    WSAerror_str(errno_val));
	}
	return;
    }
#endif

#ifdef HAVE_FDS_BITS
    svc->verbose(sd->cm, CMSelectVerbose, "select returned, rd_set started %lx, %lx, %lx, %lx, result was %lx, %lx, %lx, %lx",
	    (long) ((fd_set *) sd->fdset)->fds_bits[0],
	    (long) ((fd_set *) sd->fdset)->fds_bits[1],
	    (long) ((fd_set *) sd->fdset)->fds_bits[2],
	    (long) ((fd_set *) sd->fdset)->fds_bits[3],
	    (long) ((fd_set *) &rd_set)->fds_bits[0],
	    (long) ((fd_set *) &rd_set)->fds_bits[1],
	    (long) ((fd_set *) &rd_set)->fds_bits[2],
	    (long) ((fd_set *) &rd_set)->fds_bits[3]);
    svc->verbose(sd->cm, CMSelectVerbose, "       write_set started %lx, %lx, %lx, %lx, result was %lx, %lx, %lx, %lx",
	    (long) ((fd_set *) sd->write_set)->fds_bits[0],
	    (long) ((fd_set *) sd->write_set)->fds_bits[1],
	    (long) ((fd_set *) sd->write_set)->fds_bits[2],
	    (long) ((fd_set *) sd->write_set)->fds_bits[3],
	    (long) ((fd_set *) &wr_set)->fds_bits[0],
	    (long) ((fd_set *) &wr_set)->fds_bits[1],
	    (long) ((fd_set *) &wr_set)->fds_bits[2],
	    (long) ((fd_set *) &wr_set)->fds_bits[3]);
#endif
    /* 
     * if upon returning from select the consistency number
     * has changed, the only safe thing to do is return so 
     * that a new select can happen with consistent info.
     */
    if (sd->select_consistency_number != 
	tmp_select_consistency_number) return;
    /* 
     * Careful!  We're reading the control list here without locking!
     * Something bad *might* happen, it's just unlikely.
     */
    for(i = 0; i < res; i++) {
    	if (sd->closed) {
    		sd->server_thread =  (thr_thread_t) NULL;
    		return;
    	}
    	fd = events[i].data.fd;
    	if(events[i].events & EPOLLIN) {
    		if(sd->select_items[fd].func != NULL) {
    			svc->verbose(sd->cm, CMSelectVerbose,
    				"Running select read action on fd %d", fd);
    			sd->select_items[fd].func(sd->select_items[fd].arg1,
    				sd->select_items[fd].arg2);
    		}
    	}
    	if (sd->select_consistency_number !=
    		tmp_select_consistency_number) return;
    	if(events[i].events & EPOLLOUT) {
    		if(sd->write_items[fd].func != NULL) {
    			svc->verbose(sd->cm, CMSelectVerbose,
    					"Running select write action on fd %d", fd);
    			sd->write_items[fd].func(sd->write_items[fd].arg1,
    				sd->write_items[fd].arg2);
    		} else {
    			fprintf(stderr, "FD %d is polled, but no write item function.\n", fd);
    		}
    	}
    	if (sd->select_consistency_number !=
    		tmp_select_consistency_number) return;
    }

    if (sd->periodic_task_list != NULL) {
	/* handle periodic tasks */
	periodic_task_handle this_periodic_task = sd->periodic_task_list;
	struct timeval now;

#ifndef HAVE_WINDOWS_H
	gettimeofday(&now, NULL);
#else
	/* GSE...  No gettimeofday on windows.  
	 * Must use _ftime, get millisec time, convert to usec.  Bleh.
	 */
	struct _timeb nowb;
	_ftime(&nowb);
	now.tv_sec = nowb.time;
	now.tv_usec = nowb.millitm * 1000;
#endif
	while (this_periodic_task != NULL ) {
	    periodic_task_handle next = this_periodic_task->next;
	    if (timercmp(&now, &this_periodic_task->next_time, >)) {
		increment_time(&this_periodic_task->next_time,
			       this_periodic_task->period_sec,
			       this_periodic_task->period_usec);
		if (this_periodic_task->executing == (thr_thread_t)-1) {
		    this_periodic_task->executing = thr_thread_self();
		    DROP_CM_LOCK(svc, sd->cm);
		    this_periodic_task->func(this_periodic_task->arg1,
					     this_periodic_task->arg2);
		    ACQUIRE_CM_LOCK(svc, sd->cm);
		    next = this_periodic_task->next;
		    this_periodic_task->executing = (thr_thread_t) -1;
		    if ((this_periodic_task->period_sec == 0) &&
			(this_periodic_task->period_usec == 0)) {
		        remove_periodic_task(sd, this_periodic_task);
		    }
		}
		if (sd->closed) {
		    shutdown_wake_mechanism(sd);
		    return;
		}
	    }
	    /* 
	     * if upon returning from a handler the consistency number
	     * has changed, the only safe thing to do is return so 
	     * that a new select can happen with consistent info.
	     */
	    if (sd->select_consistency_number != 
		tmp_select_consistency_number) return;
	    this_periodic_task = next;
	}
    }
    sd->select_consistency_number++;
}

extern void
libcmepoll_LTX_add_select(CMtrans_services svc, select_data_ptr *sdp, int fd, select_list_func func, void *arg1, void *arg2)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);
    struct epoll_event ep_event;

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)sdp, (CManager)NULL);
	sd = *((select_data_ptr *)sdp);
    }
    sd->select_consistency_number++;
    if (fd > sd->sel_item_max) {
	int i;
	int size = (fd+1)*sizeof(sd->select_items[0]);
	sd->write_items = 
	    (FunctionListElement *) svc->realloc_func(sd->write_items, size);
	sd->select_items = 
	    (FunctionListElement *) svc->realloc_func(sd->select_items, size);
	if ((sd->select_items == NULL) || (sd->write_items == NULL)) {
	    perror("Realloc failed\n");
	    exit(1);
	}

	for (i = sd->sel_item_max + 1; i <= fd; i++) {
	    sd->write_items[i].func = NULL;
	    sd->write_items[i].arg1 = NULL;
	    sd->write_items[i].arg2 = NULL;
	    sd->select_items[i].func = NULL;
	    sd->select_items[i].arg1 = NULL;
	    sd->select_items[i].arg2 = NULL;
	}
	sd->sel_item_max = fd;
    }
    memset(&ep_event, 0, sizeof(ep_event));
    ep_event.events = EPOLLIN;
    ep_event.data.fd = fd;
    if(epoll_ctl(sd->epfd, EPOLL_CTL_ADD, fd, &ep_event) < 0) {
    	if(errno == EEXIST) {
    		/* This is fd is already armed for read */
    		ep_event.events = EPOLLIN | EPOLLOUT;
    		if (epoll_ctl(sd->epfd, EPOLL_CTL_MOD, fd, &ep_event) < 0) {
    			fprintf(stderr, "Something bad in %s. %d\n", __func__, errno);
    		}
    	} else {
    		fprintf(stderr, "Something bad in %s. %d\n", __func__, errno);
    	}
    }

    svc->verbose(sd->cm, CMSelectVerbose, "Adding fd %d to select read list", fd);
    sd->select_items[fd].func = func;
    sd->select_items[fd].arg1 = arg1;
    sd->select_items[fd].arg2 = arg2;
    wake_server_thread(sd);
}

extern void
libcmepoll_LTX_write_select(CMtrans_services svc, select_data_ptr *sdp, int fd, select_list_func func, void *arg1, void *arg2)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);
    struct epoll_event ep_event;

    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)sdp, (CManager)NULL);
	sd = *((select_data_ptr *)sdp);
    }
    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    memset(&ep_event, 0, sizeof(ep_event));
    sd->select_consistency_number++;
    if (fd > sd->sel_item_max) {
	int i;
	int size = (fd+1)*sizeof(sd->select_items[0]);
	sd->select_items = 
	    (FunctionListElement *) svc->realloc_func(sd->select_items, size);
	sd->write_items = 
	    (FunctionListElement *) svc->realloc_func(sd->write_items, size);
	if ((sd->select_items == NULL) || (sd->write_items == NULL)) {
	    perror("Realloc failed\n");
	    exit(1);
	}

	for (i = sd->sel_item_max + 1; i <= fd; i++) {
	    sd->write_items[i].func = NULL;
	    sd->write_items[i].arg1 = NULL;
	    sd->write_items[i].arg2 = NULL;
	    sd->select_items[i].func = NULL;
	    sd->select_items[i].arg1 = NULL;
	    sd->select_items[i].arg2 = NULL;
	}
	sd->sel_item_max = fd;
    }
    ep_event.data.fd = fd;
    if(func != NULL) {
    	ep_event.events = EPOLLOUT;
    	if(epoll_ctl(sd->epfd, EPOLL_CTL_ADD, fd, &ep_event) < 0) {
    		if(errno == EEXIST) {
    		    /* This is fd is already armed for read */
    		    ep_event.events = EPOLLIN | EPOLLOUT;
    		    if (epoll_ctl(sd->epfd, EPOLL_CTL_MOD, fd, &ep_event) < 0) {
    		    	fprintf(stderr, "Something bad in %s. %d\n", __func__, errno);
    		    }
    		} else {
    		    fprintf(stderr, "Something bad in %s. %d\n", __func__, errno);
    		}
    	}
    } else if(sd->select_items[fd].func) {
    	/* This fd should stay armed for read */
    	ep_event.events = EPOLLIN;
    	if (epoll_ctl(sd->epfd, EPOLL_CTL_MOD, fd, &ep_event) < 0) {
    	    fprintf(stderr, "Something bad in %s. %d\n", __func__, errno);
    	}
	} else {
    	if(epoll_ctl(sd->epfd, EPOLL_CTL_DEL, fd, &ep_event) < 0) {
    	    fprintf(stderr, "Something bad happened in %s. %d\n", __func__, errno);
    	}
    }

    sd->write_items[fd].func = func;
    sd->write_items[fd].arg1 = arg1;
    sd->write_items[fd].arg2 = arg2;
    wake_server_thread(sd);
}

extern periodic_task_handle
libcmepoll_LTX_add_periodic(CMtrans_services svc, select_data_ptr *sdp, int interval_sec, int interval_usec, select_list_func func, void *arg1, void *arg2)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);
    periodic_task_handle handle = malloc(sizeof(struct _periodic_task));
    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)sdp, (CManager)NULL);
	sd = *((select_data_ptr *)sdp);
    }

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    handle->period_sec = interval_sec;
    handle->period_usec = interval_usec;
    handle->executing = (thr_thread_t) -1;
#ifndef HAVE_WINDOWS_H
    gettimeofday(&handle->next_time, NULL);
#else
	/* GSE...  No gettimeofday on windows.  
	 * Must use _ftime, get millisec time, convert to usec.  Bleh.
	 */
    {
	struct _timeb nowb;
	_ftime(&nowb);
	handle->next_time.tv_sec = nowb.time;
	handle->next_time.tv_usec = nowb.millitm * 1000;
    }
#endif
    increment_time(&handle->next_time, interval_sec, interval_usec);
    handle->func = func;
    handle->arg1 = arg1;
    handle->arg2 = arg2;
    handle->next = NULL;

    if (sd->periodic_task_list == NULL) {
	sd->periodic_task_list = handle;
    } else {
	handle->next = sd->periodic_task_list;
	sd->periodic_task_list = handle;
    }
    wake_server_thread(sd);
    return handle;
}


extern periodic_task_handle
libcmepoll_LTX_add_delayed_task(CMtrans_services svc, select_data_ptr *sdp, int delay_sec, int delay_usec, select_list_func func, void *arg1, void *arg2)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);
    periodic_task_handle handle = malloc(sizeof(struct _periodic_task));
    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)sdp, (CManager)NULL);
	sd = *((select_data_ptr *)sdp);
    }

    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    handle->period_sec = 0;
    handle->period_usec = 0;
    handle->executing = (thr_thread_t) -1;
#ifndef HAVE_WINDOWS_H
    gettimeofday(&handle->next_time, NULL);
#else
    {
	/* GSE...  No gettimeofday on windows.  
	 * Must use _ftime, get millisec time, convert to usec.  Bleh.
	 */
	struct _timeb nowb;
	_ftime(&nowb);
	handle->next_time.tv_sec = nowb.time;
	handle->next_time.tv_usec = nowb.millitm * 1000;
    }
#endif
    increment_time(&handle->next_time, delay_sec, delay_usec);
    handle->func = func;
    handle->arg1 = arg1;
    handle->arg2 = arg2;
    handle->next = NULL;

    if (sd->periodic_task_list == NULL) {
	sd->periodic_task_list = handle;
    } else {
	handle->next = sd->periodic_task_list;
	sd->periodic_task_list = handle;
    }
    wake_server_thread(sd);
    return handle;
}

static int
remove_periodic_task(select_data_ptr sd, periodic_task_handle handle)
{
    periodic_task_handle list, last = NULL;
    list = sd->periodic_task_list;
    
    while(list != handle) {
	last = list;
	list = list->next;
	if (list == NULL) {
	    return 0;
	}
    }
    /* unlink task */
    if (last == NULL) {
	sd->periodic_task_list = list->next;
    } else {
	last->next = list->next;
    }
    if (handle->executing != thr_thread_self()) {
	/* someone besides us executing this ? */
        int i = 0;
	while (handle->executing != (thr_thread_t)-1) {
	    /* wait until they're done */
	    thr_thread_yield();
	    i++;
	    if (i > 1000) {
	        /* give up */
	        continue;
	    }
	}
    }
    free(handle);
    sd->select_consistency_number++;
    return 1;
}


extern void
libcmepoll_LTX_remove_periodic(CMtrans_services svc, select_data_ptr *sdp, periodic_task_handle handle)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);
    if (sd == NULL) return;
    if (remove_periodic_task(sd, handle) == 0) {
	fprintf(stderr, "Periodic task not found for removal\n");
    }
}

extern void
libcmepoll_LTX_remove_select(CMtrans_services svc, select_data_ptr *sdp, int fd)
{
    select_data_ptr sd = *((select_data_ptr *)sdp);

    struct epoll_event ep_event = {0}; // for a dumb kernel bug that we will never see

    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)sdp, (CManager)NULL);
	sd = *((select_data_ptr *)sdp);
    }
    sd->select_consistency_number++;
    if(sd->write_items[fd].func) {
    	/* this fd should stay armed for write */
    	ep_event.data.fd = fd;
    	ep_event.events = EPOLLOUT;
    	if (epoll_ctl(sd->epfd, EPOLL_CTL_MOD, fd, &ep_event) < 0) {
    	    fprintf(stderr, "Something bad happened in %s. %d\n", __func__, errno);
    	}
    } else {
    	if(epoll_ctl(sd->epfd, EPOLL_CTL_DEL, fd, &ep_event) < 0) {
    		fprintf(stderr, "Something bad happened in %s. %d\n", __func__, errno);
        }
    }

    sd->select_items[fd].func = NULL;
    sd->select_items[fd].arg1 = NULL;
    sd->select_items[fd].arg2 = NULL;
    wake_server_thread(sd);
}

static void
shutdown_wake_mechanism(select_data_ptr sd)
{
    if (sd->wake_read_fd == -1) return;
    close(sd->wake_read_fd);
    close(sd->wake_write_fd);
    sd->wake_read_fd = sd->wake_write_fd = -1;
}

static void read_wake_fd(void *fd_as_ptr, void *junk)
{
    char buffer;
    int fd = (int) (long)fd_as_ptr;
#ifdef HAVE_WINDOWS_H
    recv(fd, &buffer, 1, 0);
#else
    if (read(fd, &buffer, 1) != 1) {
	perror("wake read failed\n");
    }
#endif
}

#ifdef HAVE_WINDOWS_H
static char*
WSAerror_str(int err)
{
    switch(err) {
    case WSAEINTR: return "WSAEINTR";
    case WSAEBADF: return "WSAEBADF";
    case WSAEACCES: return "WSAEACCES";
    case WSAEFAULT: return "WSAEFAULT";
    case WSAEINVAL: return "WSAEINVAL";
    case WSAEMFILE: return "WSAEMFILE";
    case WSAEWOULDBLOCK: return "WSAEWOULDBLOCK";
    case WSAEINPROGRESS: return "WSAEINPROGRESS";
    case WSAEALREADY: return "WSAEALREADY";
    case WSAENOTSOCK: return "WSAENOTSOCK";
    case WSAEDESTADDRREQ: return "WSAEDESTADDRREQ";
    case WSAEMSGSIZE: return "WSAEMSGSIZE";
    case WSAEPROTOTYPE: return "WSAEPROTOTYPE";
    case WSAENOPROTOOPT: return "WSAENOPROTOOPT";
    case WSAEPROTONOSUPPORT: return "WSAEPROTONOSUPPORT";
    case WSAESOCKTNOSUPPORT: return "WSAESOCKTNOSUPPORT";
    case WSAEOPNOTSUPP: return "WSAEOPNOTSUPP";
    case WSAEPFNOSUPPORT: return "WSAEPFNOSUPPORT";
    case WSAEAFNOSUPPORT: return "WSAEAFNOSUPPORT";
    case WSAEADDRINUSE: return "WSAEADDRINUSE";
    case WSAEADDRNOTAVAIL: return "WSAEADDRNOTAVAIL";
    case WSAENETDOWN: return "WSAENETDOWN";
    case WSAENETUNREACH: return "WSAENETUNREACH";
    case WSAENETRESET: return "WSAENETRESET";
    case WSAECONNABORTED: return "WSAECONNABORTED";
    case WSAECONNRESET: return "WSAECONNRESET";
    case WSAENOBUFS: return "WSAENOBUFS";
    case WSAEISCONN: return "WSAEISCONN";
    case WSAENOTCONN: return "WSAENOTCONN";
    case WSAESHUTDOWN: return "WSAESHUTDOWN";
    case WSAETOOMANYREFS: return "WSAETOOMANYREFS";
    case WSAETIMEDOUT: return "WSAETIMEDOUT";
    case WSAECONNREFUSED: return "WSAECONNREFUSED";
    case WSAELOOP: return "WSAELOOP";
    case WSAENAMETOOLONG: return "WSAENAMETOOLONG";
    case WSAEHOSTDOWN: return "WSAEHOSTDOWN";
    case WSAEHOSTUNREACH: return "WSAEHOSTUNREACH";
    case WSAENOTEMPTY: return "WSAENOTEMPTY";
    case WSAEPROCLIM: return "WSAEPROCLIM";
    case WSAEUSERS: return "WSAEUSERS";
    case WSAEDQUOT: return "WSAEDQUOT";
    case WSAESTALE: return "WSAESTALE";
    case WSAEREMOTE: return "WSAEREMOTE";
    case WSAEDISCON: return "WSAEDISCON";
    case WSASYSNOTREADY: return "WSASYSNOTREADY";
    case WSAVERNOTSUPPORTED: return "WSAVERNOTSUPPORTED";
    case WSANOTINITIALISED: return "WSANOTINITIALISED";
    default: return "Unknown Winsock error";
    }
}
/*
 *  Note.  Unfortunately, the _pipe() function on WinNT 
 *  produces FDs that you can't use in select().  This ruins what we want
 *  this pipe for, which is to wake up a thread sleeping in select().
 *  So, we need to introduce a pipe function that returns two socket FDs.
 *  NT Sux.
 */

static int
pipe(int filedes[2])
{
    
    int length;
    struct sockaddr_in sock_addr;
    int sock_opt_val = 1;
    int sock1, sock2, conn_sock;
    unsigned long block = TRUE;
    int delay_value = 1;
   
    conn_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (conn_sock == SOCKET_ERROR) {
	fprintf(stderr, "Cannot open INET socket\n");
	return -1;
    }
    sock_addr.sin_family = PF_INET;
    sock_addr.sin_addr.s_addr = INADDR_ANY;
    sock_addr.sin_port = 0;
    if (bind(conn_sock, (struct sockaddr *) &sock_addr,
	     sizeof sock_addr) == SOCKET_ERROR) {
	fprintf(stderr, "Cannot bind INET socket\n");
	return -1;
    }
    length = sizeof sock_addr;
    if (getsockname(conn_sock, (struct sockaddr *) &sock_addr, &length) < 0) {
	fprintf(stderr, "Cannot get socket name\n");
	return -1;
    }
    /* begin listening for conns */
    if (listen(conn_sock, FD_SETSIZE)) {
	fprintf(stderr, "listen failed\n");
	return -1;
    }

/* send sock */
    if ((sock1 = socket(AF_INET, SOCK_STREAM, 0)) == SOCKET_ERROR) {
	return -1;
    }
    sock_addr.sin_addr.s_addr = 0x0100007f;  /* loopback */
    sock_addr.sin_family = PF_INET;
    if (ioctlsocket(sock1, FIONBIO, &block) != 0) {
	printf("ioctl failed\n");
    }
    if (connect(sock1, (struct sockaddr *) &sock_addr,
		sizeof sock_addr) == SOCKET_ERROR) {
	int err = WSAGetLastError();
	if (err != WSAEWOULDBLOCK) {
	    printf("unexpected error from connect, %s\n", WSAerror_str(err));
	}
    }

    if ((sock2 = accept(conn_sock, (struct sockaddr *) 0, (int *) 0)) == SOCKET_ERROR) {
	    int err = WSAGetLastError();
	    printf("err was %s\n", WSAerror_str(err));
    }
    
    setsockopt(sock2, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));
    {
	fd_set stXcptFDS,stWriteFDS;
	struct timeval stTimeOut;	/* for select() timeout (none) */
	int wRet;

	EVPATH_FD_ZERO((fd_set FAR*)&(stXcptFDS));
	EVPATH_FD_ZERO((fd_set FAR*)&(stWriteFDS));
	FD_SET(sock1, (fd_set FAR*)&(stWriteFDS));
	FD_SET(sock1, (fd_set FAR*)&(stXcptFDS));
	stTimeOut.tv_sec  = 10;
	stTimeOut.tv_usec = 0;
	wRet = select(-1, NULL, 
		      (fd_set FAR*)&(stWriteFDS),
		      (fd_set FAR*)&(stXcptFDS), 
		      NULL);
	if (wRet == SOCKET_ERROR) {
	    int err = WSAGetLastError();
	    printf("err was %s\n", WSAerror_str(err));
	}
    }
    setsockopt(sock1, IPPROTO_TCP, TCP_NODELAY, (char *) &delay_value,
	       sizeof(delay_value));

    filedes[0] = sock1;
    filedes[1] = sock2;
    return 0;
}
#endif

static void
setup_wake_mechanism(CMtrans_services svc, select_data_ptr *sdp)
{
    int filedes[2];

    select_data_ptr sd = *((select_data_ptr *)sdp);
    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    if (sd->wake_read_fd != -1) return;
    if (pipe(filedes) != 0) {
	perror("Pipe for wake not created.  Wake mechanism inoperative.");
	return;
    }
    sd->wake_read_fd = filedes[0];
    sd->wake_write_fd = filedes[1];
    svc->verbose(sd->cm, CMSelectVerbose, "CMSelect Adding read_wake_fd as action on fd %d",
		   sd->wake_read_fd);
    libcmepoll_LTX_add_select(svc, sdp, sd->wake_read_fd, read_wake_fd, 
			       (void*)(long)sd->wake_read_fd, NULL);
}

extern void
libcmepoll_LTX_wake_function(CMtrans_services svc, select_data_ptr *sdp)
{
    if (*sdp != NULL) {
	wake_server_thread(*sdp);
    }
}

static void
wake_server_thread(select_data_ptr sd)
{
    static char buffer = 'W';  /* doesn't matter what we write */
    if (sd->wake_write_fd != -1) {
#ifdef HAVE_WINDOWS_H
	send(sd->wake_write_fd, &buffer, 1, 0);
#else
	if (write(sd->wake_write_fd, &buffer, 1) != 1) {
	    printf("Whoops, wake write failed\n");
	}
#endif
    }
}

extern void
libcmepoll_LTX_blocking_function(CMtrans_services svc, void *client_data)
{
    select_data_ptr sd = *((select_data_ptr *)client_data);
    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)client_data, (CManager)NULL);
	sd = *((select_data_ptr *)client_data);
    }
    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    socket_select(svc, sd, -1, 0); /* no timeout */
}

extern void
libcmepoll_LTX_polling_function(CMtrans_services svc, void *client_data)
{
    select_data_ptr sd = *((select_data_ptr *)client_data);
    if (sd == NULL) {
	init_select_data(svc, (select_data_ptr*)client_data, (CManager)NULL);
	sd = *((select_data_ptr *)client_data);
    }
    if (sd->cm) {
	/* assert CM is locked */
	assert(CM_LOCKED(svc, sd->cm));
    }
    socket_select(svc, sd, 0, 0); /* no blocking, immediate timeout */
}

extern void
libcmepoll_LTX_select_initialize(CMtrans_services svc, CManager cm, void *client_data)
{
    if (*((select_data_ptr *)client_data) == NULL) {
	init_select_data(svc, (select_data_ptr*)client_data, cm);
    }
}

extern void
libcmepoll_LTX_select_shutdown(CMtrans_services svc, CManager cm, void *client_data)
{
    select_data_ptr *sdp = client_data;
    select_data_ptr sd = *sdp;

    svc->verbose(sd->cm, CMSelectVerbose, "CMSelect Shutdown task called");
    if (sd->server_thread != thr_thread_self()) {
	sd->closed = 1;
	close(sd->epfd);
	wake_server_thread(sd);
    }
}

extern void
libcmepoll_LTX_select_free(CMtrans_services svc, CManager cm, void *client_data)
{
    select_data_ptr *sdp = client_data;
    select_data_ptr sd = *sdp;

    svc->verbose(sd->cm, CMFreeVerbose, "CMSelect free task called");

    if (*((select_data_ptr *)client_data) != NULL) {
	free_epoll_data(svc, sdp);
    }
}

extern void
libcmepoll_LTX_select_stop(CMtrans_services svc, void *client_data)
{
    if (*((select_data_ptr *)client_data) != NULL) {
	(*((select_data_ptr*)client_data))->closed = 1;
    }
}

extern void
libcmepoll_init_sel_item(struct _select_item *sel_item)
{
    sel_item->add_select = (CMAddSelectFunc)libcmepoll_LTX_add_select;
    sel_item->remove_select = (CMRemoveSelectFunc)libcmepoll_LTX_remove_select;
    sel_item->write_select = (CMAddSelectFunc) libcmepoll_LTX_write_select;
    sel_item->add_periodic = (CMAddPeriodicFunc)libcmepoll_LTX_add_periodic;
    sel_item->add_delayed_task = 
	 (CMAddPeriodicFunc)libcmepoll_LTX_add_delayed_task;
    sel_item->remove_periodic = (CMRemovePeriodicFunc)libcmepoll_LTX_remove_periodic;
    sel_item->wake_function = (CMWakeSelectFunc)libcmepoll_LTX_wake_function;
    sel_item->blocking_function = (CMPollFunc)libcmepoll_LTX_blocking_function;
    sel_item->polling_function =  (CMPollFunc)libcmepoll_LTX_polling_function;
    sel_item->initialize = (SelectInitFunc)libcmepoll_LTX_select_initialize;
    sel_item->shutdown = (SelectInitFunc) libcmepoll_LTX_select_shutdown;
    sel_item->free = (SelectInitFunc) libcmepoll_LTX_select_free;
    sel_item->stop = (CMWakeSelectFunc) libcmepoll_LTX_select_stop;
}
