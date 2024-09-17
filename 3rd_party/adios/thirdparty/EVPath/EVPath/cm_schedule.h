
#ifndef __CM_SCHEDULE__H__
#define __CM_SCHEDULE__H__
/*! \file */

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef _MSC_VER
#define FD_SETSIZE 1024
#include <winsock2.h>
#endif
typedef struct _avail_period {
    struct timeval offset;
    struct timeval duration;
} *CMavail_period_ptr;

/*
 *  For network transports which support scheduling of RDMA pull requests,
 *  install a pull schedule.
 * \param cm The CManager in which to schedule pull messages.
 * \param start_time The starting time of a periodic schedule.
 * \param period The duration of a schedule period.
 * \param avail The offset and duration of times available for pulls.  
 *  This list is terminated with a {0, 0} value.
 */
extern int
CMinstall_pull_schedule(CManager cm, struct timeval *base_time, struct timeval *period, CMavail_period_ptr avail);


#ifdef	__cplusplus
}
#endif

#endif
