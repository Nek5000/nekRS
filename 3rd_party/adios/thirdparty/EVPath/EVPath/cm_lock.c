#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "config.h"

#include "ffs.h"

#include "atl.h"
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"
#undef NDEBUG
#include "assert.h"

extern void
IntCManager_lock(CManager cm, const char *file, int line)
{
    CMtrace_out(cm, CMLowLevelVerbose, "CManager Lock at \"%s\" line %d\n",
		file, line);
    thr_mutex_lock(cm->exchange_lock);
    cm->locked++;
    if (cm->locked != 1) {
	printf("CManager lock inconsistency, %d\n", cm->locked);
    }
}

extern void
IntCManager_unlock(CManager cm, const char *file, int line)
{
    CMtrace_out(cm, CMLowLevelVerbose, "CManager Unlock at \"%s\" line %d\n",
		file, line);
    cm->locked--;
    if (cm->locked != 0) {
	printf("CManager unlock inconsistency, %d\n", cm->locked);
    }
    thr_mutex_unlock(cm->exchange_lock);
}

extern int
CManager_locked(CManager cm)
{
    return cm->locked;
}


