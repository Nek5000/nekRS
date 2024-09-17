#include <stdio.h>
#include <atl.h>
#include <evpath.h>
#undef NDEBUG
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cm_internal.h>
#include <cm_transport.h>

static transport_entry
create_thread_read_transport(CManager cm, transport_entry original) {
    if (original->data_available) {
        /* first check if we've already done this */
        transport_entry *cur = cm->transports;
        struct _transport_item new_entry;
        while (*cur) {
            if (*cur != original && !strcmp((*cur)->trans_name, original->trans_name)
                && !(*cur)->data_available) {
                return *cur;
            }
            ++cur;
        }

        /* otherwise */
        new_entry = *original;
        new_entry.data_available = NULL;
        return add_transport_to_cm(cm, &new_entry);
    } else {
        return original;
    }
}

static
void* read_thread_func(void *conn_raw) {
    CMConnection conn = conn_raw;
    transport_entry trans = conn->trans;
    while (!conn->closed && !conn->failed) {
        CMDataAvailable(trans, conn);
    }
    return NULL;
}

static thr_thread_t 
thr_fork(void*(*func)(void*), void *arg)
{
    thr_thread_t new_thread = 0;
    int err = thr_thread_create(&new_thread, NULL, (void*(*)(void*))func, arg);
    if (err != 0) {
	return  (thr_thread_t) (intptr_t) NULL;
    } else {
	return (thr_thread_t) new_thread;
    }
}

void
INT_CMstart_read_thread(CMConnection conn) {
    thr_thread_t thread;
    conn->trans = create_thread_read_transport(conn->cm, conn->trans);
    thread = thr_fork(read_thread_func, (void*)conn);
    assert(thread);
    thr_thread_detach(thread);
}


