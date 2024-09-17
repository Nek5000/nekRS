#include "config.h"

#include <stdio.h>
#include <assert.h>
#include <stdarg.h>
#include <atl.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <arpa/inet.h>
#include "evpath.h"
#ifdef HAVE_WINDOWS_H
#include <windows.h>
#define drand48() (((double)rand())/((double)RAND_MAX))
#define lrand48() rand()
#define srand48(x)
#else
#include <sys/wait.h>
#endif

static int failures = 0;
static int deferred = 0;
static int quiet = 1;
static int use_remote = 0;
 
static void fail(const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vprintf(format, ap);
    va_end(ap);
    printf("\n");
    ++failures;
}

typedef struct send_request {
    int start_seqno;
    int size;
    int count;
    const char *target_contact;
    int target_stone;
} send_request, *send_request_ptr;

static FMField send_request_field_list[] = {
    {"start_seqno", "integer", sizeof(int), FMOffset(send_request_ptr, start_seqno)},
    {"size", "integer", sizeof(int), FMOffset(send_request_ptr, size)},
    {"count", "integer", sizeof(int), FMOffset(send_request_ptr, count)},
    {"target_contact", "string", sizeof(char*), FMOffset(send_request_ptr, target_contact)},
    {"target_stone", "integer", sizeof(int), FMOffset(send_request_ptr, target_stone)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec send_request_format_list[] = {
    {"send_request", send_request_field_list},
    {NULL, NULL}
};

typedef struct data_event {
    int seqno;
    char *data;
    int data_size;
} data_event, *data_event_ptr;

static FMField data_event_field_list[] = {
    {"seqno", "integer", sizeof(int), FMOffset(data_event_ptr, seqno)},
    {"data_size", "integer", sizeof(int), FMOffset(data_event_ptr, data_size)},
    {"data", "char[data_size]", sizeof(char), FMOffset(data_event_ptr, data)},
    {NULL, NULL}
};

static FMStructDescRec data_event_format_list[] = {
    {"data_event", data_event_field_list},
    {NULL, NULL}
};

typedef struct data_sender_state {
    EVstone my_stone;
    int active;
    int cur_seqno;
    int last_seqno;
    int cur_size;
    EVsource target;
    EVstone target_stone;
} ds_state, *ds_state_ptr;

static void send_some(CManager, ds_state_ptr);

#ifdef DEBUG
static void dump_request(send_request_ptr request) {
    fprintf(stderr, "from %d, %d items, size %d\n",
        request->start_seqno, request->count, request->size);
}
#endif

/* terminal handler for requesting that data be sent */
static int 
interpret_request(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    send_request_ptr request = vevent;
    (void) attrs;
    ds_state_ptr state = client_data;

#ifdef DEBUG
    dump_request(request);
#endif
    assert(!state->active);

    state->cur_seqno = request->start_seqno;
    state->last_seqno = request->start_seqno + request->count;
    state->cur_size = request->size;

#ifdef DEBUG
    fprintf(stderr, "about to freeze\n");
#endif
    EVfreeze_stone(cm, state->my_stone);
    state->active = 1;

    if (state->target_stone == -1) {
        state->target_stone = EValloc_stone(cm);
#ifdef DEBUG
        fprintf(stderr, "about to create output action for <%s>\n", request->target_contact);
#endif
        EVassoc_bridge_action(cm, state->target_stone,
            attr_list_from_string((char *) request->target_contact), request->target_stone);
        state->target = EVcreate_submit_handle(cm, state->target_stone,
            data_event_format_list);
    }

    send_some(cm, state);

    return 0;
}

static void
send_one_cb(CManager cm, EVstone ignored, void *state_raw) {
    ds_state_ptr state = state_raw; 
    (void) ignored;

    ++deferred;

    send_some(cm, state);
}
static void 
send_some(CManager cm, ds_state_ptr state) {
    char *data_data = calloc(1, state->cur_size);
    data_event data;
    int successp = 0;

#ifdef DEBUG
    fprintf(stderr, "sending starting with %d --> %d\n", state->cur_seqno, state->last_seqno);
#endif
    assert(state->active);
    while (state->cur_seqno != state->last_seqno) {
        data.data = data_data;
        data.data_size = state->cur_size;
        data.seqno = state->cur_seqno;
        successp = EVsubmit_or_wait(state->target, &data, NULL, send_one_cb, state);
        if (!successp) {
#ifdef DEBUG
            fprintf(stderr, "deferred at %d\n", state->cur_seqno);
#endif
            return;
        }
        state->cur_seqno++;
    } 
    state->active = 0;
#ifdef DEBUG
    fprintf(stderr, "about to unfreeze\n");
#endif
    EVunfreeze_stone(cm, state->my_stone); /* ready for next request */
}

static void
setup_contact(CManager cm) {
    EVstone contact;
    ds_state_ptr state = malloc(sizeof(ds_state));
    contact = EValloc_stone(cm);
    state->my_stone = contact;
    state->active = 0;
    state->target = NULL;
    state->target_stone = -1;
    assert(contact == 0);
    EVassoc_terminal_action(cm, contact, send_request_format_list, interpret_request, state);
}

static int message_count = 0;
static int last_seqno = 0;
static int
check_data(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    data_event_ptr data = vevent;
    (void) cm;
    (void) client_data;
    (void) attrs;

    if (data->seqno <= last_seqno || data->seqno % 1000 == 0) {
        fail("Bad sequence number");
    }
    last_seqno = data->seqno;
    ++message_count;

    return 0;
}

static int
subprocess_work(const char *contact_list) {
    CManager cm;
    EVstone check_stone, out_stone;
    EVsource source;
    int i;
    send_request req;
    cm = CManager_create();
    CMlisten(cm);
    check_stone = EValloc_stone(cm);
    EVassoc_terminal_action(cm, check_stone, data_event_format_list, check_data, NULL);
    out_stone = EValloc_stone(cm);
    EVassoc_bridge_action(cm, out_stone, attr_list_from_string((char *) contact_list), 0);
    source = EVcreate_submit_handle(cm, out_stone, send_request_format_list);
    req.target_contact = attr_list_to_string(CMget_contact_list(cm));
    req.target_stone = check_stone;
    req.size = 100000;
    for (i = 0; i < 20; ++i) {
        req.start_seqno = i * 1000 + 1;
        req.count = 999;
        EVsubmit(source, &req, NULL);
    }
    CMsleep(cm, 20);
    if (message_count != 999 * 20) {
        fail("Wrong number of messages (got %d, expected %d)", message_count, 999 * 20);
    }
    return failures;
}

static pid_t run_subprocess(char **args);
static pid_t subproc_proc = 0;
static void wait_for_subprocess(void);
static void do_remote_test(void) {
    CManager cm;
    attr_list contact_list;
    char *subproc_args[] = { "submit_wait_test", "-r", NULL, NULL };
    cm = CManager_create();
    setup_contact(cm);
    CMlisten(cm);
    /* CMfork_comm_thread(cm); */
    contact_list = CMget_contact_list(cm);
    assert(contact_list);
    subproc_args[2] = attr_list_to_string(contact_list);
    free_attr_list(contact_list);
    run_subprocess(subproc_args);
    CMsleep(cm, 30);
    wait_for_subprocess();
    if (deferred == 0) {
        fail("Didn't actually need to defer sending anything");
    }
    /* CManager_close(cm); */
}

char *transport = NULL;
#include "support.c"

int
main(int argc, char **argv)
{ 
    int remote_child = 0;

    PARSE_ARGS();

    srand48(getpid());
    if (remote_child) {
        use_remote = 1;
        return subprocess_work(argv[1]);
    }
    
    (void)use_remote;
    do_remote_test();
    return failures;
}

static void wait_for_subprocess(void) {
    int exit_state;
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

}
