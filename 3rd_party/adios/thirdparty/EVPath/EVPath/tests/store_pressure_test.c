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

typedef struct _complex_rec {
    double r;
    double i;
} complex, *complex_ptr;

typedef struct _nested_rec {
    complex item;
} nested, *nested_ptr;

static FMField nested_field_list[] =
{
    {"item", "complex", sizeof(complex), FMOffset(nested_ptr, item)},
    {NULL, NULL, 0, 0}
};

static FMField complex_field_list[] =
{
    {"r", "double", sizeof(double), FMOffset(complex_ptr, r)},
    {"i", "double", sizeof(double), FMOffset(complex_ptr, i)},
    {NULL, NULL, 0, 0}
};

static int seq_num = 0;
static int gen_seq_num = 0;

typedef struct _simple_rec {
    int sequence_number;
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;


static FMField simple_field_list[] =
{
    {"sequence_number", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, sequence_number)}, 
    {"integer_field", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, integer_field)},
    {"short_field", "integer",
     sizeof(short), FMOffset(simple_rec_ptr, short_field)},
    {"long_field", "integer",
     sizeof(long), FMOffset(simple_rec_ptr, long_field)},
    {"nested_field", "nested",
     sizeof(nested), FMOffset(simple_rec_ptr, nested_field)},
    {"double_field", "float",
     sizeof(double), FMOffset(simple_rec_ptr, double_field)},
    {"char_field", "char",
     sizeof(char), FMOffset(simple_rec_ptr, char_field)},
    {"scan_sum", "integer",
     sizeof(int), FMOffset(simple_rec_ptr, scan_sum)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec simple_format_list[] =
{
    {"simple", simple_field_list, sizeof(simple_rec), NULL},
    {"complex", complex_field_list, sizeof(complex), NULL},
    {"nested", nested_field_list, sizeof(nested), NULL},
    {NULL, NULL}
};

static
void 
generate_record(simple_rec_ptr event)
{
    long sum = 0;
    event->sequence_number = gen_seq_num++; 
    event->integer_field = (int) lrand48() % 100;
    sum += event->integer_field % 100;
    event->short_field = ((short) lrand48());
    sum += event->short_field % 100;
    event->long_field = ((long) lrand48());
    sum += event->long_field % 100;

    event->nested_field.item.r = drand48();
    sum += ((int) (event->nested_field.item.r * 100.0)) % 100;
    event->nested_field.item.i = drand48();
    sum += ((int) (event->nested_field.item.i * 100.0)) % 100;

    event->double_field = drand48();
    sum += ((int) (event->double_field * 100.0)) % 100;
    event->char_field = lrand48() % 128;
    sum += event->char_field;
    sum = sum % 100;
    event->scan_sum = (int) sum;
}

enum {
    CMD_FREEZE,
    CMD_UNFREEZE,
    CMD_EXPECT_EXACTLY
};

typedef struct _control_rec {
    int command;
    int argument;
} control_rec, *control_rec_ptr;

static FMField control_field_list[] =
{
    {"command", "integer",
     sizeof(int), FMOffset(control_rec_ptr, command)}, 
    {"argument", "integer",
     sizeof(int), FMOffset(control_rec_ptr, argument)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec control_format_list[] =
{
    {"contral", control_field_list},
    {NULL, NULL}
};




static int sink_stone = -1;

static int failures = 0;

static void fail(const char *format, ...) {
    va_list ap;
    va_start(ap, format);
    vprintf(format, ap);
    va_end(ap);
    printf("\n");
    ++failures;
}

static
int interpret_command(CManager cm, void *vevent, void *client_data, attr_list attrs)
{
    control_rec_ptr control = vevent;
    (void) client_data;
    (void) attrs;
    switch (control->command) {
    case CMD_FREEZE:
        EVfreeze_stone(cm, sink_stone);
        break;
    case CMD_UNFREEZE:
        EVunfreeze_stone(cm, sink_stone);
        break;
    case CMD_EXPECT_EXACTLY:
        if (seq_num != control->argument) {
            fail("Expected %d items recieved, found %d items\n", control->argument, seq_num);
        }
        break;
    default:
        fail("Unexpected control message (%d %d)\n", control->command, control->argument);
    }
    return 0;
}


int quiet = 1;

static
int dummy_handler(CManager cm, void *vevent, void *client_data, attr_list attrs) {
    (void) cm;
    (void) vevent;
    (void) client_data;
    (void) attrs;
    return 0;
}

static
int
simple_handler(CManager cm,void *vevent, void *client_data, attr_list attrs)
{
    simple_rec_ptr event = vevent;
    long sum = 0, scan_sum = 0;
    (void)cm;
    if (event->sequence_number != seq_num++) {
        printf("Sequence number %d but expected %d\n", event->sequence_number,
            seq_num - 1);
    }
    sum += event->integer_field % 100;
    sum += event->short_field % 100;
    sum += event->long_field % 100;
    sum += ((int) (event->nested_field.item.r * 100.0)) % 100;
    sum += ((int) (event->nested_field.item.i * 100.0)) % 100;
    sum += ((int) (event->double_field * 100.0)) % 100;
    sum += event->char_field;
    sum = sum % 100;
    scan_sum = event->scan_sum;
    if (sum != scan_sum) {
	printf("Received record checksum does not match. expected %d, got %d\n",
	       (int) sum, (int) scan_sum);
    }
    if ((quiet <= 0) || (sum != scan_sum)) {
	printf("In the handler, event data is :\n");
	printf("	integer_field = %d\n", event->integer_field);
	printf("	short_field = %d\n", event->short_field);
	printf("	long_field = %ld\n", event->long_field);
	printf("	double_field = %g\n", event->double_field);
	printf("	char_field = %c\n", event->char_field);
	printf("Data was received with attributes : \n");
	if (attrs) dump_attr_list(attrs);
    }
    if (client_data != NULL) {
	int tmp = *((int *) client_data);
	*((int *) client_data) = tmp + 1;
    }
    return 0;
}

static int repeat_count = 10000; 
static int use_remote = 0; /* have stall be on remote host? */
static int use_split = 1;  /* use a split stone? */
static int use_prefilter = 1; /* use a filter stone before any split stone? */
static int use_postfilter = 1; /* use a filter stone after the split stone? Irrelevant if no split stone. */
static int use_cycle = 1; /* use a non-accepting filter from a split stone to create a cycle. 
                            Irrelevant if no split stone. */
static int use_immed_filter = 1; /* use a filter immediate action (compiled ECL) rather than
                                    a callback from EVassoc_filter_action */
/* where to go when use_remote is set */
static EVstone output_stone; 
static EVstone command_output_stone;
static EVsource command_source;

static
void send_command(CManager cm, int type, int argument)
{
    control_rec cmd;
    (void)cm;
    cmd.command = type;
    cmd.argument = argument;
    EVsubmit(command_source, &cmd, NULL);
}

static EVstone setup_storage(CManager cm, EVstone sink) {
    EVstone st;
    EVaction act;

    st = EValloc_stone(cm);
    act = EVassoc_store_action(cm, st, sink, -1);

    assert(act == 0);

    printf("Storage stone is %d\n", st);

    return st;
}

enum { ACCEPT_FILTER, REJECT_FILTER };

static const char *filter_code[] = { "{return 1;}", "{return 0;}" };
static
int filter_accept(CManager cm, void *msg, void *client_data, attr_list attrs) {
    (void)cm;
    (void)msg;
    (void)client_data;
    (void)attrs;

    return 1;
}
static
int filter_reject(CManager cm, void *msg, void *client_data, attr_list attrs) {
    (void)cm;
    (void)msg;
    (void)client_data;
    (void)attrs;
    return 0;
}
static EVSimpleHandlerFunc filter_funcs[] = { filter_accept, filter_reject };

static EVstone make_dummy_filter(CManager cm, EVstone filter_stone, int filter_type, EVstone target) {
    EVaction filter_action;

    if (filter_stone == (EVstone) -1) {
        filter_stone = EValloc_stone(cm);
    }
    if (use_immed_filter) {
        char *filter;
        filter = create_filter_action_spec(simple_format_list, (char *) filter_code[filter_type]);
        filter_action = EVassoc_immediate_action(cm, filter_stone, filter, NULL);
        EVaction_set_output(cm, filter_stone, filter_action, 0, target);
    } else {
        filter_action = EVassoc_filter_action(cm, filter_stone, simple_format_list, 
                            filter_funcs[filter_type], target, NULL);
    }
    printf("A filter stone is %d\n", filter_stone);
    return filter_stone;
}

static void setup_sink(CManager cm, EVstone *ptarget, EVstone *psink) {
    EVstone other_sink;
    EVstone split_stone;
    EVstone prefilter = -1;

    *psink = use_remote ? output_stone :
        EVcreate_terminal_action(cm, simple_format_list, simple_handler, NULL);
    if (use_prefilter) {
        prefilter = EValloc_stone(cm); 
    }
    if (use_split) {
        EVstone target_list[4];
        split_stone = EValloc_stone(cm);
        other_sink = EVcreate_terminal_action(cm, simple_format_list, dummy_handler, NULL);
        target_list[0] = use_postfilter ? make_dummy_filter(cm, -1, ACCEPT_FILTER, *psink) :*psink;
        target_list[1] = other_sink;
        target_list[2] = -1;
        if (use_cycle) {
            target_list[2] = make_dummy_filter(cm, -1, REJECT_FILTER, 
                use_prefilter ? prefilter : split_stone);
            target_list[3] = -1;
        }
        EVassoc_split_action(cm, split_stone, target_list);
    }
    *ptarget = use_split ? split_stone : *psink;
    if (use_prefilter) {
        (void) make_dummy_filter(cm, prefilter, ACCEPT_FILTER, *ptarget);
        *ptarget = prefilter;
    }
}

static void fill_storage(CManager cm, EVstone storage, int count) {
    EVsource source_handle;
    source_handle = EVcreate_submit_handle(cm, storage, simple_format_list);
    while(count != 0) {
        simple_rec data;
        generate_record(&data);
        EVsubmit(source_handle, &data, NULL);
        count--;
    }
    EVfree_source(source_handle);
}

static int drain(CManager cm, EVstone storage) {
    int fail_count = 0;
    int last_count = EVstore_count(cm, storage, 0);
    printf("starting draining at %d\n", last_count);
    if (!EVstore_is_sending(cm, storage, 0)) {
        EVstore_start_send(cm, storage, 0);
    }
    for (;;) {
        int cur_count;
        cur_count = EVstore_count(cm, storage, 0);
        printf("(presleep) drained to %d\n", cur_count);
        CMsleep(cm, 1);
        cur_count = EVstore_count(cm, storage, 0);
        printf("drained to %d\n", cur_count);
        if (cur_count == last_count) {
            if (++fail_count == 4) {
                break;
            }
        } else {
            fail_count = 0;
        }
        last_count = cur_count;
    }
    if (last_count > 0 && !EVstore_is_sending(cm, storage, 0 )) {
        fail("Not sending despite being full");
    } else if (last_count == 0 && EVstore_is_sending(cm, storage, 0)) {
        fail("Still sending despite being empty");
    }
        
    return last_count;
}

static void
setup_contact(CManager cm) {
    EVstone contact;
    EVstone target = -1;
    contact = EValloc_stone(cm);
    printf("contact = %d\n", contact);
    EVassoc_terminal_action(cm, contact, control_format_list, interpret_command, NULL);
    assert(contact == 0);
    use_prefilter = use_postfilter = use_split = 0;
    setup_sink(cm, &sink_stone, &target);
    printf("target = %d\n", target);
    assert(target == 1);
}


static void run_test(CManager cm) {
    EVstone storage, sink, target;
    int last_count;

    setup_sink(cm, &target, &sink);
    storage = setup_storage(cm, target);

    if (use_remote)
        send_command(cm, CMD_FREEZE, 0);
    else
        EVfreeze_stone(cm, sink);
    
    fill_storage(cm, storage, repeat_count);
    last_count = drain(cm, storage);

    if (last_count == 0) {
        fail("No items left despite frozen sink");
    }

    if (use_remote)
        send_command(cm, CMD_UNFREEZE, 0);
    else
        EVunfreeze_stone(cm, sink);

    if (use_prefilter) {
        EVstall_stone(cm, target);
        last_count = drain(cm, storage);
        if (last_count == 0) {
            fail("No items left despite explicitly stalled filter");
        }
        EVunstall_stone(cm, target);
    }

    last_count = drain(cm, storage);

    if (last_count != 0) {
        fail("Didn't empty completely despite okay sink");
    }
}

static int
subprocess_work(const char *contact_list) {
    CManager cm;
    cm = CManager_create();
    command_output_stone = EValloc_stone(cm);
    EVassoc_bridge_action(cm, command_output_stone, 
        attr_list_from_string((char *) contact_list), 0);
    command_source = EVcreate_submit_handle(cm, command_output_stone, control_format_list);
    output_stone = EValloc_stone(cm);
    EVassoc_bridge_action(cm, output_stone,
        attr_list_from_string((char *) contact_list), 1);
    run_test(cm);
    EVfree_source(command_source);
    /* CManager_close(cm); */
    return failures;
}

static pid_t run_subprocess(char **args);
static pid_t subproc_proc = 0;
static void wait_for_subprocess(void);
static void do_remote_test(void) {
    CManager cm;
    cm = CManager_create();
    attr_list contact_list;
    char *subproc_args[] = { "store_pressure_test", "-r", NULL, NULL };
    setup_contact(cm);
    CMlisten(cm);
    /* CMfork_comm_thread(cm); */
    contact_list = CMget_contact_list(cm);
    assert(contact_list);
    subproc_args[2] = attr_list_to_string(contact_list);
    free_attr_list(contact_list);
    run_subprocess(subproc_args);
    CMsleep(cm, 20);
    wait_for_subprocess();
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

    /*
    for (i = 0; i < 2; ++i) {
        use_immed_filter = i;
        cm = CManager_create();
        run_test(cm);
        CManager_close(cm);
    }
    use_remote = 1;
    */
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
