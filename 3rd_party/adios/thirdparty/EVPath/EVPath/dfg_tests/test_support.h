#include "evpath.h"
#ifdef _MSC_VER
#define pid_t intptr_t
#include <process.h>
#endif

typedef struct _complex_rec {
    double r;
    double i;
} complex, *complex_ptr;

typedef struct _nested_rec {
    complex item;
} nested, *nested_ptr;

typedef struct _simple_rec {
    int integer_field;
    short short_field;
    long long_field;
    nested nested_field;
    double double_field;
    char char_field;
    int scan_sum;
} simple_rec, *simple_rec_ptr;

typedef struct _delay_struct {
    char **list;
    char *master_contact;
} delay_struct;


extern void(*on_exit_handler)();
extern FMStructDescRec simple_format_list[];
extern void generate_simple_record(simple_rec_ptr event);
extern int checksum_simple_record(simple_rec_ptr event, attr_list attrs, 
				  int quiet);
extern void test_fork_children(char **list, char *master_contact);
extern int wait_for_children(char **list);
extern int quiet;
extern int be_test_master(int argc, char **argv);
extern int be_test_child(int argc, char **argv);

extern void delayed_fork_children(CManager cm, char **list, char *master_contact, int delay_seconds);
