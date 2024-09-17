
#define INSN_OUT(c, i)                                                         \
    if ((((char*)c->p->cur_ip) + sizeof(virtual_insn)) >=                      \
        (char*)c->p->code_limit) {                                             \
        extend_dill_stream(c);                                                 \
    }                                                                          \
    *(virtual_insn*)c->p->cur_ip = i;                                          \
    if (c->dill_debug) {                                                       \
        printf("%p  -- ", c->p->cur_ip);                                       \
        virtual_print_insn(c, NULL, c->p->cur_ip);                             \
        printf("\n");                                                          \
    }                                                                          \
    c->p->cur_ip = ((char*)c->p->cur_ip) + sizeof(virtual_insn)

typedef struct bitv {
    short len;
    char vec[2];
} * bit_vec;

typedef struct basic_block {
    ssize_t start;
    ssize_t end;
    int label;
    int end_branch_label;
    int fall_through;
    int pred_count;
    int* pred_list;
    int succ_count;
    int* succ_list;
    bit_vec regs_used;
    bit_vec regs_defined;
    bit_vec live_at_end;
    short* reg_assigns;
    int visited;
    int is_loop_start;
    int is_loop_end;
} * basic_block;

typedef struct mach_info {
    arg_info_list arg_info;
    int bbcount;
    basic_block bblist;
    int prefix_code_start;
} * virtual_mach_info;

extern int
dill_type_of(dill_stream c, size_t vreg);
extern void
virtual_proc_start(dill_stream c,
                   char* subr_name,
                   int arg_count,
                   arg_info_list args,
                   dill_reg* arglist);
extern void
virtual_end(dill_stream c);
extern void
virtual_package_end(dill_stream c);
extern int
virtual_print_insn(dill_stream c, void* info_ptr, void* i);
