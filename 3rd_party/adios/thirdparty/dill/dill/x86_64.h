#ifndef _X86_64_H
#define _X86_64_H

#define BYTE_OUT1(c, insn1)                                                    \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 1;                              \
    } while (0)

#define BYTE_OUT2(c, insn1, insn2)                                             \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 2;                              \
    } while (0)

#define BYTE_OUT3(c, insn1, insn2, insn3)                                      \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 3;                              \
    } while (0)

#define BYTE_OUT4(c, insn1, insn2, insn3, insn4)                               \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        *(((unsigned char*)c->p->cur_ip) + 3) = (unsigned char)insn4;          \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 4;                              \
    } while (0)

#define BYTE_OUT5(c, insn1, insn2, insn3, insn4, insn5)                        \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        *(((unsigned char*)c->p->cur_ip) + 3) = (unsigned char)insn4;          \
        *(((unsigned char*)c->p->cur_ip) + 4) = (unsigned char)insn5;          \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 5;                              \
    } while (0)

#define BYTE_OUT6(c, insn1, insn2, insn3, insn4, insn5, insn6)                 \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        *(((unsigned char*)c->p->cur_ip) + 3) = (unsigned char)insn4;          \
        *(((unsigned char*)c->p->cur_ip) + 4) = (unsigned char)insn5;          \
        *(((unsigned char*)c->p->cur_ip) + 5) = (unsigned char)insn6;          \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 6;                              \
    } while (0)

#define BYTE_OUT3I(c, insn1, insn2, insn3, imm32)                              \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        *(unsigned int*)(((unsigned char*)c->p->cur_ip) + 3) =                 \
            (unsigned int)imm32;                                               \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 7;                              \
    } while (0)

#define BYTE_OUT3S(c, insn1, insn2, insn3, imm16)                              \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        *(((unsigned char*)c->p->cur_ip) + 2) = (unsigned char)insn3;          \
        *(unsigned short*)(((unsigned char*)c->p->cur_ip) + 3) =               \
            (unsigned short)imm16;                                             \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 5;                              \
    } while (0)

#define BYTE_OUT1I(c, insn1, imm32)                                            \
    do {                                                                       \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(unsigned int*)(((unsigned char*)c->p->cur_ip) + 1) =                 \
            (unsigned int)imm32;                                               \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 5;                              \
    } while (0)

#define BYTE_OUT2I(c, insn1, insn2, imm32)                                     \
    do {                                                                       \
        unsigned int tmp = (unsigned int)imm32;                                \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        memcpy((((unsigned char*)c->p->cur_ip) + 2), &tmp, 4);                 \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 6;                              \
    } while (0)

#define BYTE_OUT2II(c, insn1, insn2, imm32, imm32_2)                           \
    do {                                                                       \
        unsigned int tmp = (unsigned int)imm32;                                \
        unsigned int tmp2 = (unsigned int)imm32_2;                             \
        if (c->p->cur_ip >= c->p->code_limit) {                                \
            extend_dill_stream(c);                                             \
        }                                                                      \
        *(unsigned char*)c->p->cur_ip = (unsigned char)insn1;                  \
        *(((unsigned char*)c->p->cur_ip) + 1) = (unsigned char)insn2;          \
        memcpy((((unsigned char*)c->p->cur_ip) + 2), &tmp, 4);                 \
        memcpy((((unsigned char*)c->p->cur_ip) + 6), &tmp2, 4);                \
        if (c->dill_debug)                                                     \
            dump_cur_dill_insn(c);                                             \
        c->p->cur_ip = ((char*)c->p->cur_ip) + 10;                             \
    } while (0)

enum {
    AL = 0,
    CL,
    DL,
    BL,
    AH,
    CH,
    DH,
    BH, /* r8 */
    AX = 0,
    CX,
    DX,
    BX,
    SP,
    BP,
    SI,
    DI, /* r16 */
    EAX = 0,
    ECX,
    EDX,
    EBX,
    ESP,
    EBP,
    ESI,
    EDI, /* r32 */
    RAX = 0,
    RCX,
    RDX,
    RBX,
    RSP,
    RBP,
    RSI,
    RDI,
    R8,
    R9,
    R10,
    R11,
    R12,
    R13,
    R14,
    R15, /* r64 */
    XMM0 = 0,
    XMM1,
    XMM2,
    XMM3,
    XMM4,
    XMM5,
    XMM6,
    XMM7,
    XMM8,
    XMM9,
    XMM10,
    XMM11,
    XMM12,
    XMM13,
    XMM14,
    XMM15
};

#define REX_W 0x8
#define REX_R 0x4
#define REX_X 0x2
#define REX_B 0x1

typedef struct x86_64_mach_info {
    int act_rec_size;
    int stack_align;
    int stack_constant_offset;
    int conversion_word;
    int fcu_word;
    int save_base;
    int cur_arg_offset;
    int int_arg_count;
    int float_arg_offset;
    int float_arg_count;
    int backpatch_offset;
    int call_backpatch_offset;
    int call_stack_space;
    char pending_prefix;
    char* last_proc_ret_end;
    int varidiac_call;
    int non_var_args;
} * x86_64_mach_info;

extern int x86_64_type_align[];
extern int x86_64_type_size[];
extern void*
gen_x86_64_mach_info(dill_stream s);
extern void
x86_64_arith3(dill_stream c, int op, int commut, int dest, int src1, int src2);
extern void
x86_64_arith2(dill_stream c, int op, int subop, int dest, int src);
extern void
x86_64_mul(dill_stream c,
           int signed,
           int imm,
           int dest,
           int src1,
           IMM_TYPE src2);

extern void
x86_64_div_mod(dill_stream c, int sign, int div, int dest, int src1, int src2);
extern void
x86_64_div_modi(dill_stream c,
                int sign,
                int div,
                int dest,
                int src1,
                IMM_TYPE imm);
extern void
x86_64_arith3i(dill_stream c,
               int op,
               int commut,
               int dest,
               int src1,
               IMM_TYPE src2);
extern void
x86_64_shift(dill_stream c, int op, int junk, int dest, int src1, int src2);
extern void
x86_64_shifti(dill_stream c, int op, int junk, int dest, int src, IMM_TYPE imm);
extern void
x86_64_special(dill_stream c, special_operations type, IMM_TYPE param);
extern void
x86_64_set(dill_stream c, int r, IMM_TYPE imm);
extern void
x86_64_proc_start(dill_stream c,
                  char* subr_name,
                  int arg_count,
                  arg_info_list args,
                  dill_reg* arglist);
extern void
x86_64_end(dill_stream c);
extern void
x86_64_package_end(dill_stream c);
extern void*
x86_64_clone_code(dill_stream c, void* base, int size);
extern void
x86_64_ret(dill_stream c, int data1, int data2, int src);
extern void
x86_64_reti(dill_stream c, int data1, int data2, IMM_TYPE imm);
extern void
x86_64_retf(dill_stream c, int data1, int data2, double imm);
extern int
x86_64_getreg(dill_stream c, dill_reg* reg_p, int type, int class);
extern int
x86_64_putreg(dill_stream c, dill_reg reg, int type);
extern void
x86_64_ploadi(dill_stream c,
              int type,
              int junk,
              int dest,
              int src,
              IMM_TYPE offset);
extern void
x86_64_pload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
x86_64_pbsloadi(dill_stream c,
                int type,
                int junk,
                int dest,
                int src,
                IMM_TYPE offset);
extern void
x86_64_pbsload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
x86_64_pstorei(dill_stream c,
               int type,
               int junk,
               int dest,
               int src,
               IMM_TYPE offset);
extern void
x86_64_pstore(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
x86_64_modi(dill_stream c,
            int type,
            int junk,
            int dest,
            int src,
            IMM_TYPE offset);
extern void
x86_64_mod(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
x86_64_divi(dill_stream c,
            int type,
            int junk,
            int dest,
            int src,
            IMM_TYPE offset);
extern void
x86_64_div(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
x86_64_converti(dill_stream c,
                int from_type,
                int to_type,
                int dest,
                IMM_TYPE src);
extern void
x86_64_convert(dill_stream c, int from_type, int to_type, int dest, int src);
extern void
x86_64_mov(dill_stream c, int type, int junk, int dest, int src);
extern void
x86_64_pset(dill_stream c, int type, int junk, int dest, IMM_TYPE imm);
extern void
x86_64_setf(dill_stream c, int type, int junk, int dest, double imm);
extern void
x86_64_setp(dill_stream c, int type, int junk, int dest, void* imm);
extern void
x86_64_branch(dill_stream c, int op, int type, int src1, int src2, int label);
extern void
x86_64_branchi(dill_stream c,
               int op,
               int type,
               int src,
               IMM_TYPE imm,
               int label);
extern void
x86_64_compare(dill_stream s, int op, int type, int dest, int src1, int src2);
extern void
x86_64_comparei(dill_stream s,
                int op,
                int type,
                int dest,
                int src,
                IMM_TYPE imm);
extern void
x86_64_lea(dill_stream c, int junk, int junk1, int dest, int src, IMM_TYPE imm);
extern void
x86_64_farith(dill_stream c, int op, int typ, int dest, int src1, int src2);
extern void
x86_64_farith2(dill_stream c, int op, int typ, int dest, int src);
extern void
x86_64_bswap(dill_stream c, int op, int typ, int dest, int src);
extern void
x86_64_jump_to_label(dill_stream c, unsigned long label);
extern void
x86_64_jump_to_reg(dill_stream c, unsigned long reg);
extern void
x86_64_jump_to_imm(dill_stream c, void* imm);
extern void
x86_64_jal(dill_stream c, int return_addr_reg, int target);
extern int
x86_64_calli(dill_stream c, int type, void* xfer_address, const char* name);
extern int
x86_64_callr(dill_stream c, int type, int src);
extern void
x86_64_push(dill_stream c, int type, int reg);
extern void
x86_64_pushi(dill_stream c, int type, IMM_TYPE value);
extern void
x86_64_pushpi(dill_stream c, int type, void* value);
extern void
x86_64_pushfi(dill_stream c, int type, double value);
extern int
x86_64_local_op(dill_stream c, int flag, int val);
extern int
x86_64_local(dill_stream c, int type);
extern int
x86_64_localb(dill_stream c, int size);
extern void
x86_64_save_restore_op(dill_stream c, int save_restore, int type, int reg);
extern int
x86_64_init_disassembly_info(dill_stream c, void* ptr);
extern int
x86_64_print_insn(dill_stream c, void* info_ptr, void* insn);
extern int
x86_64_count_insn(dill_stream c, int start, int end);
extern void
x86_64_print_reg(dill_stream c, int typ, int reg);
#endif
