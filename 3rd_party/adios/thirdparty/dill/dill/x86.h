#ifndef _X86_H
#define _X86_H

extern void x86_FORM3_farith(dill_stream s, int op3, int op, int dest, int src1, int src2);
extern void x86_FORM2_farith(dill_stream s, int op3, int op, int dest, int src);

#define BYTE_OUT1(c, insn1) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+1;	 \
} while (0)

#define BYTE_OUT2(c, insn1, insn2) \
do {\
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+2;\
} while (0)

#define BYTE_OUT3(c, insn1, insn2, insn3) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+3;\
} while (0)

#define BYTE_OUT4(c, insn1, insn2, insn3, insn4) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(((unsigned char*)c->p->cur_ip) + 3)= (unsigned char)insn4;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+4;\
} while (0)

#define BYTE_OUT5(c, insn1, insn2, insn3, insn4, insn5) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(((unsigned char*)c->p->cur_ip) + 3)= (unsigned char)insn4;\
*(((unsigned char*)c->p->cur_ip) + 4)= (unsigned char)insn5;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+5;\
} while (0)

#define BYTE_OUT6(c, insn1, insn2, insn3, insn4, insn5, insn6) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(((unsigned char*)c->p->cur_ip) + 3)= (unsigned char)insn4;\
*(((unsigned char*)c->p->cur_ip) + 4)= (unsigned char)insn5;\
*(((unsigned char*)c->p->cur_ip) + 5)= (unsigned char)insn6;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+6;\
} while (0)

#define BYTE_OUT3I(c, insn1, insn2, insn3, imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 3)= (unsigned int)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+7;\
} while (0)

#define BYTE_OUT4I(c, insn1, insn2, insn3, insn4, imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(((unsigned char*)c->p->cur_ip) + 3)= (unsigned char)insn4;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 4)= (unsigned int)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+8;\
} while (0)

#define BYTE_OUT5I(c, insn1, insn2, insn3, insn4, insn5, imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(((unsigned char*)c->p->cur_ip) + 3)= (unsigned char)insn4;\
*(((unsigned char*)c->p->cur_ip) + 4)= (unsigned char)insn5;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 5)= (unsigned int)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+9;\
} while (0)

#define BYTE_OUT3S(c, insn1, insn2, insn3, imm16) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(((unsigned char*)c->p->cur_ip) + 2)= (unsigned char)insn3;\
*(unsigned short *)(((unsigned char*)c->p->cur_ip) + 3)= (unsigned short)imm16;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+5;\
} while (0)

#define BYTE_OUT1I(c, insn1, imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 1)= (unsigned int)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+5;\
} while (0)

#define BYTE_OUT1PI(c, insn1, imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(void **)(((unsigned char*)c->p->cur_ip) + 1)= (void *)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+5;\
} while (0)

#define BYTE_OUT2I(c, insn1, insn2,imm32) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 2)= (unsigned int)imm32;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+6;\
} while (0)

#define BYTE_OUT2II(c, insn1, insn2,imm32, imm32_2) \
do { \
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(unsigned char*)c->p->cur_ip = (unsigned char)insn1;\
*(((unsigned char*)c->p->cur_ip) + 1)= (unsigned char)insn2;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 2)= (unsigned int)imm32;\
*(unsigned int *)(((unsigned char*)c->p->cur_ip) + 6)= (unsigned int)imm32_2;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = ((char*)c->p->cur_ip)+10;\
} while (0)

#ifdef LINUX_KERNEL_MODULE
#undef EAX
#undef ECX
#undef EDX
#undef EBX
#undef ESP
#undef EBP
#undef ESI
#undef EDI
#endif

enum {
    AL = 0, CL, DL, BL, AH, CH, DH, BH,   /* r8 */
    AX = 0, CX, DX, BX, SP, BP, SI, DI,   /* r16 */
    EAX = 0, ECX, EDX, EBX, ESP, EBP, ESI, EDI,   /* r32 */
/* used for sse instructions */
    XMM0 = 0, XMM1, XMM2, XMM3, XMM4, XMM5, XMM6, XMM7, XMM8, XMM9, XMM10, XMM11, XMM12, XMM13, XMM14, XMM15
};

typedef struct x86_mach_info {
    int act_rec_size;
    int stack_align;
    int stack_constant_offset;
    int conversion_word;
    int fcu_word;
    int save_base;
    int cur_arg_offset;
    int backpatch_offset;
    int call_stack_adjust_point;
    char pending_prefix;
    int generate_SSE;
} *x86_mach_info;

extern int x86_type_align[];
extern int x86_type_size[];
extern void *gen_x86_mach_info();
extern void x86_arith3(dill_stream s, int op, int commut, int dest, int src1, int src2);
extern void x86_arith2(dill_stream s, int op, int subop, int dest, int src);
extern void x86_mul(dill_stream s, int signed, int imm, int dest, int src1, int src2);

extern void x86_farith(dill_stream s, int b1, int b2, int dest, int src1, int src2);
extern void x86_farith2(dill_stream s, int b1, int b2, int dest, int src);

extern void x86_div_mod(dill_stream s, int sign, int div, int dest, int src1, int src2);
extern void x86_div_modi(dill_stream s, int sign, int div, int dest, int src1, long imm);
extern void x86_arith3i(dill_stream s, int op, int commut, int dest, int src1, long src2);
extern void x86_shift(dill_stream s, int op, int junk, int dest, int src1, int src2);
extern void x86_shifti(dill_stream s, int op, int junk, int dest, int src, long imm);
extern void x86_special(dill_stream s, special_operations type, long param);
extern void x86_set(dill_stream s, int r, long imm);
extern void x86_proc_start(dill_stream s, char *subr_name, int arg_count, 
			   arg_info_list args, dill_reg *arglist);
extern void x86_end(dill_stream s);
extern void x86_package_end(dill_stream s);
extern void *x86_clone_code(dill_stream s, void *base, int size);
extern void x86_ret(dill_stream s, int data1, int data2, int src);
extern void x86_reti(dill_stream s, int data1, int data2, long imm);
extern void x86_retf(dill_stream s, int data1, int data2, double imm);
extern int x86_getreg(dill_stream s, dill_reg *reg_p, int type, int class);
extern int x86_putreg(dill_stream s, dill_reg reg, int type);
extern void
x86_ploadi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
x86_pload(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
x86_pbsloadi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
x86_pbsload(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
x86_pstorei(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
x86_pstore(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
x86_modi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
x86_mod(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
x86_divi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
x86_div(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
x86_converti(dill_stream s, int from_type, int to_type, int dest, long src);
extern void
x86_convert(dill_stream s, int from_type, int to_type, int dest, int src);
extern void
x86_mov(dill_stream s, int type, int junk, int dest, int src);
extern void
x86_pset(dill_stream s, int type, int junk, int dest, long imm);
extern void
x86_setf(dill_stream s, int type, int junk, int dest, double imm);
extern void
x86_setp(dill_stream s, int type, int junk, int dest, void *imm);
extern void
x86_branch(dill_stream s, int op, int type, int src1, int src2, int label);
extern void
x86_branchi(dill_stream s, int op, int type, int src, long imm, int label);
extern void
x86_compare(dill_stream s, int op, int type, int dest, int src1, int src2);
extern void
x86_comparei(dill_stream s, int op, int type, int dest, int src, long imm);
extern void 
x86_lea(dill_stream s, int junk, int junk1, int dest, int src, long imm);
extern void x86_bswap(dill_stream s, int op, int typ, int dest, int src);
extern void x86_jump_to_label(dill_stream s, unsigned long label);
extern void x86_jump_to_reg(dill_stream s, unsigned long reg);
extern void x86_jump_to_imm(dill_stream s, void *imm);
extern void x86_jal(dill_stream s, int return_addr_reg, int target);
extern int x86_calli(dill_stream s, int type, void *xfer_address, const char *name);
extern int x86_callr(dill_stream s, int type, int src);
extern void x86_push(dill_stream s, int type, int reg);
extern void x86_pushi(dill_stream s, int type, long value);
extern void x86_pushpi(dill_stream s, int type, void *value);
extern void x86_pushfi(dill_stream s, int type, double value);
extern int x86_local_op(dill_stream s, int flag, int val);
extern int x86_local(dill_stream s, int type);
extern int x86_localb(dill_stream s, int size);
extern void x86_save_restore_op(dill_stream s, int save_restore, int type,
				 int reg);
extern int x86_init_disassembly_info(dill_stream s, void * ptr);
extern int x86_print_insn(dill_stream s, void *info_ptr, void *insn);
extern int x86_count_insn(dill_stream s, int start, int end);
extern void x86_print_reg(dill_stream s, int typ, int reg);
#endif
