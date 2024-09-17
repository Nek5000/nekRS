#ifndef _ARM5_H
#define _ARM5_H
enum arm5_cond {EQ, NE, CS, CC, MI, PL, VS, VC, HI, LS, GE, LT, GT, LE, AL, NV};
enum arm5_opcode {AND, EOR, SUB, RSB, ADD, ADC, SBC, RSC, TST, TEQ, CMP, CMN, ORR, MOV, BIC, MVN};

extern void arm5_dproc(dill_stream s, int op, int shift_code, int dest, int src1, int src2);
extern void arm5_dproc2(dill_stream s, int op3, int op, int dest, int src);
extern void arm5_dproci(dill_stream s, int op, int shift_code, int dest, int src1, long imm);
extern void arm5_fproc(dill_stream s, int op3, int op, int dest, int src1, int src2);
extern void arm5_fproc2(dill_stream s, int op3, int op, int dest, int src);
extern void arm5_fproci(dill_stream s, int op3, int op, int dest, int src1, long imm);

enum {
    _r0,  _r1,  _r2,  _r3,  _r4,  _r5,  _r6,  _r7, /* globals */
    _r8,  _r9,  _r10,  _r11,  _r12,  _r13,  _r14,  _r15, /* globals */

    _sp = _r13,	/* stack pointer */
    _link = _r14,	/* link address */
    _pc = _r15,	/* program counter */
    _fp = _r11,
    _a1 = _r0, _a2 = _r1, _a3 = _r2, _a4 = _r3,
    _v1 = _r4, _v2 = _r5, _v3 = _r6, _v4 = _r7, _v5 = _r8, _v6 = _r9, 
    _v7 = _r10,

    /* floating point */
    _f0=0,  _f1,  _f2,  _f3,  _f4,  _f5,  _f6,  _f7
};

#define LLshift 0x10
#define LRshift 0x11
#define ARshift 0x12

typedef struct arm5_mach_info {
    int act_rec_size;
    int stack_align;
    int stack_constant_offset;
    int gp_save_offset;
    int fp_save_offset;
    int fp_save_end;
    int conversion_word;
    int cur_arg_offset;
    int save_insn_offset;
    int max_arg_size;
} *arm5_mach_info;

extern int arm5_type_align[];
extern int arm5_type_size[];
extern void *gen_arm5_mach_info(dill_stream s, int v9);
extern void arm5_set(dill_stream s, int r, long imm);
extern void arm5_proc_start(dill_stream s, char *subr_name, int arg_count, 
			 arg_info_list args, dill_reg *arglist);
extern void arm5_end(dill_stream s);
extern void *arm5_clone_code(dill_stream s, void *base, int size);
extern void arm5_ret(dill_stream s, int data1, int data2, int src);
extern void arm5_reti(dill_stream s, int data1, int data2, long imm);
extern int arm5_getreg(dill_stream s, dill_reg *reg_p, int type, int class);
extern int arm5_putreg(dill_stream s, dill_reg reg, int type);
extern void
arm5_ploadi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
arm5_pload(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
arm5_pbsloadi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
arm5_pbsload(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
arm5_pstorei(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
arm5_pstore(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
arm5_modi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
arm5_mod(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
arm5_divi(dill_stream s, int type, int junk, int dest, int src, long offset);
extern void
arm5_mul(dill_stream s, int unsign, int junk, int dest, int src1, int src2);
extern void
arm5_muli(dill_stream s, int unsign, int junk, int dest, int src, long imm);
extern void
arm5_div(dill_stream s, int type, int junk, int dest, int src1, int src2);
extern void
arm5_convert(dill_stream s, int from_type, int to_type, int dest, int src);
extern void
arm5_mov(dill_stream s, int type, int junk, int dest, int src);
extern void
arm5_pset(dill_stream s, int type, int junk, int dest, long imm);
extern void
arm5_setf(dill_stream s, int type, int junk, int dest, double imm);
extern void
arm5_branch(dill_stream s, int op, int type, int src1, int src2, int label);
extern void
arm5_branchi(dill_stream s, int op, int type, int src, long imm, int label);
extern void 
arm5_lea(dill_stream s, int junk, int junk1, int dest, int src, long imm);
extern void arm5_jump_to_label(dill_stream s, unsigned long label);
extern void arm5_jump_to_reg(dill_stream s, unsigned long reg);
extern void arm5_jump_to_imm(dill_stream s, void *imm);
extern void arm5_jal(dill_stream s, int return_addr_reg, int target);
extern int arm5_calli(dill_stream s, int type, void *xfer_address, const char *name);
extern int arm5_callr(dill_stream s, int type, int src);
extern void arm5_push(dill_stream s, int type, int reg);
extern void arm5_pushi(dill_stream s, int type, long value);
extern void arm5_pushfi(dill_stream s, int type, double value);
extern void arm5_pushpi(dill_stream s, int type, void *value);
extern int arm5_local_op(dill_stream s, int flag, int val);
extern int arm5_local(dill_stream s, int type);
extern int arm5_localb(dill_stream s, int size);
extern void arm5_save_restore_op(dill_stream s, int save_restore, int type,
				 int reg);
extern int arm5_init_disassembly_info(dill_stream s, void * ptr);
extern int arm5_print_insn(dill_stream s, void *info_ptr, void *insn);
extern int arm5_count_insn(dill_stream s, int start, int end);
extern void arm5_print_reg(dill_stream s, int typ, int reg);
#endif
