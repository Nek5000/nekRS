#ifndef _PPC64LE_H
#define _PPC64LE_H

extern void ppc64le_XOFORM_arith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void ppc64le_swap_arith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void ppc64le_log_arith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void ppc64le_imm_arith(dill_stream c, int op3, int op, int dest, int src1, long imm);
extern void ppc64le_shift_arith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void ppc64le_shiftimm_arith(dill_stream c, int op3, int op, int dest, int src1, long imm);
extern void ppc64le_logimm_arith(dill_stream c, int op3, int op, int dest, int src1, long imm);
extern void ppc64le_farith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void ppc64le_XFORM2_farith(dill_stream c, int op3, int op, int dest, int src);
extern void ppc64le_FORM2_arith(dill_stream c, int op3, int op, int dest, int src);


#define B_FORM(op, BO, BI, BD, AA, LK) ((op<<26)|(BO<<21)|(BI<<16)|(BD<<2)|(AA<<1)|LK)
#define I_FORM(op, LI, AA, LK) ((op<<26)|(LI<<2)|(AA<<1)|LK)
#define XL_FORM(op, BO, BI, BH, XO, LK) ((op<<26)|(BO<<21)|(BI<<16)|(BH<<11)|(XO<<1)|LK)
#define XO_FORM(op, RT, RA, RB, XO) ((op<<26)|(RT<<21)|(RA<<16)|(RB<<11)|(XO<<1))
#define X_FORM(op, RS, RA, RB, XO) ((op<<26)|(RS<<21)|(RA<<16)|(RB<<11)|(XO<<1))
#define XS_FORM(OP, RS, RA, SHA, XS, SHB) ((OP<<26)|(RS<<21)|(RA<<16)|((SHA)<<11)|((XS)<<2)|((SHB)<<1))
#define XX2_FORM(OP, T, RB, XOP) ((OP<<26)|(T<<21)|(RB<<11)|(XOP<<2))
#define XFX_FORM(OP, RT, SPR, XOP) ((OP<<26)|(RT<<21)|(SPR<<11)|(XOP<<1))
#define XX3_FORM(OP, T, RA, RB, XOP) ((OP<<26)|(T<<21)|(RA<<16)|(RB<<11)|(XOP<<3))
#define D_FORM(OP, RT, RA, SI) ((OP<<26)|(RT<<21)|(RA<<16)|(SI & 0xffff))
#define M_FORM(op, RS, RA, SH, MB, ME) ((op<<26)|(RS<<21)|(RA<<16)|(SH<<11)|(MB<<6)|(ME<<1))
#define MD_FORM(op, RS, RA, SH, MB, XO, SH2) ((op<<26)|(RS<<21)|(RA<<16)|((SH)<<11)|((MB)<<5)|((XO)<<2)|((SH2)<<1))
#define A_FORM(OP, FRT, FRA, FRB, XO) ((OP<<26)|(FRT<<21)|(FRA<<16)|(FRB<<11)|((XO)<<1))
#define INSN_OUT(c, insn) do {\
if (c->p->cur_ip >= c->p->code_limit) {\
   extend_dill_stream(c);\
}\
*(int*)c->p->cur_ip = (unsigned int)insn;\
if (c->dill_debug) dump_cur_dill_insn(c);\
c->p->cur_ip = (void*)(((long)c->p->cur_ip)+4);\
} while (0)\

enum {
    _gpr0 = 0,  _gpr1,  _gpr2,  _gpr3,  _gpr4,  _gpr5,  _gpr6,  _gpr7,
    _gpr8,  _gpr9,  _gpr10,  _gpr11,  _gpr12,  _gpr13,  _gpr14,  _gpr15,
    _gpr16,  _gpr17,  _gpr18,  _gpr19,  _gpr20,  _gpr21,  _gpr22,  _gpr23,
    _gpr24,  _gpr25,  _gpr26,  _gpr27,  _gpr28,  _gpr29,  _gpr30,  _gpr31,


    /* floating point */
    _fpr0 = 0,  _fpr1,  _fpr2,  _fpr3,  _fpr4,  _fpr5,  _fpr6,  _fpr7,
    _fpr8,  _fpr9,  _fpr10,  _fpr11,  _fpr12,  _fpr13,  _fpr14,  _fpr15,
    _fpr16,  _fpr17,  _fpr18,  _fpr19,  _fpr20,  _fpr21,  _fpr22,  _fpr23,
    _fpr24,  _fpr25,  _fpr26,  _fpr27,  _fpr28,  _fpr29,  _fpr30,  _fpr31,


    _sp = _gpr1,/* stack pointer */
    _lr = _gpr2,/* link register */

};

typedef struct ppc64le_mach_info {
    int act_rec_size;
    int stack_align;
    int stack_constant_offset;
    int gp_save_offset;
    int fp_save_offset;
    int fp_save_end;
    int cur_arg_offset;
    int last_int_arg;
    int last_fp_arg;
    int varidiac_call;
    int save_insn_offset;
    int start_label;
    int entry_label;
    int epilogue_label;
} *ppc64le_mach_info;

extern int ppc64le_type_align[];
extern int ppc64le_type_size[];
extern void *gen_ppc64le_mach_info(dill_stream c, int v9);
extern void ppc64le_set(dill_stream c, int r, long imm);
extern void ppc64le_proc_start(dill_stream c, char *subr_name, int arg_count, 
			 arg_info_list args, dill_reg *arglist);
extern void ppc64le_end(dill_stream c);
extern void ppc64le_package_end(dill_stream c);
extern void * ppc64le_clone_code(dill_stream c, void *new_base, int available_size);
extern void ppc64le_ret(dill_stream c, int data1, int data2, int src);
extern void ppc64le_reti(dill_stream c, int data1, int data2, long imm);
extern int ppc64le_getreg(dill_stream c, dill_reg *reg_p, int type, int class);
extern int ppc64le_putreg(dill_stream c, dill_reg reg, int type);
extern void
ppc64le_ploadi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
ppc64le_pload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
ppc64le_pbsloadi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
ppc64le_pbsload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
ppc64le_pstorei(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
ppc64le_pstore(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
ppc64le_modi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
ppc64le_mod(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
ppc64le_divi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
ppc64le_div(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
ppc64le_converti(dill_stream c, int from_type, int to_type, int dest, long src);
extern void
ppc64le_convert(dill_stream c, int from_type, int to_type, int dest, int src);
extern void
ppc64le_mov(dill_stream c, int type, int junk, int dest, int src);
extern void
ppc64le_pset(dill_stream c, int type, int junk, int dest, long imm);
extern void
ppc64le_setf(dill_stream c, int type, int junk, int dest, double imm);
extern void
ppc64le_setp(dill_stream c, int type, int junk, int dest, void *imm);
extern void
ppc64le_branch(dill_stream c, int op, int type, int src1, int src2, int label);
extern void
ppc64le_branchi(dill_stream c, int op, int type, int src, long imm, int label);
extern void
ppc64le_compare(dill_stream c, int op, int type, int dest, int src1, int src2);
extern void
ppc64le_comparei(dill_stream c, int op, int type, int dest, int src, long imm);
extern void 
ppc64le_lea(dill_stream c, int junk, int junk1, int dest, int src, long imm);
extern void ppc64le_bswap(dill_stream c, int junk, int typ, int dest, int src);
extern void ppc64le_jump_to_label(dill_stream c, unsigned long label);
extern void ppc64le_jump_to_reg(dill_stream c, unsigned long reg);
extern void ppc64le_jump_to_imm(dill_stream c, void* imm);
extern void ppc64le_jal(dill_stream c, int return_addr_reg, int target);
extern int ppc64le_calli(dill_stream c, int type, void *xfer_address, const char*name);
extern int ppc64le_callr(dill_stream c, int type, int src);
extern void ppc64le_push(dill_stream c, int type, int reg);
extern void ppc64le_pushi(dill_stream c, int type, long value);
extern void ppc64le_pushfi(dill_stream c, int type, double value);
extern void ppc64le_pushpi(dill_stream c, int type, void *value);
extern int ppc64le_local_op(dill_stream c, int flag, int val);
extern int ppc64le_local(dill_stream c, int type);
extern int ppc64le_localb(dill_stream c, int size);
extern void ppc64le_save_restore_op(dill_stream c, int save_restore, int type,
				 int reg);
extern int ppc64le_init_disassembly_info(dill_stream c, void * ptr);
extern int ppc64le_print_insn(dill_stream c, void *info_ptr, void *insn);
extern int ppc64le_count_insn(dill_stream c, int start, int end);
extern void ppc64le_print_reg(dill_stream c, int typ, int reg);
#endif
