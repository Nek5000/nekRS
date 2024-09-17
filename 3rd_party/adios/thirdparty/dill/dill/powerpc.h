#ifndef _POWERPC_H
#define _POWERPC_H
 #define OP(x)		(((x)&0x6) << 26) /* opcd field of insns */

 #define HDR(x)		((unsigned)((x)&0x3) << 30) /* hdr field of all insns */
 #define ASI(asi) 	((asi) << 5)
 #define OPF(x)		(((x)&0x1ff) << 5) /* opf field of float insns */
 #define RD(x)		(((x)&0x1f) << 25) /* destination register field */
 #define RS1(x)		(((x)&0x1f) << 14) /* rs1 field */
 #define RS2(x)		((x)&0x1f) /* rs2 field */
 #define SIMM13(x)	((x)&0x1fff) /* simm13 field */
 #define IM		0x2000
 #define COND(x) (((x)&0xf)<<25)
 #define CC(x)  (((x)&0x3)<<20)
 #define P(x)  (((x)&0x1)<<19)
/* upper 22 of 64 bits */
#define hh(x) (((unsigned long) x) >> 42)
/* bits 11-32 */
#define lm(x) hi((x) & (((unsigned long)1<<32) -1))
/* lower 10 of upper 32 */
#define hm(x) lo( ((unsigned long) x)>>32)

/* upper 22 bits of 32 */
#define hi(x) (((unsigned long)x) >> 10)
/* lower 10 bits */
#define lo(x) ((x) & ((1 << 11) - 1))

extern void powerpc_FORM3_arith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void powerpc_FORM3imm_arith(dill_stream c, int op3, int op, int dest, int src1, long imm);
extern void powerpc_FORM3_farith(dill_stream c, int op3, int op, int dest, int src1, int src2);
extern void powerpc_FORM2_farith(dill_stream c, int op3, int op, int dest, int src);


#define D_FORM(op, r1, r2, i) ((op<<26)|(r1<<21)|(r2<<16)|i)
#define XL_FORM(op, BO, BI, BH, XO, LK) ((op<<26)|(BO<<21)|(BI<<16)|(BH<<11)|(XO<<1)|LK)
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


    _sp = _gpr3,/* stack pointer */

};

typedef struct powerpc_mach_info {
    int act_rec_size;
    int stack_align;
    int stack_constant_offset;
    int gp_save_offset;
    int fp_save_offset;
    int fp_save_end;
    int conversion_word;
    int cur_arg_offset;
    int save_insn_offset;
} *powerpc_mach_info;

extern int powerpc_type_align[];
extern int powerpc_type_size[];
extern void *gen_powerpc_mach_info(dill_stream c, int v9);
extern void powerpc_set(dill_stream c, int r, long imm);
extern void powerpc_proc_start(dill_stream c, char *subr_name, int arg_count, 
			 arg_info_list args, dill_reg *arglist);
extern void powerpc_end(dill_stream c);
extern void powerpc_package_end(dill_stream c);
extern void * powerpc_clone_code(dill_stream c, void *new_base, int available_size);
extern void powerpc_ret(dill_stream c, int data1, int data2, int src);
extern void powerpc_reti(dill_stream c, int data1, int data2, long imm);
extern int powerpc_getreg(dill_stream c, dill_reg *reg_p, int type, int class);
extern int powerpc_putreg(dill_stream c, dill_reg reg, int type);
extern void
powerpc_ploadi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
powerpc_pload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
powerpc_pbsloadi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
powerpc_pbsload(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
powerpc_pstorei(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
powerpc_pstore(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
powerpc_modi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
powerpc_mod(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
powerpc_divi(dill_stream c, int type, int junk, int dest, int src, long offset);
extern void
powerpc_div(dill_stream c, int type, int junk, int dest, int src1, int src2);
extern void
powerpc_converti(dill_stream c, int from_type, int to_type, int dest, long src);
extern void
powerpc_convert(dill_stream c, int from_type, int to_type, int dest, int src);
extern void
powerpc_mov(dill_stream c, int type, int junk, int dest, int src);
extern void
powerpc_pset(dill_stream c, int type, int junk, int dest, long imm);
extern void
powerpc_setf(dill_stream c, int type, int junk, int dest, double imm);
extern void
powerpc_setp(dill_stream c, int type, int junk, int dest, void *imm);
extern void
powerpc_branch(dill_stream c, int op, int type, int src1, int src2, int label);
extern void
powerpc_branchi(dill_stream c, int op, int type, int src, long imm, int label);
extern void
powerpc_compare(dill_stream c, int op, int type, int dest, int src1, int src2);
extern void
powerpc_comparei(dill_stream c, int op, int type, int dest, int src, long imm);
extern void 
powerpc_lea(dill_stream c, int junk, int junk1, int dest, int src, long imm);
extern void powerpc_bswap(dill_stream c, int junk, int typ, int dest, int src);
extern void powerpc_jump_to_label(dill_stream c, unsigned long label);
extern void powerpc_jump_to_reg(dill_stream c, unsigned long reg);
extern void powerpc_jump_to_imm(dill_stream c, unsigned long imm);
extern void powerpc_jal(dill_stream c, int return_addr_reg, int target);
extern int powerpc_calli(dill_stream c, int type, void *xfer_address, char*name);
extern int powerpc_callr(dill_stream c, int type, int src);
extern void powerpc_push(dill_stream c, int type, int reg);
extern void powerpc_pushi(dill_stream c, int type, long value);
extern void powerpc_pushfi(dill_stream c, int type, double value);
extern void powerpc_pushpi(dill_stream c, int type, void *value);
extern int powerpc_local_op(dill_stream c, int flag, int val);
extern int powerpc_local(dill_stream c, int type);
extern int powerpc_localb(dill_stream c, int size);
extern void powerpc_save_restore_op(dill_stream c, int save_restore, int type,
				 int reg);
extern int powerpc_init_disassembly_info(dill_stream c, void * ptr);
extern int powerpc_print_insn(dill_stream c, void *info_ptr, void *insn);
extern int powerpc_count_insn(dill_stream c, int start, int end);
extern void powerpc_print_reg(dill_stream c, int typ, int reg);
#endif
