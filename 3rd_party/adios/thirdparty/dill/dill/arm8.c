#include "dill.h"
#include "dill_internal.h"
#include "arm8.h"
#include "config.h"
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#undef NDEBUG
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#define INSN_OUT(s, insn) do {\
if (s->p->cur_ip >= s->p->code_limit) {\
   extend_dill_stream(s);\
}\
*(int*)s->p->cur_ip = (unsigned int)insn;\
if (s->dill_debug) dump_cur_dill_insn(s);\
s->p->cur_ip = (void*)(((long)s->p->cur_ip)+4);\
} while (0)\

 #define COND(x)	((unsigned)((x)&0xf) << 28)
 #define CLASS(x)	(((x)&0x7) << 25)
 #define OPCODE(x)	(((x)&0xf) << 21) /* opcode field */
 #define p(x)		(((x)&0x1) << 23)
 #define D(x)		(((x)&0x1) << 22)
 #define q(x)		(((x)&0x1) << 21)
 #define r(x)		(((x)&0x1) << 20)
 #define S(x)		(((x)&0x1) << 20) /* update cond codes? */
 #define FN(x)		(((x)&0xf) << 16) /* Fn field */
 #define FD(x)		(((x)&0xf) << 12) /* Fn field */
 #define cp_num(x)	(((x)&0xf) << 8) /* cpu */
 #define N(x)		(((x)&0x1) << 7)
 #define s(x)		(((x)&0x1) << 6)
 #define M(x)		(((x)&0x1) << 5)
 #define RN(x)		(((x)&0xf) << 16) /* Rn field */
 #define RD(x)		(((x)&0xf) << 12) /* Rd field */
 #define RM(x)		(((x)&0xf) << 0) /* Rm field */
 #define SHIFTI(x,t)	((((x)&0x1f) << 7) | ((t)&0x3)<<5)
 #define SHIFTR(r,t)	((((r)&0xf) << 8) | ((t)&0x3)<<5| 1<<4)

 #define IMM(x,r)	(((x)&0xff) | ((((32-(r))>>1)&0xf)<< 8) | (1<<25)) /* simm8 field */
 #define IM		0x2000
 #define P(x)  (((x)&0x1)<<19)

#define arm8_savei(s, imm) arm8_dproci(s, SUB, 0, _sp, _sp, ar_size);
#define arm8_andi(s, dest, src, imm) arm8_dproci(s, AND, 0, dest, src, imm)
#define arm8_movi(s, dest, src) arm8_dproc(s, MOV, 0, dest, 0, src) 
#define arm8_movf(s, dest, src) arm8_fproc2(s, 0, 0, 0, dest, src)
#define arm8_movd(s, dest, src) arm8_fproc2(s, 0, 1, 0, dest, src)
#define arm8_lshi(s, dest, src,imm) arm8_dproci(s, MOV, LLshift, dest, src, imm)
#define arm8_rshi(s,dest,src,imm) arm8_dproci(s, MOV, LRshift, dest, src, imm)
#define arm8_rshai(s,dest,src,imm) arm8_dproci(s, MOV, ARshift, dest, src, imm)

#define arm8_nop(s) arm8_movi(s, _r0, _r0)
#define arm8_raw_push(s, reg) INSN_OUT(s, COND(AL)|CLASS(2)|1<<24|1<<21|0xd<<16|reg<<12|4)
#define arm8_raw_pop(s, reg) INSN_OUT(s, COND(AL)|CLASS(2)|1<<23|1<<20|0xd<<16|reg<<12|4)
#define IREG 0
#define FREG 1

#define roundup(a,b) ((a + (b-1)) & (-b))

static void
arm8_pldsti(dill_stream s, int type, int ls, int dest, int src, long offset);
static void
arm8_pldst(dill_stream s, int type, int ls, int dest, int src1, int src2);
static void int_arm8_bswap(dill_stream s, int type, int reg);

extern void 
arm8_bswap(dill_stream s, int type, int data2, int dest, int src)
{
    switch(type) {
    case DILL_L: case DILL_UL: case DILL_I: case DILL_U:
	INSN_OUT(s, COND(AL)| 0b011010111111<<16|RD(dest)|RM(src)|0xf<<8|3<<4);
	break;
    case DILL_US: case DILL_S:
	INSN_OUT(s, COND(AL)| 0b011010111111<<16|RD(dest)|RM(src)|0xf<<8|0xb<<4);
	break;
    case DILL_C: case DILL_UC:
	/* nothing to do */
	break;
    case DILL_F:
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|1<<20|FN(src>>1)|N(src)|RD(_r0)|cp_num(0xa)|1<<4);/*fmrs*/
	int_arm8_bswap(s, DILL_L, _r0);
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(dest>>1)|N(dest)|RD(_r0)|cp_num(0xa)|1<<4);/*fmsr*/
	break;
    case DILL_D: 
	INSN_OUT(s, COND(AL)| CLASS(6)|OPCODE(2)|r(1)|FN(_r1)|RD(_r0)|cp_num(0xb)|1<<4|(src>>1)|(src&1)<<5);/*fmrrd*/
	int_arm8_bswap(s, DILL_L, _r0);
	int_arm8_bswap(s, DILL_L, _r1);
	INSN_OUT(s, COND(AL)| CLASS(6)|OPCODE(2)|FN(_r0)|RD(_r1)|cp_num(0xb)|1<<4|(dest>>1)|(dest&1)<<5);/*fmdrr*/
	break;
    }
}

static void
int_arm8_bswap(dill_stream s, int type, int reg)
{
  arm8_bswap(s, type, 0, reg, reg);
}

extern void
arm8_pbsloadi(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    arm8_pldsti(s, type, 1, dest, src, offset);
    int_arm8_bswap(s, type, dest);
}


extern void
arm8_pbsload(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    arm8_pldst(s, type, 1, dest, src1, src2);
    int_arm8_bswap(s, type, dest);
}

static 
struct basic_type_info 
{   char size;
    char align;
    char reg_type;
} type_info[] = {
    { 1, 1, IREG},  /* C */
    { 1, 1, IREG},  /* UC */
    { 2, 2, IREG},  /* S */
    { 2, 2, IREG},  /* US */
    { 4, 4, IREG},  /* I */
    { 4, 4, IREG},  /* U */
    { sizeof(long), sizeof(long), IREG},  /* UL */
    { sizeof(long), sizeof(long), IREG},  /* L */
    { sizeof(char*), sizeof(char*), IREG},  /* P */
    { sizeof(float), sizeof(float), FREG},  /* F */
    { sizeof(double), 4, FREG},  /* D */
    { 0, 8, IREG}, /* V */
    { -1, 8, IREG}, /* B */
    { 4, 8, IREG}, /* EC */
};

int arm8_type_align[] = {
        1, /* C */
        1, /* UC */
        2, /* S */
        2, /* US */
        4, /* I */
        4, /* U */
        sizeof(unsigned long), /* UL */
        sizeof(long), /* L */
        sizeof(char*), /* P */
        4, /* F */
        4, /* D */
	1, /* V */
        4, /* B */
	sizeof(long), /* EC */
};

int arm8_type_size[] = {
        1, /* C */
        1, /* UC */
        2, /* S */
        2, /* US */
        4, /* I */
        4, /* U */
        sizeof(unsigned long), /* UL */
        sizeof(long), /* L */
        sizeof(char*), /* P */
        4, /* F */
        8, /* D */
	1, /* V */
        8, /* B */
        sizeof(char*), /* EC */
};

extern void arm8_dproc(s, op, shift_code, dest, src1, src2)
dill_stream s;
int op;
int shift_code;
int dest;
int src1;
int src2;
{
    int shift = 0;
    if (shift_code != 0) {
	shift_code &= 0x3;
	shift = SHIFTR(src2, shift_code);
	src2 = src1;
	src1 = 0;
	
    }
    INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(op)|S(0)|RN(src1)|RD(dest)|RM(src2)|shift);
}

extern void arm8_dproc2(s, op, fop, dest, src)
dill_stream s;
int op;
int fop;
int dest;
int src;
{
    if (op == RSB) {
	arm8_dproci(s, RSB, 0, dest, src, 0);
    } else if (op == CMN) {  /* !a */
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(CMP)|S(1)|RN(src)|RD(src)|IMM(0,0));
	INSN_OUT(s, COND(NE)|CLASS(0x0)|OPCODE(MOV)|S(0)|RN(0)|RD(dest)|IMM(0, 0));
	INSN_OUT(s, COND(EQ)|CLASS(0x0)|OPCODE(MOV)|S(0)|RN(0)|RD(dest)|IMM(1, 0));
    } else {
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(op)|S(0)|RN(src)|RD(dest)|RM(src));
    }
}

extern void arm8_negf(s, op,fd, dest, src)
dill_stream s;
int op;
int fd;
int dest;
int src;
{
    arm8_fproc2(s, op, fd, 0, dest, src);
}

extern void arm8_fproc2(s, op,fd, n, dest, src)
dill_stream s;
int op;
int fd;
int n;
int dest;
int src;
{
  INSN_OUT(s, COND(AL)|CLASS(0x7)|p(1)|D(dest&1)|q(1)|r(1)|FN(op)|N(n)|FD(dest>>1)|(0xa+fd)<<8|s(1)|M(src&1)|((src>>1)&0xf));
}

extern int
arm8_local(dill_stream s, int type)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;

    ami->act_rec_size += roundup(type_info[type].size, ami->stack_align);
    return (-ami->act_rec_size)  - 14 * 4 /* int regs to save */ 
	- 8 * 3 * 4 /* float regs to save */;
}

extern int
arm8_localb(dill_stream s, int size)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    if (size < 0) size = 0;
    ami->act_rec_size = roundup(ami->act_rec_size, size);

    ami->act_rec_size += roundup(size, ami->stack_align);
    return (-ami->act_rec_size) - 14 * 4 /* int regs to save */ 
	- 8 * 3 * 4 /* float regs to save */;
}

extern int arm8_local_op(dill_stream s, int flag, int val)
{
    int size = val;
    if (flag == 0) {
	size = type_info[val].size;
    }
    if (size < 0) size = 0;
    return arm8_localb(s, size);
}	

static int 
is_temp(int ireg)
{
    return (ireg <= _r4);  /* higher regs are saved by the callee */
}

static int 
is_ftemp(int freg)
{
    return (freg <= _f4);  /* higher regs are saved by the callee */
}

extern void
arm8_save_restore_op(dill_stream s, int save_restore, int type, int reg)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    if (save_restore == 0) { /* save */
	switch (type) {
	case DILL_D: case DILL_F:
	    if (is_ftemp(reg)) {
		arm8_pstorei(s, type, 0, reg, _fp, - 13*4 - reg * 12);
	    }
	    break;
	default:
	    if (is_temp(reg)) {
		arm8_pstorei(s, type, 0, reg, _sp, ami->gp_save_offset + (reg - _r0) * ami->stack_align);
	    }
	    break;
	}
    } else {  /* restore */
	switch (type) {
	case DILL_D: case DILL_F:
	    if (is_ftemp(reg)) {
		arm8_ploadi(s, type, 0, reg, _fp, -13*4 - reg * 12);
	    }
	    break;
	default:
	    if (is_temp(reg)) {
		arm8_ploadi(s, type, 0, reg, _sp, ami->gp_save_offset + (reg - _r0) * ami->stack_align);
	    }
	    break;
	}
    }
}	

static void
arm8_movi2f(dill_stream s, int dest, int src)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    arm8_pstorei(s, DILL_I, 0, src, _fp, ami->conversion_word);
    arm8_ploadi(s, DILL_F, 0, dest, _fp, ami->conversion_word);
}
    
static void
arm8_movf2i(dill_stream s, int dest, int src)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    arm8_pstorei(s, DILL_F, 0, src, _fp, ami->conversion_word);
    arm8_ploadi(s, DILL_I, 0, dest, _fp, ami->conversion_word);
}
    
static void
arm8_movd2i(dill_stream s, int dest, int src)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    arm8_pstorei(s, DILL_D, 0, src, _fp, ami->conversion_word);
    if (ami->stack_align == 8) {
	arm8_ploadi(s, DILL_L, 0, dest, _fp, ami->conversion_word);
    } else {
	arm8_ploadi(s, DILL_I, 0, dest, _fp, ami->conversion_word);
	arm8_ploadi(s, DILL_I, 0, dest+1, _fp, ami->conversion_word+4);
    }
}
    
static void
arm8_movi2d(dill_stream s, int dest, int src)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    if (ami->stack_align == 8) {
	arm8_pstorei(s, DILL_L, 0, src, _fp, ami->conversion_word);
    } else {
	arm8_pstorei(s, DILL_I, 0, src, _fp, ami->conversion_word);
	arm8_pstorei(s, DILL_I, 0, src+1, _fp, ami->conversion_word+4);
    }
    arm8_ploadi(s, DILL_D, 0, dest, _fp, ami->conversion_word);
}
    
extern void arm8_fproc(s, arm8_op, fd, dest, src1, src2)
dill_stream s;
int arm8_op;
int fd;
int dest;
int src1;
int src2;
{
  INSN_OUT(s, COND(AL)|CLASS(0x7)|D(dest&0x1)|p(arm8_op>>3)|q(arm8_op>>2)|r(arm8_op>>1)|s(arm8_op)|(arm8_op&0x1)<<15|((src1>>1)&0xf)<<16|((dest>>1)&0xf)<<12|(0xa+(fd))<<8|(src1&1)<<7|(src2>>1)&0xf|(src2&0x1)<<5);
}

extern void arm8_dproci(s, op, shift_code, dest, src1, imm)
dill_stream s;
int op;
int shift_code;
int dest;
int src1;
long imm;
{
    int shift = 0;
    int setcc = 0;
    if (op == CMP) setcc = 1;
    if (shift_code != 0) {
	/* must already be a mov op */
	shift_code &= 0x3;
	shift = SHIFTI(imm, shift_code);
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(op)|S(0)|RN(src1)|RD(dest)|shift|RM(src1));
	return;
    }
    if ((imm >= 0) && (imm < 256)) {
	/* arith format */
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(op)|S(setcc)|RN(src1)|RD(dest)|IMM(imm, 0));
    } else {
	arm8_set(s, _v1, imm);
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(op)|S(setcc)|RN(src1)|RD(dest)|RM(_v1));
    }
}

/*
 *	ARM stack frame organization
 *	
 *		pushed args
 *			------------  SP value at entry    FP value
 *		callee-saved int regs 
 *		pushed with stmdb	space for 14
 *		callee-saved float regs space for 8
 *
 *	 		------------  SP value after STMDB
 *		local variables
 *			------------  final SP value
 */

extern void
arm8_proc_start(dill_stream s, char *subr_name, int arg_count, arg_info_list args,
	     dill_reg *arglist)
{
    int i;
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    int cur_arg_offset = 0;
    int next_core_register = _r0;
    int next_float_register = _f0;

    /* emit start insns */
    INSN_OUT(s, 0xFF000000);
    INSN_OUT(s, 0xFF000000);
    INSN_OUT(s, 0xFF000000);
    arm8_movi(s, _r12, _sp);
    /* stmdb sp!, {r11, r12, lr, pc} */
    INSN_OUT(s, COND(AL)|CLASS(4)|1<<24/*p*/|RN(_sp)|1<<_r11|1<<_r12|1<<_link|1<<_pc);
    arm8_dproci(s, SUB, 0, _sp, _sp, 14*4 ); /* instead of write back */
    arm8_nop(s);  /* placeholder for float save */
    arm8_dproci(s, SUB, 0, _r11, _r12, 4);
    ami->save_insn_offset = (long)s->p->cur_ip - (long)s->p->code_base;
    arm8_nop(s);	/* room for largest stack adjust insn, 5 nops */
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
    ami->conversion_word = arm8_local(s, DILL_D);
    ami->conversion_word = arm8_local(s, DILL_D);
    ami->conversion_word = arm8_local(s, DILL_D);

    /* load params from regs */
    for (i = 0; i < arg_count; i++) {
        int item_size = 1;
	switch (args[i].type) {
	case DILL_D:
	    if (ami->hard_float) {
	        if (next_float_register % 2) {
		    /* double is only even regs, skip one */
		    next_float_register++;
		}
	    } else {
	        if (next_core_register % 2) {
		    /* double is only even regs, skip one */
		    next_core_register++;
		    cur_arg_offset += 4;
		}
	    }
	    item_size = 2;
	case DILL_F:
	    /* falling through */
	    if (ami->hard_float) {
	        if (next_float_register <= _f31) {
	            args[i].is_register = 1;
		    args[i].in_reg = next_float_register;
		    args[i].out_reg = next_float_register;
		} else {
		    args[i].is_register = 0;
		}
		next_float_register += ((args[i].type == DILL_D) ? 2 : 1);
		break;
	    }
	    /* if soft float, fall through, with item size at 2 for D */
	default:
	    if (next_core_register < _r4) {
		args[i].is_register = 1;
		args[i].in_reg = next_core_register;
		args[i].out_reg = next_core_register;
	    } else {
		args[i].is_register = 0;
	    }
	    next_core_register+=item_size;
	    break;
	}
	args[i].offset = cur_arg_offset;
	cur_arg_offset += roundup(type_info[(int)args[i].type].size, ami->stack_align);
    }
    
    for (i = 0; i < arg_count; i++) {
	int tmp_reg;
	//	printf("Handling arg %d, is_reg %d\n", i, args[i].is_register);
	if (args[i].is_register) {
	    /* only some moved into registers */
	    if (!dill_raw_getreg(s, &tmp_reg, args[i].type, DILL_VAR)) {
		/* not enough regs for this, store it to the stack */
		int real_offset = - args[i].offset - 4*4; 
		if (arglist != NULL) arglist[i] = -1;
		arm8_pstorei(s, DILL_I, 0, args[i].in_reg, _fp, 
				    real_offset);
		args[i].in_reg = -1;
		args[i].out_reg = -1;
		args[i].offset = real_offset;
		args[i].is_register = 0;
		continue;
	    }
	    if (args[i].is_register) {
		if ((args[i].type != DILL_F) && (args[i].type != DILL_D)) {
		    arm8_movi(s, tmp_reg, args[i].in_reg);
		} else if (args[i].type == DILL_F) {	    /* must be float */
		    if (ami->hard_float) {
		        arm8_movf(s, tmp_reg, args[i].in_reg);
		    } else {
		        int src = args[i].in_reg;
			int dest = tmp_reg;
		        INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(dest>>1)|N(dest)|RD(src)|cp_num(0xa)|1<<4);/*fmsr*/
		    }
		} else {
		    if (ami->hard_float) {
		        arm8_movd(s, tmp_reg, args[i].in_reg);
		    } else {
		        INSN_OUT(s, COND(AL)| CLASS(6)|OPCODE(2)|FN(args[i].in_reg+1)|RD(args[i].in_reg)|cp_num(0xb)|1<<4|(tmp_reg>>1)|(tmp_reg&1)<<5);/*fmdrr*/
		    }			  
		}
	    } else {
		/* general offset from fp*/
		int real_offset = args[i].offset - 3*4; 
		arm8_ploadi(s, args[i].type, 0, tmp_reg, _fp, real_offset);
	    }
	    if (arglist != NULL) arglist[i] = tmp_reg;
	    args[i].in_reg = tmp_reg;
	    args[i].is_register = 1;
	} else {
	    if (!ami->hard_float && ((args[i].type == DILL_F) || (args[i].type == DILL_D))) {
	        int tmp_reg;
		dill_raw_getreg(s, &tmp_reg, args[i].type, DILL_VAR);
		/* general offset from fp*/
		int real_offset = args[i].offset - 3*4; 
		arm8_ploadi(s, args[i].type, 0, tmp_reg, _fp, real_offset);
		if (arglist != NULL) arglist[i] = tmp_reg;
		args[i].in_reg = tmp_reg;
		args[i].is_register = 1;
	    } else {
	        /* leave it on the stack */
	        int real_offset = args[i].offset - 3*4; 
		if (arglist != NULL) arglist[i] = -1;
		args[i].in_reg = -1;
		args[i].out_reg = -1;
		args[i].offset = real_offset;
	    }
	}
    }
}


extern void
arm8_ploadi(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    arm8_pldsti(s, type, 1, dest, src, offset);
}

extern void
arm8_pload(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    arm8_pldst(s, type, 1, dest, src1, src2);
}

/* byte and whole word version */
#define ARM8_LDSTI(s,u,b,ls,rn,rd,offset) INSN_OUT(s, COND(AL)|CLASS(2)|(1<<24)|((u&1)<<23)|((b&1)<<22)|(ls&1)<<20|RN(rn)|RD(rd)|(0x7ff&offset))

/* halfword version */
#define ARM8_LDSTHI(s,u,ls,rn,rd,sh,offset) INSN_OUT(s, COND(AL)|CLASS(0)|(1<<24)|((u&1)<<23)|(1<<22)|(ls&1)<<20|RN(rn)|RD(rd)|(1<<7)|((sh&0x3)<<5)|(1<<4)|(0xf&offset)|((offset&0xf0)<<4))

/* float version */
#define ARM8_LDSTFI(s,u,fd,ls,rn,rd,offset) INSN_OUT(s, COND(AL)|CLASS(6)|(1<<24)|((u&1)<<23)|(ls&1)<<20|RN(rn)|D(rd&1)|RD(rd>>1)|((0xa+(fd))<<8)|(0xff&(offset>>2)))

extern void
arm8_pstorei(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    arm8_pldsti(s, type, 0, dest, src, offset);
}

static void
arm8_pldsti(dill_stream s, int type, int ls, int dest, int src, long offset)
{
    int u = 1;
    int max_offset;

    switch (type) {
    case DILL_S: case DILL_US: case DILL_D: case DILL_F:
	max_offset = 256;
	break;
    default:
	max_offset = 2048;
	break;
    }
    if  (((long)offset) >= max_offset || ((long)offset) < -max_offset) {
	arm8_set(s, _v1, offset);
	arm8_pldst(s, type, ls, dest, src, _v1);
	return;
    }
    if (offset < 0) {
	u = 0;
	offset = -offset;
    }

    switch (type) {
    case DILL_F:
	ARM8_LDSTFI(s, u, 0, ls, src, dest, offset);
	break;
    case DILL_D:
	ARM8_LDSTFI(s, u, 1, ls, src, dest, offset);
	break;
    case DILL_C:
    case DILL_UC:
	ARM8_LDSTI(s, u, 1, ls, src, dest, offset);
	break;
    case DILL_I:
    case DILL_U:
    case DILL_L:
    case DILL_UL:
    case DILL_P:
	ARM8_LDSTI(s, u, 0, ls, src, dest, offset);
	break;
    case DILL_S:
	if (ls == 1) { /* this is a load */
	    ARM8_LDSTHI(s,u,ls,src,dest,0x3,offset);
	    break;
	}
	/* fall through */
    case DILL_US:
	ARM8_LDSTHI(s,u,ls,src,dest,0x1,offset);
	break;
    default:
	break;
    }
}
#define ARM8_LDST(s,u,b,ls,rn,rd,rm) INSN_OUT(s, COND(AL)|CLASS(3)|(1<<24)|((u&1)<<23)|((b&1)<<22)|(ls&1)<<20|RN(rn)|RD(rd)|(0xf&rm))
#define ARM8_LDSTH(s,u,ls,rn,rd,sh,rm) INSN_OUT(s, COND(AL)|CLASS(0)|(1<<24)|((u&1)<<23)|(ls&1)<<20|RN(rn)|RD(rd)|(1<<7)|((sh&0x3)<<5)|(1<<4)|(0xf&rm))

extern void
arm8_pstore(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    arm8_pldst(s, type, 0, dest, src1, src2);
}

static void
arm8_pldst(dill_stream s, int type, int ls, int dest, int src1, int src2)
{
    switch (type) {
    case DILL_F:
	arm8_dproc(s, ADD, 0, _v1, src1, src2);
	ARM8_LDSTFI(s, 0, 0, ls, _v1, dest, 0);
	break;
    case DILL_D:
	arm8_dproc(s, ADD, 0, _v1, src1, src2);
	ARM8_LDSTFI(s, 0, 1, ls, _v1, dest, 0);
	break;
    case DILL_L: case DILL_UL: case DILL_P: case DILL_I: case DILL_U: case DILL_EC:
	ARM8_LDST(s,1,0,ls,src1,dest,src2);
	break;
    case DILL_S:
	if (ls == 1) { /* this is a load */
	    ARM8_LDSTH(s,1,ls,src1,dest,0x3,src2);
	    break;
	}
	/* fall through */
    case DILL_US:
	ARM8_LDSTH(s,1,ls,src1,dest,0x1,src2);
	break;
    case DILL_C: case DILL_UC:
	ARM8_LDST(s,1,1,ls,src1,dest,src2);
	break;
    default:
	break;
    }
}

extern int arm8_hidden_modi(int a, int b);
extern long arm8_hidden_mod(long a, long b);
extern unsigned long arm8_hidden_umod(unsigned long a, unsigned long b);
extern unsigned int arm8_hidden_umodi(unsigned int a, unsigned int b);
extern double arm8_hidden_ultod(unsigned long a);
extern float arm8_hidden_ultof(unsigned long a);
extern unsigned long arm8_hidden_dtoul(double a);
extern unsigned int arm8_hidden_dtou(double a);
extern unsigned long arm8_hidden_ftoul(float a);
extern unsigned int arm8_hidden_ftou(float a);
extern unsigned long arm8_hidden_udiv(unsigned long a, unsigned long b);
extern long arm8_hidden_div(long a, long b);

#define FIX_SRC2_POSTFIX     \
if (tmp_src2 != -1) {\
      arm8_raw_pop(s, tmp_src2);\
    }

#define FIX_SRC2_PREFIX	     \
    int tmp_src2 = -1;\
    /*  \
     * DILL lets us use ret_reg in operators in some circumstances, so src1 \
     * or src2 might be ret_reg.  This is a problem if it's src2 and we do \
     * a call.\
     */\
    if (src2 == _a1) {\
      tmp_src2 = _a4; /* arbitrary */\
      if (src1 == tmp_src2) tmp_src2 = _a3;\
      arm8_raw_push(s, tmp_src2);\
      dill_movl(s, tmp_src2, _a1);\
      src2 = tmp_src2;\
    }

extern void arm8_mod(dill_stream s, int sign, int type_long, int dest, 
		      int src1, int src2)
{
    int return_reg;
    FIX_SRC2_PREFIX;
    if (sign == 1) {
	/* signed case */
	if (type_long) {
	    return_reg = dill_scalll(s, (void*)arm8_hidden_mod, "arm8_hidden_mod", "%l%l", src1, src2);
	    dill_movl(s, dest, return_reg);
	} else {
	    return_reg = dill_scalli(s, (void*)arm8_hidden_modi, "arm8_hidden_modi", "%i%i", src1, src2);
	    dill_movi(s, dest, return_reg);
	}
    } else {
	/* unsigned case */
	if (type_long) {
	    return_reg = dill_scalll(s, (void*)arm8_hidden_umod, "arm8_hidden_umod", "%l%l", src1, src2);
	    dill_movul(s, dest, return_reg);
	} else {
	    return_reg = dill_scallu(s, (void*)arm8_hidden_umodi, "arm8_hidden_umodi", "%u%u", src1, src2);
	    dill_movu(s, dest, return_reg);
	}
    }
    FIX_SRC2_POSTFIX;
}

extern void arm8_modi(dill_stream s, int data1, int data2, int dest, int src1, 
		      long imm)
{
    arm8_set(s, _v1, imm);
    arm8_mod(s, data1, data2, dest, src1, _v1);
}

extern void arm8_div(dill_stream s, int unsign, int junk, int dest, int src1,
		      int src2)
{
    int return_reg;
    void *routine = (void*) &arm8_hidden_div;
    if (unsign) routine = (void*) &arm8_hidden_udiv;

    FIX_SRC2_PREFIX;
    return_reg = dill_scalll(s, routine, "routine", "%l%l", src1, src2);
    FIX_SRC2_POSTFIX;
    dill_movl(s, dest, return_reg);
}

#define MUL(s,A,S,Rd,Rs,Rm) INSN_OUT(s, COND(AL)|(A&1)<<21|(S&1)<<20|RN(Rd)|RD(0)|(Rs&0xf)<<8|0x90|(Rm&0xf))

extern void arm8_mul(dill_stream s, int unsign, int junk, int dest, int src1,
		      int src2)
{
    MUL(s, 0, 0, dest, src1, src2);
}

extern void arm8_muli(dill_stream s, int unsign, int junk, int dest, int src,
		      long imm)
{
    arm8_set(s, _v1, imm);
    MUL(s, 0, 0, dest, src, _v1);
}

extern void arm8_divi(dill_stream s, int unsign, int junk, int dest, int src, 
		      long imm)
{
    arm8_set(s, _v1, imm);
    arm8_div(s, unsign, junk, dest, src,	_v1);
}

extern void
arm8_mov(dill_stream s, int type, int junk, int dest, int src)
{
    if (src == dest) return;
    switch(type) {
    case DILL_D:
	arm8_movd(s, dest, src);
	break;
    case DILL_F:
	arm8_movf(s, dest, src);
	break;
    default:
	arm8_movi(s, dest, src);
    }
}


static void
arm8_saverestore_floats(dill_stream s, int saverestore)
{
    int i;
    for (i=1; i <8; i++) {
	if (dill_mustsave(&s->p->tmp_f, i)) {
	    arm8_save_restore_op(s, saverestore, DILL_D, i);
	}
    }
}

#define CONV(x,y) ((x*100)+y)
extern void
arm8_convert(dill_stream s, int from_type, int to_type, 
	      int dest, int src)
{
    from_type &= 0xf;
    to_type &= 0xf;
    switch(CONV(from_type, to_type)) {
    case CONV(DILL_I, DILL_L):
    case CONV(DILL_I, DILL_U):
    case CONV(DILL_I,DILL_UL):
    case CONV(DILL_UL,DILL_I):
    case CONV(DILL_UL,DILL_U):
    case CONV(DILL_L,DILL_U):
    case CONV(DILL_U,DILL_UL):
    case CONV(DILL_U,DILL_L):
    case CONV(DILL_L,DILL_I):
    case CONV(DILL_UL,DILL_L):
    case CONV(DILL_L,DILL_UL):
    case CONV(DILL_P,DILL_UL):
    case CONV(DILL_UL,DILL_P):
    case CONV(DILL_U,DILL_I):
	if(src == dest) return;
	arm8_movi(s, dest,src);
	break;
    case CONV(DILL_F,DILL_D):
        arm8_fproc2(s, 0b0111, 0, 1, dest, src);  /* fcvtds */
	break;
    case CONV(DILL_F,DILL_L):
    case CONV(DILL_F,DILL_I):
    case CONV(DILL_F,DILL_S):
    case CONV(DILL_F,DILL_C):
	arm8_fproc2(s, 0b1101, 0, 1, src, src);  /* ftosis */
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|1<<20|FN(src>>1)|N(src)|RD(dest)|cp_num(0xa)|1<<4);/*fmrs*/
	break;
    case CONV(DILL_F,DILL_U):
    case CONV(DILL_F,DILL_US):
    case CONV(DILL_F,DILL_UC):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scallu(s, (void*)arm8_hidden_ftou, "arm8_hidden_ftou", "%f", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_UL, 0, dest, ret);
	}
	break;
	/* fallthrough */
    case CONV(DILL_F,DILL_UL):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scallul(s, (void*)arm8_hidden_ftoul, "arm8_hidden_ftoul", "%f", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_UL, 0, dest, ret);
	}
	break;
    case CONV(DILL_D,DILL_F):
	arm8_fproc2(s, 0b0111, 1, 1, dest, src);  /* fcvtds */
	break;
    case CONV(DILL_D,DILL_L):
    case CONV(DILL_D,DILL_I):
    case CONV(DILL_D,DILL_S):
    case CONV(DILL_D,DILL_C):
	arm8_fproc2(s, 0b1101, 1, 1, src, src);  /* ftosid */
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|1<<20|FN(src>>1)|N(src)|RD(dest)|cp_num(0xa)|1<<4);/*fmsr*/
	break;
    case CONV(DILL_D,DILL_U):
    case CONV(DILL_D,DILL_US):
    case CONV(DILL_D,DILL_UC):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scallu(s, (void*)arm8_hidden_dtou, "arm8_hidden_dtou", "%d", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_U, 0, dest, ret);
	}
	break;
    case CONV(DILL_D,DILL_UL):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scallul(s, (void*)arm8_hidden_dtoul, "arm8_hidden_dtoul", "%d", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_UL, 0, dest, ret);
	}
	break;
    case CONV(DILL_I,DILL_D):
    case CONV(DILL_L,DILL_D):
    case CONV(DILL_C,DILL_D):
    case CONV(DILL_S,DILL_D):
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(dest>>1)|N(dest)|RD(src)|cp_num(0xa)|1<<4);/*fmsr*/
	arm8_fproc2(s, 0x8, 1, 1, dest, dest);  /* fsitod */
	break;
    case CONV(DILL_U,DILL_D):
    case CONV(DILL_UL,DILL_D): 
    case CONV(DILL_UC,DILL_D):
    case CONV(DILL_US,DILL_D):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scalld(s, (void*)arm8_hidden_ultod, "arm8_hidden_ultod", "%l", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_D, 0, dest, ret);
	}
	break;
    case CONV(DILL_I,DILL_F):
    case CONV(DILL_L,DILL_F):
    case CONV(DILL_C,DILL_F):
    case CONV(DILL_S,DILL_F):
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(dest>>1)|N(dest)|RD(src)|cp_num(0xa)|1<<4);/*fmsr*/
	arm8_fproc2(s, 0x8, 0, 1, dest, dest);  /* fsitos */
	break;
    case CONV(DILL_U,DILL_F):
    case CONV(DILL_UL,DILL_F):
    case CONV(DILL_UC,DILL_F):
    case CONV(DILL_US,DILL_F):
        {
	    int ret;
	    arm8_saverestore_floats(s, 0);
	    ret = dill_scallf(s, (void*)arm8_hidden_ultof, "arm8_hidden_ultof", "%l", src);
	    arm8_saverestore_floats(s, 1);
	    arm8_mov(s, DILL_D, 0, dest, ret);
	}
	break;
    case CONV(DILL_C,DILL_UL):
    case CONV(DILL_C,DILL_L):
    case CONV(DILL_C,DILL_I):
    case CONV(DILL_C,DILL_U):
    case CONV(DILL_C, DILL_S):
    case CONV(DILL_S, DILL_C):
    case CONV(DILL_US, DILL_C):
	arm8_lshi(s, dest, src, 24);
	arm8_rshai(s, dest, dest, 24);
	break;
    case CONV(DILL_I, DILL_C):
    case CONV(DILL_U, DILL_C):
    case CONV(DILL_L, DILL_C):
    case CONV(DILL_UL, DILL_C):
    case CONV(DILL_C, DILL_UC):
    case CONV(DILL_I, DILL_UC):
    case CONV(DILL_S, DILL_UC):
    case CONV(DILL_U, DILL_UC):
    case CONV(DILL_L, DILL_UC):
    case CONV(DILL_UL, DILL_UC):
    case CONV(DILL_US, DILL_UC):
	arm8_andi(s, dest, src, 0xff);
	break;
    case CONV(DILL_S,DILL_L):
    case CONV(DILL_S,DILL_UL):
    case CONV(DILL_S,DILL_I):
    case CONV(DILL_S,DILL_U):
    case CONV(DILL_S,DILL_US):
    case CONV(DILL_US,DILL_S):
	arm8_lshi(s, dest, src, 16);
	arm8_rshai(s, dest, dest, 16);
	break;
    case CONV(DILL_C, DILL_US):
	/* signext24 - lsh24, rsha24, trunc 16 */
	arm8_lshi(s, dest, src, 24);
	arm8_rshai(s, dest, dest, 24);
	arm8_andi(s, dest, dest, 0xffff);
	break;
    case CONV(DILL_US,DILL_I):
    case CONV(DILL_US,DILL_L):
    case CONV(DILL_US,DILL_U):
    case CONV(DILL_US,DILL_UL):
    case CONV(DILL_I,DILL_S):
    case CONV(DILL_U,DILL_S):
    case CONV(DILL_L,DILL_S):
    case CONV(DILL_UL,DILL_S):
    case CONV(DILL_I,DILL_US):
    case CONV(DILL_U,DILL_US):
    case CONV(DILL_L,DILL_US):
    case CONV(DILL_UL,DILL_US):
	arm8_lshi(s, dest, src, 16);
	arm8_rshi(s, dest, dest, 16);
	break;
    default:
	printf("Unknown case in arm convert %d\n", CONV(from_type,to_type));
    }
}

static signed char op_conds[] = {
    EQ, /* dill_beq_code */  /* signed */
    GE, /* dill_bge_code */
    GT, /* dill_bgt_code */
    LE, /* dill_ble_code */
    LT, /* dill_blt_code */
    NE, /* dill_bne_code */

    EQ, /* dill_beq_code */  /* unsigned */
    CS, /* dill_bge_code */
    HI, /* dill_bgt_code */ 
    LS, /* dill_ble_code */
    CC, /* dill_blt_code */
    NE, /* dill_bne_code */
};

#define CMF 0x4
extern void
arm8_branch(dill_stream s, int op, int type, int src1, int src2, int label)
{
    switch(type) {
    case DILL_D:
    case DILL_F:
	arm8_fproc2(s, 0b0100, (type == DILL_D), 0, src1, src2);  /* fcmps */
	INSN_OUT(s, COND(AL)|0b111011110001<<16|0b1111101000010000); /* fmstat */
	dill_mark_branch_location(s, label);
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x5)|/*disp */0);/* b*/
	break;
	break;
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    default:
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(CMP)|S(1)|RN(src1)|RD(0)|RM(src2));
	dill_mark_branch_location(s, label);
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x5)|/*disp */0);/* b*/
    }
    /*    arm8_nop(s);*/
}
extern void
arm8_compare(dill_stream s, int op, int type, int dest, int src1, int src2)
{
    int setcc = 0;
    if (op == CMP) setcc = 1;
    arm8_set(s, dest, 0);
    switch(type) {
    case DILL_D:
    case DILL_F:
	arm8_fproc2(s, 0b0100, (type == DILL_D), 0, src1, src2);  /* fcmps */
	INSN_OUT(s, COND(AL)|0b111011110001<<16|0b1111101000010000); /* fmstat */
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x0)|OPCODE(MOV)|S(setcc)|RN(src1)|RD(dest)|IMM(1, 0));
	break;
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    default:
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(CMP)|S(1)|RN(src1)|RD(0)|RM(src2));
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x0)|OPCODE(MOV)|S(setcc)|RD(dest)|IMM(1, 0));
    }
    /*    arm8_nop(s);*/
}

extern void 
arm8_jal(dill_stream s, int return_addr_reg, int target)
{

}

extern void 
arm8_jump_to_label(dill_stream s, unsigned long label)
{
    dill_mark_branch_location(s, label);
    INSN_OUT(s, COND(AL)|CLASS(5)|(1<<24)/*link*/);
}

extern void arm8_jump_to_reg(dill_stream s, unsigned long reg)
{
    arm8_dproc(s, MOV, 0, _link, _pc, _pc);
    arm8_dproc(s, MOV, 0, _pc, reg, reg);
}

extern void arm8_jump_to_imm(dill_stream s, void * imm)
{

}

static void internal_push(dill_stream s, int type, int immediate, 
			  void *value_ptr)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    struct arg_info arg;
    int real_offset;

    arg.is_immediate = immediate;
    switch(type) {
    case DILL_C: case DILL_S:  case DILL_I: case DILL_L:
	arg.type = DILL_L;
	break;
    case DILL_UC: case DILL_US: case DILL_U: case DILL_UL:
	arg.type = DILL_UL;
	break;
    default:
	arg.type = type;
    }
	
    switch(arg.type) {
    case DILL_C: case DILL_S:  case DILL_I: case DILL_L:
    case DILL_UC: case DILL_US: case DILL_U: case DILL_UL: case DILL_P:
        if (ami->next_core_register < _r4) {
	    arg.is_register = 1;
	    arg.in_reg = ami->next_core_register;
	    arg.out_reg = ami->next_core_register;
	    ami->next_core_register++;
	} else {
	    arg.is_register = 0;
	}
	break;
    case DILL_D:
      if (ami->varidiac_call || (ami->hard_float == 0)) {
	    if (ami->next_core_register % 2) {
		/* double is only even regs, skip one */
		ami->next_core_register++;
		ami->cur_arg_offset += 4;
	    }
	} else {
	    if (ami->next_float_register % 2) {
		/* double is only even regs, skip one */
		ami->next_float_register++;
		ami->cur_arg_offset += 4;
	    }
	}
	/* falling through */
    case DILL_F: 
        if (ami->varidiac_call || (ami->hard_float == 0)) {
	    if (ami->next_core_register < _r4) {
		arg.is_register = 1;
		arg.in_reg = ami->next_core_register;
		arg.out_reg = ami->next_core_register;
		ami->next_core_register++;
		if (arg.type == DILL_D) ami->next_core_register++;/* two, or split */
	    } else {
		arg.is_register = 0;
	    }
	} else {
	    if (ami->next_float_register < _f16) {
		arg.is_register = 1;
		arg.in_reg = ami->next_float_register;
		arg.out_reg = ami->next_float_register;
		ami->next_float_register++;
	    } else {
		arg.is_register = 0;
	    }
	}
	break;
    default:
	assert(0);
    }
    arg.offset = ami->cur_arg_offset;
    ami->cur_arg_offset += 
	roundup(type_info[(int)arg.type].size, ami->stack_align);
    real_offset = arg.offset - 4*4; /* first 16 bytes in regs */

    if (ami->cur_arg_offset > ami->max_arg_size) {
	ami->max_arg_size = ami->cur_arg_offset;
    }
    if (arg.is_register == 0) {
	/* store it on the stack only */
	if (arg.is_immediate) {
	    if (type != DILL_D) {
		if (type == DILL_F) {
		    union {
		      float f;
		      int i;
		    } u;
		    u.f = (float) *(double*)value_ptr;
		    arm8_set(s, _v1, u.i);
		} else {
		    arm8_set(s, _v1, *(long*)value_ptr);
		}
		arm8_pstorei(s, arg.type, 0, _v1, _sp, real_offset);
	    } else {
		arm8_set(s, _v1, *(int*)value_ptr);
		arm8_pstorei(s, DILL_I, 0, _v1, _sp, real_offset);
		arm8_set(s, _v1, *(((int*)value_ptr)+1));
		arm8_pstorei(s, DILL_I, 0, _v1, _sp, real_offset+4);
	    }		
	} else {
	    arm8_pstorei(s, arg.type, 0, *(int*)value_ptr, _sp, real_offset);
	}
    } else {
	if ((type != DILL_F) && (type != DILL_D)) {
	    if (arg.is_immediate) {
		arm8_set(s, arg.out_reg, *(long*)value_ptr);
	    } else {
		arm8_mov(s, type, 0, arg.out_reg, *(int*) value_ptr);
	    }
	} else {
	    if (ami->varidiac_call || (ami->hard_float == 0)) {
		union {
		    float f;
		    int i;
		} a;
		union {
		    double d;
		    long l;
		    int i[2];
		} b;
		arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
		if (arg.is_immediate) {
		    if (type == DILL_F) {
			a.f = *(double*)value_ptr;
			arm8_set(s, arg.out_reg, a.i);
		    } else {
			b.d =  *(double*)value_ptr;
			arm8_set(s, arg.out_reg, b.i[0]);
			arm8_set(s, arg.out_reg+1, b.i[1]);
		    }
		} else {
		    switch(type) {
		    case DILL_D:
			INSN_OUT(s, COND(AL)| 0b11000101<<20|FN(arg.out_reg+1)|RD(arg.out_reg)|cp_num(0xb)|1<<4|(*(int*)value_ptr)>>1);/*fmsrrd*/
			break;
		    case DILL_F:
			INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|1<<20|FN(*(int*)value_ptr>>1)|N(*(int*)value_ptr)|RD(arg.out_reg)|cp_num(0xa)|1<<4);/*fmrs*/
			arm8_movf(s, arg.out_reg, *(int*)value_ptr);
			break;
		    default:
			assert(0);
		    }
		}
	    } else {
		if (arg.is_immediate) {
		    arm8_setf(s, type, 0, arg.out_reg, *(double*)value_ptr);
		} else {
		    switch(type) {
		    case DILL_D:
			arm8_movd(s, arg.out_reg, *(int*)value_ptr);
			break;
		    case DILL_F:
			arm8_movf(s, arg.out_reg, *(int*)value_ptr);
			break;
		    default:
			assert(0);
		    }
		}
	    }
	}
    }
}

static void push_init(dill_stream s)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    ami->cur_arg_offset = 0;
    ami->next_core_register = _r0;
    ami->next_float_register = _f0;
    ami->varidiac_call = 0;
}

extern void arm8_push(dill_stream s, int type, int reg)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    if ((type == DILL_V) && (reg <= -1)) {
	push_init(s);
	if (reg <= -2) {
	    ami->varidiac_call = 1;
	}
    } else {
	internal_push(s, type, 0, &reg);
    }
}

extern void arm8_pushi(dill_stream s, int type, long value)
{
    internal_push(s, type, 1, &value);
}

extern void arm8_pushfi(dill_stream s, int type, double value)
{
    internal_push(s, type, 1, &value);
}

extern void arm8_pushpi(dill_stream s, int type, void *value)
{
    internal_push(s, type, 1, &value);
}

extern int arm8_calli(dill_stream s, int type, void *xfer_address, const char *name)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    int caller_side_ret_reg = _a1;
    (void) name;
    /* save temporary registers */
#ifdef HOST_ARM8
    dill_mark_call_location(s, name, xfer_address);
    INSN_OUT(s, COND(AL)|CLASS(5)|(1<<24)/*link*/);
#else
    /* killing r12 here */
    arm8_set(s, _r12, (long) xfer_address);
    /* blx (register) */
    INSN_OUT(s, 0xe12fff30|RM(_r12));
#endif
    /* restore temporary registers */
    if ((type == DILL_D) || (type == DILL_F)) {
        if (ami->hard_float == 0) {
	    /* move _r0 to _f0 */
	    if (type == DILL_F) {
	        int src = _r0;
		int dest = _f0;
		INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(dest>>1)|N(dest)|RD(src)|cp_num(0xa)|1<<4);/*fmsr*/
	    } else {
	      INSN_OUT(s, COND(AL)| CLASS(6)|OPCODE(2)|FN(_r1)|RD(_r0)|cp_num(0xb)|1<<4|(_f0>>1)|(_f0&1)<<5);/*fmdrr*/
	    }
	}
	caller_side_ret_reg = _f0;
    }
    push_init(s);
    return caller_side_ret_reg;
}

extern int arm8_callr(dill_stream s, int type, int src)
{
    int caller_side_ret_reg = _a1;

    arm8_dproc(s, MOV, 0, _link, _pc, _pc);
    arm8_dproc(s, MOV, 0, _pc, src, src);

    /* restore temporary registers */
    if ((type == DILL_D) || (type == DILL_F)) {
	caller_side_ret_reg = _f0;
    }
    push_init(s);
    return caller_side_ret_reg;
}

extern void
arm8_branchi(dill_stream s, int op, int type, int src, long imm, int label)
{
    switch(type) {
    case DILL_F:
    case DILL_D:
	fprintf(stderr, "Shouldn't happen\n");
	break;
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    default:
	arm8_dproci(s, CMP, 0, 0/*dest*/, src, imm);
	dill_mark_branch_location(s, label);
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x5)|/*disp */0);/* b*/
    }
}

extern void
arm8_comparei(dill_stream s, int op, int type, int dest, int src, long imm)
{
    int setcc = 0;
    if (op == CMP) setcc = 1;
    arm8_set(s, dest, 0);
    switch(type) {
    case DILL_F:
    case DILL_D:
	fprintf(stderr, "Shouldn't happen\n");
	break;
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    default:
	arm8_dproci(s, CMP, 0, 0/*dest*/, src, imm);
	INSN_OUT(s, COND(op_conds[op])|CLASS(0x0)|OPCODE(MOV)|S(setcc)|RD(dest)|IMM(1, 0));
    }
}

static void
arm8_simple_ret(dill_stream s)
{
    dill_mark_ret_location(s);
    INSN_OUT(s, COND(AL)|CLASS(4)|1<<24/*p*/|1<<20/*l*/|RN(_r11)|1<<_r11|1<<_sp|1<<_pc);
    arm8_nop(s);  /* ldmea may slide back here if we have to restore floats */
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
    arm8_nop(s);
}

extern void arm8_ret(dill_stream s, int data1, int data2, int src)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    switch (data1) {
    case DILL_C:
    case DILL_UC:
    case DILL_S:
    case DILL_US:
    case DILL_I:
    case DILL_U:
    case DILL_L:
    case DILL_UL:
    case DILL_P:
	if (src != _a1) arm8_movi(s, _a1, src);
	break;
    case DILL_F:
        if (ami->hard_float) {
	    if (src != _f0) arm8_movf(s, _f0, src);
	} else {
	    INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|1<<20|FN(src>>1)|N(src)|RD(_r0)|cp_num(0xa)|1<<4);/*fmrs*/
	}
	break;
    case DILL_D:
        if (ami->hard_float) {
	    if (src != _f0) arm8_movd(s, _f0, src);
	} else {
	    INSN_OUT(s, COND(AL)| CLASS(6)|OPCODE(2)|r(1)|FN(_r1)|RD(_r0)|cp_num(0xb)|1<<4|(src>>1)|(src&1)<<5);/*fmrrd*/
	}
	break;
    }
    arm8_simple_ret(s);
}

extern void arm8_retf(dill_stream s, int data1, int data2, double imm)
{
    union {
	float f;
	int i;
    } a;
    union {
	double d;
	int i[2];
	long l;
    } b;
  
    switch(data1) {
    case DILL_F:
	a.f = imm;
	arm8_set(s, _r0, a.i);
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(_f0>>1)|N(_f0)|RD(_r0)|cp_num(0xa)|1<<4);/*fmsr*/
	break;
    case DILL_D:
	b.d = imm;
	arm8_set(s, _r0, b.i[0]);
	arm8_set(s, _r1, b.i[1]);
	INSN_OUT(s, COND(AL)| CLASS(7)|OPCODE(0)|0<<20|FN(_f0>>1)|N(_f0)|RD(_r0)|cp_num(0xa)|1<<4);/*fmsr*/
	break;
    }
}

extern void arm8_reti(dill_stream s, int data1, int data2, long imm)
{
    switch (data1) {
    case DILL_C:
    case DILL_UC:
    case DILL_S:
    case DILL_US:
    case DILL_I:
    case DILL_U:
    case DILL_L:
    case DILL_UL:
    case DILL_P:
	arm8_set(s, _a1, imm);
	break;
    case DILL_F:
    case DILL_D:
	break;/* no return immediate of floats */
    }
    arm8_simple_ret(s);
}

static void
arm8_data_link(dill_stream s)
{
  /*    struct branch_table *t = &s->p->branch_table;
    int i;
    for (i=0; i < t->data_mark_count; i++) {
	int label = t->data_marks[i].label;
	void *label_addr = t->label_locs[label] + (char*)s->p->code_base;
	*t->data_marks[i].addr = (long) label_addr;
	}*/
}

static void
arm8_branch_link(dill_stream s)
{
    struct branch_table *t = &s->p->branch_table;
    int i;

    for(i=0; i< t->branch_count; i++) {
	int label = t->branch_locs[i].label;
	int label_offset = t->label_locs[label] - t->branch_locs[i].loc;
	int *branch_addr = (int*)((char *)s->p->code_base + 
				  t->branch_locs[i].loc);
	/* compensate for arm PC lookahead */
	label_offset = label_offset - 8;
        /* div addr diff by 4 for arm offset value */
	label_offset = label_offset >> 2;  
	*branch_addr &= 0xff000000;
	*branch_addr |= (label_offset & 0xffffff);
    }
}

/*
 * on ARM, we need a procedure linkage table to manage 
 * calls to DLLs in an address range that is typically more than 26 bits
 * away from malloc'd memory.  We emit a PLT that is basically a set_reg, 
 * then jump_through_reg for each routine.  Later, during call linkage, 
 * we'll call to the PLT entry rather than directly to the routine.
 * If creating a package, emit a PLT entry for every routine because we 
 * don't know what the offsets will be when the code is stitched together.
 */
extern void
arm8_PLT_emit(dill_stream s, int package)
{
    call_t *t = &s->p->call_table;
    int i;

    for(i=0; i< t->call_count; i++) {
	int *call_addr = (int*) ((unsigned long)s->p->code_base + 
				 t->call_locs[i].loc);
	long call_offset = (unsigned long)t->call_locs[i].xfer_addr - 
	    (unsigned long)call_addr;
        /* div addr diff by 4 for arm offset value */
	call_offset = call_offset >> 2;
	call_offset = call_offset >> 24;
#ifdef HOST_ARM8
	if (((call_offset != 0) && (call_offset != -1)) || package) {
	    t->call_locs[i].mach_info = (void*)
		((long)s->p->cur_ip - (long)s->p->code_base);
	    arm8_pldsti(s, DILL_P, 1, _v1, _pc, 0);
	    arm8_movi(s, _pc, _v1);
	    arm8_nop(s);  /* place where addr will be loaded */
	}
#endif
    }
}

extern void arm8_rt_call_link(char *code, call_t *t);

static void
arm8_call_link(dill_stream s)
{
    arm8_rt_call_link(s->p->code_base, &s->p->call_table);
}


#ifndef CLEAR_CACHE_DEFINED
extern void __clear_cache(void *, void *);
#endif

static void
arm8_flush(void *base, void *limit)
{
#if defined(HOST_ARM8) || defined(HOST_ARM7)
  __clear_cache(base, limit);
#endif
}    

static void
arm8_emit_save(dill_stream s)
{
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    void *save_ip = s->p->cur_ip;
    int ar_size = ami->act_rec_size + ami->max_arg_size;
    int float_count = 0;
    int int_count = 3;  /* fp, ip, lr */

    int reg;
    int mask = 0;
    ret_t *t = &s->p->ret_table;
    int i;

    ar_size += 14 * 4 /* int save */ + 8 * 3 * 4 /* float save */;
    
    ar_size = roundup(ar_size, 8);
    
    switch(ami->max_arg_size) {
    case 0: case 4:
	mask |= 1<<_a3;
    case 8:
	mask |= 1<<_a3;
    case 12:
	mask |= 1<<_a4;
    default:
	/* save nothing */
	break;
    }
    mask |= 1<< _v1;
    for (reg = _v2; reg <= _v7; reg++) {
	if (dill_wasused(&s->p->tmp_i, reg)) {
	    mask |= (1<<reg);
	    int_count++;
	}
    }
    for (reg = _f16; reg <= _f30; reg++) {
	if (dill_wasused(&s->p->tmp_f, reg)) {
	    float_count = reg - _f16 + 2;
	}
    }
    s->p->cur_ip = (char*)s->p->code_base + ami->save_insn_offset - 16;
    INSN_OUT(s, COND(AL)|CLASS(4)|1<<24/*p*/|RN(_sp)| mask|1<<_r11|1<<_r12|1<<_link|1<<_pc;);
    s->p->cur_ip = (char*)s->p->code_base + ami->save_insn_offset - 8;
    if (float_count > 0) {
	INSN_OUT(s, COND(AL)|CLASS(6)|1<<24|1<<21|RN(_sp)|(_f8)<<12|0b1011<<8|float_count); /*sfm*/
    } else {
	arm8_nop(s);
    }
    s->p->cur_ip = (char*)s->p->code_base + ami->save_insn_offset;
    arm8_savei(s, -ar_size);

    for(i=0; i< t->ret_count; i++) {
	s->p->cur_ip = (char*)((char *)s->p->code_base + t->ret_locs[i]);
	arm8_dproci(s, ADD, 0, _sp, _sp, ar_size);
	if (float_count > 0) {
	    int offset = 4 * 12 + 14*4 - 4;
	    INSN_OUT(s, COND(AL)|CLASS(6)|1<<23|1<<21|1<<20|RN(_sp)|(_f8)<<12|0b1011<<8|float_count); /*lfm*/
	    arm8_dproci(s, ADD, 0, _sp, _sp, 4*14);

	} else {
	    arm8_dproci(s, ADD, 0, _sp, _sp, ar_size + 4*14);
	}
        INSN_OUT(s, COND(AL)|CLASS(4)|1<<24/*p*/|1<<20/*l*/|RN(_r11)|1<<_r11|1<<_sp|1<<_pc|mask);
    }
    s->p->fp = (char*)s->p->code_base + 12; /* skip 3 swinv */
    s->p->cur_ip = save_ip;
}
    
extern void
arm8_end(s)
dill_stream s;
{
    arm8_nop(s);
    arm8_simple_ret(s);
    arm8_PLT_emit(s, 0);   /* must be done before linking */
    arm8_branch_link(s);
    arm8_call_link(s);
    arm8_data_link(s);
    arm8_emit_save(s);
    arm8_flush(s->p->code_base, s->p->code_limit);
}

extern void
arm8_package_end(s)
dill_stream s;
{
    arm8_nop(s);
    arm8_simple_ret(s);
    arm8_PLT_emit(s, 1);
    arm8_branch_link(s);
    arm8_emit_save(s);
}

extern void *
arm8_clone_code(s, new_base, available_size)
dill_stream s;
void *new_base;
int available_size;
{
    int size = dill_code_size(s);
    if (available_size < size) {
	return NULL;
    }
    void *old_base = s->p->code_base;
    void *native_base = s->p->code_base;
    if (native_base == NULL) native_base = s->p->native.code_base;
    memcpy(new_base, native_base, size);
    s->p->code_base = new_base;
    s->p->cur_ip = new_base + size;
    s->p->fp = new_base;
    arm8_branch_link(s);
    arm8_call_link(s);
    arm8_data_link(s);
    arm8_flush(new_base, (void*)((long)new_base + size));
    s->p->code_base = old_base;
    s->p->cur_ip = old_base + size;
    s->p->fp = old_base;
    while (*(int*)new_base == 0xFF000000) {
	/* skip UNIMPs */
	new_base = (void*)((long) new_base + 4);
    }
    return new_base;
}

extern void
arm8_pset(dill_stream s, int type, int junk, int dest, long imm)
{
    arm8_set(s, dest, imm);
}	

extern void
arm8_setp(dill_stream s, int type, int junk, int dest, void* imm)
{
    arm8_pset(s, DILL_L, 0, dest, (long)imm);
}	

extern void
arm8_setf(dill_stream s, int type, int junk, int dest, double imm)
{
    union {
	float f;
	int i;
    } a;
    union {
	double d;
	long l;
	int i[2];
    } b;
    arm8_mach_info ami = (arm8_mach_info) s->p->mach_info;
    if (type == DILL_F) {
	a.f = (float) imm;
	arm8_set(s, _v1, a.i);
	arm8_movi2f(s, dest, _v1);
    } else {
	b.d = imm;
	arm8_set(s, _v1, b.i[0]);
	arm8_pstorei(s, DILL_I, 0, _v1, _fp, ami->conversion_word);
	arm8_set(s, _v1, b.i[1]);
	arm8_pstorei(s, DILL_I, 0, _v1, _fp, ami->conversion_word+4);
	arm8_ploadi(s, DILL_D, 0, dest, _fp, ami->conversion_word);
    }
}	


extern void
arm8_set(s, r, val)
dill_stream s;
int r;
long val;
{
#if defined(HOST_ARM7)
    /* movw */
    INSN_OUT(s, COND(AL)|3<<24|((val>>12)&0xf)<<16| RD(r)| val&0xfff);
    if ((val & 0xffff0000) != 0) {
	int imm = (val >> 16) & 0xffff;
	/* movt */
	INSN_OUT(s, COND(AL)|3<<24|4<<20|((imm>>12)&0xf)<<16| RD(r)| imm&0xfff);
    }
#else
    arm8_dproci(s, MOV, 0, r, 0, val & 0xff);
    if ((val & 0xff00) != 0) {
	int imm = (val >> 8) & 0xff;
	/* or in the byte */
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(ORR)|S(0)|RN(r)|RD(r)|IMM(imm, 8));
    }
    if ((val & 0xff0000) != 0) {
	int imm = (val >> 16) & 0xff;
	/* or in the byte */
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(ORR)|S(0)|RN(r)|RD(r)|IMM(imm, 16));
    }
    if ((val & 0xff000000) != 0) {
	int imm = (val >> 24) & 0xff;
	/* or in the byte */
	INSN_OUT(s, COND(AL)|CLASS(0x0)|OPCODE(ORR)|S(0)|RN(r)|RD(r)|IMM(imm, 24));
    }
#endif
}

#define bit_R(x) ((unsigned long)1<<x)

extern void
arm8_reg_init(dill_stream s)
{
    s->p->var_i.init_avail[0] = 0;
    s->p->var_i.members[0] = s->p->var_i.init_avail[0];
    s->p->tmp_i.init_avail[0] = (bit_R(_v2)|bit_R(_v3)|bit_R(_v4)|
				 bit_R(_v5)|bit_R(_v6)|bit_R(_v7));
    s->p->tmp_i.members[0] = s->p->tmp_i.init_avail[0] | bit_R(_v1) |
	(bit_R(_a1)|bit_R(_a2)|bit_R(_a3)|bit_R(_a4));
    s->p->var_f.init_avail[0] = 0;
    s->p->var_f.members[0] = s->p->var_f.init_avail[0];
    /* in reality, there are 32 single precision regs, overlapping 16 
     * double-precision regs.  DILL isn't quite smart enough to handle 
     * overlapping registers, so we allocate only the even ones.  _f2 if 
     * used as a double-precision corresponds to ARM register D1.
     */
    s->p->tmp_f.init_avail[0] = (bit_R(_f16)|bit_R(_f18)|bit_R(_f20)|bit_R(_f22)|
				 bit_R(_f24)|bit_R(_f26)|bit_R(_f28)|bit_R(_f30));
    s->p->tmp_f.members[0] = s->p->tmp_f.init_avail[0];
}

extern void*
gen_arm8_mach_info(s, v9)
dill_stream s;
int v9;
{
    arm8_mach_info ami = malloc(sizeof(*ami));
    if (s->p->mach_info != NULL) {
	free(s->p->mach_info);
	s->p->mach_info = NULL;
	s->p->native.mach_info = NULL;
    }
    arm8_reg_init(s);
    ami->act_rec_size = 0;
    ami->conversion_word = 0;
    ami->gp_save_offset = 0;
    ami->cur_arg_offset = 0;
    ami->next_core_register = _r0;
    ami->next_float_register = _f0;
    ami->stack_align = 4;
    ami->stack_constant_offset = 0;
    ami->fp_save_offset = ami->gp_save_offset + 8 * ami->stack_align;
    ami->fp_save_end = ami->fp_save_offset + 8 * 8;
    ami->max_arg_size = 0;
    ami->hard_float = ARM_HARD_FLOAT;
    return ami;
}

#if defined(HAVE_DIS_ASM_H) && !defined(NO_DISASSEMBLER)
/* GENERIC BINUTILS DISASSEMBLER */
#include "dis-asm.h"

#define MAXLENGTH (1<<23) /* Max length of function that can be disassembled */

extern int
arm8_init_disassembly_info(dill_stream s, void * ptr)
{
    struct disassemble_info *i = ptr;
#ifdef INIT_DISASSEMBLE_INFO_THREE_ARG
    INIT_DISASSEMBLE_INFO(*i, stdout,fprintf);
    i->endian = BFD_ENDIAN_BIG;
#else
    INIT_DISASSEMBLE_INFO(*i, stdout);
#endif
#ifdef bfd_mach_arm8_5
    i->mach = bfd_mach_arm8_5;
#elif defined (bfd_mach_arm8_4)
    i->mach = bfd_mach_arm8_4;
#elif defined (bfd_mach_arm8_3)
    i->mach = bfd_mach_arm8_3;
#endif
    if (s->p->code_base != NULL) {
	i->buffer = (bfd_byte *)s->p->code_base;
	i->buffer_vma = (bfd_vma)s->p->code_base;
    } else {
	i->buffer = (bfd_byte *)s->p->native.code_base;
	i->buffer_vma = (bfd_vma)s->p->native.code_base;
    }
    i->buffer_length = MAXLENGTH;
#ifdef HAVE_PRINT_INSN_ARM
    return 1;
#elif defined(HAVE_PRINT_INSN_LITTLE_ARM)
    return 1;
#else
    return 0;
#endif
}

extern int
arm8_print_insn(dill_stream s, void *info_ptr, void *insn)
{
#ifdef HAVE_PRINT_INSN_ARM
    return print_insn_arm((unsigned long) insn, (disassemble_info*)info_ptr);
#elif defined(HAVE_PRINT_INSN_LITTLE_ARM)
    return print_insn_little_arm((unsigned long) insn, (disassemble_info*)info_ptr);
#else
    return 0;
#endif
}
#else
extern int
arm8_init_disassembly_info(dill_stream s, void * ptr){return 0;}
extern int arm8_print_insn(dill_stream s, void *info_ptr, void *insn){return 0;}
#endif

extern void
arm8_print_reg(dill_stream s, int typ, int reg)
{
    switch(typ) {
    case DILL_C: case DILL_UC:
    case DILL_S: case DILL_US:
    case DILL_I: case DILL_U: case DILL_L: case DILL_UL:
	if (reg == _sp) {
	    printf("sp");
	    return;
	} else if (reg == _link) {
	    printf("link");
	    return;
	} else if (reg == _pc) {
	    printf("pc");
	    return;
	} else if (reg == _fp) {
	    printf("fp");
	    return;
	} else if (reg <= _r3) {
	    printf("r%d(a%d)\n", reg, reg +1);
	    return;
	} else if (reg <= _r10) {
	    printf("r%d(v%d)\n", reg, reg - 3);
	    return;
	}
	break;
    case DILL_F: case DILL_D:
	printf("F%d", reg);
	return;
    }
    printf("NoReg(%d)", reg);
}

extern int
arm8_count_insn(dill_stream s, int start, int end)
{
    return (end - start)>>2;
}
