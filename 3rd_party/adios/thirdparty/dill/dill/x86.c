#include "config.h"

#ifdef LINUX_KERNEL_MODULE
#ifndef MODULE
#define MODULE
#endif
#ifndef __KERNEL__
#define __KERNEL__
#endif
#include <linux/kernel.h>
#include <linux/module.h>
#endif

#undef NDEBUG
#include "assert.h"

#ifndef LINUX_KERNEL_MODULE
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef STDC_HEADERS
#include <stdlib.h>
#endif
#include <string.h>
#else
#include "kdill.h"
#define fprintf(fmt, args...)	printk(args)
#define printf	printk
#define malloc (void *)DAllocMM
#define free(a) DFreeMM((addrs_t) a)
#endif

/*
 * GANEV: note that we have to include "x86.h" _after_ including
 * "kdill.h" because it needs to #undef and then re-#define a few
 * symbols (namely, EAX, EBX, etc. defined in <asm/ptrace.h>
 */
#include "dill.h"
#include "dill_internal.h"
#include "x86.h"

#define MOV32 0x89
#define Mod(x) ((x)<<6)
#define RegOp(x) ((x)<<3)
#define RM(x) (x)
#define ModRM(mod,reg,rm) (Mod(mod)|RegOp(reg)|RM(rm))
#define SIB(scale, index, base) (((scale)<<6) | ((index)<<3) | (base))

#define INSN_OUT(s, i) printf("Bad opout, line %d\n", __LINE__)
#define x86_savei(s, imm) 
#define x86_seti(s, r, val) BYTE_OUT1I(s, 0xb8 + r, val);
#define x86_rshai(s, dest, src, imm) if (dest !=src) x86_movi(s, dest, src); BYTE_OUT3(s, 0xc1, ModRM(0x3, 0x7, dest), imm & 0xff);
#define x86_rshi(s, dest, src, imm) if (dest !=src) x86_movi(s, dest, src); BYTE_OUT3(s, 0xc1, ModRM(0x3, 0x5, dest), imm & 0xff);
#define x86_lshi(s, d, src, imm) if (d !=src) x86_movi(s, d, src); BYTE_OUT3(s, 0xc1, ModRM(0x3, 0x4, d), imm & 0xff);
#define x86_andi(s, dest, src, imm) x86_arith3i(s, 0x4, 0x2, dest, src, imm)
#define x86_movi(s, dest, src) 	do { if (src != dest) BYTE_OUT2(s, MOV32, ModRM(0x3, src, dest)); } while (0)

#define x86_push_reg(s, src) BYTE_OUT1(s, 0x50+src)
#define x86_pop_reg(s, src) BYTE_OUT1(s, 0x58+src)
#define x86_nop(c) BYTE_OUT1(s, 0x90)

#define IREG 0
#define FREG 1

#define _temp_reg EAX
#define _frame_reg EBP
#define roundup(a,b) ((a + (b-1)) & (-b))

extern void
x86_sse_ploadi(dill_stream s, int type, int junk, int dest, int src, long offset);

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
    { sizeof(double), sizeof(double), FREG},  /* D */
    { 0, 8, IREG}, /* V */
    { -1, 8, IREG}, /* B */
    { sizeof(long), sizeof(long), IREG}, /* EC */
};

int x86_type_align[] = {
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
        4, /* V */
        4, /* B */
        sizeof(char*), /* EC */
};

int x86_type_size[] = {
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
        4, /* V */
        0, /* B */
        sizeof(char*), /* EC */
};

static void x86_pmov(dill_stream s, int typ, int dest, int src);
#define x86_movd(s, dest, src)  x86_pmov(s, DILL_D, dest, src)
#define x86_movf(s, dest, src)  x86_pmov(s, DILL_F, dest, src)

static void
x86_pmov(dill_stream s, int typ, int dest, int src)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    if (src != dest) {
	switch (typ) {
	case DILL_D: case DILL_F:
	    if (smi->generate_SSE) {
		BYTE_OUT3(s, 0x0f, 0x29, ModRM(0x3, src, dest));
	    }
	    break;
	default:
	    BYTE_OUT2(s, MOV32, ModRM(0x3, src, dest));
	    break;
	}
    }
}

extern int
x86_local(dill_stream s, int type)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    s->p->used_frame++;
    smi->act_rec_size += roundup(type_info[type].size, smi->stack_align);
    return (-smi->act_rec_size) + smi->stack_constant_offset;
}

extern int
x86_localb(dill_stream s, int size)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    s->p->used_frame++;
    smi->act_rec_size = roundup(smi->act_rec_size, size);

    smi->act_rec_size += roundup(size, smi->stack_align);
    return (-smi->act_rec_size) + smi->stack_constant_offset;
}

extern int x86_local_op(dill_stream s, int flag, int val)
{
    int size = val;
    if (flag == 0) {
	size = type_info[val].size;
    }
    return x86_localb(s, size);
}	

#define BEGIN_FLOAT_SAVE 32
extern void
x86_save_restore_op(dill_stream s, int save_restore, int type, int reg)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    int offset;
    switch (type) {
    case DILL_D: case DILL_F:
	offset = reg * 8 + BEGIN_FLOAT_SAVE;
	break;
    default:
	offset = reg * smi->stack_align;
	break;
    }
    if (save_restore == 0) { /* save */
	x86_pstorei(s, type, 0, reg, _frame_reg, smi->save_base + offset + smi->stack_constant_offset);
    } else {  /* restore */
	x86_ploadi(s, type, 0, reg, _frame_reg, smi->save_base + offset + smi->stack_constant_offset);
    }
    s->p->used_frame++;
}	

static void
save_required_regs(dill_stream s, int force)
{
    if (s->p->call_table.call_count != 0) force++;

    if (force) {
	dill_andii(s, ESP, ESP, -16); /* make sure it's multiple of 16 */
    }
    /* callee is supposed to save these */
    if (force || dill_wasused(&s->p->var_i, EDI) || dill_wasused(&s->p->tmp_i, EDI)) {
	x86_push_reg(s, EDI);
    }
    if (force || dill_wasused(&s->p->var_i, ESI) || dill_wasused(&s->p->tmp_i, ESI)) {
	x86_push_reg(s, ESI);
    }
    if (force || dill_wasused(&s->p->var_i, EBX) || dill_wasused(&s->p->tmp_i, EBX)) {
	x86_push_reg(s, EBX);
    }
    if (force) {
	x86_push_reg(s, EBX);
    }
}

static int 
generate_prefix_code(dill_stream s, int force, int ar_size )
{
    int end, start = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    int i;
    arg_info_list args = s->p->c_param_args;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    if ((s->p->c_param_count > 0) || s->p->used_frame || (s->p->call_table.call_count != 0) || force) {
	x86_push_reg(s, EBP);
	x86_movi(s, EBP, ESP);
    }

    if (force || s->p->used_frame) {
	/* do local space reservation */
	dill_subii(s, ESP, ESP, ar_size);
    }

    save_required_regs(s, force);

    for (i = 0; i < s->p->c_param_count; i++) {
	if (args[i].is_register) {
	    if ((args[i].type != DILL_F) && (args[i].type != DILL_D)) {
		x86_ploadi(s, DILL_I, 0, args[i].in_reg, EBP, args[i].offset);
	    } else {
		if (smi->generate_SSE) {
		    x86_ploadi(s, args[i].type, 0, args[i].in_reg, EBP, args[i].offset);
		    dill_dealloc_specific(s, args[i].in_reg, args[i].type, DILL_TEMP);
		}
	    }
	}
    }

    end = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    return end - start;
}

static void *last_proc_ret_end = 0;

extern void
x86_proc_start(dill_stream s, char *subr_name, int arg_count, arg_info_list args,
	     dill_reg *arglist)
{
    int i, fp_arg_count, int_arg_count;
    int cur_arg_offset = 0;

    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    smi->pending_prefix = 0;

    /* leave some space */ x86_local(s, DILL_D);
    smi->conversion_word = x86_local(s, DILL_D);
    smi->fcu_word = x86_local(s, DILL_I);
    smi->save_base = x86_localb(s, 8 * 4 + /* floats */ 8 * 8);
    s->p->used_frame = 0;  /* don't count our own use of the frame */

    cur_arg_offset = 8;
    fp_arg_count = 0;
    int_arg_count = 0;
    for (i = 0; i < arg_count; i++) {
	/* at most 2 args moved into registers */
	if (((args[i].type == DILL_F) || (args[i].type == DILL_D)) &&
	    smi->generate_SSE && (fp_arg_count < 2)) {
	    /* put up to two float params in XMM reg set */
	    int reg;
	    fp_arg_count++;
	    args[i].is_register = dill_raw_getreg(s, &reg, args[i].type, DILL_VAR);
	    args[i].in_reg = reg;
	    args[i].out_reg = -1;
	} else if ((args[i].type != DILL_F) && (args[i].type != DILL_D) && 
		   (int_arg_count < 2)) {
	    int reg;
	    int_arg_count++;
	    args[i].is_register = dill_raw_getreg(s, &reg, DILL_L, DILL_VAR);
	    args[i].in_reg = reg;
	    args[i].out_reg = -1;
	} else {
	    args[i].is_register = 0;
	}
	if ((args[i].is_register)  && (arglist != NULL)) 
	    arglist[i] = args[i].in_reg;

	args[i].offset = cur_arg_offset;
	cur_arg_offset += roundup(type_info[(int)args[i].type].size, smi->stack_align);
    }

    /* make local space reservation constant big so we have a word to patch */
    /* use the nop op code so that if we don't use all of it we get nops */
    (void) generate_prefix_code(s, 1 /* force */, 0x90909090);

    smi->backpatch_offset = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);

    last_proc_ret_end = 0;
}

static void
x86_proc_ret(dill_stream s) 
{
    int force = 0;
    if (s->p->call_table.call_count != 0) force++;

    /* last thing was a ret and there's not a label here, we don't need another ret */
    if ((last_proc_ret_end == s->p->cur_ip) && !dill_is_label_mark(s))
	return;

    if (force) {
	x86_pop_reg(s, EBX);
    }
    if (force || dill_wasused(&s->p->var_i, EBX) || dill_wasused(&s->p->tmp_i, EBX)) {
	x86_pop_reg(s, EBX);
    }
    if (force || dill_wasused(&s->p->var_i, ESI) || dill_wasused(&s->p->tmp_i, ESI)) {
	x86_pop_reg(s, ESI);
    }
    if (force || dill_wasused(&s->p->var_i, EDI) || dill_wasused(&s->p->tmp_i, EDI)) {
	x86_pop_reg(s, EDI);
    }
    if ((s->p->c_param_count > 0) || s->p->used_frame || (s->p->call_table.call_count != 0)) {
	x86_movi(s, ESP, EBP);
	x86_pop_reg(s, EBP);
    }
    BYTE_OUT1(s, 0xc3);
    last_proc_ret_end = s->p->cur_ip;
}

static unsigned char ld_opcodes[] = {
    0x8a, /* DILL_C */
    0x8a, /* DILL_UC */
    0x8b, /* DILL_S */
    0x8b, /* DILL_US */
    0x8b, /* DILL_I */
    0x8b, /* DILL_U */
    0x8b, /* DILL_L */
    0x8b, /* DILL_UL */
    0x8b, /* DILL_P */
    0xd9, /* DILL_F */
    0xdd, /* DILL_D */
    0x00, /* DILL_V */
    0x00, /* DILL_B */
    0x8b, /* DILL_EC */
};

static void x86_clear(s, dest)
dill_stream s;
int dest;
{
    BYTE_OUT2(s, 0x33, ModRM(0x3, dest, dest));  /* xor dest, dest */
}

extern void
x86_ploadi(dill_stream s, int type, int force_8087, int dest, int src, long offset)
{
    unsigned char opcode = ld_opcodes[type];
    int tmp_dest = dest;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((smi->generate_SSE && !force_8087) &&
	((type == DILL_F) || (type == DILL_D))) {
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    tmp_dest = dest;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    tmp_dest = dest;
	    break;
	}
	if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT5(s, opcode, 0x0f, 0x10, ModRM(0x1, tmp_dest, src), offset & 0xff);
	} else {
	    BYTE_OUT4I(s, opcode, 0x0f, 0x10, ModRM(0x2, tmp_dest, src), offset);
	}
	return;
    }
    switch (type) {
    case DILL_F:
	opcode = 0xd9;
	tmp_dest = dest = 0;
	break;
    case DILL_D:
	opcode = 0xdd;
	tmp_dest = dest = 0;
	break;
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, load to eax */
	    tmp_dest = EAX;
	}
	if (type == DILL_UC) {
	    /* clear high bytes */
	    if (src == tmp_dest) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	break;
    case DILL_S: case DILL_US:
	if (type == DILL_US) {
	    /* clear high bytes */
	    if (src == tmp_dest) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	BYTE_OUT1(s, 0x66);
	break;
    case DILL_L: case DILL_UL: case DILL_P:
    /* fall through */
    default:
	break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    if (src == ESP) {
	if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT4(s, opcode, ModRM(0x1, tmp_dest, ESP), SIB(0, 4, ESP), offset & 0xff);
	} else {
	    BYTE_OUT3I(s, opcode, ModRM(0x2, tmp_dest, src), SIB(0, 4, ESP), offset);
	}
    } else {
	if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT3(s, opcode, ModRM(0x1, tmp_dest, src), offset & 0xff);
	} else {
	    BYTE_OUT2I(s, opcode, ModRM(0x2, tmp_dest, src), offset);
	}
    }
    switch(type){
    case DILL_C:
	x86_lshi(s, dest, tmp_dest, 24);
	x86_rshi(s, dest, dest, 24);
	break;
    case DILL_S:
	x86_lshi(s, dest, tmp_dest, 16);
	x86_rshi(s, dest, dest, 16);
	break;
    case DILL_UC: case DILL_US:
	if (dest != tmp_dest)
	    x86_movi(s, dest, tmp_dest);
	break;
    }
}

extern void
x86_sse_ploadi(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    unsigned char opcode = ld_opcodes[type];
    int tmp_dest = dest;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((type == DILL_F) || (type == DILL_D)) {
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    tmp_dest = dest;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    tmp_dest = dest;
	    break;
	}
	if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT5(s, opcode, 0x0f, 0x10, ModRM(0x1, tmp_dest, src), offset & 0xff);
	} else {
	    BYTE_OUT4I(s, opcode, 0x0f, 0x10, ModRM(0x2, tmp_dest, src), offset);
	}
	return;
    }
    switch (type) {
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, load to eax */
	    tmp_dest = EAX;
	}
	if (type == DILL_UC) {
	    /* clear high bytes */
	    if (src == tmp_dest) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	break;
    case DILL_S: case DILL_US:
	if (type == DILL_US) {
	    /* clear high bytes */
	    if (src == tmp_dest) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	BYTE_OUT1(s, 0x66);
	break;
    case DILL_L: case DILL_UL: case DILL_P:
    /* fall through */
    default:
	break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    if (((long)offset <= 127) && ((long)offset > -128)) {
	BYTE_OUT3(s, opcode, ModRM(0x1, tmp_dest, src), offset & 0xff);
    } else {
	BYTE_OUT2I(s, opcode, ModRM(0x2, tmp_dest, src), offset);
    }
    switch(type){
    case DILL_C:
	x86_lshi(s, dest, tmp_dest, 24);
	x86_rshi(s, dest, dest, 24);
	break;
    case DILL_S:
	x86_lshi(s, dest, tmp_dest, 16);
	x86_rshi(s, dest, dest, 16);
	break;
    case DILL_UC: case DILL_US:
	if (dest != tmp_dest)
	    x86_movi(s, dest, tmp_dest);
	break;
    }
}

extern void
x86_pload(dill_stream s, int type, int force_8087, int dest, int src1, int src2)
{
    unsigned char opcode = ld_opcodes[type];
    int tmp_dest = dest;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((smi->generate_SSE && !force_8087) &&
	((type == DILL_F) || (type == DILL_D))) {
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    tmp_dest = dest;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    tmp_dest = dest;
	    break;
	}
	BYTE_OUT5(s, opcode, 0x0f, 0x10, ModRM(0x0, dest, 0x4), SIB(0, src1, src2));
	return;
    }
    switch(type){
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, load to eax */
	    tmp_dest = EAX;
	}
	if (type == DILL_UC) {
	    /* clear high bytes */
	    if ((src1 == tmp_dest) || (src2 == tmp_dest)) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	break;
    case DILL_S: case DILL_US:
	if (type == DILL_US) {
	    /* clear high bytes */
	    if ((src1 == tmp_dest) || (src2 == tmp_dest)) {
		/* don't destroy source */
		tmp_dest = EAX;
	    }
	    x86_clear(s, tmp_dest);
	}
	BYTE_OUT1(s, 0x66);
	break;
    case DILL_F:
	if (smi->pending_prefix != 0) {
	    BYTE_OUT1(s, smi->pending_prefix);
	    smi->pending_prefix = 0;
	}
	BYTE_OUT3(s, 0xd9, ModRM(0x0, 0x0, 0x4), SIB(0, src1, src2));
	return;
	break;
    case DILL_D:
	if (smi->pending_prefix != 0) {
	    BYTE_OUT1(s, smi->pending_prefix);
	    smi->pending_prefix = 0;
	}
	BYTE_OUT3(s, 0xdd, ModRM(0x0, 0x0, 0x4), SIB(0, src1, src2));
	return;
	break;
    case DILL_L: case DILL_UL: case DILL_P:
	if (smi->stack_align == 4) {
	    type = DILL_I;
	}
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    BYTE_OUT3(s, opcode, ModRM(0x0, tmp_dest, 0x4), SIB(0,src1,src2));
    switch(type){
    case DILL_C:
	x86_lshi(s, dest, tmp_dest, 24);
	x86_rshi(s, dest, dest, 24);
	break;
    case DILL_S:
	x86_lshi(s, dest, tmp_dest, 16);
	x86_rshi(s, dest, dest, 16);
	break;
    case DILL_UC: case DILL_US:
	if (dest != tmp_dest)
	    x86_movi(s, dest, tmp_dest);
	break;
    }
}

extern void
x86_pbsloadi(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    int ltype = type;
    switch (type) {
    case DILL_D:
	ltype = DILL_L;
	dest = EAX;
	break;
    case DILL_F:
	ltype = DILL_I;
	dest = EAX;
	break;
    }
    x86_ploadi(s, ltype, junk, dest, src, offset);
    switch(type) {
    case DILL_F:
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word);
	x86_ploadi(s, DILL_F, 0, 0, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case DILL_D:
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word+4);
	x86_ploadi(s, ltype, junk, dest, src, offset + 4);
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word);
	x86_ploadi(s, DILL_D, 0, 0, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case DILL_L: case DILL_UL: case DILL_P: case DILL_I: case DILL_U:
	BYTE_OUT2(s, 0x0f, 0xc8 + dest);   /* byteswap dest */
	break;
    case DILL_S: case DILL_US:
	x86_bswap(s, 0, DILL_S, dest, dest);
/*	BYTE_OUT2(s, 0x86, (0xc0 | (dest << 3)) | (dest+4));*/ /* xchange dest */
	break;
    case DILL_C: case DILL_UC:
	break;
    }
}


extern void
x86_pbsload(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    int fdest = dest;
    int ltype = type;
    switch (type) {
    case DILL_D:
	ltype = DILL_L;
	fdest = dest;
	dest = EAX;
	break;
    case DILL_F:
	ltype = DILL_I;
	fdest = dest;
	dest = EAX;
	break;
    }
    x86_pload(s, ltype, junk, dest, src1, src2);
    switch(type) {
    case DILL_F:
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word);
	x86_ploadi(s, DILL_F, 0, 0, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case DILL_D:
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word+4);
	dill_addi(s, EAX, src1, src2);
	x86_ploadi(s, ltype, junk, dest, EAX, 4);
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	x86_pstorei(s, DILL_I, 0, EAX, _frame_reg, smi->conversion_word);
	x86_ploadi(s, DILL_D, 0, 0, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case DILL_L: case DILL_UL: case DILL_P: case DILL_I: case DILL_U:
	BYTE_OUT2(s, 0x0f, 0xc8 + dest);   /* byteswap dest */
	break;
    case DILL_S: case DILL_US:
	BYTE_OUT2(s, 0x86, (0xc0 | (dest << 3)) | (dest+4)); /* xchange dest */
	break;
    case DILL_C: case DILL_UC:
	break;
    }
}

static unsigned char st_opcodes[] = {
    0x88, /* DILL_C */
    0x88, /* DILL_UC */
    0x89, /* DILL_S */
    0x89, /* DILL_US */
    0x89, /* DILL_I */
    0x89, /* DILL_U */
    0x89, /* DILL_L */
    0x89, /* DILL_UL */
    0x89, /* DILL_P */
    0xd9, /* DILL_F */
    0xdd, /* DILL_D */
    0x00, /* DILL_V */
    0x00, /* DILL_B */
    0x89, /* DILL_EC */
};
extern void
x86_pstorei(dill_stream s, int type, int force_8087, int dest, int src, long offset)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((smi->generate_SSE && !force_8087) &&
	((type == DILL_F) || (type == DILL_D))) {
	int opcode = 0;
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    break;
	}
	if (src != ESP) {
	    if (((long)offset <= 127) && ((long)offset > -128)) {
	        BYTE_OUT5(s, opcode, 0x0f, 0x11, ModRM(0x1, dest, src), offset & 0xff);
	    } else {
	        BYTE_OUT4I(s, opcode, 0x0f, 0x11, ModRM(0x2, dest, src), offset);
	    }
	} else {
	    if (((long)offset <= 127) && ((long)offset > -128)) {
	        BYTE_OUT6(s, opcode, 0x0f, 0x11, ModRM(0x1, dest, 0x4), 0x24, offset & 0xff);
	    } else {
	        BYTE_OUT5I(s, opcode, 0x0f, 0x11, ModRM(0x2, dest, 0x4), 0x24, offset);
	    }
	}
	return;
    }
    switch(type) {
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, move to eax */
	    x86_movi(s, EAX, dest);
	    dest = EAX;
	}
	break;
    }
    switch (type) {
    case DILL_F:
    case DILL_D:
	dest = 3;
	break;
    case DILL_S: case DILL_US:
	BYTE_OUT1(s, 0x66);
	break;
    default:
	break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    if (src != ESP) {
        if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT3(s, st_opcodes[type], ModRM(0x1, dest, src), offset & 0xff);
        } else {
	    BYTE_OUT2I(s, st_opcodes[type], ModRM(0x2, dest, src), offset);
        }
    } else {
        /* can't use ModRM straight for SP */
        if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT4(s, st_opcodes[type], ModRM(0x1, dest, 0x4), 0x24, offset &0xff);
	} else {
	    BYTE_OUT3I(s, st_opcodes[type], ModRM(0x2, dest, 0x4), 0x24, offset &0xff);
	}
    }
}

extern void
x86_sse_pstorei(dill_stream s, int type, int junk, int dest, int src, long offset)
{
    unsigned char opcode = st_opcodes[type];
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((type == DILL_F) || (type == DILL_D)) {
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    break;
	}
	if (((long)offset <= 127) && ((long)offset > -128)) {
	    BYTE_OUT5(s, opcode, 0x0f, 0x11, ModRM(0x1, dest, src), offset & 0xff);
	} else {
	    BYTE_OUT4I(s, opcode, 0x0f, 0x11, ModRM(0x2, dest, src), offset);
	}
	return;
    }
    switch(type) {
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, move to eax */
	    x86_movi(s, EAX, dest);
	    dest = EAX;
	}
	break;
    }
    switch (type) {
    case DILL_F:
    case DILL_D:
	dest = 3;
	break;
    case DILL_S: case DILL_US:
	BYTE_OUT1(s, 0x66);
	break;
    default:
	break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    if (((long)offset <= 127) && ((long)offset > -128)) {
	BYTE_OUT3(s, st_opcodes[type], ModRM(0x1, dest, src), offset & 0xff);
    } else {
	BYTE_OUT2I(s, st_opcodes[type], ModRM(0x2, dest, src), offset);
    }
}

extern void
x86_pstore(dill_stream s, int type, int force_8087, int dest, int src1, int src2)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((smi->generate_SSE && !force_8087) &&
	((type == DILL_F) || (type == DILL_D))) {
	int opcode = 0;
	switch (type) {
	case DILL_F:
	    opcode = 0xf3;
	    break;
	case DILL_D:
	    opcode = 0xf2;
	    break;
	}
	BYTE_OUT5(s, opcode, 0x0f, 0x11, ModRM(0x0, dest, 0x4), SIB(0, src1, src2));
	return;
    }
    switch(type) {
    case DILL_C:
    case DILL_UC:
	if (dest >= ESP) {
	    /* movb doesn't work for big regs, move to eax */
	    x86_movi(s, EAX, dest);
	    dest = EAX;
	}
	break;
    }
    switch (type) {
    case DILL_F:
    case DILL_D:
	dest = 3;
	break;
    case DILL_S: case DILL_US:
	BYTE_OUT1(s, 0x66);
	break;
    default:
	break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
	smi->pending_prefix = 0;
    }
    BYTE_OUT3(s, st_opcodes[type], ModRM(0x0, dest, 0x4), SIB(0,src1,src2));
}

extern long dill_x86_hidden_mod(long a, long b);
extern long dill_x86_hidden_umod(unsigned long a, unsigned long b);

extern void x86_mod(dill_stream s, int data1, int data2, int dest, int src1, 
		      int src2)
{
    int return_reg;
    if (data1 == 1) {
	/* signed case */
	return_reg = dill_scalll(s, (void*)dill_x86_hidden_mod, "dill_x86_hidden_mod", "%l%l", src1, src2);
	dill_movl(s, dest, return_reg);
    } else {
	/* unsigned case */
	return_reg = dill_scalll(s, (void*)dill_x86_hidden_umod, "dill_x86_hidden_umod", "%l%l", src1, src2);
	dill_movul(s, dest, return_reg);
    }
}

extern void x86_modi(dill_stream s, int data1, int data2, int dest, int src1, 
		      long imm)
{
    x86_seti(s, _temp_reg, imm);
    x86_mod(s, data1, data2, dest, src1, _temp_reg);
}

extern void x86_div(dill_stream s, int op3, int op, int dest, int src1, 
		      int src2)
{
}

extern void x86_divi(dill_stream s, int op3, int op, int dest, int src, 
		      long imm)
{
}

extern void
x86_mov(dill_stream s, int type, int junk, int dest, int src)
{
    x86_pmov(s, type, dest, src);
}

static int b1_sse[] = {0x58, 0x5c, 0x59, 0x5e};
static int b1[] = {0xde, 0xde, 0xde, 0xde};
static int b2[] = {0xc1, 0xe9, 0xc9, 0xf9};

extern void x86_farith(s, op, typ, dest, src1, src2)
dill_stream s;
int op;
int typ;
int dest;
int src1;
int src2;
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    
    if (smi->generate_SSE) {
	int pre = 0xf3;
	if (src2 == dest) {
	    /* protect src2 from overwrite (below) */
	    x86_movd(s, XMM0, src2);
	    src2 = XMM0;
	}
	if (src1 != dest) x86_movd(s, dest, src1);
	if (typ == DILL_D) pre = 0xf2;
	BYTE_OUT4(s, pre, 0x0f, b1_sse[op], ModRM(0x3, dest, src2));
    } else {
	/* 8087 operations */
	BYTE_OUT2(s, b1[op], b2[op]);
    }
}

extern void x86_farith2(s, b1, typ, dest, src)
dill_stream s;
int b1;
int typ;
int dest;
int src;
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    if (smi->generate_SSE) {
	/* this is fneg */
	int op = 0xf3;
    
	/* clear dest */
	{
	    BYTE_OUT3(s, 0x0f, 0x57, ModRM(0x3, dest, dest));
	}
	if (typ == DILL_D) op = 0xf2;
	BYTE_OUT4(s, op, 0x0f, 0x5c, ModRM(0x3, dest, src));
    } else {
	BYTE_OUT2(s, 0xd9, 0xe0);
    }
}

extern void x86_arith3(s, op, commut, dest, src1, src2)
dill_stream s;
int op;
int commut;
int dest;
int src1;
int src2;
{
    if (commut && (dest == src1)) {
	BYTE_OUT2(s, op, ModRM(0x3, dest, src2));
    } else if (commut && (dest == src2)) {
	BYTE_OUT2(s, op, ModRM(0x3, dest, src1));
    } else if (dest == src2) {
	assert(op == 0x2b);	/* must be subtract */
	BYTE_OUT2(s, 0xf7, ModRM(0x3, 0x3, dest));   /* neg src2/dest */
	BYTE_OUT2(s, 0x03, ModRM(0x3, dest, src1));  /* add src1, dest */
    } else {
	BYTE_OUT2(s, MOV32, ModRM(0x3, src1, dest));
	BYTE_OUT2(s, op, ModRM(0x3, dest, src2));
    }
}

extern void x86_arith2(s, op, subop, dest, src)
dill_stream s;
int op;
int subop;
int dest;
int src;
{
    if (op == 0) {
	int tmp_dest = dest;
	/* must be not */
	if (dest >= ESP) {
	    tmp_dest = EAX;
	}
	BYTE_OUT3(s, 0x83, ModRM(0x3, 0x7, src), 0);  /* cmp */
	x86_seti(s, tmp_dest, 0);
	BYTE_OUT3(s, 0x0f, 0x94, ModRM(0x3, src, tmp_dest));  /* sete dest */
	if (tmp_dest != dest) {
	    x86_movi(s, dest, tmp_dest);
	}	    
    } else {
	if (src != dest) {
	    BYTE_OUT2(s, MOV32, ModRM(0x3, src, dest));
	}
	BYTE_OUT2(s, op, ModRM(0x3, subop, dest));
    }
}

extern void x86_bswap(s, junk, typ, dest, src)
dill_stream s;
int junk;
int typ;
int dest;
int src;
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    switch(typ) {
    case DILL_F: 
	BYTE_OUT4(s, 0x66, 0x0f, 0x7e, ModRM(0x3, src, EAX));
	x86_bswap(s, 0, DILL_I, EAX, EAX);
	BYTE_OUT4(s, 0x66, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
	break;
    case DILL_D:
	if (src != dest) x86_movd(s, dest, src);
	BYTE_OUT4(s, 0x0f, 0xc6, ModRM(0x3, dest, dest), 0xb1); /* SHUFPS */
	x86_pstorei(s, DILL_D, 0, dest, _frame_reg, smi->conversion_word);
	x86_pbsloadi(s, DILL_I, 0, EAX,  _frame_reg, smi->conversion_word);
	x86_pstorei(s, DILL_I, 0, EAX,  _frame_reg, smi->conversion_word);
	x86_pbsloadi(s, DILL_I, 0, EAX,  _frame_reg, smi->conversion_word + 4);
	x86_pstorei(s, DILL_I, 0, EAX,  _frame_reg, smi->conversion_word + 4);
	x86_ploadi(s, DILL_D, 0, dest, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case DILL_L: case DILL_UL: case DILL_P: 
    case DILL_I: case DILL_U:
	if (dest != src) {
	    x86_movi(s, dest, src);
	}
	BYTE_OUT2(s, 0x0f, 0xc8 + (dest&0x7));   /* byteswap dest */
	break;
    case DILL_S: case DILL_US:
	if (dest != src) {
	    x86_movi(s, dest, src);
	}
	/* byteswap 32 bits and shift down 16 */
	BYTE_OUT2(s, 0x0f, 0xc8 + (dest&0x7));   /* byteswap dest */
	x86_rshi(s, dest, dest, 16);
	break;
    case DILL_C: case DILL_UC:
	break;
    }
}

extern void x86_mul(s, sign, imm, dest, src1, src2)
dill_stream s;
int sign;
int imm;
int dest;
int src1;
int src2;
{
    /* make src1 be EAX */
    if (dest != EAX) {
	x86_push_reg(s, EAX);
    }
    if (dest != EDX) {
	x86_push_reg(s, EDX);
    }
	
    if ((src2 == EAX) && !imm){
	int tmp = src2;
	src1 = src2;
	src2 = tmp;
    }
    if (src1 != EAX) {
	x86_movi(s, EAX, src1);
    }
    if (imm == 0) {
	BYTE_OUT2(s, 0xf7, ModRM(0x3, sign ? 0x5 : 0x4, src2));
    } else {
	/* src2 is really immediate */
	if (sign) {
	    BYTE_OUT2I(s, 0x69, ModRM(0x3, 0, EAX), src2);
	} else {
	    x86_seti(s, EDX, src2);
	    BYTE_OUT2(s, 0xf7, ModRM(0x3, 0x4, EDX));
	}
    }
    if (dest != EDX) {
	x86_pop_reg(s, EDX);
    }
    if (dest != EAX) {
	x86_movi(s, dest, EAX);
	x86_pop_reg(s, EAX);
    }
}

extern void x86_div_modi(s, sign, div, dest, src1, imm)
dill_stream s;
int sign;
int div;
int dest;
int src1;
long imm;
{
    x86_push_reg(s, EBP);
    x86_seti(s, EBP, imm);
    x86_div_mod(s, sign, div, dest, src1, EBP);
    x86_pop_reg(s, EBP);
}

extern void x86_div_mod(s, sign, div, dest, src1, src2)
dill_stream s;
int sign;
int div;
int dest;
int src1;
int src2;
{
    int tmp_src2 = src2;

    /* make src1 be EAX */
    if (dest != EAX) {
	x86_push_reg(s, EAX);
    }
    if (dest != EDX) {
	x86_push_reg(s, EDX);
    }
	
    if (src1 != EAX) {
	x86_movi(s, EAX, src1);
    }
    if (src2 == EDX) {
	tmp_src2 = EBP;
	x86_push_reg(s, EBP);
	x86_movi(s, EBP, src2);
    }
    if (sign) {
	x86_rshai(s, EDX, EAX, 31);
    } else {
	x86_seti(s, EDX, 0);
    }
    BYTE_OUT2(s, 0xf7, ModRM(0x3, sign ? 0x7 : 0x6, tmp_src2));
    if (src2 == EDX) {
	x86_pop_reg(s, EBP);
    }
    if (div && (dest != EAX)) {
	x86_movi(s, dest, EAX);
    }
    if (!div && (dest != EDX)) {
	x86_movi(s, dest, EDX);
    }
    if (dest != EDX) {
	x86_pop_reg(s, EDX);
    }
    if (dest != EAX) {
	x86_pop_reg(s, EAX);
    }
	
}

static int group1_eax_op[] = {
    0x05 /* add */,
    0x0d /* or */,
    0x15 /* adc */,
    0x1d /* sbb */,
    0x25 /* and */,
    0x2d /* sub */,
    0x35 /* xor */,
    0x3d /* cmp */
};

extern void x86_arith3i(s, op, junk, dest, src, imm)
dill_stream s;
int op;
int junk;
int dest;
int src;
long imm;
{
    if (dest != src) {
	BYTE_OUT2(s, MOV32, ModRM(0x3, src, dest));
    }
    if ((imm <= 127) && (imm > -128)) {
	BYTE_OUT3(s, 0x83, ModRM(0x3, op, dest), imm & 0xff);
	return;
    }
    if (dest == EAX) {
	BYTE_OUT1I(s, group1_eax_op[op], imm);
    } else {
	BYTE_OUT2I(s, 0x81, ModRM(0x3, op, dest), imm);
    }
}

extern void x86_shift(s, op, junk, dest, src1, src2)
dill_stream s;
int op;
int junk;
int dest;
int src1;
int src2;
{
    int tmp_dest = dest;
    if ((dest == ECX) || (dest == src2)) {
	x86_push_reg(s, EAX);
	tmp_dest = EAX;
    }
    if (tmp_dest != src1) {
	x86_movi(s, tmp_dest, src1);
    }
    if (src2 != ECX) {
	x86_push_reg(s, ECX);
	x86_movi(s, ECX, src2);
    }
    BYTE_OUT2(s, 0xd3, ModRM(0x3, op, tmp_dest));
    if (src2 != ECX) {
	x86_pop_reg(s, ECX);
    }
    if ((dest == ECX) || (dest == src2)) {
	x86_movi(s, dest, tmp_dest);
	x86_pop_reg(s, EAX);
    }
}

extern void x86_shifti(s, op, junk, dest, src, imm)
dill_stream s;
int op;
int junk;
int dest;
int src;
long imm;
{
    if (dest != src) {
	x86_movi(s, dest, src);
    }
    BYTE_OUT3(s, 0xc1, ModRM(0x3, op, dest), imm & 0xff);
}

#define CONV(x,y) ((x*100)+y)
extern void
x86_convert(dill_stream s, int from_type, int to_type, 
	      int dest, int src)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    switch(CONV(from_type, to_type)) {
    case CONV(DILL_I, DILL_U):
    case CONV(DILL_I,DILL_UL):
    case CONV(DILL_I, DILL_L):
    case CONV(DILL_U,DILL_I):
    case CONV(DILL_U,DILL_UL):
    case CONV(DILL_U,DILL_L):
    case CONV(DILL_UL,DILL_I):
    case CONV(DILL_UL,DILL_U):
    case CONV(DILL_UL,DILL_L):
    case CONV(DILL_L,DILL_U):
    case CONV(DILL_L,DILL_UL):
    case CONV(DILL_L,DILL_I):
    case CONV(DILL_P,DILL_UL):
    case CONV(DILL_UL,DILL_P):
	if(src == dest) return;
	x86_movi(s, dest,src);
	break;
    case CONV(DILL_D,DILL_F):
    case CONV(DILL_F,DILL_D):
	if (smi->generate_SSE) {
	    BYTE_OUT4(s, (from_type == DILL_D) ? 0xf2 : 0xf3, 0xf, 0x5a, ModRM(0x3, dest, src));
	}
	break;
    case CONV(DILL_F,DILL_I):
    case CONV(DILL_F,DILL_L):
    case CONV(DILL_F,DILL_S):
    case CONV(DILL_F,DILL_C):
    case CONV(DILL_D,DILL_I):
    case CONV(DILL_D,DILL_S):
    case CONV(DILL_D,DILL_C):
    case CONV(DILL_D,DILL_L):
	if (smi->generate_SSE) {
	    BYTE_OUT4(s, (from_type == DILL_D) ? 0xf2 : 0xf3, 0xf, 0x2c, ModRM(0x3, dest, src));
	    break;
	}
	/* fall through */
    case CONV(DILL_F,DILL_U):
    case CONV(DILL_F,DILL_UL):
    case CONV(DILL_F,DILL_US):
    case CONV(DILL_F,DILL_UC):
    case CONV(DILL_D,DILL_U):
    case CONV(DILL_D,DILL_UL):
    case CONV(DILL_D,DILL_US):
    case CONV(DILL_D,DILL_UC):
	if (smi->generate_SSE) {
	    /* NO direct SSE 32-bit unsigned conversion */
	    /* load operand from mmx to 8087 stack */
	    x86_pstorei(s, from_type, 0, src,  _frame_reg, smi->conversion_word);
	    x86_ploadi(s, from_type, 1 /* force 8087 */, 0, _frame_reg, smi->conversion_word);
	}
	/* use dest for float control word */
	/* fstcw (store control word) */
	BYTE_OUT3(s, 0xd9, ModRM(0x1, 0x7, _frame_reg), smi->fcu_word);
	x86_push_reg(s, EAX);
	/* movsw eax, fcu(ebx)   Load control word to reg */
	BYTE_OUT3I(s, 0x66, 0x8b, ModRM(0x2, EAX, _frame_reg), smi->fcu_word);
	/* modify control word OR it with 0x60 to set RC to 11 */
	BYTE_OUT1(s, 0x66);
	BYTE_OUT3(s, 0x0d,0x0, 0xc);
	/* movsw fcu(ebx), eax   store mod control word */
	BYTE_OUT3I(s, 0x66, 0x89, ModRM(0x2, EAX, _frame_reg), smi->fcu_word +2);
	/* fldcw (load modified control word) */
	BYTE_OUT3(s, 0xd9, ModRM(0x1, 0x5, _frame_reg), smi->fcu_word + 2);
	x86_pop_reg(s, EAX);
	/* fistpl (do the conversion) */
	switch(to_type){
	case DILL_U: case DILL_UL:
	case DILL_UC: case DILL_US:
	    BYTE_OUT3(s, 0xdf, ModRM(0x1, 0x7, _frame_reg), smi->conversion_word);
	    break;
	case DILL_I: case DILL_L:
	case DILL_C: case DILL_S:
	    BYTE_OUT3(s, 0xdb, ModRM(0x1, 0x3, _frame_reg), smi->conversion_word);
	    break;
	}
	/* fldcw (restore original) */
	BYTE_OUT3(s, 0xd9, ModRM(0x1, 0x5, _frame_reg), smi->fcu_word);
	x86_ploadi(s, DILL_I, 0, dest, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	break;
    case CONV(DILL_I,DILL_D):
    case CONV(DILL_L,DILL_D):
    case CONV(DILL_I,DILL_F):
    case CONV(DILL_L,DILL_F):
    case CONV(DILL_C,DILL_D):
    case CONV(DILL_C,DILL_F):
    case CONV(DILL_S,DILL_D):
    case CONV(DILL_S,DILL_F):
	x86_pstorei(s, DILL_I, 0, src, _frame_reg, smi->conversion_word);
	s->p->used_frame++;
	if (smi->generate_SSE) {
	    	BYTE_OUT4(s, (to_type == DILL_D) ? 0xf2 : 0xf3, 0xf, 0x2a, ModRM(0x3, dest, src));
	} else {
	    BYTE_OUT3(s, 0xdb, ModRM(0x1, 0x0, _frame_reg), smi->conversion_word);
	}
	break;
    case CONV(DILL_U,DILL_D):
    case CONV(DILL_UL,DILL_D):
    case CONV(DILL_UC,DILL_D):
    case CONV(DILL_US,DILL_D):
    case CONV(DILL_UC,DILL_F):
    case CONV(DILL_US,DILL_F):
    case CONV(DILL_U,DILL_F):
    case CONV(DILL_UL,DILL_F):
	if (smi->generate_SSE) {
	    x86_pstorei(s, DILL_I, 0, src, _frame_reg, smi->conversion_word);
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word + 4, 0);
	    BYTE_OUT3(s, 0xdf, ModRM(0x1, 0x5, _frame_reg), smi->conversion_word);
	    x86_pstorei(s, to_type, 1 /* force 8087 */, 0,  _frame_reg, smi->conversion_word);
	    x86_ploadi(s, to_type, 0, dest, _frame_reg, smi->conversion_word);
	    s->p->used_frame++;
	} else {
	    x86_pstorei(s, DILL_I, 0, src, _frame_reg, smi->conversion_word);
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word + 4, 0);
	    BYTE_OUT3(s, 0xdf, ModRM(0x1, 0x5, _frame_reg), smi->conversion_word);
	    s->p->used_frame++;
	}
	break;
    case CONV(DILL_C,DILL_I):
    case CONV(DILL_C,DILL_L):
    case CONV(DILL_C,DILL_U):
    case CONV(DILL_C,DILL_UL):
    case CONV(DILL_C, DILL_S):
    case CONV(DILL_S, DILL_C):
    case CONV(DILL_US, DILL_C):
	/* signext24 - lsh24, rsha24 */
	x86_lshi(s, dest, src, 24);
	x86_rshai(s, dest, dest, 24);
	break;
    case CONV(DILL_I, DILL_UC):
    case CONV(DILL_S, DILL_UC):
    case CONV(DILL_L, DILL_UC):
    case CONV(DILL_I, DILL_C):
    case CONV(DILL_U, DILL_C):
    case CONV(DILL_L, DILL_C):
    case CONV(DILL_UL, DILL_C):
    case CONV(DILL_C, DILL_UC):
    case CONV(DILL_U, DILL_UC):
    case CONV(DILL_US, DILL_UC):
    case CONV(DILL_UL, DILL_UC):
	x86_andi(s, dest, src, 0xff);
	break;
    case CONV(DILL_S,DILL_I):
    case CONV(DILL_S,DILL_L):
    case CONV(DILL_S,DILL_U):
    case CONV(DILL_S,DILL_UL):
    case CONV(DILL_S,DILL_US):
    case CONV(DILL_US,DILL_S):
	/* signext16 - lsh16, rsha16 */
	x86_lshi(s, dest, src, 16);
	x86_rshai(s, dest, dest, 16);
	break;
    case CONV(DILL_C, DILL_US):
	/* signext24 - lsh24, rsha24, trunc 16 */
	x86_lshi(s, dest, src, 24);
	x86_rshai(s, dest, dest, 24);
	x86_andi(s, dest, dest, 0xffff);
	break;
    case CONV(DILL_US,DILL_I):
    case CONV(DILL_US,DILL_L):
    case CONV(DILL_US,DILL_U):
    case CONV(DILL_US,DILL_UL):
    case CONV(DILL_I, DILL_S):
    case CONV(DILL_U, DILL_S):
    case CONV(DILL_L, DILL_S):
    case CONV(DILL_UL, DILL_S):
    case CONV(DILL_I, DILL_US):
    case CONV(DILL_U, DILL_US):
    case CONV(DILL_L, DILL_US):
    case CONV(DILL_UL, DILL_US):
	/* zero uppper 16 - lsh16, rsh16 */
	x86_lshi(s, dest, src, 16);
	x86_rshi(s, dest, dest, 16);
	break;
    default:
	printf("Unknown case in x86 convert %d\n", CONV(from_type,to_type));
    }
}

static unsigned char set_op_conds[] = {
    0x94, /* dill_beq_code */  /* signed */
    0x9d, /* dill_bge_code */
    0x9F, /* dill_bgt_code */
    0x9e, /* dill_ble_code */
    0x9c, /* dill_blt_code */
    0x95, /* dill_bne_code */

    0x94, /* dill_beq_code */  /* unsigned */
    0x93, /* dill_bge_code */
    0x97, /* dill_bgt_code */ 
    0x96, /* dill_ble_code */
    0x92, /* dill_blt_code */
    0x95, /* dill_bne_code */
};

static unsigned char fop_conds[] = {
    0x94, /* dill_beq_code */   /* z = 1*/
    0x96, /* dill_bge_code */   /* jna */
    0x92, /* dill_bgt_code */
    0x93, /* dill_ble_code */   /* c = 0 */
    0x97, /* dill_blt_code */   
    0x95, /* dill_bne_code */
};

extern void
x86_compare(dill_stream s, int op, int type, int dest, int src1, int src2)
{
    int op_cond;
    switch(type) {
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    }

    if (type == DILL_D) {
	BYTE_OUT4(s, 0x66, 0x0f, 0x2e, ModRM(0x3, src2, src1));
	op_cond = fop_conds[op];
    } else if (type == DILL_F) {
	BYTE_OUT3(s, 0x0f, 0x2e, ModRM(0x3, src2, src1));
	op_cond = fop_conds[op];
    } else {
	BYTE_OUT2(s, 0x39, ModRM(0x3, src2, src1));
	op_cond = set_op_conds[op];
    }
    BYTE_OUT3(s, 0x0f, op_cond, ModRM(0x3, EAX,EAX));	  /* sete */
    BYTE_OUT3(s, 0x0f, 0xb6, ModRM(0x3, dest, EAX));	/* movzbl */
}

extern void
x86_comparei(dill_stream s, int op, int type, int dest, int src, long imm)
{
    switch(type) {
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    }
    BYTE_OUT2I(s, 0x81, ModRM(0x3, 0x7, src), imm);  /* cmp */
    BYTE_OUT3(s, 0x0f, set_op_conds[op],ModRM(0x3, EAX,EAX));	  /* sete */
    BYTE_OUT3(s, 0x0f, 0xb6, ModRM(0x3, dest, EAX));	/* movzbl */
}

static unsigned char br_op_conds[] = {
    0x84, /* dill_beq_code */  /* signed */
    0x8d, /* dill_bge_code */
    0x8F, /* dill_bgt_code */
    0x8e, /* dill_ble_code */
    0x8c, /* dill_blt_code */
    0x85, /* dill_bne_code */

    0x84, /* dill_beq_code */  /* unsigned */
    0x83, /* dill_bge_code */
    0x87, /* dill_bgt_code */ 
    0x86, /* dill_ble_code */
    0x82, /* dill_blt_code */
    0x85, /* dill_bne_code */
};

static unsigned char br_fop_conds[] = {
    0x84, /* dill_beq_code */   /* z = 1*/
    0x86, /* dill_bge_code */
    0x82, /* dill_bgt_code */
    0x83, /* dill_ble_code */   /* c = 0 */
    0x87, /* dill_blt_code */   
    0x85, /* dill_bne_code */
};

extern void
x86_branch(dill_stream s, int op, int type, int src1, int src2, int label)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    switch(type) {
    case DILL_F:
    case DILL_D:
	if (smi->generate_SSE) {
	    if (type == DILL_D) {
		BYTE_OUT4(s, 0x66, 0x0f, 0x2e, ModRM(0x3, src2, src1));
	    } else if (type == DILL_F) {
		BYTE_OUT3(s, 0x0f, 0x2e, ModRM(0x3, src2, src1));
	    }
	} else {
	    BYTE_OUT2(s, 0xde, 0xd9);   /* fcompp */
	    BYTE_OUT2(s, 0xdf, 0xe0);   /* fnstsw ax */
	    BYTE_OUT1(s, 0x9e);	    /* sahf */
	}
	dill_mark_branch_location(s, label);
	BYTE_OUT2I(s, 0x0f, br_fop_conds[op], 0);
	break;
    case DILL_US:
    case DILL_UC:
    case DILL_U:
    case DILL_UL:
	op += 6; /* second set of codes */
	/* fall through */
    default:
	BYTE_OUT2(s, 0x39, ModRM(0x3, src2, src1));
	dill_mark_branch_location(s, label);
	BYTE_OUT2I(s, 0x0f, br_op_conds[op], 0);
    }
    x86_nop(c);
}

extern void 
x86_jump_to_label(dill_stream s, unsigned long label)
{
    dill_mark_branch_location(s, label);
    BYTE_OUT1I(s, 0xe9, 0);
}

extern void x86_jump_to_reg(dill_stream s, unsigned long reg)
{
    BYTE_OUT2(s, 0xff, ModRM(0x3, 0x4, reg));
}

extern void x86_jump_to_imm(dill_stream s, void *imm)
{
    x86_setp(s, DILL_P, 0, EAX, imm);
    BYTE_OUT2(s, 0xff, ModRM(0x3, 0x4, EAX));
}

extern void 
x86_jal(dill_stream s, int return_addr_reg, int target)
{
/* jump, source addr to return_addr_reg */
}

extern void 
x86_special(dill_stream s, special_operations type, long param)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    switch (type) {
    case DILL_NOP:
	break;
    case DILL_SEGMENTED_FOLLOWS:
	switch(param) {
	case DILL_X86_CS_PREFIX:
	    smi->pending_prefix = 0x2e;
	    break;
	case DILL_X86_SS_PREFIX:
	    smi->pending_prefix = 0x36;
	    break;
	case DILL_X86_DS_PREFIX:
	    smi->pending_prefix = 0x3e;
	    break;
	case DILL_X86_ES_PREFIX:
	    smi->pending_prefix = 0x26;
	    break;
	case DILL_X86_FS_PREFIX:
	    smi->pending_prefix = 0x64;
	    break;
	case DILL_X86_GS_PREFIX:
	    smi->pending_prefix = 0x65;
	    break;
	default:
	    fprintf(stderr, "Unknown x86 segment prefix!\n");
	}
	break;
    }
}

static void internal_push(dill_stream s, int type, int immediate, 
			  void *value_ptr)
{
    if (!immediate) {
	x86_push_reg(s, *(int*)value_ptr);
    } else {
	BYTE_OUT1I(s, 0x68, *(int*)value_ptr);
    }
}

static void push_init(dill_stream s)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    smi->cur_arg_offset = 0;
    smi->call_stack_adjust_point = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    /* reserve space.  This sub will be overwritten */
    dill_subii(s, ESP, ESP, 60);
}

extern void x86_push(dill_stream s, int type, int reg)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if ((type == DILL_V) && (reg <= -1)) {
	push_init(s);
    } else if ((type == DILL_F) || (type == DILL_D)) {
        x86_pstorei(s, DILL_D, 0, reg, ESP, smi->cur_arg_offset);
	smi->cur_arg_offset += 8;
    } else {
	x86_pstorei(s, DILL_I, 0, reg, ESP, smi->cur_arg_offset);
	smi->cur_arg_offset += 4;
    }
}

extern void x86_pushi(dill_stream s, int type, long value)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    x86_seti(s, EAX, value);
    x86_pstorei(s, DILL_I, 0, EAX, ESP, smi->cur_arg_offset);
    smi->cur_arg_offset += 4;
}

extern void x86_pushfi(dill_stream s, int type, double value)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    internal_push(s, type, 1, &value);
    smi->cur_arg_offset += 8;
}

extern void x86_pushpi(dill_stream s, int type, void *value)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    x86_seti(s, EAX, value);
    x86_pstorei(s, DILL_I, 0, EAX, ESP, smi->cur_arg_offset);
    smi->cur_arg_offset += 4;
}

extern int x86_calli(dill_stream s, int type, void *xfer_address, const char *name)
{
    int caller_side_ret_reg = EAX;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    /* save temporary registers */
    dill_mark_call_location(s, name, xfer_address);
    BYTE_OUT1I(s, 0xe8, 0);
    /* restore temporary registers */
    if ((type == DILL_D) || (type == DILL_F)) {
        if (smi->generate_SSE) {
	    x86_pstorei(s, type, 1 /* force 8087 */, 0,  _frame_reg, smi->conversion_word);
	    x86_ploadi(s, type, 0, XMM0, _frame_reg, smi->conversion_word);
	    caller_side_ret_reg = XMM0;
	    s->p->used_frame++;
	}
    }
    /* make stack adjust be multiple of 16 */
    smi->cur_arg_offset = ((smi->cur_arg_offset + 15) & -16);
    int tmp = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    
    s->p->cur_ip = (char*)s->p->code_base + smi->call_stack_adjust_point;
    /* backpatch in stack adjust at front end of arg load */
    dill_subii(s, ESP, ESP, smi->cur_arg_offset);
    /* restore IP */
    s->p->cur_ip = (char*)s->p->code_base + tmp;
    
    dill_addii(s, ESP, ESP, smi->cur_arg_offset);
    return caller_side_ret_reg;
}

extern int x86_callr(dill_stream s, int type, int src)
{
    int caller_side_ret_reg = EAX;
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;

    /* save temporary registers */
    /* call through reg */
    BYTE_OUT2(s, 0xff, ModRM(0x3, 0x2, src));
    /* restore temporary registers */
    if ((type == DILL_D) || (type == DILL_F)) {
/*	caller_side_ret_reg = _f0;*/
    }
    /* make stack adjust be multiple of 16 */
    smi->cur_arg_offset = ((smi->cur_arg_offset + 15) & -16);
    int tmp = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    
    s->p->cur_ip = (char*)s->p->code_base + smi->call_stack_adjust_point;
    /* backpatch in stack adjust at front end of arg load */
    dill_subii(s, ESP, ESP, smi->cur_arg_offset);
    /* restore IP */
    s->p->cur_ip = (char*)s->p->code_base + tmp;
    
    dill_addii(s, ESP, ESP, smi->cur_arg_offset);
    return caller_side_ret_reg;
}

extern void
x86_branchi(dill_stream s, int op, int type, int src, long imm, int label)
{
    switch(type) {
    case DILL_F:
    case DILL_D:
	fprintf(stderr, "Shouldn't happen\n");
	break;
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
/*
	switch(op) {
	case dill_bge_code: {
	    imm = imm-1;
	    op = dill_bgt_code;
	    break;
	}
	case dill_blt_code: {
	    imm = imm-1;
	    op = dill_ble_code;
	    break;
	}
	}
*/
	op += 6; /* second set of codes */
	/* fall through */
    default:
/*	BYTE_OUT2(s, 0x39, ModRM(0x3, src2, src1));*/
	BYTE_OUT2I(s, 0x81, ModRM(0x3, 0x7, src), imm);  /* cmp */
	dill_mark_branch_location(s, label);
	BYTE_OUT2I(s, 0x0f, br_op_conds[op], 0);
    }
}

extern void x86_ret(dill_stream s, int data1, int data2, int src)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
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
	if (src != EAX) x86_movi(s, EAX, src);
	break;
    case DILL_F:
	if (smi->generate_SSE) {
	    /* store XMM reg to conversion word, reload to ST(0) */
	    x86_pstorei(s, DILL_F, 0, src, _frame_reg, smi->conversion_word);
	    x86_ploadi(s, DILL_F, 1 /* force 8087*/, src, 
		       _frame_reg, smi->conversion_word);
	    s->p->used_frame++;
	}
	break;
    case DILL_D:
	if (smi->generate_SSE) {
	    /* store XMM reg to conversion word, reload to ST(0) */
	    x86_pstorei(s, DILL_D, 0, src, _frame_reg, smi->conversion_word);
	    x86_ploadi(s, DILL_D, 1 /* force 8087*/, src, 
		       _frame_reg, smi->conversion_word);
	    s->p->used_frame++;
	}
	break;
    }
    x86_proc_ret(s);
}

extern void x86_retf(dill_stream s, int data1, int data2, double imm)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    if (smi->generate_SSE) {
	switch(data1) {
	case DILL_F:
	    if (smi->generate_SSE) {
		x86_setf(s, DILL_F, 0, XMM0, imm);
		x86_pstorei(s, DILL_F, 0, XMM0, _frame_reg, smi->conversion_word);
		x86_ploadi(s, DILL_F, 1 /* force 8087*/, 0, 
			   _frame_reg, smi->conversion_word);
		s->p->used_frame++;
	    }
	    break;
	case DILL_D:
	    if (smi->generate_SSE) {
		x86_setf(s, DILL_D, 0, XMM0, imm);
		x86_pstorei(s, DILL_D, 0, XMM0, _frame_reg, smi->conversion_word);
		x86_ploadi(s, DILL_D, 1 /* force 8087*/, 0, 
			   _frame_reg, smi->conversion_word);
		s->p->used_frame++;
	    }
	    break;
	}
    }
    x86_proc_ret(s);
}

extern void x86_reti(dill_stream s, int data1, int data2, long imm)
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
	x86_seti(s, EAX, imm);
	break;
/*    case DILL_F:
	x86_setf(s, _f0, imm);
    case DILL_D:
    x86_setd(s, _f0, imm);*/
    }
    x86_proc_ret(s);
}

extern void x86_rt_call_link(char *code, call_t *t);

static void
x86_call_link(dill_stream s)
{
    x86_rt_call_link(s->p->code_base, &s->p->call_table);
}

static void
x86_data_link(dill_stream s)
{
/*    struct branch_table *t = &s->p->branch_table;
    int i;
    for (i=0; i < t->data_mark_count; i++) {
	int label = t->data_marks[i].label;
	void *label_addr = t->label_locs[label] + (char*)s->p->code_base;
	*t->data_marks[i].addr = label_addr;
	}*/
}

static void
x86_branch_link(dill_stream s)
{
    struct branch_table *t = &s->p->branch_table;
    int i;

    for(i=0; i< t->branch_count; i++) {
	int label = t->branch_locs[i].label;
	int label_offset = t->label_locs[label] - t->branch_locs[i].loc;
	char *branch_addr = (char*)((char *)s->p->code_base + 
				  t->branch_locs[i].loc);
	int offset;
	if (*branch_addr == 0x0f) {
	    branch_addr +=2; /* conditional */
	    offset = 6;
	} else {
	    branch_addr++; /* unconditional */
	    offset = 5;
	}
	*(int*)branch_addr = label_offset - offset;
    }
}

static void
x86_emit_save(dill_stream s)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    void *save_ip = s->p->cur_ip;
    int ar_size = smi->act_rec_size;
    int prefix_size;
    ar_size = roundup(ar_size, 8) + 16;

    s->p->cur_ip = (char*)s->p->code_base;

    prefix_size = generate_prefix_code(s, 0, ar_size);

    s->p->cur_ip = (char*)s->p->code_base - prefix_size + smi->backpatch_offset;
    s->p->fp = (char*)s->p->cur_ip;

    /* re-generate in the right place */
    prefix_size = generate_prefix_code(s, 0, ar_size);

    s->p->cur_ip = save_ip;
}
    
extern void
x86_end(s)
dill_stream s;
{
    x86_proc_ret(s);
    x86_branch_link(s);
    x86_call_link(s);
    x86_data_link(s);
    x86_emit_save(s);
}

extern void
x86_package_end(dill_stream s)
{
    x86_proc_ret(s);
    x86_branch_link(s);
    x86_emit_save(s);
    /* at this point, code segment is finalized */
    
}

extern void *
x86_clone_code(s, new_base, available_size)
dill_stream s;
void *new_base;
int available_size;
{
    int size = dill_code_size(s);
    void *old_base = s->p->code_base;
    void *native_base = s->p->code_base;
    if (available_size < size) {
	return NULL;
    }
    if (native_base == NULL) native_base = s->p->native.code_base;
    memcpy(new_base, native_base, size);
    s->p->code_base = new_base;
    s->p->cur_ip = (char*)new_base + size;
    s->p->fp = new_base;
    x86_branch_link(s);
    x86_call_link(s);
    x86_data_link(s);
    s->p->code_base = old_base;
    s->p->cur_ip = (char*)old_base + size;
    s->p->fp = old_base;
    return new_base;
}

extern void
x86_pset(dill_stream s, int type, int junk, int dest, long imm)
{
    x86_seti(s, dest, imm);
}	

extern void
x86_setp(dill_stream s, int type, int junk, int r, void *val) 
{
    BYTE_OUT1PI(s, 0xb8 + r, val);
}

extern void
x86_setf(dill_stream s, int type, int junk, int dest, double imm)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    union {
	float f;
	int i;
    } a;
    union {
	double d;
	int i[2];
    } b;

    if (smi->generate_SSE) {
	if (type == DILL_F) {
	    a.f = imm;
	    x86_seti(s, EAX, a.i);
	    BYTE_OUT4(s, 0x66, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
	} else {
	    b.d = imm;
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word, b.i[0]);
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word + 4, b.i[1]);

	    x86_ploadi(s, DILL_D, 0, dest, 
		       _frame_reg, smi->conversion_word);
	    s->p->used_frame++;
	}
    } else {
	if (type == DILL_F) {
	    a.f = (float) imm;
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word, a.i);
	    /* flds */
	    BYTE_OUT3(s, 0xd9, ModRM(0x1, 0x0, _frame_reg), smi->conversion_word);
	    s->p->used_frame++;
	} else {
	    b.d = imm;
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word, b.i[0]);
	    BYTE_OUT3I(s, 0xc7, ModRM(0x1, 0x0, _frame_reg), 
		       smi->conversion_word + 4, b.i[1]);
	    /* fldd */
	    BYTE_OUT3(s, 0xdd, ModRM(0x1, 0x0, _frame_reg), smi->conversion_word);
	    s->p->used_frame++;
	}
    }

}	

#define bit_R(x) (1<<x)

extern void
x86_reg_init(dill_stream s, x86_mach_info smi)
{
    s->p->var_i.init_avail[0] = (bit_R(EBX)|bit_R(ESI)|bit_R(EDI));
    s->p->var_i.members[0] = s->p->var_i.init_avail[0];
    s->p->tmp_i.init_avail[0] = (bit_R(EDX)|bit_R(ECX));
    s->p->tmp_i.members[0] = s->p->tmp_i.init_avail[0] | bit_R(EAX);
    if (smi->generate_SSE) {
	s->p->tmp_f.init_avail[0] = (bit_R(XMM1)|bit_R(XMM2)|bit_R(XMM3)|bit_R(XMM4)|bit_R(XMM5)|bit_R(XMM6)|bit_R(XMM7));
	s->p->var_f.init_avail[0] = 0;
    } else {
	s->p->tmp_f.init_avail[0] = 0;
	s->p->var_f.init_avail[0] = 0;
    }
    s->p->var_f.members[0] = (bit_R(XMM0)|s->p->var_f.init_avail[0]);
    s->p->tmp_f.members[0] = s->p->tmp_f.init_avail[0];
}

#define bit_mmx (1<<23)
#define bit_sse (1<<25)

extern void*
gen_x86_mach_info(s)
dill_stream s;
{
    static int host_supports_SSE = -1;
    x86_mach_info smi = malloc(sizeof(*smi));
    if (s->p->mach_info != NULL) {
		free(s->p->mach_info);
		s->p->mach_info = NULL;
		s->p->native.mach_info = NULL;
    }
    if (host_supports_SSE == -1) {
#if defined(HOST_X86) && defined(_MSC_VER)
#define _SSE2_FEATURE_BIT 0x04000000
	int dwFeature = 0;
	_asm {
	    push ebx
		push ecx
		push edx
		
		//get the Feature bits
		xor eax, eax
	        mov eax, 1
		cpuid
		mov dwFeature, edx
		
		pop edx
		pop ecx
		pop ebx
		}
	if (dwFeature & _SSE2_FEATURE_BIT) {
	    host_supports_SSE = 1;
	}
#elif defined(HOST_X86) && defined(__GNUC__)
	int fl1, fl2;
	/* Invoke CPUID(1), return %edx; caller can examine bits to
	   determine what's supported.  */
	__asm__ ("pushl %%ecx; pushl %%ebx; cpuid; popl %%ebx; popl %%ecx"
		 : "=d" (fl2), "=a" (fl1) : "1" (1) : "cc");
	host_supports_SSE = ((fl2 & bit_sse) == bit_sse);
	host_supports_SSE |= ((fl2 & bit_mmx) == bit_mmx);
#else
	host_supports_SSE = 1;
#endif
	if (getenv("DILL_NO_SSE")) host_supports_SSE = 0;
	if (s->dill_debug) {
	    if (host_supports_SSE) {
		printf("x86 SSE code generated\n");
	    } else {
		printf("x86 8087 code generated\n");
	    }
	}
    }
    smi->generate_SSE = host_supports_SSE;
    x86_reg_init(s, smi);
    smi->act_rec_size = 0;
    smi->stack_align = 4; /* 8 for x86v9 */
    smi->stack_constant_offset = 0; /* 2047 for x86v9 */
    smi->conversion_word = 0;
    smi->fcu_word = 0;
    smi->save_base = 0;
    smi->backpatch_offset = 0;
    return smi;
}
#if defined(HAVE_DIS_ASM_H) && !defined(NO_DISASSEMBLER)
/* GENERIC BINUTILS DISASSEMBLER */
#include "dis-asm.h"

#define MAXLENGTH (1<<23) /* Max length of function that can be disassembled */

#ifdef LINUX_KERNEL_MODULE
extern int
kfprintf(FILE *file, const char *fmt, ...)
{
    static char printk_buf[1024];
    int val;
    va_list ap;
    va_start(ap, fmt);
    val = vsnprintf(printk_buf, sizeof(printk_buf), fmt, ap);
    printk("%s", printk_buf);
    va_end(ap);
    return val;
}

#undef  stdout
#define stdout (FILE*)0
#define FPRINTF_FUNCTION kfprintf
#else
#define FPRINTF_FUNCTION fprintf
#endif
extern int
x86_init_disassembly_info(dill_stream s, void * ptr)
{
    struct disassemble_info *i = ptr;
#ifdef INIT_DISASSEMBLE_INFO_THREE_ARG
    INIT_DISASSEMBLE_INFO(*i, stdout, FPRINTF_FUNCTION);
    i->endian = BFD_ENDIAN_LITTLE;
#else
    INIT_DISASSEMBLE_INFO(*i, stdout);
#endif
    i->mach = bfd_mach_i386_i386;
    if (s->p->code_base != NULL) {
	i->buffer = (bfd_byte *)s->p->code_base;
	i->buffer_vma = (bfd_vma)s->p->code_base;
    } else {
	i->buffer = (bfd_byte *)s->p->native.code_base;
	i->buffer_vma = (bfd_vma)s->p->native.code_base;
    }
    i->buffer_length = MAXLENGTH;
#ifdef HAVE_PRINT_INSN_I386
    return 1;
#else
    return 0;
#endif
}

extern int
x86_print_insn(dill_stream s, void *info_ptr, void *insn)
{
#ifdef HAVE_PRINT_INSN_I386
    return print_insn_i386((bfd_vma)insn, info_ptr);
#else
    return 0;
#endif
}

extern void null_func(){}
extern int
x86_count_insn(dill_stream s, int start, int end)
{
#ifdef HAVE_PRINT_INSN_I386
    struct disassemble_info i;
    int count;
    char *insn_ptr;
#ifdef INIT_DISASSEMBLE_INFO_THREE_ARG
    INIT_DISASSEMBLE_INFO(i, stdout, (fprintf_ftype) null_func);
    i.endian = BFD_ENDIAN_LITTLE;
#else
    INIT_DISASSEMBLE_INFO(i, stdout);
#endif
    i.mach = bfd_mach_i386_i386;
    if (s->p->code_base != NULL) {
	i.buffer = (bfd_byte *)s->p->code_base;
	i.buffer_vma = (bfd_vma)s->p->code_base;
    } else {
	i.buffer = (bfd_byte *)s->p->native.code_base;
	i.buffer_vma = (bfd_vma)s->p->native.code_base;
    }
    i.buffer_length = MAXLENGTH;
    count = 0;
    insn_ptr = (char*)i.buffer + start;
    while((long)insn_ptr < (long)i.buffer + end) {
	insn_ptr += print_insn_i386((bfd_vma)insn_ptr, &i);
	count++;
    }
    return count;
#else
   /* no print insn, just return the buffer length */
    return end - start;
#endif
}
#else
extern int
x86_count_insn(dill_stream s, int start, int end)
{   /* no print insn, just return the buffer length */
    return end - start;
}
extern int
x86_init_disassembly_info(dill_stream s, void * ptr){return 0;}
extern int x86_print_insn(dill_stream s, void *info_ptr, void *insn){return 0;}
#endif

static char *char_regs[] = {"AL", "CL", "DL", "BL", "AH", "CH", "DH", "BH"};
static char *short_regs[] = {"AX", "CX", "DX", "BX", "SP", "BP", "SI", "DI"};
static char *int_regs[] = {"EAX", "ECX", "EDX", "EBX", "ESP", "EBP", "ESI", "EDI"};

extern void
x86_print_reg(dill_stream s, int typ, int reg)
{
    x86_mach_info smi = (x86_mach_info) s->p->mach_info;
    switch(typ) {
    case DILL_C: case DILL_UC:
	if ((reg > 0) && (reg < sizeof(char_regs)/sizeof(char_regs[0]))) {
	    printf("%s", char_regs[reg]);
	    return;
	}
	break;
    case DILL_S: case DILL_US:
	if ((reg > 0) && (reg < sizeof(short_regs)/sizeof(short_regs[0]))) {
	    printf("%s", short_regs[reg]);
	    return;
	}
	break;
/*    case DILL_C: case DILL_UC: case DILL_S: case DILL_US:*/
    case DILL_I: case DILL_U: case DILL_L: case DILL_UL:
	if ((reg > 0) && (reg < sizeof(int_regs)/sizeof(int_regs[0]))) {
	    printf("%s", int_regs[reg]);
	    return;
	}
	break;
    case DILL_F: case DILL_D:
	if (smi->generate_SSE) {
	    printf("xmm%d", reg);
	} else {
	    printf("Fstack");
	}
	return;
    }
    printf("NoReg(%d)", reg);
}
	
