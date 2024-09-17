#include "config.h"

#undef NDEBUG
#include "assert.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dill.h"
#include "dill_internal.h"
#include "x86_64.h"

#define MOV32 0x89
#define Mod(x) ((x) << 6)
#define RegOp(x) ((7 & (x)) << 3)
#define RM(x) (7 & (x))
#define ModRM(mod, reg, rm) (Mod(mod) | RegOp(reg) | RM(rm))
#define SIB(scale, index, base)                                                \
    (((scale) << 6) | ((0x7 & (index)) << 3) | (0x7 & (base)))

#define INSN_OUT(s, i) printf("Bad opout, line %d\n", __LINE__)
#define x86_64_savei(s, imm)
#define x86_64_andi(s, dest, src, imm)                                         \
    x86_64_arith3i(s, 0x4, 0x2, dest, src, imm)

static void
x86_64_push_reg(dill_stream s, int src);
static void
x86_64_pop_reg(dill_stream s, int src);
#define x86_64_nop(s) BYTE_OUT1(s, 0x90)

static char* char_regs[] = {"AL", "CL", "DL", "BL", "AH", "CH", "DH", "BH"};
static char* short_regs[] = {"AX", "CX", "DX", "BX", "SP", "BP", "SI", "DI"};
static char* int_regs[] = {"EAX", "ECX", "EDX", "EBX",
                           "ESP", "EBP", "ESI", "EDI"};
char* long_regs[] = {"RAX", "RCX", "RDX", "RBX", "RSP", "RBP", "RSI", "RDI"};
char* float_regs[] = {"RAX", "RCX", "RDX", "RBX", "RSP", "RBP", "RSI", "RDI"};

#define IREG 0
#define FREG 1

#define _temp_reg EAX
#define _frame_reg EBP
#define roundup(a, b) ((a + (b - 1)) & (-b))
static void
x86_64_pmov(dill_stream s, int typ, int dest, int src);

static struct basic_type_info {
    char size;
    char align;
    char reg_type;
} type_info[] = {
    {1, 1, IREG},                                 /* C */
    {1, 1, IREG},                                 /* UC */
    {2, 2, IREG},                                 /* S */
    {2, 2, IREG},                                 /* US */
    {4, 4, IREG},                                 /* I */
    {4, 4, IREG},                                 /* U */
    {sizeof(uintptr_t), sizeof(uintptr_t), IREG}, /* UL */
    {sizeof(intptr_t), sizeof(intptr_t), IREG},   /* L */
    {sizeof(char*), sizeof(char*), IREG},         /* P */
    {sizeof(float), sizeof(float), FREG},         /* F */
    {sizeof(double), sizeof(double), FREG},       /* D */
    {0, 8, IREG},                                 /* V */
    {-1, 8, IREG},                                /* B */
    {sizeof(void*), sizeof(void*), IREG},         /* EC */
};

int x86_64_type_align[] = {
    1,                 /* C */
    1,                 /* UC */
    2,                 /* S */
    2,                 /* US */
    4,                 /* I */
    4,                 /* U */
    sizeof(uintptr_t), /* UL */
    sizeof(intptr_t),  /* L */
    sizeof(char*),     /* P */
    4,                 /* F */
    4,                 /* D */
    4,                 /* V */
    4,                 /* B */
    sizeof(char*),     /* EC */
};

int x86_64_type_size[] = {
    1,                 /* C */
    1,                 /* UC */
    2,                 /* S */
    2,                 /* US */
    4,                 /* I */
    4,                 /* U */
    sizeof(uintptr_t), /* UL */
    sizeof(intptr_t),  /* L */
    sizeof(char*),     /* P */
    4,                 /* F */
    8,                 /* D */
    4,                 /* V */
    0,                 /* B */
    sizeof(char*),     /* EC */
};

static void
BYTE_OUT2R(dill_stream s, int rex, int insn1, int insn2)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        *(tmp_ip + 2) = (unsigned char)insn2;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 2;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R(dill_stream s, int rex, int insn1)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 1;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT2IR(dill_stream s, int rex, int insn1, int insn2, int imm32)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        int tmp = imm32;
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        *(tmp_ip + 2) = (unsigned char)insn2;
        memcpy(tmp_ip + 3, &tmp, 4);
    } else {
        int tmp = imm32;
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        memcpy(tmp_ip + 2, &tmp, 4);
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 6;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT3IR(dill_stream s, int rex, int insn1, int insn2, int insn3, int imm32)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        int tmp = imm32;
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        memcpy(tmp_ip + 4, &tmp, 4);
    } else {
        int tmp = imm32;
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        memcpy(tmp_ip + 3, &tmp, 4);
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 7;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1IR(dill_stream s, int rex, int insn1, int imm32)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        int tmp = imm32;
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        memcpy(tmp_ip + 2, &tmp, 4);
    } else {
        int tmp = imm32;
        *(tmp_ip) = (unsigned char)insn1;
        memcpy(tmp_ip + 1, &tmp, 4);
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 5;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1LR(dill_stream s, int rex, int insn1, intptr_t imm64)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        intptr_t tmp = imm64;
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        memcpy(tmp_ip + 2, &tmp, 8);
    } else {
        intptr_t tmp = imm64;
        *(tmp_ip) = (unsigned char)insn1;
        memcpy(tmp_ip + 1, &tmp, 8);
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 9;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT3R(dill_stream s, int rex, int insn1, int insn2, int insn3)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 3;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT4R(dill_stream s, int rex, int insn1, int insn2, int insn3, int insn4)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *tmp_ip = (unsigned char)rex | 0x40;
        *(tmp_ip + 1) = (unsigned char)insn1;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 4;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R3I(dill_stream s,
             int insn1,
             int rex,
             int insn2,
             int insn3,
             int insn4,
             int imm)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *(tmp_ip + 0) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)rex | 0x40;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
        *((int*)(tmp_ip + 5)) = imm;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
        *((int*)(tmp_ip + 4)) = imm;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 8;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R4I(dill_stream s,
             int insn1,
             int rex,
             int insn2,
             int insn3,
             int insn4,
             int insn5,
             int imm)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *(tmp_ip + 0) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)rex | 0x40;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
        *(tmp_ip + 5) = (unsigned char)insn5;
        *((int*)(tmp_ip + 6)) = imm;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
        *(tmp_ip + 4) = (unsigned char)insn5;
        *((int*)(tmp_ip + 5)) = imm;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 9;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R3(dill_stream s, int insn1, int rex, int insn2, int insn3, int insn4)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *(tmp_ip + 0) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)rex | 0x40;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 4;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R4(dill_stream s,
            int insn1,
            int rex,
            int insn2,
            int insn3,
            int insn4,
            int insn5)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *(tmp_ip + 0) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)rex | 0x40;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
        *(tmp_ip + 5) = (unsigned char)insn5;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
        *(tmp_ip + 4) = (unsigned char)insn5;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 5;
    if (rex != 0)
        s->p->cur_ip++;
}

static void
BYTE_OUT1R5(dill_stream s,
            int insn1,
            int rex,
            int insn2,
            int insn3,
            int insn4,
            int insn5,
            int insn6)
{
    unsigned char* tmp_ip;
    if (s->p->cur_ip >= s->p->code_limit) {
        extend_dill_stream(s);
    }
    tmp_ip = (unsigned char*)s->p->cur_ip;
    if (rex != 0) {
        *(tmp_ip + 0) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)rex | 0x40;
        *(tmp_ip + 2) = (unsigned char)insn2;
        *(tmp_ip + 3) = (unsigned char)insn3;
        *(tmp_ip + 4) = (unsigned char)insn4;
        *(tmp_ip + 5) = (unsigned char)insn5;
        *(tmp_ip + 6) = (unsigned char)insn6;
    } else {
        *(tmp_ip) = (unsigned char)insn1;
        *(tmp_ip + 1) = (unsigned char)insn2;
        *(tmp_ip + 2) = (unsigned char)insn3;
        *(tmp_ip + 3) = (unsigned char)insn4;
        *(tmp_ip + 4) = (unsigned char)insn5;
        *(tmp_ip + 5) = (unsigned char)insn6;
    }
    if (s->dill_debug)
        dump_cur_dill_insn(s);
    s->p->cur_ip = ((char*)s->p->cur_ip) + 6;
    if (rex != 0)
        s->p->cur_ip++;
}

#define x86_64_movi(s, dest, src) x86_64_pmov(s, DILL_I, dest, src)
#define x86_64_movl(s, dest, src) x86_64_pmov(s, DILL_L, dest, src)
#define x86_64_movd(s, dest, src) x86_64_pmov(s, DILL_D, dest, src)
#define x86_64_movf(s, dest, src) x86_64_pmov(s, DILL_F, dest, src)

static void
x86_64_push_reg(dill_stream s, int src)
{
    int rex = 0;
    if (src > RDI)
        rex |= REX_B;
    BYTE_OUT1R(s, rex, 0x50 + (src & 0x7));
}

static void
x86_64_pop_reg(dill_stream s, int src)
{
    int rex = 0;
    if (src > RDI)
        rex |= REX_B;
    BYTE_OUT1R(s, rex, 0x58 + (src & 0x7));
}

static void
x86_64_rshai(dill_stream s, int dest, int src, int imm)
{
    int rex = REX_W;
    if (dest != src) {
        x86_64_movl(s, dest, src);
    }
    if (dest > RDI)
        rex |= REX_B;
    BYTE_OUT3R(s, rex, 0xc1, ModRM(0x3, 0x7, dest), imm & 0xff);
}

static void
x86_64_rshi(dill_stream s, int dest, int src, int imm)
{
    int rex = REX_W;
    if (dest != src) {
        x86_64_movl(s, dest, src);
    }
    if (dest > RDI)
        rex |= REX_B;
    BYTE_OUT3R(s, rex, 0xc1, ModRM(0x3, 0x5, dest), imm & 0xff);
}

static void
x86_64_lshi(dill_stream s, int dest, int src, int imm)
{
    int rex = REX_W;
    if (dest != src) {
        x86_64_movl(s, dest, src);
    }
    if (dest > RDI)
        rex |= REX_B;
    BYTE_OUT3R(s, rex, 0xc1, ModRM(0x3, 0x4, dest), imm & 0xff);
}

static void
x86_64_pmov(dill_stream s, int typ, int dest, int src)
{
    int rex = 0;
    if ((typ == DILL_L) || (typ == DILL_UL) || (typ == DILL_P))
        rex = REX_W;
    if (src > RDI)
        rex |= REX_R;
    if (dest > RDI)
        rex |= REX_B;
    if (src != dest) {
        switch (typ) {
        case DILL_D:
        case DILL_F:
            //	The 0x29 form of this instruction is not supported by
            // valgrind
            //      the 0x28 version reverses the operands.
            //	    BYTE_OUT1(s, 0x66);
            //	    BYTE_OUT3R(s, rex, 0x0f, 0x29, ModRM(0x3, src,
            // dest));
            rex = 0;
            if (dest > RDI)
                rex |= REX_R;
            if (src > RDI)
                rex |= REX_B;
            BYTE_OUT1(s, 0x66);
            BYTE_OUT3R(s, rex, 0x0f, 0x28, ModRM(0x3, dest, src));
            break;
        default:
            BYTE_OUT2R(s, rex, MOV32, ModRM(0x3, src, dest));
            break;
        }
    }
}

static void
x86_64_sxmov(dill_stream s, int typ, int src, int dest)
{
    int rex = REX_W;
    if (src > RDI)
        rex |= REX_R;
    if (dest > RDI)
        rex |= REX_B;
    if (src != dest) {
        switch (typ) {
        case DILL_C:
            BYTE_OUT3R(s, rex, 0x0f, 0xbe, ModRM(0x3, src, dest));
            break;
        case DILL_UC:
            BYTE_OUT3R(s, rex, 0x0f, 0xb6, ModRM(0x3, src, dest));
            break;
        case DILL_S:
            BYTE_OUT3R(s, rex, 0x0f, 0xbf, ModRM(0x3, src, dest));
            break;
        case DILL_US:
            BYTE_OUT3R(s, rex, 0x0f, 0xb7, ModRM(0x3, src, dest));
            break;
        }
    }
}

extern void
x86_64_farith2(dill_stream s, int b1, int typ, int dest, int src)
{
    /* this is fneg */
    int rex = 0;
    int op = 0xf3;
    if (src > XMM7)
        rex |= REX_B;
    if (dest > XMM7)
        rex |= REX_R;

    /* clear dest */
    {
        int rex1 = 0;
        if (dest > XMM7)
            rex1 = REX_R | REX_B;
        BYTE_OUT3R(s, rex1, 0x0f, 0x57,
                   ModRM(0x3, dest, dest));  // GSE really rex1?  Late fix.
    }
    if (typ == DILL_D)
        op = 0xf2;
    BYTE_OUT1R3(s, op, rex, 0x0f, 0x5c, ModRM(0x3, dest, src));
}

extern void
x86_64_farith(dill_stream s, int b1, int typ, int dest, int src1, int src2)
{
    int rex = 0;
    int op = 0xf3;

    if (src2 == dest) {
        /* protect src2 from overwrite (below) */
        x86_64_movd(s, XMM0, src2);
        src2 = XMM0;
    }
    if (src2 > XMM7)
        rex |= REX_B;
    if (dest > XMM7)
        rex |= REX_R;
    if (src1 != dest)
        x86_64_movd(s, dest, src1);
    if (typ == DILL_D)
        op = 0xf2;
    BYTE_OUT1R3(s, op, rex, 0x0f, b1, ModRM(0x3, dest, src2));
}

static void
x86_64_seti(dill_stream s, int r, int val)
{
    int rex = 0;
    if (r > RDI)
        rex |= REX_B;
    BYTE_OUT1IR(s, rex, 0xb8 + (0x7 & r), val);
}

static void
x86_64_setl(dill_stream s, int r, IMM_TYPE val)
{
    int rex = REX_W;
    if (r > RDI)
        rex |= REX_B;
    BYTE_OUT1LR(s, rex, 0xb8 + ((0x7) & r), val);
}

extern void
x86_64_setp(dill_stream s, int type, int junk, int r, void* val)
{
    int rex = REX_W;
    union {
        intptr_t l;
        void* a;
    } a;
    a.a = val;
    if (r > RDI)
        rex |= REX_B;
    BYTE_OUT1LR(s, rex, 0xb8 + ((0x7) & r), a.l);
}

extern int
x86_64_local(dill_stream s, int type)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;

    s->p->used_frame++;
    smi->act_rec_size += roundup(type_info[type].size, smi->stack_align);
    return (-smi->act_rec_size) + smi->stack_constant_offset;
}

extern int
x86_64_localb(dill_stream s, int size)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    s->p->used_frame++;
    smi->act_rec_size = roundup(smi->act_rec_size, size);

    smi->act_rec_size += roundup(size, smi->stack_align);
    return (-smi->act_rec_size) + smi->stack_constant_offset;
}

extern int
x86_64_local_op(dill_stream s, int flag, int val)
{
    int size = val;
    if (flag == 0) {
        size = type_info[val].size;
    }
    return x86_64_localb(s, size);
}

#define BEGIN_FLOAT_SAVE 64
extern void
x86_64_save_restore_op(dill_stream s, int save_restore, int type, int reg)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int offset = 0;
    switch (type) {
    case DILL_D:
    case DILL_F:
        offset = reg * smi->stack_align + BEGIN_FLOAT_SAVE;
        break;
    default:
        if ((reg == RBX) || ((reg >= R12) && (reg <= R15))) {
            /* RBX, R12-R15 inclusive are callee-saved, caller doesn't save
             */
            return;
        }
        switch (reg) {
        case RDI:
            offset = 0;
            break;
        case RSI:
            offset = 8;
            break;
        case RDX:
            offset = 16;
            break;
        case RCX:
            offset = 24;
            break;
        case R8:
            offset = 32;
            break;
        case R9:
            offset = 40;
            break;
        case R10:
            offset = 48;
            break;
        case R11:
            offset = 56;
            break;
        default:
            printf("x86_64 save_restore called with int reg %d\n", reg);
            assert(0); /* shouldn't happen */
        }
    }
    if (save_restore == 0) { /* save */
        x86_64_pstorei(s, type, 0, reg, _frame_reg, smi->save_base + offset);
    } else { /* restore */
        x86_64_ploadi(s, type, 0, reg, _frame_reg, smi->save_base + offset);
    }
    s->p->used_frame++;
}

/*
 *     This is linux convention.  Windows is different.
 *
 * register   use  			       preserved across function calls
 *  %rax temporary register; with variable arguments
 *        passes information about the
 *	  number of SSE registers used; 1st return
 *	  register				            	No
 *  %rbx callee-saved register; optionally used
 *        as base pointer				    	Yes
 *  %rcx used to pass 4th integer argument to
 *	  functions					    	No
 *  %rdx used to pass 3rd argument to functions;
 *	  2nd return register				    	No
 *  %rsp stack pointer 					    	Yes
 *  %rbp callee-saved register; optionally used
 *       as frame pointer				    	Yes
 *  %rsi used to pass 2nd argument to functions 		No
 *  %rdi used to pass 1st argument to functions			No
 *  %r8  used to pass 5th argument to functions 		No
 *  %r9  used to pass 6th argument to functions			No
 *  %r10 temporary register, used for passing a
 *       function?s static chain pointer 			No
 * %r11 temporary register					No
 * %r12-r15    callee-saved registers				Yes
 * %xmm0-%xmm1 used to pass and return floating point
 * 			arguments 				No
 * %xmm2-%xmm7 used to pass floating point arguments		No
 * %xmm8-%xmm15 temporary registers 				No
 * %mmx0-%mmx7 temporary registers 				No
 * %st0 temporary register; used to return long double args	No
 * %st1 temporary registers; used to return long double args 	No
 * %st2-%st7 temporary registers 				No
 * %fs Reserved for system use (as thread specific data register)

Figure 3.33: Register Save Area   (modified -- GSE )
Register Offset
%rdi 0
%rsi 8
%rdx 16
%rcx 24
%r8 32
%r9 40
%r10  48
%r11  56
%xmm0 64	BEGIN_FLOAT_SAVE
%xmm1 68
. . .
%xmm15 288

*/
#ifndef USE_WINDOWS_CALLS
static int arg_regs[] = {RDI, RSI, RDX, RCX, R8, R9};
#else
static int arg_regs[] = {RCX, RDX, R8, R9};
#endif

static void
save_required_regs(dill_stream s, int force)
{
    if (s->p->call_table.call_count != 0)
        force++;

#ifdef OSX
    if (force) {
        dill_andii(s, ESP, ESP, -16); /* make sure it's multiple of 16 */
    }
#endif

    /* callee is supposed to save these */
    if (force || dill_wasused(&s->p->var_i, EBX) ||
        dill_wasused(&s->p->tmp_i, EBX)) {
        x86_64_push_reg(s, EBX);
    }
    if (force || dill_wasused(&s->p->var_i, R12) ||
        dill_wasused(&s->p->tmp_i, R12)) {
        x86_64_push_reg(s, R12);
    }
    if (force || dill_wasused(&s->p->var_i, R13) ||
        dill_wasused(&s->p->tmp_i, R13)) {
        x86_64_push_reg(s, R13);
    }
    if (force || dill_wasused(&s->p->var_i, R14) ||
        dill_wasused(&s->p->tmp_i, R14)) {
        x86_64_push_reg(s, R14);
    }
    if (force || dill_wasused(&s->p->var_i, R15) ||
        dill_wasused(&s->p->tmp_i, R15)) {
        x86_64_push_reg(s, R15);
    }
}

static int
generate_prefix_code(dill_stream s, int force, int ar_size, dill_reg* arglist)
{
    int end, start = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    int i, int_arg_count, float_arg_count;
    arg_info_list args = s->p->c_param_args;

    if ((s->p->c_param_count > 0) || s->p->used_frame ||
        (s->p->call_table.call_count != 0) || force) {
        x86_64_push_reg(s, EBP);
        x86_64_movl(s, EBP, ESP);
    }

    if (force || s->p->used_frame || (s->p->call_table.call_count != 0)) {
        /* do local space reservation */
        dill_subli(s, ESP, ESP, ar_size);
    }

    save_required_regs(s, force);

    int_arg_count = float_arg_count = 0; /* reset */
    for (i = 0; i < s->p->c_param_count; i++) {
        int came_in_a_reg = 0;
        switch (args[i].type) {
        case DILL_D:
        case DILL_F:
#ifdef USE_WINDOWS_CALLS
            came_in_a_reg = (float_arg_count + int_arg_count) < 4;
            float_arg_count++;
#else
            came_in_a_reg = (float_arg_count++ < 8);
#endif
            break;
        default:
#ifdef USE_WINDOWS_CALLS
            came_in_a_reg = (float_arg_count + int_arg_count) < 4;
            int_arg_count++;
#else
            came_in_a_reg = (int_arg_count++ < 6);
#endif
            break;
        }
        if (came_in_a_reg && (args[i].in_reg == -1)) {
            /* not enough regs for this, store it to the stack */
            int real_offset = args[i].offset;
            x86_64_pstorei(s, args[i].type, 0, arg_regs[int_arg_count - 1], EBP,
                           real_offset);
            continue;
        }
        if (came_in_a_reg && args[i].is_register) {
            switch (args[i].type) {
            case DILL_D:
            case DILL_F:
#ifdef USE_WINDOWS_CALLS
                x86_64_movd(s, args[i].in_reg,
                            float_arg_count + int_arg_count - 1);
#else
                x86_64_movd(s, args[i].in_reg, float_arg_count - 1);
#endif
                break;
            case DILL_UC:
            case DILL_C:
            case DILL_US:
            case DILL_S:
#ifdef USE_WINDOWS_CALLS
                x86_64_sxmov(s, args[i].type, args[i].in_reg,
                             arg_regs[float_arg_count + int_arg_count - 1]);
#else
                x86_64_sxmov(s, args[i].type, args[i].in_reg,
                             arg_regs[int_arg_count - 1]);
#endif
                break;
            default:
#ifdef USE_WINDOWS_CALLS
                x86_64_movl(s, args[i].in_reg,
                            arg_regs[int_arg_count + float_arg_count - 1]);
#else
                x86_64_movl(s, args[i].in_reg, arg_regs[int_arg_count - 1]);
#endif
                break;
            }
            continue;
        }
        if (args[i].is_register) {
            /* general offset from fp*/
            int real_offset = args[i].offset;
            x86_64_ploadi(s, args[i].type, 0, args[i].in_reg, EBP, real_offset);
        }
    }

    end = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    return end - start;
}

extern void
x86_64_proc_start(dill_stream s,
                  char* subr_name,
                  int arg_count,
                  arg_info_list args,
                  dill_reg* arglist)
{
    int i, int_arg_count, float_arg_count;
    int cur_arg_offset = 0;

    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    smi->pending_prefix = 0;
    smi->last_proc_ret_end = 0;

    /* leave some space */ x86_64_local(s, DILL_D);
    smi->conversion_word = x86_64_local(s, DILL_D);
    smi->fcu_word = x86_64_local(s, DILL_I);
    smi->save_base = x86_64_localb(s, BEGIN_FLOAT_SAVE + 16 * 16);
    s->p->used_frame = 0;

    cur_arg_offset = 16;
    int_arg_count = 0;
    float_arg_count = 0;
    for (i = 0; i < arg_count; i++) {
        args[i].in_reg = -1;
        args[i].out_reg = -1;
        args[i].is_register = 0;
        if (arglist != NULL)
            arglist[i] = -1;
#ifdef USE_WINDOWS_CALLS
        if ((args[i].type != DILL_F) && (args[i].type != DILL_D) &&
            (int_arg_count + float_arg_count < 4)) {
            cur_arg_offset +=
                roundup(type_info[(int)args[i].type].size, smi->stack_align);
#else
        if ((args[i].type != DILL_F) && (args[i].type != DILL_D) &&
            (int_arg_count < 6)) {
#endif

            args[i].offset = smi->save_base + int_arg_count * 8;
            int_arg_count++;
#ifdef USE_WINDOWS_CALLS
            if (int_arg_count + float_arg_count <= 4) {
#else
            if (int_arg_count <= 4) {
#endif
                dill_reg tmp_reg;
                if (dill_raw_getreg(s, &tmp_reg, args[i].type, DILL_VAR)) {
                    args[i].in_reg = tmp_reg;
                    if (arglist != NULL)
                        arglist[i] = tmp_reg;
                    args[i].is_register = 1;
                }
            } else {
#ifdef USE_WINDOWS_CALLS
                if (int_arg_count + float_arg_count <=
                    (sizeof(arg_regs) / sizeof(arg_regs[0]))) {
#else
                if (int_arg_count <= (sizeof(arg_regs) / sizeof(arg_regs[0]))) {
#endif
                    args[i].is_register = 1;
                    args[i].in_reg = arg_regs[int_arg_count - 1];
                }
            }
#ifdef USE_WINDOWS_CALLS
        } else if (((args[i].type == DILL_F) || (args[i].type == DILL_D)) &&
                   (int_arg_count + float_arg_count < 4)) {
            cur_arg_offset +=
                roundup(type_info[(int)args[i].type].size, smi->stack_align);
#else
        } else if (((args[i].type == DILL_F) || (args[i].type == DILL_D)) &&
                   (float_arg_count < 8)) {
#endif

            args[i].offset =
                smi->save_base + BEGIN_FLOAT_SAVE + float_arg_count * 8;
            float_arg_count++;
            if (float_arg_count <= 4) {
                dill_reg tmp_reg;
                if (dill_raw_getreg(s, &tmp_reg, args[i].type, DILL_VAR)) {
                    args[i].in_reg = tmp_reg;
                    if (arglist != NULL)
                        arglist[i] = tmp_reg;
                    args[i].is_register = 1;
                }
            }
        } else {
            args[i].offset = cur_arg_offset;
            cur_arg_offset +=
                roundup(type_info[(int)args[i].type].size, smi->stack_align);
        }
    }
    /* make local space reservation constant big so we have a word to patch */
    /* use the nop op code so that if we don't use all of it we get nops */
    (void)generate_prefix_code(s, 1 /* force */, 0x909090, arglist);

    smi->backpatch_offset = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
}

static void
x86_64_proc_ret(dill_stream s)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int force = 0;
    if (s->p->call_table.call_count != 0)
        force++;

    /* last thing was a ret and there's not a label here, we don't need another
     * ret */
    if ((smi->last_proc_ret_end == s->p->cur_ip) && !dill_is_label_mark(s))
        return;

    if (force || dill_wasused(&s->p->var_i, R15) ||
        dill_wasused(&s->p->tmp_i, R15)) {
        x86_64_pop_reg(s, R15);
    }
    if (force || dill_wasused(&s->p->var_i, R14) ||
        dill_wasused(&s->p->tmp_i, R14)) {
        x86_64_pop_reg(s, R14);
    }
    if (force || dill_wasused(&s->p->var_i, R13) ||
        dill_wasused(&s->p->tmp_i, R13)) {
        x86_64_pop_reg(s, R13);
    }
    if (force || dill_wasused(&s->p->var_i, R12) ||
        dill_wasused(&s->p->tmp_i, R12)) {
        x86_64_pop_reg(s, R12);
    }
    if (force || dill_wasused(&s->p->var_i, EBX) ||
        dill_wasused(&s->p->tmp_i, EBX)) {
        x86_64_pop_reg(s, EBX);
    }
    if ((s->p->c_param_count > 0) || s->p->used_frame ||
        (s->p->call_table.call_count != 0)) {
        x86_64_movl(s, ESP, EBP);
        x86_64_pop_reg(s, EBP);
    }
    BYTE_OUT1(s, 0xc3);
    smi->last_proc_ret_end = s->p->cur_ip;
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

static void
x86_64_clear(dill_stream s, int dest)
{
    int rex = REX_W;
    if (dest > RDI)
        rex |= REX_B | REX_R;
    BYTE_OUT2R(s, rex, 0x33, ModRM(0x3, dest, dest)); /* xor dest, dest */
}

extern void
x86_64_ploadi(dill_stream s,
              int type,
              int junk,
              int dest,
              int src,
              IMM_TYPE offset)
{
    unsigned char opcode = ld_opcodes[type];
    int tmp_dest = dest;
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int rex = 0;
    int float_op = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
    }
    switch (type) {
    case DILL_F:
        float_op = 0xf3;
        break;
    case DILL_D:
        float_op = 0xf2;
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
            x86_64_clear(s, tmp_dest);
        }
        break;
    case DILL_S:
    case DILL_US:
        if (type == DILL_US) {
            /* clear high bytes */
            if (src == tmp_dest) {
                /* don't destroy source */
                tmp_dest = EAX;
            }
            x86_64_clear(s, tmp_dest);
        }
        BYTE_OUT1(s, 0x66);
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
    /* fall through */
    default:
        break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
        smi->pending_prefix = 0;
    }
    if (src > RDI)
        rex |= REX_B;
    if (tmp_dest > RDI)
        rex |= REX_R;
    if (((0x7 & src) == RSP) &&
        (((offset & 0xffffffff80000000) == 0) ||
         ((offset & 0xffffffff80000000) == 0xffffffff80000000))) {
        /* must use SIB because ModRM has a discontinuity */
        if (offset == 0) {
            if (float_op != 0) {
                BYTE_OUT1R4(s, float_op, rex, 0x0f, 0x10,
                            ModRM(0x0, tmp_dest, 0x4), SIB(0, 4, src));
            } else {
                BYTE_OUT3R(s, rex, opcode, ModRM(0x0, tmp_dest, 0x4),
                           SIB(0, 4, src));
            }
        } else if (((intptr_t)offset <= 127) && ((intptr_t)offset > -128)) {
            if (float_op != 0) {
                BYTE_OUT1R5(s, float_op, rex, 0x0f, 0x10,
                            ModRM(0x1, tmp_dest, src), SIB(0, ESP, 0x4),
                            offset & 0xff);
            } else {
                BYTE_OUT4R(s, rex, opcode, ModRM(0x1, tmp_dest, src),
                           SIB(0, ESP, 0x4), offset & 0xff);
            }
        } else {
            if (float_op != 0) {
                BYTE_OUT1R4I(s, float_op, rex, 0x0f, 0x10,
                             ModRM(0x2, tmp_dest, 0x4), SIB(0, 4, src),
                             (int)offset);
            } else {
                BYTE_OUT3IR(s, rex, opcode, ModRM(0x2, tmp_dest, 0x4),
                            SIB(0, 4, src), (int)offset);
            }
        }
    } else {
        if ((offset == 0) &&
            ((src & 0x7) != 5)) { /* avoid discontinuity in ModRM */
            if (float_op != 0) {
                BYTE_OUT1R3(s, float_op, rex, 0x0f, 0x10,
                            ModRM(0x0, tmp_dest, src));
            } else {
                BYTE_OUT2R(s, rex, opcode, ModRM(0x0, tmp_dest, src));
            }
        } else if (((intptr_t)offset <= 127) && ((intptr_t)offset > -128)) {
            if (float_op != 0) {
                BYTE_OUT1R4(s, float_op, rex, 0x0f, 0x10,
                            ModRM(0x1, tmp_dest, src), offset & 0xff);
            } else {
                BYTE_OUT3R(s, rex, opcode, ModRM(0x1, tmp_dest, src),
                           offset & 0xff);
            }
        } else if (((offset & 0xffffffff80000000) == 0) ||
                   ((offset & 0xffffffff80000000) == 0xffffffff80000000)) {
            if (float_op != 0) {
                BYTE_OUT1R3I(s, float_op, rex, 0x0f, 0x10,
                             ModRM(0x2, tmp_dest, src), (int)offset);
            } else {
                BYTE_OUT2IR(s, rex, opcode, ModRM(0x2, tmp_dest, src),
                            (int)offset);
            }
        } else {
            /* really big offset */
            if (dest != src) {
                /* use dest since it'll be destroyed by the load */
                x86_64_setl(s, dest, offset);
                x86_64_pload(s, type, 0, dest, src, dest);
            } else {
                /* destroy src, but since it's the same as dest, it's lost
                 * anyway */
                x86_64_arith3i(s, 0, DILL_L, src, src, offset);
                x86_64_ploadi(s, type, 0, dest, src, 0);
            }
        }
    }
    switch (type) {
    case DILL_C:
        x86_64_lshi(s, dest, tmp_dest, 56);
        x86_64_rshai(s, dest, dest, 56);
        break;
    case DILL_S:
        x86_64_lshi(s, dest, tmp_dest, 48);
        x86_64_rshai(s, dest, dest, 48);
        break;
    case DILL_UC:
    case DILL_US:
        if (dest != tmp_dest)
            x86_64_movi(s, dest, tmp_dest);
        break;
    }
}

extern void
x86_64_pload(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    unsigned char opcode = ld_opcodes[type];
    int tmp_dest = dest;
    int clear_high_bytes_after_load = 0;

    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int rex = 0;
    int float_op = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
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
            if ((src1 == tmp_dest) || (src2 == tmp_dest)) {
                /* don't destroy source */
                tmp_dest = EAX;
            }
            if ((src1 == EAX) || (src2 == EAX)) {
                clear_high_bytes_after_load = 1;
            } else {
                x86_64_clear(s, tmp_dest);
            }
        }
        break;
    case DILL_S:
    case DILL_US:
        if (type == DILL_US) {
            /* clear high bytes */
            if ((src1 == tmp_dest) || (src2 == tmp_dest)) {
                /* don't destroy source */
                tmp_dest = EAX;
            }
            if ((src1 == EAX) || (src2 == EAX)) {
                clear_high_bytes_after_load = 1;
            } else {
                x86_64_clear(s, tmp_dest);
            }
        }
        BYTE_OUT1(s, 0x66);
        break;
    case DILL_F:
        float_op = 0xf3;
        break;
    case DILL_D:
        float_op = 0xf2;
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        if (smi->stack_align == 4) {
            type = DILL_I;
        }
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
        smi->pending_prefix = 0;
    }
    if ((src2 & 0x7) == 0x5) { /* avoid SIB discontinuity */
        int tmp = src1;
        src1 = src2;
        src2 = tmp;
    }
    if (src1 > RDI)
        rex |= REX_X;
    if (src2 > RDI)
        rex |= REX_B;
    if (tmp_dest > RDI)
        rex |= REX_R;
    if (float_op != 0) {
        BYTE_OUT1R4(s, float_op, rex, 0x0f, 0x10, ModRM(0x0, tmp_dest, 0x4),
                    SIB(0, src1, src2));
    } else {
        BYTE_OUT3R(s, rex, opcode, ModRM(0x0, tmp_dest, 0x4),
                   SIB(0, src1, src2));
    }
    switch (type) {
    case DILL_C:
        x86_64_lshi(s, dest, tmp_dest, 56);
        x86_64_rshai(s, dest, dest, 56);
        break;
    case DILL_S:
        x86_64_lshi(s, dest, tmp_dest, 48);
        x86_64_rshai(s, dest, dest, 48);
        break;
    case DILL_UC:
    case DILL_US:
        if (dest != tmp_dest) {
            if (clear_high_bytes_after_load) {
                int bits = 56;
                if (type == DILL_US)
                    bits = 48;
                x86_64_lshi(s, dest, tmp_dest, bits);
                x86_64_rshi(s, dest, dest, bits);
            } else {
                x86_64_movi(s, dest, tmp_dest);
            }
        }
        break;
    }
}

extern void
x86_64_pbsloadi(dill_stream s,
                int type,
                int junk,
                int dest,
                int src,
                IMM_TYPE offset)
{
    int rex = 0;
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
    x86_64_ploadi(s, ltype, junk, dest, src, offset);
    if (dest > RDI)
        rex |= REX_B;
    switch (type) {
    case DILL_F:
        rex = 0;
        if (fdest > RDI)
            rex |= REX_R;
        x86_64_bswap(s, 0, DILL_I, EAX, EAX);
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, fdest, EAX));
        break;
    case DILL_D:
        rex = REX_W;
        if (fdest > RDI)
            rex |= REX_R;
        x86_64_bswap(s, 0, DILL_L, EAX, EAX);
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, fdest, EAX));
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        rex |= REX_W;
        /* falling through */
    case DILL_I:
    case DILL_U:
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        break;
    case DILL_S:
    case DILL_US:
        /* byteswap 32 bits and shift down 16 */
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        x86_64_rshi(s, dest, dest, 16);
        break;
    case DILL_C:
    case DILL_UC:
        break;
    }
}

extern void
x86_64_pbsload(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    int rex = 0;
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
    if (dest > RDI)
        rex |= REX_B;
    x86_64_pload(s, ltype, junk, dest, src1, src2);
    switch (type) {
    case DILL_F:
        rex = 0;
        x86_64_bswap(s, 0, DILL_I, EAX, EAX);
        if (fdest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, fdest, EAX));
        break;
    case DILL_D:
        rex = REX_W;
        if (fdest > RDI)
            rex |= REX_R;
        x86_64_bswap(s, 0, DILL_L, EAX, EAX);
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, fdest, EAX));
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        rex |= REX_W;
    case DILL_I:
    case DILL_U:
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        break;
    case DILL_S:
    case DILL_US:
        /* byteswap 32 bits and shift down 16 */
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        x86_64_rshi(s, dest, dest, 16);
        break;
    case DILL_C:
    case DILL_UC:
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
x86_64_pstorei(dill_stream s,
               int type,
               int junk,
               int dest,
               int src,
               IMM_TYPE offset)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int rex = 0, float_op = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
    }
    switch (type) {
    case DILL_C:
    case DILL_UC:
        if (dest >= ESP) {
            /* movb doesn't work for big regs, move to eax */
            x86_64_movi(s, EAX, dest);
            dest = EAX;
        }
        break;
    }
    switch (type) {
    case DILL_F:
        float_op = 0xf3;
        break;
    case DILL_D:
        float_op = 0xf2;
        break;
    case DILL_S:
    case DILL_US:
        BYTE_OUT1(s, 0x66);
        break;
    default:
        break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
        smi->pending_prefix = 0;
    }
    if (dest > RDI)
        rex |= REX_R;
    if (src > RDI)
        rex |= REX_B;
    if (((src & 0x7) == ESP) &&
        (((offset & 0xffffffff80000000) == 0) ||
         ((offset & 0xffffffff80000000) == 0xffffffff80000000))) {
        if (((intptr_t)offset <= 127) && ((intptr_t)offset > -128)) {
            if (float_op != 0) {
                BYTE_OUT1R5(s, float_op, rex, 0x0f, 0x11, ModRM(0x1, dest, src),
                            SIB(0x0, ESP, 0x4), offset & 0xff);
            } else {
                BYTE_OUT4R(s, rex, st_opcodes[type], ModRM(0x1, dest, src),
                           SIB(0x0, ESP, 0x4), offset & 0xff);
            }
        } else {
            if (float_op != 0) {
                BYTE_OUT1R4I(s, float_op, rex, 0x0f, 0x11,
                             ModRM(0x2, dest, 0x4), SIB(0, 4, src),
                             (int)offset);
            } else {
                BYTE_OUT3IR(s, rex, st_opcodes[type], ModRM(0x2, dest, 0x4),
                            SIB(0, 4, src), (int)offset);
            }
        }
    } else {
        if ((offset == 0) &&
            ((src & 0x7) != 5)) { /* avoid discontinuity in ModRM */
            if (float_op != 0) {
                BYTE_OUT1R3(s, float_op, rex, 0x0f, 0x11,
                            ModRM(0x0, dest, src));
            } else {
                BYTE_OUT2R(s, rex, st_opcodes[type], ModRM(0x0, dest, src));
            }
        } else if (((intptr_t)offset <= 127) && ((intptr_t)offset > -128)) {
            if (float_op != 0) {
                BYTE_OUT1R4(s, float_op, rex, 0x0f, 0x11, ModRM(0x1, dest, src),
                            offset & 0xff);
            } else {
                BYTE_OUT3R(s, rex, st_opcodes[type], ModRM(0x1, dest, src),
                           offset & 0xff);
            }
        } else if (((offset & 0xffffffff80000000) == 0) ||
                   ((offset & 0xffffffff80000000) == 0xffffffff80000000)) {
            /* safe INT offset using only low 31 bits */
            if (float_op != 0) {
                BYTE_OUT1R3I(s, float_op, rex, 0x0f, 0x11,
                             ModRM(0x2, dest, src), (int)offset);
            } else {
                BYTE_OUT2IR(s, rex, st_opcodes[type], ModRM(0x2, dest, src),
                            (int)offset);
            }
        } else {
            x86_64_push_reg(s, src);
            x86_64_arith3i(s, 0, DILL_L, src, src, offset);
            x86_64_pstorei(s, type, 0, dest, src, 0);
            x86_64_pop_reg(s, src);
        }
    }
}

extern void
x86_64_pstore(dill_stream s, int type, int junk, int dest, int src1, int src2)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int rex = 0, float_op = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
    }
    switch (type) {
    case DILL_F:
        float_op = 0xf3;
        break;
    case DILL_D:
        float_op = 0xf2;
        break;
    case DILL_S:
    case DILL_US:
        BYTE_OUT1(s, 0x66);
        break;
    default:
        break;
    }
    if (smi->pending_prefix != 0) {
        BYTE_OUT1(s, smi->pending_prefix);
        smi->pending_prefix = 0;
    }
    if ((src2 & 0x7) == 0x5) { /* avoid SIB discontinuity */
        int tmp = src1;
        src1 = src2;
        src2 = tmp;
    }
    if (src1 > RDI)
        rex |= REX_X;
    if (src2 > RDI)
        rex |= REX_B;
    if (dest > RDI)
        rex |= REX_R;
    if (float_op != 0) {
        BYTE_OUT1R4(s, float_op, rex, 0x0f, 0x11, ModRM(0x0, dest, 0x4),
                    SIB(0, src1, src2));
    } else {
        BYTE_OUT3R(s, rex, st_opcodes[type], ModRM(0x0, dest, 0x4),
                   SIB(0, src1, src2));
    }
}

extern double
dill_x86_64_hidden_ULtoD(size_t a);
extern size_t
dill_x86_64_hidden_DtoUL(double a);

extern void
x86_64_div(dill_stream s, int op3, int op, int dest, int src1, int src2)
{
}

extern void
x86_64_divi(dill_stream s, int op3, int op, int dest, int src, IMM_TYPE imm)
{
}

extern void
x86_64_mov(dill_stream s, int type, int junk, int dest, int src)
{
    x86_64_pmov(s, type, dest, src);
}

extern void
x86_64_arith3(dill_stream s, int op, int typ, int dest, int src1, int src2)
{
    int commut = (op != 0x2b); /* subtract not commutative */
    int rex = 0;
    if ((typ == DILL_L) || (typ == DILL_UL) || (typ == DILL_P)) {
        rex = REX_W;
    }

    if (commut && (dest == src1)) {
        if (dest > RDI)
            rex |= REX_R;
        if (src2 > RDI)
            rex |= REX_B;
        BYTE_OUT2R(s, rex, op, ModRM(0x3, dest, src2));
    } else if (commut && (dest == src2)) {
        if (dest > RDI)
            rex |= REX_R;
        if (src1 > RDI)
            rex |= REX_B;
        BYTE_OUT2R(s, rex, op, ModRM(0x3, dest, src1));
    } else if (dest == src2) {
        assert(op == 0x2b); /* must be subtract */
        int rex1 = rex, rex2 = rex;
        if (dest > RDI)
            rex1 |= REX_B;
        BYTE_OUT2R(s, rex1, 0xf7, ModRM(0x3, 0x3, dest)); /* neg src2/dest */
        if (dest > RDI)
            rex2 |= REX_R;
        if (src1 > RDI)
            rex2 |= REX_B;
        BYTE_OUT2R(s, rex2, 0x03, ModRM(0x3, dest, src1)); /* add src1, dest */
    } else {
        int rex1 = rex, rex2 = rex;
        if (dest > RDI)
            rex1 |= REX_B;
        if (src1 > RDI)
            rex1 |= REX_R;
        BYTE_OUT2R(s, rex1, MOV32, ModRM(0x3, src1, dest));
        if (src2 > RDI)
            rex2 |= REX_B;
        if (dest > RDI)
            rex2 |= REX_R;
        BYTE_OUT2R(s, rex2, op, ModRM(0x3, dest, src2));
    }
}

extern void
x86_64_arith2(dill_stream s, int op, int subop, int dest, int src)
{
    if (op == 0) {
        int tmp_dest = dest;
        int rex1 = REX_W, rex2 = REX_W;
        /* must be not */
        if (dest >= ESP) {
            tmp_dest = EAX;
        }
        if (src > RDI)
            rex1 |= REX_B;
        BYTE_OUT3R(s, rex1, 0x83, ModRM(0x3, 0x7, src), 0); /* cmp */
        x86_64_seti(s, tmp_dest, 0);
        if (src > RDI)
            rex2 |= REX_R;
        if (tmp_dest > RDI)
            rex2 |= REX_B;
        BYTE_OUT3R(s, rex2, 0x0f, 0x94,
                   ModRM(0x3, src, tmp_dest)); /* sete dest */
        if (tmp_dest != dest) {
            x86_64_movi(s, dest, tmp_dest);
        }
    } else {
        int rex = REX_W;
        if (src != dest) {
            int rex1 = rex;
            if (dest > RDI)
                rex1 |= REX_B;
            if (src > RDI)
                rex1 |= REX_R;
            BYTE_OUT2R(s, rex1, MOV32, ModRM(0x3, src, dest));
        }
        if (dest > RDI)
            rex |= REX_B;
        BYTE_OUT2R(s, rex, op, ModRM(0x3, subop, dest));
    }
}

extern void
x86_64_bswap(dill_stream s, int junk, int typ, int dest, int src)
{
    int rex = 0;
    switch (typ) {
    case DILL_F:
        rex = 0;
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x7e, ModRM(0x3, src, EAX));
        x86_64_bswap(s, 0, DILL_I, EAX, EAX);
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
        break;
    case DILL_D:
        rex = REX_W;
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x7e, ModRM(0x3, src, EAX));
        x86_64_bswap(s, 0, DILL_L, EAX, EAX);
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        rex |= REX_W;
    case DILL_I:
    case DILL_U:
        if (dest > RDI)
            rex |= REX_B;
        if (dest != src) {
            x86_64_movl(s, dest, src);
        }
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        break;
    case DILL_S:
    case DILL_US:
        if (dest > RDI)
            rex |= REX_B;
        if (dest != src) {
            x86_64_movl(s, dest, src);
        }
        /* byteswap 32 bits and shift down 16 */
        BYTE_OUT2R(s, rex, 0x0f, 0xc8 + (dest & 0x7)); /* byteswap dest */
        x86_64_rshi(s, dest, dest, 16);
        break;
    case DILL_C:
    case DILL_UC:
        break;
    }
}

extern void
x86_64_mul(dill_stream s,
           int sign,
           int imm,
           int dest,
           int src1,
           IMM_TYPE imm_src2)
{
    int src2 = (int)imm_src2;
    int rex = REX_W;
    /* make src1 be EAX */
    if (dest != EAX) {
        x86_64_push_reg(s, EAX);
    }
    if (dest != EDX) {
        x86_64_push_reg(s, EDX);
    }

    if ((src2 == EAX) && !imm) {
        int tmp = src2;
        src1 = src2;
        src2 = tmp;
    }
    if (src1 != EAX) {
        x86_64_movl(s, EAX, src1);
    }
    if (imm == 0) {
        if (src2 > RDI)
            rex |= REX_B;
        BYTE_OUT2R(s, rex, 0xf7, ModRM(0x3, sign ? 0x5 : 0x4, src2));
    } else {
        /* src2 is really immediate */
        if (sign && ((imm_src2 & 0xffffffff80000000) == 0)) {
            BYTE_OUT2IR(s, rex, 0x69, ModRM(0x3, 0, EAX), (int)imm_src2);
        } else {
            x86_64_setl(s, EDX, imm_src2);
            BYTE_OUT2R(s, rex, 0xf7, ModRM(0x3, 0x4, EDX));
        }
    }
    if (dest != EDX) {
        x86_64_pop_reg(s, EDX);
    }
    if (dest != EAX) {
        x86_64_movl(s, dest, EAX);
        x86_64_pop_reg(s, EAX);
    }
}

extern void
x86_64_div_modi(dill_stream s,
                int div,
                int type,
                int dest,
                int src1,
                IMM_TYPE imm)
{
    x86_64_push_reg(s, EBP);
    x86_64_setl(s, EBP, imm);
    x86_64_div_mod(s, div, type, dest, src1, EBP);
    x86_64_pop_reg(s, EBP);
}

extern void
x86_64_div_mod(dill_stream s, int div, int type, int dest, int src1, int src2)
{
    int tmp_src2 = src2;
    int rex = 0;
    int sign = ((type == DILL_I) || (type == DILL_L));

    if ((type == DILL_UL) || (type == DILL_L)) {
        rex = REX_W;
    }
    /*  Well, this is awkward.  divi must use EDX and EAX for src1, destroying
     * both */

    /* save eax and edx, (unless one is the dest) */
    if (dest != EAX) {
        x86_64_push_reg(s, EAX);
    }
    if (dest != EDX) {
        x86_64_push_reg(s, EDX);
    }
    if ((src2 == EAX) || (src2 == EDX)) {
        /* we'll move src2 to EBP */
        tmp_src2 = EBP;
        x86_64_push_reg(s, EBP);
        x86_64_movl(s, EBP, src2);
    }
    if (src1 != EAX) {
        x86_64_movl(s, EAX, src1);
    }
    if (type == DILL_I) {
        BYTE_OUT1(s, 0x99); /* cltd  (altern  cdq) */
                            /*	rex |= REX_W;*/
    } else if (sign) {
        x86_64_rshai(s, EDX, EAX, 63);
        /*	x86_64_movl(s, EDX, EAX);*/
    } else {
        x86_64_clear(s, EDX);
    }
    if (tmp_src2 > RDI)
        rex |= REX_B;
    BYTE_OUT2R(s, rex, 0xf7, ModRM(0x3, sign ? 0x7 : 0x6, tmp_src2));
    if ((src2 == EDX) || (src2 == EAX)) {
        x86_64_pop_reg(s, EBP);
    }
    if (div && (dest != EAX)) {
        x86_64_movl(s, dest, EAX);
    }
    if (!div && (dest != EDX)) {
        x86_64_movl(s, dest, EDX);
    }
    if (dest != EDX) {
        x86_64_pop_reg(s, EDX);
    }
    if (dest != EAX) {
        x86_64_pop_reg(s, EAX);
    }
}

static int group1_eax_op[] = {
    0x05 /* add */, 0x0d /* or */,  0x15 /* adc */, 0x1d /* sbb */,
    0x25 /* and */, 0x2d /* sub */, 0x35 /* xor */, 0x3d /* cmp */
};

extern void
x86_64_arith3i(dill_stream s, int op, int typ, int dest, int src, IMM_TYPE imm)
{
    int rex = 0;
    if ((typ == DILL_L) || (typ == DILL_UL) || (typ == DILL_P)) {
        rex = REX_W;
    }

    if (dest != src) {
        int rex1 = rex;
        if (dest > RDI)
            rex1 |= REX_B;
        if (src > RDI)
            rex1 |= REX_R;
        BYTE_OUT2R(s, rex1, MOV32, ModRM(0x3, src, dest));
    }
    if ((imm <= 127) && (imm > -128)) {
        int rex1 = rex;
        if (dest > RDI)
            rex1 |= REX_B;
        BYTE_OUT3R(s, rex1, 0x83, ModRM(0x3, op, dest), imm & 0xff);
        return;
    }
    if ((0xffffffff80000000 & imm) == 0) {
        /* int-sized immediate */
        if (dest == EAX) {
            BYTE_OUT1IR(s, rex, group1_eax_op[op], (int)imm);
        } else {
            int rex2 = rex;
            if (dest > RDI)
                rex2 |= REX_B;
            BYTE_OUT2IR(s, rex2, 0x81, ModRM(0x3, op, dest), (int)imm);
        }
    } else {
        int tmp_reg = dest;
        if (src == dest) {
            if (src == EAX) {
                tmp_reg = R11;
            } else {
                tmp_reg = EAX;
            }
            x86_64_push_reg(s, tmp_reg);
        }
        x86_64_pset(s, typ, 0, tmp_reg, imm);
        switch (op) {
        case 0:        /* add */
            op = 0x03; /* reg to reg */
            break;
        case 5:        /* sub */
            op = 0x2b; /* reg to reg */
            break;
        case 4:        /* and */
            op = 0x23; /* reg to reg */
            break;
        case 1:        /* or  */
            op = 0x0b; /* reg to reg */
            break;
        case 6:        /* xor */
            op = 0x33; /* reg to reg */
            break;
        }
        x86_64_arith3(s, op, typ, dest, src, tmp_reg);
        if (src == dest) {
            x86_64_pop_reg(s, tmp_reg);
        }
    }
}

extern void
x86_64_shift(dill_stream s, int op, int type, int dest, int src1, int src2)
{
    int tmp_dest = dest;
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
    }
    if ((dest == ECX) || (dest == src2)) {
        x86_64_push_reg(s, EAX);
        tmp_dest = EAX;
    }
    if (tmp_dest != src1) {
        x86_64_movl(s, tmp_dest, src1);
    }
    if (src2 != ECX) {
        x86_64_push_reg(s, ECX);
        x86_64_movl(s, ECX, src2);
    }
    if (tmp_dest > RDI)
        rex |= REX_B;
    BYTE_OUT2R(s, rex, 0xd3, ModRM(0x3, op, tmp_dest));
    if (src2 != ECX) {
        x86_64_pop_reg(s, ECX);
    }
    if ((dest == ECX) || (dest == src2)) {
        x86_64_movl(s, dest, tmp_dest);
        x86_64_pop_reg(s, EAX);
    }
}

extern void
x86_64_shifti(dill_stream s, int op, int type, int dest, int src, IMM_TYPE imm)
{
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P)) {
        rex = REX_W;
    }
    if (dest != src) {
        x86_64_movl(s, dest, src);
    }
    if (dest > RDI)
        rex |= REX_B;
    BYTE_OUT3R(s, rex, 0xc1, ModRM(0x3, op, dest), imm & 0xff);
}

#define CONV(x, y) ((x * 100) + y)
extern void
x86_64_convert(dill_stream s, int from_type, int to_type, int dest, int src)
{
    switch (CONV(from_type, to_type)) {
    case CONV(DILL_I, DILL_L):
    case CONV(DILL_L, DILL_I):
    case CONV(DILL_I, DILL_UL):
        /* signext24 - lsh32, rsha32 */
        x86_64_lshi(s, dest, src, 32);
        x86_64_rshai(s, dest, dest, 32);
        break;
    case CONV(DILL_U, DILL_I):
    case CONV(DILL_U, DILL_UL):
    case CONV(DILL_U, DILL_L):
    case CONV(DILL_I, DILL_U):
    case CONV(DILL_UL, DILL_U):
    case CONV(DILL_L, DILL_U):
        /* clear upper 32 - lsh32, rsh32 */
        x86_64_lshi(s, dest, src, 32);
        x86_64_rshi(s, dest, dest, 32);
        break;
    case CONV(DILL_UL, DILL_I):
    case CONV(DILL_UL, DILL_L):
    case CONV(DILL_L, DILL_UL):
    case CONV(DILL_P, DILL_UL):
    case CONV(DILL_UL, DILL_P):
        if (src == dest)
            return;
        x86_64_movl(s, dest, src);
        break;
    case CONV(DILL_D, DILL_F):
    case CONV(DILL_F, DILL_D):
        /* cvtss2sd */
        /* cvtsd2ss */
        {
            int rex = 0;
            /* cvts{d,s}2s{s,d} */
            if (src > RDI)
                rex |= REX_B;
            if (dest > RDI)
                rex |= REX_R;
            BYTE_OUT1R3(s, (from_type == DILL_D) ? 0xf2 : 0xf3, rex, 0xf, 0x5a,
                        ModRM(0x3, dest, src));
            break;
        }
    case CONV(DILL_F, DILL_U):
    case CONV(DILL_F, DILL_UL):
    case CONV(DILL_F, DILL_US):
    case CONV(DILL_F, DILL_UC):
    case CONV(DILL_D, DILL_U):
    case CONV(DILL_D, DILL_UL):
    case CONV(DILL_D, DILL_US):
    case CONV(DILL_D, DILL_UC): {
        int return_reg;
        if (from_type == DILL_F) {
            x86_64_convert(s, DILL_F, DILL_D, XMM0, src);
            src = XMM0;
        }
        return_reg = dill_scallul(s, (void*)dill_x86_64_hidden_DtoUL,
                                  "dill_x86_64_hidden_DtoUL", "%d", src);
        if (to_type == DILL_U) {
            x86_64_lshi(s, dest, return_reg, 32);
            x86_64_rshi(s, dest, dest, 32);
        } else {
            x86_64_pmov(s, to_type, dest, return_reg);
        }
        break;
    }
    case CONV(DILL_F, DILL_I):
    case CONV(DILL_F, DILL_L):
    case CONV(DILL_F, DILL_S):
    case CONV(DILL_F, DILL_C):
    case CONV(DILL_D, DILL_I):
    case CONV(DILL_D, DILL_L):
    case CONV(DILL_D, DILL_S):
    case CONV(DILL_D, DILL_C): {
        int rex = 0;
        /* cvttsd2si */
        if (to_type == DILL_L)
            rex = REX_W;
        if (src > RDI)
            rex |= REX_B;
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, (from_type == DILL_D) ? 0xf2 : 0xf3, rex, 0xf, 0x2c,
                    ModRM(0x3, dest, src));
        break;
    }
    case CONV(DILL_I, DILL_D):
    case CONV(DILL_I, DILL_F):
    case CONV(DILL_S, DILL_D):
    case CONV(DILL_S, DILL_F):
    case CONV(DILL_C, DILL_D):
    case CONV(DILL_C, DILL_F): {
        int rex = 0;
        /* cvtsi2s{s,d} */
        if (src > RDI)
            rex |= REX_B;
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, (to_type == DILL_D) ? 0xf2 : 0xf3, rex, 0xf, 0x2a,
                    ModRM(0x3, dest, src));
        break;
    }
    case CONV(DILL_L, DILL_F):
    case CONV(DILL_L, DILL_D):
    case CONV(DILL_U, DILL_D):
    case CONV(DILL_U, DILL_F):
    case CONV(DILL_US, DILL_D):
    case CONV(DILL_US, DILL_F):
    case CONV(DILL_UC, DILL_D):
    case CONV(DILL_UC, DILL_F): {
        int rex = REX_W;
        /* cvtsi2s{s,d} */
        if (src > RDI)
            rex |= REX_B;
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, (to_type == DILL_D) ? 0xf2 : 0xf3, rex, 0xf, 0x2a,
                    ModRM(0x3, dest, src));
        break;
    }
    case CONV(DILL_UL, DILL_D):
    case CONV(DILL_UL, DILL_F): {
        int return_reg = dill_scalld(s, (void*)dill_x86_64_hidden_ULtoD,
                                     "dill_x86_64_hidden_ULtoD", "%l", src);
        if (to_type == DILL_F) {
            BYTE_OUT1R3(s, 0xf2, 0, 0xf, 0x5a,
                        ModRM(0x3, return_reg, return_reg));
        }
        x86_64_pmov(s, to_type, dest, return_reg);
        break;
    }
    case CONV(DILL_C, DILL_I):
    case CONV(DILL_C, DILL_L):
    case CONV(DILL_C, DILL_U):
    case CONV(DILL_C, DILL_UL):
    case CONV(DILL_C, DILL_S):
    case CONV(DILL_S, DILL_C):
    case CONV(DILL_US, DILL_C):
        /* signext56 - lsh56, rsha56 */
        x86_64_lshi(s, dest, src, 56);
        x86_64_rshai(s, dest, dest, 56);
        break;
    case CONV(DILL_C, DILL_US):
        /* signext56 - lsh56, rsha56 */
        x86_64_lshi(s, dest, src, 56);
        x86_64_rshai(s, dest, dest, 56);
        x86_64_andi(s, dest, dest, 0xffff);
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
        x86_64_andi(s, dest, src, 0xff);
        break;
    case CONV(DILL_S, DILL_I):
    case CONV(DILL_S, DILL_L):
    case CONV(DILL_S, DILL_U):
    case CONV(DILL_S, DILL_UL):
    case CONV(DILL_S, DILL_US):
    case CONV(DILL_US, DILL_S):
        /* signext48 - lsh48, rsha48 */
        x86_64_lshi(s, dest, src, 48);
        x86_64_rshai(s, dest, dest, 48);
        break;
    case CONV(DILL_US, DILL_I):
    case CONV(DILL_US, DILL_L):
    case CONV(DILL_US, DILL_U):
    case CONV(DILL_US, DILL_UL):
    case CONV(DILL_I, DILL_S):
    case CONV(DILL_U, DILL_S):
    case CONV(DILL_L, DILL_S):
    case CONV(DILL_UL, DILL_S):
    case CONV(DILL_I, DILL_US):
    case CONV(DILL_U, DILL_US):
    case CONV(DILL_L, DILL_US):
    case CONV(DILL_UL, DILL_US):
        /* zero uppper 48 - lsh48, rsh48 */
        x86_64_lshi(s, dest, src, 48);
        x86_64_rshi(s, dest, dest, 48);
        break;
    default:
        printf("Unknown case in x86_64 convert %d\n", CONV(from_type, to_type));
    }
}

static unsigned char op_conds[] = {
    0x84,
    /* dill_beq_code */ /* signed */
    0x8d,               /* dill_bge_code */
    0x8F,               /* dill_bgt_code */
    0x8e,               /* dill_ble_code */
    0x8c,               /* dill_blt_code */
    0x85,               /* dill_bne_code */

    0x84,
    /* dill_beq_code */ /* unsigned */
    0x83,               /* dill_bge_code */
    0x87,               /* dill_bgt_code */
    0x86,               /* dill_ble_code */
    0x82,               /* dill_blt_code */
    0x85,               /* dill_bne_code */
};

static unsigned char fop_conds[] = {
    0x84,
    /* dill_beq_code */ /* z = 1*/
    0x86,
    /* dill_bge_code */ /* jna */
    0x82,               /* dill_bgt_code */
    0x83,
    /* dill_ble_code */ /* c = 0 */
    0x87,               /* dill_blt_code */
    0x85,               /* dill_bne_code */
};

extern void
x86_64_branch(dill_stream s, int op, int type, int src1, int src2, int label)
{
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P))
        rex |= REX_W;

    switch (type) {
    case DILL_U:
    case DILL_UL:
    case DILL_US:
    case DILL_UC:
        op += 6; /* second set of codes */
                 /* fall through */
    default:
        if (src1 > RDI)
            rex |= REX_B;
        if (src2 > RDI)
            rex |= REX_R;
        if (type == DILL_D) {
            BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x2e, ModRM(0x3, src2, src1));
        } else if (type == DILL_F) {
            BYTE_OUT3R(s, rex, 0x0f, 0x2e, ModRM(0x3, src2, src1));
        } else {
            BYTE_OUT2R(s, rex, 0x39, ModRM(0x3, src2, src1));
        }
        dill_mark_branch_location(s, label);
        if ((type == DILL_D) || (type == DILL_F)) {
            BYTE_OUT2I(s, 0x0f, fop_conds[op], 0);
        } else {
            BYTE_OUT2I(s, 0x0f, op_conds[op], 0);
        }
    }
    x86_64_nop(s);
}

static unsigned char set_op_conds[] = {
    0x94,
    /* dill_beq_code */ /* signed */
    0x9d,               /* dill_bge_code */
    0x9F,               /* dill_bgt_code */
    0x9e,               /* dill_ble_code */
    0x9c,               /* dill_blt_code */
    0x95,               /* dill_bne_code */

    0x94,
    /* dill_beq_code */ /* unsigned */
    0x93,               /* dill_bge_code */
    0x97,               /* dill_bgt_code */
    0x96,               /* dill_ble_code */
    0x92,               /* dill_blt_code */
    0x95,               /* dill_bne_code */

    0x94,
    /* dill_beq_code */ /* floating */ /* z = 1*/
    0x96,
    /* dill_bge_code */ /* jna */
    0x92,               /* dill_bgt_code */
    0x93,
    /* dill_ble_code */ /* c = 0 */
    0x97,               /* dill_blt_code */
    0x95,               /* dill_bne_code */
};

extern void
x86_64_compare(dill_stream s, int op, int type, int dest, int src1, int src2)
{
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P))
        rex |= REX_W;

    switch (type) {
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
        op += 6; /* second set of codes */
        /* fall through */
        break;
    case DILL_F:
    case DILL_D:
        op += 12; /* third set of codes */
    }
    if (src1 > RDI)
        rex |= REX_B;
    if (src2 > RDI)
        rex |= REX_R;
    if (type == DILL_D) {
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x2e, ModRM(0x3, src2, src1));
    } else if (type == DILL_F) {
        BYTE_OUT3R(s, rex, 0x0f, 0x2e, ModRM(0x3, src2, src1));
    } else {
        BYTE_OUT2R(s, rex, 0x39, ModRM(0x3, src2, src1)); /* compare */
    }
    BYTE_OUT3(s, 0x0f, set_op_conds[op], ModRM(0x3, EAX, EAX)); /* sete */
    rex = 0;
    if (dest > RDI)
        rex = REX_R;
    BYTE_OUT3R(s, rex, 0x0f, 0xb6, ModRM(0x3, dest, EAX)); /* movzbl */
}

extern void
x86_64_comparei(dill_stream s,
                int op,
                int type,
                int dest,
                int src,
                IMM_TYPE imm)
{
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P))
        rex |= REX_W;

    switch (type) {
    case DILL_F:
    case DILL_D:
        fprintf(stderr, "Shouldn't happen\n");
        break;
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
        op += 6; /* second set of codes */
                 /* fall through */
    }
    if (src > RDI)
        rex |= REX_B;
    if (imm < 0xffffffff) {
        BYTE_OUT2IR(s, rex, 0x81, ModRM(0x3, 0x7, src), (int)imm); /* cmp */
    } else {
        x86_64_setl(s, EAX, imm);
        BYTE_OUT2R(s, rex, 0x39, ModRM(0x3, EAX, src));
    }
    BYTE_OUT3(s, 0x0f, set_op_conds[op], ModRM(0x3, EAX, EAX)); /* sete */
    BYTE_OUT3(s, 0x0f, 0xb6, ModRM(0x3, dest, EAX));            /* movzbl */
    rex = 0;
    if (dest > RDI)
        rex = REX_R;
    BYTE_OUT3R(s, rex, 0x0f, 0xb6, ModRM(0x3, dest, EAX)); /* movzbl */
}

extern void
x86_64_jump_to_label(dill_stream s, unsigned long label)
{
    dill_mark_branch_location(s, label);
    BYTE_OUT1I(s, 0xe9, 0);
}

extern void
x86_64_jump_to_reg(dill_stream s, unsigned long reg)
{
    int rex = 0;
    if (reg > RDI)
        rex |= REX_B;
    BYTE_OUT2R(s, rex, 0xff, ModRM(0x3, 0x4, reg));
}

extern void
x86_64_jump_to_imm(dill_stream s, void* imm)
{
    x86_64_seti(s, EAX, (int)(intptr_t)imm);
    BYTE_OUT2(s, 0xff, ModRM(0x3, 0x4, EAX));
}

extern void
x86_64_jal(dill_stream s, int return_addr_reg, int target)
{
    /* jump, source addr to return_addr_reg */
}

extern void
x86_64_special(dill_stream s, special_operations type, intptr_t param)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;

    switch (type) {
    case DILL_NOP:
        break;
    case DILL_SEGMENTED_FOLLOWS:
        switch (param) {
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

static void
internal_push(dill_stream s, int type, int immediate, void* value_ptr)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    struct arg_info arg;

    arg.is_immediate = immediate;
    switch (type) {
    case DILL_C:
    case DILL_S:
    case DILL_I:
    case DILL_L:
    case DILL_P:
    case DILL_EC:
        arg.type = DILL_L;
        break;
    case DILL_UC:
    case DILL_US:
    case DILL_U:
    case DILL_UL:
        arg.type = DILL_UL;
        break;
    case DILL_D:
    case DILL_F:
        arg.type = type;
        break;
    default:
        assert(0);
    }

    if (smi->varidiac_call &&
        (smi->int_arg_count + smi->float_arg_count >= smi->non_var_args)) {
        if (type == DILL_F) {
            arg.type = DILL_D;
            if (immediate) {
                /* value_ptr is already pointing to a double, no special case */
                type = DILL_D;
            }
        }
    }
    if ((arg.type != DILL_D) && (arg.type != DILL_F)) {
#ifndef USE_WINDOWS_CALLS
        if (smi->int_arg_count < sizeof(arg_regs) / sizeof(arg_regs[0])) {
#else
        if (smi->float_arg_count + smi->int_arg_count <
            sizeof(arg_regs) / sizeof(arg_regs[0])) {
#endif
            arg.is_register = 1;
            arg.in_reg = arg.out_reg = arg_regs[smi->int_arg_count];
#ifdef USE_WINDOWS_CALLS
            arg.in_reg = arg.out_reg =
                arg_regs[smi->int_arg_count + smi->float_arg_count];
            arg.offset = smi->cur_arg_offset;
            smi->cur_arg_offset +=
                roundup(type_info[(int)arg.type].size, smi->stack_align);
#endif
            smi->int_arg_count++;
        } else {
            arg.is_register = 0;
            arg.offset = smi->cur_arg_offset;
            smi->cur_arg_offset +=
                roundup(type_info[(int)arg.type].size, smi->stack_align);
        }
    } else {
#ifndef USE_WINDOWS_CALLS
        if (smi->float_arg_count < 8) {
#else
        if (smi->float_arg_count + smi->int_arg_count < 4) {
#endif
            arg.is_register = 1;
            arg.in_reg = arg.out_reg = smi->float_arg_count;
#ifdef USE_WINDOWS_CALLS
            arg.in_reg = arg.out_reg =
                smi->float_arg_count + smi->int_arg_count;
            arg.offset = smi->cur_arg_offset;
            smi->cur_arg_offset +=
                roundup(type_info[(int)arg.type].size, smi->stack_align);
#endif
            smi->float_arg_count++;
        } else {
            arg.is_register = 0;
            arg.offset = smi->cur_arg_offset;
            smi->cur_arg_offset +=
                roundup(type_info[(int)arg.type].size, smi->stack_align);
        }
    }
    if (arg.is_register == 0) {
        if (arg.offset == 0) {
            smi->call_backpatch_offset =
                (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
            dill_subli(s, ESP, ESP, 0x70909090); /* tentative for backpatch */
            smi->call_stack_space = 128;
        }
        /* store it on the stack only */
        if (arg.is_immediate) {
            int arg_type = arg.type;
            if (type == DILL_F) {
                union {
                    float f;
                    int i;
                } a;
                a.f = (float)*(double*)value_ptr;
                x86_64_setl(s, EAX, a.i);
                arg_type = DILL_I;
            } else {
                x86_64_setl(s, EAX, *(intptr_t*)value_ptr);
                arg_type = DILL_L;
            }
            x86_64_pstorei(s, arg_type, 0, EAX, ESP, arg.offset);
        } else {
            /* need to handle DILL_F upconvert to DILL_D here? */
            x86_64_pstorei(s, arg.type, 0, *(int*)value_ptr, ESP,
                           (IMM_TYPE)arg.offset);
        }
    } else {
        if ((type != DILL_F) && (type != DILL_D)) {
            if (arg.is_immediate) {
                x86_64_setl(s, arg.out_reg, *(intptr_t*)value_ptr);
            } else {
                x86_64_pmov(s, arg.type, arg.out_reg, *(int*)value_ptr);
            }
        } else {
            /* float types handled here */
            if (arg.is_immediate) {
                if ((type == DILL_F) || (type == DILL_D)) {
                    /* set appropriate register */
#ifdef USE_WINDOWS_CALLS
                    x86_64_setf(s, type, 0, arg.out_reg, *(double*)value_ptr);
                    x86_64_setl(s, arg_regs[arg.out_reg],
                                *(IMM_TYPE*)value_ptr);
#else
                    x86_64_setf(s, type, 0, arg.out_reg, *(double*)value_ptr);
#endif
                } else {
                    x86_64_setl(s, arg.out_reg, *(intptr_t*)value_ptr);
                }
            } else {
                /* move to the appropriate float reg */
                if ((type == DILL_F) && (arg.type == DILL_D)) {
                    /* special case for upconverting varidiac args */
                    x86_64_convert(s, DILL_F, DILL_D, arg.out_reg,
                                   *(int*)value_ptr);
                } else {
                    x86_64_mov(s, type, 0, arg.out_reg, *(int*)value_ptr);
                }
#ifdef USE_WINDOWS_CALLS
                {
                    int rex = REX_W;
                    if (arg_regs[arg.out_reg] > RDI)
                        rex |= REX_B;
                    if (arg.out_reg > RDI)
                        rex |= REX_R;
                    // move float reg to corresponding integer
                    BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x7e,
                                ModRM(0x3, arg.out_reg, arg_regs[arg.out_reg]));
                }
#endif
            }
        }
    }
}

static void
push_init(dill_stream s)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    smi->varidiac_call = 0;
    smi->cur_arg_offset = 0;
    smi->int_arg_count = 0;
    smi->float_arg_count = 0;
    smi->float_arg_offset = 0;
    smi->call_stack_space = 0;
}

extern void
x86_64_push(dill_stream s, int type, int reg)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    if ((type == DILL_V) && (reg <= -1)) {
        push_init(s);
        if (reg <= -2) {
            smi->varidiac_call = 1;
            smi->non_var_args = -(reg + 2);
        }
    } else {
        internal_push(s, type, 0, &reg);
    }
}

extern void
x86_64_pushi(dill_stream s, int type, IMM_TYPE value)
{
    internal_push(s, type, 1, &value);
}

extern void
x86_64_pushpi(dill_stream s, int type, void* value)
{
    internal_push(s, type, 1, &value);
}

extern void
x86_64_pushfi(dill_stream s, int type, double value)
{
    internal_push(s, type, 1, &value);
}

extern int
x86_64_calli(dill_stream s, int type, void* xfer_address, const char* name)
{
    int rex = REX_W;
    int i;
    int tmp_call_reg = R11;
    if (tmp_call_reg > RDI)
        rex |= REX_B;

    /* save temporary registers */
    for (i = XMM8; i < XMM15; i += 1) {
        if (dill_mustsave(&s->p->tmp_f, i)) {
            x86_64_save_restore_op(s, 0, DILL_D, i);
        }
    }

    /* save temporary registers */
    dill_mark_call_location(s, name, xfer_address);
    BYTE_OUT1LR(s, rex, 0xb8 + (0x7 & tmp_call_reg), 0); /* setl */
    int ret_reg = x86_64_callr(s, type, R11);

    /* restore temporary registers */
    for (i = XMM8; i < XMM15; i += 1) {
        if (dill_mustsave(&s->p->tmp_f, i)) {
            x86_64_save_restore_op(s, 1, DILL_D, i);
        }
    }
    return ret_reg;
}

extern int
x86_64_callr(dill_stream s, int type, int src)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    int caller_side_ret_reg = EAX;
    int rex = REX_W;
    if (src > RDI)
        rex |= REX_B;

    /* save temporary registers */
    /* call through reg */
    x86_64_setl(s, EAX, smi->float_arg_count);
    BYTE_OUT2R(s, rex, 0xff, ModRM(0x3, 0x2, src));
    /* restore temporary registers */
    if ((type == DILL_D) || (type == DILL_F)) {
        /*	caller_side_ret_reg = _f0;*/
    }
    if (smi->call_stack_space != 0) {
        void* save_ip = s->p->cur_ip;
        int call_stack_size = roundup(smi->call_stack_space, 8) + 16;

        /* backpatch space reservation */
        s->p->cur_ip = (char*)s->p->code_base + smi->call_backpatch_offset;
        dill_subli(s, ESP, ESP, call_stack_size);

        s->p->cur_ip = save_ip;
        /* undo arg space reservation */
        dill_addli(s, ESP, ESP, call_stack_size);
    }
    return caller_side_ret_reg;
}

extern void
x86_64_branchi(dill_stream s,
               int op,
               int type,
               int src,
               IMM_TYPE imm,
               int label)
{
    int rex = 0;
    if ((type == DILL_L) || (type == DILL_UL) || (type == DILL_P))
        rex |= REX_W;

    switch (type) {
    case DILL_F:
    case DILL_D:
        fprintf(stderr, "Shouldn't happen\n");
        break;
    case DILL_U:
    case DILL_UL:
    case DILL_US:
    case DILL_UC:
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
        if (src > RDI)
            rex |= REX_B;
        if (((uintptr_t)imm) < 0x7fffffff) {
            BYTE_OUT2IR(s, rex, 0x81, ModRM(0x3, 0x7, src), (int)imm); /* cmp */
        } else {
            x86_64_setl(s, EAX, imm);
            BYTE_OUT2R(s, rex, 0x39, ModRM(0x3, EAX, src));
        }
        dill_mark_branch_location(s, label);
        BYTE_OUT2I(s, 0x0f, op_conds[op], 0);
    }
}

extern void
x86_64_ret(dill_stream s, int data1, int data2, int src)
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
        if (src != EAX)
            x86_64_movl(s, EAX, src);
        break;
    case DILL_F:
        if (src != XMM0)
            x86_64_movf(s, XMM0, src);
        break;
    case DILL_D:
        if (src != XMM0)
            x86_64_movd(s, XMM0, src);
        break;
    }
    x86_64_proc_ret(s);
}

extern void
x86_64_retf(dill_stream s, int data1, int data2, double imm)
{
    switch (data1) {
    case DILL_F:
        x86_64_setf(s, DILL_F, 0, XMM0, imm);
        break;
    case DILL_D:
        x86_64_setf(s, DILL_D, 0, XMM0, imm);
        break;
    }
}

extern void
x86_64_reti(dill_stream s, int data1, int data2, IMM_TYPE imm)
{
    switch (data1) {
    case DILL_C:
    case DILL_UC:
    case DILL_S:
    case DILL_US:
    case DILL_I:
    case DILL_U:
        x86_64_seti(s, EAX, (int)imm);
        break;
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        x86_64_setl(s, EAX, imm);
        break;
    }
    x86_64_proc_ret(s);
}

static void
x86_64_data_link(dill_stream s)
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
x86_64_branch_link(dill_stream s)
{
    struct branch_table* t = &s->p->branch_table;
    int i;

    for (i = 0; i < t->branch_count; i++) {
        int label = t->branch_locs[i].label;
        int label_offset = t->label_locs[label] - t->branch_locs[i].loc;
        char* branch_addr =
            (char*)((char*)s->p->code_base + t->branch_locs[i].loc);
        int offset;
        if (*branch_addr == 0x0f) {
            branch_addr += 2; /* conditional */
            offset = 6;
        } else {
            branch_addr++; /* unconditional */
            offset = 5;
        }
        int tmp = label_offset - offset;
        memcpy(branch_addr, &tmp, 4);
    }
}

extern void
x86_64_rt_call_link(char* code, call_t* t);

static void
x86_64_call_link(dill_stream s)
{
    call_t* t = &s->p->call_table;

    x86_64_rt_call_link(s->p->code_base, t);
}

static void
x86_64_emit_save(dill_stream s)
{
    x86_64_mach_info smi = (x86_64_mach_info)s->p->mach_info;
    void* save_ip = s->p->cur_ip;
    int ar_size = smi->act_rec_size;
    int prefix_size;
    ar_size = roundup(ar_size, 16) + 8 + 32;

    s->p->cur_ip = (char*)s->p->code_base;

    prefix_size = generate_prefix_code(s, 0, ar_size, NULL);

    s->p->cur_ip = (char*)s->p->code_base - prefix_size + smi->backpatch_offset;
    s->p->fp = (char*)s->p->cur_ip;

    /* re-generate in the right place */
    if (prefix_size != generate_prefix_code(s, 0, ar_size, NULL)) {
        printf("2nd generation different than first\n");
    }

    s->p->cur_ip = save_ip;
}

#ifdef USE_VIRTUAL_PROTECT
#include <windows.h>
#include <intrin.h>
#include <memoryapi.h>
#endif

static void
x86_64_flush(void* base, void* limit)
{
#if defined(HOST_X86_64)
    {
        volatile void* ptr = base;

        /* flush every 8 bytes of preallocated insn stream. */
        while ((char*)ptr < (char*)limit) {
#ifndef _MSC_VER
#ifdef __x86_64__
            asm volatile("clflush (%0)" : /* */ : "r"(ptr));
#endif
#else
            _mm_clflush((const void*)ptr);
#endif
            ptr = (char*)ptr + 8;
        }
#ifndef _MSC_VER
        asm volatile("nop");
        asm volatile("nop");
        asm volatile("nop");
        asm volatile("nop");
        asm volatile("nop");
#endif
    }
#endif
#ifdef USE_VIRTUAL_PROTECT
    int result;
    DWORD dummy;
    size_t size = ((intptr_t)limit - (intptr_t)base);
    result = VirtualProtect(base, size, PAGE_EXECUTE_READWRITE, &dummy);
#endif
}
extern void
x86_64_end(dill_stream s)
{
    x86_64_proc_ret(s);
    x86_64_branch_link(s);
    x86_64_call_link(s);
    x86_64_data_link(s);
    x86_64_emit_save(s);
    x86_64_flush(s->p->code_base, s->p->code_limit);
}

extern void
x86_64_package_end(dill_stream s)
{
    x86_64_proc_ret(s);
    x86_64_branch_link(s);
    x86_64_emit_save(s);
    x86_64_flush(s->p->code_base, s->p->code_limit);
}

extern void*
x86_64_clone_code(dill_stream s, void* new_base, int available_size)
{
    int size = dill_code_size(s);
    void* old_base = s->p->code_base;
    void* native_base = s->p->code_base;
    if (available_size < size) {
        return NULL;
    }
    if (native_base == NULL)
        native_base = s->p->native.code_base;
    memcpy(new_base, native_base, size);
    s->p->code_base = new_base;
    s->p->cur_ip = (char*)new_base + size;
    s->p->fp = new_base;
    x86_64_branch_link(s);
    x86_64_call_link(s);
    x86_64_data_link(s);
    s->p->code_base = old_base;
    s->p->cur_ip = (char*)old_base + size;
    s->p->fp = old_base;
    x86_64_flush(new_base, (char*)new_base + size);

    return new_base;
}

extern void
x86_64_pset(dill_stream s, int type, int junk, int dest, IMM_TYPE imm)
{
    switch (type) {
    case DILL_L:
    case DILL_UL:
    case DILL_P:
        x86_64_setl(s, dest, imm);
        break;
    default:
        x86_64_seti(s, dest, (int)imm);
        break;
    }
    s->p->used_frame++;
}

extern void
x86_64_setf(dill_stream s, int type, int junk, int dest, double imm)
{
    union {
        float f;
        int i;
    } a;
    union {
        double d;
        int i[2];
        intptr_t l;
    } b;
    if (type == DILL_F) {
        int rex = 0;
        a.f = (float)imm;
        x86_64_seti(s, EAX, a.i);
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
    } else {
        int rex = REX_W;
        b.d = imm;
        x86_64_setl(s, EAX, b.l);
        if (dest > RDI)
            rex |= REX_R;
        BYTE_OUT1R3(s, 0x66, rex, 0x0f, 0x6e, ModRM(0x3, dest, EAX));
    }
}

#define bit_R(x) (1 << x)

extern void
x86_64_reg_init(dill_stream s)
{
    s->p->var_i.init_avail[0] =
        (bit_R(EBX) | bit_R(R12) | bit_R(R13) | bit_R(R14) | bit_R(R15));
    s->p->var_i.members[0] = s->p->var_i.init_avail[0];
    s->p->tmp_i.init_avail[0] = (bit_R(R11));
    s->p->tmp_i.members[0] = s->p->tmp_i.init_avail[0] | bit_R(EAX);
    s->p->var_f.init_avail[0] = 0;
    s->p->var_f.members[0] = s->p->var_f.init_avail[0];
    s->p->tmp_f.init_avail[0] =
        (bit_R(XMM8) | bit_R(XMM9) | bit_R(XMM10) | bit_R(XMM11) |
         bit_R(XMM12) | bit_R(XMM13) | bit_R(XMM14) | bit_R(XMM15));
    s->p->tmp_f.members[0] =
        s->p->tmp_f.init_avail[0] | bit_R(XMM0) | bit_R(XMM1) | bit_R(XMM2) |
        bit_R(XMM3) | bit_R(XMM4) | bit_R(XMM5) | bit_R(XMM6) | bit_R(XMM7);
}

extern void*
gen_x86_64_mach_info(dill_stream s)
{
    x86_64_mach_info smi = malloc(sizeof(*smi));
    if (s->p->mach_info != NULL) {
        free(s->p->mach_info);
        s->p->mach_info = NULL;
        s->p->native.mach_info = NULL;
    }
    x86_64_reg_init(s);
    smi->act_rec_size = 0;
    smi->stack_align = 8;           /* 8 for x86_64 */
    smi->stack_constant_offset = 0; /* 2047 for x86_64v9 */
    smi->conversion_word = 0;
    smi->fcu_word = 0;
    smi->save_base = 0;
    smi->backpatch_offset = 0;
    return smi;
}
#if defined(HAVE_DIS_ASM_H) && !defined(NO_DISASSEMBLER)
/* GENERIC BINUTILS DISASSEMBLER */
#include "dis-asm.h"

#define MAXLENGTH                                                              \
    (1 << 23) /* Max length of function that can be disassembled */

#ifdef LINUX_KERNEL_MODULE
extern int
kfprintf(FILE* file, const char* fmt, ...)
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

#undef stdout
#define stdout (FILE*)0
#define FPRINTF_FUNCTION kfprintf
#else
#define FPRINTF_FUNCTION fprintf
#endif
extern int
x86_64_init_disassembly_info(dill_stream s, void* ptr)
{
    struct disassemble_info* i = ptr;
#ifdef INIT_DISASSEMBLE_INFO_THREE_ARG
    INIT_DISASSEMBLE_INFO(*i, stdout, FPRINTF_FUNCTION);
    i->endian = BFD_ENDIAN_LITTLE;
#else
    INIT_DISASSEMBLE_INFO(*i, stdout);
#endif
    i->mach = bfd_mach_x86_64;
    i->disassembler_options = "x86-64";
    if (s->p->code_base != NULL) {
        i->buffer = (bfd_byte*)s->p->code_base;
        i->buffer_vma = (bfd_vma)s->p->code_base;
    } else {
        i->buffer = (bfd_byte*)s->p->native.code_base;
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
x86_64_print_insn(dill_stream s, void* info_ptr, void* insn)
{
#ifdef HAVE_PRINT_INSN_I386
    return print_insn_i386((bfd_vma)insn, info_ptr);
#else
    return 0;
#endif
}

extern void
null_func()
{
}
extern int
x86_64_count_insn(dill_stream s, int start, int end)
{
#ifdef HAVE_PRINT_INSN_I386
    struct disassemble_info i;
    int count;
    char* insn_ptr;
#ifdef INIT_DISASSEMBLE_INFO_THREE_ARG
    INIT_DISASSEMBLE_INFO(i, stdout, (fprintf_ftype)null_func);
    i.endian = BFD_ENDIAN_LITTLE;
#else
    INIT_DISASSEMBLE_INFO(i, stdout);
#endif
    i.mach = bfd_mach_i386_i386;
    if (s->p->code_base != NULL) {
        i.buffer = (bfd_byte*)s->p->code_base;
        i.buffer_vma = (bfd_vma)s->p->code_base;
    } else {
        i.buffer = (bfd_byte*)s->p->native.code_base;
        i.buffer_vma = (bfd_vma)s->p->native.code_base;
    }
    i.buffer_length = MAXLENGTH;
    count = 0;
    insn_ptr = (char*)(i.buffer + start);
    while ((intptr_t)insn_ptr < (intptr_t)i.buffer + end) {
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
x86_64_count_insn(dill_stream s, int start, int end)
{ /* no print insn, just return the buffer length */
    return end - start;
}
extern int
x86_64_init_disassembly_info(dill_stream s, void* ptr)
{
    return 1;
}
unsigned int
x86_64_disassemble(unsigned char* bytes,
                   unsigned int max,
                   int offset,
                   char* output);
extern int
x86_64_print_insn(dill_stream s, void* info_ptr, void* insn)
{
    char out[128] = "";
    int ret = x86_64_disassemble(insn, sizeof(out), 0, out);
    printf("%s", out);
    return ret;
}
#endif

extern void
x86_64_print_reg(dill_stream s, int typ, int reg)
{
    switch (typ) {
    case DILL_C:
    case DILL_UC:
        if ((reg > 0) && (reg < sizeof(char_regs) / sizeof(char_regs[0]))) {
            printf("%s", char_regs[reg]);
            return;
        }
        break;
    case DILL_S:
    case DILL_US:
        if ((reg > 0) && (reg < sizeof(short_regs) / sizeof(short_regs[0]))) {
            printf("%s", short_regs[reg]);
            return;
        }
        break;
        /*    case DILL_C: case DILL_UC: case DILL_S: case DILL_US:*/
    case DILL_I:
    case DILL_U:
    case DILL_L:
    case DILL_UL:
        if ((reg > 0) && (reg < sizeof(int_regs) / sizeof(int_regs[0]))) {
            printf("%s", int_regs[reg]);
            return;
        }
        break;
    case DILL_F:
    case DILL_D:
        printf("Fstack");
        return;
    }
    printf("NoReg(%d)", reg);
}
