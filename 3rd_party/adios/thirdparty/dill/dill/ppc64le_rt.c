#include "config.h"
#include "dill.h"
#include "dill_internal.h"
#include "sys/mman.h"
#ifdef HAVE_MEMORY_H
#include "memory.h"
#endif
#include "ppc64le.h"
#include <string.h>
#include <stdio.h>

#define ppc64le_ori(s, dest, src, imm) 	INSN_OUT(s, D_FORM(24, dest, src, lo16(imm)));
#define ppc64le_andi(s, dest, src, imm) INSN_OUT(s, D_FORM(28, dest, src, lo16(imm)));
#define ppc64le_or(s, dest, src1, src2) INSN_OUT(s, X_FORM(31,src1, dest, src2, 444))
#define ppc64le_lshi(s, dest, src,imm) INSN_OUT(s, MD_FORM(30,dest,src,imm & 0x1f,63-imm, 0, imm>>5));

#define lo16(im) (((long)im) & 0xffff)
#define hi16(im) ((((long)im) & 0xffffffff) >> 16)
#define lo32(im) (((long)im) & 0xfffffffff)
#define hi32(im) ((((long)im) >> 32) &  0xfffffffff)

extern unsigned long dill_ppc64le_hidden_ftoul(float a)
{ return (unsigned long) a; }
extern unsigned int dill_ppc64le_hidden_ftou(float a)
{ long l = (long) a;  unsigned int ret = l & 0xffffffff; return ret; }
extern unsigned int dill_ppc64le_hidden_dtou(double a)
{ long l = (long) a;  unsigned int ret = l & 0xffffffff; return ret; }

static xfer_entry ppc64le_xfer_recs[] = {
    {"dill_ppc64le_hidden_ftoul", (void*)dill_ppc64le_hidden_ftoul},
    {"dill_ppc64le_hidden_ftou", (void*)dill_ppc64le_hidden_ftou},
    {"dill_ppc64le_hidden_dtou", (void*)dill_ppc64le_hidden_dtou},
    {(char*)0, (void*)0}};


extern void
ppc64le_rt_call_link(char *code, call_t *t, int force_plt)
{
    int i;

    for(i=0; i< t->call_count; i++) {
	int *call_addr = (int*) (code + t->call_locs[i].loc);
	long xfer_addr = (unsigned long)t->call_locs[i].xfer_addr;
	call_addr[0] = D_FORM(15, _gpr12, _gpr0, hi16(hi32(xfer_addr)));
	call_addr[1] = D_FORM(24, _gpr12, _gpr12, lo16(hi32(xfer_addr)));
	call_addr[2] = MD_FORM(30, _gpr12, _gpr12, (32&0x1f), ((((63-32)&0x1f) << 1) | ((63-32)>>5)), 1, (32>>5));
	call_addr[3] = D_FORM(25, _gpr12, _gpr12, hi16(lo32(xfer_addr)));
	call_addr[4] = D_FORM(24, _gpr12, _gpr12, lo16(lo32(xfer_addr)));
    }
}

static void
ppc64le_flush(void *base, void *limit)
{
#if defined(HOST_PPC64LE)
    {
	volatile void *ptr = base;

#ifdef __GNUC__
	/* flush every 8 bytes of preallocated insn stream. */
	while((char*)ptr < (char*) limit) {
	    asm volatile ("dcbst 0, %0" : /* */ : "r" (ptr));
	    ptr = (char *)ptr + 128;
	}
	asm volatile("sync");
	while((char*)ptr < (char*) limit) {
	    asm volatile ("icbi 0, %0" : /* */ : "r" (ptr));
	    ptr = (char *)ptr + 128;
	}
	asm volatile("isync");
#else
	while((char*)ptr < (char*) limit) {
	    asm("dcbst 0, %r0");
	    ptr = (char *)ptr + 128;
	}
	asm ("sync");
	while((char*)ptr < (char*) limit) {
	    asm("icbi 0, %r0");
	    ptr = (char *)ptr + 128;
	}
	asm ("isync");
#endif
    }
#endif
}    

extern char *
ppc64le_package_stitch(char *code, call_t *t, dill_pkg pkg)
{
    char *tmp = code;
    dill_lookup_xfer_addrs(t, &ppc64le_xfer_recs[0]);
#ifdef USE_MMAP_CODE_SEG
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif
    tmp = (void*)mmap(0, pkg->code_size,
		      PROT_EXEC | PROT_READ | PROT_WRITE, 
		      MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
    memcpy(tmp, code, pkg->code_size);
#endif
    ppc64le_rt_call_link(tmp, t, 0);
    ppc64le_flush(tmp, tmp + pkg->code_size);
    return tmp + pkg->entry_offset;
}
