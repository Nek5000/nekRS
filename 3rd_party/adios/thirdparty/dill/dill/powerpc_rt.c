#include "config.h"
#include "dill.h"
#include "dill_internal.h"
#include "powerpc.h"

extern long dill_powerpc_hidden_modi(int a, int b)
{ return a % b; }
extern long dill_powerpc_hidden_mod(long a, long b)
{ return a % b; }
extern unsigned long dill_powerpc_hidden_umod(unsigned long a, unsigned long b)
{ return a % b; }
extern unsigned int dill_powerpc_hidden_umodi(unsigned int a, unsigned int b)
{ return a % b; }
extern double dill_powerpc_hidden_ultod(unsigned long a)
{ return (double) a; }
extern float dill_powerpc_hidden_ultof(unsigned long a)
{ return (float) a; }
extern unsigned long dill_powerpc_hidden_dtoul(double a)
{ return (unsigned long) a; }
extern unsigned int dill_powerpc_hidden_dtou(double a)
{ return (unsigned int) a; }
extern unsigned long dill_powerpc_hidden_ftoul(float a)
{ return (unsigned long) a; }
extern unsigned int dill_powerpc_hidden_ftou(float a)
{ return (unsigned int) a; }
extern long dill_powerpc_hidden_udiv(unsigned long a, unsigned long b)
{ return a / b; }

static xfer_entry powerpc_xfer_recs[] = {
    {"dill_powerpc_hidden_modi", (void*)dill_powerpc_hidden_modi},
    {"dill_powerpc_hidden_mod", (void*)dill_powerpc_hidden_mod},
    {"dill_powerpc_hidden_umodi", (void*)dill_powerpc_hidden_umodi},
    {"dill_powerpc_hidden_umod", (void*)dill_powerpc_hidden_umod},
    {"dill_powerpc_hidden_ultod", (void*)dill_powerpc_hidden_ultod},
    {"dill_powerpc_hidden_ultof", (void*)dill_powerpc_hidden_ultof},
    {"dill_powerpc_hidden_dtoul", (void*)dill_powerpc_hidden_dtoul},
    {"dill_powerpc_hidden_ftoul", (void*)dill_powerpc_hidden_ftoul},
    {"dill_powerpc_hidden_dtou", (void*)dill_powerpc_hidden_dtou},
    {"dill_powerpc_hidden_ftou", (void*)dill_powerpc_hidden_ftou},
    {"dill_powerpc_hidden_udiv", (void*)dill_powerpc_hidden_udiv},
    {(char*)0, (void*)0}};

#define powerpc_sethi(r, imm) (HDR(0)|RD(r)|OP2(0x4)|imm)
#define powerpc_ori(dest, src, imm) (HDR(0x2)|OP3(0x2)|IM|RD(dest)|RS1(src)|imm)
#define powerpc_or(dest, src1, src2) (HDR(0x2)|OP3(0x2)|RD(dest)|RS1(src1)|RS2(src2))
#define powerpc_sllx(dest, src1,imm) (HDR(0x2)|OP3(0x25)|RD(dest)|RS1(src1)|IM|SIMM13(imm)|(1<<12) )

extern void
powerpc_rt_call_link(char *code, call_t *t, int force_plt)
{
    int i;

    for(i=0; i< t->call_count; i++) {
	int *call_addr = (int*) (code + t->call_locs[i].loc);
	if (!force_plt && (t->call_locs[i].mach_info == (void*)0)) {
	    /* no PLT */
	    int call_offset = (unsigned long)t->call_locs[i].xfer_addr -
		(unsigned long)call_addr;
	
	    /* div addr diff by 4 for powerpc offset value */
	    call_offset = call_offset >> 2;
	    *call_addr &= 0xc0000000;
	    *call_addr |= (call_offset & 0x3fffffff);
	} else {
	    int *plt = (int*)(code + t->call_locs[i].loc);
	    union {
		void *xa;
		long xl;
	    }u;
	    u.xa = t->call_locs[i].xfer_addr;
	    plt[0] = powerpc_sethi(_g1, hh(u.xl));
	    plt[1] = powerpc_ori(_g1, _g1, hm(u.xl));
	    plt[2] = powerpc_sethi(_l0, lm(u.xl));
	    plt[3] = powerpc_ori(_l0, _l0, lo(u.xl));
	    plt[4] = powerpc_sllx(_g1, _g1, 32);
	    plt[5] = powerpc_or(_g1, _g1, _l0);
	}
    }
}

static void
powerpc_flush(void *base, void *limit)
{
#if defined(HOST_POWERPC) || defined(HOST_POWERPCV9)
    {
	volatile void *ptr = base;

#ifdef __GNUC__
	/* flush every 8 bytes of preallocated insn stream. */
	while((char*)ptr < (char*) limit) {
	  //	    asm volatile ("iflush %0" : /* */ : "r" (ptr));
	    ptr = (char *)ptr + 8;
	}
	asm volatile("nop");
	asm volatile("nop");
	asm volatile("nop");
	asm volatile("nop");
	asm volatile("nop");
#else
	int nbytes = (char*)limit - (char*)base;
	for(; nbytes > 0;nbytes -= 8) {
	    asm("add %i0, 8, %i0");
	    asm ("iflush %i0");
	}

	asm ("nop");
	asm ("nop");
	asm ("nop");
	asm ("nop");
	asm ("nop");
#endif
#ifdef USE_MEMBAR
	asm("membar #Sync");
#endif
    }
#endif
}    

extern char *
powerpc_package_stitch(char *code, call_t *t, dill_pkg pkg)
{
    int force_plt = 0;
    dill_lookup_xfer_addrs(t, &powerpc_xfer_recs[0]);
#if defined(HOST_POWERPCV9)
    force_plt = 1;
#endif
    powerpc_rt_call_link(code, t, force_plt);
    powerpc_flush(code, code + pkg->code_size);
    return code + pkg->entry_offset;
}
