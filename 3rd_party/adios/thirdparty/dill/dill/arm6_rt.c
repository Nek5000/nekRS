#include "config.h"
#include <string.h>
#include "dill.h"
#include "dill_internal.h"
#include "sys/mman.h"
#ifdef HAVE_MEMORY_H
#include "memory.h"
#endif
#include "arm6.h"
#include <string.h>

extern int arm6_hidden_modi(int a, int b)
{ return a % b; }
extern long arm6_hidden_mod(long a, long b)
{ return a % b; }
extern unsigned long arm6_hidden_umod(unsigned long a, unsigned long b)
{ return a % b; }
extern unsigned int arm6_hidden_umodi(unsigned int a, unsigned int b)
{ return a % b; }
extern double arm6_hidden_ultod(unsigned long a)
{ return (double) a; }
extern float arm6_hidden_ultof(unsigned long a)
{ return (float) a; }
extern unsigned long arm6_hidden_dtoul(double a)
{ return (unsigned long) a; }
extern unsigned int arm6_hidden_dtou(double a)
{ return (unsigned int) a; }
extern unsigned long arm6_hidden_ftoul(float a)
{ return (unsigned long) a; }
extern unsigned int arm6_hidden_ftou(float a)
{ return (unsigned int) a; }
extern unsigned long arm6_hidden_udiv(unsigned long a, unsigned long b)
{ return a / b; }
extern long arm6_hidden_div(long a, long b)
{ return a / b; }

#define COND(x)	((unsigned)((x)&0xf) << 28)
#define CLASS(x)	(((x)&0x7) << 25)

static xfer_entry arm6_xfer_recs[] = {
    {"arm6_hidden_modi", arm6_hidden_modi},
    {"arm6_hidden_mod", arm6_hidden_mod},
    {"arm6_hidden_umod", arm6_hidden_umod},
    {"arm6_hidden_umodi", arm6_hidden_umodi},
    {"arm6_hidden_ultod", arm6_hidden_ultod},
    {"arm6_hidden_ultof", arm6_hidden_ultof},
    {"arm6_hidden_dtoul", arm6_hidden_dtoul},
    {"arm6_hidden_dtou", arm6_hidden_dtou},
    {"arm6_hidden_ftoul", arm6_hidden_ftoul},
    {"arm6_hidden_ftou", arm6_hidden_ftou},
    {"arm6_hidden_udiv", arm6_hidden_udiv},
    {"arm6_hidden_div", arm6_hidden_div},
    {(char*)0, (void*)0}};

static void
arm6_rt_set_PLT_locs(call_t* t, dill_pkg pkg)
{
    /* 
     * Must set mach_info with PLT entry location.  
     * One entry per call, stacked at the end of the code block 
     */
    int i;
    int PLT_offset = pkg->code_size - 3 * 4;  /* 3 insn per PLT entry */
    for (i=t->call_count-1; i>=0; i--) {
        t->call_locs[i].mach_info = (void*) (intptr_t)PLT_offset;
	PLT_offset -= 3 * 4;
    }
}

extern void
arm6_rt_call_link(char *code, call_t *t)
{
    int i;
    
    for(i=0; i< t->call_count; i++) {
	int *call_addr = (int*) ((unsigned long)code + 
				 t->call_locs[i].loc);
	if (t->call_locs[i].mach_info == NULL) {
	    /* no PLT */
	    int call_offset = (unsigned long)t->call_locs[i].xfer_addr -
		(unsigned long)call_addr;
	    int thumb_target = (((unsigned long)t->call_locs[i].xfer_addr & 0x1) == 0x1);
	    int bit1;
	    /* compensate for arm PC lookahead */
	    call_offset = call_offset - 8;
	    /* div addr diff by 4 for arm offset value */
	    bit1 = (call_offset & 0x2) >> 1;
	    call_offset = call_offset >> 2;
	    *call_addr &= 0xff000000;
	    *call_addr |= (call_offset & 0xffffff);
	    if (thumb_target) {
	      *call_addr &= 0x00ffffff; /*  kill top bit */
	      *call_addr |= (COND(0xf)|CLASS(5)|bit1<<24);  /* blx */
	    }
	} else {
	    /* call through PLT */
	    unsigned long PLT_addr = (unsigned long)code + 
				      (unsigned long)t->call_locs[i].mach_info;
	    int call_offset = PLT_addr - (unsigned long)call_addr;

	    /* compensate for arm PC lookahead */
	    call_offset = call_offset - 8;
	    call_offset = call_offset >> 2;
	    *call_addr &= 0xff000000;
	    *call_addr |= (call_offset & 0x00ffffff);
	    PLT_addr += 8;
	    *(unsigned long*)PLT_addr = (unsigned long)t->call_locs[i].xfer_addr;
	}
    }
}

#ifndef CLEAR_CACHE_DEFINED
extern void __clear_cache(void *, void *);
#endif

static void
arm6_flush(void *base, void *limit)
{
#if defined(HOST_ARM6) || defined(HOST_ARM7)
    __clear_cache(base, limit);
#endif
}    

extern char *
arm6_package_stitch(char *code, call_t *t, dill_pkg pkg)
{
    char *tmp = code;
    dill_lookup_xfer_addrs(t, &arm6_xfer_recs[0]);
    arm6_rt_set_PLT_locs(t, pkg);
#ifdef USE_MMAP_CODE_SEG
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif
    tmp = (void*)mmap(0, pkg->code_size,
		      PROT_EXEC | PROT_READ | PROT_WRITE, 
		      MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
    memcpy(tmp, code, pkg->code_size);
#endif
    arm6_rt_call_link(tmp, t);
    arm6_flush(code, tmp+pkg->code_size);
    return tmp + pkg->entry_offset;
}

