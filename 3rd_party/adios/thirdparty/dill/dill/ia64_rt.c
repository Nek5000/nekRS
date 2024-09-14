#include "config.h"
#ifdef USE_MMAP_CODE_SEG
#include <sys/mman.h>
#endif
#include "string.h"
#include "dill.h"
#include "dill_internal.h"
#include "x86.h"

extern long dill_ia64_hidden_mod(long a, long b)
{ return a % b; }
extern long dill_ia64_hidden_umod(unsigned long a, unsigned long b)
{ return a % b; }
extern long dill_ia64_hidden_div(long a, long b)
{ return a / b; }
extern long dill_ia64_hidden_udiv(unsigned long a, unsigned long b)
{ return a / b; }
extern double dill_ia64_hidden_fdiv(double a, double b)
{ return a / b; }

static xfer_entry ia64_xfer_recs[] = {
    {"dill_ia64_hidden_mod", dill_ia64_hidden_mod},
    {"dill_ia64_hidden_umod", dill_ia64_hidden_umod},
    {"dill_ia64_hidden_div", dill_ia64_hidden_div},
    {"dill_ia64_hidden_udiv", dill_ia64_hidden_udiv},
    {"dill_ia64_hidden_fdiv", dill_ia64_hidden_fdiv},
    {(char*)0, (void*)0}};

extern void
ia64_rt_call_link(char *code, call_t *t)
{
    int i;

    for(i=0; i< t->call_count; i++) {
	int *call_addr = (int*) (code + t->call_locs[i].loc + 1);
	int call_offset = (unsigned long)t->call_locs[i].xfer_addr - 
	    (unsigned long)((char*)call_addr + 4);  /* add len of call insn */

	*call_addr = call_offset;
    }
}

extern char *
ia64_package_stitch(char *code, call_t *t, dill_pkg pkg)
{
    char *tmp = code;
    dill_lookup_xfer_addrs(t, &ia64_xfer_recs[0]);
    ia64_rt_call_link(code, t);
#ifdef USE_MMAP_CODE_SEG
    tmp = (void*)mmap(0, pkg->code_size,
		      PROT_EXEC | PROT_READ | PROT_WRITE, 
		      MAP_ANONYMOUS|MAP_PRIVATE, -1, 0);
    memcpy(tmp, code, pkg->code_size);
#endif
    *((void**)tmp) = (char*)tmp + 16;
    *((void**)(tmp+8)) = 0;
    return tmp;
}

