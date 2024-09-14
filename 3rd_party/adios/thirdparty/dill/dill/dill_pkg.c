#include "config.h"

#ifdef STDC_HEADERS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#include "dill.h"
#include "dill_internal.h"

static void
unpack_package(char* package, call_t* t, char** code_p)
{
    int count;
    int pkg_size = sizeof(struct dill_pkg_1);
    char* p;
    struct dill_pkg_1* pkg = (struct dill_pkg_1*)package;

    pkg_size = (pkg_size + 7) & -8; /* round up to mod 8 */
    p = package + pkg_size;
    if (pkg->magic != 0xbeef)
        printf("Not valid package\n");
    if (pkg->pkg_version != 1)
        printf("Not valid package version\n");
    t->call_alloc = t->call_count = pkg->symbol_count;
    t->call_locs = malloc(sizeof(t->call_locs[0]) * pkg->symbol_count);
    memset(t->call_locs, 0, sizeof(t->call_locs[0]) * pkg->symbol_count);
    for (count = 0; count < pkg->symbol_count; count++) {
        size_t call_len;
        t->call_locs[count].loc = *((int*)p);
        t->call_locs[count].xfer_name = (p + sizeof(int));
        call_len = sizeof(int) + strlen(t->call_locs[count].xfer_name) + 1;
        call_len = (call_len + 7) & (size_t)-8; /* round up to mod 8 */
        p += call_len;
    }
    *code_p = p;
}

extern void*
dill_package_entry(char* package)
{
    struct dill_pkg_1* pkg = (struct dill_pkg_1*)package;
    return package + pkg->code_offset + pkg->entry_offset;
}

extern char*
sparc_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
x86_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
x86_64_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
arm5_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
arm6_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
arm8_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
ia64_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
powerpc_package_stitch(char* code, call_t* t, dill_pkg pkg);
extern char*
ppc64le_package_stitch(char* code, call_t* t, dill_pkg pkg);

extern void
dill_lookup_xfer_addrs(call_t* t, xfer_entry* x)
{
    int i;
    for (i = 0; i < t->call_count; i++) {
        xfer_entry* e = x;
        while (e->xfer_name != NULL) {
            if (strcmp(e->xfer_name, t->call_locs[i].xfer_name) == 0) {
                t->call_locs[i].xfer_addr = e->xfer_addr;
            }
            e++;
        }
    }
}

extern dill_exec_handle
dill_package_stitch(char* pkg, dill_extern_entry* extra_externs)
{
    dill_exec_handle handle = malloc(sizeof(*handle));
    char* code;
    call_t t;
    unpack_package(pkg, &t, &code);
    if (extra_externs) {
        dill_lookup_xfer_addrs(&t, (xfer_entry*)extra_externs);
    }
#if defined(HOST_X86) && !defined(DILL_IGNORE_NATIVE)
    char* p = x86_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_SPARC) && !defined(DILL_IGNORE_NATIVE)
    char* p = sparc_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_SPARCV9) && !defined(DILL_IGNORE_NATIVE)
    char* p = sparc_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_X86_64) && !defined(DILL_IGNORE_NATIVE)
    char* p = x86_64_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_IA64) && !defined(DILL_IGNORE_NATIVE)
    char* p = ia64_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_ARM5) && !defined(DILL_IGNORE_NATIVE)
    char* p = arm5_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_ARM6) && !defined(DILL_IGNORE_NATIVE)
    char* p = arm6_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_ARM7) && !defined(DILL_IGNORE_NATIVE)
    char* p = arm6_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_ARM8) && !defined(DILL_IGNORE_NATIVE)
    char* p = arm8_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_POWERPC) && !defined(DILL_IGNORE_NATIVE)
    char* p = powerpc_package_stitch(code, &t, (dill_pkg)pkg);
#elif defined(HOST_PPC64LE) && !defined(DILL_IGNORE_NATIVE)
    char* p = ppc64le_package_stitch(code, &t, (dill_pkg)pkg);
#else
    char* p = NULL;
#endif
    free(t.call_locs);
    handle->fp = (void (*)())p;
    return handle;
}
