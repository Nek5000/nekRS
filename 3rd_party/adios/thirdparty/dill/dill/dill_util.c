#include "config.h"

#ifndef NATIVE_ARCH
#define NATIVE_ARCH "UNKNOWN"
#endif

#include "dill.h"

#include <ctype.h>
#include <stdio.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdint.h>
#include <string.h>
#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <unistd.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#ifdef USE_MMAP_CODE_SEG
#include <sys/mman.h>
#endif
#include "dill_internal.h"

static char* DILL_version = "DILL Version " DILL_VERSION "\n";

void
DILLprint_version()
{
    printf("%s", DILL_version);
}

extern void
init_code_block(dill_stream s);
static void
free_code_blocks(dill_stream s);

#define END_OF_CODE_BUFFER 60

static void
reset_regset(reg_set* regs)
{
    *(regs->avail) = *(regs->init_avail);
    memset(&regs->used, 0, sizeof(regs->used));
    memset(&regs->mustsave, 0, sizeof(regs->mustsave));
}

static void
dill_register_init(dill_stream s)
{
    reset_regset(&s->p->var_i);
    reset_regset(&s->p->tmp_i);
    reset_regset(&s->p->var_f);
    reset_regset(&s->p->tmp_f);
}

static void
dill_branch_init(dill_stream s)
{
    struct branch_table* t = &s->p->branch_table;
    int i;

    t->next_label = 0;
    for (i = 0; i < t->max_alloc; i++) {
        t->label_locs[i] = -1;
        t->label_name[i] = NULL;
    }
    t->branch_count = 0;
    t->data_segment_size = 0;
}

static void
dill_call_init(dill_stream s)
{
    struct call_table* t = &s->p->call_table;

    t->call_count = 0;
}

static void
dill_ret_init(dill_stream s)
{
    struct ret_table* t = &s->p->ret_table;

    t->ret_count = 0;
}
#ifdef BUILD_EMULATOR
#include "ffi.h"
void
free_emulator_handler_bits(dill_exec_handle handle);
#endif

static void
reset_context(dill_stream s)
{
    s->p->mach_reset(s);
    s->p->cur_ip = s->p->code_base;
    dill_register_init(s);
    dill_branch_init(s);
    dill_call_init(s);
    dill_ret_init(s);
    if (s->p->emu_args)
        free(s->p->emu_args);
    if (s->p->cifp)
        free(s->p->cifp);
#ifdef BUILD_EMULATOR
    if (s->p->closure)
        ffi_closure_free(s->p->closure);
#endif
    s->p->emu_args = NULL;
    s->p->cifp = NULL;
    s->p->closure = NULL;
}

extern void
dill_sparc_init(dill_stream s);
extern void
dill_sparcv9_init(dill_stream s);
extern void
dill_x86_init(dill_stream s);
extern void
dill_x86_64_init(dill_stream s);
extern void
dill_arm5_init(dill_stream s);
extern void
dill_arm6_init(dill_stream s);
extern void
dill_arm7_init(dill_stream s);
extern void
dill_arm8_init(dill_stream s);
extern void
dill_ppc64le_init(dill_stream s);
extern void
dill_ia64_init(dill_stream s);
#if defined(EMULATION_ONLY)
static void
null_init(dill_stream s)
{
    extern int virtual_type_align[];
    extern int virtual_type_size[];
    s->j = malloc(sizeof(*s->j));
    memset(s->j, 0, sizeof(*s->j));
    s->j->type_size = virtual_type_size;
    s->j->type_align = virtual_type_align;
}
#endif

static int
set_mach_reset(dill_stream s, char* arch)
{
#if (defined(MULTI_TARGET) || defined(HOST_SPARC) || defined(HOST_SPARCV9)) && \
    !defined(DILL_IGNORE_NATIVE)
    if (strcmp(arch, "sparc") == 0) {
        s->p->mach_reset = dill_sparc_init;
        return 1;
    } else if (strcmp(arch, "sparcv9") == 0) {
        s->p->mach_reset = dill_sparcv9_init;
        return 1;
    } else
#endif
#if (defined(MULTI_TARGET) || defined(HOST_X86)) && !defined(DILL_IGNORE_NATIVE)
        if (strcmp(arch, "x86") == 0) {
        s->p->mach_reset = dill_x86_init;
        return 1;
    } else
#endif
#if (defined(MULTI_TARGET) || defined(HOST_X86_64)) &&                         \
    !defined(DILL_IGNORE_NATIVE)
        if (strcmp(arch, "x86_64") == 0) {
        s->p->mach_reset = dill_x86_64_init;
        return 1;
    } else
#endif
#if (defined(MULTI_TARGET) || defined(HOST_IA64)) &&                           \
    !defined(DILL_IGNORE_NATIVE)
        if (strcmp(arch, "ia64") == 0) {
        s->p->mach_reset = dill_ia64_init;
        return 1;
    } else
#endif
#if (defined(MULTI_TARGET) || defined(HOST_ARM5)) &&                           \
    !defined(DILL_IGNORE_NATIVE)
        if (strcmp(arch, "arm5") == 0) {
        s->p->mach_reset = dill_arm5_init;
        return 1;
    }
#endif
#if (defined(MULTI_TARGET) || defined(HOST_ARM6) || defined(HOST_ARM7)) &&     \
    !defined(DILL_IGNORE_NATIVE)
    if (strcmp(arch, "arm6") == 0) {
        s->p->mach_reset = dill_arm6_init;
        return 1;
    }
#endif
#if (defined(MULTI_TARGET) || defined(HOST_ARM7)) &&                           \
    !defined(DILL_IGNORE_NATIVE)
    if (strcmp(arch, "arm7") == 0) {
        s->p->mach_reset = dill_arm6_init;
        return 1;
    }
#endif
#if (defined(MULTI_TARGET) || defined(HOST_ARM8)) &&                           \
    !defined(DILL_IGNORE_NATIVE)
    if (strcmp(arch, "arm8") == 0) {
        s->p->mach_reset = dill_arm8_init;
        return 1;
    }
#endif
#if (defined(MULTI_TARGET) || defined(HOST_PPC64LE)) &&                        \
    !defined(DILL_IGNORE_NATIVE)
    if (strcmp(arch, "ppc64le") == 0) {
        s->p->mach_reset = dill_ppc64le_init;
        return 1;
    }
#endif
#if defined(EMULATION_ONLY) || defined(DILL_IGNORE_NATIVE)
    s->p->mach_reset = null_init;
    return 1;
#endif
    return 0;
}

extern void
dill_free_stream(dill_stream s)
{
    if (s->p->branch_table.label_locs)
        free(s->p->branch_table.label_locs);
    if (s->p->branch_table.label_name) {
        int i = 0;
        for (i = 0; i < s->p->branch_table.max_alloc; i++) {
            if (s->p->branch_table.label_name[i])
                free(s->p->branch_table.label_name[i]);
        }
        free(s->p->branch_table.label_name);
    }
    if (s->p->branch_table.branch_locs)
        free(s->p->branch_table.branch_locs);
    if (s->p->branch_table.data_segment)
        free(s->p->branch_table.data_segment);
    free(s->p->call_table.call_locs);
    free(s->p->ret_table.ret_locs);
    free(s->p->c_param_regs);
    free(s->p->c_param_args);
    free(s->p->c_param_structs);
    free_code_blocks(s);
    if (s->p->mach_info) {
        if ((s->p->mach_info != s->p->virtual.mach_info) &&
            (s->p->mach_info != s->p->native.mach_info))
            free(s->p->mach_info);
    }
    if (s->p->vregs)
        free(s->p->vregs);
    if (s->p->virtual.mach_info)
        free(s->p->virtual.mach_info);
    if (s->p->native.mach_info)
        free(s->p->native.mach_info);
    if (s->p->emu_args)
        free(s->p->emu_args);
    if (s->p->cifp)
        free(s->p->cifp);
#ifdef BUILD_EMULATOR
    if (s->p->closure)
        ffi_closure_free(s->p->closure);
#endif
    free(s->p);
    s->p = NULL;
    free(s);
}

extern void
DILLprint_version();

extern dill_stream
dill_cross_init(char* arch)
{
    dill_stream s = (dill_stream)malloc(sizeof(struct dill_stream_s));
    char* env = getenv("DILL_DEBUG");
    struct branch_table* bt;
    struct call_table* ct;
    struct ret_table* rt;
    memset(s, 0, sizeof(struct dill_stream_s));
    s->p = (private_ctx)malloc(sizeof(struct dill_private_ctx));
    memset(s->p, 0, sizeof(struct dill_private_ctx));
    if (env == NULL) {
        s->dill_debug = 0;
    } else {
        s->dill_debug = 1;
        DILLprint_version();
    }
    s->p->mach_info = NULL;
    if (!set_mach_reset(s, arch)) {
        fprintf(stderr, "DILL support for architecture %s not found.\n", arch);
        free(s->p);
        free(s);
        return NULL;
    }
    init_code_block(s);
    s->p->cur_ip = s->p->code_base;
    bt = &s->p->branch_table;
    bt->max_alloc = 1;
    bt->label_locs = malloc(sizeof(bt->label_locs[0]));
    bt->label_name = malloc(sizeof(bt->label_name[0]));
    bt->branch_alloc = 1;
    bt->branch_locs = malloc(sizeof(bt->branch_locs[0]));
    bt->data_segment_size = 0;
    bt->data_segment = malloc(1);
    ct = &s->p->call_table;
    ct->call_alloc = 1;
    ct->call_count = 0;
    ct->call_locs = malloc(sizeof(ct->call_locs[0]));
    rt = &s->p->ret_table;
    rt->ret_alloc = 1;
    rt->ret_count = 0;
    rt->ret_locs = malloc(sizeof(rt->ret_locs[0]));
    s->p->c_param_count = 0;
    s->p->c_param_regs = NULL;
    s->p->c_param_args = NULL;
    s->p->c_param_structs = NULL;
    s->p->vreg_count = 0;
    s->p->vregs = malloc(1);
    reset_context(s);
    s->p->native.mach_jump = s->j;
    s->p->native.mach_reset = s->p->mach_reset;
    s->p->native.mach_info = NULL;
    s->p->native.code_base = NULL;
    s->p->native.cur_ip = s->p->cur_ip;
    s->p->native.code_limit = s->p->code_limit;
    s->p->virtual.mach_jump = NULL;
    s->p->virtual.mach_reset = NULL;
    s->p->virtual.mach_info = NULL;
    s->p->virtual.code_base = NULL;
    s->p->virtual.cur_ip = NULL;
    s->p->virtual.code_limit = NULL;
    s->p->unavail_called = 0;
    return s;
}

extern void
dill_virtual_init(dill_stream s);

extern dill_stream
dill_create_stream()
{
    dill_stream s;
    s = dill_cross_init(NATIVE_ARCH);
    s->p->native_mach_reset = s->p->mach_reset;
    s->p->native.code_base = s->p->code_base;
    s->p->native.mach_info = s->p->mach_info;
    s->p->mach_reset = dill_virtual_init;
    init_code_block(s);
    s->p->cur_ip = s->p->code_base;
    s->p->mach_reset = dill_virtual_init;
    s->p->mach_info = NULL;

    s->p->virtual.mach_jump = s->j;
    s->p->virtual.mach_reset = s->p->mach_reset;
    s->p->virtual.mach_info = s->p->mach_info;
    s->p->virtual.code_base = s->p->code_base;
    s->p->virtual.cur_ip = s->p->cur_ip;
    s->p->virtual.code_limit = s->p->code_limit;
    return s;
}

#ifndef EMULATION_ONLY
extern void
dill_native_dcg()
{
}
#endif

extern dill_stream
dill_create_raw_stream()
{
    return dill_cross_init(NATIVE_ARCH);
}

static void
extend_params(dill_stream s, int argno)
{
    int i;
    if (s->p->c_param_count == 0) {
        s->p->c_param_regs =
            malloc(sizeof(s->p->c_param_regs[0]) * (argno + 2));
        s->p->c_param_args =
            malloc(sizeof(s->p->c_param_args[0]) * (argno + 2));
        s->p->c_param_structs =
            malloc(sizeof(s->p->c_param_structs[0]) * (argno + 2));
    } else if (s->p->c_param_count <= (argno + 1)) {
        s->p->c_param_regs = realloc(
            s->p->c_param_regs, sizeof(s->p->c_param_regs[0]) * (argno + 2));
        s->p->c_param_args = realloc(
            s->p->c_param_args, sizeof(s->p->c_param_args[0]) * (argno + 2));
        s->p->c_param_structs =
            realloc(s->p->c_param_structs,
                    sizeof(s->p->c_param_structs[0]) * (argno + 2));
    }
    for (i = s->p->c_param_count; i <= argno; i++) {
        s->p->c_param_regs[i] = NULL;
        s->p->c_param_args[i].type = DILL_V;
        s->p->c_param_args[i].is_register = 0;
        s->p->c_param_args[i].is_immediate = 0;
        s->p->c_param_args[i].in_reg = 0;
        s->p->c_param_args[i].out_reg = 0;
        s->p->c_param_args[i].offset = 0;
        s->p->c_param_structs[i] = NULL;
    }
    s->p->c_param_count = (argno + 1);
}

extern dill_reg
dill_vparam(dill_stream s, int param_no)
{
    return param_no;
}

extern int
dill_param_reg(dill_stream s, int argno)
{
    if (argno >= s->p->c_param_count) {
        printf("Warning, dill_param_reg requested param %d, largest is %d\n",
               argno, s->p->c_param_count - 1);
        return -1;
    }
    return s->p->c_param_args[argno].in_reg;
}

extern void
dill_param_alloc(dill_stream s, int argno, int type, dill_reg* reg_p)
{
    extend_params(s, argno);
    s->p->c_param_regs[argno] = reg_p;
    s->p->c_param_args[argno].type = type;
}

extern void
dill_param_struct_alloc(dill_stream s,
                        int argno,
                        int type,
                        dill_parameter_type* struct_p)
{
    extend_params(s, argno);
    s->p->c_param_structs[argno] = struct_p;
    s->p->c_param_args[argno].type = type;
}

extern void
dill_start_simple_proc(dill_stream s, const char* subr_name, int ret_type)
{
    int i;
    if (!s->p->unavail_called)
        reset_context(s);
    s->p->ret_type = ret_type;
    s->p->unavail_called = 0;
    (s->j->proc_start)(s, (char*)subr_name, s->p->c_param_count,
                       s->p->c_param_args, NULL);
    for (i = 0; i < s->p->c_param_count; i++) {
        if (s->p->c_param_regs[i] != NULL) {
            *s->p->c_param_regs[i] = s->p->c_param_args[i].in_reg;
        }
        if (s->p->c_param_structs[i] != NULL) {
            s->p->c_param_structs[i]->is_register =
                s->p->c_param_args[i].is_register;
            s->p->c_param_structs[i]->reg = s->p->c_param_args[i].in_reg;
            s->p->c_param_structs[i]->offset = s->p->c_param_args[i].offset;
        }
    }
    s->p->c_param_count = 0;
    if (s->p->c_param_regs) {
        free(s->p->c_param_regs);
        s->p->c_param_regs = NULL;
    }
    if (s->p->c_param_args) {
        free(s->p->c_param_args);
        s->p->c_param_args = NULL;
    }
    if (s->p->c_param_structs) {
        free(s->p->c_param_structs);
        s->p->c_param_structs = NULL;
    }
}

extern void*
dill_get_fp(dill_exec_handle h)
{
    return (void*)h->fp;
}

extern dill_exec_handle
dill_finalize(dill_stream s)
{
    dill_exec_handle handle = malloc(sizeof(*handle));
    (s->j->end)(s);
    s->p->save_param_count = s->p->c_param_count;
    s->p->c_param_count = 0;
    memset(handle, 0, sizeof(*handle));
    handle->fp = (void (*)())s->p->fp;
    handle->ref_count = 1;
    handle->size = 0;
    return handle;
}

extern dill_exec_handle
dill_get_handle(dill_stream s)
{
    char* native_base = s->p->native.code_base;
    dill_exec_handle handle = malloc(sizeof(*handle));
    intptr_t size = (intptr_t)s->p->native.code_limit -
                    (intptr_t)s->p->native.code_base + END_OF_CODE_BUFFER;
    s->p->native.code_base = NULL;
    if (native_base == 0) {
        native_base = s->p->code_base;
        size = (intptr_t)s->p->code_limit - (intptr_t)s->p->code_base +
               END_OF_CODE_BUFFER;
        s->p->code_base = NULL;
    }
    handle->fp = (void (*)())s->p->fp;
    handle->ref_count = 1;
    handle->size = (int)size;
    handle->code_base = native_base;
    /* below used for libffi emulation */
    handle->emu_args = s->p->emu_args;
    handle->cifp = s->p->cifp;
    handle->closure = s->p->closure;
    s->p->emu_args = NULL;
    s->p->cifp = NULL;
    s->p->closure = NULL;
    return handle;
}

extern void
dill_free_handle(dill_exec_handle handle)
{
    handle->ref_count--;
    if (handle->ref_count > 0)
        return;
    if (handle->size != 0) {
        if (handle->code_base) {
#ifdef USE_MMAP_CODE_SEG
            if (munmap(handle->code_base, handle->size) == -1)
                perror("unmap 1");
#else
            free(handle->code_base);
#endif
        }
    }
    handle->code_base = NULL;
    handle->size = 0;
#ifdef BUILD_EMULATOR
    free_emulator_handler_bits(handle);
#endif
    free(handle);
}

extern void
dill_ref_handle(dill_exec_handle handle)
{
    handle->ref_count++;
}

static char*
dill_build_package(dill_stream s, int* pkg_len)
{
    struct call_table* t = &s->p->call_table;
    int pkg_size = sizeof(struct dill_pkg_1);
    char* pkg;
    int i;

    pkg_size = (pkg_size + 7) & -8; /* round up to mod 8 */
    pkg = malloc(pkg_size);
    memset(pkg, 0, pkg_size);
    ((struct dill_pkg_1*)pkg)->magic = 0xbeef;
    ((struct dill_pkg_1*)pkg)->pkg_version = 1;
    ((struct dill_pkg_1*)pkg)->symbol_count = t->call_count;
    for (i = 0; i < t->call_count; i++) {
        int call_len =
            (int)(sizeof(int) + strlen(t->call_locs[i].xfer_name) + 1);
        char* call_loc;

        call_len = (call_len + 7) & -8; /* round up to mod 8 */
        pkg = realloc(pkg, pkg_size + call_len);
        call_loc = pkg + pkg_size;
        pkg_size += call_len;
        *((int*)call_loc) = t->call_locs[i].loc;
        *((int*)(call_loc + call_len - 4)) = 0; /* zero last bit */
        strcpy(call_loc + 4, t->call_locs[i].xfer_name);
    }
    pkg = realloc(pkg, pkg_size + dill_code_size(s));
    ((struct dill_pkg_1*)pkg)->code_size = dill_code_size(s);
    ((struct dill_pkg_1*)pkg)->code_offset = pkg_size;
    memcpy(pkg + pkg_size, s->p->code_base, dill_code_size(s));
    pkg_size += dill_code_size(s);
    *pkg_len = pkg_size;
    ((struct dill_pkg_1*)pkg)->entry_offset =
        (short)((char*)s->p->fp - (char*)s->p->code_base);
    return pkg;
}

extern char*
dill_finalize_package(dill_stream s, int* pkg_len)
{
    (s->j->package_end)(s);
    s->p->save_param_count = s->p->c_param_count;
    s->p->c_param_count = 0;
    char* p = dill_build_package(s, pkg_len);
    return p;
}

extern int
dill_code_size(dill_stream s)
{
    char* native_base = s->p->native.code_base;
    if (native_base == 0) {
        native_base = s->p->code_base;
    }
    return (int)((char*)s->p->cur_ip - native_base);
}

extern void*
dill_clone_code(dill_stream s, void* new_base, int size)
{
    return (s->j->clone_code)(s, new_base, size);
}

extern void*
dill_take_code(dill_stream s)
{
    void* ret = s->p->code_base;
    s->p->code_base = NULL;
    s->p->native.code_base = NULL;
    return ret;
}

extern int
dill_alloc_label(dill_stream s, char* name)
{
    struct branch_table* t = &s->p->branch_table;
    if (t->next_label == t->max_alloc) {
        t->max_alloc++;
        t->label_locs = realloc(t->label_locs, sizeof(int) * t->max_alloc);
        t->label_name = realloc(t->label_name, sizeof(char*) * t->max_alloc);
    }
    t->label_locs[t->next_label] = -1;
    t->label_name[t->next_label] = NULL;
    if (name)
        t->label_name[t->next_label] = strdup(name);
    return t->next_label++;
}

extern void
dill_mark_label(dill_stream s, int label)
{
    struct branch_table* t = &s->p->branch_table;
    int label_loc = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    t->label_locs[label] = label_loc;
    if (s->j->mark_label)
        (s->j->mark_label)(s, 0, 0, 0, 0, label);
    if (s->dill_debug) {
        if (t->label_name[label] == NULL) {
            printf("L%d:\n", label);
        } else {
            printf("L%d<%s>:\n", label, t->label_name[label]);
        }
    }
}

extern int
dill_is_label_mark(dill_stream s)
{
    struct branch_table* t = &s->p->branch_table;
    int i;
    int current_loc = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    for (i = 0; i < t->max_alloc; i++) {
        if (t->label_locs[0] == current_loc)
            return 1;
    }
    return 0;
}

extern void
dill_mark_branch_location(dill_stream s, int label)
{
    struct branch_table* t = &s->p->branch_table;
    int branch_loc = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);

    if (t->branch_count == t->branch_alloc) {
        t->branch_alloc++;
        t->branch_locs = realloc(
            t->branch_locs, sizeof(struct branch_location) * t->branch_alloc);
    }
    t->branch_locs[t->branch_count].label = label;
    t->branch_locs[t->branch_count].loc = branch_loc;
    t->branch_count++;
}

extern void
dill_mark_ret_location(dill_stream s)
{
    struct ret_table* t = &s->p->ret_table;
    int ret_loc = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);

    if (t->ret_count == t->ret_alloc) {
        t->ret_alloc++;
        t->ret_locs = realloc(t->ret_locs, sizeof(int) * t->ret_alloc);
    }
    t->ret_locs[t->ret_count] = ret_loc;
    t->ret_count++;
}

extern void
dill_mark_call_location(dill_stream s,
                        const char* xfer_name,
                        void* xfer_address)
{
    struct call_table* t = &s->p->call_table;
    int call_loc = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);

    if (t->call_count == t->call_alloc) {
        t->call_alloc++;
        t->call_locs =
            realloc(t->call_locs, sizeof(struct call_location) * t->call_alloc);
    }
    t->call_locs[t->call_count].loc = call_loc;
    t->call_locs[t->call_count].xfer_addr = xfer_address;
    t->call_locs[t->call_count].xfer_name = xfer_name;
    t->call_locs[t->call_count].mach_info = NULL;
    t->call_count++;
}

extern int
dill_add_const(dill_stream s, void* addr, int size)
{
    struct branch_table* t = &s->p->branch_table;
    int offset = t->data_segment_size;
    t->data_segment = realloc(t->data_segment, offset + size);
    memcpy(t->data_segment + offset, addr, size);
    t->data_segment_size += size;
    return offset;
}

static arg_info_list
translate_arg_str(const char* string, int* count)
{
    int cnt = 0;
    arg_info_list list = malloc(sizeof(list[0]));

    while (string && (*string != 0)) {
        if (*string != '%') {
            fprintf(stderr, "invalid format, expected %%, got \"%c\"\n",
                    *string);
            return NULL;
        }
        string++;
        list[cnt].is_register = 0;
        list[cnt].is_immediate = (isupper((int)(*string)) != 0);
        list[cnt].in_reg = 0;
        list[cnt].out_reg = 0;
        list[cnt].offset = 0;
        list[cnt].type = 0;
        switch (tolower(*string)) {
        case 'u':
            string++;
            switch (tolower(*string)) {
            case 'l':
                list[cnt].type = DILL_UL;
                string++;
                break;
            case 's':
                list[cnt].type = DILL_US;
                string++;
                break;
            case 'c':
                list[cnt].type = DILL_UC;
                string++;
                break;
            case '%':
            case '\0':
                list[cnt].type = DILL_U;
                break;
            default:
                fprintf(stderr,
                        "invalid format, unexpect char \"%c\" after %%u\n",
                        *string);
            }
            break;
        case 'p':
            list[cnt].type = DILL_P;
            string++;
            break;
        case 'l':
            list[cnt].type = DILL_L;
            string++;
            break;
        case 'i':
            list[cnt].type = DILL_I;
            string++;
            break;
        case 's':
            list[cnt].type = DILL_S;
            string++;
            break;
        case 'c':
            list[cnt].type = DILL_C;
            string++;
            break;
        case 'b':
            list[cnt].type = DILL_B;
            string++;
            break;
        case 'f':
            list[cnt].type = DILL_F;
            string++;
            break;
        case 'd':
            list[cnt].type = DILL_D;
            string++;
            break;
        case 'e':
            string++;
            if (tolower(*string) == 'c') {
                if (cnt == 0) {
                    list[cnt].type = DILL_EC;
                    string++;
                    break;
                } else {
                    fprintf(stderr, "%%ec format must be first format\n");
                }
            } else {
                fprintf(stderr,
                        "invalid format, unexpect char \"%c\" after %%e\n",
                        *string);
            }
            string++;
            break;
        default:
            fprintf(stderr, "invalid format, unexpect char \"%c\" after %%\n",
                    *string);
            string++;
        }
        cnt++;
        list = realloc(list, sizeof(list[0]) * (cnt + 1));
    }
    *count = cnt;
    return list;
}

extern void
dill_start_proc(dill_stream s, char* name, int ret_type, char* arg_str)
{
    int arg_count = 0;
    arg_info_list args;
    if (!s->p->unavail_called)
        reset_context(s);
    s->p->c_param_count = 0;
    s->p->ret_type = ret_type;
    if (s->p->c_param_args != NULL) {
        free(s->p->c_param_args);
        s->p->c_param_args = NULL;
    }
    args = s->p->c_param_args = translate_arg_str(arg_str, &arg_count);
    s->p->c_param_count = arg_count;
    (s->j->proc_start)(s, name, arg_count, args, NULL);
}

void
init_code_block(dill_stream s)
{
    static unsigned long size = INIT_CODE_SIZE;
#ifdef USE_MMAP_CODE_SEG
#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif
    static long ps = -1;
    if (ps == -1) {
        ps = (getpagesize());
    }
    if (ps > size)
        size = ps;
    s->p->code_base = (void*)mmap(0, 4096 /*INIT_CODE_SIZE*/,
                                  PROT_EXEC | PROT_READ | PROT_WRITE,
                                  MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
    if (s->p->code_base == (void*)-1)
        perror("mmap");
#else
    s->p->code_base = (void*)malloc(size);
#endif
    s->p->code_limit = ((char*)s->p->code_base) + size - END_OF_CODE_BUFFER;
}

static void
free_code_blocks(dill_stream s)
{
#ifdef USE_MMAP_CODE_SEG
    if (s->p->code_base) {
        int size =
            (long)s->p->code_limit - (long)s->p->code_base + END_OF_CODE_BUFFER;
        if (munmap(s->p->code_base, size) == -1)
            perror("unmap 1");
    }
    if (s->p->virtual.code_base &&
        (s->p->virtual.code_base != s->p->code_base)) {
        int vsize = (long)s->p->virtual.code_limit -
                    (long)s->p->virtual.code_base + END_OF_CODE_BUFFER;
        if (munmap(s->p->code_base, vsize) == -1)
            perror("unmap v");
    }
    if (s->p->native.code_base && (s->p->native.code_base != s->p->code_base)) {
        int nsize = (long)s->p->native.code_limit -
                    (long)s->p->native.code_base + END_OF_CODE_BUFFER;
        if (munmap(s->p->code_base, nsize) == -1)
            perror("unmap n");
    }
#else
    if (s->p->code_base)
        free(s->p->code_base);
    if (s->p->virtual.code_base && (s->p->virtual.code_base != s->p->code_base))
        free(s->p->virtual.code_base);
    if (s->p->native.code_base && (s->p->native.code_base != s->p->code_base))
        free(s->p->native.code_base);
#endif
}

extern void
extend_dill_stream(dill_stream s)
{
    int size = (int)((char*)s->p->code_limit - (char*)s->p->code_base +
                     END_OF_CODE_BUFFER);
    int cur_ip = (int)((char*)s->p->cur_ip - (char*)s->p->code_base);
    int new_size = size * 2;
#ifdef USE_MMAP_CODE_SEG
    {
        void* old = s->p->code_base;
        void* new = mmap(0, new_size, PROT_EXEC | PROT_READ | PROT_WRITE,
                         MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
        if (new == (void*)-1)
            perror("mmap1");
        memcpy(new, old, size);
        s->p->code_base = new;
        if (munmap(old, size) == -1)
            perror("munmap exp");
    }
#else
    s->p->code_base = realloc(s->p->code_base, new_size);
#endif
    s->p->cur_ip = ((char*)s->p->code_base) + cur_ip;
    s->p->code_limit = ((char*)s->p->code_base) + new_size - END_OF_CODE_BUFFER;
}

extern jmp_table
alloc_dill_jump_table()
{
    fprintf(stderr,
            "Don't use alloc_dill_jump_table().  Rewrite to use "
            "DECLARE_JUMP_TABLE() and FILL_JUMP_STRUCTURE()\n");
    return NULL;
}

#define BIT_ON(i, bitvec)                                                      \
    (((((long)1) << reg) & bitvec[0]) == (((long)1) << reg))
#define SET_BIT(i, bitvec) (bitvec[0] |= (((long)1) << reg))
#define RESET_BIT(i, bitvec) (bitvec[0] &= (~(((long)1) << reg)))

void
dill_markused(dill_stream s, int type, int reg)
{
    switch (type) {
    case DILL_D:
    case DILL_F:
        SET_BIT(reg, s->p->var_f.used);
        SET_BIT(reg, s->p->tmp_f.used);
        break;
    default:
        SET_BIT(reg, s->p->var_i.used);
        SET_BIT(reg, s->p->tmp_i.used);
        break;
    }
}

int
dill_mustsave(reg_set* regs, int reg)
{
    return BIT_ON(reg, regs->mustsave);
}

int
dill_wasused(reg_set* regs, int reg)
{
    return BIT_ON(reg, regs->used);
}

static int
reg_alloc(reg_set* regs)
{
    int reg;
    if (regs->avail[0] == 0) {
        return -1;
    }
    reg = 0;
    while (!BIT_ON(reg, regs->avail)) {
        reg++;
    }
    RESET_BIT(reg, regs->avail);
    return reg;
}

extern void
dill_alloc_specific(dill_stream s, dill_reg reg, int type, int class)
{
    switch (type) {
    case DILL_D:
    case DILL_F:
        if (class == DILL_VAR) {
            SET_BIT(reg, s->p->var_f.avail);
            RESET_BIT(reg, s->p->var_f.mustsave);
        } else {
            SET_BIT(reg, s->p->tmp_f.avail);
        }
        break;
    default:
        if (class == DILL_VAR) {
            SET_BIT(reg, s->p->var_i.avail);
            RESET_BIT(reg, s->p->var_i.mustsave);
        } else {
            SET_BIT(reg, s->p->tmp_i.avail);
        }
        break;
    }
}

extern void
dill_dealloc_specific(dill_stream s, dill_reg reg, int type, int class)
{
    switch (type) {
    case DILL_D:
    case DILL_F:
        if (class == DILL_VAR) {
            RESET_BIT(reg, s->p->var_f.avail);
            SET_BIT(reg, s->p->var_f.mustsave);
        } else {
            RESET_BIT(reg, s->p->tmp_f.avail);
        }
        break;
    default:
        if (class == DILL_VAR) {
            RESET_BIT(reg, s->p->var_i.avail);
            SET_BIT(reg, s->p->var_i.mustsave);
        } else {
            RESET_BIT(reg, s->p->tmp_i.avail);
        }
        break;
    }
}

extern void
dill_raw_unavailreg(dill_stream s, int type, dill_reg reg)
{
    if (s->p->unavail_called == 0) {
        reset_context(s);
        s->p->unavail_called = 1;
    }
    switch (type) {
    case DILL_D:
    case DILL_F:
        RESET_BIT(reg, s->p->var_f.avail);
        RESET_BIT(reg, s->p->tmp_f.avail);
        break;
    default:
        RESET_BIT(reg, s->p->var_i.avail);
        RESET_BIT(reg, s->p->tmp_i.avail);
        break;
    }
}

extern void
dill_raw_availreg(dill_stream s, int type, dill_reg reg)
{
    switch (type) {
    case DILL_D:
    case DILL_F:
        if (BIT_ON(reg, s->p->tmp_f.members)) {
            SET_BIT(reg, s->p->tmp_f.avail);
        } else if (BIT_ON(reg, s->p->var_f.members)) {
            SET_BIT(reg, s->p->var_f.avail);
        } else {
            printf("mk avail not in set error %d\n", reg);
        }
        break;
    default:
        if (BIT_ON(reg, s->p->tmp_i.members)) {
            SET_BIT(reg, s->p->tmp_i.avail);
        } else if (BIT_ON(reg, s->p->var_i.members)) {
            SET_BIT(reg, s->p->var_i.avail);
        } else {
            printf("mk avail not in set error %d\n", reg);
        }
        break;
    }
}

extern int
dill_getreg(dill_stream s, int typ)
{
    s->p->vregs =
        realloc(s->p->vregs, (s->p->vreg_count + 1) * sizeof(vreg_info));
    s->p->vregs[s->p->vreg_count].typ = typ;
    s->p->vregs[s->p->vreg_count].use_info.use_count = 0;
    s->p->vregs[s->p->vreg_count].use_info.def_count = 0;
    s->p->vregs[s->p->vreg_count].offset = 0xdeadbeef;
    return ((s->p->vreg_count++) + 100);
}

extern int
dill_getvblock(dill_stream s, int size)
{
    s->p->vregs =
        realloc(s->p->vregs, (s->p->vreg_count + 1) * sizeof(vreg_info));
    s->p->vregs[s->p->vreg_count].typ = DILL_B;
    s->p->vregs[s->p->vreg_count].offset = size;
    s->p->vregs[s->p->vreg_count].use_info.use_count = 0;
    s->p->vregs[s->p->vreg_count].use_info.def_count = 0;
    return ((s->p->vreg_count++) + 100);
}

extern int
dill_raw_getreg(dill_stream s, dill_reg* reg_p, int type, int class)
{
    int reg = -1;
    switch (type) {
    case DILL_D:
    case DILL_F:
        if (class == DILL_VAR) {
            if ((reg = reg_alloc(&s->p->var_f)) == -1) {
                reg = reg_alloc(&s->p->tmp_f);
            }
            if (reg != -1) {
                SET_BIT(reg, s->p->tmp_f.mustsave);
                SET_BIT(reg, s->p->tmp_f.used);
            }
            *reg_p = reg;
            return (reg != -1);
        } else {
            if ((reg = reg_alloc(&s->p->tmp_f)) == -1) {
                reg = reg_alloc(&s->p->var_f);
            }
            if (reg != -1) {
                SET_BIT(reg, s->p->tmp_f.used);
            }
            *reg_p = reg;
            return (reg != -1);
        }
    default:
        if (class == DILL_VAR) {
            if ((reg = reg_alloc(&s->p->var_i)) == -1) {
                reg = reg_alloc(&s->p->tmp_i);
                if (reg != -1) {
                    SET_BIT(reg, s->p->tmp_i.mustsave);
                    SET_BIT(reg, s->p->tmp_i.used);
                }
            } else {
                SET_BIT(reg, s->p->var_i.used);
            }
            *reg_p = reg;
            return (reg != -1);
        } else {
            if ((reg = reg_alloc(&s->p->tmp_i)) == -1) {
                reg = reg_alloc(&s->p->var_i);
            }
            if (reg != -1) {
                SET_BIT(reg, s->p->tmp_i.used);
            }
            *reg_p = reg;
            return (reg != -1);
        }
    }
}

extern void
dill_raw_putreg(dill_stream s, dill_reg reg, int type)
{
    switch (type) {
    case DILL_F:
    case DILL_D:
        if (BIT_ON(reg, s->p->tmp_f.members)) {
            dill_alloc_specific(s, reg, type, DILL_TEMP);
        } else if (BIT_ON(reg, s->p->var_f.members)) {
            dill_alloc_specific(s, reg, type, DILL_VAR);
        } else {
            /*
             * special case, if the put reg is invalid (-1) and
             * the member set is empty, give no warning.
             */
            if ((s->p->var_f.members[0] != 0) || (reg != -1)) {
                printf("Putreg not in set error %d\n", reg);
            }
        }
        break;
    default:
        if (BIT_ON(reg, s->p->tmp_i.members)) {
            dill_alloc_specific(s, reg, type, DILL_TEMP);
        } else if (BIT_ON(reg, s->p->var_i.members)) {
            dill_alloc_specific(s, reg, type, DILL_VAR);
        } else {
            printf("Putreg not in set error %d\n", reg);
        }
        break;
    }
}

extern int
dill_do_reverse_vararg_push(dill_stream s)
{
    if (s->j->do_reverse_push) {
        s->p->doing_reverse_push = 1;
        return 1;
    }
    return 0;
}

extern void
dill_end_vararg_push(dill_stream s)
{
    s->p->doing_reverse_push = 0;
}

extern int
dill_type_of(dill_stream s, int vreg)
{
    if (vreg >= 100) {
        return s->p->vregs[vreg - 100].typ;
    } else {
        return s->p->c_param_args[vreg].type;
    }
}

typedef union {
    unsigned int u;
    int i;
    long l;
    unsigned long ul;
    void* p;
    float f;
    double d;
} type_union;

static void
do_vararg_push(dill_stream s, const char* arg_str, va_list ap)
{
    int i, arg_count;
    int reverse = 0;
    arg_info_list args = translate_arg_str(arg_str, &arg_count);
    type_union value[256];

    dill_push_init(s);

    reverse = dill_do_reverse_vararg_push(s);

    for (i = 0; i < arg_count; i++) {
        if (!args[i].is_immediate) {
            value[i].i = va_arg(ap, int);
        } else {
            switch (args[i].type) {
            case DILL_UC:
            case DILL_US:
            case DILL_U:
                value[i].u = va_arg(ap, unsigned);
                break;
            case DILL_C:
            case DILL_S:
            case DILL_I:
                value[i].i = va_arg(ap, int);
                break;
            case DILL_L:
                value[i].l = va_arg(ap, long);
                break;
            case DILL_UL:
                value[i].ul = va_arg(ap, unsigned long);
                break;
            case DILL_P:
                value[i].p = va_arg(ap, void*);
                break;
            case DILL_F:
                value[i].f = (float)va_arg(ap, double);
                break;
            case DILL_D:
                value[i].d = va_arg(ap, double);
                break;
            default:
                fprintf(stderr, "do_push_args: unknown type\n");
                exit(1);
            }
        }
    }

    /* push all arguments */
    for (i = 0; i < arg_count; i++) {
        int arg = i;
        if (reverse) {
            arg = arg_count - i - 1;
        }
        if (!args[arg].is_immediate) {
            dill_push_arg(s, args[arg].type, value[arg].i);
        } else {
            switch (args[arg].type) {
            case DILL_UC:
            case DILL_US:
            case DILL_U:
                dill_push_argui(s, value[arg].u);
                break;
            case DILL_C:
            case DILL_S:
            case DILL_I:
                dill_push_argii(s, value[arg].i);
                break;
            case DILL_L:
                dill_push_argli(s, value[arg].l);
                break;
            case DILL_UL:
                dill_push_arguli(s, value[arg].ul);
                break;
            case DILL_P:
                dill_push_argpi(s, value[arg].p);
                break;
            case DILL_F:
                dill_push_argfi(s, value[arg].f);
                break;
            case DILL_D:
                dill_push_argdi(s, value[arg].d);
                break;
            default:
                fprintf(stderr, "do_push_args: unknown type\n");
                exit(1);
            }
        }
    }
    free(args);
}

void
dill_scallv(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    va_list ap;

    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    (void)s->j->calli(s, DILL_V, ptr, (char*)name);
    va_end(ap);
}

int
dill_scalli(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_I, ptr, (char*)name);
    va_end(ap);
    return ret_reg;
}

int
dill_scallu(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_U, ptr, (char*)name);
    va_end(ap);
    return ret_reg;
}

int
dill_scallp(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_P, ptr, (char*)name);
    va_end(ap);
    return ret_reg;
}

int
dill_scallul(dill_stream s,
             void* ptr,
             const char* name,
             const char* arg_str,
             ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_UL, ptr, (char*)name);
    va_end(ap);
    return ret_reg;
}

int
dill_scalll(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_L, ptr, name);
    va_end(ap);
    return ret_reg;
}

int
dill_scallf(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_F, ptr, name);
    va_end(ap);
    return ret_reg;
}

int
dill_scalld(dill_stream s,
            void* ptr,
            const char* name,
            const char* arg_str,
            ...)
{
    int ret_reg;
    va_list ap;
    va_start(ap, arg_str);
    do_vararg_push(s, arg_str, ap);
    ret_reg = s->j->calli(s, DILL_D, ptr, name);
    va_end(ap);
    return ret_reg;
}

extern void
dill_error(char* msg)
{
    printf("%s", msg);
}

extern void
dill_pbr(dill_stream s,
         int op_type,
         int data_type,
         dill_reg src1,
         dill_reg src2,
         int label)
{
    int index;
    if ((op_type < dill_eq_code) || (op_type > dill_ne_code)) {
        printf("Bad op type in dill_pbr\n");
    }
    index = data_type + 11 * op_type;
    (s->j->jmp_b)[index](s, s->j->b_data[index].data1,
                         s->j->b_data[index].data2, src1, src2, label);
}

extern void
dill_pcompare(dill_stream s,
              int op_type,
              int data_type,
              dill_reg dest,
              dill_reg src1,
              dill_reg src2)
{
    int index;
    if ((op_type < dill_eq_code) || (op_type > dill_ne_code)) {
        printf("Bad op type in dill_pcompare\n");
    }
    index = data_type + 11 * op_type;
    (s->j->jmp_c)[index](s, s->j->c_data[index].data1,
                         s->j->c_data[index].data2, dest, src1, src2);
}

extern void
dill_dump_reg(dill_stream s, int typ, int reg)
{
    s->j->print_reg(s, typ, reg);
}

#if !defined(HAVE_DIS_ASM_H) || defined(NO_DISASSEMBLER)
struct disassemble_info {
    void* junk;
};
#endif

#if defined(HAVE_DIS_ASM_H) && !defined(NO_DISASSEMBLER)
/* GENERIC BINUTILS DISASSEMBLER */
/* include some things from libbfd so we don't have to have it all */
#include "dis-asm.h"
#ifndef NO_INLINED_BFD_PROCS
#if defined(BFD64)
bfd_vma bfd_getl64(addr)
#ifdef BFD_BYTE
    register const bfd_byte* addr;
#else
    const void* addr;
#endif
{
#if SIZEOF_LONG == 8
    register const unsigned char* addrp = addr;
    bfd_vma low, high;
    high =
        (((((((addrp[7] << 8) | addrp[6]) << 8) | addrp[5]) << 8) | addrp[4]));

    low = ((((((((bfd_vma)addrp[3] << 8) | addrp[2]) << 8) | addrp[1]) << 8) |
            addrp[0]));

    return high << 32 | low;
#else
    return 0;
#endif
}
#endif

bfd_vma bfd_getb32(addr)
#ifdef BFD_BYTE
    register const bfd_byte* addr;
#else
    const void* addr;
#endif
{
    register const unsigned char* addrp = addr;
    return (((((bfd_vma)addrp[0] << 8) | addrp[1]) << 8) | addrp[2]) << 8 |
           addrp[3];
}

bfd_vma bfd_getl32(addr)
#ifdef BFD_BYTE
    register const bfd_byte* addr;
#else
    const void* addr;
#endif
{
    register const unsigned char* addrp = addr;
    return (((((bfd_vma)addrp[3] << 8) | addrp[2]) << 8) | addrp[1]) << 8 |
           addrp[0];
}

bfd_vma bfd_getb16(addr)
#ifdef BFD_BYTE
    register const bfd_byte* addr;
#else
    const void* addr;
#endif
{
    register const unsigned char* addrp = addr;
    return (addrp[0] << 8) | addrp[1];
}

bfd_vma bfd_getl16(addr)
#ifdef BFD_BYTE
    register const bfd_byte* addr;
#else
    const void* addr;
#endif
{
    register const unsigned char* addrp = addr;
    return (addrp[1] << 8) | addrp[0];
}

enum bfd_architecture bfd_get_arch(abfd) bfd* abfd;
{
    return abfd->arch_info->arch;
}

void* xmalloc(size) int size;
{
    return (void*)malloc(size);
}

char* xstrdup(s) const char* s;
{
    int len;
    char* ret;

    len = strlen(s);
    ret = xmalloc(len + 1);
    strcpy(ret, s);
    return ret;
}
#endif

extern void
dump_dill_insn(dill_stream s, void* p)
{
    struct disassemble_info info;
    int l;

    if (s->j->init_disassembly(s, &info) == 0)
        return;

    printf("%lx  - %x - ", (unsigned long)p, (unsigned int)*(int*)p);
    l = s->j->print_insn(s, &info, p);
    printf("\n");
}

extern void
dump_cur_dill_insn(dill_stream s)
{
    dump_dill_insn(s, s->p->cur_ip);
}

#else
extern void
dump_cur_dill_insn(dill_stream s)
{
    struct disassemble_info info;
    s->j->init_disassembly(s, &info);

    void* p = s->p->cur_ip;
    printf("%p  - %x - ", p, (unsigned)*(int*)p);
    (void)s->j->print_insn(s, &info, (void*)p);
    printf("\n");
}
#endif

extern void
dill_dump(dill_stream s)
{
    struct disassemble_info info;
    void* base = s->p->code_base;
    int native_missing = 0;

    if ((base != s->p->virtual.code_base) &&
        (s->p->virtual.code_base != NULL) &&
        (s->p->virtual.mach_jump != NULL)) {
        /* Section to dump virtual code base *after* dcg completion */
        /* only do this if we're not currently in the middle of virtual
         * generation */
        void* code_limit = s->p->virtual.cur_ip;
        base = s->p->virtual.code_base;
        s->p->virtual.mach_jump->init_disassembly(s, &info);
        void* p;
        int l;
        int insn_count = 0;
        printf("\nDILL virtual instruction stream\n\n");
        for (p = base; p < code_limit;) {
            printf("%p  - %x - ", p, (unsigned)*(int*)p);
            l = s->p->virtual.mach_jump->print_insn(s, &info, (void*)p);
            printf("\n");
            if (l == -1)
                return;
            p = (char*)p + l;
            insn_count++;
        }
        printf("\nDumped %d virtual instructions\n\n", insn_count);
    }
#if defined(NO_DISASSEMBLER)
    native_missing = 1;
#endif
    base = s->p->code_base;
    if (base == NULL) {
        base = s->p->native.code_base;
    }
    if (base == NULL) {
        printf("No code to dump\n");
        return;
    }
    /* if ((s->j != s->p->virtual.mach_jump) && native_missing) { */
    /*     printf("No native disassembler available\n"); */
    /*     return; */
    /* } */
    (void)native_missing;
    if (s->j->init_disassembly(s, &info) == 0) {
        printf("No native disassembler available\n");
    } else {
        void* p;
        int l;
        int insn_count = 0;
        if ((s->j != s->p->virtual.mach_jump) && (s->p->fp != NULL))
            base = s->p->fp;
        for (p = base; (char*)p < s->p->cur_ip;) {
            int i;
            struct branch_table* t = &s->p->branch_table;
            for (i = 0; i < t->next_label; i++) {
                if (t->label_locs[i] == ((char*)p - (char*)base)) {
                    printf("L%d:\n", i);
                }
            }
            if (p == s->p->fp) {
                printf("Function entry point:\n");
            }
            printf("%p  - %x - ", p, (unsigned)*(int*)p);
            l = s->j->print_insn(s, &info, (void*)p);
            printf("\n");
            if (l <= 0)
                return;
            p = (char*)p + l;
            insn_count++;
        }
        printf("\nDumped %d instructions\n\n", insn_count);
    }
}

#undef malloc
#undef realloc

extern void*
dill_malloc(size_t size)
{
    void* tmp = malloc(size);
    if (tmp)
        return tmp;
    fprintf(stderr, "Dill out of memory, exiting\n");
    exit(1);
}

extern void*
dill_realloc(void* ptr, size_t size)
{
    void* tmp = realloc(ptr, size);
    if (tmp)
        return tmp;
    fprintf(stderr, "Dill out of memory, exiting\n");
    exit(1);
}
