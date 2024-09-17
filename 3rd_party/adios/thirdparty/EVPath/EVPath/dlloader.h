
#ifndef _MSC_VER
#include <dlfcn.h>
#else
#define RTLD_GLOBAL 0x100 /* do not hide entries in this module */
#define RTLD_LOCAL  0x000 /* hide entries in this module */

#define RTLD_LAZY   0x000 /* accept unresolved externs */
#define RTLD_NOW    0x001 /* abort if module has unresolved externs */
#endif
#define lt_dlopen(x) CMdlopen(cm, x, 0)
#define lt_dladdsearchdir(x) CMdladdsearchdir(x)
#define lt_dlsym(x, y) CMdlsym(x, y)
#define lt_dlhandle void*
#define MODULE_EXT CMAKE_SHARED_MODULE_SUFFIX
extern void CMdladdsearchdir(char *dir);
extern void* CMdlopen(void *CMTrace_file, char *library, int mode);
extern void CMdlclose(void *handle);
extern void CMdlclearsearchlist();
extern void* CMdlsym(void *handle, char *symbol);
extern void CMset_dlopen_verbose(int verbose);
