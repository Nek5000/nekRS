#include "config.h"
#ifndef _MSC_VER
#include <dlfcn.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "dlloader.h"

static char **search_list = NULL;
#ifdef _MSC_VER
#include <windows.h>

static struct {
    long lasterror;
    const char* err_rutin;
} var = {
    0,
    NULL
};

void* dlopen(const char* filename, int flags)
{
    HINSTANCE hInst;

    hInst = LoadLibrary(filename);
    if (hInst == NULL) {
	var.lasterror = GetLastError();
	var.err_rutin = "dlopen";
    }
    return hInst;
}

int dlclose(void* handle)
{
    BOOL ok;
    int rc = 0;

    ok = FreeLibrary((HINSTANCE)handle);
    if (!ok) {
	var.lasterror = GetLastError();
	var.err_rutin = "dlclose";
	rc = -1;
    }
    return rc;
}

void* dlsym(void* handle, const char* name)
{
    FARPROC fp;

    fp = GetProcAddress((HINSTANCE)handle, name);
    if (!fp) {
	var.lasterror = GetLastError();
	var.err_rutin = "dlsym";
    }
    return (void*)(intptr_t)fp;
}

const char* dlerror(void)
{
    static char errstr[88];

    if (var.lasterror) {
	sprintf(errstr, "%s error #%ld", var.err_rutin, var.lasterror);
	return errstr;
    }
    else {
	return NULL;
    }
}
#endif

void
CMdladdsearchdir(char *string)
{
    int count = 0;
    if (search_list == NULL) {
	search_list = malloc(2*sizeof(char*));
    } else {
	while(search_list[count] != NULL) count++;
	search_list = realloc(search_list, (count+2)*sizeof(char*));
    }
    search_list[count] = strdup(string);
    search_list[count+1] = NULL;
}

typedef struct {
    void *dlopen_handle;
    char *lib_prefix;
} *dlhandle;

static int dlopen_verbose = -1;

void
CMset_dlopen_verbose(int verbose)
{
    dlopen_verbose = verbose;
}
#undef dlopen

void *
CMdlopen(void *CMTrace_filev, char *in_lib, int mode)
{
#if NO_DYNAMIC_LINKING
    return NULL;
#else
    dlhandle dlh;
    void *handle = NULL;
    char *tmp;
    char *lib;
    FILE *CMTrace_file = (FILE*)CMTrace_filev;
    if (dlopen_verbose == -1) {
	dlopen_verbose = (getenv("CMTransportVerbose") != NULL);
    }
    tmp = strrchr(in_lib, '.');
    if (dlopen_verbose) fprintf(CMTrace_file, "Trying to dlopen %s\n", in_lib);
    if (tmp && (strcmp(tmp, ".la") == 0)) {
	/* can't open .la files */
	lib = malloc(strlen(in_lib) + strlen(MODULE_EXT) + 8);
	strcpy(lib, in_lib);
	strcpy(strrchr(lib, '.'), MODULE_EXT);
	if (dlopen_verbose) fprintf(CMTrace_file, "Dlopen module name replaced, now %s\n", lib);
    } else {
	lib = strdup(in_lib);
    }
    char **list = search_list;
    while(list && (list[0] != NULL)) {
        char *tmp = malloc(strlen(list[0]) + strlen(lib) + 2);
	sprintf(tmp, "%s/%s", list[0], lib);
	handle = dlopen(tmp, RTLD_LAZY);
	const char *err = dlerror();
	if (dlopen_verbose) {
	    if (err) {
		fprintf(CMTrace_file, "Failed to dlopen %s, error is %s\n", tmp, err);
	    } else {
		fprintf(CMTrace_file, "DLopen of %s succeeded\n", tmp);
	    }
	}
	free(tmp);
 	list++;
	if (handle) list = NULL; // fall out
    }
    if (!handle) {
        handle = dlopen(lib, RTLD_LAZY);
	const char *err = dlerror();
	if (dlopen_verbose) {
	    if (err) {
		fprintf(CMTrace_file, "Failed to dlopen %s, error is %s\n", tmp, err);
	    } else {
		fprintf(CMTrace_file, "DLopen of %s succeeded\n", tmp);
	    }
	}
    }
    if (!handle) return NULL;
    dlh = malloc(sizeof(*dlh));
    tmp = strrchr(lib, '/'); /* find name start */
    if (!tmp) tmp = lib;

    char *cm_lib_prefix;
    if(strlen(CM_LIBRARY_PREFIX) > 0 &&
      (cm_lib_prefix = strstr(tmp, CM_LIBRARY_PREFIX))) {
      dlh->lib_prefix = malloc(strlen(tmp) - strlen(CM_LIBRARY_PREFIX) + 4);
      strncpy(dlh->lib_prefix, tmp, cm_lib_prefix-tmp);
      strcpy(dlh->lib_prefix + (cm_lib_prefix - tmp),
          cm_lib_prefix + strlen(CM_LIBRARY_PREFIX));
    }
    else {
      dlh->lib_prefix = malloc(strlen(tmp) + 4);
      strcpy(dlh->lib_prefix, tmp);
    }
    tmp = strrchr(dlh->lib_prefix, '.');
    strcpy(tmp, "_LTX_");  /* kill postfix, add _LTX_ */
    dlh->dlopen_handle = handle;
    free(lib);
    return (void*)dlh;
#endif
}

void*
CMdlsym(void *vdlh, char *sym)
{
#if NO_DYNAMIC_LINKING
    return NULL;
#else
    dlhandle dlh = (dlhandle)vdlh;
    char *tmp = malloc(strlen(sym) + strlen(dlh->lib_prefix) + 3);
    void *sym_val;
    strcpy(tmp, dlh->lib_prefix);
    strcat(tmp, sym);
    sym_val = dlsym(dlh->dlopen_handle, tmp);
    if (!sym_val) {
      // try with lib prefix
      char *tmp2 = malloc(strlen(tmp) + 4);
      strcpy(tmp2, "lib");
      strcat(tmp2, tmp);
      sym_val = dlsym(dlh->dlopen_handle, tmp2);
      free(tmp2);
    }
    free(tmp);
    if (!sym_val) 
	sym_val = dlsym(dlh->dlopen_handle, sym);
    return sym_val;
#endif
}

void
CMdlclose(void *vdlh)
{
#if NO_DYNAMIC_LINKING
    return;
#else
    dlhandle dlh = (dlhandle)vdlh;
#ifdef ACTUALL_DO_DLCLOSE
    dlclose(dlh->dlopen_handle);
#endif
    free(dlh->lib_prefix);
    free(dlh);
#endif
}

void
CMdlclearsearchlist()
{
    int i = 0;
    while(search_list[i]) {
        free(search_list[i]);
	i++;
    }
    free(search_list);
}
