#include "config.h"
#ifndef MODULE

#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdio.h>
#ifdef STDC_HEADERS
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <errno.h>
#else
#include "kernel/kcm.h"
#include "kernel/cm_kernel.h"
#include "kernel/library.h"
#endif
#include "atl.h"
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"


extern void EVfprint_version(FILE* out);
extern void CMset_dlopen_verbose(int verbose);

int CMtrace_val[CMLastTraceType] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
int CMtrace_timing = 0;
int CMtrace_PID = 0;

static int CMTrace_file_num = -1;

extern void INT_CMTrace_file_id(int ID)
{
    CMTrace_file_num = ID;
}

extern int CMtrace_init(CManager cm, CMTraceType trace_type)
{
    int i, trace = 0;
    char *str;
    CMtrace_val[0] = 0;
    CMtrace_val[EVWarning] = 1;  /* default on */
    CMtrace_val[CMControlVerbose] = (getenv("CMControlVerbose") != NULL);
    CMtrace_val[CMConnectionVerbose] = (getenv("CMConnectionVerbose") != NULL);
    CMtrace_val[CMDataVerbose] = (getenv("CMDataVerbose") != NULL);
    CMtrace_val[CMTransportVerbose] = (getenv("CMTransportVerbose") != NULL);
    CMtrace_val[CMFormatVerbose] = (getenv("CMFormatVerbose") != NULL);
    CMtrace_val[CMFreeVerbose] = (getenv("CMFreeVerbose") != NULL);
    CMtrace_val[CMAttrVerbose] = (getenv("CMAttrVerbose") != NULL);
    CMtrace_val[CMBufferVerbose] = (getenv("CMBufferVerbose") != NULL);
    CMtrace_val[EVerbose] = (getenv("EVerbose") != NULL);
    CMtrace_val[CMSelectVerbose] = (getenv("CMSelectVerbose") != NULL);    
    CMtrace_val[EVdfgVerbose] = (getenv("EVdfgVerbose") != NULL);
    CMtrace_timing = (getenv("CMTraceTiming") != NULL);
    CMtrace_PID = (getenv("CMTracePID") != NULL);
    if ((str = getenv("EVWarning")) != NULL) {
	sscanf(str, "%d", &CMtrace_val[EVWarning]);
    }
    if (getenv("CMVerbose") != NULL) {
	int j;
	for (j=0; j<CMLastTraceType; j++)
	    CMtrace_val[j] = 1;
    }
    /* for low level verbose, value overrides general CMVerbose */
    CMtrace_val[CMLowLevelVerbose] = (getenv("CMLowLevelVerbose") != NULL);

    if (getenv("CMTraceFile") != NULL) {
	CMTrace_file_num = getpid();
    }
    if (CMTrace_file_num != -1) {
	char name[40];
	static int cm_count = 0;
	if (cm_count == 0) {
	    sprintf(name, "CMTrace_output.%d", (int)CMTrace_file_num);
	} else {
	    sprintf(name, "CMTrace_output.%d_%d", (int)CMTrace_file_num, cm_count);
	}
	cm_count++;
	cm->CMTrace_file = fopen(name, "w");
	if (cm->CMTrace_file == NULL) {
	    printf("Failed to open trace file %s\n", name);
	    cm->CMTrace_file = stdout;
	} else {
	    fprintf(cm->CMTrace_file, "Trace flags set : \n");
	    if (CMtrace_val[CMAlwaysTrace]) fprintf(cm->CMTrace_file, "CMAlwaysTrace, ");
	    if (CMtrace_val[CMControlVerbose]) fprintf(cm->CMTrace_file, "CMControlVerbose, ");
	    if (CMtrace_val[CMConnectionVerbose]) fprintf(cm->CMTrace_file, "CMConnectionVerbose, ");
	    if (CMtrace_val[CMLowLevelVerbose]) fprintf(cm->CMTrace_file, "CMLowLevelVerbose, ");
	    if (CMtrace_val[CMDataVerbose]) fprintf(cm->CMTrace_file, "CMDataVerbose, ");
	    if (CMtrace_val[CMTransportVerbose]) fprintf(cm->CMTrace_file, "CMTransportVerbose, ");
	    if (CMtrace_val[CMFormatVerbose]) fprintf(cm->CMTrace_file, "CMFormatVerbose, ");
	    if (CMtrace_val[CMFreeVerbose]) fprintf(cm->CMTrace_file, "CMFreeVerbose, ");
	    if (CMtrace_val[CMAttrVerbose]) fprintf(cm->CMTrace_file, "CMAttrVerbose, ");
	    if (CMtrace_val[CMBufferVerbose]) fprintf(cm->CMTrace_file, "CMBufferVerbose, ");
	    if (CMtrace_val[EVerbose]) fprintf(cm->CMTrace_file, "EVerbose, ");
	    if (CMtrace_val[EVWarning]) fprintf(cm->CMTrace_file, "EVWarning, ");
	    if (CMtrace_val[CMSelectVerbose]) fprintf(cm->CMTrace_file, "CMSelectVerbose, ");
	    if (CMtrace_val[EVdfgVerbose]) fprintf(cm->CMTrace_file, "EVdfgVerbose, ");
	    fprintf(cm->CMTrace_file, "\n");
	}
    } else {
	cm->CMTrace_file = stdout;
    }
    for (i = 0; i < sizeof(CMtrace_val)/sizeof(CMtrace_val[0]); i++) {
	if (i!=EVWarning) trace |= CMtrace_val[i];
    }
    if (CMtrace_val[CMTransportVerbose]) {
	CMset_dlopen_verbose(1);
    }

    if (trace != 0) {
	EVfprint_version(cm->CMTrace_file);
    }
    fflush(cm->CMTrace_file);
    return CMtrace_val[trace_type];
}

/*extern int
CMtrace_on(CManager cm, CMTraceType trace_type)
{
    if (CMtrace_val[0] == -1) {
	CMtrace_init(cm);
    }

    return CMtrace_val[trace_type];
    }*/

 /*extern void
CMtrace_out(CManager cm, CMTraceType trace_type, char *format, ...)
{
#ifndef MODULE
    va_list ap;

    if (CMtrace_on(cm, trace_type)) {
	if (CMtrace_on(cm, CMLowLevelVerbose)) {
	    printf("P%lxT%lx - ", (long) getpid(), (long)thr_thread_self());
	}
#ifdef STDC_HEADERS
	va_start(ap, format);
#else
	va_start(ap);
#endif
	vfprintf(cm->CMTrace_file, format, ap);
	va_end(ap);
	fprintf(cm->CMTrace_file, "\n");
    }
#endif
}
 */
extern void
CMtransport_verbose(CManager cm, CMTraceType trace, const char *format, ...)
{
#ifndef MODULE
    va_list ap;
    if (CMtrace_on(cm, trace)) {
        if (CMtrace_PID) {
            fprintf(cm->CMTrace_file, "P%lxT%lx - ", (long) getpid(), (long)thr_thread_self());
        }
        if (CMtrace_timing) {
            TRACE_TIME_DECL;
            TRACE_TIME_GET;
            fprintf(cm->CMTrace_file, TRACE_TIME_PRINTDETAILS);
        }
#ifdef STDC_HEADERS
	va_start(ap, format);
#else
	va_start(ap);
#endif
	vfprintf(cm->CMTrace_file, format, ap);
	va_end(ap);
	(void)cm;
	fprintf(cm->CMTrace_file, "\n");
    }
#endif
}

extern void
CMtransport_trace(CManager cm, const char *format, ...)
{
#ifndef MODULE
    va_list ap;
    if (CMtrace_on(cm, CMTransportVerbose)) {
        if (CMtrace_PID) {
            fprintf(cm->CMTrace_file, "P%lxT%lx - ", (long) getpid(), (long)thr_thread_self());
        }
        if (CMtrace_timing) {
            TRACE_TIME_DECL;
            TRACE_TIME_GET;
            fprintf(cm->CMTrace_file, TRACE_TIME_PRINTDETAILS);
        }
#ifdef STDC_HEADERS
	va_start(ap, format);
#else
	va_start(ap);
#endif
	vfprintf(cm->CMTrace_file, format, ap);
	va_end(ap);
	(void)cm;
	fprintf(cm->CMTrace_file, "\n");
    }
#endif
}

extern attr_list 
CMint_create_attr_list(CManager cm, char *file, int line)
{
    attr_list list = create_attr_list();
    (void)cm;
    CMtrace_out(cm, CMAttrVerbose, "Creating attr list %p at %s:%d\n", 
		list, file, line);
    return list;
}

extern void 
CMint_free_attr_list(CManager cm, attr_list l, char *file, int line)
{
    int count = attr_list_ref_count(l);
    (void)cm;
    CMtrace_out(cm, CMAttrVerbose, "Freeing attr list %p at %s:%d, ref count was %d\n", 
		l, file, line, count);
    free_attr_list(l);
}


extern attr_list 
CMint_add_ref_attr_list(CManager cm, attr_list l, char *file, int line)
{
    int count;
    (void)cm;
    if (l == NULL) return NULL;
    count = attr_list_ref_count(l);
    CMtrace_out(cm, CMAttrVerbose, "Adding ref attr list %p at %s:%d, ref count now %d\n", 
		l, file, line, count+1);
    return add_ref_attr_list(l);
}

extern attr_list 
CMint_attr_copy_list(CManager cm, attr_list l, char *file, int line)
{
    attr_list ret = attr_copy_list(l);
    (void)cm;
    CMtrace_out(cm, CMAttrVerbose, "Copy attr list %p at %s:%d, new list %p\n", 
		l, file, line, ret);
    return ret;
}

extern void
CMint_attr_merge_lists(CManager cm, attr_list l1, attr_list l2, 
		       char *file, int line)
{
    (void)cm;
    (void)file;
    (void)line;
    attr_merge_lists(l1, l2);
}

extern attr_list 
CMint_decode_attr_from_xmit(CManager cm, void * buf, char *file, int line)
{
    attr_list l = decode_attr_from_xmit(buf);
    (void)cm;
    CMtrace_out(cm, CMAttrVerbose, "decode attr list from xmit at %s:%d, new list %p\n", 
		file, line, l);
    return l;
}
#undef realloc
#undef malloc

extern void*
INT_CMrealloc(void *ptr, size_t size)
{
    void *tmp = realloc(ptr, size);
    if ((tmp == 0) && (size != 0)) {
	printf("Realloc failed on ptr %p, size %zd\n", ptr, size);
	perror("realloc");
    }
    return tmp;
}

extern void*
INT_CMmalloc(size_t size)
{
    void* tmp = malloc(size);
    if ((tmp == 0) && (size != 0)) {
	printf("Malloc failed on size %zd\n", size);
	perror("malloc");
    }
    return tmp;
}

extern void
INT_CMfree(void *ptr)
{
    free(ptr);
}

