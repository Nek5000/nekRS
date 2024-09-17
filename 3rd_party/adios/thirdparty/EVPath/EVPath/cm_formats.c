#include "config.h"

#ifndef MODULE
#include <stdio.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>

#else
#include "kernel/kcm.h"
#include "kernel/cm_kernel.h"
#include "kernel/library.h"
#undef NDEBUG
#include "assert.h"

#endif
#include "ffs.h"
#include "atl.h"
#include "evpath.h"
#include "chr_time.h"
#include "cm_internal.h"

static CMFormat add_format_to_cm(CManager cm, CMFormat format);

CMFormat
INT_CMlookup_format(CManager cm, FMStructDescList format_list)
{
    int i;
    for (i=0; i< cm->reg_format_count; i++) {
	if (cm->reg_formats[i]->format_list == format_list) {
	    return cm->reg_formats[i];
	}
    }
    return NULL;
}

CMFormat
INT_CMregister_simple_format(CManager cm, char *format_name, FMFieldList field_list, int struct_size)
{
  FMStructDescRec *format_list = malloc(sizeof(*format_list) * 2);
  /* NOTE:  This is a memory leak as we do not, by design, deallocate the format_list */
  format_list[0].format_name = format_name;
  format_list[0].field_list = field_list;
  format_list[0].struct_size = struct_size;
  format_list[0].opt_info = NULL;
  format_list[1].format_name = NULL;
  format_list[1].field_list = NULL;
  format_list[1].struct_size = 0;
  format_list[1].opt_info = NULL;
  return INT_CMregister_format(cm, &format_list[0]);
}

CMFormat
INT_CMregister_format(CManager cm, FMStructDescList format_list)
{
    CMFormat format;

    if ((format_list == NULL) || (cm == NULL)) 
	return NULL;

    format = INT_CMmalloc(sizeof(struct _CMFormat));
    
    format->cm = cm;
    format->format_name = INT_CMmalloc(strlen(format_list[0].format_name) + 1);
    strcpy(format->format_name, format_list[0].format_name);
    format->fmformat = NULL;
    format->format_list_addr = format_list;
    format->handler = (CMHandlerFunc) NULL;
    format->client_data = NULL;
    format->format_list = format_list;
    format->registration_pending = 1;

    format = add_format_to_cm(cm, format);
    return format;
}

#ifdef COMPAT
void *INT_CMcreate_compat_info(format, xform_code, len_p)
CMFormat format;
char *xform_code;
int *len_p;
{
    if(!format) return NULL;
    if(!format->format) CMcomplete_format_registration(format, 0);
    if(!format->format) return NULL;
    return create_compat_info(format->format, xform_code, len_p);
}
#endif

extern void 
CMcomplete_format_registration(CMFormat format, int lock)
{
    
    CManager cm = format->cm;
    FMContext c = FMContext_from_FFS(format->cm->FFScontext);
    int i;
    format->fmformat = register_data_format(c, format->format_list);
    format->ffsformat = FFSset_fixed_target(format->cm->FFScontext, 
					    format->format_list);
	    
    for(i=0; i < cm->in_format_count; i++) {
	if (cm->in_formats[i].format == format->ffsformat) {
	    format->fmformat = NULL;
	    return;
	}
    }
    cm->in_formats = INT_CMrealloc(cm->in_formats,
			       sizeof(struct _CMincoming_format) *
			       (cm->in_format_count + 1));
    cm->in_formats[cm->in_format_count].format = format->ffsformat;
    cm->in_formats[cm->in_format_count].handler = format->handler;
    cm->in_formats[cm->in_format_count].client_data = format->client_data;
    cm->in_formats[cm->in_format_count].older_format = NULL;
    cm->in_formats[cm->in_format_count].f2_format = format;
    cm->in_formats[cm->in_format_count].f1_struct_size = 0;
    cm->in_formats[cm->in_format_count].code = NULL;
    cm->in_formats[cm->in_format_count].local_iocontext = NULL;
    cm->in_format_count++;

    if (format->fmformat == NULL) {
	fprintf(stderr, "Format registration failed for format \"%s\"\n",
		format->format_name);
        free_CMFormat (format);
	if (lock) CManager_unlock(format->cm);
	return;
    }
    format->registration_pending = 0;
}

static CMFormat
add_format_to_cm(CManager cm, CMFormat format)
{
    char *format_name = format->format_name;
    int insert_before = 0, i;

    i = 0;
    for (i=0; i< cm->reg_format_count; i++) {
	int order = strcmp(format_name, cm->reg_formats[i]->format_name);
	if (order < 0) {
	    insert_before = i;
	    break;
	} else if (order == 0) {
	    /* have the same name */
	    FMformat_order suborder;
	    if (cm->reg_formats[i]->registration_pending) {
		CMcomplete_format_registration(cm->reg_formats[i], 0);
	    }
	    if (format->registration_pending) {
		CMcomplete_format_registration(format, 0);
		if (format->registration_pending) {
		    int i = 0;
		    /* if still pending, it's a duplicate format */
		    for(i=0; i < cm->in_format_count; i++) {
			if (cm->in_formats[i].format == format->ffsformat) {
			    free(format->format_name);
			    free(format);
			    return cm->in_formats[i].f2_format;
			}
		    }
		    printf("Gack, duplicate format, but didn't find it\n");
		    return NULL;
		}
	    }
	    suborder = FMformat_cmp(format->fmformat, 
				    cm->reg_formats[i]->fmformat);
	    if ((suborder == Format_Greater) || 
		(suborder == Format_Incompatible)) {
		insert_before = i;
		break;
	    } else if (suborder == Format_Equal) {
                insert_before = i;
	    }
	}
    }
    if (i == cm->reg_format_count) {
	insert_before = i;
    }
    cm->reg_formats = INT_CMrealloc(cm->reg_formats, sizeof(CMFormat) * 
				(cm->reg_format_count + 1));
    for (i = cm->reg_format_count; i > insert_before; i--) {
	/* move this up */
	cm->reg_formats[i] = cm->reg_formats[i-1];
    }
    cm->reg_formats[insert_before] = format;
    cm->reg_format_count++;
    return format;
}

extern void
free_CMFormat(CMFormat format)
{
    INT_CMfree(format);
}

extern CMincoming_format_list
CMidentify_CMformat(CManager cm, FFSTypeHandle format)
{
    int i;
    char *format_name = name_of_FMformat(FMFormat_of_original(format));
    FMStructDescList native_format_list;

    for (i=0; i< cm->reg_format_count; i++) {
	int order = strcmp(format_name, cm->reg_formats[i]->format_name);
	if (order < 0) {
	    return NULL;
	} else if (order == 0) {
	    /* 
	     *  we found a registered format with the same name as the 
	     *  incoming record, is it compatible? 
	     */
	    if (cm->reg_formats[i]->registration_pending) {
		CMcomplete_format_registration(cm->reg_formats[i], 0);
	    }
	    if (cm->reg_formats[i]->fmformat == NULL) {
		continue;
	    }
	    switch(FMformat_cmp(FMFormat_of_original(format), 
				cm->reg_formats[i]->fmformat)) {
	    case Format_Equal:
	    case Format_Greater:
		/* 
		 * if the incoming format has the same or more fields
		 * as the registered format, we're cool
		 */
		break;
	    case Format_Less:
	    case Format_Incompatible:
		/* 
		 * if the incoming format has fewer fields than the 
		 * registered handler requires or they're incompatible
		 * keep looking;
		 */
		 continue;
	    }
	    break;
	}
    }

    if (i >= cm->reg_format_count) return NULL;
    native_format_list = cm->reg_formats[i]->format_list;

    establish_conversion(cm->FFScontext, format, native_format_list);
    cm->in_formats = INT_CMrealloc(cm->in_formats, 
			       sizeof(struct _CMincoming_format) * 
			       (cm->in_format_count + 1));
    cm->in_formats[cm->in_format_count].format = format;
    cm->in_formats[cm->in_format_count].handler = 
	cm->reg_formats[i]->handler;
    cm->in_formats[cm->in_format_count].client_data = 
	cm->reg_formats[i]->client_data;
    return &cm->in_formats[cm->in_format_count++];
}

extern FFSTypeHandle
INT_CMget_format_IOcontext(CManager cm, FFSContext context, void *buffer)
{
    FFSTypeHandle ret;
    (void)cm;
    ret = FFSTypeHandle_from_encode(context, buffer);
    return ret;
}

extern FFSTypeHandle
INT_CMget_format_app_IOcontext(CManager cm, FFSContext context, void *buffer,
			       void *app_context)
{
    FFSTypeHandle ret;
    (void) cm;
    (void) app_context;
    ret = FFSTypeHandle_from_encode(context, buffer);
    return ret;
}

extern void
INT_CMset_conversion_IOcontext(CManager cm, FFSContext context,
			       FFSTypeHandle format,
			       FMStructDescList format_list)
{
    (void)cm;
    establish_conversion(context, format, format_list);
}

extern int CMself_hosted_formats;

static void
preload_pbio_format(CMConnection conn, FMFormat ioformat)
{
    CMtrace_out(conn->cm, CMFormatVerbose, 
		"CMpbio preloading format %s on connection %p\n", 
		name_of_FMformat(ioformat), conn);
    if (CMpbio_send_format_preload(ioformat, conn) != 1) {
	if (!conn->closed) 
	    CMtrace_out(conn->cm, CMFormatVerbose, "CMpbio preload failed\n");
	return;
    }
#ifndef MODULE
    if (CMtrace_on(conn->cm, CMFormatVerbose)) {
	int junk;
	fprintf(conn->cm->CMTrace_file, "CMpbio Preload is format ");
	fprint_server_ID(conn->cm->CMTrace_file, (unsigned char *)get_server_ID_FMformat(ioformat, &junk));
	fprintf(conn->cm->CMTrace_file, "\n");
    }
#endif
}

extern void
CMformat_preload(CMConnection conn, CMFormat format)
{
    int load_count = 0;
    FMFormat *loaded_list = conn->preloaded_formats;
    int my_FFSserver_ID = conn->cm->FFSserver_identifier;
    int remote_FFSserver_ID = conn->remote_format_server_ID;
    int preload = 0;

    if (my_FFSserver_ID == -1) preload = 1;   /* we're self hosting formats */
    if (remote_FFSserver_ID == -1) preload = 1;   /* they're self hosting formats */
    if (remote_FFSserver_ID == 0) preload = 1; /* 0 is the unset state, we don't know their server, preload to be safe */
    if (remote_FFSserver_ID != my_FFSserver_ID) preload = 1;  /* if we're both using a format server, but not the same one, preload */

    if (preload == 0) return;

    if (conn->closed) return;

    /* if (CMtrace_on(conn->cm, CMFormatVerbose)) { */
    /* 	int junk; */
    /* 	fprintf(conn->cm->CMTrace_file, "Considering CMpbio Preload of format "); */
    /* 	fprint_server_ID(conn->cm->CMTrace_file, (unsigned char *)get_server_ID_FMformat(format->fmformat, &junk)); */
    /* 	fprintf(conn->cm->CMTrace_file, " value %p\n", format); */
    /* } */
    while (loaded_list && (*loaded_list != NULL)) {
      if (*loaded_list == format->fmformat) return;
	loaded_list++;
	load_count++;
    }
    
    preload_pbio_format(conn, format->fmformat);

    if (conn->preloaded_formats == NULL) {
	loaded_list = malloc(2*sizeof(*loaded_list));
    } else {
	loaded_list = realloc(conn->preloaded_formats, 
			      sizeof(*loaded_list) * (load_count + 2));
    }
    loaded_list[load_count] = format->fmformat;
    loaded_list[load_count+1] = NULL;
    conn->preloaded_formats = loaded_list;
}
