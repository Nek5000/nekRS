#include "config.h"
#include <stdio.h>
#include <string.h>
#undef NDEBUG
#include <assert.h>
#ifdef HAVE_NETDB_H
#include <netdb.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdlib.h>
#ifdef HAVE_WINDOWS_H
#define FD_SETSIZE 1024
#include <winsock2.h>
#define __ANSI_CPP__
#else
#include <netinet/in.h>
#include <arpa/inet.h>
#endif
#include <ffs.h>
#include <atl.h>
#include "evpath.h"
#include "cm_internal.h"
#include "cm_transport.h"
#ifdef HAVE_COD_H
#include "cod.h"
#endif

/* 
 * Creates the context of conversion for the incoming 'format',
 * caches it and returns a pointer to the cached context
 * */
extern CMincoming_format_list
CMidentify_rollbackCMformat(CManager cm, char *data_buffer)
{
    int i;
    int nearest_format = -1;
    FMcompat_formats older_format = NULL;
    FFSTypeHandle format;

    for (i = 0; i < cm->reg_format_count; i++) {
	if (cm->reg_formats[i]->registration_pending)
	    CMcomplete_format_registration(cm->reg_formats[i], 0);

    }
    format = FFS_target_from_encode(cm->FFScontext, data_buffer);

    cm->in_formats = INT_CMrealloc(cm->in_formats,
			       sizeof(struct _CMincoming_format) *
			       (cm->in_format_count + 1));
    cm->in_formats[cm->in_format_count].format = format;
    cm->in_formats[cm->in_format_count].handler =
	cm->reg_formats[nearest_format]->handler;
    cm->in_formats[cm->in_format_count].client_data =
	cm->reg_formats[nearest_format]->client_data;
    cm->in_formats[cm->in_format_count].older_format = older_format;
    cm->in_formats[cm->in_format_count].f2_format =
	cm->reg_formats[nearest_format];
    cm->in_formats[cm->in_format_count].f1_struct_size = 0;
    cm->in_formats[cm->in_format_count].code = NULL;
    cm->in_formats[cm->in_format_count].local_iocontext = NULL;
    return &cm->in_formats[cm->in_format_count++];
}

/* 
 * Creates the conversion context as specified by "cm_format".
 * If no rollback to older format is required, set the normal conversion
 * using set_conversion_FMcontext, else use FMlocalize_conv and 
 * FMlocalize_register_conv pair to localize and register wire format 
 * and set conversion to native format.
 * Then invoke gen_rollback_code to inorder to generate the rollback 
 * conversion code.
 * */
extern void
CMcreate_conversion(CManager cm, CMincoming_format_list cm_format)
{
    FMStructDescList native_format_list = cm_format->f2_format->format_list;
    FFSTypeHandle format = cm_format->format;
    FFSContext context = cm->FFScontext;
    FMcompat_formats older_format = cm_format->older_format;

    if (!older_format) {
	establish_conversion(context, format,
			     native_format_list);
	return;
    }
#ifdef EVOL
    FMStructDescList local_f0_formats, local_f1_formats;
    local_f0_formats = IOlocalize_conv(context, format);
    local_f1_formats = IOlocalize_register_conv(cm_format->local_iocontext,
						older_format->prior_format,
						native_field_list,
						native_subformat_list,
						&cm_format->
						local_prior_format,
						&cm_format->
						f1_struct_size);
    cm_format->code =
	gen_rollback_code(local_f0_formats, local_f1_formats,
			  older_format->xform_code);
    free_FMStructDescList(local_f0_formats);
    free_FMStructDescList(local_f1_formats);
#endif
}

/* 
 * Read the cached context and do the conversion of data
 * The conversion take place in several steps:
 * 1. Convert the localized data (*decode_buff) to a different format 
 *              that the client can understand using the generated ecode
 * 2. Do pbio encoding of the buffer returned by ecode funtion
 * 3. Decode the above encoded buffer using the conversion set in 
 *              build_conversion
 * 
 * Returns 1 on success, 0 otherwise
 * '*decode_buff' is appropriately set to the decoded buffer which can
 * now be passed to the registered handler.
 * '*cm_decode_buffer' is set so that the above buffer can be freed when the 
 * application handler returns.
 * */
#ifdef EVOL
extern int
process_old_format_data(CManager cm, CMincoming_format_list cm_format,
			char **decode_buff, CMbuffer * cm_decode_buffer)
{
    int encoded_buf_size;
    char *encoded_buf, *f1_buffer, *decode_buffer = *decode_buff;
    CMbuffer cm_f1_buf, cm_decode_buf = *cm_decode_buffer;
    IOFormat f1_format = cm_format->local_prior_format;
    IOContext local_iocontext = cm_format->local_iocontext;
    int (*func)();

    assert(cm_format->code && cm_format->f1_struct_size);
    cm_f1_buf = cm_get_data_buf(cm, cm_format->f1_struct_size);
    f1_buffer = cm_f1_buf->buffer;

    memset(f1_buffer, 0, cm_format->f1_struct_size);
    func = (int(*)())cm_format->code->func;
    func(decode_buffer, f1_buffer);

    assert(has_conversion_IOformat(f1_format));

    encoded_buf = encode_IOcontext_buffer(local_iocontext, f1_format,
					  f1_buffer, &encoded_buf_size);
    if (cm_f1_buf) {
	/* TODO: Need to free ecl f1_buffer IOfree_var_rec_elements
	 * segment faults */
	// IOfree_var_rec_elements((IOFile)local_iocontext, f1_format,
	// f1_buffer);
	cm_return_data_buf(cm, cm_f1_buf);
	cm_f1_buf = NULL;
    }
    if (cm_decode_buf) {
	cm_return_data_buf(cm, cm_decode_buf);
	cm_decode_buf = NULL;
    }
    if (decode_in_place_possible(f1_format)) {
	if (!decode_in_place_IOcontext(local_iocontext, encoded_buf,
				       (void **) (long) &decode_buffer)) {
	    printf("Decode failed at line:%d\n", __LINE__);
	    {
		CManager cm = NULL;
		printf("Stopped at line:%d\n", __LINE__);
		printf("%p\n", cm->IOcontext);
	    }
	    return 0;		/* FAILED */
	}
    } else {
	int decoded_length =
	    this_IOrecord_length(local_iocontext, encoded_buf,
				 encoded_buf_size);
	cm_decode_buf = cm_get_data_buf(cm, decoded_length);
	decode_buffer = cm_decode_buf->buffer;
	decode_to_buffer_IOcontext(local_iocontext, encoded_buf,
				   decode_buffer);
    }
    *cm_decode_buffer = cm_decode_buf;
    *decode_buff = decode_buffer;
    return 1;			/* SUCCESS */
}
#endif
