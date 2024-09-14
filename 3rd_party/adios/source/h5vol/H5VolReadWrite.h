#ifndef ADIOS_VOL_WRITER_H
#define ADIOS_VOL_WRITER_H

#define H5F_FRIEND /*suppress error about including H5Fpkg   */

#include "H5VolError.h"
#include "H5VolUtil.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <adios2_c.h>

#include "H5Vol_def.h"

//
// VL definition
//
static herr_t H5VL_adios2_init(hid_t vipl_id) { return 0; }

static herr_t H5VL_adios2_term(void)
{
    gExitADIOS2();
    return 0;
}

extern herr_t H5VL_adios2_begin_read_step(const char *);
extern herr_t H5VL_adios2_begin_write_step(const char *);

extern herr_t H5VL_adios2_beginstep(const char *engine_name, adios2_step_mode m);

extern herr_t H5VL_adios2_endstep(const char *engine_nane);

static herr_t H5VL_adios2_introspect_get_cap_flags(const void *info, uint64_t *cap_flags);

static herr_t H5VL_adios2_introspect_opt_query(void *obj, H5VL_subclass_t cls, int opt_type,
                                               uint64_t *supported)
{
    *supported = 0;
    return 0;
}

static herr_t H5VL_adios2_datatype_close(void *dt, hid_t H5_ATTR_UNUSED dxpl_id,
                                         void H5_ATTR_UNUSED **req)
{
    return 0;
}

//
// Define H5VOL functions
//
static const H5VL_class_t H5VL_adios2_def = {
    H5VL_VERSION, /* Version # of connector, needed for v1.13 */
    (H5VL_class_value_t)H5VL_ADIOS2_VALUE,
    H5VL_ADIOS2_NAME,      /* name */
    H5VL_ADIOS2_VERSION,   /* version of this vol, not as important  */
    H5VL_ADIOS2_CAP_FLAGS, /* Capability flags for connector */
    H5VL_adios2_init,      /* initialize */
    H5VL_adios2_term,      /* terminate */
    {
        /* info_cls */
        (size_t)0, /* info size    */
        NULL,      /* info copy    */
        NULL,      /* info compare */
        NULL,      /* info free    */
        NULL,      /* info to str  */
        NULL       /* str to info  */
    },
    {
        NULL, /* get_object   */
        NULL, /* get_wrap_ctx */
        NULL, /* wrap_object  */
        NULL, /* unwrap_object */
        NULL  /* free_wrap_ctx */
    },
    {H5VL_adios2_attr_create, H5VL_adios2_attr_open, H5VL_adios2_attr_read, H5VL_adios2_attr_write,
     H5VL_adios2_attr_get, H5VL_adios2_attr_specific,
     NULL, // H5VL_adios2_attr_optional,
     H5VL_adios2_attr_close},
    {
        /* dataset_cls */
        H5VL_adios2_dataset_create, H5VL_adios2_dataset_open, H5VL_adios2_dataset_read,
        H5VL_adios2_dataset_write, H5VL_adios2_dataset_get, /* get properties*/
        NULL,                                               // H5VL_adios2_dataset_specific
        NULL,                                               // optional
        H5VL_adios2_dataset_close                           /* close */
    },
    {
        /* datatype_cls */
        NULL,                      // H5VL_adios2_datatype_commit,           /* commit */
        NULL,                      // H5VL_adios2_datatype_open,             /* open */
        NULL,                      // H5VL_adios2_datatype_get, /* get_size */
        NULL,                      // H5VL_adios2_datatype_specific,
        NULL,                      // H5VL_adios2_datatype_optional,
        H5VL_adios2_datatype_close /* close */
    },
    {H5VL_adios2_file_create, H5VL_adios2_file_open,
     NULL, // H5VL_adios2_file_get,
     H5VL_adios2_file_specific,
     NULL, // H5VL_adios2_file_optional,
     H5VL_adios2_file_close},
    {H5VL_adios2_group_create, H5VL_adios2_group_open, H5VL_adios2_group_get, NULL, NULL,
     H5VL_adios2_group_close},
    {
        NULL, // H5VL_adios2_link_create,
        NULL, // H5VL_adios2_link_copy,
        NULL, // H5VL_adios2_link_move,
        H5VL_adios2_link_get, H5VL_adios2_link_specific,
        NULL // H5VL_adios2_link_remove
    },
    {
        H5VL_adios2_object_open,
        NULL, // H5VL_adios2_object_copy,
        H5VL_adios2_object_get,
        NULL, // H5VL_adios2_object_specific,
        NULL  // H5VL_adios2_object_optional,
    },
    {
        /* introspect_cls */
        NULL, // H5VL_pass_through_introspect_get_conn_cls,  /* get_conn_cls */
        H5VL_adios2_introspect_get_cap_flags, /* get_cap_flags */
        H5VL_adios2_introspect_opt_query,     /* opt_query */
    },
    {
        /* request_cls */
        NULL, /* wait */
        NULL, /* notify */
        NULL, /* cancel */
        NULL, /* specific */
        NULL, /* optional */
        NULL  /* free */
    },
    {
        /* blob_cls */
        NULL, /* put */
        NULL, /* get */
        NULL, /* specific */
        NULL  /* optional */
    },
    {
        /* token_cls */
        NULL, /* cmp */
        NULL, /* to_str */
        NULL  /* from_str */
    },
    NULL /*/optional*/
};

static herr_t H5VL_adios2_introspect_get_cap_flags(const void *info, uint64_t *cap_flags)
{
    *cap_flags = H5VL_adios2_def.cap_flags;
    return 0;
}

#endif // ADIOS_VOL_WRITER_H
