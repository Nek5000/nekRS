#include "H5Vol_def.h"

//
// NOTE: this is called from H5F.c when a new file or  trunc file is asked
//       so no need to check flags here. if do need to, use & not ==
//
void *H5VL_adios2_file_create(const char *name, unsigned flags, hid_t fcpl_id, hid_t fapl_id,
                              hid_t dxpl_id, void **req)
{
    gInitADIOS2(fapl_id);
    if (flags & H5F_ACC_TRUNC)
    {
        H5VL_FileDef_t *handle = gADIOS2CreateFile(name);
        return gFileToVolObj(handle);
    }

    if (flags & H5F_ACC_EXCL)
    {
        H5VL_FileDef_t *handle = gADIOS2OpenFile(name);
        if (NULL == handle->m_Engine)
        {
            H5VL_FileDef_t *handle = gADIOS2CreateFile(name);
            return gFileToVolObj(handle);
        }
        else
        {
            gADIOS2CloseFile(handle); // f is freed in this call
        }
        // exists,
    }
    return NULL;
}

void *H5VL_adios2_file_open(const char *name, unsigned flags, hid_t fapl_id, hid_t dxpl_id,
                            void **req)
{
    gInitADIOS2(fapl_id);
    H5VL_FileDef_t *handle = gADIOS2OpenFile(name);
    return gFileToVolObj(handle);
}

herr_t H5VL_adios2_file_specific(void *file, H5VL_file_specific_args_t *args, hid_t dxpl_id,
                                 void **req)
{
    //
    // This function is called after H5Fopen/create. Do not remove
    //
    return 0;
}

herr_t H5VL_adios2_file_close(void *file, hid_t dxpl_id, void **req)
{
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)file;
    REQUIRE_SUCC((ROOT == vol->m_ObjType), -1);
    H5VL_FileDef_t *f = (H5VL_FileDef_t *)(vol->m_ObjPtr);
    gADIOS2CloseFile(f); // f is freed in this call
    gFreeVol(vol);
    return 0;
}
