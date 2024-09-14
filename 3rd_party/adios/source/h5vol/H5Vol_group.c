#ifndef __H5VOL_GROUP_FUNC
#define __H5VOL_GROUP_FUNC

#include "H5Vol_def.h"

void *H5VL_adios2_group_create(void *obj, const H5VL_loc_params_t *loc_params, const char *name,
                               hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id, hid_t dxpl_id,
                               void **req)
{
    REQUIRE_NOT_NULL_ERR(obj, NULL);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    // NOTE: flaky here as I do not check whether this group was create before.
    //       good faith on users
    if ((H5I_GROUP == loc_params->obj_type) || (H5I_FILE == loc_params->obj_type))
    {
        H5VL_GroupDef_t *grp = gCreateGroupDef(name);
        return gGroupToVolObj(grp, vol);
    }
    return NULL;
}

herr_t H5VL_adios2_group_close(void *obj, hid_t dxpl_id, void **req)
{
    if (NULL == obj)
        return 0;

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;
    H5VL_GroupDef_t *grp = (H5VL_GroupDef_t *)(vol->m_ObjPtr);
    SAFE_FREE(grp->m_Name);

    SAFE_FREE(grp);

    gFreeVol(vol);
    return 0;
}

void *H5VL_adios2_group_open(void *obj, const H5VL_loc_params_t *loc_params, const char *name,
                             hid_t gapl_id, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(obj, NULL);
    REQUIRE_NOT_NULL_ERR(loc_params, NULL);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    if ((H5I_GROUP == loc_params->obj_type) | (H5I_FILE == loc_params->obj_type))
    {
        H5VL_GroupDef_t *grp = gCreateGroupDef(name);
        return gGroupToVolObj(grp, vol);
    }

    return NULL;
}

herr_t H5VL_adios2_group_get(void *obj, H5VL_group_get_args_t *args, hid_t H5_ATTR_UNUSED dxpl_id,
                             void H5_ATTR_UNUSED **req)
{
    REQUIRE_NOT_NULL_ERR(obj, -1);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    switch (args->op_type)
    {
    case H5VL_GROUP_GET_INFO: {
        const H5VL_loc_params_t *loc_params = &args->args.get_info.loc_params;
        H5G_info_t *group_info = args->args.get_info.ginfo;

        if (loc_params->type == H5VL_OBJECT_BY_SELF)
        {
            gLoadContent(vol);
            gLoadSubGroups(vol);
            group_info->storage_type = H5G_STORAGE_TYPE_COMPACT;
            // group_info->nlinks       = vol->m_NumAttrs + vol->m_NumVars +
            // vol->m_NumSubGroups;
            group_info->nlinks = vol->m_NumVars + vol->m_NumSubGroups;
            return 0;
        }
    }
    default:
        break;
    }
    return -1;
}

#endif
