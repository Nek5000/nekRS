#ifndef __H5VOL_OBJECT_FUNC
#define __H5VOL_OBJECT_FUNC

#include "H5Vol_def.h"

void *H5VL_adios2_object_open(void *obj, const H5VL_loc_params_t *loc_params,
                              H5I_type_t *opened_type, hid_t H5_ATTR_UNUSED dxpl_id,
                              void H5_ATTR_UNUSED **req)
{
    REQUIRE_NOT_NULL_ERR(loc_params, NULL);
    REQUIRE_NOT_NULL_ERR(obj, NULL);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    // can not open anything under attr
    if (ATTR == vol->m_ObjType)
        return NULL;

    switch (loc_params->type)
    {
    case H5VL_OBJECT_BY_NAME: {
        const char *obj_name = loc_params->loc_data.loc_by_name.name;
        adios2_attribute *attr = gLocateAttrFrom(vol, obj_name);

        if (NULL != attr)
        {
            H5VL_AttrDef_t *attrDef = gCreateAttrDef(obj_name, -1, -1);
            attrDef->m_Attribute = attr;
            gLoadAttrDef(attrDef);
            *opened_type = H5I_ATTR;
            return gAttrToVolObj(attrDef, vol);
        }

        H5VL_ObjDef_t *varObj = gGetVarObjDef(obj_name, vol);
        if (NULL != varObj)
        {
            *opened_type = H5I_DATASET;
            gLoadContent(varObj);
            return varObj;
        }

        // Good faith: assume to be a valid GROUP
        H5VL_GroupDef_t *group = gCreateGroupDef(obj_name);
        H5VL_ObjDef_t *result = (H5VL_ObjDef_t *)gGroupToVolObj(group, vol);

        gLoadContent(result);
        *opened_type = H5I_GROUP;
        return result;
    }
    default:
        break;
    }

    return NULL;
}

herr_t H5VL_adios2_object_get(void *obj, const H5VL_loc_params_t *loc_params,
                              H5VL_object_get_args_t *args, hid_t H5_ATTR_UNUSED dxpl_id,
                              void H5_ATTR_UNUSED **req)
{
    REQUIRE_NOT_NULL_ERR(loc_params, -1);
    REQUIRE_NOT_NULL_ERR(obj, -1);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    switch (args->op_type)
    {
    case H5VL_OBJECT_GET_INFO: {
        H5O_info2_t *oinfo = args->args.get_info.oinfo;
        if (loc_params->type == H5VL_OBJECT_BY_SELF)
        {
            oinfo->fileno = 1;
            oinfo->num_attrs = vol->m_NumAttrs;
            if (GROUP == vol->m_ObjType)
                oinfo->type = H5O_TYPE_GROUP;
            else if (VAR == vol->m_ObjType)
                oinfo->type = H5O_TYPE_DATASET;
            else if (ATTR == vol->m_ObjType)
                oinfo->type = H5O_TYPE_UNKNOWN; // no attr type
            else
                oinfo->type = H5O_TYPE_GROUP; // no file type

            return 0;
        }
        else if (loc_params->type == H5VL_OBJECT_BY_IDX)
        {
            oinfo->fileno = 1;
            int idx = loc_params->loc_data.loc_by_idx.n;
            if ((GROUP == vol->m_ObjType) || (ROOT == vol->m_ObjType))
            {
                gLoadContent(vol);
                gLoadSubGroups(vol);

                if (idx < (vol->m_NumVars))
                    oinfo->type = H5O_TYPE_DATASET;
                else
                    oinfo->type = H5O_TYPE_GROUP;
                return 0;
            }
        }
        // We do not support Get_INFO_BY_NAME or BY_IDX
    }

    case H5VL_OBJECT_GET_NAME:
    case H5VL_OBJECT_GET_TYPE:
    case H5VL_OBJECT_GET_FILE:
    default:
        break;
    }
    return -1;
}

#endif
