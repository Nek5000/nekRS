#ifndef __H5VOL_ATTR_FUNC
#define __H5VOL_ATTR_FUNC

#include "H5Vol_def.h"

///
// attribute handlers
//
void *H5VL_adios2_attr_create(void *obj, const H5VL_loc_params_t *loc_params, const char *name,
                              hid_t type_id, hid_t space_id, hid_t acpl_id, hid_t aapl_id,
                              hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(loc_params, NULL);
    REQUIRE_NOT_NULL_ERR(obj, NULL);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    H5VL_AttrDef_t *attrDef = gCreateAttrDef(name, type_id, space_id);
    return gAttrToVolObj(attrDef, vol);
}

void *H5VL_adios2_attr_open(void *obj, const H5VL_loc_params_t *loc_params, const char *name,
                            hid_t aapl_id, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(loc_params, NULL);
    REQUIRE_NOT_NULL_ERR(obj, NULL);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;
    adios2_attribute *attr = gLocateAttrFrom(vol, name);

    H5VL_AttrDef_t *attrDef = NULL;
    if (NULL == attr)
    {
        if ('/' == name[0])
        {
            SHOW_ERROR_MSG("H5VL_ADIOS2: Error: No such ATTRIBUTE: [%s] in file\n ", name);
            return NULL;
        }
        else
        {
            char withSlash[strlen(name) + 2];
            snprintf(withSlash, sizeof(withSlash), "/%s", name);
            withSlash[strlen(name) + 1] = '\0';
            attr = gLocateAttrFrom(vol, withSlash);
            if (NULL == attr)
            {
                SHOW_ERROR_MSG("H5VL_ADIOS2: Error: No such ATTRIBUTE: [%s] "
                               "found in file\n ",
                               withSlash);
                return NULL;
            }
            attrDef = gCreateAttrDef(withSlash, -1, -1);
        }
    }
    else
        attrDef = gCreateAttrDef(name, -1, -1);

    // H5VL_AttrDef_t *attrDef = gCreateAttrDef(name, -1, -1);
    attrDef->m_Attribute = attr;

    gLoadAttrDef(attrDef);

    return gAttrToVolObj(attrDef, vol);
}

herr_t H5VL_adios2_attr_read(void *attrObj, hid_t mem_type_id, void *buf, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(attrObj, -1);

    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)attrObj;

    H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)(vol->m_ObjPtr);
    adios2_attribute *attr = attrDef->m_Attribute;

    if (NULL == attr)
        return -1;

    if ((!attrDef->m_IsScalar) && (H5Tget_class(mem_type_id) == H5T_STRING))
    {
        htri_t isVariableSize = H5Tis_variable_str(mem_type_id);
        if (!isVariableSize)
        {
            size_t strSize = H5Tget_size(mem_type_id);
            // buf is char* instead of char**, adios2 c api expects char**

            char **result = (char **)(malloc(sizeof(char *) * (int)(attrDef->m_Size)));
            size_t k = 0;
            for (k = 0; k < attrDef->m_Size; k++)
                result[k] = (char *)(malloc(strSize));

            adios2_attribute_data(result, &(attrDef->m_Size), attr);

            char *out = (char *)buf;

            if (strlen(out) > 0)
            {
                buf = result;
                return 0;
            }
            // from openPMD api.. using char*[n*max_length] to hold result
            for (k = 0; k < attrDef->m_Size; k++)
            {
                strncpy(out + k * strSize, result[k], strSize);
                out[k * strSize + strlen(result[k])] = '\0';
                free(result[k]);
            }

            free(result);
            return 0;
        }
    }

    adios2_attribute_data(buf, &(attrDef->m_Size), attr);

    return 0;
}

herr_t H5VL_adios2_attr_write(void *attr, hid_t mem_type_id, const void *buf, hid_t dxpl_id,
                              void **req)
{
    REQUIRE_NOT_NULL_ERR(attr, -1);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)attr;

    H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)(vol->m_ObjPtr);
    attrDef->m_Data = (void *)buf;
    gADIOS2CreateAttr(vol->m_FileIO, attrDef, vol->m_Path);

    return 0;
}

void GetFromAttribute(void *attrObj, hid_t *ret_id, H5VL_attr_get_t get_type)
{
    H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)attrObj;

    switch (get_type)
    {
    case H5VL_ATTR_GET_SPACE: {
        *ret_id = H5Scopy(attrDef->m_SpaceID);
        break;
    }
    case H5VL_ATTR_GET_TYPE: {
        *ret_id = H5Tcopy(attrDef->m_TypeID);
        break;
    }
    default:
        break;
    }
}

herr_t H5VL_adios2_attr_get(void *obj, H5VL_attr_get_args_t *args, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(obj, -1);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    switch (args->op_type)
    {
    case H5VL_ATTR_GET_SPACE: {
        hid_t *ret_id = (hid_t *)args->args.get_space.space_id;
        GetFromAttribute((vol->m_ObjPtr), ret_id, args->op_type);
        return 0;
    }
    case H5VL_ATTR_GET_TYPE: {
        hid_t *ret_id = (hid_t *)args->args.get_type.type_id;
        GetFromAttribute((vol->m_ObjPtr), ret_id, args->op_type);
        return 0;
    }
    default:
        break;
    }

    switch (args->op_type)
    {
    case H5VL_ATTR_GET_NAME: {
        char *buf = args->args.get_name.buf;
        size_t *ret_val = (size_t *)args->args.get_name.attr_name_len;

        const H5VL_loc_params_t *loc_params = &args->args.get_name.loc_params;
        REQUIRE_NOT_NULL_ERR(loc_params, -1);

        if (H5VL_OBJECT_BY_SELF == loc_params->type)
        {
            H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)(vol->m_ObjPtr);
            *ret_val = strlen(attrDef->m_Name);
            if (buf)
            {
                memcpy(buf, attrDef->m_Name, *ret_val);
            }
        }
        else if (H5VL_OBJECT_BY_IDX == loc_params->type)
        {
            // The  number of attrs is from H5Oget_info(), then iterate each by
            // calling H5Aget_name_by_idx, to reach here
            *ret_val = gGetNameOfNthAttr(vol, loc_params->loc_data.loc_by_idx.n, buf);
        }
        return 0;
    }
    default:
        break;
    }

    ADIOS_VOL_NOT_SUPPORTED_ERR("UNABLE TO SUPPORT feature at attr_get\n");
    return -1;
}

herr_t H5VL_adios2_attr_close(void *attr, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(attr, -1);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)attr;

    H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)(vol->m_ObjPtr);
    free(attrDef->m_Name);
    if (-1 != attrDef->m_SpaceID)
        H5Sclose(attrDef->m_SpaceID);
    free(attrDef);
    attrDef = NULL;

    gFreeVol(vol);
    vol = NULL;
    return 0;
}

herr_t H5VL_adios2_attr_specific(void *obj, const H5VL_loc_params_t *loc_params,
                                 H5VL_attr_specific_args_t *args, hid_t dxpl_id, void **req)
{
    REQUIRE_NOT_NULL_ERR(obj, -1);
    H5VL_ObjDef_t *vol = (H5VL_ObjDef_t *)obj;

    const char *attr_name;
    adios2_attribute *attr;

    switch (args->op_type)
    {
    case H5VL_ATTR_DELETE: {
        attr_name = (const char *)args->args.del.name;
        attr = gLocateAttrFrom(vol, attr_name);

        if (NULL != attr)
        {
            if (NULL == vol->m_Path)
                gADIOS2RemoveAttr(vol->m_FileIO, attr_name);
            else
            {
                char fullPath[strlen(vol->m_Path) + 4 + strlen(attr_name)];
                gGenerateFullPath(fullPath, vol->m_Path, attr_name);
                gADIOS2RemoveAttr(vol->m_FileIO, fullPath);
            }
            return 0;
        }
    }
    case H5VL_ATTR_EXISTS: {
        hbool_t *ret = args->args.exists.exists;

        attr_name = (const char *)args->args.exists.name;
        attr = gLocateAttrFrom(vol, attr_name);

        if (NULL == attr)
        {
            *ret = 0;
        }
        else
        {
            *ret = 1;
        }

        return 0;
    }
    default:
        break;
    }
    return 0;
}

#endif
