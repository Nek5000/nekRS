#include "H5VLpublic.h"
#include "hdf5.h"
#define H5S_FRIEND // suppress error for H5Spkg
#define H5O_FRIEND // suppress error for H5Opkg

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include "H5Epublic.h"
#include "H5VolReadWrite.h"

#include "H5VolUtil.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h> // sleep
#elif defined HAVE_WINDOWS_H
#include <windows.h>
#define sleep(x) Sleep(1000 * (x))
#endif

// these are in h5private.h
#define SUCCEED 1
#define FAIL 0

static hid_t H5VL_ADIOS_g = -1;

static adios2_adios *m_ADIOS2 = NULL;
// static adios2_io *m_IO = NULL;
static int m_MPIRank = 0;

#define RANK_ZERO_MSG(...)                                                                         \
    {                                                                                              \
        if (0 == m_MPIRank)                                                                        \
        {                                                                                          \
            fprintf(stderr, "## VOL info:");                                                       \
            fprintf(stderr, __VA_ARGS__);                                                          \
            fflush(stderr);                                                                        \
        }                                                                                          \
    }

void gGenerateFullPath(char *fullPath, const char *parentPath, const char *name)
{
    size_t ps = strlen(parentPath);
    size_t ns = strlen(name);
    size_t length = ps;
    bool startsWithDotSlash = ((ns > 1) && ('/' == name[1]) && ('.' == name[0]));

    if ('/' == parentPath[ps - 1])
    {
        if (startsWithDotSlash)
        {
            sprintf(fullPath, "%s%s", parentPath, name + 2);
            length += ns - 2;
        }
        else
        {
            sprintf(fullPath, "%s%s", parentPath, name);
            length += ns;
        }
    }
    else
    {
        if (startsWithDotSlash)
        {
            sprintf(fullPath, "%s/%s", parentPath, name + 2);
            length += 1 + ns - 2;
        }
        else
        {
            sprintf(fullPath, "%s/%s", parentPath, name);
            length += 1 + ns;
        }
    }
    fullPath[length] = '\0';
}

herr_t H5VL_adios2_begin_read_step(const char *filename)
{
    return H5VL_adios2_beginstep(filename, adios2_step_mode_read);
}

herr_t H5VL_adios2_begin_write_step(const char *filename)
{
    return H5VL_adios2_beginstep(filename, adios2_step_mode_append);
}

herr_t H5VL_adios2_beginstep(const char *filename, adios2_step_mode m)
{
    adios2_io *m_IO = adios2_at_io(m_ADIOS2, filename);
    adios2_engine *currEngine = adios2_get_engine(m_IO, filename);
    if (NULL == currEngine)
        return -1;

    adios2_step_status status;
    adios2_begin_step(currEngine, m, 0.0, &status);

    if (adios2_step_status_end_of_stream == status)
    {
        RANK_ZERO_MSG("..end_of_stream \n");
        return -1;
    }
    else if (adios2_step_status_not_ready == status)
    {
        RANK_ZERO_MSG(".. not ready \n");
        while (adios2_step_status_not_ready == status)
        {
            sleep(1);
            adios2_begin_step(currEngine, m, 0.0, &status);
        }
        RANK_ZERO_MSG("... other status \n");
        if (adios2_step_status_ok == status)
        {
            return 0;
        }
        return -1;
    }
    else if (adios2_step_status_ok == status)
    {
        RANK_ZERO_MSG(".. stream ready \n");
        return 0;
    }
    return -1;
}

herr_t H5VL_adios2_endstep(const char *filename)
{
    adios2_io *m_IO = adios2_at_io(m_ADIOS2, filename);
    adios2_engine *engine = adios2_get_engine(m_IO, filename);
    if (NULL == engine)
        return -1;

    adios2_end_step(engine);
    return 0;
}

void gInitADIOS2(hid_t acc_tpl)
{
    if (NULL != m_ADIOS2)
        return;

#if ADIOS2_USE_MPI
    int flag = 0;
    MPI_Initialized(&flag);
    if (!flag)
    {
        RANK_ZERO_MSG("H5VL_ADIOS2 WARNING: MPI is not initialized, will use "
                      "Serial ADIOS\n");
        m_ADIOS2 = adios2_init_serial();
    }
    else
    {
        MPI_Comm comm = MPI_COMM_WORLD;
        if (H5Pget_driver(acc_tpl) == H5FD_MPIO)
        {
            MPI_Info info;
            H5Pget_fapl_mpio(acc_tpl, &comm, &info);
            // int rank;
        }
        MPI_Comm_rank(comm, &m_MPIRank);
        m_ADIOS2 = adios2_init_mpi(comm);
    }
#else
    m_ADIOS2 = adios2_init_serial();
#endif
    REQUIRE_NOT_NULL(m_ADIOS2);
}

void gExitADIOS2()
{
    if (NULL == m_ADIOS2)
        return;
    adios2_remove_all_ios(m_ADIOS2);
    adios2_finalize(m_ADIOS2);
    m_ADIOS2 = NULL;
}

void loadPath(H5VL_ObjDef_t *result, const char *name, H5VL_ObjDef_t *parent)
{
    if (NULL == parent->m_Path)
    {
        result->m_Path = (char *)(SAFE_CALLOC(strlen(name) + 1, sizeof(char)));
        strcpy(result->m_Path, name);
        result->m_Path[strlen(name)] = '\0';
    }
    else if ((strlen(parent->m_Path) == 1) && (parent->m_Path[0] == '/'))
    {
        result->m_Path = (char *)(SAFE_CALLOC(strlen(name) + 2, sizeof(char)));
        sprintf(result->m_Path, "%s%s", parent->m_Path, name);
        result->m_Path[strlen(name) + 1] = '\0';
    }
    else
    {
        size_t size = strlen(name) + strlen(parent->m_Path) + 1;
        if (parent->m_Path[strlen(parent->m_Path) - 1] == '/')
        {
            size = size - 1;
            result->m_Path = (char *)(SAFE_CALLOC(size + 1, sizeof(char)));
            sprintf(result->m_Path, "%s%s", parent->m_Path, name);
        }
        else
        {
            result->m_Path = (char *)(SAFE_CALLOC(size + 1, sizeof(char)));
            sprintf(result->m_Path, "%s/%s", parent->m_Path, name);
        }
        // result->m_Path=
        // (char*)(SAFE_CALLOC(strlen(name)+4+strlen(parent->m_Path),
        // sizeof(char)));
        result->m_Path[size] = '\0';
    }
}

H5VL_ObjDef_t *initVolObj(H5VL_ObjDef_t *parent)
{
    H5VL_ObjDef_t *result = (H5VL_ObjDef_t *)SAFE_CALLOC(1, sizeof(H5VL_ObjDef_t));
    result->m_Parent = parent;

    result->m_ObjPtr = NULL;
    result->m_Path = NULL;
    result->m_SubGroupNames = NULL;

    result->m_NumAttrs = 0;
    result->m_NumVars = 0;
    result->m_NumSubGroups = 0;

    result->m_Vars = NULL;
    result->m_Attrs = NULL;

    if (NULL == parent)
        result->m_FileIO = NULL;
    else
        result->m_FileIO = parent->m_FileIO;

    return result;
}

void gFreeVol(H5VL_ObjDef_t *vol)
{
    if (NULL == vol)
        return;

    if (NULL != vol->m_Vars)
    {
        SAFE_FREE(vol->m_Vars);
    }

    if (NULL != vol->m_Attrs)
    {
        SAFE_FREE(vol->m_Attrs);
    }

    if (NULL != vol->m_SubGroupNames)
    {
        int i;
        for (i = 0; i < vol->m_NumSubGroups; i++)
            SAFE_FREE(vol->m_SubGroupNames[i]);

        SAFE_FREE(vol->m_SubGroupNames);
    }

    SAFE_FREE(vol->m_Path);
    SAFE_FREE(vol);
    vol = NULL;
}

void *gAttrToVolObj(H5VL_AttrDef_t *attr, H5VL_ObjDef_t *parent)
{
    H5VL_ObjDef_t *result = initVolObj(parent);

    result->m_ObjType = ATTR;
    result->m_ObjPtr = attr;

    loadPath(result, attr->m_Name, parent);
    return result;
}

void *gGroupToVolObj(H5VL_GroupDef_t *group, H5VL_ObjDef_t *parent)
{
    H5VL_ObjDef_t *result = initVolObj(parent);
    result->m_ObjType = GROUP;
    result->m_ObjPtr = group;

    loadPath(result, group->m_Name, parent);

    return result;
}

void *gVarToVolObj(H5VL_VarDef_t *var, H5VL_ObjDef_t *parent)
{
    H5VL_ObjDef_t *result = initVolObj(parent);
    result->m_ObjType = VAR;
    result->m_ObjPtr = var;

    loadPath(result, var->m_Name, parent);
    return result;
}

void *gFileToVolObj(H5VL_FileDef_t *f)
{
    H5VL_ObjDef_t *result = initVolObj(NULL);
    result->m_ObjType = ROOT;
    result->m_ObjPtr = f;

    result->m_FileIO = f->m_IO;
    return result;
}

adios2_attribute *gLocateAttrFrom(H5VL_ObjDef_t *owner, const char *attrName)
{
    if (NULL == owner)
        return NULL;

    if (ROOT == owner->m_ObjType)
        return adios2_inquire_attribute(owner->m_FileIO, attrName);

    if ((GROUP == owner->m_ObjType) || (VAR == owner->m_ObjType))
    {
        size_t ss = strlen(owner->m_Path);
        if ('/' == (owner->m_Path)[ss - 1])
        {
            char fullPath[strlen(owner->m_Path) + 4 + strlen(attrName)];
            sprintf(fullPath, "%s%s", owner->m_Path, attrName);
            return adios2_inquire_attribute(owner->m_FileIO, fullPath);
        }
        else
        {
            char fullPath[strlen(owner->m_Path) + 4 + strlen(attrName)];
            sprintf(fullPath, "%s/%s", owner->m_Path, attrName);
            return adios2_inquire_attribute(owner->m_FileIO, fullPath);
        }
    }
    return NULL;
}

//
// returns 0 is obj_name is not under owner
// returns 1 is obj_name is under owner
//
htri_t gExistsUnderGrp(H5VL_ObjDef_t *owner, const char *obj_name)
{
    if (NULL == owner)
        return 0;

    if (ROOT == owner->m_ObjType)
    {
        if (NULL != adios2_inquire_attribute(owner->m_FileIO, obj_name))
            return 1;
        if (NULL != adios2_inquire_variable(owner->m_FileIO, obj_name))
            return 1;
        return 0;
    }

    if (GROUP != owner->m_ObjType)
        return 0;

    char fullPath[strlen(owner->m_Path) + 4 + strlen(obj_name)];
    sprintf(fullPath, "%s/%s", owner->m_Path, obj_name);

    if (NULL != adios2_inquire_attribute(owner->m_FileIO, fullPath))
        return 1;

    if (NULL != adios2_inquire_variable(owner->m_FileIO, fullPath))
        return 1;

    return 0;
}

//
// return:
//    <true>   if obj_name is found and deleted.
//    <falses> otherwise
//
bool gRemoveUnderGrp(H5VL_ObjDef_t *owner, const char *obj_name)
{
    if (NULL == owner)
        return 0;

    adios2_bool result;

    if (ROOT == owner->m_ObjType)
    {
        if (adios2_error_none == adios2_remove_attribute(&result, owner->m_FileIO, obj_name))
            if (adios2_true == result)
                return true;
        if (adios2_error_none == adios2_remove_variable(&result, owner->m_FileIO, obj_name))
            if (adios2_true == result)
                return true;
        return false;
    }

    if (GROUP != owner->m_ObjType)
        return false;

    char fullPath[strlen(owner->m_Path) + 4 + strlen(obj_name)];
    gGenerateFullPath(fullPath, owner->m_Path, obj_name);
    // sprintf(fullPath, "%s/%s", owner->m_Path, obj_name);

    if (adios2_error_none == adios2_remove_attribute(&result, owner->m_FileIO, fullPath))
        if (adios2_true == result)
            return true;
    if (adios2_error_none == adios2_remove_variable(&result, owner->m_FileIO, fullPath))
        if (adios2_true == result)
            return true;

#ifdef NEVER
    return false;
#else
    printf("\n......... NOTE: unable to remove GROUP %s \n\n", fullPath);
    return true;
#endif
}

void gLoadAttrDef(H5VL_AttrDef_t *attrDef)
{
    adios2_attribute *attr = attrDef->m_Attribute;
    if (NULL == attr)
        return;

    adios2_bool boolVal;
    adios2_attribute_is_value(&boolVal, attr);
    if (adios2_true == boolVal)
    {
        attrDef->m_SpaceID = H5Screate(H5S_SCALAR);
        attrDef->m_Size = 1;
    }
    else
    {
        attrDef->m_IsScalar = false;
        adios2_attribute_size(&(attrDef->m_Size), attr);
        attrDef->m_SpaceID = H5Screate(H5S_SIMPLE);
        // NOTE: adios only supports 1 dimension
        hsize_t ss = (hsize_t)(attrDef->m_Size);
        H5Sset_extent_simple(attrDef->m_SpaceID, 1, &ss, NULL);
    }

    adios2_type adiosType;
    adios2_attribute_type(&adiosType, attr);
    attrDef->m_TypeID = gUtilHDF5Type(adiosType);

    /*
    if (adios2_type_string == adiosType) {
      char typeNameStr[30];
      size_t typeNameSize;
      adios2_attribute_type_string (typeNameStr, &typeNameSize, attr);
      // need a way to get actual string value size, not type sie
      //printf(" string type= %s, name=%s, size= %lu\n", attrDef->m_Name,
    typeNameStr, typeNameSize);
    }
    */
}

H5VL_ObjDef_t *gGetVarObjDef(const char *name, H5VL_ObjDef_t *vol)
{
    if (vol->m_ObjType == ROOT)
    { // no vol->m_Path
        if ((strlen(name) == 1) && (name[0] == '/'))
        {
            // root group, not var
            return NULL;
        }
        H5VL_FileDef_t *handle = (H5VL_FileDef_t *)(vol->m_ObjPtr);
        adios2_variable *var = gADIOS2InqVar(vol->m_FileIO, name);
        if (NULL == var)
        {
            if (name[strlen(name) - 1] != '/')
            {
                SHOW_ERROR_MSG("H5VL_ADIOS2: Error: No such variable: %s in file\n ", name);
                return NULL;
            }

            char *n = (char *)(SAFE_CALLOC(strlen(name) + 1, sizeof(char)));
            strcpy(n, name);
            n[strlen(name) - 1] = '\0';
            var = gADIOS2InqVar(vol->m_FileIO, n);
            SAFE_FREE(n);
            if (NULL == var)
            {
                return NULL;
            }
        }
        H5VL_VarDef_t *varDef = gCreateVarDef(name, handle->m_Engine, var, -1, -1);
        return gVarToVolObj(varDef, vol);
    }

    char fullPath[strlen(vol->m_Path) + 4 + strlen(name)];
    gGenerateFullPath(fullPath, vol->m_Path, name);

    if ('/' == name[strlen(name) - 1])
    {
        fullPath[strlen(fullPath) - 1] = '\0';
    }

    adios2_variable *var = gADIOS2InqVar(vol->m_FileIO, fullPath);
    if (NULL == var)
    {
        SHOW_ERROR_MSG("H5VL_ADIOS2: Error: No such variable:: %s in file\n ", fullPath);
        return NULL;
    }

    H5VL_ObjDef_t *curr = vol;
    while (curr)
    {
        if (curr->m_Parent == NULL)
        {
            H5VL_VarDef_t *varDef = gCreateVarDef(
                fullPath, ((H5VL_FileDef_t *)(curr->m_ObjPtr))->m_Engine, var, -1, -1);

            return gVarToVolObj(varDef, vol);
        }
        else
        {
            curr = curr->m_Parent;
        }
    }

    return NULL;
}

H5VL_AttrDef_t *gCreateAttrDef(const char *name, hid_t type_id, hid_t space_id)
{
    H5VL_AttrDef_t *attrDef = (H5VL_AttrDef_t *)SAFE_CALLOC(1, sizeof(H5VL_AttrDef_t));
    attrDef->m_Attribute = NULL;
    attrDef->m_Size = 0;

    attrDef->m_Name = (char *)SAFE_CALLOC(strlen(name) + 1, sizeof(char));
    sprintf(attrDef->m_Name, "%s", name);
    attrDef->m_Name[strlen(name)] = '\0';

    attrDef->m_TypeID = type_id;

    attrDef->m_IsScalar = true;

    if (space_id != -1) // most likely from H5Dcreate
        attrDef->m_SpaceID = H5Scopy(space_id);
    else
        attrDef->m_SpaceID = -1;

    return attrDef;
}

H5VL_VarDef_t *gCreateVarDef(const char *name, adios2_engine *engine, adios2_variable *var,
                             hid_t space_id, hid_t type_id)
{
    if ((-1 == type_id) && (NULL == var))
    {
        printf("UNABLE to create var with unknown var _and_ unknown type");
        return NULL;
    }
    H5VL_VarDef_t *varDef = (H5VL_VarDef_t *)SAFE_CALLOC(1, sizeof(H5VL_VarDef_t));
    varDef->m_Name = (char *)SAFE_CALLOC(strlen(name) + 1, sizeof(char));
    sprintf(varDef->m_Name, "%s", name);

    varDef->m_Engine = engine;
    varDef->m_Variable = var;
    varDef->m_DimCount = (size_t)-1; // default: unknown
    varDef->m_TypeID = -1;           // default: unknown
    varDef->m_Data = NULL;

    if (space_id != -1) // most likely from H5Dcreate
    {
        varDef->m_ShapeID = H5Scopy(space_id);
        varDef->m_DimCount = H5Sget_simple_extent_ndims(space_id);
    }
    else
    { // likely from H5Dopen, so get space info from adios var:
        REQUIRE_NOT_NULL_ERR(var, NULL);
        size_t nDims;
        if (adios2_error_none != adios2_variable_ndims(&nDims, var))
        {
            SAFE_FREE(varDef);
            return NULL;
        }

        varDef->m_DimCount = nDims;

        size_t shape[nDims];
        if (adios2_error_none != adios2_variable_shape(shape, var))
        {
            SAFE_FREE(varDef);
            return NULL;
        }

        hid_t filespace = H5Screate_simple(nDims, (hsize_t *)shape, NULL);
        varDef->m_ShapeID = filespace;
    }

    if (type_id != -1)
        varDef->m_TypeID = type_id;
    else
    {
        adios2_type adiosType;
        adios2_variable_type(&adiosType, var);
        varDef->m_TypeID = gUtilHDF5Type(adiosType);

        /*
        if (adios2_type_string == adiosType) {
          char typeNameStr[30];
          size_t typeNameSize;
          adios2_variable_type_string (typeNameStr, &typeNameSize, var);
          printf(" var type= %s\n", typeNameStr);
        }
        */
    }

    return varDef;
}

H5VL_GroupDef_t *gCreateGroupDef(const char *name)
{
    H5VL_GroupDef_t *grp = (H5VL_GroupDef_t *)SAFE_CALLOC(1, sizeof(H5VL_GroupDef_t));

    grp->m_Name = (char *)SAFE_CALLOC(strlen(name) + 1, sizeof(char));
    sprintf(grp->m_Name, "%s", name);

    return grp;
}

size_t gGetBranchNameLength(H5VL_ObjDef_t *vol, size_t namelen)
{
    if (vol->m_Path != NULL)
        if ('/' == (vol->m_Path)[strlen(vol->m_Path) - 1])
            return namelen - strlen(vol->m_Path); // minus parent path
        else
            return namelen - strlen(vol->m_Path) - 1; // minus parent path & seperator
    else
        return namelen;
}

void gGetBranchName(H5VL_ObjDef_t *vol, const char *fullPath, char *name)
{
    size_t namelen = strlen(fullPath);

    if ('/' == (vol->m_Path)[strlen(vol->m_Path) - 1])
        // minus parent path & seperator
        strncpy(name, fullPath + strlen(vol->m_Path), namelen - strlen(vol->m_Path));
    else
        strncpy(name, fullPath + strlen(vol->m_Path) + 1, namelen - strlen(vol->m_Path) - 1);
}

//
//  returns length of name
//
size_t gGetNameOfNthAttr(H5VL_ObjDef_t *vol, uint32_t idx, char *name)
{
    gLoadContent(vol);
    if (0 == vol->m_NumAttrs)
        return 0;

    if (idx >= vol->m_NumAttrs)
        return 0;

    adios2_attribute *curr = (vol->m_Attrs)[idx];

    size_t namelen;
    adios2_attribute_name(NULL, &namelen, curr);

    if ((NULL != name) && (NULL == vol->m_Path))
        adios2_attribute_name(name, &namelen, curr);
    else if (NULL != name)
    {
        char fullPath[namelen + 1];
        adios2_attribute_name(fullPath, &namelen, curr);
        fullPath[namelen] = '\0';

        gGetBranchName(vol, fullPath, name);
        // printf(".... [%s] vs [%s]\n", fullPath, name);
    }

    return gGetBranchNameLength(vol, namelen);
}

//
//  returns length of name
//  of either dataset or subgroup
//  called from: H5Gget_info, then H5Gget_objname_by_idx
//  (which calls H5Lget_name.. )
//
size_t gGetNameOfNthItem(H5VL_ObjDef_t *vol, uint32_t idx, char *name)
{
    gLoadContent(vol);

    /*
    if (idx < vol->m_NumAttrs)
      return gGetNameOfNthAttr(vol, idx, name);

    uint32_t vIdx = idx - vol->m_NumAttrs;
    if (vIdx < vol->m_NumVars)
    */
    if (idx < vol->m_NumVars)
    {
        adios2_variable *curr = (vol->m_Vars)[idx];

        size_t namelen;
        adios2_variable_name(NULL, &namelen, curr);

        if ((NULL != name) && (NULL == vol->m_Path))
            adios2_variable_name(name, &namelen, curr);
        else if (NULL != name)
        {
            char fullPath[namelen + 1];
            adios2_variable_name(fullPath, &namelen, curr);
            fullPath[namelen] = '\0';

            gGetBranchName(vol, fullPath, name);
            // printf(".... [%s] vs [%s]\n", fullPath, name);
        }
        return gGetBranchNameLength(vol, namelen);
    }

    // now try subgroups
    if (0 == vol->m_NumSubGroups)
        return 0;

    uint32_t vIdx = idx - vol->m_NumVars;
    if (vIdx >= vol->m_NumSubGroups)
        return 0;

    char *n = vol->m_SubGroupNames[vIdx];
    if (NULL != name)
    {
        sprintf(name, "%s", n);
    }
    return strlen(n);
}

//
// load attributes and vars that belongs to this obj
//
void gLoadContent(H5VL_ObjDef_t *obj)
{
    if (0 < (obj->m_NumVars + obj->m_NumAttrs))
        return;

    size_t nvars, nattrs;
    {
        if ((ROOT == obj->m_ObjType) || (GROUP == obj->m_ObjType))
        {
            adios2_variable **adios_vars;
            adios2_inquire_group_variables(&adios_vars, obj->m_Path, &nvars, obj->m_FileIO);

            obj->m_NumVars = nvars;
            if (nvars > 0)
                obj->m_Vars = adios_vars;
        }
    }

    {
        if (ATTR != obj->m_ObjType)
        {
            adios2_attribute **adios_attrs;
            adios2_inquire_group_attributes(&adios_attrs, obj->m_Path, &nattrs, obj->m_FileIO);

            obj->m_NumAttrs = nattrs;
            if (nattrs > 0)
                obj->m_Attrs = adios_attrs;
        }
    }
}

//
//
void gLoadSubGroups(H5VL_ObjDef_t *obj)
{
    // NOTE: subgroups will be supported when ADIOS supports hierarchy properly
    if ((ROOT == obj->m_ObjType) || (GROUP == obj->m_ObjType))
    {
        if (obj->m_NumSubGroups > 0)
            return;

        adios2_inquire_subgroups(&(obj->m_SubGroupNames), obj->m_Path, &(obj->m_NumSubGroups),
                                 obj->m_FileIO);
        return;
    }
}

void gChooseEngine(adios2_io *io)
{
    const char *engineType = getenv("ADIOS2_ENGINE");

    if (engineType != NULL)
    {
        if (0 == m_MPIRank)
            printf("  ADIOS2 will apply engine: %s\n", engineType);
        adios2_set_engine(io, engineType);
    }
    else
    {
        adios2_set_engine(io, "BP4");
    }
}

H5VL_FileDef_t *gADIOS2CreateFile(const char *name)
{
    H5VL_FileDef_t *handle = NULL;
    handle = (H5VL_FileDef_t *)SAFE_CALLOC(1, sizeof(H5VL_FileDef_t));

    // if (0 == m_MPIRank)
    // printf("===========> creating: file:%s\n", name);

    handle->m_IO = adios2_declare_io(m_ADIOS2, name);

    if (NULL == handle->m_IO)
        handle->m_IO = adios2_at_io(m_ADIOS2, name);

    if (NULL == handle->m_IO)
    {
        SAFE_FREE(handle);
        return NULL;
    }

    adios2_set_parameter(handle->m_IO, "Profile", "Off");
    gChooseEngine(handle->m_IO);
    handle->m_Engine = adios2_open(handle->m_IO, name, adios2_mode_write);

    return handle;
}

H5VL_FileDef_t *gADIOS2OpenFile(const char *name)
{
    H5VL_FileDef_t *handle = NULL;
    handle = (H5VL_FileDef_t *)SAFE_CALLOC(1, sizeof(H5VL_FileDef_t));

    // printf("===========> opening: file:%s\n", name);

    handle->m_IO = adios2_declare_io(m_ADIOS2, name);

    gChooseEngine(handle->m_IO);
    handle->m_Engine = adios2_open(handle->m_IO, name, adios2_mode_read);

    char engineType[10];
    size_t engineTypeSize;
    adios2_engine_get_type(engineType, &engineTypeSize, handle->m_Engine);
    printf("==> engine type:%s", engineType);
    if ((engineType[0] == 'B') && (engineType[2] == '5'))
    { // BP5's IO is empty until beginStep is called.
      // since there is no concept for steps in HDF5, sufficient to call
      // begin/end steps here once.
        H5VL_adios2_begin_read_step(name);
        H5VL_adios2_endstep(name);
    }
    return handle;
}

void gADIOS2CloseFile(H5VL_FileDef_t *handle)
{
    if (NULL == handle)
        return;

    if (NULL != handle->m_Engine)
        adios2_close(handle->m_Engine);

    SAFE_FREE(handle);
}

adios2_variable *gADIOS2InqVar(adios2_io *io, const char *name)
{
    return adios2_inquire_variable(io, name);
}

herr_t gADIOS2ReadVar(H5VL_VarDef_t *varDef)
{
    REQUIRE_NOT_NULL_ERR(varDef, -1);
    REQUIRE_NOT_NULL_ERR(varDef->m_Variable, -1);

    int varDim = varDef->m_DimCount;
    if (varDim < 0)
        return -1;

    if (varDim > 0)
    {
        size_t start[varDim], count[varDim];
        if (H5VL_CODE_FAIL == gUtilADIOS2GetBlockInfo(varDef->m_HyperSlabID, start, count, varDim))
            return -1;

        adios2_set_selection(varDef->m_Variable, varDef->m_DimCount, start, count);

        if (varDef->m_MemSpaceID > 0)
        {
            RANK_ZERO_MSG("\n## No memory space is supported for reading in ADIOS...\n");
        }
    }
    adios2_get(varDef->m_Engine, varDef->m_Variable, varDef->m_Data, adios2_mode_sync);

    return 0;
}

adios2_variable *gADIOS2DefineVar(adios2_io *io, H5VL_VarDef_t *varDef)
{
    adios2_variable *variable = adios2_inquire_variable(io, varDef->m_Name);
    if (NULL == variable)
    {
        adios2_type varType = gUtilADIOS2Type(varDef->m_TypeID);

        size_t varDim = 1; // gUtilGetSpaceInfo(;

        varDim = gUtilADIOS2GetDim(varDef->m_ShapeID);

        if (0 == varDim)
        { //  scalar
            variable = adios2_define_variable(io, varDef->m_Name, varType, varDim, NULL, NULL, NULL,
                                              adios2_constant_dims_true);
        }
        else
        {
            varDef->m_DimCount = varDim;

            size_t shape[varDim];
            gUtilADIOS2GetShape(varDef->m_ShapeID, shape, varDim);

            variable = adios2_define_variable(io, varDef->m_Name, varType, varDim, shape, NULL,
                                              NULL, adios2_constant_dims_false);
        }
    }

    return variable;
}

adios2_variable *gADIOS2CreateVar(adios2_io *io, H5VL_VarDef_t *varDef)
{
    REQUIRE_NOT_NULL_ERR(varDef, NULL);

    adios2_variable *variable = adios2_inquire_variable(io, varDef->m_Name);
    size_t varDim = 1;

    if (NULL == variable)
    {
        adios2_type varType = gUtilADIOS2Type(varDef->m_TypeID);
        if (adios2_type_unknown == varType)
        {
            SHOW_ERROR_MSG("... ERROR! Unsupported type. Cannot identify var type.  %s\n",
                           varDef->m_Name);
            return NULL;
        }

        varDim = gUtilADIOS2GetDim(varDef->m_ShapeID);

        if (0 == varDim)
        { //  scalar
            variable = adios2_define_variable(io, varDef->m_Name, varType, varDim, NULL, NULL, NULL,
                                              adios2_constant_dims_true);
        }
        else
        {
            varDef->m_DimCount = varDim;

            size_t shape[varDim];
            gUtilADIOS2GetShape(varDef->m_ShapeID, shape, varDim);

            size_t start[varDim], count[varDim];
            if (H5VL_CODE_FAIL ==
                gUtilADIOS2GetBlockInfo(varDef->m_HyperSlabID, start, count, varDim))
                return NULL;

            variable = adios2_define_variable(io, varDef->m_Name, varType, varDim, shape, start,
                                              count, adios2_constant_dims_true);
        }
    }

    if (NULL != varDef->m_Data)
    {
        varDim = gUtilADIOS2GetDim(varDef->m_ShapeID);
        if (varDim > 0)
        {
            size_t start[varDim], count[varDim];
            if (H5VL_CODE_FAIL ==
                gUtilADIOS2GetBlockInfo(varDef->m_HyperSlabID, start, count, varDim))
                return NULL;
            adios2_set_selection(variable, varDim, start, count);

            if ((varDef->m_MemSpaceID > 0) && (varDef->m_MemSpaceID != varDef->m_ShapeID))
            {
                RANK_ZERO_MSG("\n## No support of memory space for writing in ADIOS.\n");
            }
        }
        adios2_put(varDef->m_Engine, variable, varDef->m_Data, adios2_mode_sync);
    }

    return variable;
}

adios2_attribute *gADIOS2CreateAttr(adios2_io *io, H5VL_AttrDef_t *input, const char *fullPath)
{
    adios2_type attrType = gUtilADIOS2Type(input->m_TypeID);

    if (adios2_type_unknown == attrType)
    {
        SHOW_ERROR_MSG("... ERROR Unsupported type. Cannot create attr %s\n", fullPath);
        return NULL;
    }

    size_t attrDim = 0;

    if (NULL != adios2_inquire_attribute(io, fullPath))
    {
        /*
        printf("... adios2 attribute %s is already created.\n", fullPath);
        return NULL;
        */
        gADIOS2RemoveAttr(io, fullPath);
    }

    if (gUtilADIOS2IsScalar(input->m_SpaceID))
    {
        return adios2_define_attribute(io, fullPath, attrType, input->m_Data);
    }
    else
    {
        attrDim = gUtilADIOS2GetDim(input->m_SpaceID);

        if (1 != attrDim)
        {
            printf("Unable to support 2+D arrays  in ADIOS2 attributes. Use "
                   "Vars instead.");
            return NULL;
        }

        size_t shape[attrDim];
        gUtilADIOS2GetShape(input->m_SpaceID, shape, attrDim);

        if (adios2_type_string == attrType)
        {
            size_t strSize = H5Tget_size(input->m_TypeID);
            htri_t isVariableSize = H5Tis_variable_str(input->m_TypeID);
            /*printf("attr: %s is string array, is variable str? %d, "
                   "H5Tget_size=%lu\n",
                   input->m_Name, isVariableSize, strSize);
            */
            if (isVariableSize)
                return adios2_define_attribute_array(io, fullPath, attrType, (input->m_Data),
                                                     shape[0]);
            else
            {
                int i;

                char *arrayOfStr[shape[0]];
                for (i = 0; i < shape[0]; i++)
                {
                    arrayOfStr[i] = malloc(sizeof(char) * strSize + 1);
                    strncpy(arrayOfStr[i], (char *)(input->m_Data) + strSize * i, strSize);
                    arrayOfStr[i][strSize] = '\0';

                    // printf(".... output attr: %d, [%s]", i, arrayOfStr[i]);
                }
                adios2_attribute *result =
                    adios2_define_attribute_array(io, fullPath, attrType, arrayOfStr, shape[0]);
                for (i = 0; i < shape[0]; i++)
                    free(arrayOfStr[i]);
                return result;
            }
        }
        else
            // return adios2_define_attribute_array(io, fullPath, attrType,
            // &(input->m_Data), shape[0]);
            return adios2_define_attribute_array(io, fullPath, attrType, (input->m_Data), shape[0]);
    }
}

adios2_attribute *gADIOS2InqAttr(adios2_io *io, const char *name)
{
    return adios2_inquire_attribute(io, name);
}

bool gADIOS2RemoveAttr(adios2_io *io, const char *name)
{
    adios2_bool result;
    adios2_remove_attribute(&result, io, name);
    return (adios2_true == result);
}

hid_t H5VL_adios_register(void)
{
    if (H5I_VOL != H5Iget_type(H5VL_ADIOS_g))
    {

        H5VL_ADIOS_g = H5VLregister_connector(&H5VL_adios2_def, H5P_DEFAULT);
        if (H5VL_ADIOS_g <= 0)
        {
            SHOW_ERROR_MSG("  [ECP ADIOS VOL ERROR] with H5VLregister_connector\n");
            return -1;
        }
    }

    return H5VL_ADIOS_g;
}

/*
 */
herr_t H5P_unset_adios() { return H5VLunregister_connector(H5VL_ADIOS_g); }
