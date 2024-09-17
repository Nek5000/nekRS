#ifndef __H5VOL_DEFS
#define __H5VOL_DEFS

#include <hdf5.h>

#include "H5VolError.h"
#include "H5VolUtil.h"
#include <adios2_c.h>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define H5_ATTR_UNUSED /**/

// Define the values for the first entries of
// VOL struct H5VL_class_t
#define H5VL_ADIOS2_NAME "ADIOS2_VOL"
#define H5VL_ADIOS2_VALUE 511 /* VOL connector ID */
#define H5VL_ADIOS2_VERSION 0
#define H5VL_ADIOS2_CAP_FLAGS                                                                      \
    (H5VL_CAP_FLAG_FILE_BASIC | H5VL_CAP_FLAG_GROUP_BASIC | H5VL_CAP_FLAG_DATASET_BASIC |          \
     H5VL_CAP_FLAG_ATTR_BASIC | H5VL_CAP_FLAG_OBJECT_BASIC)

typedef struct H5VL_ADIOS2_t
{
    adios2_engine *m_Engine;
    adios2_io *m_IO;
} H5VL_FileDef_t;

typedef struct H5VL_VARDef_t
{
    char *m_Name;
    hid_t m_ShapeID;
    hid_t m_TypeID;
    hid_t m_MemSpaceID;
    hid_t m_HyperSlabID;
    hid_t m_PropertyID;
    void *m_Data; // used for both write & read so not going to be const
    adios2_engine *m_Engine;
    adios2_variable *m_Variable; // mainly for read reuse
    size_t m_DimCount;
} H5VL_VarDef_t;

typedef struct H5VL_ATTRDef_t
{
    char *m_Name;
    hid_t m_SpaceID;
    hid_t m_TypeID;
    const void *m_Data;
    bool m_IsScalar;
    size_t m_Size;                 // num of elements
    adios2_attribute *m_Attribute; // for read reuse
} H5VL_AttrDef_t;

typedef struct H5VL_GROUPDef_t
{
    char *m_Name;
} H5VL_GroupDef_t;

//
// if (m_Parent == NULL)  => FILE
// otherwise type is the non null pointer
//
// m_Parent is of type < H5VL_ObjDef_t >
//
enum H5VL_PTR_TYPES
{
    ROOT = 3,
    GROUP = 2,
    VAR = 1,
    ATTR = 0
};

typedef struct H5VL_ObjDef
{
    void *m_ObjPtr;
    void *m_Parent;
    char *m_Path; // NULL if is ROOT.
    enum H5VL_PTR_TYPES m_ObjType;

    size_t m_NumVars;
    adios2_variable **m_Vars;

    size_t m_NumAttrs;
    adios2_attribute **m_Attrs;

    size_t m_NumSubGroups;
    char **m_SubGroupNames;

    adios2_io *m_FileIO;
} H5VL_ObjDef_t;

#ifdef __cplusplus
extern "C" {
#endif
//
// file functions
//
extern void *H5VL_adios2_file_create(const char *name, unsigned flags, hid_t fcpl_id, hid_t fapl_id,
                                     hid_t dxpl_id, void **req);

extern void *H5VL_adios2_file_open(const char *name, unsigned flags, hid_t fapl_id, hid_t dxpl_id,
                                   void **req);

extern herr_t H5VL_adios2_file_specific(void *file, H5VL_file_specific_args_t *args, hid_t dxpl_id,
                                        void **req);

extern herr_t H5VL_adios2_file_close(void *file, hid_t dxpl_id, void **req);

//
// attr functions
//
extern void *H5VL_adios2_attr_create(void *obj, const H5VL_loc_params_t *loc_params,
                                     const char *name, hid_t type_id, hid_t space_id, hid_t acpl_id,
                                     hid_t aapl_id, hid_t dxpl_id, void **req);

extern void *H5VL_adios2_attr_open(void *obj, const H5VL_loc_params_t *loc_params, const char *name,
                                   hid_t aapl_id, hid_t dxpl_id, void **req);

extern herr_t H5VL_adios2_attr_read(void *attrObj, hid_t mem_type_id, void *buf, hid_t dxpl_id,
                                    void **req);

extern herr_t H5VL_adios2_attr_write(void *attr, hid_t mem_type_id, const void *buf, hid_t dxpl_id,
                                     void **req);

extern herr_t H5VL_adios2_attr_get(void *obj, H5VL_attr_get_args_t *args, hid_t dxpl_id,
                                   void **req);

extern herr_t H5VL_adios2_attr_close(void *attr, hid_t dxpl_id, void **req);

extern herr_t H5VL_adios2_attr_specific(void *obj, const H5VL_loc_params_t *loc_params,
                                        H5VL_attr_specific_args_t *args, hid_t dxpl_id, void **req);

//
// object functions:
//
extern void *H5VL_adios2_object_open(void *obj, const H5VL_loc_params_t *loc_params,
                                     H5I_type_t *opened_type, hid_t H5_ATTR_UNUSED dxpl_id,
                                     void H5_ATTR_UNUSED **req);

extern herr_t H5VL_adios2_object_get(void *obj, const H5VL_loc_params_t *loc_params,
                                     H5VL_object_get_args_t *args, hid_t H5_ATTR_UNUSED dxpl_id,
                                     void H5_ATTR_UNUSED **req);

// dataset functions:
extern void *H5VL_adios2_dataset_create(void *obj, const H5VL_loc_params_t *loc_params,
                                        const char *name, hid_t lcpl_id, hid_t type_id,
                                        hid_t space_id, hid_t dcpl_id, hid_t dapl_id, hid_t dxpl_id,
                                        void **req);

extern void *H5VL_adios2_dataset_open(void *obj, const H5VL_loc_params_t *loc_params,
                                      const char *name, hid_t dapl_id, hid_t dxpl_id, void **req);

extern herr_t H5VL_adios2_dataset_read(size_t count, void *dset[], hid_t mem_type_id[],
                                       hid_t mem_space_id[], hid_t file_space_id[], hid_t dxpl_id,
                                       void *buf[], void **req);

extern herr_t H5VL_adios2_dataset_get(void *dset, H5VL_dataset_get_args_t *args, hid_t dxpl_id,
                                      void **req);

extern herr_t H5VL_adios2_dataset_write(size_t count, void *dset[], hid_t mem_type_id[],
                                        hid_t mem_space_id[], hid_t file_space_id[], hid_t dxpl_id,
                                        const void *buf[], void **req);

extern herr_t H5VL_adios2_dataset_close(void *dset, hid_t dxpl_id, void **req);

//
// link functions:
//
extern herr_t H5VL_adios2_link_specific(void *obj, const H5VL_loc_params_t *loc_params,
                                        H5VL_link_specific_args_t *args,
                                        hid_t H5_ATTR_UNUSED dxpl_id, void H5_ATTR_UNUSED **req);

extern herr_t H5VL_adios2_link_get(void *obj, const H5VL_loc_params_t *loc_params,
                                   H5VL_link_get_args_t *args, hid_t H5_ATTR_UNUSED dxpl_id,
                                   void H5_ATTR_UNUSED **req);

//
// group functions:
//
extern void *H5VL_adios2_group_create(void *obj, const H5VL_loc_params_t *loc_params,
                                      const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id,
                                      hid_t dxpl_id, void **req);

extern herr_t H5VL_adios2_group_close(void *obj, hid_t dxpl_id, void **req);

extern void *H5VL_adios2_group_open(void *obj, const H5VL_loc_params_t *loc_params,
                                    const char *name, hid_t gapl_id, hid_t dxpl_id, void **req);

extern herr_t H5VL_adios2_group_get(void *obj, H5VL_group_get_args_t *args,
                                    hid_t H5_ATTR_UNUSED dxpl_id, void H5_ATTR_UNUSED **req);

//
// general definitions:
//
extern void *gAttrToVolObj(H5VL_AttrDef_t *attr, H5VL_ObjDef_t *parent);

extern void *gFileToVolObj(H5VL_FileDef_t *f);

extern void *gGroupToVolObj(H5VL_GroupDef_t *group, H5VL_ObjDef_t *parent);

extern void *gVarToVolObj(H5VL_VarDef_t *var, H5VL_ObjDef_t *parent);

extern void gFreeVol(H5VL_ObjDef_t *vol);

extern adios2_attribute *gLocateAttrFrom(H5VL_ObjDef_t *owner, const char *name);
extern htri_t gExistsUnderGrp(H5VL_ObjDef_t *owner, const char *obj_name);
extern bool gRemoveUnderGrp(H5VL_ObjDef_t *owner, const char *obj_name);

extern void gLoadAttrDef(H5VL_AttrDef_t *attrDef);
extern void gLoadContent(H5VL_ObjDef_t *obj);
extern void gLoadSubGroups(H5VL_ObjDef_t *obj);

extern size_t gGetNameOfNthAttr(H5VL_ObjDef_t *obj, uint32_t idx, char *name);
extern size_t gGetNameOfNthItem(H5VL_ObjDef_t *obj, uint32_t idx, char *name);

extern H5VL_ObjDef_t *gGetVarObjDef(const char *fullPath, H5VL_ObjDef_t *vol);
extern H5VL_VarDef_t *gCreateVarDef(const char *name, adios2_engine *engine, adios2_variable *var,
                                    hid_t space_id, hid_t type_id);

extern H5VL_AttrDef_t *gCreateAttrDef(const char *name, hid_t type_id, hid_t space_id);
extern H5VL_GroupDef_t *gCreateGroupDef(const char *name);

extern void gGenerateFullPath(char *fullPath, const char *parentPath, const char *name);

/*
 */

extern void gInitADIOS2(hid_t acc_tpl);
extern void gExitADIOS2();

// extern adios2_attribute *gADIOS2CreateAttr(H5VL_AttrDef_t *input);
extern adios2_attribute *gADIOS2CreateAttr(adios2_io *io, H5VL_AttrDef_t *input,
                                           const char *fullPath);
extern adios2_attribute *gADIOS2InqAttr(adios2_io *io, const char *name);
extern bool gADIOS2RemoveAttr(adios2_io *io, const char *name);

extern H5VL_FileDef_t *gADIOS2CreateFile(const char *name);
extern H5VL_FileDef_t *gADIOS2OpenFile(const char *name);
extern void gADIOS2CloseFile(H5VL_FileDef_t *handle);

extern adios2_variable *gADIOS2CreateVar(adios2_io *io, H5VL_VarDef_t *var);
extern adios2_variable *gADIOS2DefineVar(adios2_io *io, H5VL_VarDef_t *varDef);
extern adios2_variable *gADIOS2InqVar(adios2_io *io, const char *name);
extern herr_t gADIOS2ReadVar(H5VL_VarDef_t *var);

#ifdef __cplusplus
}
#endif

#endif
