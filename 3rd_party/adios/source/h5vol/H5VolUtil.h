#ifndef H5VL_UTIL
#define H5VL_UTIL

#include <adios2_c.h>
#include <hdf5.h>

#define H5VL_CODE_SUCC 1
#define H5VL_CODE_FAIL 0

hid_t gUtilHDF5Type(adios2_type adios2Type);

adios2_type gUtilADIOS2Type(hid_t h5Type);

/*
int gUtilADIOS2IsNone(hid_t space_id)
{
// note: using npts is better than type
// for example, the following line will get type SIMPLE instead of NON
// memspace = H5Screate_simple(DIM, count, NULL); H5Sselect_none(memspace);

  hsize_t npts = H5Sget_select_npoints(space_id);
  H5S_class_t stype = H5Sget_simple_extent_type(space_id);
  printf(" type = %d npts=%ld\n", stype, npts);
  if (0 == npts)
    return H5VL_CODE_SUCC;
  if ((H5S_NULL == stype) || (H5S_NO_CLASS == stype))
    return H5VL_CODE_SUCC;

  return H5VL_CODE_FAIL;
}
*/

int gUtilADIOS2IsScalar(hid_t space_id);
int gUtilADIOS2GetDim(hid_t space_id);

//
// h5 uses hsize_t for dimensions (unsigned long long)
// adios uses size_t
//
void gUtilConvert(hsize_t *fromH5, size_t *to, size_t ndims);

int gUtilADIOS2GetShape(hid_t space_id, size_t *shape, size_t ndims);

int gUtilADIOS2GetBlockInfo(hid_t hyperSlab_id, size_t *start, size_t *count, hsize_t ndims);

#endif
