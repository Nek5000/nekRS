#ifndef ADIOS_VOL_H
#define ADIOS_VOL_H

#include "H5VolReadWrite.h"

#ifdef __cplusplus
extern "C" {
#endif

// extern void H5VL_ADIOS2_set(hid_t fapl);
void H5VL_ADIOS2_set(hid_t fapl);
void H5VL_ADIOS2_unset();

#ifdef __cplusplus
}
#endif
#endif
