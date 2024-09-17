/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Info.h
 *
 *  Created on: June 22, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "Info.h"
#include "adios2/common/ADIOSConfig.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

const int adios2_version_major = ADIOS2_VERSION_MAJOR;
const int adios2_version_minor = ADIOS2_VERSION_MINOR;
const int adios2_version_patch = ADIOS2_VERSION_PATCH;
const char adios2_version_str[] = ADIOS2_VERSION_STR;

static const char *const engines[] = {"BP3",
                                      "BP4",
                                      "BP5",
#ifdef ADIOS2_HAVE_HDF5
                                      "HDF5",
#endif
#ifdef ADIOS2_HAVE_SST
                                      "SST",
#endif
#ifdef ADIOS2_HAVE_MPI
                                      "SSC",
#endif
#ifdef ADIOS2_HAVE_DataMan
                                      "DataMan",
#endif
#ifdef ADIOS2_HAVE_DataSpaces
                                      "DataSpaces",
#endif
                                      "Inline",
#ifdef ADIOS2_HAVE_DAOS
                                      "DAOS",
#endif
#ifdef ADIOS2_HAVE_MHS
                                      "MHS",
#endif
#ifdef ADIOS2_HAVE_CATALYST
                                      "ParaViewADIOSInSituEngine",
#endif
                                      "Null",
                                      "Skeleton",
                                      nullptr};

void adios2_available_engines(size_t *nentries, const char *const **list)
{
    *nentries = (sizeof(engines) / sizeof(const char *)) - 1;
    *list = engines;
}

static const char *const operators[] = {
#ifdef ADIOS2_HAVE_BZIP2
    "BZip2",
#endif
#ifdef ADIOS2_HAVE_BLOSC2
    "Blosc",
#endif
#ifdef ADIOS2_HAVE_MGARD
    "MGARD",
    "MGARDPlus",
#endif
#ifdef ADIOS2_HAVE_MGARD_MDR
    "MDR",
#endif
#ifdef ADIOS2_HAVE_SZ
    "SZ",
#endif
#ifdef ADIOS2_HAVE_ZFP
    "ZFP",
#endif
#ifdef ADIOS2_HAVE_PNG
    "PNG",
#endif
#ifdef ADIOS2_HAVE_SIRIUS
    "Sirius",
#endif
#ifdef ADIOS2_HAVE_LIBPRESSIO
    "libpressio",
#ifdef ADIOS2_HAVE_SODIUM
    "Sodium plugin",
#endif
#endif
    nullptr};

void adios2_available_operators(size_t *nentries, const char *const **list)
{
    *nentries = (sizeof(operators) / sizeof(const char *)) - 1;
    *list = operators;
}

const char *adios2_feature_list[] = {ADIOS2_FEATURE_LIST};

void adios2_available_features(size_t *nentries, const char *const **list)
{
    *nentries = (sizeof(adios2_feature_list) / sizeof(const char *)) - 1;
    *list = adios2_feature_list;
}

#ifdef __cplusplus
} // end extern C
#endif
