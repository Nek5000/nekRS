/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Info.h
 *
 *  Created on: June 22, 2023
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#ifndef ADIOS2_INFO_H_
#define ADIOS2_INFO_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

extern const int adios2_version_major;
extern const int adios2_version_minor;
extern const int adios2_version_patch;
extern const char adios2_version_str[];

/** Return the list of available Engines in the installed adios2 library */
void adios2_available_engines(size_t *nentries, const char *const **list);

/** Return the list of available Engines in the installed adios2 library */
void adios2_available_operators(size_t *nentries, const char *const **list);

/** Return the list of available features in the installed adios2 library */
void adios2_available_features(size_t *nentries, const char *const **list);

#ifdef __cplusplus
} // end extern C
#endif

#endif /* ADIOS2_INFO_H_ */
