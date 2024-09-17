/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ds.h
 *
 *  Created on: Jan 4, 2019
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */

#ifndef SOURCE_ADIOS2_TOOLKIT_DATASPACES_DS_H_
#define SOURCE_ADIOS2_TOOLKIT_DATASPACES_DS_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "adios2/toolkit/dataspaces/DSpacesConfig.h"
#include "ds_data.h"
#include <stddef.h>
#ifdef HAVE_DSPACES2
#include <dspaces.h>
#endif /* HAVE_DSPACES2 */

typedef struct adios_ds_data_struct DsData;

extern int adios_dataspaces_init(void *comm, DsData *md);
extern int adios_read_dataspaces_init(void *comm, DsData *data);

#ifdef HAVE_DSPACES2
extern dspaces_client_t *get_client_handle();
#endif /* HAVE_DSPACES2 */
extern void adios_dataspaces_open(char *fname, DsData *md);
extern void globals_adios_set_dataspaces_connected_from_reader();
extern void globals_adios_set_dataspaces_disconnected_from_reader();
extern void globals_adios_set_dataspaces_connected_from_writer();
extern void globals_adios_set_dataspaces_disconnected_from_writer();
extern int globals_adios_is_dataspaces_connected(); // from any
extern int globals_adios_is_dataspaces_connected_from_reader();
extern int globals_adios_is_dataspaces_connected_from_writer();
extern int globals_adios_is_dataspaces_connected_from_both();
extern void globals_adios_set_application_id(int id);
extern int globals_adios_get_application_id(int *was_set);

#ifdef __cplusplus
}
#endif

#endif /* SOURCE_ADIOS2_TOOLKIT_DATASPACES_DS_H_ */
