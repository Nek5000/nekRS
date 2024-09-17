/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ds_writer.c
 *
 *  Created on: Jan 4, 2019
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */
#include <assert.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/toolkit/dataspaces/DSpacesConfig.h"
#ifdef HAVE_DSPACES2
#include "dspaces.h"
#else
#include "dataspaces.h"
#endif /* HAVE_DSPACES2 */
#include "ds.h"
#include <mpi.h>

#define MAX_DS_NAMELEN 128

#ifdef HAVE_DSPACES2
static dspaces_client_t client = NULL;

dspaces_client_t *get_client_handle() { return (&client); }
#endif /* HAVE_DSPACES2 */

static int globals_rank;
static int globals_adios_appid = -1;
static int globals_adios_was_set = 0;
void globals_adios_set_application_id(int id)
{
    globals_adios_appid = id;
    globals_adios_was_set = 1;
}

int globals_adios_get_application_id(int *was_set)
{
    *was_set = globals_adios_was_set;
    return globals_adios_appid;
}

enum DATASPACES_CONNECTION
{
    dataspaces_disconnected = 0,
    dataspaces_connected_from_reader = 1,
    dataspaces_connected_from_writer = 2,
    dataspaces_connected_from_both = 3
};
static enum DATASPACES_CONNECTION globals_adios_connected_to_dataspaces = dataspaces_disconnected;

void globals_adios_set_dataspaces_connected_from_reader()
{
    if (globals_adios_connected_to_dataspaces == dataspaces_disconnected)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_reader;
    else if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_writer)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_both;
}
void globals_adios_set_dataspaces_disconnected_from_reader()
{
    if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_reader)
    {
        globals_adios_connected_to_dataspaces = dataspaces_disconnected;
#ifdef HAVE_DSPACES2
        if (globals_rank == 0)
        {
            dspaces_kill(client);
        }
#endif
    }
    else if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_both)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_writer;
}
void globals_adios_set_dataspaces_connected_from_writer()
{
    if (globals_adios_connected_to_dataspaces == dataspaces_disconnected)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_writer;
    else if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_reader)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_both;
}
void globals_adios_set_dataspaces_disconnected_from_writer()
{
    if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_writer)
    {
        globals_adios_connected_to_dataspaces = dataspaces_disconnected;
#ifdef HAVE_DSPACES2
        if (globals_rank == 0)
        {
            dspaces_kill(client);
        }
#endif
    }
    else if (globals_adios_connected_to_dataspaces == dataspaces_connected_from_both)
        globals_adios_connected_to_dataspaces = dataspaces_connected_from_reader;
}
int globals_adios_is_dataspaces_connected()
{
    return (globals_adios_connected_to_dataspaces != dataspaces_disconnected);
}
int globals_adios_is_dataspaces_connected_from_reader()
{
    return (globals_adios_connected_to_dataspaces == dataspaces_connected_from_reader ||
            globals_adios_connected_to_dataspaces == dataspaces_connected_from_both);
}
int globals_adios_is_dataspaces_connected_from_writer()
{
    return (globals_adios_connected_to_dataspaces == dataspaces_connected_from_writer ||
            globals_adios_connected_to_dataspaces == dataspaces_connected_from_both);
}
int globals_adios_is_dataspaces_connected_from_both()
{
    return (globals_adios_connected_to_dataspaces == dataspaces_connected_from_both);
}

int adios_dataspaces_init(void *comm, DsData *md)
{
    int nproc;
    int rank, err;
    if (!globals_adios_is_dataspaces_connected())
    {

        MPI_Comm_rank(*(MPI_Comm *)comm, &rank);
        MPI_Comm_size(*(MPI_Comm *)comm, &nproc);
        globals_rank = rank;

        if (md->appid == 0)
        {
            md->appid = 1;
            fprintf(stderr, "AppID not found in xml file. Setting it as 1 for "
                            "the Writer in DATASPACES\n");
        }

#ifdef HAVE_DSPACES2
        err = dspaces_init(rank, &client);
        if (err != dspaces_SUCCESS)
        {
            fprintf(stderr, "Failed to connect to the DataSpaces server.\n");
            return (err);
        }
#else
        err = dspaces_init(nproc, md->appid, (MPI_Comm *)comm, NULL);
        if (err < 0)
        {
            fprintf(stderr, "Failed to connect with DATASPACES\n");
            return err;
        }
#endif /* HAVE_DSPACES2 */
        md->rank = rank;
        md->mpi_comm = *(MPI_Comm *)comm;
    }
    globals_adios_set_dataspaces_connected_from_writer();
    return 0;
}

void adios_dataspaces_open(char *fname, DsData *md)
{
    // if we have MPI and a communicator, we can get the exact size of this
    // application that we need to tell DATASPACES
    MPI_Comm_rank(md->mpi_comm, &(md->rank));
    MPI_Comm_size(md->mpi_comm, &(md->peers));
#ifndef HAVE_DSPACES2
    dspaces_lock_on_write(fname, &md->mpi_comm);
#endif /* HAVE_DSPACES2 */
}

int adios_read_dataspaces_init(void *comm, DsData *md)
{
    int nproc;
    int rank, err;
    if (!globals_adios_is_dataspaces_connected())
    {

        MPI_Comm_rank(*(MPI_Comm *)comm, &rank);
        MPI_Comm_size(*(MPI_Comm *)comm, &nproc);
        globals_rank = rank;

        if (md->appid == 0)
        {
            md->appid = 2;
            fprintf(stderr, "AppID not found in xml file. Setting it as 2 for "
                            "the Reader in DATASPACES\n");
        }
#ifdef HAVE_DSPACES2
        err = dspaces_init(rank, &client);
        if (err != dspaces_SUCCESS)
        {
            fprintf(stderr, "Failed to connect to the DataSpaces server.\n");
            return (err);
        }
#else
        err = dspaces_init(nproc, md->appid, (MPI_Comm *)comm, NULL);
        if (err < 0)
        {
            fprintf(stderr, "Failed to connect with DATASPACES\n");
            return err;
        }
#endif /* HAVE_DSPACES2 */
        md->rank = rank;
        md->mpi_comm = *(MPI_Comm *)comm;
    }
    globals_adios_set_dataspaces_connected_from_reader();
    return 0;
}
