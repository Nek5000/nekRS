/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ds_data.h
 *
 *  Created on: Dec 12, 2018
 *      Author: Pradeep Subedi
 *      		pradeep.subedi@rutgers.edu
 */

#ifndef _DS_DATA_H_
#define _DS_DATA_H_

#ifndef _SYS_TYPES_H_
#include <sys/types.h>
#endif
#include "mpi.h"

#define MAX_DS_NDIM 10

struct adios_ds_data_struct
{
    int rank;          // dataspaces rank or MPI rank if MPI is available
    int peers;         // from xml parameter or group communicator
    int appid;         // from xml parameter or 1
    int n_writes;      // how many times adios_write has been called
    MPI_Comm mpi_comm; // for use in open..close
};

#endif /* _DS_DATA_H_ */
