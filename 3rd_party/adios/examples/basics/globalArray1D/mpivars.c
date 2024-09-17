/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Created by Dmitry Ganyushin ganyushindi@ornl.gov
 */
#include "mpivars.h"
#include <mpi.h>

int rank, nproc;
int wrank, wnproc;
MPI_Comm app_comm;

void init_mpi(int color, int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    /* World comm spans all applications started with the same mpirun command */
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wnproc);

    /* Have to split and create a 'world' communicator for this app only
     color must be unique for each application*/
    MPI_Comm_split(MPI_COMM_WORLD, color, wrank, &app_comm);
    MPI_Comm_rank(app_comm, &rank);
    MPI_Comm_size(app_comm, &nproc);
    return;
}

void finalize_mpi() { MPI_Finalize(); }
