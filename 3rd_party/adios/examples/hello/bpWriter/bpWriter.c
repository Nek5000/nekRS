/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * bpWriter.c : C bindings version of bpWriter.cpp
 *
 *  Created on: Aug 8, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include <stdio.h>  // printf
#include <stdlib.h> // malloc, free, exit

#include <adios2_c.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

void check_error(const int error)
{
    if (error)
    {
        printf("ERROR: %d\n", error);
        exit(error);
    }
}

void check_handler(const void *handler, const char *message)
{
    if (handler == NULL)
    {
        printf("ERROR: invalid %s handler \n", message);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
    int rank, size;

#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif

    adios2_error errio;
    // application input, data in heap
    const size_t Nx = 10;
    float *myFloats;
    myFloats = malloc(sizeof(float) * Nx);

    unsigned int i;
    for (i = 0; i < Nx; ++i)
    {
        myFloats[i] = (float)i;
    }

#if ADIOS2_USE_MPI
    adios2_adios *adios = adios2_init_mpi(MPI_COMM_WORLD);
#else
    adios2_adios *adios = adios2_init_serial();
#endif

    adios2_step_status err;
    check_handler(adios, "adios");

    adios2_io *io = adios2_declare_io(adios, "BPFile_Write");
    check_handler(io, "io");

    // dims are allocated in the stack
    size_t shape[1];
    shape[0] = (size_t)size * Nx;

    size_t start[1];
    start[0] = (size_t)rank * Nx;

    size_t count[1];
    count[0] = Nx;

    adios2_variable *variable = adios2_define_variable(io, "bpFloats", adios2_type_float, 1, shape,
                                                       start, count, adios2_constant_dims_true);
    check_handler(variable, "variable");

    adios2_engine *engine = adios2_open(io, "myVector_c.bp", adios2_mode_write);
    check_handler(engine, "engine");

    adios2_begin_step(engine, adios2_step_mode_append, 0.0f, &err);

    errio = adios2_put(engine, variable, myFloats, adios2_mode_deferred);
    check_error(errio);

    adios2_end_step(engine);

    errio = adios2_close(engine);
    check_error(errio);

    // deallocate adios
    errio = adios2_finalize(adios);
    check_error(errio);

    free(myFloats);

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
