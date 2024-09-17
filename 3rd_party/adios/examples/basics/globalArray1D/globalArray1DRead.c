/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Created by Dmitry Ganyushin ganyushindi@ornl.gov
 */
#include "decomp.h"
#include "mpivars.h"
#include <adios2_c.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void reader(adios2_adios *adios)
{
    int step;
    float *g = NULL;
    const char *streamname = "adios2-global-array-1d-c.bp";
    adios2_step_status err;

    size_t shape[1], start[1], count[1];

    adios2_io *io = adios2_declare_io(adios, "input");

    adios2_engine *engine = adios2_open(io, streamname, adios2_mode_read);
    step = 0;
    while (1)
    {
        adios2_begin_step(engine, adios2_step_mode_read, 10.0, &err);
        if (err == adios2_step_status_end_of_stream)
        {
            break;
        }
        if (err != adios2_step_status_ok)
        {
            printf("Unexpected status when calling adios2_begin_step() for step %d", step);
            break;
        }
        adios2_variable *var_g = adios2_inquire_variable(io, "GlobalArray");
        if (step == 0)
        {
            /* fixed_shape is allocated in the next call*/
            adios2_variable_shape(shape, var_g);
            decomp_1d(shape, start, count);
            g = malloc(count[0] * sizeof(float));

            printf("Read plan rank = %d global shape = %zu local count = %zu offset = %zu\n", rank,
                   shape[0], count[0], start[0]);
        }

        adios2_variable *var = adios2_inquire_variable(io, "GlobalArray");
        if (!var)
        {
            printf("ERROR: Variable 'GlobalArray' was not found in step %d", step);
            break;
        }

        adios2_set_selection(var, 1, start, count);
        // Initiate reading data in default/deferred mode: data is available after end_step
        adios2_get(engine, var, g, adios2_mode_deferred);

        adios2_end_step(engine);
        step++;
    }
    // Close the output
    adios2_close(engine);

    if (g != NULL)
    {
        free(g);
    }
}

int main(int argc, char *argv[])
{
#if ADIOS2_USE_MPI
    init_mpi(123, argc, argv);
#endif

    {
#if ADIOS2_USE_MPI

        adios2_adios *adios = adios2_init_mpi(MPI_COMM_WORLD);
#else
        adios2_adios *adios = adios2_init();
#endif
        reader(adios);
        adios2_finalize(adios);
    }

#if ADIOS2_USE_MPI
    finalize_mpi();
#endif

    return 0;
}
