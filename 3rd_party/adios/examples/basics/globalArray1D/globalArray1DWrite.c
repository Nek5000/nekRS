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

void writer(adios2_adios *adios)
{
    int step, i;
    float *g;
    const int mincount = 2;
    const int maxcount = 5;
    const int numsteps = 5;
    adios2_step_status err;

    size_t shape[1], start[1], count[1];

    /* Application variables
     g = 1D distributed array,
     global shape and per-process size is fixed */
    count[0] = get_random(mincount, maxcount, rank);
    g = malloc((size_t)count[0] * sizeof(float));
    gather_decomp_1d(count, shape, start);

    adios2_io *io = adios2_declare_io(adios, "output");

    adios2_variable *var_g = adios2_define_variable(io, "GlobalArray", adios2_type_float, 1, shape,
                                                    start, count, adios2_constant_dims_true);

    adios2_engine *engine = adios2_open(io, "adios2-global-array-1d-c.bp", adios2_mode_write);
    printf("Decomp rank = %d global shape = %zu local count = %zu  offset = %zu\n", rank, shape[0],
           count[0], start[0]);
    for (step = 0; step < numsteps; step++)
    {
        for (i = 0; i < count[0]; i++)
        {
            g[i] = (float)(rank + (step + 1) / 100.0);
        }

        adios2_begin_step(engine, adios2_step_mode_append, 10.0f, &err);
        adios2_put(engine, var_g, g, adios2_mode_deferred);
        adios2_end_step(engine);
    }
    // Close the output
    adios2_close(engine);
    free(g);

    if (rank == 0)
    {
        printf("Try the following: \n");
        printf("  bpls -la adios2-global-array-1d-c.bp GlobalArray -d -n %zu \n", shape[0]);
        printf("  bpls -la adios2-global-array-1d-c.bp GlobalArray -d -t -n %zu \n ", shape[0]);
        printf("  mpirun -n 2 ./adios2_basics_globalArray1DRead_c \n");
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

        writer(adios);
        adios2_finalize(adios);
    }

#if ADIOS2_USE_MPI
    finalize_mpi();
#endif

    return 0;
}
