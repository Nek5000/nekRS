/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestUtilsCWriter
 *
 *  Created on: Sept 19, 2018
 *      Author: Norbert Podhorszki
 */

#include <stdio.h>
#include <string.h>

#include <adios2_c.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include "../SmallTestData_c.h"

int main(int argc, char *argv[])
{
    int rank = 0;
    int nproc = 1;
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    adios2_adios *adiosH = adios2_init_mpi(MPI_COMM_WORLD);
#else
    adios2_adios *adiosH = adios2_init_serial();
#endif

    char engineName[32] = "BPFile";
    if (argc > 1)
    {
        strncpy(engineName, argv[1], sizeof(engineName) - 1);
    }

    // IO
    adios2_io *ioH = adios2_declare_io(adiosH, "CArrayTypes");
    // Set engine parameters
    adios2_set_engine(ioH, engineName);
    adios2_set_parameter(ioH, "ProfileUnits", "Microseconds");
    adios2_set_parameter(ioH, "Threads", "1");

    // 1D shape
    size_t shape[1], start[1], count[1];
    shape[0] = nproc * data_Nx;
    start[0] = rank * data_Nx;
    count[0] = data_Nx;

    // 2D shape
    size_t shape2[2], start2[2], count2[2];
    shape2[0] = nproc * d2_Nx;
    shape2[1] = d2_Ny;
    start2[0] = rank * d2_Nx;
    start2[1] = 0;
    count2[0] = d2_Nx;
    count2[1] = d2_Ny;

    // Define variables in ioH
    adios2_define_variable(ioH, "nproc", adios2_type_int32_t, 0, NULL, NULL, NULL,
                           adios2_constant_dims_true);

    adios2_define_variable(ioH, "varI8", adios2_type_int8_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varI16", adios2_type_int16_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varI32", adios2_type_int32_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varI64", adios2_type_int64_t, 1, shape, start, count,
                           adios2_constant_dims_true);

    adios2_define_variable(ioH, "varU8", adios2_type_uint8_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varU16", adios2_type_uint16_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varU32", adios2_type_uint32_t, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varU64", adios2_type_uint64_t, 1, shape, start, count,
                           adios2_constant_dims_true);

    adios2_define_variable(ioH, "varR32", adios2_type_float, 1, shape, start, count,
                           adios2_constant_dims_true);
    adios2_define_variable(ioH, "varR64", adios2_type_double, 1, shape, start, count,
                           adios2_constant_dims_true);

    adios2_define_variable(ioH, "R64_2d", adios2_type_double, 2, shape2, start2, count2,
                           adios2_constant_dims_true);

    // Define attributes in ioH
    adios2_define_attribute(ioH, "name", adios2_type_string, "TestUtilsCWrite");
    adios2_define_attribute_array(ioH, "strarray", adios2_type_string, strarray,
                                  sizeof(strarray) / sizeof(char *));
    adios2_define_attribute(ioH, "nwriters", adios2_type_int32_t, &nproc);
    unsigned short shape2D[2] = {(unsigned short)d2_Nx, (unsigned short)d2_Ny};
    adios2_define_attribute_array(ioH, "shape2D", adios2_type_uint16_t, shape2D, 2);
    adios2_define_attribute(ioH, "aI8", adios2_type_int8_t, data_I8);
    adios2_define_attribute(ioH, "aI16", adios2_type_int16_t, data_I16);
    adios2_define_attribute(ioH, "aI32", adios2_type_int32_t, data_I32);
    adios2_define_attribute(ioH, "aU8", adios2_type_uint8_t, data_U8);
    adios2_define_attribute(ioH, "aU16", adios2_type_uint16_t, data_U16);
    adios2_define_attribute(ioH, "aU32", adios2_type_uint32_t, data_U32);
    adios2_define_attribute(ioH, "aR32", adios2_type_float, data_R32);
    adios2_define_attribute(ioH, "aR64", adios2_type_double, data_R64);

    // inquire variables
    adios2_variable *varNproc = adios2_inquire_variable(ioH, "nproc");
    adios2_variable *varI8 = adios2_inquire_variable(ioH, "varI8");
    adios2_variable *varI16 = adios2_inquire_variable(ioH, "varI16");
    adios2_variable *varI32 = adios2_inquire_variable(ioH, "varI32");
    adios2_variable *varI64 = adios2_inquire_variable(ioH, "varI64");
    adios2_variable *varU8 = adios2_inquire_variable(ioH, "varU8");
    adios2_variable *varU16 = adios2_inquire_variable(ioH, "varU16");
    adios2_variable *varU32 = adios2_inquire_variable(ioH, "varU32");
    adios2_variable *varU64 = adios2_inquire_variable(ioH, "varU64");
    adios2_variable *varR32 = adios2_inquire_variable(ioH, "varR32");
    adios2_variable *varR64 = adios2_inquire_variable(ioH, "varR64");
    adios2_variable *R64_2d = adios2_inquire_variable(ioH, "R64_2d");

    adios2_engine *engineH = adios2_open(ioH, "TestUtilsCWriter.bp", adios2_mode_write);

    adios2_put(engineH, varNproc, &nproc, adios2_mode_deferred);
    adios2_put(engineH, varI8, data_I8, adios2_mode_deferred);
    adios2_put(engineH, varI16, data_I16, adios2_mode_deferred);
    adios2_put(engineH, varI32, data_I32, adios2_mode_deferred);
    adios2_put(engineH, varI64, data_I64, adios2_mode_deferred);

    adios2_put(engineH, varU8, data_U8, adios2_mode_deferred);
    adios2_put(engineH, varU16, data_U16, adios2_mode_deferred);
    adios2_put(engineH, varU32, data_U32, adios2_mode_deferred);
    adios2_put(engineH, varU64, data_U64, adios2_mode_deferred);

    adios2_put(engineH, varR32, data_R32, adios2_mode_deferred);
    adios2_put(engineH, varR64, data_R64, adios2_mode_deferred);

    adios2_put(engineH, R64_2d, d2_R64, adios2_mode_deferred);
    adios2_close(engineH);

    // deallocate adiosH
    adios2_finalize(adiosH);

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
