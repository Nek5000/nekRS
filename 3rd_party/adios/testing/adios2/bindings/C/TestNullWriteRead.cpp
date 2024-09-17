/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */

#include <adios2_c.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include <cstring>
#include <gtest/gtest.h>

class NullWriteReadTests_C_API : public ::testing::Test
{
public:
    NullWriteReadTests_C_API() = default;
};

//******************************************************************************
// 1D 1x8 test data
//******************************************************************************

// ADIOS2 BP write, native ADIOS1 read
TEST_F(NullWriteReadTests_C_API, NullWriteRead1D8)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 3;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using BP

#if ADIOS2_USE_MPI
    adios2_adios *adios = adios2_init_mpi(MPI_COMM_WORLD);
    const std::string fname("NullWriteRead1D8_c_MPI.bp");
#else
    adios2_adios *adios = adios2_init_serial();
    const std::string fname("NullWriteRead1D8_c.bp");
#endif
    {
        adios2_io *io = adios2_declare_io(adios, "WriteNull");
        const std::vector<std::size_t> shape{static_cast<size_t>(Nx * mpiSize)};
        const std::vector<std::size_t> start{static_cast<size_t>(Nx * mpiRank)};
        const std::vector<std::size_t> count{Nx};
        std::vector<double> data(Nx, 1.0);

        adios2_variable *var =
            adios2_define_variable(io, "r64", adios2_type_double, 1, shape.data(), start.data(),
                                   count.data(), adios2_constant_dims_true);
        adios2_set_engine(io, "NULL");

        adios2_engine *nullWriter = adios2_open(io, fname.data(), adios2_mode_write);

        adios2_step_status status;

        for (size_t step = 0; step < NSteps; ++step)
        {
            adios2_begin_step(nullWriter, adios2_step_mode_append, -1., &status);
            adios2_put(nullWriter, var, data.data(), adios2_mode_deferred);
            adios2_perform_puts(nullWriter);
            adios2_end_step(nullWriter);
        }

        // Close the file
        adios2_close(nullWriter);
    }

    {
        adios2_io *io = adios2_declare_io(adios, "ReadNull");
        adios2_set_engine(io, "null");

        adios2_engine *nullReader = adios2_open(io, fname.data(), adios2_mode_read);

        std::vector<double> R64;

        adios2_step_status status;
        size_t currentStep;

        for (size_t t = 0; t < NSteps; ++t)
        {
            adios2_begin_step(nullReader, adios2_step_mode_read, -1., &status);
            EXPECT_EQ(status, adios2_step_status_end_of_stream);

            adios2_variable *var = adios2_inquire_variable(io, "r64");
            EXPECT_EQ(var, nullptr);

            adios2_current_step(&currentStep, nullReader);
            EXPECT_EQ(currentStep, t);

            auto ret = adios2_get(nullReader, var, nullptr, adios2_mode_deferred);
            EXPECT_EQ(ret, 1);
            adios2_perform_gets(nullReader);
            adios2_end_step(nullReader);
        }
        adios2_close(nullReader);
    }
    adios2_finalize(adios);
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
