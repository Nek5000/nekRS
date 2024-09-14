/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TestBPWriteAggregateReadLocal.cpp
 *
 *  Created on: Jan 11, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include <adios2_c.h>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

#include <gtest/gtest.h>

#include <numeric> //std::iota
#include <thread>

void LocalAggregate1D(const std::string substreams)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
#if ADIOS2_USE_MPI
    const std::string fname("LocalAggregate1D_" + substreams + "_MPI.bp");
#else
    const std::string fname("LocalAggregate1D_" + substreams + ".bp");
#endif
    int mpiRank = 0, mpiSize = 1;
    // Number of steps
    constexpr size_t NSteps = 5;
    constexpr size_t Nx = 20;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    std::vector<int32_t> inumbers(NSteps * Nx);
    std::vector<float> fnumbers(NSteps * Nx);

    adios2_adios *adiosH = adios2_init_mpi(MPI_COMM_WORLD);

    // writer
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "TestIO");

        size_t count[1];
        count[0] = Nx;

        adios2_variable *variNumbers = adios2_define_variable(
            ioH, "ints", adios2_type_int32_t, 1, NULL, NULL, count, adios2_constant_dims_true);

        adios2_variable *varfNumbers = adios2_define_variable(
            ioH, "floats", adios2_type_float, 1, NULL, NULL, count, adios2_constant_dims_true);

        // TODO adios2_set_parameter(ioH, "CollectiveMetadata", "Off");
        adios2_set_parameter(ioH, "Profile", "Off");

        if (mpiSize > 1)
        {
            adios2_set_parameter(ioH, "substreams", substreams.c_str());
        }

        adios2_engine *bpWriter = adios2_open(ioH, fname.c_str(), adios2_mode_write);

        adios2_step_status step_status;
        for (size_t i = 0; i < NSteps; ++i)
        {
            adios2_begin_step(bpWriter, adios2_step_mode_read, -1., &step_status);

            std::iota(inumbers.begin() + i * Nx, inumbers.begin() + i * Nx + Nx, mpiRank);

            const float randomStart = static_cast<float>(rand() % mpiSize);
            std::iota(fnumbers.begin() + i * Nx, fnumbers.begin() + i * Nx + Nx, randomStart);

            //            if (mpiRank % 3 == 0)
            //            {
            adios2_put(bpWriter, variNumbers, &inumbers[i * Nx], adios2_mode_sync);
            //}

            //            if (mpiRank % 3 == 1)
            //            {
            adios2_put(bpWriter, varfNumbers, &fnumbers[i * Nx], adios2_mode_sync);
            //}
            adios2_end_step(bpWriter);
        }

        adios2_close(bpWriter);
        bpWriter = NULL;
    }

    // Reader TODO, might need to generate metadata file
    // if (false)
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "Reader");
        adios2_engine *bpReader = adios2_open(ioH, fname.c_str(), adios2_mode_read);

        adios2_step_status step_status;
        while (true)
        {
            adios2_begin_step(bpReader, adios2_step_mode_read, -1, &step_status);

            if (step_status == adios2_step_status_end_of_stream)
            {
                break;
            }

            size_t currentStep;
            adios2_current_step(&currentStep, bpReader);

            adios2_variable *varInts = adios2_inquire_variable(ioH, "ints");
            adios2_variable *varFloats = adios2_inquire_variable(ioH, "floats");

            //            if (mpiRank % 3 == 0)
            //            {
            EXPECT_NE(varInts, nullptr);
            std::vector<int32_t> inVarInts(Nx);

            adios2_set_block_selection(varInts, mpiRank);
            adios2_get(bpReader, varInts, inVarInts.data(), adios2_mode_sync);

            for (size_t i = 0; i < Nx; ++i)
            {
                ASSERT_EQ(inVarInts[i], inumbers[currentStep * Nx + i]);
            }
            //}

            //            if (mpiRank % 3 == 1)
            //            {
            EXPECT_NE(varFloats, nullptr);
            std::vector<float> inVarFloats(Nx);

            adios2_set_block_selection(varFloats, mpiRank);
            adios2_get(bpReader, varFloats, inVarFloats.data(), adios2_mode_sync);

            for (size_t i = 0; i < Nx; ++i)
            {
                ASSERT_EQ(inVarFloats[i], fnumbers[currentStep * Nx + i]);
            }
            //}

            adios2_end_step(bpReader);
        }

        adios2_close(bpReader);
        bpReader = NULL;
    }

    adios2_finalize(adiosH);
}

void LocalAggregate1DBlock0(const std::string substreams)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
#if ADIOS2_USE_MPI
    const std::string fname("LocalAggregate1DSubFile_" + substreams + "_MPI.bp");
#else
    const std::string fname("LocalAggregate1DSubFile_" + substreams + ".bp");
#endif
    int mpiRank = 0, mpiSize = 1;
    // Number of steps
    constexpr size_t NSteps = 5;
    constexpr size_t Nx = 20;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

    std::vector<int32_t> inumbers(NSteps * Nx);
    std::vector<float> fnumbers(NSteps * Nx);

    adios2_adios *adiosH = adios2_init_mpi(MPI_COMM_WORLD);

    // writer
    {
        adios2_io *ioH = adios2_declare_io(adiosH, "TestIO");
        adios2_set_engine(ioH, "BP3");

        size_t count[1];
        count[0] = Nx;

        adios2_variable *variNumbers = adios2_define_variable(
            ioH, "ints", adios2_type_int32_t, 1, NULL, NULL, count, adios2_constant_dims_true);

        adios2_variable *varfNumbers = adios2_define_variable(
            ioH, "floats", adios2_type_float, 1, NULL, NULL, count, adios2_constant_dims_true);

        // adios2_set_parameter(ioH, "CollectiveMetadata", "Off");
        adios2_set_parameter(ioH, "Profile", "Off");

        if (mpiSize > 1)
        {
            adios2_set_parameter(ioH, "substreams", substreams.c_str());
        }

        adios2_engine *bpWriter = adios2_open(ioH, fname.c_str(), adios2_mode_write);

        adios2_step_status step_status;
        for (size_t i = 0; i < NSteps; ++i)
        {
            adios2_begin_step(bpWriter, adios2_step_mode_read, -1., &step_status);

            std::iota(inumbers.begin() + i * Nx, inumbers.begin() + i * Nx + Nx, mpiRank);

            const float randomStart = static_cast<float>(rand() % mpiSize);
            std::iota(fnumbers.begin() + i * Nx, fnumbers.begin() + i * Nx + Nx, randomStart);

            //            if (mpiRank % 3 == 0)
            //            {
            adios2_put(bpWriter, variNumbers, &inumbers[i * Nx], adios2_mode_sync);
            //}

            //            if (mpiRank % 3 == 1)
            //            {
            adios2_put(bpWriter, varfNumbers, &fnumbers[i * Nx], adios2_mode_sync);
            //}
            adios2_end_step(bpWriter);
        }

        adios2_close(bpWriter);
        bpWriter = NULL;
    }

    {
        // read block zero from collective, then compare to subfile read
        adios2_io *ioH = adios2_declare_io(adiosH, "Reader");
        adios2_set_engine(ioH, "File");

        adios2_engine *bpReader = adios2_open(ioH, fname.c_str(), adios2_mode_read);

        // subfile read
        adios2_io *ioH0 = adios2_declare_io(adiosH, "Reader0");
        adios2_set_engine(ioH0, "File");

        const std::string fnameBP0 = fname + ".dir/" + fname + ".0";

        adios2_engine *bpReader0 = adios2_open(ioH0, fnameBP0.c_str(), adios2_mode_read);

        adios2_step_status step_status;
        while (true)
        {
            adios2_begin_step(bpReader, adios2_step_mode_read, -1, &step_status);
            adios2_begin_step(bpReader0, adios2_step_mode_read, -1, &step_status);

            if (step_status == adios2_step_status_end_of_stream)
            {
                break;
            }

            size_t currentStep;
            adios2_current_step(&currentStep, bpReader0);

            adios2_variable *varInts = adios2_inquire_variable(ioH, "ints");
            adios2_variable *varFloats = adios2_inquire_variable(ioH, "floats");

            adios2_variable *varInts0 = adios2_inquire_variable(ioH0, "ints");
            adios2_variable *varFloats0 = adios2_inquire_variable(ioH0, "floats");

            EXPECT_NE(varInts, nullptr);
            EXPECT_NE(varInts0, nullptr);

            std::vector<int32_t> inVarInts(Nx);
            std::vector<int32_t> inVarInts0(Nx);

            adios2_set_block_selection(varInts, 0);
            adios2_get(bpReader, varInts, inVarInts.data(), adios2_mode_sync);

            adios2_set_block_selection(varInts0, 0);
            adios2_get(bpReader0, varInts0, inVarInts0.data(), adios2_mode_sync);

            for (size_t i = 0; i < Nx; ++i)
            {
                ASSERT_EQ(inVarInts0[i], inVarInts[i]);
            }

            if (mpiRank == 0)
            {
                for (size_t i = 0; i < Nx; ++i)
                {
                    ASSERT_EQ(inVarInts0[i], inumbers[currentStep * Nx + i]);
                }
            }

            // floats
            EXPECT_NE(varFloats, nullptr);
            EXPECT_NE(varFloats0, nullptr);

            std::vector<float> inVarFloats(Nx);
            std::vector<float> inVarFloats0(Nx);

            adios2_set_block_selection(varFloats, 0);
            adios2_get(bpReader, varFloats, inVarFloats.data(), adios2_mode_sync);

            adios2_set_block_selection(varFloats0, 0);
            adios2_get(bpReader0, varFloats0, inVarFloats0.data(), adios2_mode_sync);

            for (size_t i = 0; i < Nx; ++i)
            {
                ASSERT_EQ(inVarFloats[i], inVarFloats[i]);
            }

            if (mpiRank == 0)
            {
                for (size_t i = 0; i < Nx; ++i)
                {
                    ASSERT_EQ(inVarFloats0[i], fnumbers[currentStep * Nx + i]);
                }
            }

            adios2_end_step(bpReader0);
            adios2_end_step(bpReader);
        }

        adios2_close(bpReader0);
        adios2_close(bpReader);
        bpReader0 = NULL;
        bpReader = NULL;
    }

    adios2_finalize(adiosH);
}

class BPWriteAggregateReadLocalTest : public ::testing::TestWithParam<std::string>
{
public:
    BPWriteAggregateReadLocalTest() = default;

    virtual void SetUp() {}
    virtual void TearDown() {}
};

TEST_P(BPWriteAggregateReadLocalTest, Aggregate1D) { LocalAggregate1D(GetParam()); }

TEST_P(BPWriteAggregateReadLocalTest, Aggregate1DBlock0) { LocalAggregate1DBlock0(GetParam()); }

INSTANTIATE_TEST_SUITE_P(Substreams, BPWriteAggregateReadLocalTest,
                         ::testing::Values("1", "2", "3", "4"));

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
