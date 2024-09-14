/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <algorithm> //std::min_element, std::max_element
#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

class BPReadFBlocks : public ::testing::Test
{
public:
    BPReadFBlocks() = default;
};

// ADIOS2 Fortran BP write, ADIOS2 C++ read
TEST_F(BPReadFBlocks, FHeatMap2D)
{
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(;
#endif
    // read
    {
        adios2::IO io = adios.DeclareIO("ReadF");
        adios2::Engine bpReader = io.Open("HeatMap2D_f.bp", adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            adios2::Variable<float> var_r32 = io.InquireVariable<float>("temperatures_r4");

            EXPECT_TRUE(var_r32);
            EXPECT_EQ(var_r32.Shape().size(), 2);

            const std::vector<adios2::Variable<float>::Info> r32Blocks =
                bpReader.BlocksInfo(var_r32, bpReader.CurrentStep());

            EXPECT_EQ(r32Blocks.size(), static_cast<size_t>(mpiSize));

            for (size_t i = 0; i < r32Blocks.size(); ++i)
            {
                if (i == static_cast<size_t>(mpiRank))
                {
                    EXPECT_EQ(r32Blocks[i].Start[0],
                              r32Blocks[i].Count[0] * static_cast<size_t>(mpiRank));

                    EXPECT_EQ(r32Blocks[i].Start[1], 0);
                    EXPECT_EQ(r32Blocks[i].Count[1], var_r32.Shape()[1]);
                    EXPECT_TRUE(r32Blocks[i].IsReverseDims);
                }
            }
            bpReader.EndStep();
        }
    }
}

TEST_F(BPReadFBlocks, FHeatMap3D)
{
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    // read
    {
        adios2::IO io = adios.DeclareIO("ReadF");
        adios2::Engine bpReader = io.Open("HeatMap3D_f.bp", adios2::Mode::Read);

        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            adios2::Variable<float> var_r32 = io.InquireVariable<float>("temperatures_r4");

            EXPECT_TRUE(var_r32);
            EXPECT_EQ(var_r32.Shape().size(), 3);

            const std::vector<adios2::Variable<float>::Info> r32Blocks =
                bpReader.BlocksInfo(var_r32, bpReader.CurrentStep());

            EXPECT_EQ(r32Blocks.size(), static_cast<size_t>(mpiSize));

            for (size_t i = 0; i < r32Blocks.size(); ++i)
            {
                if (i == static_cast<size_t>(mpiRank))
                {
                    EXPECT_EQ(r32Blocks[i].Start[0],
                              r32Blocks[i].Count[0] * static_cast<size_t>(mpiRank));

                    EXPECT_EQ(r32Blocks[i].Start[1], 0);
                    EXPECT_EQ(r32Blocks[i].Count[1], var_r32.Shape()[1]);

                    EXPECT_EQ(r32Blocks[i].Start[2], 0);
                    EXPECT_EQ(r32Blocks[i].Count[2], var_r32.Shape()[2]);
                    EXPECT_EQ(r32Blocks[i].IsReverseDims, true);
                }
            }
            bpReader.EndStep();
        }
    }
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
