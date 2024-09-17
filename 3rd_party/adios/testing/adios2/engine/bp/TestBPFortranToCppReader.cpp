/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPFortranToCppRead : public ::testing::Test
{
public:
    BPFortranToCppRead() = default;

    SmallTestData m_TestData;
};

TEST_F(BPFortranToCppRead, ADIOS2BPFortranToCppRead)
{

    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("FortranToCpp_MPI.bp");
#else
    const std::string fname("FortranToCpp.bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    adios2::IO io = adios.DeclareIO("ReadIO");
    if (!engineName.empty())
    {
        // io.SetEngine(engineName);
    }

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    int32_t inx;
    int step = 0;
    while (true)
    {
        auto status = bpReader.BeginStep();
        if (status != adios2::StepStatus::OK)
        {
            break;
        }
        const size_t currentStep = bpReader.CurrentStep();
        EXPECT_EQ(currentStep, static_cast<size_t>(step));

        // Global value
        if (step == 0)
        {
            auto var_GlobalValue = io.InquireVariable<int32_t>("inx");
            EXPECT_TRUE(var_GlobalValue);
            EXPECT_EQ(var_GlobalValue.ShapeID(), adios2::ShapeID::GlobalValue);
            EXPECT_EQ(var_GlobalValue.Shape().size(), 0);
            bpReader.Get(var_GlobalValue, inx, adios2::Mode::Sync);
            EXPECT_EQ(inx, 10);
        }

        std::vector<int64_t> a;

        // 1D local array
        auto var_localarray_1D = io.InquireVariable<int64_t>("localarray_1D");
        EXPECT_TRUE(var_localarray_1D);
        EXPECT_EQ(var_localarray_1D.ShapeID(), adios2::ShapeID::LocalArray);
        EXPECT_EQ(var_localarray_1D.Shape().size(), 0);
        var_localarray_1D.SetBlockSelection(mpiRank);
        a.clear();
        bpReader.Get(var_localarray_1D, a, adios2::Mode::Sync);
        EXPECT_EQ(a.size(), static_cast<size_t>(inx));
        for (int i = 0; i < inx; ++i)
        {
            EXPECT_DOUBLE_EQ(a[i], (mpiRank + 1) * 1000 + (step + 1) * 100 + (i + 1));
        }

        // 2D local array
        auto var_localarray_2D = io.InquireVariable<int64_t>("localarray_2D");
        EXPECT_TRUE(var_localarray_2D);
        EXPECT_EQ(var_localarray_2D.ShapeID(), adios2::ShapeID::LocalArray);
        EXPECT_EQ(var_localarray_2D.Shape().size(), 0);
        var_localarray_2D.SetBlockSelection(mpiRank);
        a.clear();
        bpReader.Get(var_localarray_2D, a, adios2::Mode::Sync);
        EXPECT_EQ(a.size(), static_cast<size_t>(inx * 2));
        for (int i = 0; i < inx; ++i)
        {
            EXPECT_DOUBLE_EQ(a[i], (mpiRank + 1) * 1000 + (step + 1) * 100 + (i + 1));
            EXPECT_DOUBLE_EQ(a[i + inx], (mpiRank + 1) * 1000 + (step + 1) * 100 + (inx + i + 1));
        }

        // 1D changing local array
        auto var_localarray_1D_changing = io.InquireVariable<int64_t>("localarray_1D_changing");
        EXPECT_TRUE(var_localarray_1D_changing);
        EXPECT_EQ(var_localarray_1D_changing.ShapeID(), adios2::ShapeID::LocalArray);
        EXPECT_EQ(var_localarray_1D_changing.Shape().size(), 0);
        var_localarray_1D_changing.SetBlockSelection(mpiRank);
        a.clear();
        bpReader.Get(var_localarray_1D_changing, a, adios2::Mode::Sync);
        EXPECT_EQ(a.size(), static_cast<size_t>(step + 1));
        for (int i = 0; i < step + 1; ++i)
        {
            EXPECT_DOUBLE_EQ(a[i], (mpiRank + 1) * 1000 + (step + 1) * 100 + (i + 1) + (step + 1));
        }

        bpReader.EndStep();
        ++step;
    }

    EXPECT_EQ(step, 3);
    bpReader.Close();
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
