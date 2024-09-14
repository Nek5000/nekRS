/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>
#include <tuple>

#include <adios2.h>
#include <adios2/common/ADIOSConfig.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

using ParamType = std::tuple<std::string, adios2::Params, double>;

class BPChangingShapeWithinStep : public ::testing::TestWithParam<ParamType>
{
public:
    BPChangingShapeWithinStep() = default;
    virtual void SetUp(){};
    virtual void TearDown(){};
};

TEST_P(BPChangingShapeWithinStep, MultiBlock)
{
    // Write multiple blocks and change shape in between
    // At read, the last shape should be used not the first one
    // This test guarantees that one can change the variable shape
    // until EndStep()

    auto operatorName = std::get<0>(GetParam());
    auto params = std::get<1>(GetParam());
    double epsilon = std::get<2>(GetParam());

    const int nsteps = 2;
    const std::vector<int> nblocks = {2, 3};
    const int N = 16384; // size of one block (should be big enough to compress)
    EXPECT_EQ(nsteps, nblocks.size());
    int rank = 0, nproc = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    const std::string fname("BPChangingShapeMultiblock_" + operatorName + "_MPI.bp");
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    const std::string fname("BPChangingShapeMultiblock_" + operatorName + ".bp");
    adios2::ADIOS adios;
#endif

    // Writer
    {
        adios2::IO outIO = adios.DeclareIO("Output");

        if (!engineName.empty())
        {
            outIO.SetEngine(engineName);
        }

        adios2::Engine writer = outIO.Open(fname, adios2::Mode::Write);

        const size_t dim0 = static_cast<size_t>(nproc);
        const size_t off0 = static_cast<size_t>(rank);
        auto var = outIO.DefineVariable<double>("v", {dim0, 1}, {off0, 0}, {1, 1});

        if (operatorName != "none")
        {
            auto op = adios.DefineOperator("compressor", operatorName, params);
            var.AddOperation(op);
        }

        if (!rank)
        {
            std::cout << "Writing to :" << fname;
            if (operatorName != "none")
            {
                std::cout << " with operator " << operatorName;
            }
            std::cout << std::endl;
        }
        for (size_t i = 0; i < nsteps; i++)
        {
            writer.BeginStep();

            double value = static_cast<double>(rank) + static_cast<double>(i + 1) / 10.0;

            for (size_t j = 0; j < static_cast<size_t>(nblocks[i]); j++)
            {
                std::vector<double> data(N, value);
                var.SetShape({dim0, (j + 1) * N});
                var.SetSelection({{off0, j * N}, {1, N}});

                if (!rank)
                {
                    std::cout << "Step " << i << " block " << j << " shape (" << var.Shape()[0]
                              << ", " << var.Shape()[1] << ")"
                              << " value = " << value << " data[] = " << data[0] << " .. "
                              << data[N - 1] << std::endl;
                }

                writer.Put(var, data.data(), adios2::Mode::Sync);
                value += 0.01;
            }
            writer.EndStep();
        }
        writer.Close();
    }

    // Reader with streaming
    {
        adios2::IO inIO = adios.DeclareIO("Input");

        if (!engineName.empty())
        {
            inIO.SetEngine(engineName);
        }
        adios2::Engine reader = inIO.Open(fname, adios2::Mode::Read);

        if (!rank)
        {
            std::cout << "Reading as stream with BeginStep/EndStep:" << std::endl;
        }

        int step = 0;
        while (true)
        {
            adios2::StepStatus status = reader.BeginStep(adios2::StepMode::Read);

            if (status != adios2::StepStatus::OK)
            {
                break;
            }

            size_t expected_shape = N * nblocks[step];

            auto var = inIO.InquireVariable<double>("v");
            EXPECT_TRUE(var);

            if (!rank)
            {

                std::cout << "Step " << step << " shape (" << var.Shape()[0] << ", "
                          << var.Shape()[1] << ")" << std::endl;
            }

            EXPECT_EQ(var.Shape()[0], nproc);
            EXPECT_EQ(var.Shape()[1], expected_shape);

            var.SetSelection({{0, 0}, {static_cast<size_t>(nproc), expected_shape}});

            // Check data on rank 0
            if (!rank)
            {
                std::vector<double> data(nproc * expected_shape);
                reader.Get(var, data.data());

                reader.PerformGets();

                for (int i = 0; i < nproc; i++)
                {
                    double value = static_cast<double>(i) + static_cast<double>(step + 1) / 10.0;

                    for (int j = 0; j < nblocks[step]; j++)
                    {
                        EXPECT_LE(fabs(data[(i * nblocks[step] + j) * N] - value), epsilon);
                        value += 0.01;
                    }
                }
            }

            reader.EndStep();
            ++step;
        }
        reader.Close();
    }

    // Reader with file reading
    {
        adios2::IO inIO = adios.DeclareIO("InputFile");

        if (!engineName.empty())
        {
            inIO.SetEngine(engineName);
        }
        adios2::Engine reader = inIO.Open(fname, adios2::Mode::ReadRandomAccess);

        if (!rank)
        {
            std::cout << "Reading as file with SetStepSelection:" << std::endl;
        }

        auto var = inIO.InquireVariable<double>("v");
        EXPECT_TRUE(var);
        for (int step = 0; step < nsteps; step++)
        {
            var.SetStepSelection({step, 1});
            if (!rank)
            {
                std::cout << "Step " << step << " shape (" << var.Shape()[0] << ", "
                          << var.Shape()[1] << ")" << std::endl;
            }
            size_t expected_shape = N * nblocks[step];
            EXPECT_EQ(var.Shape()[0], nproc);
            EXPECT_EQ(var.Shape()[1], expected_shape);

            var.SetSelection({{0, 0}, {static_cast<size_t>(nproc), expected_shape}});

            std::vector<double> data(nproc * expected_shape);
            reader.Get(var, data.data());
            reader.PerformGets();

            for (int i = 0; i < nproc; i++)
            {
                double value = static_cast<double>(i) + static_cast<double>(step + 1) / 10.0;

                for (int j = 0; j < nblocks[step]; j++)
                {
                    EXPECT_LE(fabs(data[(i * nblocks[step] + j) * N] - value), epsilon);
                    value += 0.01;
                }
            }

            EXPECT_THROW(reader.EndStep(), std::logic_error);
        }
    }
}

adios2::Params paccuracy = {{"accuracy", "0.001"}};
adios2::Params pblosc = {{"clevel", "9"}};

std::vector<ParamType> p = {{"none", paccuracy, 0.0}
#ifdef ADIOS2_HAVE_BLOSC2
                            ,
                            {"blosc", pblosc, 0.0}
#endif
#ifdef ADIOS2_HAVE_MGARD
                            ,
                            {"mgard", paccuracy, 0.001}
#endif
#ifdef ADIOS2_HAVE_ZFP
#ifndef ADIOS2_HAVE_ZFP_CUDA // only test on CPU
                            ,
                            {"zfp", paccuracy, 0.001}
#endif
#endif

};

INSTANTIATE_TEST_SUITE_P(Multiblock, BPChangingShapeWithinStep, ::testing::ValuesIn(p));

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
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
