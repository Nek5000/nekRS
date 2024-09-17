/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <iostream> //std::cout
#include <numeric>  //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

class BPWriteReadMGARDMDR : public ::testing::TestWithParam<std::string>
{
public:
    BPWriteReadMGARDMDR() = default;
    virtual void SetUp(){};
    virtual void TearDown(){};
};

TEST_F(BPWriteReadMGARDMDR, BPWRMGARD1D)
{
    // Refactor a dataset with MDR, then
    // read back with various accuracies

    int mpiRank = 0, mpiSize = 1;
    const size_t Nx = 30000; // 100k minimum data size for MDR
    const size_t NSteps = 1;

    std::vector<float> r32s(Nx);
    std::vector<double> r64s(Nx);

    const double pi = 3.141592653589793238462643383279;
    const double twopi = 2 * pi;
    const double value = (2 * pi) / double(Nx);
    for (size_t x = 0; x < Nx; ++x)
    {
        r64s[x] = std::cos(value * x);
        r32s[x] = (float)r64s[x];
    }

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWRMGARDMDR1D_MPI.bp");
#else
    const std::string fname("BPWRMGARDMDR1D.bp");
#endif

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif
    {
        adios2::IO io = adios.DeclareIO("WriteIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
        const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
        const adios2::Dims count{Nx};

        auto var_r32 = io.DefineVariable<float>("r32", shape, start, count, adios2::ConstantDims);
        auto var_r64 = io.DefineVariable<double>("r64", shape, start, count, adios2::ConstantDims);
        auto var_raw64 =
            io.DefineVariable<double>("raw64", shape, start, count, adios2::ConstantDims);

        // add operations
        adios2::Operator mgardOp = adios.DefineOperator("mgardCompressor", adios2::ops::MDR);

        var_r32.AddOperation(mgardOp, {});
        var_r64.AddOperation(mgardOp, {});

        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            bpWriter.BeginStep();
            bpWriter.Put<float>("r32", r32s.data());
            bpWriter.Put<double>("r64", r64s.data());
            bpWriter.Put<double>("raw64", r64s.data());
            bpWriter.EndStep();
        }

        bpWriter.Close();
    }

    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }

        adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);
        size_t t = 0;
        while (bpReader.BeginStep() == adios2::StepStatus::OK)
        {
            auto var_r32 = io.InquireVariable<float>("r32");
            EXPECT_TRUE(var_r32);
            ASSERT_EQ(var_r32.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r32.Steps(), NSteps);
            ASSERT_EQ(var_r32.Shape()[0], mpiSize * Nx);

            auto var_r64 = io.InquireVariable<double>("r64");
            EXPECT_TRUE(var_r64);
            ASSERT_EQ(var_r64.ShapeID(), adios2::ShapeID::GlobalArray);
            ASSERT_EQ(var_r64.Steps(), NSteps);
            ASSERT_EQ(var_r64.Shape()[0], mpiSize * Nx);

            const adios2::Dims start{mpiRank * Nx};
            const adios2::Dims count{Nx};
            const adios2::Box<adios2::Dims> sel(start, count);
            var_r32.SetSelection(sel);
            var_r64.SetSelection(sel);

            std::vector<double> errors = {0.0, 0.0000001, 0.0001, 0.1};
            for (auto error : errors)
            {
                std::cout << "===== Read with error = " << error << std::endl;

                std::vector<float> read32s;
                std::vector<double> read64s;

                adios2::Accuracy accuracyRequested = {
                    error, std::numeric_limits<double>::infinity(), false};
                var_r32.SetAccuracy(accuracyRequested);
                var_r64.SetAccuracy(accuracyRequested);
                bpReader.Get(var_r32, read32s, adios2::Mode::Sync);
                bpReader.Get(var_r64, read64s, adios2::Mode::Sync);

                auto accuracyGot32 = var_r32.GetAccuracy();
                assert(accuracyGot32.error <=
                       std::max(accuracyRequested.error,
                                adios2::ops::mdr::FLOAT_ROUNDING_ERROR_LIMIT));
                assert(accuracyGot32.norm == accuracyRequested.norm);
                assert(accuracyGot32.relative == accuracyRequested.relative);

                auto accuracyGot64 = var_r64.GetAccuracy();
                assert(accuracyGot64.error <=
                       std::max(accuracyRequested.error,
                                adios2::ops::mdr::DOUBLE_ROUNDING_ERROR_LIMIT));
                assert(accuracyGot64.norm == accuracyRequested.norm);
                assert(accuracyGot64.relative == accuracyRequested.relative);

                double maxDiff = 0.0, relativeMaxDiff = 0.0;
                size_t maxDiffPos = 0;

                for (size_t i = 0; i < Nx; ++i)
                {
                    double diff = std::abs(r64s[i] - read64s[i]);
                    if (diff > maxDiff)
                    {
                        maxDiff = diff;
                        maxDiffPos = i;
                    }
                }

                auto r64s_Max = std::max_element(r64s.begin(), r64s.end());
                relativeMaxDiff = maxDiff / *r64s_Max;
                std::cout << "double array: Relative Max Diff " << relativeMaxDiff << " Max Diff "
                          << maxDiff << " requested error " << error << " result error "
                          << accuracyGot64.error << " max diff pos " << maxDiffPos << " orig value "
                          << r64s[maxDiffPos] << " read value " << read64s[maxDiffPos] << "\n";
                ASSERT_LE(maxDiff, accuracyGot64.error);

                for (size_t i = 0; i < Nx; ++i)
                {
                    double diff = std::abs(r32s[i] - read32s[i]);
                    if (diff > maxDiff)
                    {
                        maxDiff = diff;
                        maxDiffPos = i;
                    }
                }

                auto r32s_Max = std::max_element(r32s.begin(), r32s.end());
                relativeMaxDiff = maxDiff / *r32s_Max;

                std::cout << "float array: Relative Max Diff " << relativeMaxDiff << " Max Diff "
                          << maxDiff << " requested error " << error << " result error "
                          << accuracyGot32.error << " max diff pos " << maxDiffPos << " orig value "
                          << r32s[maxDiffPos] << " read value " << read32s[maxDiffPos] << "\n";
                ASSERT_LE(maxDiff, accuracyGot32.error);
            }
            bpReader.EndStep();
            ++t;
        }

        EXPECT_EQ(t, NSteps);

        bpReader.Close();
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
