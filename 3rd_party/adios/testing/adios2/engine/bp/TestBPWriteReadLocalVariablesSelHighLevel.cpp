/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#include "../SmallTestData.h"

std::string engineName; // comes from command line

class BPWriteReadLocalVariablesSelHighLevel : public ::testing::Test
{
public:
    BPWriteReadLocalVariablesSelHighLevel() = default;

    SmallTestData m_TestData;
};

TEST_F(BPWriteReadLocalVariablesSelHighLevel, BPWriteReadLocal1DSel)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadLocal1DSelHighLevel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal1DSelHighLevel.bp");
#endif

    // Write test data using BP
    // writer
    {
        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Nx};

#if ADIOS2_USE_MPI
        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t step32 = static_cast<int32_t>(step);
            // global value
            oStream.write("stepsGlobalValue", step32, adios2::GlobalValue);
            oStream.write("stepsGlobalValueString", std::to_string(step), adios2::GlobalValue);

            // local value
            const int32_t localValue = static_cast<int32_t>(mpiRank + step);
            oStream.write("ranksLocalValue", localValue, adios2::LocalValue);
            oStream.write("ranksLocalValueString", std::to_string(localValue), adios2::LocalValue);

            // local arrays
            oStream.write("i8", currentTestData.I8.data(), shape, start, count);
            oStream.write("i16", currentTestData.I16.data(), shape, start, count);
            oStream.write("i32", currentTestData.I32.data(), shape, start, count);
            oStream.write("i64", currentTestData.I64.data(), shape, start, count);
            oStream.write("u8", currentTestData.U8.data(), shape, start, count);
            oStream.write("u16", currentTestData.U16.data(), shape, start, count);
            oStream.write("u32", currentTestData.U32.data(), shape, start, count);
            oStream.write("u64", currentTestData.U64.data(), shape, start, count);
            oStream.write("r32", currentTestData.R32.data(), shape, start, count);
            oStream.write("r64", currentTestData.R64.data(), shape, start, count);
            oStream.write("cr32", currentTestData.CR32.data(), shape, start, count);
            oStream.write("cr64", currentTestData.CR64.data(), shape, start, count);
            oStream.end_step();
        }
        oStream.close();

#if ADIOS2_USE_MPI
        adios2::fstream iStream(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream iStream(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        for (adios2::fstep iStep; adios2::getstep(iStream, iStep);)
        {
            const size_t currentStep = iStep.current_step();
            EXPECT_EQ(currentStep, t);

            const std::vector<int32_t> stepsGlobalValue = iStep.read<int32_t>("stepsGlobalValue");
            EXPECT_EQ(stepsGlobalValue.front(), t);

            const std::vector<std::string> stepsGlobalValueString =
                iStep.read<std::string>("stepsGlobalValueString");
            EXPECT_EQ(stepsGlobalValueString.front(), std::to_string(t));

            const std::vector<int32_t> ranksGlobalArray = iStep.read<int32_t>("ranksLocalValue");

            const std::vector<std::string> ranksGlobalArrayString =
                iStep.read<std::string>("ranksLocalValueString");

            // loop blocks
            for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
            {
                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                EXPECT_EQ(ranksGlobalArray[b], b + t);
                EXPECT_EQ(ranksGlobalArrayString[b], std::to_string(b + t));

                // loop selections
                for (size_t s = 0; s < Nx; ++s)
                {
                    const adios2::Dims start = {s};
                    const adios2::Dims count = {Nx - s};

                    const auto I8 = iStep.read<int8_t>("i8", start, count, b);
                    const auto I16 = iStep.read<int16_t>("i16", start, count, b);
                    const auto I32 = iStep.read<int32_t>("i32", start, count, b);
                    const auto I64 = iStep.read<int64_t>("i64", start, count, b);
                    const auto U8 = iStep.read<uint8_t>("u8", start, count, b);
                    const auto U16 = iStep.read<uint16_t>("u16", start, count, b);
                    const auto U32 = iStep.read<uint32_t>("u32", start, count, b);
                    const auto U64 = iStep.read<uint64_t>("u64", start, count, b);
                    const auto R32 = iStep.read<float>("r32", start, count, b);
                    const auto R64 = iStep.read<double>("r64", start, count, b);

                    const auto CR32 = iStep.read<std::complex<float>>("cr32", start, count, b);
                    const auto CR64 = iStep.read<std::complex<double>>("cr64", start, count, b);

                    for (size_t i = s; i < Nx; ++i)
                    {
                        std::stringstream ss;
                        ss << "t=" << t << " s=" << s << " rank=" << mpiRank;
                        std::string msg = ss.str();

                        EXPECT_EQ(I8[i - s], currentTestData.I8[i]) << msg;
                        EXPECT_EQ(I16[i - s], currentTestData.I16[i]) << msg;
                        EXPECT_EQ(I32[i - s], currentTestData.I32[i]) << msg;
                        EXPECT_EQ(I64[i - s], currentTestData.I64[i]) << msg;
                        EXPECT_EQ(U8[i - s], currentTestData.U8[i]) << msg;
                        EXPECT_EQ(U16[i - s], currentTestData.U16[i]) << msg;
                        EXPECT_EQ(U32[i - s], currentTestData.U32[i]) << msg;
                        EXPECT_EQ(U64[i - s], currentTestData.U64[i]) << msg;
                        EXPECT_EQ(R32[i - s], currentTestData.R32[i]) << msg;
                        EXPECT_EQ(R64[i - s], currentTestData.R64[i]) << msg;
                        EXPECT_EQ(CR32[i - s], currentTestData.CR32[i]) << msg;
                        EXPECT_EQ(CR64[i - s], currentTestData.CR64[i]) << msg;
                    } // element
                }     // selection
            }         // blockID
            ++t;
        } // step
    }
}

TEST_F(BPWriteReadLocalVariablesSelHighLevel, BPWriteReadLocal2D2x4Sel)
{

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 4;
    const size_t Ny = 2;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadLocal2D2x4SelHighLevel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal2D2x4SelHighLevel.bp");
#endif

    // Write test data using BP
    // writer
    {
        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Ny, Nx};

#if ADIOS2_USE_MPI
        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t step32 = static_cast<int32_t>(step);
            // global value
            oStream.write("stepsGlobalValue", step32, adios2::GlobalValue);
            oStream.write("stepsGlobalValueString", std::to_string(step), adios2::GlobalValue);

            // local value
            const int32_t localValue = static_cast<int32_t>(mpiRank + step);
            oStream.write("ranksLocalValue", localValue, adios2::LocalValue);
            oStream.write("ranksLocalValueString", std::to_string(localValue), adios2::LocalValue);

            // local arrays
            oStream.write("i8", currentTestData.I8.data(), shape, start, count);
            oStream.write("i16", currentTestData.I16.data(), shape, start, count);
            oStream.write("i32", currentTestData.I32.data(), shape, start, count);
            oStream.write("i64", currentTestData.I64.data(), shape, start, count);
            oStream.write("u8", currentTestData.U8.data(), shape, start, count);
            oStream.write("u16", currentTestData.U16.data(), shape, start, count);
            oStream.write("u32", currentTestData.U32.data(), shape, start, count);
            oStream.write("u64", currentTestData.U64.data(), shape, start, count);
            oStream.write("r32", currentTestData.R32.data(), shape, start, count);
            oStream.write("r64", currentTestData.R64.data(), shape, start, count);
            oStream.write("cr32", currentTestData.CR32.data(), shape, start, count);
            oStream.write("cr64", currentTestData.CR64.data(), shape, start, count);
            oStream.end_step();
        }
        oStream.close();

#if ADIOS2_USE_MPI
        adios2::fstream iStream(fname, adios2::fstream::in, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream iStream(fname, adios2::fstream::in, engineName);
#endif

        size_t t = 0;
        for (adios2::fstep iStep; adios2::getstep(iStream, iStep);)
        {
            const size_t currentStep = iStep.current_step();
            EXPECT_EQ(currentStep, t);

            const std::vector<int32_t> stepsGlobalValue = iStep.read<int32_t>("stepsGlobalValue");
            EXPECT_EQ(stepsGlobalValue.front(), t);

            const std::vector<std::string> stepsGlobalValueString =
                iStep.read<std::string>("stepsGlobalValueString");
            EXPECT_EQ(stepsGlobalValueString.front(), std::to_string(t));

            const std::vector<int32_t> ranksGlobalArray = iStep.read<int32_t>("ranksLocalValue");

            const std::vector<std::string> ranksGlobalArrayString =
                iStep.read<std::string>("ranksLocalValueString");

            // loop blocks
            for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
            {
                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(currentStep), static_cast<int>(b), mpiSize);

                EXPECT_EQ(ranksGlobalArray[b], b + t);
                EXPECT_EQ(ranksGlobalArrayString[b], std::to_string(b + t));

                // loop selections
                for (size_t j = 0; j < Ny; ++j)
                {
                    for (size_t i = 0; i < Nx; ++i)
                    {
                        const adios2::Dims start = {j, i};
                        const adios2::Dims count = {Ny - j, Nx - i};

                        const auto I8 = iStep.read<int8_t>("i8", start, count, b);
                        const auto I16 = iStep.read<int16_t>("i16", start, count, b);
                        const auto I32 = iStep.read<int32_t>("i32", start, count, b);
                        const auto I64 = iStep.read<int64_t>("i64", start, count, b);
                        const auto U8 = iStep.read<uint8_t>("u8", start, count, b);
                        const auto U16 = iStep.read<uint16_t>("u16", start, count, b);
                        const auto U32 = iStep.read<uint32_t>("u32", start, count, b);
                        const auto U64 = iStep.read<uint64_t>("u64", start, count, b);
                        const auto R32 = iStep.read<float>("r32", start, count, b);
                        const auto R64 = iStep.read<double>("r64", start, count, b);

                        const auto CR32 = iStep.read<std::complex<float>>("cr32", start, count, b);
                        const auto CR64 = iStep.read<std::complex<double>>("cr64", start, count, b);

                        // element loop
                        for (size_t q = j; q < Ny; ++q)
                        {
                            for (size_t p = i; p < Nx; ++p)
                            {
                                std::stringstream ss;
                                ss << "t=" << t << " q=" << q << " p=" << p << " rank=" << mpiRank;
                                std::string msg = ss.str();

                                const size_t indexSel = (q - j) * (Nx - i) + (p - i);
                                const size_t indexBlock = q * Nx + p;

                                EXPECT_EQ(I8[indexSel], currentTestData.I8[indexBlock]) << msg;

                                EXPECT_EQ(I16[indexSel], currentTestData.I16[indexBlock]) << msg;

                                EXPECT_EQ(I32[indexSel], currentTestData.I32[indexBlock]) << msg;

                                EXPECT_EQ(I64[indexSel], currentTestData.I64[indexBlock]) << msg;

                                EXPECT_EQ(U8[indexSel], currentTestData.U8[indexBlock]) << msg;

                                EXPECT_EQ(U16[indexSel], currentTestData.U16[indexBlock]) << msg;

                                EXPECT_EQ(U32[indexSel], currentTestData.U32[indexBlock]) << msg;

                                EXPECT_EQ(U64[indexSel], currentTestData.U64[indexBlock]) << msg;

                                EXPECT_EQ(R32[indexSel], currentTestData.R32[indexBlock]) << msg;

                                EXPECT_EQ(R64[indexSel], currentTestData.R64[indexBlock]) << msg;

                                EXPECT_EQ(CR32[indexSel], currentTestData.CR32[indexBlock]) << msg;

                                EXPECT_EQ(CR64[indexSel], currentTestData.CR64[indexBlock]) << msg;
                            } // q index
                        }     // p index
                    }         // j sel
                }             // i sel
            }                 // blockID
            ++t;
        } // step
    }
}

TEST_F(BPWriteReadLocalVariablesSelHighLevel, BPWriteReadLocal1DAllStepsSel)
{
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

    int mpiRank = 0, mpiSize = 1;
    // Number of rows
    const size_t Nx = 8;

    // Number of steps
    const size_t NSteps = 5;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    const std::string fname("BPWriteReadLocal1DAllStepsSelHighLevel_MPI.bp");
#else
    const std::string fname("BPWriteReadLocal1DAllStepsSelHighLevel.bp");
#endif

    // Write test data using BP
    // writer
    {
        const adios2::Dims shape{};
        const adios2::Dims start{};
        const adios2::Dims count{Nx};

#if ADIOS2_USE_MPI
        adios2::fstream oStream(fname, adios2::fstream::out, MPI_COMM_WORLD, engineName);
#else
        adios2::fstream oStream(fname, adios2::fstream::out, engineName);
#endif

        for (size_t step = 0; step < NSteps; ++step)
        {
            SmallTestData currentTestData =
                generateNewSmallTestData(m_TestData, static_cast<int>(step), mpiRank, mpiSize);

            const int32_t step32 = static_cast<int32_t>(step);
            // global value
            oStream.write("stepsGlobalValue", step32, adios2::GlobalValue);
            oStream.write("stepsGlobalValueString", std::to_string(step), adios2::GlobalValue);

            // local value
            const int32_t localValue = static_cast<int32_t>(mpiRank + step);
            oStream.write("ranksLocalValue", localValue, adios2::LocalValue);
            oStream.write("ranksLocalValueString", std::to_string(localValue), adios2::LocalValue);

            // local arrays
            oStream.write("i8", currentTestData.I8.data(), shape, start, count);
            oStream.write("i16", currentTestData.I16.data(), shape, start, count);
            oStream.write("i32", currentTestData.I32.data(), shape, start, count);
            oStream.write("i64", currentTestData.I64.data(), shape, start, count);
            oStream.write("u8", currentTestData.U8.data(), shape, start, count);
            oStream.write("u16", currentTestData.U16.data(), shape, start, count);
            oStream.write("u32", currentTestData.U32.data(), shape, start, count);
            oStream.write("u64", currentTestData.U64.data(), shape, start, count);
            oStream.write("r32", currentTestData.R32.data(), shape, start, count);
            oStream.write("r64", currentTestData.R64.data(), shape, start, count);
            oStream.write("cr32", currentTestData.CR32.data(), shape, start, count);
            oStream.write("cr64", currentTestData.CR64.data(), shape, start, count);
            oStream.end_step();
        }
        oStream.close();
    }
#if ADIOS2_USE_MPI
    adios2::fstream iStream(fname, adios2::fstream::in_random_access, MPI_COMM_WORLD, engineName);
#else
    adios2::fstream iStream(fname, adios2::fstream::in_random_access, engineName);
#endif

    const size_t stepStart = 0;
    const size_t stepCount = NSteps;

    const std::vector<int32_t> stepsGlobalValue =
        iStream.read<int32_t>("stepsGlobalValue", stepStart, stepCount);

    const std::vector<std::string> stepsGlobalValueString =
        iStream.read<std::string>("stepsGlobalValueString", stepStart, stepCount);

    const std::vector<int32_t> ranksGlobalArray =
        iStream.read<int32_t>("ranksLocalValue", stepStart, stepCount);

    const std::vector<std::string> ranksGlobalArrayString =
        iStream.read<std::string>("ranksLocalValueString", stepStart, stepCount);

    // loop blocks
    for (size_t b = 0; b < static_cast<size_t>(mpiSize); ++b)
    {
        // loop selections
        for (size_t i = 0; i < Nx; ++i)
        {
            const adios2::Dims start = {i};
            const adios2::Dims count = {Nx - i};

            const auto I8 = iStream.read<int8_t>("i8", start, count, stepStart, stepCount, b);

            const auto I16 = iStream.read<int16_t>("i16", start, count, stepStart, stepCount, b);
            const auto I32 = iStream.read<int32_t>("i32", start, count, stepStart, stepCount, b);
            const auto I64 = iStream.read<int64_t>("i64", start, count, stepStart, stepCount, b);
            const auto U8 = iStream.read<uint8_t>("u8", start, count, stepStart, stepCount, b);
            const auto U16 = iStream.read<uint16_t>("u16", start, count, stepStart, stepCount, b);
            const auto U32 = iStream.read<uint32_t>("u32", start, count, stepStart, stepCount, b);
            const auto U64 = iStream.read<uint64_t>("u64", start, count, stepStart, stepCount, b);
            const auto R32 = iStream.read<float>("r32", start, count, stepStart, stepCount, b);
            const auto R64 = iStream.read<double>("r64", start, count, stepStart, stepCount, b);

            const auto CR32 =
                iStream.read<std::complex<float>>("cr32", start, count, stepStart, stepCount, b);
            const auto CR64 =
                iStream.read<std::complex<double>>("cr64", start, count, stepStart, stepCount, b);

            for (size_t s = 0; s < NSteps; ++s)
            {
                SmallTestData currentTestData = generateNewSmallTestData(
                    m_TestData, static_cast<int>(s), static_cast<int>(b), mpiSize);

                for (size_t j = i; j < Nx; ++j)
                {
                    std::stringstream ss;
                    ss << "t=" << s << " j=" << j << " rank=" << mpiRank;
                    std::string msg = ss.str();

                    ASSERT_EQ(I8[s * (Nx - i) + j - i], currentTestData.I8[j]) << msg;
                    ASSERT_EQ(I16[s * (Nx - i) + j - i], currentTestData.I16[j]) << msg;
                    EXPECT_EQ(I32[s * (Nx - i) + j - i], currentTestData.I32[j]) << msg;
                    EXPECT_EQ(I64[s * (Nx - i) + j - i], currentTestData.I64[j]) << msg;

                    EXPECT_EQ(U8[s * (Nx - i) + j - i], currentTestData.U8[j]) << msg;
                    EXPECT_EQ(U16[s * (Nx - i) + j - i], currentTestData.U16[j]) << msg;
                    EXPECT_EQ(U32[s * (Nx - i) + j - i], currentTestData.U32[j]) << msg;
                    EXPECT_EQ(U64[s * (Nx - i) + j - i], currentTestData.U64[j]) << msg;
                    EXPECT_EQ(R32[s * (Nx - i) + j - i], currentTestData.R32[j]) << msg;
                    EXPECT_EQ(R64[s * (Nx - i) + j - i], currentTestData.R64[j]) << msg;
                    EXPECT_EQ(CR32[s * (Nx - i) + j - i], currentTestData.CR32[j]) << msg;
                    EXPECT_EQ(CR64[s * (Nx - i) + j - i], currentTestData.CR64[j]) << msg;
                }
            }

        } // selection
    }     // blockID
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
