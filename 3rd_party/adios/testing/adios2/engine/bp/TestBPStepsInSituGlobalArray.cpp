/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <array>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

// Number of elements per process
const std::size_t Nx = 10;
using DataArray = std::array<int32_t, Nx>;

class BPStepsInSituGlobalArray : public ::testing::Test
{
protected:
    BPStepsInSituGlobalArray() = default;

    const DataArray I32 = {{512, 513, -510, 515, -508, 517, -506, 519, -504, 521}};

    DataArray GenerateData(const int step, const int rank, const int size)
    {
        DataArray d;
        int j = rank + 1 + step * size;
        for (size_t i = 0; i < d.size(); ++i)
        {
            d[i] = I32[i] + j;
        }
        return d;
    }

    std::string ArrayToString(const int32_t *data, const size_t nelems)
    {
        std::stringstream ss;
        ss << "[";
        for (size_t i = 0; i < nelems; ++i)
        {
            ss << data[i];
            if (i < nelems - 1)
            {
                ss << " ";
            }
        }
        ss << "]";
        return ss.str();
    }
};

enum class ReadMode
{
    ReadGlobal,
    ReadBlocks
};

std::string ReadModeToString(const ReadMode r)
{
    switch (r)
    {
    case ReadMode::ReadGlobal:
        return "ReadGlobal";
    case ReadMode::ReadBlocks:
        return "ReadBlocks";
    }
    return "unknown";
}

enum class Act
{
    Write,
    Read
};

const std::vector<std::vector<Act>> Schedules = {
    {Act::Write, Act::Write, Act::Write, Act::Read, Act::Read, Act::Read},
    {Act::Write, Act::Write, Act::Read, Act::Write, Act::Read, Act::Read},
    {Act::Write, Act::Write, Act::Read, Act::Read, Act::Write, Act::Read},
    {Act::Write, Act::Read, Act::Write, Act::Write, Act::Read, Act::Read},
    {Act::Write, Act::Read, Act::Write, Act::Read, Act::Write, Act::Read}};

std::string ScheduleToString(const std::vector<Act> &schedule)
{
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < schedule.size(); ++i)
    {
        if (schedule[i] == Act::Write)
        {
            ss << "Write";
        }
        else if (schedule[i] == Act::Read)
        {
            ss << "Read";
        }
        if (i < schedule.size() - 1)
        {
            ss << " ";
        }
    }
    ss << "]";
    return ss.str();
}

class BPStepsInSituGlobalArrayReaders
: public BPStepsInSituGlobalArray,
  public ::testing::WithParamInterface<std::tuple<size_t, ReadMode>>
{
protected:
    const std::vector<Act> &GetSchedule() { return Schedules[std::get<0>(GetParam())]; };
    ReadMode GetReadMode() { return std::get<1>(GetParam()); };
    size_t GetScheduleID() { return std::get<0>(GetParam()); };
};

// Basic case: Variable written every step
TEST_P(BPStepsInSituGlobalArrayReaders, EveryStep)
{
    const std::vector<Act> &schedule = GetSchedule();
    const ReadMode readMode = GetReadMode();
    using namespace std;
    std::string BaseName = engineName;
    const std::string fname_prefix = "BPStepsInSituGlobalArray.EveryStep." +
                                     std::to_string(GetScheduleID()) + "." +
                                     ReadModeToString(readMode) + "." + engineName;
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    std::vector<DataArray> m_TestData;
    adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
    adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
    adios2::Dims count{static_cast<unsigned int>(Nx)};

    std::string fname;
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    fname = fname_prefix + BaseName + ".MPI.bp";
#else
    adios2::ADIOS adios;
    fname = fname_prefix + BaseName + ".Serial.bp";
#endif

    if (!mpiRank)
    {
        std::cout << "Test with Schedule " << GetScheduleID() << " " << ScheduleToString(schedule)
                  << " Read Mode " << ReadModeToString(readMode) << std::endl;
    }

    // Start writer
    adios2::IO iow = adios.DeclareIO("Write");
    if (!engineName.empty())
    {
        iow.SetEngine(engineName);
    }
    adios2::Engine writer = iow.Open(fname, adios2::Mode::Write);
    EXPECT_TRUE(writer);
    auto var_i32 = iow.DefineVariable<int32_t>("i32", shape, start, count);
    EXPECT_TRUE(var_i32);

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Start reader
    adios2::IO ior = adios.DeclareIO("Read");
    if (!engineName.empty())
    {
        ior.SetEngine(engineName);
    }
    ior.SetParameter("OpenTimeoutSecs", "10.0");
    adios2::Engine reader = ior.Open(fname, adios2::Mode::Read);
    EXPECT_TRUE(reader);

    int stepsWritten = 0;
    int stepsRead = 0;

    for (const auto act : schedule)
    {
        if (act == Act::Write)
        {
            // Write test data using ADIOS2
            if (!mpiRank)
            {
                std::cout << "Write step " << stepsWritten << std::endl;
            }

            // Generate test data for each process uniquely
            m_TestData.push_back(GenerateData(stepsWritten, mpiRank, mpiSize));
            std::cout << "Rank " << mpiRank << " write step " << stepsWritten << ": "
                      << ArrayToString(m_TestData[stepsWritten].data(), Nx) << std::endl;
            writer.BeginStep();
            writer.Put(var_i32, m_TestData[stepsWritten].data());
            writer.EndStep();
            ++stepsWritten;
        }
        else if (act == Act::Read)
        {
            if (readMode == ReadMode::ReadGlobal)
            {
                /// Read back data with global selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Global selection" << std::endl;
                }

                reader.BeginStep();
                auto var_i32 = ior.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                // EXPECT_EQ(var_i32.Steps(), 1);
                EXPECT_EQ(var_i32.StepsStart(), 0);
                size_t start = static_cast<size_t>(mpiRank) * Nx;
                var_i32.SetSelection({{start}, {Nx}});
                DataArray d;
                reader.Get(var_i32, d.data(), adios2::Mode::Sync);
                std::cout << "Rank " << mpiRank << " read step " << stepsRead << ": "
                          << ArrayToString(d.data(), Nx) << std::endl;
                for (size_t i = 0; i < Nx; ++i)
                {
                    EXPECT_EQ(d[i], m_TestData[stepsRead][i]);
                }
                reader.EndStep();
            }
            else if (readMode == ReadMode::ReadBlocks)
            {
                /// Read back data with block selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Block selection" << std::endl;
                }

                reader.BeginStep();
                auto var_i32 = ior.InquireVariable<int32_t>("i32");
                EXPECT_TRUE(var_i32);
                // EXPECT_EQ(var_i32.Steps(), 1);
                EXPECT_EQ(var_i32.StepsStart(), 0);
                size_t blockID = static_cast<size_t>(mpiRank);
                var_i32.SetBlockSelection(blockID);
                DataArray d;
                reader.Get(var_i32, d.data(), adios2::Mode::Sync);
                std::cout << "Rank " << mpiRank << " read step " << stepsRead << " block "
                          << blockID << ": " << ArrayToString(d.data(), Nx) << std::endl;
                auto start = var_i32.Start();
                auto count = var_i32.Count();
                EXPECT_EQ(start[0], mpiRank * Nx);
                EXPECT_EQ(count[0], 1 * Nx);
                for (size_t i = 0; i < Nx; ++i)
                {
                    EXPECT_EQ(d[i], m_TestData[stepsRead][i]);
                }
                reader.EndStep();
            }
            ++stepsRead;
        }
#if ADIOS2_USE_MPI
        std::flush(std::cout);
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    writer.Close();
    reader.Close();
}

// A new variable is created and written every step
TEST_P(BPStepsInSituGlobalArrayReaders, NewVarPerStep)
{
    const std::vector<Act> &schedule = GetSchedule();
    const ReadMode readMode = GetReadMode();
    using namespace std;
    std::string BaseName = engineName;
    const std::string fname_prefix = "BPStepsInSituGlobalArray.NewVarPerStep." +
                                     std::to_string(GetScheduleID()) + "." +
                                     ReadModeToString(readMode) + "." + engineName;
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    std::vector<DataArray> m_TestData;
    adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
    adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
    adios2::Dims count{static_cast<unsigned int>(Nx)};

    std::string fname;
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    fname = fname_prefix + BaseName + ".MPI.bp";
#else
    adios2::ADIOS adios;
    fname = fname_prefix + BaseName + ".Serial.bp";
#endif

    auto lf_VarName = [](int step) -> std::string { return "i32_" + std::to_string(step); };

    std::cout << "Test with Schedule " << GetScheduleID() << " " << ScheduleToString(schedule)
              << " Read Mode " << ReadModeToString(readMode) << std::endl;

    // Start writer
    adios2::IO iow = adios.DeclareIO("Write");
    if (!engineName.empty())
    {
        iow.SetEngine(engineName);
    }
    adios2::Engine writer = iow.Open(fname, adios2::Mode::Write);
    EXPECT_TRUE(writer);
    auto var_i32 = iow.DefineVariable<int32_t>("i32", shape, start, count);
    EXPECT_TRUE(var_i32);

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Start reader
    adios2::IO ior = adios.DeclareIO("Read");
    if (!engineName.empty())
    {
        ior.SetEngine(engineName);
    }
    ior.SetParameter("OpenTimeoutSecs", "10.0");
    adios2::Engine reader = ior.Open(fname, adios2::Mode::Read);
    EXPECT_TRUE(reader);

    int stepsWritten = 0;
    int stepsRead = 0;

    for (const auto act : schedule)
    {
        if (act == Act::Write)
        {
            if (!mpiRank)
            {
                std::cout << "Write step " << stepsWritten << std::endl;
            }

            const std::string varName = lf_VarName(stepsWritten);
            auto var = iow.DefineVariable<int32_t>(varName, shape, start, count);
            // Generate test data for each process uniquely
            m_TestData.push_back(GenerateData(stepsWritten, mpiRank, mpiSize));
            std::cout << "Rank " << mpiRank << " write step " << stepsWritten << " var " << varName
                      << ": " << ArrayToString(m_TestData[stepsWritten].data(), Nx) << std::endl;
            writer.BeginStep();
            writer.Put(var, m_TestData[stepsWritten].data());
            writer.EndStep();
            ++stepsWritten;
        }
        else if (act == Act::Read)
        {
            if (readMode == ReadMode::ReadGlobal)
            {
                /// Read back data with global selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Global selection" << std::endl;
                }

                reader.BeginStep();
                const std::string varName = lf_VarName(stepsRead);
                auto var = ior.InquireVariable<int32_t>(varName);
                EXPECT_TRUE(var);
                EXPECT_EQ(var.Steps(), 1);
                EXPECT_EQ(var.StepsStart(), 0);
                size_t start = static_cast<size_t>(mpiRank) * Nx;
                var.SetSelection({{start}, {Nx}});
                DataArray d;
                reader.Get(var, d.data(), adios2::Mode::Sync);
                std::cout << "Rank " << mpiRank << " read var " << varName << ": "
                          << ArrayToString(d.data(), Nx) << std::endl;

                for (size_t i = 0; i < Nx; ++i)
                {
                    EXPECT_EQ(d[i], m_TestData[stepsRead][i]);
                }
                reader.EndStep();
            }
            else if (readMode == ReadMode::ReadBlocks)
            {
                /// Read back data with block selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Block selection" << std::endl;
                }

                reader.BeginStep();
                const std::string varName = lf_VarName(stepsRead);
                auto var = ior.InquireVariable<int32_t>(varName);
                EXPECT_TRUE(var);
                EXPECT_EQ(var.Steps(), 1);
                EXPECT_EQ(var.StepsStart(), 0);
                size_t blockID = static_cast<size_t>(mpiRank);
                var.SetBlockSelection(blockID);
                DataArray d;
                reader.Get(var, d.data(), adios2::Mode::Sync);
                std::cout << "Rank " << mpiRank << " read step " << stepsRead << " block "
                          << blockID << ": " << ArrayToString(d.data(), Nx) << std::endl;
                auto start = var.Start();
                auto count = var.Count();
                EXPECT_EQ(start[0], mpiRank * Nx);
                EXPECT_EQ(count[0], 1 * Nx);
                for (size_t i = 0; i < Nx; ++i)
                {
                    EXPECT_EQ(d[i], m_TestData[stepsRead][i]);
                }
                reader.EndStep();
            }
            ++stepsRead;
        }
#if ADIOS2_USE_MPI
        std::flush(std::cout);
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    writer.Close();
    reader.Close();
}

INSTANTIATE_TEST_SUITE_P(BPStepsInSituGlobalArray, BPStepsInSituGlobalArrayReaders,
                         ::testing::Values(std::make_tuple(0, ReadMode::ReadGlobal),
                                           std::make_tuple(1, ReadMode::ReadGlobal),
                                           std::make_tuple(2, ReadMode::ReadGlobal),
                                           std::make_tuple(3, ReadMode::ReadGlobal),
                                           std::make_tuple(4, ReadMode::ReadGlobal),
                                           std::make_tuple(0, ReadMode::ReadBlocks),
                                           std::make_tuple(1, ReadMode::ReadBlocks),
                                           std::make_tuple(2, ReadMode::ReadBlocks),
                                           std::make_tuple(3, ReadMode::ReadBlocks),
                                           std::make_tuple(4, ReadMode::ReadBlocks)));

class BPStepsInSituGlobalArrayParameters
: public BPStepsInSituGlobalArray,
  public ::testing::WithParamInterface<std::tuple<size_t, size_t, ReadMode>>
{
protected:
    const std::vector<Act> &GetSchedule() { return Schedules[std::get<0>(GetParam())]; };
    size_t GetOddity() { return std::get<1>(GetParam()); };
    ReadMode GetReadMode() { return std::get<2>(GetParam()); };
    size_t GetScheduleID() { return std::get<0>(GetParam()); };
};

// A variable written every other step either from step 0 (EVEN) or from
// step 1 (ODD)
TEST_P(BPStepsInSituGlobalArrayParameters, EveryOtherStep)
{
    const std::vector<Act> &schedule = GetSchedule();
    const std::size_t Oddity = GetOddity();
    const ReadMode readMode = GetReadMode();
    using namespace std;
    std::string BaseName = engineName;
    const std::string fname_prefix =
        "BPStepsInSituGlobalArray.EveryOtherStep.Schedule" + std::to_string(GetScheduleID()) +
        ".Oddity" + std::to_string(Oddity) + "." + ReadModeToString(readMode) + "." + engineName;
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    std::vector<DataArray> m_TestData;
    adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
    adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
    adios2::Dims count{static_cast<unsigned int>(Nx)};

    std::string fname;
#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    fname = fname_prefix + BaseName + ".MPI.bp";
#else
    adios2::ADIOS adios;
    fname = fname_prefix + BaseName + ".Serial.bp";

#endif

    if (!mpiRank)
    {
        std::cout << "Test with Schedule " << GetScheduleID() << " " << ScheduleToString(schedule)
                  << " Oddity " << Oddity << " Read Mode " << ReadModeToString(readMode)
                  << std::endl;
    }

    // Start writer
    adios2::IO iow = adios.DeclareIO("Write");
    if (!engineName.empty())
    {
        iow.SetEngine(engineName);
    }
    adios2::Engine writer = iow.Open(fname, adios2::Mode::Write);
    EXPECT_TRUE(writer);
    auto var_i32 = iow.DefineVariable<int32_t>("i32", shape, start, count);
    EXPECT_TRUE(var_i32);
    auto var_step = iow.DefineVariable<int>("step");
    EXPECT_TRUE(var_step);

#if ADIOS2_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Start reader
    adios2::IO ior = adios.DeclareIO("Read");
    if (!engineName.empty())
    {
        ior.SetEngine(engineName);
    }
    ior.SetParameter("OpenTimeoutSecs", "10.0");
    adios2::Engine reader = ior.Open(fname, adios2::Mode::Read);
    EXPECT_TRUE(reader);

    int stepsWritten = 0;
    int stepsRead = 0;
    int varStepsWritten = 0;
    int varStepsRead = 0;
    for (const auto act : schedule)
    {
        if (act == Act::Write)
        {
            if (!mpiRank)
            {
                std::cout << "Write step " << stepsWritten << std::endl;
            }

            // Generate test data for each process uniquely
            writer.BeginStep();
            writer.Put(var_step, stepsWritten);
            if (stepsWritten % 2 == static_cast<int>(Oddity))
            {
                m_TestData.push_back(GenerateData(stepsWritten, mpiRank, mpiSize));
                std::cout << "Rank " << mpiRank << " at step " << stepsWritten << " write var-step "
                          << varStepsWritten << ": "
                          << ArrayToString(m_TestData[varStepsWritten].data(), Nx) << std::endl;
                writer.Put(var_i32, m_TestData[varStepsWritten].data());
                ++varStepsWritten;
            }
            writer.EndStep();
            ++stepsWritten;
        }
        else if (act == Act::Read)
        {

            if (readMode == ReadMode::ReadGlobal)
            {
                /// Read back data with global selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Global selection" << std::endl;
                }

                reader.BeginStep();
                if (stepsRead % 2 == static_cast<int>(Oddity))
                {
                    auto var_i32 = ior.InquireVariable<int32_t>("i32");
                    EXPECT_TRUE(var_i32);
                    // EXPECT_EQ(var_i32.Steps(), 1);
                    EXPECT_EQ(var_i32.StepsStart(), 0);
                    size_t start = static_cast<size_t>(mpiRank) * Nx;
                    var_i32.SetSelection({{start}, {Nx}});
                    DataArray d;
                    reader.Get(var_i32, d.data(), adios2::Mode::Sync);
                    std::cout << "Rank " << mpiRank << " read at step " << stepsRead << " var-step "
                              << varStepsRead << ": " << ArrayToString(d.data(), Nx) << std::endl;

                    for (size_t i = 0; i < Nx; ++i)
                    {
                        EXPECT_EQ(d[i], m_TestData[varStepsRead][i]);
                    }
                    ++varStepsRead;
                }
                reader.EndStep();
            }
            else if (readMode == ReadMode::ReadBlocks)
            {
                /// Read back data with block selection
                if (!mpiRank)
                {
                    std::cout << "Read step " << stepsRead << " with Block selection" << std::endl;
                }
                reader.BeginStep();
                if (stepsRead % 2 == static_cast<int>(Oddity))
                {
                    auto var_i32 = ior.InquireVariable<int32_t>("i32");
                    EXPECT_TRUE(var_i32);
                    // EXPECT_EQ(var_i32.Steps(), 1);
                    EXPECT_EQ(var_i32.StepsStart(), 0);
                    size_t blockID = static_cast<size_t>(mpiRank);
                    var_i32.SetBlockSelection(blockID);
                    DataArray d;
                    reader.Get(var_i32, d.data(), adios2::Mode::Sync);
                    std::cout << "Rank " << mpiRank << " read at step " << stepsRead << " var-step "
                              << varStepsRead << " block " << blockID << ": "
                              << ArrayToString(d.data(), Nx) << std::endl;
                    auto start = var_i32.Start();
                    auto count = var_i32.Count();
                    EXPECT_EQ(start[0], mpiRank * Nx);
                    EXPECT_EQ(count[0], 1 * Nx);
                    for (size_t i = 0; i < Nx; ++i)
                    {
                        EXPECT_EQ(d[i], m_TestData[varStepsRead][i]);
                    }
                    ++varStepsRead;
                }
                reader.EndStep();
            }
            ++stepsRead;
        }

#if ADIOS2_USE_MPI
        std::flush(std::cout);
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    writer.Close();
    reader.Close();
}

INSTANTIATE_TEST_SUITE_P(
    BPStepsInSituGlobalArray, BPStepsInSituGlobalArrayParameters,
    ::testing::Values(
        std::make_tuple(0, 0, ReadMode::ReadGlobal), std::make_tuple(0, 0, ReadMode::ReadBlocks),
        std::make_tuple(0, 1, ReadMode::ReadGlobal), std::make_tuple(0, 1, ReadMode::ReadBlocks),
        std::make_tuple(1, 0, ReadMode::ReadGlobal), std::make_tuple(1, 0, ReadMode::ReadBlocks),
        std::make_tuple(1, 1, ReadMode::ReadGlobal), std::make_tuple(1, 1, ReadMode::ReadBlocks),
        std::make_tuple(2, 0, ReadMode::ReadGlobal), std::make_tuple(2, 0, ReadMode::ReadBlocks),
        std::make_tuple(2, 1, ReadMode::ReadGlobal), std::make_tuple(2, 1, ReadMode::ReadBlocks),
        std::make_tuple(3, 0, ReadMode::ReadGlobal), std::make_tuple(3, 0, ReadMode::ReadBlocks),
        std::make_tuple(3, 1, ReadMode::ReadGlobal), std::make_tuple(3, 1, ReadMode::ReadBlocks),
        std::make_tuple(4, 0, ReadMode::ReadGlobal), std::make_tuple(4, 0, ReadMode::ReadBlocks),
        std::make_tuple(4, 1, ReadMode::ReadGlobal), std::make_tuple(4, 1, ReadMode::ReadBlocks)));
//******************************************************************************
// main
//******************************************************************************

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
