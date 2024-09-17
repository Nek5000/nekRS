/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <fstream>
#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

// std::string engineName; // comes from command line

struct QueryTestData
{
    std::vector<int> m_IntData;
    std::vector<double> m_DoubleData;
};

void WriteXmlQuery1D(const std::string &queryFile, const std::string &ioName,
                     const std::string &varName)
{
    std::ofstream file(queryFile.c_str());
    file << "<adios-query>" << std::endl;
    file << " <io name=\"" << ioName << "\">" << std::endl;
    file << "   <var name=\"" << varName << "\">" << std::endl;
    file << "      <boundingbox  start=\"5\" count=\"80\"/>" << std::endl;
    file << "       <op value=\"OR\">" << std::endl;
    file << "         <range  compare=\"GT\" value=\"100.6\"/>" << std::endl;
    file << "         <range  compare=\"LT\" value=\"-0.17\"/>" << std::endl;
    file << "         <op value=\"AND\">" << std::endl;
    file << "            <range  compare=\"LT\" value=\"11.9\"/>" << std::endl;
    file << "            <range  compare=\"GT\" value=\"2.8\"/>" << std::endl;
    file << "         </op>" << std::endl;
    file << "       </op>" << std::endl;
    file << "   </var>" << std::endl;
    file << " </io>" << std::endl;
    file << "</adios-query>" << std::endl;
    file.close();
}

void LoadTestData(QueryTestData &input, int step, int rank, int dataSize)
{
    input.m_IntData.clear();
    input.m_IntData.reserve(dataSize);
    input.m_DoubleData.clear();
    input.m_DoubleData.reserve(dataSize);

    for (auto i = 0; i < dataSize; i++)
    {
        input.m_IntData.push_back(step * 1000 + rank * 100 + i);
        input.m_DoubleData.push_back(step * 100 + rank * 10 + i * 0.01);
    }
}

class BPQueryTest : public ::testing::Test
{
public:
    BPQueryTest() = default;

    void WriteFile(const std::string &fname, adios2::ADIOS &adios, const std::string &engineName);
    void QueryDoubleVar(const std::string &fname, adios2::ADIOS &adios,
                        const std::string &engineName);
    void QueryIntVar(const std::string &fname, adios2::ADIOS &adios, const std::string &engineName);

    QueryTestData m_TestData;

    // Number of rows
    const size_t Nx = 100;
    // Number of steps
    const size_t NSteps = 3;

    int mpiRank = 0, mpiSize = 1;
};

void BPQueryTest::QueryIntVar(const std::string &fname, adios2::ADIOS &adios,
                              const std::string &engineName)
{
#if ADIOS2_USE_MPI
    std::string ioName = "IOQueryTestInt_MPI" + engineName;
#else
    std::string ioName = "IOQueryTestInt" + engineName;
#endif
    adios2::IO io = adios.DeclareIO(ioName.c_str());

    if (!engineName.empty())
    {
        io.SetEngine(engineName);
    }

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    EXPECT_EQ(bpReader.Steps(), NSteps);

    std::string queryFile = "./" + ioName + "test.xml"; //"./test.xml";
    std::cout << ioName << std::endl;
    WriteXmlQuery1D(queryFile, ioName, "intV");

    std::vector<size_t> rr;
    if (engineName.compare("BP4") == 0)
        rr = {2, 1, 1};
    else
        rr = {1, 1, 1};

    while (bpReader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::QueryWorker w = adios2::QueryWorker(queryFile, bpReader);
        std::vector<adios2::Box<adios2::Dims>> touched_blocks;
        w.GetResultCoverage(touched_blocks);
        ASSERT_EQ(touched_blocks.size(), rr[bpReader.CurrentStep()]);
        bpReader.EndStep();
    }
    bpReader.Close();
}

void BPQueryTest::QueryDoubleVar(const std::string &fname, adios2::ADIOS &adios,
                                 const std::string &engineName)
{
    std::string ioName = "IOQueryTestDouble" + engineName;
    adios2::IO io = adios.DeclareIO(ioName.c_str());

    if (!engineName.empty())
    {
        io.SetEngine(engineName);
    }

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

    EXPECT_EQ(bpReader.Steps(), NSteps);

    // std::string queryFile = "./.test.xml";
    std::string queryFile = "./" + ioName + "test.xml";
    WriteXmlQuery1D(queryFile, ioName, "doubleV");

    std::vector<size_t> rr; //= {0,9,9};
    if (engineName.compare("BP4") == 0)
        rr = {0, 3, 1};
    else
        rr = {0, 1, 1};
    while (bpReader.BeginStep() == adios2::StepStatus::OK)
    {
        adios2::QueryWorker w = adios2::QueryWorker(queryFile, bpReader);
        std::vector<adios2::Box<adios2::Dims>> touched_blocks;
        w.GetResultCoverage(touched_blocks);
        ASSERT_EQ(touched_blocks.size(), rr[bpReader.CurrentStep()]);
        bpReader.EndStep();
    }
    bpReader.Close();
}

void BPQueryTest::WriteFile(const std::string &fname, adios2::ADIOS &adios,
                            const std::string &engineName)
{

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

    // Write test data using BP

    {
        adios2::IO io = adios.DeclareIO("TestQueryIOWriter");

        // Declare 1D variables (NumOfProcesses * Nx)
        // The local process' part (start, count) can be defined now or later
        // before Write().
        {
            const adios2::Dims shape{static_cast<size_t>(Nx * mpiSize)};
            const adios2::Dims start{static_cast<size_t>(Nx * mpiRank)};
            const adios2::Dims count{Nx};

            auto var_i32 = io.DefineVariable<int32_t>("intV", shape, start, count);
            auto var_r64 = io.DefineVariable<double>("doubleV", shape, start, count);
            EXPECT_TRUE(var_i32);
            EXPECT_TRUE(var_r64);
        }

        if (!engineName.empty())
        {
            io.SetEngine(engineName);
        }
        else
        {
            // Create the BP Engine
            io.SetEngine("BPFile");
        }

        if (engineName.compare("BP4") == 0)
        {
            io.SetParameters("statslevel=1");
            io.SetParameters("statsblocksize=10");
        }
        io.AddTransport("file");

        // QUESTION: It seems that BPFilterWriter cannot overwrite existing
        // files
        // Ex. if you tune Nx and NSteps, the test would fail. But if you clear
        // the cache in
        // ${adios2Build}/testing/adios2/engine/bp/ADIOS2BPWriteADIOS1Read1D800.bp.dir,
        // then it works
        adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

        for (size_t step = 0; step < NSteps; ++step)
        {
            // Generate test data for each process uniquely
            LoadTestData(m_TestData, static_cast<int>(step), mpiRank, static_cast<int>(Nx));

            auto var_i32 = io.InquireVariable<int32_t>("intV");
            auto var_r64 = io.InquireVariable<double>("doubleV");

            adios2::Box<adios2::Dims> sel({mpiRank * Nx}, {Nx});

            var_i32.SetSelection(sel);
            var_r64.SetSelection(sel);

            bpWriter.BeginStep();
            bpWriter.Put(var_i32, m_TestData.m_IntData.data());
            bpWriter.Put(var_r64, m_TestData.m_DoubleData.data());

            bpWriter.EndStep();
        }

        // Close the file
        bpWriter.Close();
    }
}
//******************************************************************************
// 1D  test data
//******************************************************************************

TEST_F(BPQueryTest, BP5)
{
    std::string engineName = "BP5";
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
    const std::string fname(engineName + "Query1D_MPI.bp");
#else
    adios2::ADIOS adios;
    const std::string fname(engineName + "Query1D.bp");
#endif

    WriteFile(fname, adios, engineName);

    if (mpiSize == 1)
    {
        QueryDoubleVar(fname, adios, engineName);
        QueryIntVar(fname, adios, engineName);
    }
}

//******************************************************************************
// 2D  test data
//******************************************************************************

TEST_F(BPQueryTest, BP4)
{
    std::string engineName = "BP4";
    // Each process would write a 1x8 array and all processes would
    // form a mpiSize * Nx 1D array
    const std::string fname(engineName + "Query1D.bp");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(MPI_COMM_WORLD);
#else
    adios2::ADIOS adios;
#endif

    WriteFile(fname, adios, engineName);

    if (mpiSize == 1)
    {
        QueryDoubleVar(fname, adios, engineName);
        QueryIntVar(fname, adios, engineName);
    }
}

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    if (argc > 1)
    {
        // engineName = std::string(argv[1]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
