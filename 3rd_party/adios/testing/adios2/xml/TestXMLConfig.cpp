#include <cstdint>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#define str_helper(X) #X
#define str(X) str_helper(X)

class XMLConfigTest : public ::testing::Test
{
public:
    XMLConfigTest() : configDir(str(XML_CONFIG_DIR)) {}

    // protected:
    // virtual void SetUp() { }

    // virtual void TearDown() { }
    std::string configDir;
};

TEST_F(XMLConfigTest, TwoIOs)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "config1.xml");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(configFile, MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(configFile);
#endif

    // must be declared at least once
    EXPECT_THROW(adios2::IO io = adios.AtIO("Test IO 1"); (void)io, std::invalid_argument);

    EXPECT_NO_THROW({
        adios2::IO io = adios.DeclareIO("Test IO 1");
        const adios2::Params params = io.Parameters();
        ASSERT_EQ(params.size(), 5);
        EXPECT_THROW((void)params.at("DoesNotExist"), std::out_of_range);
        EXPECT_EQ(params.at("Threads"), "1");
        EXPECT_EQ(params.at("ProfileUnits"), "Microseconds");
        EXPECT_EQ(params.at("MaxBufferSize"), "20Mb");
        EXPECT_EQ(params.at("InitialBufferSize"), "1Mb");
        EXPECT_EQ(params.at("BufferGrowthFactor"), "2");
        adios2::Engine engine = io.Open("Test BP Writer 1", adios2::Mode::Write);
        engine.Close();
    });
    EXPECT_NO_THROW(adios2::IO io = adios.AtIO("Test IO 1"); (void)io);

    EXPECT_THROW(adios2::IO io = adios.AtIO("Test IO 2"); (void)io, std::invalid_argument);
    EXPECT_NO_THROW({
        adios2::IO io = adios.DeclareIO("Test IO 2");
        const adios2::Params params = io.Parameters();
        ASSERT_EQ(params.size(), 0);
    });
    EXPECT_NO_THROW(adios.AtIO("Test IO 2"));

    // double declaring
    EXPECT_THROW(adios.DeclareIO("Test IO 1"), std::invalid_argument);
    EXPECT_THROW(adios.DeclareIO("Test IO 2"), std::invalid_argument);
}

TEST_F(XMLConfigTest, TwoEnginesException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "config2.xml");

#if ADIOS2_USE_MPI
    EXPECT_THROW(adios2::ADIOS adios(configFile, MPI_COMM_WORLD), std::invalid_argument);
#else
    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
#endif
}

TEST_F(XMLConfigTest, OpTypeException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpTypeException.xml");

#if ADIOS2_USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        EXPECT_THROW(adios2::ADIOS adios(configFile, MPI_COMM_SELF), std::invalid_argument);
    }
#else
    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
#endif
}

TEST_F(XMLConfigTest, OpNullException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpNullException.xml");

#if ADIOS2_USE_MPI
    EXPECT_THROW(adios2::ADIOS adios(configFile, MPI_COMM_WORLD), std::invalid_argument);
#else
    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
#endif
}

TEST_F(XMLConfigTest, OpNoneException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpNoneException.xml");

#if ADIOS2_USE_MPI
    EXPECT_THROW(adios2::ADIOS adios(configFile, MPI_COMM_WORLD), std::invalid_argument);
#else
    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
#endif
}

TEST_F(XMLConfigTest, RemoveIO)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configRemoveIO.xml");

#if ADIOS2_USE_MPI
    adios2::ADIOS adios(configFile, MPI_COMM_WORLD);
#else
    adios2::ADIOS adios(configFile);
#endif

    adios2::IO io;
    adios2::Engine engine;

    std::string io_name_ = "checkpoint";
    for (int c = 0; c < 3; c++)
    {
        io = adios.DeclareIO(io_name_);
        std::string filename = "test.bp";
        engine = io.Open(filename, adios2::Mode::Write);
        EXPECT_TRUE(io.EngineType() == "BP4");
        engine.Close();
        adios.RemoveIO(io_name_);
    }
}

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
