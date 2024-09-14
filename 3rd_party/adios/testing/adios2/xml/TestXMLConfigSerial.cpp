#include <cstdint>

#include <iostream>
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

#define str_helper(X) #X
#define str(X) str_helper(X)

class XMLConfigSerialTest : public ::testing::Test
{
public:
    XMLConfigSerialTest() : configDir(str(XML_CONFIG_DIR)) {}

    // protected:
    // virtual void SetUp() { }

    // virtual void TearDown() { }
    std::string configDir;
};

TEST_F(XMLConfigSerialTest, TwoIOs)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "config1.xml");

    adios2::ADIOS adios(configFile);

    // must be declared at least once
    EXPECT_THROW(adios2::IO io = adios.AtIO("Test IO 1"), std::invalid_argument);

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
    EXPECT_NO_THROW(adios2::IO io = adios.AtIO("Test IO 1"));

    EXPECT_THROW(adios2::IO io = adios.AtIO("Test IO 2"), std::invalid_argument);
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

TEST_F(XMLConfigSerialTest, TwoEnginesException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "config2.xml");

    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
}

TEST_F(XMLConfigSerialTest, OpTypeException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpTypeException.xml");

    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
}

TEST_F(XMLConfigSerialTest, OpNullException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpNullException.xml");

    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
}

TEST_F(XMLConfigSerialTest, OpNoneException)
{
    const std::string configFile(configDir + std::string(&adios2::PathSeparator, 1) +
                                 "configOpNoneException.xml");

    EXPECT_THROW(adios2::ADIOS adios(configFile), std::invalid_argument);
}

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}
