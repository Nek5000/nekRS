/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <array>
#include <stdexcept>
#include <tuple>

#include <adios2.h>

#include <gtest/gtest.h>

class BufferTest
: public ::testing::TestWithParam<std::tuple<std::string, std::string, std::string, std::string>>
{
};

TEST_P(BufferTest, WriteRead)
{
    const std::string &transportWriteLibrary = std::get<0>(GetParam());
    const std::string &transportWriteBuffer = std::get<1>(GetParam());
    const std::string &transportReadLibrary = std::get<2>(GetParam());
    const std::string &transportReadBuffer = std::get<3>(GetParam());

    const std::string fname("FileBufferTest_" + transportWriteLibrary + "_" + transportWriteBuffer +
                            "_" + transportReadLibrary + "_" + transportReadBuffer + ".bp");

    std::array<double, 100> dataOrig;
    for (size_t i = 0; i < 100; ++i)
    {
        dataOrig[i] = static_cast<double>(i);
    }

    adios2::ADIOS adios;
    {
        adios2::IO io = adios.DeclareIO("TestIO");

        io.SetEngine("BP4");
        const size_t transportID = io.AddTransport("file");
        io.SetTransportParameter(transportID, "Library", transportWriteLibrary);
        io.SetTransportParameter(transportID, "Buffer", transportWriteBuffer);

        auto var = io.DefineVariable<double>("var", {100}, {0}, {100});
        adios2::Engine writer = io.Open(fname, adios2::Mode::Write);

        writer.Put(var, dataOrig.data());

        // Close the file
        writer.Close();
    }

    std::array<double, 100> dataRead;
    {
        adios2::IO io = adios.DeclareIO("ReadIO");

        io.SetEngine("BP4");
        const size_t transportID = io.AddTransport("file");
        io.SetTransportParameter(transportID, "Library", transportReadLibrary);
        io.SetTransportParameter(transportID, "Buffer", transportReadBuffer);

        adios2::Engine reader = io.Open(fname, adios2::Mode::Read);
        auto var = io.InquireVariable<double>("var");
        ASSERT_EQ(var.Shape().size(), 1);
        ASSERT_EQ(var.Shape()[0], 100);

        reader.Get(var, dataRead.data());
        reader.Close();
    }

    for (size_t i = 0; i < 100; ++i)
    {
        ASSERT_EQ(dataOrig[i], dataRead[i]);
    }
}

#ifdef __unix__
INSTANTIATE_TEST_SUITE_P(TransportTests, BufferTest,
                         ::testing::Values(std::make_tuple("fstream", "true", "posix", "false"),
                                           std::make_tuple("fstream", "false", "posix", "false"),
                                           std::make_tuple("posix", "false", "stdio", "true"),
                                           std::make_tuple("posix", "false", "stdio", "false"),
                                           std::make_tuple("posix", "false", "fstream", "true"),
                                           std::make_tuple("posix", "false", "fstream", "false"),
                                           std::make_tuple("posix", "false", "posix", "false"),
                                           std::make_tuple("stdio", "true", "posix", "false"),
                                           std::make_tuple("stdio", "false", "posix", "false"),

                                           std::make_tuple("stdio", "true", "stdio", "true"),
                                           std::make_tuple("stdio", "true", "stdio", "false"),
                                           std::make_tuple("stdio", "false", "stdio", "true"),
                                           std::make_tuple("stdio", "false", "stdio", "false"),
                                           std::make_tuple("stdio", "true", "fstream", "true"),
                                           std::make_tuple("stdio", "true", "fstream", "false"),
                                           std::make_tuple("stdio", "false", "fstream", "true"),
                                           std::make_tuple("stdio", "false", "fstream", "false"),

                                           std::make_tuple("fstream", "true", "stdio", "true"),
                                           std::make_tuple("fstream", "true", "stdio", "false"),
                                           std::make_tuple("fstream", "false", "stdio", "true"),
                                           std::make_tuple("fstream", "false", "stdio", "false"),
                                           std::make_tuple("fstream", "true", "fstream", "true"),
                                           std::make_tuple("fstream", "true", "fstream", "false"),
                                           std::make_tuple("fstream", "false", "fstream", "true"),
                                           std::make_tuple("fstream", "false", "fstream",
                                                           "false")));
#else
INSTANTIATE_TEST_SUITE_P(TransportTests, BufferTest,
                         ::testing::Values(std::make_tuple("stdio", "true", "stdio", "true"),
                                           std::make_tuple("stdio", "true", "stdio", "false"),
                                           std::make_tuple("stdio", "false", "stdio", "true"),
                                           std::make_tuple("stdio", "false", "stdio", "false"),
                                           std::make_tuple("stdio", "true", "fstream", "true"),
                                           std::make_tuple("stdio", "true", "fstream", "false"),
                                           std::make_tuple("stdio", "false", "fstream", "true"),
                                           std::make_tuple("stdio", "false", "fstream", "false"),

                                           std::make_tuple("fstream", "true", "stdio", "true"),
                                           std::make_tuple("fstream", "true", "stdio", "false"),
                                           std::make_tuple("fstream", "false", "stdio", "true"),
                                           std::make_tuple("fstream", "false", "stdio", "false"),
                                           std::make_tuple("fstream", "true", "fstream", "true"),
                                           std::make_tuple("fstream", "true", "fstream", "false"),
                                           std::make_tuple("fstream", "false", "fstream", "true"),
                                           std::make_tuple("fstream", "false", "fstream",
                                                           "false")));
#endif

int main(int argc, char **argv)
{
    int result;
    ::testing::InitGoogleTest(&argc, argv);

    result = RUN_ALL_TESTS();

    return result;
}
