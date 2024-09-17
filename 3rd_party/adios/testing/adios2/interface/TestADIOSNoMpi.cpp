
#include <adios2.h>

#include <gtest/gtest.h>

/**
 * Basic test to see that writing and reading a simple BP4 file works without
 * initializing MPI
 */
TEST(ADIOSInterface, ADIOSNoMpi)
{
    adios2::ADIOS adios;
    {
        adios2::IO io = adios.DeclareIO("TestIOWrite");
        adios2::Engine engine = io.Open("test.bp", adios2::Mode::Write);
        engine.Close();
    }
    {
        adios2::IO io = adios.DeclareIO("TestIORead");
        adios2::Engine engine = io.Open("test.bp", adios2::Mode::Read);
        engine.Close();
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
