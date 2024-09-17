#include <adios2.h>
#include <fstream>
#include <gtest/gtest.h>

TEST(ADIOS2ReadNonBPFile, ShortBadFile)
{
    std::ofstream of{"foo.bp"};
    of << "Hello, world!\n";
    of.close();

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("Test");
    EXPECT_THROW(io.Open("foo.bp", adios2::Mode::Read), std::exception);
    std::remove("foo.bp");
}

TEST(ADIOS2ReadNonBPFile, LongBadFile)
{
    std::ofstream of{"foo.bp"};
    of << "hellohellohellohellohellohellohellohellohellohellohellohellohello\n";
    of.close();

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("Test");
    // auto bpReader = io.Open("foo.bp", adios2::Mode::Read);
    EXPECT_THROW(io.Open("foo.bp", adios2::Mode::Read), std::exception);
    std::remove("foo.bp");
}

TEST(ADIOS2ReadNonBPFile, LongBadFileWithValidEndianness)
{
    std::ofstream of("foo.bp", std::ios::out | std::ios::binary);
    char a = 0x00;
    for (size_t i = 0; i < 100; ++i)
    {
        of.write(static_cast<char *>(&a), 1);
    }
    of.close();

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("Test");
    EXPECT_THROW(io.Open("foo.bp", adios2::Mode::Read), std::exception);
    std::remove("foo.bp");
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
