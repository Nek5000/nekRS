/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <iostream>
#include <numeric> //std::iota
#include <stdexcept>

#include <adios2.h>

#include <gtest/gtest.h>

std::string engineName;       // comes from command line
std::string engineParameters; // comes from command line

class BP4ReadOldCompressed : public ::testing::Test
{
public:
    BP4ReadOldCompressed() = default;
};

class BP4ReadOldCompressedP : public BP4ReadOldCompressed,
                              public ::testing::WithParamInterface<std::tuple<std::string, bool>>
{
protected:
    std::string GetFileName() { return std::get<0>(GetParam()); };
    bool GetSuccess() { return std::get<1>(GetParam()); };
};

// Read an old file (BP4 written with ADIOS 2.7.1) and compressed by an old
// library. Backward compatible reader should work properly if success is
// expected, otherwise a runtime error should be thrown.
TEST_P(BP4ReadOldCompressedP, Read)
{
    const std::string fname = GetFileName();
    const bool expectSuccess = GetSuccess();
    std::cout << "\nRead old file " << fname << " and expect it ";
    if (expectSuccess)
    {
        std::cout << "to succeed\n\n";
    }
    else
    {
        std::cout << "to throw a runtime error\n\n";
    }

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("ReadIO");

    if (!engineName.empty())
    {
        io.SetEngine(engineName);
    }
    if (!engineParameters.empty())
    {
        io.SetParameters(engineParameters);
    }

    adios2::Engine bpReader = io.Open(fname, adios2::Mode::ReadRandomAccess);

    auto var_a = io.InquireVariable<double>("a");
    EXPECT_TRUE(var_a);
    ASSERT_EQ(var_a.ShapeID(), adios2::ShapeID::GlobalArray);
    ASSERT_EQ(var_a.Steps(), 1);
    ASSERT_EQ(var_a.Shape().size(), 2);

    std::vector<double> R64;

    const adios2::Dims start{0, 0};
    const adios2::Dims count = var_a.Shape();
    const adios2::Box<adios2::Dims> sel(start, count);
    var_a.SetSelection(sel);
    var_a.SetStepSelection({0, 1});
    // The actual read function fails if there is no backward compatible
    // decompression
    if (expectSuccess)
    {
        bpReader.Get(var_a, R64, adios2::Mode::Sync);
        ASSERT_EQ(R64.size(), var_a.Shape()[0] * var_a.Shape()[1]);
    }
    else
    {
        ASSERT_THROW(bpReader.Get(var_a, R64, adios2::Mode::Sync), std::runtime_error);
    }
    bpReader.Close();
}

std::vector<std::tuple<std::string, bool>> GenerateParameters()
{
    std::vector<std::tuple<std::string, bool>> v;
    v.push_back(std::make_tuple("notcompressed_271.bp", true));
#ifdef ADIOS2_HAVE_BLOSC
    v.push_back(std::make_tuple("compressed_blosc_271.bp", true));
#endif
#ifdef ADIOS2_HAVE_BZIP2
    v.push_back(std::make_tuple("compressed_bzip2_271.bp", false));
#endif
    return v;
}

INSTANTIATE_TEST_SUITE_P(BP4ReadOldCompressed, BP4ReadOldCompressedP,
                         ::testing::ValuesIn(GenerateParameters()));

//******************************************************************************
// main
//******************************************************************************

int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    MPI_Init(nullptr, nullptr);
#endif

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    if (argc > 2)
    {
        engineParameters = std::string(argv[2]);
    }
    result = RUN_ALL_TESTS();

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
