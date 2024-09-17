/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <adios2.h>
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

std::string engineName;

class ADIOSInquireVariableException : public ::testing::Test
{
public:
    ADIOSInquireVariableException() = default;
};

TEST_F(ADIOSInquireVariableException, Read)
{

    // Number of steps
    const std::size_t NSteps = 5;

    int rank = 0, size = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::string filename = "ADIOSInquireVariableException";
#else
    std::string filename = "ADIOSInquireVariableException";
#endif

    // Write test data using BP
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io_w = adios.DeclareIO("Test");
        std::cout << "Setting engine name to " << engineName << std::endl;
        io_w.SetEngine(engineName);

        adios2::Engine writer = io_w.Open(filename, adios2::Mode::Write);
        const std::size_t Nx = 10;
        const adios2::Dims shape = {size * Nx};
        const adios2::Dims start = {rank * Nx};
        const adios2::Dims count = {Nx};

        auto var1 = io_w.DefineVariable<int32_t>("variable1", shape, start, count);

        for (size_t step = 0; step < NSteps; ++step)
        {
            std::vector<int32_t> Ints(10, static_cast<int32_t>(step));
            writer.BeginStep();
            writer.Put(var1, Ints.data());
            writer.EndStep();
        }
        writer.Close();

        adios2::IO io_r = adios.DeclareIO("TestIO");
        io_r.SetParameters({{"StreamReader", "true"}});
        auto reader = io_r.Open(filename, adios2::Mode::Read);
        auto var = io_r.InquireVariable<int32_t>("variable1");
        std::cout << "Engine type is " << reader.Type() << std::endl;
        for (size_t step = 0; step < NSteps; step++)
        {
            reader.BeginStep();
            std::vector<int32_t> myInts;
            std::cout << "in step " << step << std::endl;
            EXPECT_THROW(reader.Get<int32_t>(var, myInts, adios2::Mode::Sync);
                         , std::invalid_argument);
            reader.EndStep();
        }
        reader.Close();
    }
}
int main(int argc, char **argv)
{
#if ADIOS2_USE_MPI
    int provided;

    // MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
    MPI_Init_thread(nullptr, nullptr, MPI_THREAD_MULTIPLE, &provided);
#endif
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }

    int result = RUN_ALL_TESTS();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
