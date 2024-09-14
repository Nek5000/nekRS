/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>

#include <limits>
#include <stdexcept>

#include <adios2.h>
#include <adios2/common/ADIOSTypes.h>

#include <gtest/gtest.h>

std::string engineName; // comes from command line

class ADIOSReadDirectIOTest : public ::testing::Test
{
public:
    ADIOSReadDirectIOTest() = default;
};

TEST_F(ADIOSReadDirectIOTest, BufferResize)
{
    /* Test proper buffer sizes when one chunk cannot hold all variables,
       and the last chunck is resized back. It should be properly aligned
       to not cause any problems at writing that chunk.
    */

    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    std::string filename = "ADIOSDirectIO_MPI.bp";
#else
    std::string filename = "ADIOSDirectIO.bp";
#endif

    // Write test data using BP
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
        ioWrite.SetEngine(engineName);
        ioWrite.SetParameter("DirectIO", "true");
        ioWrite.SetParameter("DirectIOAlignOffset", "4096");
        ioWrite.SetParameter("DirectIOAlignBuffer", "4096");
        ioWrite.SetParameter("StripeSize", "9999");
        ioWrite.SetParameter("BufferChunkSize", "7111");
        // BufferChunkSize should be adjusted to 2*4096 by engine
        // StripeSize should be adjusted to 3*4096 by engine

        adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
        // Number of elements per process
        const std::size_t Nx = 2000;
        // two variables fit in one Chunk but not three
        adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
        adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
        adios2::Dims count{static_cast<unsigned int>(Nx)};

        auto var0 = ioWrite.DefineVariable<char>("var0", shape, start, count);
        auto var1 = ioWrite.DefineVariable<char>("var1", shape, start, count);
        auto var2 = ioWrite.DefineVariable<char>("var2", shape, start, count);

        std::vector<char> a0(Nx, 'a');
        std::vector<char> a1(Nx, 'b');
        std::vector<char> a2(Nx, 'c');

        engine.BeginStep();
        engine.Put(var0, a0.data());
        engine.Put(var1, a1.data());
        engine.Put(var2, a2.data());
        engine.EndStep();
        engine.Close();

#if ADIOS2_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        adios2::IO ioRead = adios.DeclareIO("TestIORead");
        ioRead.SetEngine(engineName);
        adios2::Engine engine_s = ioRead.Open(filename, adios2::Mode::Read);
        EXPECT_TRUE(engine_s);
        try
        {
            engine_s.BeginStep();
            adios2::Variable<char> var0 = ioRead.InquireVariable<char>("var0");
            adios2::Variable<char> var1 = ioRead.InquireVariable<char>("var1");
            adios2::Variable<char> var2 = ioRead.InquireVariable<char>("var2");

            EXPECT_TRUE(var0);
            EXPECT_TRUE(var1);
            EXPECT_TRUE(var2);

            std::vector<char> res0;
            var0.SetSelection({{Nx * mpiRank}, {Nx}});
            engine_s.Get<char>(var0, res0, adios2::Mode::Sync);
            EXPECT_EQ(res0, a0);

            std::vector<char> res1;
            var1.SetSelection({{Nx * mpiRank}, {Nx}});
            engine_s.Get<char>(var1, res1, adios2::Mode::Sync);
            EXPECT_EQ(res1, a1);

            std::vector<char> res2;
            var2.SetSelection({{Nx * mpiRank}, {Nx}});
            engine_s.Get<char>(var2, res2, adios2::Mode::Sync);
            EXPECT_EQ(res2, a2);

            engine_s.EndStep();
        }
        catch (std::exception &e)
        {
            std::cout << "Exception " << e.what() << std::endl;
        }
        engine_s.Close();
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
