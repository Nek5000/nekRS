//
// Created by ganyush on 4/23/21.
//
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

class ADIOSInquireDefineTest : public ::testing::Test
{
public:
    ADIOSInquireDefineTest() = default;
};

TEST_F(ADIOSInquireDefineTest, Read)
{

    // Number of steps
    const int32_t NSteps = 5;
    int mpiRank = 0, mpiSize = 1;

#if ADIOS2_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    std::string filename = "ADIOSInquireDefine_MPI.bp";
#else
    std::string filename = "ADIOSInquireDefine.bp";
#endif

    // Write test data using BP
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        // Number of elements per process
        const std::size_t Nx = 10;
        std::vector<int32_t> Ints0 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        std::vector<int32_t> Ints1 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        std::vector<int32_t> Ints2 = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
        std::vector<int32_t> Ints3 = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
        std::vector<int32_t> Ints4 = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
        try
        {
            adios2::IO ioWrite = adios.DeclareIO("TestIOWrite");
            adios2::Engine engine = ioWrite.Open(filename, adios2::Mode::Write);
            adios2::Dims shape{static_cast<unsigned int>(mpiSize * Nx)};
            adios2::Dims start{static_cast<unsigned int>(mpiRank * Nx)};
            adios2::Dims count{static_cast<unsigned int>(Nx)};

            auto var_g = ioWrite.DefineVariable<int32_t>("global_variable", shape, start, count);
            for (auto step = 0; step < NSteps; ++step)
            {
                engine.BeginStep();
                if (step == 0)
                {
                    engine.Put(var_g, Ints0.data());
                    auto var0 = ioWrite.DefineVariable<int32_t>("variable0", shape, start, count);
                    engine.Put(var0, Ints0.data(), adios2::Mode::Deferred);
                    if (!ioWrite.InquireVariable<int>("variable0"))
                    {
                        auto var0 =
                            ioWrite.DefineVariable<int32_t>("variable0", shape, start, count);
                        engine.Put(var0, Ints1.data(), adios2::Mode::Deferred);
                    }
                }
                else if (step == 1)
                {
                    engine.Put(var_g, Ints1.data());
                    if (!ioWrite.InquireVariable<int>("variable1"))
                    {
                        auto var1 =
                            ioWrite.DefineVariable<int32_t>("variable1", shape, start, count);
                        engine.Put(var1, Ints1.data(), adios2::Mode::Deferred);
                    }
                }
                else if (step == 2)
                {
                    engine.Put(var_g, Ints2.data());
                    if (!ioWrite.InquireVariable<int>("variable2"))
                    {
                        auto var2 =
                            ioWrite.DefineVariable<int32_t>("variable2", shape, start, count);
                        engine.Put(var2, Ints2.data(), adios2::Mode::Deferred);
                    }
                }
                else if (step == 3)
                {
                    engine.Put(var_g, Ints3.data());
                    if (!ioWrite.InquireVariable<int>("variable3"))
                    {
                        auto var3 =
                            ioWrite.DefineVariable<int32_t>("variable3", shape, start, count);
                        engine.Put(var3, Ints3.data(), adios2::Mode::Deferred);
                    }
                }
                else if (step == 4)
                {
                    engine.Put(var_g, Ints4.data());
                }
                engine.EndStep();
            }
            engine.Close();
        }
        catch (std::exception &e)
        {
            std::cout << "Exception " << e.what() << std::endl;
            EXPECT_TRUE(false);
        }

#if ADIOS2_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        try
        {
            adios2::IO ioRead = adios.DeclareIO("TestIORead");
            ioRead.SetEngine("BPFile");
            adios2::Engine engine_s = ioRead.Open(filename, adios2::Mode::Read);
            EXPECT_TRUE(engine_s);

            adios2::Variable<int> var_gr;
            for (auto step = 0; step < NSteps; step++)
            {
                engine_s.BeginStep();
                var_gr = ioRead.InquireVariable<int>("global_variable");
                auto var0r = ioRead.InquireVariable<int32_t>("variable0");
                auto var1r = ioRead.InquireVariable<int32_t>("variable1");
                auto var2r = ioRead.InquireVariable<int32_t>("variable2");
                auto var3r = ioRead.InquireVariable<int32_t>("variable3");

                if (step == 0)
                {
                    EXPECT_TRUE(var0r);
                    EXPECT_FALSE(var1r);
                    EXPECT_FALSE(var2r);
                    EXPECT_FALSE(var3r);
                    EXPECT_TRUE(var_gr);
                }
                else if (step == 1)
                {
                    EXPECT_FALSE(var0r);
                    EXPECT_TRUE(var1r);
                    EXPECT_FALSE(var2r);
                    EXPECT_FALSE(var3r);
                    EXPECT_TRUE(var_gr);
                }
                else if (step == 2)
                {
                    EXPECT_FALSE(var0r);
                    EXPECT_FALSE(var1r);
                    EXPECT_TRUE(var2r);
                    EXPECT_FALSE(var3r);
                    EXPECT_TRUE(var_gr);
                }
                else if (step == 3)
                {
                    EXPECT_FALSE(var0r);
                    EXPECT_FALSE(var1r);
                    EXPECT_FALSE(var2r);
                    EXPECT_TRUE(var3r);
                    EXPECT_TRUE(var_gr);
                }
                else if (step == 4)
                {
                    EXPECT_TRUE(var_gr);
                }
                if (var0r)
                {
                    std::vector<int32_t> res;
                    var0r.SetSelection({{Nx * mpiRank}, {Nx}});
                    engine_s.Get<int32_t>(var0r, res, adios2::Mode::Sync);
                    EXPECT_EQ(res, Ints0);
                }
                if (var1r)
                {
                    std::vector<int32_t> res;
                    var1r.SetSelection({{Nx * mpiRank}, {Nx}});
                    engine_s.Get<int32_t>(var1r, res, adios2::Mode::Sync);
                    EXPECT_EQ(res, Ints1);
                }
                if (var2r)
                {
                    std::vector<int32_t> res;
                    var2r.SetSelection({{Nx * mpiRank}, {Nx}});
                    engine_s.Get<int32_t>(var2r, res, adios2::Mode::Sync);
                    EXPECT_EQ(res, Ints2);
                }
                if (var3r)
                {
                    std::vector<int32_t> res;
                    var3r.SetSelection({{Nx * mpiRank}, {Nx}});
                    engine_s.Get<int32_t>(var3r, res, adios2::Mode::Sync);
                    EXPECT_EQ(res, Ints3);
                }
                if (var_gr)
                {
                    std::vector<int32_t> res;
                    var_gr.SetSelection({{Nx * mpiRank}, {Nx}});
                    engine_s.Get<int32_t>(var_gr, res, adios2::Mode::Sync);
                    EXPECT_EQ(res, std::vector<int32_t>(10, step));
                }
                engine_s.EndStep();
            }
            engine_s.Close();
        }
        catch (std::exception &e)
        {
            std::cout << "Exception " << e.what() << std::endl;
            EXPECT_TRUE(false);
        }
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
    int result = RUN_ALL_TESTS();
#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
