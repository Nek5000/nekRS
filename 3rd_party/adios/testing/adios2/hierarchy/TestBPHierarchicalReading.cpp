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

class ADIOSHierarchicalReadVariableTest : public ::testing::Test
{
public:
    ADIOSHierarchicalReadVariableTest() = default;
};

TEST_F(ADIOSHierarchicalReadVariableTest, Read)
{
    // Number of steps
    const std::size_t NSteps = 2;

    long unsigned int rank, size;

#if ADIOS2_USE_MPI
    std::string filename = "ADIOSHierarchicalReadVariable." + engineName + "_MPI.bp";

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    std::string filename = "ADIOSHierarchicalReadVariable." + engineName + ".bp";

    rank = 0;
    size = 1;
#endif

    // Write test data using BP
    {
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif
        adios2::IO io = adios.DeclareIO("TestIO");
        io.SetEngine("BPFile");

        io.AddTransport("file");
        adios2::Engine engine = io.Open(filename, adios2::Mode::Write);
        const std::size_t Nx = 10;
        const adios2::Dims shape = {size * Nx};
        const adios2::Dims start = {rank * Nx};
        const adios2::Dims count = {Nx};

        auto var1 = io.DefineVariable<int32_t>("group1/group2/group3/group4/variable1", shape,
                                               start, count);
        auto var2 = io.DefineVariable<int32_t>("group1/group2/group3/group4/variable2", shape,
                                               start, count);
        auto var3 = io.DefineVariable<int32_t>("group1/group2/group3/group4/variable3", shape,
                                               start, count);
        auto var4 = io.DefineVariable<int32_t>("group1/group2/group3/group4/variable4", shape,
                                               start, count);
        auto var5 = io.DefineVariable<int32_t>("group1/group2/group3/group4/variable5", shape,
                                               start, count);
        auto var6 = io.DefineVariable<int32_t>("variable6", shape, start, count);
        std::vector<int32_t> Ints = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

        for (size_t step = 0; step < NSteps; ++step)
        {
            engine.BeginStep();

            engine.Put(var1, Ints.data());
            engine.Put(var2, Ints.data());
            engine.Put(var3, Ints.data());
            engine.Put(var4, Ints.data());
            engine.Put(var5, Ints.data());
            engine.Put(var6, Ints.data());

            engine.EndStep();
        }
        auto EngineName = engine.Type();

        engine.Close();

        if (EngineName == "BP5Writer")
        {
            engine = io.Open(filename, adios2::Mode::Read);
            EXPECT_THROW(engine.BeginStep(), std::logic_error);
        }
        adios2::IO io2 = adios.DeclareIO("ReadIO");
        io2.SetEngine("BPFile");

        io2.AddTransport("file");
        engine = io2.Open(filename, adios2::Mode::Read);

        for (size_t step = 0; step < NSteps; step++)
        {
            engine.BeginStep();
            auto g = io2.InquireGroup('/');
            auto res = g.AvailableGroups();
            EXPECT_EQ(res[0], "group1");
            res = g.AvailableVariables();
            EXPECT_EQ(res[0], "variable6");
            g.setPath("group1/group2");
            res = g.AvailableGroups();
            EXPECT_EQ(res[0], "group3");
            g.setPath("group1/group2/group3");
            res = g.AvailableGroups();
            EXPECT_EQ(res[0], "group4");
            g.setPath("group1/group2/group3/group4");
            res = g.AvailableGroups();
            EXPECT_EQ(res.size(), 0);
            res = g.AvailableVariables();
            EXPECT_EQ(res.size(), 5);
            res = g.AvailableAttributes();
            EXPECT_EQ(res.size(), 0);

            g = io2.InquireGroup('/');
            auto var = g.InquireVariable<int32_t>("variable6");
            EXPECT_TRUE(var);
            if (var)
            {
                std::vector<int32_t> myInts;
                var.SetSelection({{Nx * rank}, {Nx}});
                engine.Get<int32_t>(var, myInts, adios2::Mode::Sync);
                EXPECT_EQ(Ints, myInts);
            }
            g.setPath("group1/group2/group3/group4");
            var = g.InquireVariable<int32_t>("variable1");
            EXPECT_TRUE(var);
            if (var)
            {
                std::vector<int32_t> myInts;
                var.SetSelection({{Nx * rank}, {Nx}});
                engine.Get<int32_t>(var, myInts, adios2::Mode::Sync);

                EXPECT_EQ(Ints, myInts);
            }
            engine.EndStep();
        }
        engine.Close();
    }
}
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    if (argc > 1)
    {
        engineName = std::string(argv[1]);
    }
    return RUN_ALL_TESTS();
}
