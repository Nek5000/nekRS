/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cstdint>
#include <cstring>
#include <ctime>

#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <adios2.h>

#include <gtest/gtest.h>

#include "TestData.h"

#include "ParseArgs.h"
int StartStep = 0;
int EndStep = 4 * 4 * 4 * 4 * 4; // all 4 possibilities for all 5 variables
int SmallSize = 100;

// ADIOS2  write
TEST(CommonWriteTest, ADIOS2CommonWrite)
{
    adios2::ADIOS adios;

    adios2::IO io = adios.DeclareIO("TestIO");

    adios2::Dims big_shape{static_cast<unsigned int>(DataSize)};
    adios2::Dims big_start{static_cast<unsigned int>(0)};
    adios2::Dims big_count{static_cast<unsigned int>(DataSize)};

    adios2::Dims small_shape{static_cast<unsigned int>(SmallSize)};
    adios2::Dims small_start{static_cast<unsigned int>(0)};
    adios2::Dims small_count{static_cast<unsigned int>(SmallSize)};

    std::vector<adios2::Variable<double>> vars(5);
    vars[0] = io.DefineVariable<double>("big1", big_shape, big_start, big_count);
    vars[1] = io.DefineVariable<double>("small1", small_shape, small_start, small_count);
    vars[2] = io.DefineVariable<double>("big2", big_shape, big_start, big_count);
    vars[3] = io.DefineVariable<double>("small2", small_shape, small_start, small_count);
    vars[4] = io.DefineVariable<double>("big3", big_shape, big_start, big_count);

    std::vector<std::vector<double>> data(5);
    for (int i = 0; i < 5; i++)
    {
        int size = DataSize;
        if ((i == 1) || (i == 3))
        {
            size = SmallSize;
        }
        data[i] = std::vector<double>(size);
    }

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Write);

    /*
     *   write possibilities:
     *		Don't write
     *          Sync   - always destroy data afterwards
     *		Deferred
     *		Deferred with immediate PerformPuts() or PerformDataWrite() -
     *Destroy all prior data
     *
     */
    for (int step = StartStep; step < EndStep; ++step)
    {
        int mask = step;
        engine.BeginStep();

        std::cout << "Begin Write Step " << step << " writing vars : ";
        for (int j = 0; j < 5; j++)
        {
            std::fill(data[j].begin(), data[j].end(), (double)j + 1.0);
        }
        for (int j = 0; j < 5; j++)
        {
            adios2::Mode write_mode;
            char c;
            int this_var_mask = (mask & 0x3);
            mask >>= 2;
            switch (this_var_mask)
            {
            case 0:
                continue;
            case 1:
                write_mode = adios2::Mode::Sync;
                c = 's';
                break;
            case 2:
            case 3:
                write_mode = adios2::Mode::Deferred;
                c = 'd';
                break;
            }
            std::cout << j << c << " ";
            engine.Put(vars[j], data[j].data(), write_mode);
            if (this_var_mask == 1)
            {
                std::fill(data[j].begin(), data[j].end(), -100.0);
            }
            else if (this_var_mask == 3)
            {
                if (Flush)
                {
                    std::cout << "PDW ";
                    engine.PerformDataWrite();
                }
                else
                {
                    std::cout << "P ";
                    engine.PerformPuts();
                }
                for (int k = 0; k <= j; k++)
                    std::fill(data[k].begin(), data[k].end(), -100.0);
            }
        }
        std::cout << std::endl;
        engine.EndStep();
    }

    // Close the file
    engine.Close();
}

// ADIOS2  write
TEST(CommonWriteTest, ADIOS2CommonRead)
{
    adios2::ADIOS adios;

    adios2::IO io = adios.DeclareIO("TestIO");

    std::vector<std::vector<double>> data(5);
    for (int i = 0; i < 5; i++)
    {
        int size = DataSize;
        if ((i == 1) || (i == 3))
        {
            size = SmallSize;
        }
        data[i] = (std::vector<double>(size));
    }

    // Create the Engine
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    adios2::Engine engine = io.Open(fname, adios2::Mode::Read);
    EXPECT_TRUE(engine);

    /*
     *   write possibilities:
     *		Don't write
     *          Sync   - always destroy data afterwards
     *		Deferred
     *		Deferred with immediate PerformPuts()  -  Destroy all prior data
     *
     */
    for (int step = StartStep; step < EndStep; ++step)
    {
        EXPECT_TRUE(engine.BeginStep() == adios2::StepStatus::OK);

        std::vector<adios2::Variable<double>> vars(5);
        vars[0] = io.InquireVariable<double>("big1");
        vars[1] = io.InquireVariable<double>("small1");
        vars[2] = io.InquireVariable<double>("big2");
        vars[3] = io.InquireVariable<double>("small2");
        vars[4] = io.InquireVariable<double>("big3");

        std::vector<bool> var_present(vars.size());
        for (size_t i = 0; i < data.size(); i++)
            std::fill(data[i].begin(), data[i].end(), -200.0);
        std::cout << "Variables Read in TS " << step << ": ";
        for (int j = 0; j < 5; j++)
        {
            var_present[j] = (bool)vars[j];
            if (vars[j])
            {
                std::cout << " " << j;
                engine.Get(vars[j], data[j].data());
            }
        }
        std::cout << std::endl;
        engine.EndStep();
        for (int j = 0; j < 5; j++)
        {
            if (var_present[j])
            {
                for (std::size_t index = 0; index < data[j].size(); ++index)
                {
                    EXPECT_EQ(data[j][index], j + 1.0)
                        << "Data isn't correct, for " << vars[j].Name() << "[" << index << "]";
                }
            }
        }
    }

    // Close the file
    engine.Close();
}

int main(int argc, char **argv)
{

    int result;
    ::testing::InitGoogleTest(&argc, argv);

    DelayMS = 0; // zero for common writer

    ParseArgs(argc, argv);

    result = RUN_ALL_TESTS();

    return result;
}
