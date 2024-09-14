#include <adios2.h>
#include <array>
#include <condition_variable>
#include <future>
#include <iostream>
#include <numeric>
#include <thread>

#include <gtest/gtest.h>

#include "ParseArgs.h"

using dt = long long;

int value_errors = 0;

std::mutex StdOutMtx;

int Read(std::string BaseName, int ID)
{
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("IO");

    {
        std::lock_guard<std::mutex> guard(StdOutMtx);
        std::cout << "Reader: engine = " << engine << std::endl;
    }
    io.SetEngine(engine);
    io.SetParameters(engineParams);

    {
        std::lock_guard<std::mutex> guard(StdOutMtx);
        std::cout << "Reader: call Open" << std::endl;
    }

    try
    {
        std::string FName = BaseName + std::to_string(ID);
        adios2::Engine Reader = io.Open(FName, adios2::Mode::Read);
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Reader: passed Open" << std::endl;
        }
        std::array<dt, 1000> ar;

        auto status = Reader.BeginStep();
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Reader: passed BeginStep";
        }
        if (status == adios2::StepStatus::EndOfStream)
        {
            {
                std::lock_guard<std::mutex> guard(StdOutMtx);
                std::cout << " with EndOfStream " << std::endl;
            }
            return false;
        }
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << " with success " << std::endl;
        }

        adios2::Variable<dt> var = io.InquireVariable<dt>("data");
        Reader.Get(var, ar.begin());
        Reader.EndStep();
        dt expect = 0;
        for (auto &val : ar)
        {
            if (val != expect)
            {
                value_errors++;
            }
            expect++;
        }

        Reader.Close();
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Reader got " << expect << " values, " << value_errors << " were incorrect"
                      << std::endl;
        }
    }
    catch (std::exception &e)
    {
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Reader: Exception: " << e.what() << std::endl;
        }
        return false;
    }
    return true;
}

bool Write(std::string BaseName, int ID)
{
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("IO");
    io.SetEngine(engine);
    io.SetParameters(engineParams);
    {
        std::lock_guard<std::mutex> guard(StdOutMtx);
        std::cout << "Writer: engine = " << engine << std::endl;
    }
    auto var = io.DefineVariable<dt>("data", adios2::Dims{100, 10}, adios2::Dims{0, 0},
                                     adios2::Dims{100, 10});

    std::array<dt, 1000> ar;
    std::iota(ar.begin(), ar.end(), 0);

    try
    {
        std::string FName = BaseName + std::to_string(ID);
        adios2::Engine Writer = io.Open(FName, adios2::Mode::Write);

        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Writer completed Open() " << std::endl;
        }
        Writer.BeginStep();
        Writer.Put<dt>(var, ar.begin());
        Writer.EndStep();
        Writer.Close();
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Writer completed Close() " << std::endl;
        }
    }
    catch (std::exception &e)
    {
        {
            std::lock_guard<std::mutex> guard(StdOutMtx);
            std::cout << "Writer: Exception: " << e.what() << std::endl;
        }
        return false;
    }
    return true;
}

class TestThreads : public ::testing::Test
{
public:
    TestThreads() = default;
};

TEST_F(TestThreads, Basic)
{
    using namespace std;
    std::string BaseName = engine + "_P" + std::to_string(getpid());
    auto read_fut = std::async(std::launch::async, Read, BaseName, 0);
    auto write_fut = std::async(std::launch::async, Write, BaseName, 0);
    bool reader_success = read_fut.get();
    bool writer_success = write_fut.get();
    EXPECT_TRUE(reader_success);
    EXPECT_TRUE(writer_success);
    EXPECT_EQ(value_errors, 0) << "We got " << value_errors << " erroneous values at the reader";
}

//  This test tries to push up to the limits to see if we're leaking FDs, but it
//  runs slowly, commenting it out until needed.
//
// TE S T _F(TestThreads, Repeated)
// {
//     auto high_write_fut = std::async(std::launch::async, Write, 0);
//     for (int i = 0; i < 1024; i++)
//     {
//         using namespace std;
//         std::string BaseName = std::to_string(_getpid());
//         auto read_fut = std::async(std::launch::async, Read, BaseName, i +
//         1); auto write_fut = std::async(std::launch::async, Write, BaseName,
//         i + 1); bool reader_success = read_fut.get(); bool writer_success =
//         write_fut.get(); EXPECT_TRUE(reader_success);
//         EXPECT_TRUE(writer_success);
//         EXPECT_EQ(value_errors, 0)
//             << "We got " << value_errors << " erroneous values at the
//             reader";
//         std::cout << "finished pair " << i << std::endl;
//     }
//     auto high_read_fut = std::async(std::launch::async, Read, 0);
//     bool reader_success = high_read_fut.get();
//     bool writer_success = high_write_fut.get();
//     EXPECT_TRUE(reader_success);
//     EXPECT_TRUE(writer_success);
//     EXPECT_EQ(value_errors, 0);
// }

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    ParseArgs(argc, argv);

    int result;
    result = RUN_ALL_TESTS();
    return result;
}
