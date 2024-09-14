/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cmath>
#include <cstdint>
#include <cstring>

#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>

#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"
#include <adios2.h>
#include <adios2/common/ADIOSTypes.h>
#include <adios2/helper/adiosCommDummy.h>
#include <adios2/toolkit/transport/file/FilePOSIX.h>

#include <gtest/gtest.h>

namespace adios2
{
namespace format
{

TEST(FileTransport, FailOnEOF)
{
    {
        std::vector<uint8_t> b(256, 0xef);
        helper::Comm comm = helper::CommDummy();
        std::unique_ptr<transport::FilePOSIX> w =
            std::unique_ptr<transport::FilePOSIX>(new transport::FilePOSIX(comm));

        w->Open("FailOnEOF", Mode::Write);
        w->Write((char *)b.data(), b.size());
        w->Close();
    }
    {
        std::vector<uint8_t> b(256);
        helper::Comm comm = helper::CommDummy();
        std::unique_ptr<transport::FilePOSIX> r =
            std::unique_ptr<transport::FilePOSIX>(new transport::FilePOSIX(comm));

        r->Open("FailOnEOF", Mode::Read);
        Params p = {{"FailOnEOF", "true"}};
        r->SetParameters(p);
        EXPECT_THROW(r->Read((char *)b.data(), b.size() * 2), std::ios_base::failure);
        r->Close();
    }
}

TEST(FileTransport, WaitForData)
{
    constexpr int size = 256;
    std::vector<uint8_t> b(size, 0xef);
    helper::Comm comm = helper::CommDummy();
    std::unique_ptr<transport::FilePOSIX> w =
        std::unique_ptr<transport::FilePOSIX>(new transport::FilePOSIX(comm));

    w->Open("FailOnEOF", Mode::Write);
    w->Write((char *)b.data(), b.size());
    {
        auto lf_WriteMore = [&](const transport::FilePOSIX *) {
            std::vector<uint8_t> b2(size, 0xfe);
            std::this_thread::sleep_for(std::chrono::seconds(2));
            w->Write((char *)b2.data(), size);
            std::cout << "Wrote data" << std::endl;
        };

        // write more data soon
        auto h = std::async(std::launch::async, lf_WriteMore, w.get());

        std::vector<uint8_t> b(size * 2);
        helper::Comm comm = helper::CommDummy();
        std::unique_ptr<transport::FilePOSIX> r =
            std::unique_ptr<transport::FilePOSIX>(new transport::FilePOSIX(comm));

        r->Open("FailOnEOF", Mode::Read);
        r->Read((char *)b.data(), size * 2);
        ASSERT_EQ(b[0], 0xef);
        ASSERT_EQ(b[size], 0xfe);
        r->Close();
    }
    w->Close();
}
}
}

int main(int argc, char **argv)
{

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}
