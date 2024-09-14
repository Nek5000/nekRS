/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 */
#include <cmath>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <limits>
#include <stdexcept>

#include <adios2.h>
#include <adios2/common/ADIOSTypes.h>
#include <adios2/helper/adiosMath.h>

#include <gtest/gtest.h>

template <typename T>
void printVector(const std::vector<T> &v)
{
    const size_t ndim = v.size();
    std::cout << "{";
    for (size_t d = 0; d < ndim; ++d)
    {
        std::cout << v[d];
        if (d < ndim - 1)
        {
            std::cout << ", ";
        }
        if (d > 0 && !(d % 20))
        {
            std::cout << "\n        ";
        }
    }
    std::cout << "}";
}

void printTestInfo(const adios2::Dims &count, const size_t blockSize, const size_t sourceLine)
{
    const size_t ndim = count.size();
    std::cout << "\nTest " << ndim << "-D array ";
    printVector(count);
    std::cout << " divide by blockSize " << blockSize << " at line " << sourceLine << std::endl;
}

void printBlockDivisionInfo(const adios2::helper::BlockDivisionInfo &info)
{
    std::cout << "    nblocks = " << info.NBlocks;        // << std::endl;
    std::cout << " subblockSize = " << info.SubBlockSize; // << std::endl;
    std::cout << " div = ";
    printVector(info.Div);
    // std::cout << std::endl;
    std::cout << " rem = ";
    printVector(info.Rem);
    // std::cout << std::endl;
    std::cout << " reverseDivProduct = ";
    printVector(info.ReverseDivProduct);
    std::cout << std::endl;
}

void printBlock(const adios2::Box<adios2::Dims> &block, const int blockID)
{
    std::cout << "        block " << blockID << " start = ";
    printVector(block.first);
    std::cout << "  count = ";
    printVector(block.second);
    std::cout << std::endl;
}

template <typename T>
void printMinMax(const std::vector<T> &expected_minmax, const T &expected_bmin,
                 const T &expected_bmax)
{
    std::cout << "    expected  bmin = " << expected_bmin << "  bmax = " << expected_bmax;
    std::cout << "  min-max = \n       ";
    printVector(expected_minmax);
    std::cout << std::endl;
}

void printBlock(const adios2::Box<adios2::Dims> &block, const int blockID,
                const std::vector<uint16_t> &expected_start,
                const std::vector<uint16_t> &expected_count)
{
    std::cout << "        block " << blockID << " start = ";
    printVector(block.first);
    std::cout << "  count = ";
    printVector(block.second);
    std::cout << "  expected start = ";
    printVector(expected_start);
    std::cout << "  count = ";
    printVector(expected_count);
    std::cout << std::endl;
}

void assert_block(const adios2::Box<adios2::Dims> &block, const int blockID,
                  const std::vector<uint16_t> &expected_start,
                  const std::vector<uint16_t> &expected_count)
{
    printBlock(block, blockID, expected_start, expected_count);
    const adios2::Dims &block_start = block.first;
    const adios2::Dims &block_count = block.second;
    ASSERT_EQ(block_start.size(), expected_start.size());
    const size_t nstart = expected_start.size();
    for (size_t d = 0; d < nstart; ++d)
    {
        ASSERT_EQ(block_start[d], expected_start[d]);
    }
    ASSERT_EQ(block_count.size(), expected_count.size());
    for (size_t d = 0; d < nstart; ++d)
    {
        ASSERT_EQ(block_count[d], expected_count[d]);
    }
}

template <typename T>
void assert_minmax(const size_t nBlocks, const T &bmin, const T &bmax, const std::vector<T> &minmax,
                   const T &expected_bmin, const T &expected_bmax,
                   const std::vector<T> &expected_minmax)
{
    printMinMax(expected_minmax, expected_bmin, expected_bmax);
    ASSERT_EQ(bmin, expected_bmin);
    ASSERT_EQ(bmax, expected_bmax);

    const size_t nMinMax = expected_minmax.size();
    ASSERT_EQ(nMinMax, 2 * nBlocks);
    for (size_t i = 0; i < nMinMax; ++i)
    {
        ASSERT_EQ(minmax[i], expected_minmax[i]);
    }
}

TEST(ADIOS2MinMaxs, ADIOS2MinMaxs_1D_100)
{
    std::vector<int> data(100);
    for (int i = 0; i < static_cast<int>(data.size()); ++i)
    {
        data[i] = i;
    }
    size_t blockSize;
    adios2::Dims count(1, data.size());
    adios2::Box<adios2::Dims> block;
    // int blockID;
    std::vector<int> MinMaxs;
    int Min, Max;

    {
        blockSize = 100;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {100});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 99});
    }

    {
        blockSize = 50;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {50});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {50}, {50});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 49, 50, 99});
    }

    {
        blockSize = 17;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {17});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {17}, {17});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {34}, {17});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {51}, {17});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {68}, {16});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {84}, {16});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99,
                      {0, 16, 17, 33, 34, 50, 51, 67, 68, 83, 84, 99});
    }

    {
        blockSize = 170;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {100});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 99});
    }

    {
        blockSize = 1;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 100);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 1);
        ASSERT_EQ(subBlockInfo.Div[0], 100);
        ASSERT_EQ(subBlockInfo.Rem.size(), 1);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {1});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1}, {1});

        blockID = 77;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {77}, {1});

        blockID = 99;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {99}, {1});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        std::vector<int> mm(200);
        for (int j = 0; j < 100; j++)
        {
            mm[2 * j] = j;
            mm[2 * j + 1] = j;
        }
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, mm);
    }
}

TEST(ADIOS2MinMaxs, ADIOS2MinMaxs_2D_10x10)
{
    std::vector<int> data(100);
    for (int i = 0; i < static_cast<int>(data.size()); ++i)
    {
        data[i] = i;
    }
    size_t blockSize;
    adios2::Dims count(2, 10);

    adios2::Box<adios2::Dims> block;
    // int blockID;
    std::vector<int> MinMaxs;
    int Min, Max;

    {
        blockSize = 100;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {10, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 99});
    }

    {
        blockSize = 50;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {5, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {5, 0}, {5, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 49, 50, 99});
    }

    {
        blockSize = 5;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 20);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 2);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {1, 5});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 5}, {1, 5});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0}, {1, 5});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 5}, {1, 5});

        blockID = 18;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0}, {1, 5});

        blockID = 19;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 5}, {1, 5});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99,
                      {0,  4,  5,  9,  10, 14, 15, 19, 20, 24, 25, 29, 30, 34,
                       35, 39, 40, 44, 45, 49, 50, 54, 55, 59, 60, 64, 65, 69,
                       70, 74, 75, 79, 80, 84, 85, 89, 90, 94, 95, 99});
    }

    {
        blockSize = 17;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 6);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 6);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 4);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {2, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {2, 0}, {2, 10});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {4, 0}, {2, 10});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {6, 0}, {2, 10});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {8, 0}, {1, 10});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0}, {1, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99,
                      {0, 19, 20, 39, 40, 59, 60, 79, 80, 89, 90, 99});
    }

    {
        blockSize = 170;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 1);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 1);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {10, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, {0, 99});
    }

    {
        blockSize = 1;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 100);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 10);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 10);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {1, 1});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1}, {1, 1});

        blockID = 9;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 9}, {1, 1});

        blockID = 10;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0}, {1, 1});

        blockID = 77;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {7, 7}, {1, 1});

        blockID = 99;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 9}, {1, 1});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        std::vector<int> mm(200);
        for (int j = 0; j < 100; j++)
        {
            mm[2 * j] = j;
            mm[2 * j + 1] = j;
        }
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99, mm);
    }

    {
        blockSize = 3;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 30);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 3);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 3);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {1, 4});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 4}, {1, 3});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 7}, {1, 3});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0}, {1, 4});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 4}, {1, 3});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 7}, {1, 3});

        blockID = 27;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0}, {1, 4});

        blockID = 28;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 4}, {1, 3});

        blockID = 29;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 7}, {1, 3}); */

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 99,
                      {0,  3,  4,  6,  7,  9,  10, 13, 14, 16, 17, 19, 20, 23, 24,
                       26, 27, 29, 30, 33, 34, 36, 37, 39, 40, 43, 44, 46, 47, 49,
                       50, 53, 54, 56, 57, 59, 60, 63, 64, 66, 67, 69, 70, 73, 74,
                       76, 77, 79, 80, 83, 84, 86, 87, 89, 90, 93, 94, 96, 97, 99});
    }
}

TEST(ADIOS2MinMaxs, ADIOS2MinMaxs_3D_10x10x10)
{
    std::vector<int> data(1000);
    for (int i = 0; i < static_cast<int>(data.size()); ++i)
    {
        data[i] = i;
    }
    size_t blockSize;
    adios2::Dims count(3, 10);

    adios2::Box<adios2::Dims> block;
    // int blockID;
    std::vector<int> MinMaxs;
    int Min, Max;

    {
        blockSize = 100;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 10);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 10, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 10, 10});

        blockID = 9;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0, 0}, {1, 10, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999,
                      {0,   99,  100, 199, 200, 299, 300, 399, 400, 499,
                       500, 599, 600, 699, 700, 799, 800, 899, 900, 999});
    }

    {
        blockSize = 500;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 2);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 2);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {5, 10, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {5, 0, 0}, {5, 10, 10});*/
        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999, {0, 499, 500, 999});
    }

    {
        blockSize = 50;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 20);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 2);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 5, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 5, 0}, {1, 5, 10});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 5, 10});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 5, 0}, {1, 5, 10});

        blockID = 18;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0, 0}, {1, 5, 10});

        blockID = 19;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 5, 0}, {1, 5, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999,
                      {0,   49,  50,  99,  100, 149, 150, 199, 200, 249, 250, 299, 300, 349,
                       350, 399, 400, 449, 450, 499, 500, 549, 550, 599, 600, 649, 650, 699,
                       700, 749, 750, 799, 800, 849, 850, 899, 900, 949, 950, 999});
    }

    {
        blockSize = 17;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 50);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 5);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 5);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 2, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 2, 0}, {1, 2, 10});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 4, 0}, {1, 2, 10});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 8, 0}, {1, 2, 10});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 2, 10});

        blockID = 48;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 6, 0}, {1, 2, 10});

        blockID = 49;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 8, 0}, {1, 2, 10});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        std::vector<int> mm(100);
        for (int j = 0; j < 50; j++)
        {
            mm[2 * j] = 20 * j;
            mm[2 * j + 1] = 20 * j + 19;
        }
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999, mm);
    }

    {
        blockSize = 170;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 6);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 6);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 4);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {2, 10, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {2, 0, 0}, {2, 10, 10});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {4, 0, 0}, {2, 10, 10});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {6, 0, 0}, {2, 10, 10});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {8, 0, 0}, {1, 10, 10});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0, 0}, {1, 10, 10}); */

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999,
                      {0, 199, 200, 399, 400, 599, 600, 799, 800, 899, 900, 999});
    }

    {
        blockSize = 1;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 1000);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 10);
        ASSERT_EQ(subBlockInfo.Div[2], 10);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 100);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 10);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 1, 1});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 1}, {1, 1, 1});

        blockID = 10;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 0}, {1, 1, 1});

        blockID = 17;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 7}, {1, 1, 1});

        blockID = 536;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {5, 3, 6}, {1, 1, 1});

        blockID = 999;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 9, 9}, {1, 1, 1});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        std::vector<int> mm(2000);
        for (int j = 0; j < 1000; j++)
        {
            mm[2 * j] = j;
            mm[2 * j + 1] = j;
        }
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999, mm);
    }

    {
        blockSize = 3;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 300);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 10);
        ASSERT_EQ(subBlockInfo.Div[1], 10);
        ASSERT_EQ(subBlockInfo.Div[2], 3);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 30);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 1, 4});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 4}, {1, 1, 3});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 7}, {1, 1, 3});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 0}, {1, 1, 4});

        blockID = 29;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 9, 7}, {1, 1, 3});

        blockID = 30;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 1, 4});

        blockID = 60;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {2, 0, 0}, {1, 1, 4});

        blockID = 180;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {6, 0, 0}, {1, 1, 4});

        blockID = 217;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {7, 2, 4}, {1, 1, 3});

        blockID = 299;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 9, 7}, {1, 1, 3});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        std::vector<int> mm(600);
        for (int i = 0; i < 10; i++)
        {
            for (int j = 0; j < 10; j++)
            {
                int idx = i * 10 + j;
                int minstart = 100 * i + 10 * j;
                int k = 6 * idx;
                mm[k++] = minstart;
                mm[k++] = minstart + 3;
                mm[k++] = minstart + 4;
                mm[k++] = minstart + 6;
                mm[k++] = minstart + 7;
                mm[k++] = minstart + 9;
            }
        }
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 999, mm);
    }
}

TEST(ADIOS2MinMaxs, ADIOS2MinMaxs_3D_24x24x48)
{
    std::vector<int> data(27648);
    for (int i = 0; i < static_cast<int>(data.size()); ++i)
    {
        data[i] = i;
    }
    size_t blockSize;
    adios2::Dims count({24, 24, 48});

    adios2::Box<adios2::Dims> block;
    // int blockID;

    std::vector<int> MinMaxs;
    int Min, Max;

    {
        blockSize = 5000;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 6);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 6);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Div[2], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {4, 24, 48});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {4, 0, 0}, {4, 24, 48});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {20, 0, 0}, {4, 24, 48});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 27647,
                      {0, 4607, 4608, 9215, 9216, 13823, 13824, 18431, 18432, 23039, 23040, 27647});
    }
}

TEST(ADIOS2MinMaxs, ADIOS2MinMaxs_4D_3x2x5x4)
{
    std::vector<int> data(120);
    for (int i = 0; i < static_cast<int>(data.size()); ++i)
    {
        data[i] = i;
    }
    size_t blockSize;
    adios2::Dims count({3, 2, 5, 4});

    adios2::Box<adios2::Dims> block;
    // int blockID;

    std::vector<int> MinMaxs;
    int Min, Max;

    {
        blockSize = 6;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 18);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 4);
        ASSERT_EQ(subBlockInfo.Div[0], 3);
        ASSERT_EQ(subBlockInfo.Div[1], 2);
        ASSERT_EQ(subBlockInfo.Div[2], 3);
        ASSERT_EQ(subBlockInfo.Div[3], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 4);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 2);
        ASSERT_EQ(subBlockInfo.Rem[3], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 4);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 6);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[3], 1);

        /*blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0, 0}, {1, 1, 2, 4});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 2, 0}, {1, 1, 2, 4});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 4, 0}, {1, 1, 1, 4});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 0, 0}, {1, 1, 2, 4});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 2, 0}, {1, 1, 2, 4});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 4, 0}, {1, 1, 1, 4});

        blockID = 6;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0, 0}, {1, 1, 2, 4});

        blockID = 17;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {2, 1, 4, 0}, {1, 1, 1, 4});*/

        adios2::helper::GetMinMaxSubblocks(data.data(), count, subBlockInfo, MinMaxs, Min, Max, 1);
        assert_minmax(subBlockInfo.NBlocks, Min, Max, MinMaxs, 0, 119,
                      {0,  7,  8,  15, 16, 19, 20,  27,  28,  35,  36,  39,
                       40, 47, 48, 55, 56, 59, 60,  67,  68,  75,  76,  79,
                       80, 87, 88, 95, 96, 99, 100, 107, 108, 115, 116, 119});
    }
}

int main(int argc, char **argv)
{

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}
