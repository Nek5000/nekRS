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

void printBlock(const adios2::Box<adios2::Dims> &block, const unsigned int blockID)
{
    std::cout << "        block " << blockID << " start = ";
    printVector(block.first);
    std::cout << "  count = ";
    printVector(block.second);
    std::cout << std::endl;
}

void printBlock(const adios2::Box<adios2::Dims> &block, const unsigned int blockID,
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

void assert_block(const adios2::Box<adios2::Dims> &block, const unsigned int blockID,
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

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_1D_100)
{
    size_t blockSize;
    adios2::Dims count(1, 100);

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

    {
        blockSize = 100;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 1);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 1);
        ASSERT_EQ(subBlockInfo.Div[0], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 1);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {100});
    }

    {
        blockSize = 50;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 2);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 1);
        ASSERT_EQ(subBlockInfo.Div[0], 2);
        ASSERT_EQ(subBlockInfo.Rem.size(), 1);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {50});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {50}, {50});
    }

    {
        blockSize = 17;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 6);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 1);
        ASSERT_EQ(subBlockInfo.Div[0], 6);
        ASSERT_EQ(subBlockInfo.Rem.size(), 1);
        ASSERT_EQ(subBlockInfo.Rem[0], 4);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);

        blockID = 0;
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
        assert_block(block, blockID, {84}, {16});
    }

    {
        blockSize = 170;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 1);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 1);
        ASSERT_EQ(subBlockInfo.Div[0], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 1);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0}, {100});
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

        blockID = 0;
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
        assert_block(block, blockID, {99}, {1});
    }
}

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_2D_10x10)
{
    size_t blockSize;
    adios2::Dims count(2, 10);

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

    {
        blockSize = 100;
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

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {10, 10});
    }

    {
        blockSize = 50;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 2);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 2);
        ASSERT_EQ(subBlockInfo.Div[0], 2);
        ASSERT_EQ(subBlockInfo.Div[1], 1);
        ASSERT_EQ(subBlockInfo.Rem.size(), 2);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 1);

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {5, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {5, 0}, {5, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 5}, {1, 5});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 0}, {1, 10});
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

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0}, {10, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 9}, {1, 1});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 7}, {1, 3});
    }
}

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_3D_10x10x10)
{
    size_t blockSize;
    adios2::Dims count(3, 10);

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

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

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 10, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 10, 10});

        blockID = 9;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {9, 0, 0}, {1, 10, 10});
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

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {5, 10, 10});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {5, 0, 0}, {5, 10, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 5, 0}, {1, 5, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 8, 0}, {1, 2, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 0, 0}, {1, 10, 10});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 9, 9}, {1, 1, 1});
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

        blockID = 0;
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
        assert_block(block, blockID, {9, 9, 7}, {1, 1, 3});
    }
}

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_3D_24x24x48)
{
    size_t blockSize;
    adios2::Dims count({24, 24, 48});

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

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

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {4, 24, 48});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {4, 0, 0}, {4, 24, 48});

        blockID = 5;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {20, 0, 0}, {4, 24, 48});
    }
}

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_3D_3x2x5)
{
    size_t blockSize;
    adios2::Dims count({3, 2, 5});

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

    {
        blockSize = 2;
        printTestInfo(count, blockSize, __LINE__);

        struct adios2::helper::BlockDivisionInfo subBlockInfo = adios2::helper::DivideBlock(
            count, blockSize, adios2::helper::BlockDivisionMethod::Contiguous);
        printBlockDivisionInfo(subBlockInfo);
        ASSERT_EQ(subBlockInfo.NBlocks, 12);
        ASSERT_EQ(subBlockInfo.SubBlockSize, blockSize);
        ASSERT_EQ(subBlockInfo.Div.size(), 3);
        ASSERT_EQ(subBlockInfo.Div[0], 3);
        ASSERT_EQ(subBlockInfo.Div[1], 2);
        ASSERT_EQ(subBlockInfo.Div[2], 2);
        ASSERT_EQ(subBlockInfo.Rem.size(), 3);
        ASSERT_EQ(subBlockInfo.Rem[0], 0);
        ASSERT_EQ(subBlockInfo.Rem[1], 0);
        ASSERT_EQ(subBlockInfo.Rem[2], 1);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct.size(), 3);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[0], 4);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[1], 2);
        ASSERT_EQ(subBlockInfo.ReverseDivProduct[2], 1);

        blockID = 0;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 0}, {1, 1, 3});

        blockID = 1;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 0, 3}, {1, 1, 2});

        blockID = 2;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 0}, {1, 1, 3});

        blockID = 3;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {0, 1, 3}, {1, 1, 2});

        blockID = 4;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {1, 0, 0}, {1, 1, 3});

        blockID = 11;
        block = adios2::helper::GetSubBlock(count, subBlockInfo, blockID);
        assert_block(block, blockID, {2, 1, 3}, {1, 1, 2});
    }
}

TEST(ADIOS2DivideBlock, ADIOS2DivideBlock_4D_3x2x5x4)
{
    size_t blockSize;
    adios2::Dims count({3, 2, 5, 4});

    adios2::Box<adios2::Dims> block;
    unsigned int blockID;

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

        blockID = 0;
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
        assert_block(block, blockID, {2, 1, 4, 0}, {1, 1, 1, 4});
    }
}

int main(int argc, char **argv)
{

    int result;
    ::testing::InitGoogleTest(&argc, argv);
    result = RUN_ALL_TESTS();

    return result;
}
