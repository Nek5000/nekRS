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
#include <adios2/toolkit/format/buffer/chunk/ChunkV.h>

#include <gtest/gtest.h>

namespace adios2
{
namespace format
{

static void print_array(const std::string name, const uint8_t *ptr, const size_t size)
{
    std::cout << name << " ptr = " << (void *)ptr << "= [ ";
    for (size_t n = 0; n < size; ++n)
    {
        std::cout << (unsigned int)ptr[n] << " ";
    }
    std::cout << "]\n";
}

TEST(ChunkV, AlignmentToMemBlockSize)
{
    /* Align chunks to a certain byte size, here 32 */
    size_t MemBlockSize = 32;
    size_t ChunkSize = 137;
    ChunkV b = ChunkV("test", false, 1, MemBlockSize, ChunkSize);

    /* Align "allocation" to a certain byte size, like type size, here 8 */
    size_t AlignnmentSize = 8;
    /* Make allocation not aligned to either AlignmentSize nor the
     * MemBlockSize*/
    size_t AllocSize = 51;

    {
        adios2::format::BufferV::BufferPos pos = b.Allocate(AllocSize, AlignnmentSize);
        ASSERT_EQ(pos.bufferIdx, 0);
        ASSERT_EQ(pos.globalPos, 0);
        ASSERT_EQ(pos.posInBuffer, 0);
        uint8_t *ptr = reinterpret_cast<uint8_t *>(b.GetPtr(pos.bufferIdx, pos.posInBuffer));
        for (uint8_t n = 0; n < AllocSize; ++n)
        {
            ptr[n] = n;
        }
        print_array("array 1", ptr, AllocSize);

        std::vector<core::iovec> vec = b.DataVec();
        ASSERT_EQ(vec.size(), 1);
        /* First block is alone in first chunk, so only real length (51) should
         * be here */
        ASSERT_EQ(vec[0].iov_len, 51);
    }
    {
        /* second block should fit in existing chunk but aligned to 56 bytes
           (first AllocSize aligned to AlignmentSize is 56) */
        adios2::format::BufferV::BufferPos pos = b.Allocate(AllocSize, AlignnmentSize);
        ASSERT_EQ(pos.bufferIdx, 0);
        ASSERT_EQ(pos.globalPos, 56);
        ASSERT_EQ(pos.posInBuffer, 56);
        uint8_t *ptr = reinterpret_cast<uint8_t *>(b.GetPtr(pos.bufferIdx, pos.posInBuffer));
        for (uint8_t n = 0; n < AllocSize; ++n)
        {
            ptr[n] = n;
        }
        print_array("array 2", ptr, AllocSize);

        std::vector<core::iovec> vec = b.DataVec();
        ASSERT_EQ(vec.size(), 1);
        /* Seconf block is aligned to 56, so real length should be 56+51 */
        ASSERT_EQ(vec[0].iov_len, 56 + 51);
    }
    {
        /* - third block should NOT fit in existing chunk
           - first chunk should be 128 bytes long (aligned to 32)
           - globalPos should be 128 (0 in this second chunk)
        */
        adios2::format::BufferV::BufferPos pos = b.Allocate(AllocSize, AlignnmentSize);
        ASSERT_EQ(pos.bufferIdx, 1);
        ASSERT_EQ(pos.globalPos, 128);
        ASSERT_EQ(pos.posInBuffer, 0);
        uint8_t *ptr = reinterpret_cast<uint8_t *>(b.GetPtr(pos.bufferIdx, pos.posInBuffer));
        for (uint8_t n = 0; n < AllocSize; ++n)
        {
            ptr[n] = n;
        }
        print_array("array 3", ptr, AllocSize);

        std::vector<core::iovec> vec = b.DataVec();
        ASSERT_EQ(vec.size(), 2);
        /* First chunk must be finalized and aligned to MemBlockSize, so should
         * be 128 now */
        ASSERT_EQ(vec[0].iov_len, 128);
        /* Third block is alone in second chunk, so only real length (51) should
         * be here */
        ASSERT_EQ(vec[1].iov_len, 51);

        const uint8_t *chunk0 = reinterpret_cast<const uint8_t *>(vec[0].iov_base);
        const uint8_t *chunk1 = reinterpret_cast<const uint8_t *>(vec[1].iov_base);
        /* first array in chunk0[0..AllocSize-1]*/
        ASSERT_EQ(chunk0[0], 0);
        ASSERT_EQ(chunk0[AllocSize - 1], AllocSize - 1);
        /* alignment fill for 5 bytes*/
        ASSERT_EQ(chunk0[51], 0);
        ASSERT_EQ(chunk0[52], 0);
        ASSERT_EQ(chunk0[53], 0);
        ASSERT_EQ(chunk0[54], 0);
        ASSERT_EQ(chunk0[55], 0);

        /* second array in chunk0[56..107]*/
        ASSERT_EQ(chunk0[56], 0);
        ASSERT_EQ(chunk0[56 + AllocSize - 1], AllocSize - 1);
        /* alignment fill for 5 bytes*/
        ASSERT_EQ(chunk0[107], 0);
        ASSERT_EQ(chunk0[108], 0);
        ASSERT_EQ(chunk0[109], 0);
        ASSERT_EQ(chunk0[110], 0);
        ASSERT_EQ(chunk0[111], 0);
        /* aligment fill to MemBlockSize for 16 bytes */
        ASSERT_EQ(chunk0[112], 0);
        ASSERT_EQ(chunk0[127], 0);

        /* third array in chunk1[0..AllocSize-1]*/
        ASSERT_EQ(chunk1[1], 1);
        ASSERT_EQ(chunk1[AllocSize - 1], AllocSize - 1);
    }
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
