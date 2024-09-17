/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferV.cpp
 *
 */

#include "ChunkV.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/toolkit/format/buffer/BufferV.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <stddef.h> // max_align_t
#include <string.h>

namespace adios2
{
namespace format
{

ChunkV::ChunkV(const std::string type, const bool AlwaysCopy, const size_t MemAlign,
               const size_t MemBlockSize, const size_t ChunkSize)
: BufferV(type, AlwaysCopy, MemAlign, MemBlockSize), m_ChunkSize(ChunkSize)
{
}

ChunkV::~ChunkV()
{
    for (const auto &Chunk : m_Chunks)
    {
        free(Chunk.AllocatedPtr);
    }
}

size_t ChunkV::ChunkAlloc(Chunk &v, const size_t size)
{
    // try to alloc/realloc a buffer to requested size
    // first, size must be aligned with block size
    size_t actualsize = size;
    size_t rem = size % m_MemBlockSize;
    if (rem)
    {
        actualsize = actualsize + (m_MemBlockSize - rem);
    }

    // align usable buffer to m_MemAlign bytes
    void *b = realloc(v.AllocatedPtr, actualsize + m_MemAlign - 1);
    if (b)
    {
        if (b != v.AllocatedPtr)
        {
            v.AllocatedPtr = b;
            size_t p = (size_t)v.AllocatedPtr;
            v.Ptr = (char *)((p + m_MemAlign - 1) & ~(m_MemAlign - 1));
        }
        v.Size = actualsize;
        return actualsize;
    }
    else
    {
        std::cout << "ADIOS2 ERROR: Cannot (re)allocate " << actualsize
                  << " bytes for a chunk in ChunkV. "
                     "Continue buffering with chunk size "
                  << v.Size / 1048576 << " MB" << std::endl;
        return 0;
    }
}

size_t ChunkV::AddToVec(const size_t size, const void *buf, size_t align, bool CopyReqd,
                        MemorySpace MemSpace)
{
    AlignBuffer(align); // may call back AddToVec recursively
    size_t retOffset = CurOffset;

    if (size == 0)
    {
        return CurOffset;
    }

    /* std::cout << "  AddToVec: size = " << size << " copy = " << CopyReqd
              << std::endl; */
    if (!CopyReqd && !m_AlwaysCopy)
    {
        /*std::cout << "    add external = " << size << std::endl;*/
        // just add buf to internal version of output vector
        VecEntry entry = {true, buf, 0, size};
        DataV.push_back(entry);
    }
    else
    {
        // we can possibly append this entry to the last if the last was
        // internal
        bool AppendPossible =
            DataV.size() && !DataV.back().External &&
            (m_TailChunk->Ptr + m_TailChunkPos - DataV.back().Size == DataV.back().Base);

        if (AppendPossible && (m_TailChunkPos + size > m_ChunkSize))
        {
            // No room in current chunk, close it out
            // realloc down to used size (helpful?) and set size in array
            /*std::cout << "    downsize ptr = " << m_Chunks.back().Ptr
                      << " to size = " << m_TailChunkPos << std::endl;*/
            Chunk &c = m_Chunks.back();
            size_t actualsize = ChunkAlloc(c, m_TailChunkPos);
            size_t alignment = actualsize - m_TailChunkPos;
            if (alignment)
            {
                auto p = c.Ptr + m_TailChunkPos;
                std::fill(p, p + alignment, 0);
            }
            retOffset += alignment;
            // Update entry in DataV as size and potentiall ptr has changed
            // Learned from sanitizer: downsizing realloc still may change pointer
            VecEntry &dv = DataV.back();
            dv.Size = actualsize;
            dv.Base = c.Ptr;
            m_TailChunkPos = 0;
            m_TailChunk = nullptr;
            AppendPossible = false;
        }
        if (AppendPossible)
        {
            // We can use current chunk, just append the data;
            /*std::cout << "    append to current chunk ptr = "
                      << m_TailChunk->Ptr << " at pos = " << m_TailChunkPos
                      << " add size = " << size << std::endl; */
            CopyDataToBuffer(size, buf, m_TailChunkPos, MemSpace);
            DataV.back().Size += size;
            m_TailChunkPos += size;
        }
        else
        {
            // We need a new chunk, get the larger of size or m_ChunkSize
            /*  std::cout << "    get a new chunk size = " << m_ChunkSize
                      << std::endl; */
            size_t NewSize = m_ChunkSize;
            if (size > m_ChunkSize)
                NewSize = size;
            Chunk c{nullptr, nullptr, 0};
            ChunkAlloc(c, NewSize);
            m_Chunks.push_back(c);
            m_TailChunk = &m_Chunks.back();
            CopyDataToBuffer(size, buf, 0, MemSpace);
            m_TailChunkPos = size;
            VecEntry entry = {false, m_TailChunk->Ptr, 0, size};
            DataV.push_back(entry);
        }
    }
    CurOffset = retOffset + size;
    return retOffset;
}

void ChunkV::CopyDataToBuffer(const size_t size, const void *buf, size_t pos, MemorySpace MemSpace)
{
#ifdef ADIOS2_HAVE_GPU_SUPPORT
    if (MemSpace == MemorySpace::GPU)
    {
        helper::CopyFromGPUToBuffer(m_TailChunk->Ptr, pos, buf, MemSpace, size);
        return;
    }
#endif
    memcpy(m_TailChunk->Ptr + pos, buf, size);
}

BufferV::BufferPos ChunkV::Allocate(const size_t size, size_t align)
{
    /*std::cout << "  Allocate: size = " << size << " align = " << align
              << std::endl;*/
    if (size == 0)
    {
        return BufferPos(-1, 0, CurOffset);
    }

    AlignBuffer(align);

    // we can possibly append this entry to the last if the last was
    // internal
    bool AppendPossible =
        DataV.size() && !DataV.back().External &&
        (m_TailChunk->Ptr + m_TailChunkPos - DataV.back().Size == DataV.back().Base);

    if (AppendPossible && (m_TailChunkPos + size > m_ChunkSize))
    {
        // No room in current chunk, close it out
        // realloc down to used size (helpful?) and set size in array
        Chunk &c = m_Chunks.back();
        size_t actualsize = ChunkAlloc(c, m_TailChunkPos);
        size_t alignment = actualsize - m_TailChunkPos;
        if (alignment)
        {
            auto p = c.Ptr + m_TailChunkPos;
            std::fill(p, p + alignment, 0);
            CurOffset += alignment;
        }
        // Update entry in DataV as size and potentiall ptr has changed
        // Learned from sanitizer: downsizing realloc still may change pointer
        VecEntry &dv = DataV.back();
        dv.Size = actualsize;
        dv.Base = c.Ptr;
        m_TailChunkPos = 0;
        m_TailChunk = nullptr;
        AppendPossible = false;
    }

    size_t bufferPos = 0;
    if (AppendPossible)
    {
        // We can use current chunk, just append the data;
        bufferPos = m_TailChunkPos;
        DataV.back().Size += size;
        m_TailChunkPos += size;
    }
    else
    {
        // We need a new chunk, get the larger of size or m_ChunkSize
        size_t NewSize = m_ChunkSize;
        if (size > m_ChunkSize)
            NewSize = size;
        Chunk c{nullptr, nullptr, 0};
        ChunkAlloc(c, NewSize);
        m_Chunks.push_back(c);
        m_TailChunk = &m_Chunks.back();
        bufferPos = 0;
        m_TailChunkPos = size;
        VecEntry entry = {false, m_TailChunk->Ptr, 0, size};
        DataV.push_back(entry);
    }

    BufferPos bp(static_cast<int>(DataV.size() - 1), bufferPos, CurOffset);
    // valid ptr anytime <-- DataV[idx] + bufferPos;

    CurOffset += size;

    return bp;
}

void ChunkV::DownsizeLastAlloc(const size_t oldSize, const size_t newSize)
{
    DataV.back().Size -= (oldSize - newSize);
    CurOffset -= (oldSize - newSize);
    m_TailChunkPos -= (oldSize - newSize);
}

void *ChunkV::GetPtr(int bufferIdx, size_t posInBuffer)
{
    if (bufferIdx == -1)
    {
        return nullptr;
    }
    else if (static_cast<size_t>(bufferIdx) > DataV.size() || DataV[bufferIdx].External)
    {
        helper::Throw<std::invalid_argument>(
            "Toolkit", "format::ChunkV", "GetPtr",
            "ChunkV::GetPtr(" + std::to_string(bufferIdx) + ", " + std::to_string(posInBuffer) +
                ") refers to a non-existing or deferred memory chunk.");
        return nullptr;
    }
    else
    {
        return (void *)((char *)DataV[bufferIdx].Base + posInBuffer);
    }
}

void *ChunkV::GetPtr(size_t OverallPosInBuffer)
{
    int bufferIdx = 0;
    if (DataV.size() == 0)
        return nullptr;
    while (DataV[bufferIdx].Size <= OverallPosInBuffer)
    {
        OverallPosInBuffer -= DataV[bufferIdx].Size;
        bufferIdx++;
        if (static_cast<size_t>(bufferIdx) > DataV.size())
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "format::ChunkV", "GetPtr",
                "ChunkV::GetPtr(" + std::to_string(OverallPosInBuffer) +
                    ") refers to a non-existing or deferred memory chunk.");
            return nullptr;
        }
    }
    if (DataV[bufferIdx].External)
        return ((char *)DataV[bufferIdx].External) + OverallPosInBuffer;

    return (void *)((char *)DataV[bufferIdx].Base + OverallPosInBuffer);
}

std::vector<core::iovec> ChunkV::DataVec() noexcept
{
    std::vector<core::iovec> iov(DataV.size());
    for (std::size_t i = 0; i < DataV.size(); ++i)
    {
        // For ChunkV, all entries in DataV are actual iov entries.
        iov[i].iov_base = DataV[i].Base;
        iov[i].iov_len = DataV[i].Size;
    }
    return iov;
}

} // end namespace format
} // end namespace adios2
