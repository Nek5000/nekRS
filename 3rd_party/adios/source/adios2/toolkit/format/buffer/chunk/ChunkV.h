/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_MALLOC_CHUNKV_H_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_MALLOC_CHUNKV_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/CoreTypes.h"

#include "adios2/toolkit/format/buffer/BufferV.h"

namespace adios2
{
namespace format
{

class ChunkV : public BufferV
{
public:
    uint64_t Size() noexcept;

    const size_t m_ChunkSize;

    ChunkV(const std::string type, const bool AlwaysCopy = false, const size_t MemAlign = 1,
           const size_t MemBlockSize = 1, const size_t ChunkSize = DefaultBufferChunkSize);
    virtual ~ChunkV();

    virtual std::vector<core::iovec> DataVec() noexcept;

    virtual size_t AddToVec(const size_t size, const void *buf, size_t align, bool CopyReqd,
                            MemorySpace MemSpace = MemorySpace::Host);

    virtual BufferPos Allocate(const size_t size, size_t align);
    virtual void DownsizeLastAlloc(const size_t oldSize, const size_t newSize);

    virtual void *GetPtr(int bufferIdx, size_t posInBuffer);
    virtual void *GetPtr(size_t OverallPosInBuffer);

    void CopyDataToBuffer(const size_t size, const void *buf, size_t pos, MemorySpace MemSpace);

private:
    struct Chunk
    {
        char *Ptr;          // aligned, do not free
        void *AllocatedPtr; // original ptr, free this
        size_t Size;
    };

    std::vector<Chunk> m_Chunks;
    size_t m_TailChunkPos = 0;
    Chunk *m_TailChunk = nullptr;

    // allocator function, doing aligned allocation of memory
    // return true if (re)allocation is successful
    // on failure, VecEntry is unmodified
    size_t ChunkAlloc(Chunk &v, const size_t size);
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_MALLOC_MALLOCV_H_ */
