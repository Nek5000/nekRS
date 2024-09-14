/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferMalloc.h
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_H_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_H_

#include "adios2/toolkit/format/buffer/Buffer.h"

#include "adios2/common/ADIOSMacros.h"

namespace adios2
{
namespace format
{

class BufferMalloc : public Buffer
{
public:
    BufferMalloc();
    ~BufferMalloc();

    char *Data() noexcept final;
    const char *Data() const noexcept final;

    void Resize(const size_t size, const std::string hint) final;

    void Reset(const bool resetAbsolutePosition, const bool zeroInitialize) final;

    size_t GetAvailableSize() const final;

    template <class T>
    size_t Align() const noexcept;

    void Delete();

    size_t DebugGetSize() const;
    size_t Size() const;

private:
    size_t m_size = 0;
    char *m_buffer = NULL;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_H_ */
