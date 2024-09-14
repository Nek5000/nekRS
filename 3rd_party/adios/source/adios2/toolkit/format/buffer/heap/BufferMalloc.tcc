/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferMalloc.tcc
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_TCC_

#include "BufferMalloc.h"

#include <memory>

#ifdef _WIN32
#pragma warning(disable : 4146) // Windows complains about unsigned minus
#endif

namespace adios2
{
namespace format
{

template <class T>
size_t BufferMalloc::Align() const noexcept
{
    // std::align implementation from llvm libc++
    // needed due to bug in gcc 4.8
    auto lf_align = [](const size_t alignment, const size_t size, void *&ptr, size_t &space) {
        if (size <= space)
        {
            const char *p1 = static_cast<char *>(ptr);
            const char *p2 = reinterpret_cast<char *>(
                reinterpret_cast<size_t>(p1 + (alignment - 1)) & -alignment);
            const size_t d = static_cast<size_t>(p2 - p1);
            if (d <= space - size)
            {
                space -= d;
            }
        }
    };

    void *currentAddress =
        reinterpret_cast<void *>(const_cast<char *>((char *)Data() + m_Position));
    size_t size = GetAvailableSize();
    lf_align(alignof(T), sizeof(T), currentAddress, size);
    return GetAvailableSize() - size;
}

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERMALLOC_TCC_ */
