/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferMalloc.cpp
 *
 */

#include "BufferMalloc.h"
#include "BufferMalloc.tcc"
#include "adios2/helper/adiosLog.h"
#include "time.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace adios2
{
namespace format
{

BufferMalloc::BufferMalloc() : Buffer("BufferMalloc") {}
BufferMalloc::~BufferMalloc()
{
    if (m_size)
        free(m_buffer);
}

char *BufferMalloc::Data() noexcept { return (char *)m_buffer; }

const char *BufferMalloc::Data() const noexcept { return (const char *)m_buffer; }

void BufferMalloc::Resize(const size_t size, const std::string hint)
{
    if (size == 0)
        return;
    if (m_size == 0)
    {
        m_buffer = (char *)malloc(size);
        if (m_buffer == NULL)
            helper::ThrowNested<std::runtime_error>(
                "Toolkit::Format", "buffer::heap::BufferMalloc", "BufferSystemV",
                "buffer overflow when resizing to " + std::to_string(size) + " bytes, " + hint);
    }
    else
    {
        char *tmp = (char *)realloc(m_buffer, size);
        if (tmp == NULL)
            helper::ThrowNested<std::runtime_error>(
                "Toolkit::Format", "buffer::heap::BufferMalloc", "BufferSystemV",
                "buffer overflow when resizing to " + std::to_string(size) + " bytes, " + hint);
        m_buffer = tmp;
    }
    m_size = size;
}

void BufferMalloc::Reset(const bool resetAbsolutePosition, const bool zeroInitialize)
{
    m_Position = 0;
    if (resetAbsolutePosition)
    {
        m_AbsolutePosition = 0;
    }
    if (zeroInitialize)
    {
        //        std::fill(m_Buffer.begin(), m_Buffer.end(), 0);
    }
    else
    {
        // just zero out the first and last 1kb
        const size_t bufsize = m_size;
        size_t s = (bufsize < 1024 ? bufsize : 1024);
        std::fill_n(m_buffer, s, 0);
        if (bufsize > 1024)
        {
            size_t pos = bufsize - 1024;
            if (pos < 1024)
            {
                pos = 1024;
            }
            s = bufsize - pos;
            std::fill_n(m_buffer + bufsize - pos, s, 0);
        }
    }
}

size_t BufferMalloc::GetAvailableSize() const { return m_size - m_Position; }

void BufferMalloc::Delete()
{
    if (m_size)
        free(m_buffer);
    m_size = 0;
}

size_t BufferMalloc::DebugGetSize() const { return m_size; };

size_t BufferMalloc::Size() const { return m_size; };

} // end namespace format
} // end namespace adios2
