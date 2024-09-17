/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferSTL.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BufferSTL.h"
#include "BufferSTL.tcc"
#include "adios2/helper/adiosLog.h"
#include <cstdlib>
#include <cstring>

namespace adios2
{
namespace format
{

BufferSTL::BufferSTL() : Buffer("BufferSTL") {}

char *BufferSTL::Data() noexcept { return m_Buffer.data(); }

const char *BufferSTL::Data() const noexcept { return m_Buffer.data(); }

void BufferSTL::Resize(const size_t size, const std::string hint)
{
    try
    {
        // doing this will effectively replace the STL GNU default power of 2
        // reallocation.
        m_Buffer.reserve(size);
        // must initialize memory (secure)
        m_Buffer.resize(size, '\0');
    }
    catch (...)
    {
        // catch a bad_alloc
        helper::ThrowNested<std::runtime_error>(
            "Toolkit::Format", "buffer::heap::BufferSTL", "BufferSystemV",
            "buffer overflow when resizing to " + std::to_string(size) + " bytes, " + hint);
    }
}

void BufferSTL::Reset(const bool resetAbsolutePosition, const bool zeroInitialize)
{
    m_Position = 0;
    if (resetAbsolutePosition)
    {
        m_AbsolutePosition = 0;
    }
    if (zeroInitialize)
    {
        std::fill(m_Buffer.begin(), m_Buffer.end(), 0);
    }
    else
    {
        // just zero out the first and last 1kb
        const size_t bufsize = m_Buffer.size();
        size_t s = (bufsize < 1024 ? bufsize : 1024);
        std::fill_n(m_Buffer.begin(), s, 0);
        if (bufsize > 1024)
        {
            size_t pos = bufsize - 1024;
            if (pos < 1024)
            {
                pos = 1024;
            }
            s = bufsize - pos;
            std::fill_n(next(m_Buffer.begin(), pos), s, 0);
        }
    }
}

size_t BufferSTL::GetAvailableSize() const { return m_Buffer.size() - m_Position; }

#define declare_template_instantiation(T) template size_t BufferSTL::Align<T>() const noexcept;

ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

void BufferSTL::Delete() { std::vector<char>().swap(m_Buffer); }

size_t BufferSTL::DebugGetSize() const { return m_Buffer.size(); };

} // end namespace format
} // end namespace adios2
