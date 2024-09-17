/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferV.cpp
 *
 */

#include "BufferV.h"
#include <assert.h>
#include <stddef.h> // max_align_t
#include <string.h>

namespace adios2
{
namespace format
{

BufferV::BufferV(const std::string type, const bool AlwaysCopy, const size_t MemAlign,
                 const size_t MemBlockSize)
: m_Type(type), m_MemAlign(MemAlign), m_MemBlockSize(MemBlockSize), m_AlwaysCopy(AlwaysCopy)
{
}

BufferV::~BufferV() {}

void BufferV::Reset()
{
    CurOffset = 0;
    m_internalPos = 0;
    DataV.clear();
}

uint64_t BufferV::Size() noexcept { return CurOffset; }

void BufferV::AlignBuffer(const size_t align)
{
    size_t badAlign = CurOffset % align;
    if (badAlign)
    {
        size_t addAlign = align - badAlign;
        if (addAlign > zero.size())
        {
            zero.resize(addAlign);
        }
        AddToVec(addAlign, zero.data(), 1, true);
    }
}

} // end namespace format
} // end namespace adios2
