/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Buffer.cpp
 *
 *  Created on: Jul 9, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Buffer.h"
#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace format
{

Buffer::Buffer(const std::string type, const size_t fixedSize)
: m_Type(type), m_FixedSize(fixedSize)
{
}

void Buffer::Resize(const size_t size, const std::string hint)
{
    helper::Throw<std::invalid_argument>("Toolkit", "format::Buffer", "Resize",
                                         "buffer memory of type " + m_Type + " can't call Resize " +
                                             hint);
}

void Buffer::Reset(const bool resetAbsolutePosition, const bool zeroInitialize)
{
    helper::Throw<std::invalid_argument>("Toolkit", "format::Buffer", "Reset",
                                         "buffer memory of type " + m_Type + " can't call Reset");
}

char *Buffer::Data() noexcept { return nullptr; }

const char *Buffer::Data() const noexcept { return nullptr; }

size_t Buffer::GetAvailableSize() const
{
    if (m_FixedSize > 0 && m_FixedSize >= m_Position)
    {
        return m_FixedSize - m_Position;
    }
    return 0;
}

void Buffer::Delete()
{
    helper::Throw<std::invalid_argument>("Toolkit", "format::Buffer", "Delete",
                                         "buffer memory of type " + m_Type + " can't call Delete");
}

} // end namespace format
} // end namespace adios2
