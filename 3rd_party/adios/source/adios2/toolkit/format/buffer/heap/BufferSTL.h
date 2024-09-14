/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferSTL.h
 *
 *  Created on: Sep 26, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERSTL_H_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERSTL_H_

#include "adios2/toolkit/format/buffer/Buffer.h"

#include "adios2/common/ADIOSMacros.h"

namespace adios2
{
namespace format
{

class BufferSTL : public Buffer
{
public:
    std::vector<char> m_Buffer;

    BufferSTL();
    ~BufferSTL() = default;

    char *Data() noexcept final;
    const char *Data() const noexcept final;

    void Resize(const size_t size, const std::string hint) final;

    void Reset(const bool resetAbsolutePosition, const bool zeroInitialize) final;

    size_t GetAvailableSize() const final;

    template <class T>
    size_t Align() const noexcept;

    void Delete();

    size_t DebugGetSize() const;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_HEAP_BUFFERSTL_H_ */
