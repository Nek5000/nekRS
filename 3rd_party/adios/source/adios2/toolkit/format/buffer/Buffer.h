/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Buffer.h : abstract class for all buffer memory types
 *
 *  Created on: Jul 9, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_BUFFER_H_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_BUFFER_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{
namespace format
{

class Buffer
{
public:
    const std::string m_Type;

    /** if 0: buffer can be extended, if >0: buffer has a fixed size */
    const size_t m_FixedSize = 0;

    size_t m_Position = 0;
    size_t m_AbsolutePosition = 0;

    Buffer(const std::string type, const size_t fixedSize = 0);

    virtual ~Buffer() = default;

    virtual char *Data() noexcept;
    virtual const char *Data() const noexcept;

    virtual void Resize(const size_t size, const std::string hint);

    /**
     * Reset the buffer m_Position to zero
     * @param resetAbsolutePosition true: reset m_AbsolutePosition = 0
     * @param zeroInitialize populate current buffer contents with '\0'
     */
    virtual void Reset(const bool resetAbsolutePosition, const bool zeroInitialize);

    virtual size_t GetAvailableSize() const;

    /** Deletes buffer memory manually */
    virtual void Delete();
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_BUFFER_H_ */
