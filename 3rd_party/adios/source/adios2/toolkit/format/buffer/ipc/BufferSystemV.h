/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferSystemV.h
 *
 *  Created on: Jul 9, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BUFFER_IPC_BUFFERSYSTEMV_H_
#define ADIOS2_TOOLKIT_FORMAT_BUFFER_IPC_BUFFERSYSTEMV_H_

#include "adios2/toolkit/format/buffer/Buffer.h"

namespace adios2
{
namespace format
{

class BufferSystemV : public Buffer
{
public:
    BufferSystemV(const size_t fixedSize, const std::string &name, const unsigned int projectID,
                  const bool remove);

    ~BufferSystemV();

    char *Data() noexcept final;

    const char *Data() const noexcept final;

    void Reset(const bool resetAbsolutePosition, const bool zeroInitialize) final;

private:
    /** shared memory segment ID from shmget */
    int m_ShmID = -1;

    /** pointer to shared memory segment */
    char *m_Data = nullptr;

    /** false: make it persistent, true: remove with destructor */
    const bool m_Remove;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_FORMAT_BUFFER_IPC_BUFFERSYSTEMV_H_ */
