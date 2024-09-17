/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Transport.h defines Transport abstract base class
 *
 *  Created on: Oct 6, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_TRANSPORT_H_
#define ADIOS2_TOOLKIT_TRANSPORT_TRANSPORT_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/profiling/iochrono/IOChrono.h"

namespace adios2
{
/** Abstract Base class for transports */
class Transport
{

public:
    const std::string m_Type;          ///< transport type from derived class
    const std::string m_Library;       ///< library implementation (POSIX, Mdtm, etc.)
    std::string m_Name;                ///< from Open, unique identifier (e.g. filename)
    Mode m_OpenMode = Mode::Undefined; ///< at Open from ADIOSTypes.h
    bool m_IsOpen = false;             ///< true: open for communication, false: unreachable
    helper::Comm const &m_Comm;        ///< current multi-process communicator
    profiling::IOChrono m_Profiler;    ///< profiles Open, Write/Read, Close

    struct Status
    {
        size_t Bytes;
        bool Running;
        bool Successful;
        // TODO add more thing...time?
    };

    /**
     * Base constructor that all derived classes pass
     * @param type from derived class
     * @param comm passed to m_Comm
     */
    Transport(const std::string type, const std::string library, helper::Comm const &comm);

    virtual ~Transport() = default;

    void InitProfiler(const Mode openMode, const TimeUnit timeUnit);

    /**
     * Opens transport, possibly asynchronously, required before SetBuffer,
     * Write, Read, Flush, Close
     * @param name
     * @param openMode
     * @param async
     */
    virtual void Open(const std::string &name, const Mode openMode, const bool async = false,
                      const bool directio = false) = 0;

    /**
     * Opens transport, possibly asynchronously, in a chain to avoid
     * DOS attack to the file system.
     * Required before SetBuffer, Write, Read, Flush, Close
     * @param name
     * @param openMode
     * @param chainComm
     * @param async
     */
    virtual void OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                           const bool async = false, const bool directio = false);

    /**
     * If OS buffered (FILE* or fstream), sets the buffer size
     * @param buffer address for OS buffering
     * @param size of OS buffer
     */
    virtual void SetBuffer(char *buffer, size_t size);

    /**
     * Set transport parameters
     * @param parameters transport parameters to set
     */
    virtual void SetParameters(const Params &parameters);

    /**
     * Writes to transport. Note that size is non-const due to the nature of
     * underlying transport libraries
     * @param buffer raw data to be written
     * @param size number of bytes to be written
     * @param start starting position for writing (to allow rewind), if not
     * passed then start at current stream position
     */
    virtual void Write(const char *buffer, size_t size, size_t start = MaxSizeT) = 0;

    /**
     * Writes to transport, writev version. Note that size is non-const due to
     * the nature of underlying transport libraries
     * @param iovec array pointer
     * @param iovcnt number of entries
     * @param start starting position for writing (to allow rewind), if not
     * passed then start at current stream position
     */
    virtual void WriteV(const core::iovec *iov, const int iovcnt, size_t start = MaxSizeT);

    /**
     * Reads from transport "size" bytes from a certain position. Note that size
     * and position and non-const due to the nature of underlying transport
     * libraries
     * @param buffer raw data pointer to put the read bytes (must be
     * preallocated)
     * @param size number of bytes to be read
     * @param start starting position for read, if not passed then start at
     * current stream position
     */
    virtual void Read(char *buffer, size_t size, size_t start = MaxSizeT) = 0;

    /**
     * Returns the size of current data in transport
     * @return size as size_t
     */
    virtual size_t GetSize();

    /** flushes current contents to physical medium without closing */
    virtual void Flush();

    /** closes current file, after this file becomes unreachable */
    virtual void Close() = 0;

    /** deletes current file, after this file becomes unreachable */
    virtual void Delete() = 0;

    virtual void SeekToEnd() = 0;

    virtual void SeekToBegin() = 0;

    virtual void Seek(const size_t start = MaxSizeT) = 0;

    virtual void Truncate(const size_t length) = 0;

    virtual void MkDir(const std::string &fileName) = 0;

protected:
    void ProfilerWriteBytes(size_t bytes) noexcept;
    void ProfilerStart(const std::string process) noexcept;

    void ProfilerStop(const std::string process) noexcept;

    virtual void CheckName() const;
};

} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_TRANSPORT_TRANSPORT_H_ */
