/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileDescriptor.h wrapper of POSIX library functions for file I/O
 *
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_FILE_FILEDESCRIPTOR_H_
#define ADIOS2_TOOLKIT_TRANSPORT_FILE_FILEDESCRIPTOR_H_

#include <future> //std::async, std::future

#include "adios2/common/ADIOSConfig.h"
#include "adios2/toolkit/transport/Transport.h"

namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File descriptor transport using the POSIX IO library */
class FilePOSIX : public Transport
{

public:
    FilePOSIX(helper::Comm const &comm);

    ~FilePOSIX();

    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) final;

    void OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                   const bool async = false, const bool directio = false) final;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) final;

#ifdef REALLY_WANT_WRITEV
    /* Actual writev() function, inactive for now */
    void WriteV(const core::iovec *iov, const int iovcnt, size_t start = MaxSizeT) final;
#endif

    void Read(char *buffer, size_t size, size_t start = MaxSizeT) final;

    size_t GetSize() final;

    /** Does nothing, each write is supposed to flush */
    void Flush() final;

    void Close() final;

    void Delete() final;

    void SeekToEnd() final;

    void SeekToBegin() final;

    void Seek(const size_t start = MaxSizeT) final;

    void Truncate(const size_t length) final;

    void MkDir(const std::string &fileName) final;

    void SetParameters(const Params &params) final;

private:
    /** POSIX file handle returned by Open */
    int m_FileDescriptor = -1;
    int m_Errno = 0;
    bool m_FailOnEOF = false; // default to false for historic reasons
    bool m_IsOpening = false;
    std::future<int> m_OpenFuture;
    bool m_DirectIO = false;

    /**
     * Check if m_FileDescriptor is -1 after an operation
     * @param hint exception message
     */
    void CheckFile(const std::string hint) const;
    void WaitForOpen();
    std::string SysErrMsg() const;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TRANSPORT_FILE_FILEDESCRIPTOR_H_ */
