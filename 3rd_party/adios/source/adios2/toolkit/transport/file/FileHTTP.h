/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileDescriptor.h wrapper of POSIX library functions for file I/O
 *
 *  Created on: Jul 7, 2023
 *      Author: Dmitry Ganyushin  ganyushin@gmail.com
 */
#ifndef ADIOS2_FILEHTTP_H
#define ADIOS2_FILEHTTP_H

#include "../Transport.h"
#include "adios2/common/ADIOSConfig.h"
#include <netinet/in.h>

namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File descriptor transport using the POSIX IO library */
class FileHTTP : public Transport
{

public:
    FileHTTP(helper::Comm const &comm);

    ~FileHTTP();

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

private:
    /** POSIX file handle returned by Open */
    int m_socketFileDescriptor = -1;
    int m_Errno = 0;
    bool m_IsOpening = false;
    /* if filename is very lomg, we can get lout from array boundaries */
    std::string request_template = "GET %s HTTP/1.1\r\nHost: %s\r\nRange: bytes=%d-%d\r\n\r\n";
    std::string m_hostname = "localhost";
    int m_server_port = 9999;
    struct sockaddr_in sockaddr_in;
    /* protocol number */
    int m_p_proto;

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

#endif // ADIOS2_FILEHTTP_H
