/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileDescriptor.h wrapper of POSIX library functions for file I/O
 *
 *  Created on: Jul 7, 2023
 *      Author: Dmitry Ganyushin  ganyushin@gmail.com
 */
#include "FileHTTP.h"
#include <sys/socket.h>

#include <arpa/inet.h>
#include <cstring>
#include <netdb.h>
#include <unistd.h>
namespace adios2
{
namespace transport
{

FileHTTP::FileHTTP(helper::Comm const &comm) : Transport("File", "HTTP", comm) {}

FileHTTP::~FileHTTP()
{
    if (m_IsOpen)
    {
    }
}

void FileHTTP::WaitForOpen()
{
    if (m_IsOpening)
    {
        m_IsOpening = false;
        CheckFile("couldn't open file " + m_Name + ", in call to POSIX open");
        m_IsOpen = true;
    }
}

void FileHTTP::Open(const std::string &name, const Mode openMode, const bool async,
                    const bool directio)
{
    struct protoent *protoent;
    struct hostent *hostent;
    in_addr_t in_addr;

    m_Name = name;
    /* Build the socket. */
    protoent = getprotobyname("tcp");
    if (protoent == NULL)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Open",
                                              "cannot make getprotobyname");
    }
    m_p_proto = protoent->p_proto;

    /* get from parameters. Where the proxy should run*/

    /* Build the address. */
    hostent = gethostbyname(m_hostname.c_str());
    if (hostent == NULL)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Open",
                                              "error: gethostbyname " + m_hostname);
    }
    in_addr = inet_addr(inet_ntoa(*(struct in_addr *)*(hostent->h_addr_list)));
    if (in_addr == (in_addr_t)-1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Open",
                                              "error: inet_addr " +
                                                  std::string(*(hostent->h_addr_list)));
    }
    sockaddr_in.sin_addr.s_addr = in_addr;
    sockaddr_in.sin_family = AF_INET;
    sockaddr_in.sin_port = htons(m_server_port);

    return;
}

void FileHTTP::OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                         const bool async, const bool directio)
{
    return;
}

void FileHTTP::Write(const char *buffer, size_t size, size_t start) { return; }

#ifdef REALLY_WANT_WRITEV
void FilePOSIX::WriteV(const core::iovec *iov, const int iovcnt, size_t start)
{
    auto lf_Write = [&](const core::iovec *iov, const int iovcnt) {
        ProfilerStart("write");
        errno = 0;
        size_t nBytesExpected = 0;
        for (int i = 0; i < iovcnt; ++i)
        {
            nBytesExpected += iov[i].iov_len;
        }
        const iovec *v = reinterpret_cast<const iovec *>(iov);
        const auto ret = writev(m_FileDescriptor, v, iovcnt);
        m_Errno = errno;
        ProfilerStop("write");

        size_t written;
        if (ret == -1)
        {
            if (errno != EINTR)
            {
                helper::Throw<std::ios_base::failure>(
                    "Toolkit", "transport::file::FilePOSIX", "WriteV",
                    "couldn't write to file " + m_Name + " " + SysErrMsg());
            }
            written = 0;
        }
        else
        {
            written = static_cast<size_t>(ret);
        }

        ProfilerWriteBytes(written);

        if (written < nBytesExpected)
        {
            /* Fall back to write calls with individual buffers */
            // find where the writing has ended
            int c = 0;
            size_t n = 0;
            size_t pos = 0;
            while (n < written)
            {
                if (n + iov[c].iov_len <= written)
                {
                    n += iov[c].iov_len;
                    ++c;
                }
                else
                {
                    pos = written - n;
                    n = written;
                }
            }

            // write the rest one by one
            Write(static_cast<const char *>(iov[c].iov_base) + pos, iov[c].iov_len - pos);
            for (; c < iovcnt; ++c)
            {
                Write(static_cast<const char *>(iov[c].iov_base), iov[c].iov_len);
            }
        }
    };

    WaitForOpen();
    if (start != MaxSizeT)
    {
        errno = 0;
        const auto newPosition = lseek(m_FileDescriptor, start, SEEK_SET);
        m_Errno = errno;

        if (static_cast<size_t>(newPosition) != start)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "WriteV",
                                                  "couldn't move to start position " +
                                                      std::to_string(start) + " in file " + m_Name +
                                                      " " + SysErrMsg());
        }
    }

    int cntTotal = 0;
    while (cntTotal < iovcnt)
    {
        int cnt = iovcnt - cntTotal;
        if (cnt > 8)
        {
            cnt = 8;
        }
        lf_Write(iov + cntTotal, cnt);
        cntTotal += cnt;
    }
}
#endif

void FileHTTP::Read(char *buffer, size_t size, size_t start)
{
    /* not using BUFSIZ, the server might use another value for that */
    const size_t BUF_SIZE = 8192;
    enum CONSTEXPR
    {
        MAX_REQUEST_LEN = 1024
    };
    char request[MAX_REQUEST_LEN] = {'\0'};
    int request_len = snprintf(request, MAX_REQUEST_LEN, request_template.c_str(), m_Name.c_str(),
                               m_hostname.c_str(), start, start + size - 1);
    if (request_len >= MAX_REQUEST_LEN)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Read",
                                              "request length too long:  " +
                                                  std::to_string(request_len));
    }
    m_socketFileDescriptor = socket(AF_INET, SOCK_STREAM, m_p_proto);
    if (m_socketFileDescriptor == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Read",
                                              "cannot open socket");
    }
    /* Actually connect. */
    if (connect(m_socketFileDescriptor, (struct sockaddr *)&sockaddr_in, sizeof(sockaddr_in)) == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Read",
                                              "cannot connect");
    }
    /* Send HTTP request. */
    int nbytes_total = 0;
    while (nbytes_total < request_len)
    {
        int nbytes_last =
            write(m_socketFileDescriptor, request + nbytes_total, request_len - nbytes_total);
        if (nbytes_last == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Read",
                                                  "cannot send request");
        }
        nbytes_total += nbytes_last;
    }

    /* Read the response. */
    size_t bytes_recd = 0;
    while (bytes_recd < size)
    {
        nbytes_total = read(m_socketFileDescriptor, buffer + bytes_recd,
                            std::min(size - bytes_recd, BUF_SIZE));
        if (nbytes_total == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "Read",
                                                  "cannot get response");
        }
        bytes_recd += nbytes_total;
    }

    close(m_socketFileDescriptor);
    return;
}

size_t FileHTTP::GetSize()
{
    char request_template[] = "GET %s HTTP/1.1\r\nHost: %s\r\nContent-Length: bytes\r\n\r\n";
    enum CONSTEXPR
    {
        MAX_REQUEST_LEN = 1024,
        BUF_SIZE = 128
    };
    char request[MAX_REQUEST_LEN];
    char buffer[BUF_SIZE] = {'\0'};
    int request_len =
        snprintf(request, MAX_REQUEST_LEN, request_template, m_Name.c_str(), m_hostname.c_str());
    if (request_len >= MAX_REQUEST_LEN)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "GetSize",
                                              "request length too long:  " +
                                                  std::to_string(request_len));
    }
    m_socketFileDescriptor = socket(AF_INET, SOCK_STREAM, m_p_proto);
    if (m_socketFileDescriptor == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "GetSize",
                                              "cannot bind socket");
    }
    /* Actually connect. */
    if (connect(m_socketFileDescriptor, (struct sockaddr *)&sockaddr_in, sizeof(sockaddr_in)) == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "GetSize",
                                              "cannot connect");
    }
    /* Send HTTP request. */
    int nbytes_total = 0;
    while (nbytes_total < request_len)
    {
        int nbytes_last =
            write(m_socketFileDescriptor, request + nbytes_total, request_len - nbytes_total);
        if (nbytes_last == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "GetSize",
                                                  "sending request failed");
        }
        nbytes_total += nbytes_last;
    }

    /* Read the response. */
    while ((nbytes_total = read(m_socketFileDescriptor, buffer, BUF_SIZE)) > 0)
        ;
    if (nbytes_total == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileHTTP", "GetSize",
                                              "receiving response failed");
    }
    close(m_socketFileDescriptor);
    size_t result = atoi(buffer);

    return result;
}

void FileHTTP::Flush()
{
    /* Turn this off now because BP3/BP4 calls manager Flush and this syncing
     * slows down IO performance */
}

void FileHTTP::Close() { return; }

void FileHTTP::Delete() { return; }

void FileHTTP::CheckFile(const std::string hint) const
{
    if (m_socketFileDescriptor == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "CheckFile",
                                              hint + SysErrMsg());
    }
}

std::string FileHTTP::SysErrMsg() const
{
    return std::string(": errno = " + std::to_string(m_Errno) + ": " + strerror(m_Errno));
}

void FileHTTP::SeekToEnd() { return; }

void FileHTTP::SeekToBegin() { return; }

void FileHTTP::Seek(const size_t start) { return; }

void FileHTTP::Truncate(const size_t length) { return; }

void FileHTTP::MkDir(const std::string &fileName) { return; }

} // end namespace transport
} // end namespace adios2
