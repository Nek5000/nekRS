/*
 *
 * FilePOSIX.cpp file I/O using POSIX I/O library
 *
 */
#include "FilePOSIX.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"

#ifdef ADIOS2_HAVE_O_DIRECT
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif

#include <cstdio>      // remove
#include <cstring>     // strerror
#include <errno.h>     // errno
#include <fcntl.h>     // open
#include <sys/stat.h>  // open, fstat
#include <sys/types.h> // open
#include <thread>
#ifndef _MSC_VER
#include <unistd.h> // write, close, ftruncate
#define O_BINARY 0
#else
#include <io.h>
#define close _close
#define open _open
#define lseek(a, b, c) _lseek(a, (long)b, c)
#define write(a, b, c) _write(a, b, (unsigned int)c)
#define read(a, b, c) _read(a, b, (unsigned int)c)
#define ftruncate _chsize
#endif

/// \cond EXCLUDE_FROM_DOXYGEN
#include <ios> //std::ios_base::failure
/// \endcond

namespace adios2
{
namespace transport
{

FilePOSIX::FilePOSIX(helper::Comm const &comm) : Transport("File", "POSIX", comm) {}

FilePOSIX::~FilePOSIX()
{
    if (m_IsOpen)
    {
        close(m_FileDescriptor);
    }
}

void FilePOSIX::WaitForOpen()
{
    if (m_IsOpening)
    {
        if (m_OpenFuture.valid())
        {
            m_FileDescriptor = m_OpenFuture.get();
        }
        m_IsOpening = false;
        CheckFile("couldn't open file " + m_Name + ", in call to POSIX open");
        m_IsOpen = true;
    }
}

static int __GetOpenFlag(const int flag, const bool directio)
{
#ifdef ADIOS2_HAVE_O_DIRECT
    if (directio)
    {
        return flag | O_DIRECT;
    }
    else
#endif
    {
        return flag;
    }
}

void FilePOSIX::Open(const std::string &name, const Mode openMode, const bool async,
                     const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name, const bool directio) -> int {
        ProfilerStart("open");
        errno = 0;
        int flag = __GetOpenFlag(O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, directio);
        int FD = open(m_Name.c_str(), flag, 0666);
        m_Errno = errno;
        ProfilerStop("open");
        return FD;
    };

    m_Name = name;
    CheckName();
    m_DirectIO = directio;
    m_OpenMode = openMode;
    switch (m_OpenMode)
    {

    case Mode::Write:
        if (async)
        {
            m_IsOpening = true;
            m_OpenFuture = std::async(std::launch::async, lf_AsyncOpenWrite, name, directio);
        }
        else
        {
            ProfilerStart("open");
            errno = 0;
            m_FileDescriptor =
                open(m_Name.c_str(),
                     __GetOpenFlag(O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, directio), 0666);
            m_Errno = errno;
            ProfilerStop("open");
        }
        break;

    case Mode::Append:
        ProfilerStart("open");
        errno = 0;
        m_FileDescriptor =
            open(m_Name.c_str(), __GetOpenFlag(O_RDWR | O_CREAT | O_BINARY, directio), 0777);
        lseek(m_FileDescriptor, 0, SEEK_END);
        m_Errno = errno;
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        errno = 0;
        m_FileDescriptor = open(m_Name.c_str(), O_RDONLY | O_BINARY);
        m_Errno = errno;
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to POSIX open");
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name + ", in call to POSIX open");
        m_IsOpen = true;
    }
}

void FilePOSIX::OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                          const bool async, const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name, const bool directio) -> int {
        ProfilerStart("open");
        errno = 0;
        int flag = __GetOpenFlag(O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, directio);
        int FD = open(m_Name.c_str(), flag, 0666);
        m_Errno = errno;
        ProfilerStop("open");
        return FD;
    };

    int token = 1;
    m_Name = name;
    CheckName();

    if (chainComm.Rank() > 0)
    {
        chainComm.Recv(&token, 1, chainComm.Rank() - 1, 0, "Chain token in FilePOSIX::OpenChain");
    }

    m_DirectIO = directio;
    m_OpenMode = openMode;
    switch (m_OpenMode)
    {

    case Mode::Write:
        if (async && chainComm.Size() == 1)
        {
            // only when process is a single writer, can create the file
            // asynchronously, otherwise other processes are waiting on it
            m_IsOpening = true;
            m_OpenFuture = std::async(std::launch::async, lf_AsyncOpenWrite, name, directio);
        }
        else
        {
            ProfilerStart("open");
            errno = 0;
            if (chainComm.Rank() == 0)
            {
                m_FileDescriptor =
                    open(m_Name.c_str(),
                         __GetOpenFlag(O_WRONLY | O_CREAT | O_TRUNC | O_BINARY, directio), 0666);
            }
            else
            {
                m_FileDescriptor =
                    open(m_Name.c_str(), __GetOpenFlag(O_WRONLY | O_BINARY, directio), 0666);
                lseek(m_FileDescriptor, 0, SEEK_SET);
            }
            m_Errno = errno;
            ProfilerStop("open");
        }
        break;

    case Mode::Append:
        ProfilerStart("open");
        errno = 0;
        if (chainComm.Rank() == 0)
        {
            m_FileDescriptor =
                open(m_Name.c_str(), __GetOpenFlag(O_RDWR | O_CREAT | O_BINARY, directio), 0666);
        }
        else
        {
            m_FileDescriptor = open(m_Name.c_str(), __GetOpenFlag(O_RDWR | O_BINARY, directio));
        }
        lseek(m_FileDescriptor, 0, SEEK_END);
        m_Errno = errno;
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        errno = 0;
        m_FileDescriptor = open(m_Name.c_str(), O_RDONLY | O_BINARY);
        m_Errno = errno;
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to POSIX open");
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name + ", in call to POSIX open");
        m_IsOpen = true;
    }

    if (chainComm.Rank() < chainComm.Size() - 1)
    {
        chainComm.Isend(&token, 1, chainComm.Rank() + 1, 0,
                        "Sending Chain token in FilePOSIX::OpenChain");
    }
}

void FilePOSIX::Write(const char *buffer, size_t size, size_t start)
{
    auto lf_Write = [&](const char *buffer, size_t size) {
        while (size > 0)
        {
            ProfilerStart("write");
            ProfilerWriteBytes(size);
            errno = 0;
            const auto writtenSize = write(m_FileDescriptor, buffer, size);
            m_Errno = errno;
            ProfilerStop("write");

            if (writtenSize == -1)
            {
                if (errno == EINTR)
                {
                    continue;
                }

                helper::Throw<std::ios_base::failure>(
                    "Toolkit", "transport::file::FilePOSIX", "Write",
                    "couldn't write to file " + m_Name + " " + SysErrMsg());
            }

            buffer += writtenSize;
            size -= writtenSize;
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
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "Write",
                                                  "couldn't move to start position " +
                                                      std::to_string(start) + " in file " + m_Name +
                                                      " " + SysErrMsg());
        }
    }
    else
    {
        const auto pos = lseek(m_FileDescriptor, 0, SEEK_CUR);
        start = static_cast<size_t>(pos);
    }

    if (size > DefaultMaxFileBatchSize)
    {
        const size_t batches = size / DefaultMaxFileBatchSize;
        const size_t remainder = size % DefaultMaxFileBatchSize;

        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            lf_Write(&buffer[position], DefaultMaxFileBatchSize);
            position += DefaultMaxFileBatchSize;
        }
        lf_Write(&buffer[position], remainder);
    }
    else
    {
        lf_Write(buffer, size);
    }
}

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

void FilePOSIX::Read(char *buffer, size_t size, size_t start)
{
    auto lf_Read = [&](char *buffer, size_t size) {
        size_t backoff_ns = 20;
        while (size > 0)
        {
            ProfilerStart("read");
            errno = 0;
            const auto readSize = read(m_FileDescriptor, buffer, size);
            m_Errno = errno;
            ProfilerStop("read");

            if (readSize == -1)
            {
                if (errno == EINTR)
                {
                    continue;
                }

                helper::Throw<std::ios_base::failure>(
                    "Toolkit", "transport::file::FilePOSIX", "Read",
                    "couldn't read from file " + m_Name + " " + SysErrMsg());
            }
            else if (readSize == 0)
            {
                if (m_FailOnEOF)
                {
                    // we got an EOF on data that *should* be present,
                    // but maybe we've got filesystem consistency
                    // issues.  We'll wait, but if the backoff time
                    // reaches 30 seconds (nearly 45 seconds total
                    // wait time) and we still don't have data, treat
                    // this as a real failure and throw an exception.
                    std::this_thread::sleep_for(std::chrono::nanoseconds(backoff_ns));
                    backoff_ns *= 2;
                    if (std::chrono::nanoseconds(backoff_ns) > std::chrono::seconds(30))
                        helper::Throw<std::ios_base::failure>(
                            "Toolkit", "transport::file::FilePOSIX", "Read",
                            "Read past end of file on " + m_Name + " trying to read " +
                                std::to_string(size) + " bytes " + SysErrMsg());
                }
                else
                {
                    // read past EOF, but we're to wait for data.  Exponential backoff with a limit
                    // of .5 sec (500,000,000 nanosec)
                    std::this_thread::sleep_for(std::chrono::nanoseconds(backoff_ns));
                    constexpr size_t backoff_limit = 500 * 1000 * 1000;
                    backoff_ns *= 2;
                    if (backoff_ns > backoff_limit)
                        backoff_ns = backoff_limit;
                }
            }

            buffer += readSize;
            size -= readSize;
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
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "Read",
                                                  "couldn't move to start position " +
                                                      std::to_string(start) + " in file " + m_Name +
                                                      " " + SysErrMsg());
        }
    }

    if (size > DefaultMaxFileBatchSize)
    {
        const size_t batches = size / DefaultMaxFileBatchSize;
        const size_t remainder = size % DefaultMaxFileBatchSize;

        size_t position = 0;
        for (size_t b = 0; b < batches; ++b)
        {
            lf_Read(&buffer[position], DefaultMaxFileBatchSize);
            position += DefaultMaxFileBatchSize;
        }
        lf_Read(&buffer[position], remainder);
    }
    else
    {
        lf_Read(buffer, size);
    }
}

size_t FilePOSIX::GetSize()
{
    struct stat fileStat;
    WaitForOpen();
    errno = 0;
    if (fstat(m_FileDescriptor, &fileStat) == -1)
    {
        m_Errno = errno;
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "GetSize",
                                              "couldn't get size of file " + m_Name + SysErrMsg());
    }
    m_Errno = errno;
    return static_cast<size_t>(fileStat.st_size);
}

void FilePOSIX::Flush()
{
    /* Turn this off now because BP3/BP4 calls manager Flush and this syncing
     * slows down IO performance */
#if 0
#if (_POSIX_C_SOURCE >= 199309L || _XOPEN_SOURCE >= 500)
    fdatasync(m_FileDescriptor);
#else
    fsync(m_FileDescriptor)
#endif
#endif
}

void FilePOSIX::Close()
{
    WaitForOpen();
    ProfilerStart("close");
    errno = 0;
    const int status = close(m_FileDescriptor);
    m_Errno = errno;
    ProfilerStop("close");

    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "Close",
                                              "couldn't close file " + m_Name + " " + SysErrMsg());
    }

    m_IsOpen = false;
}

void FilePOSIX::Delete()
{
    WaitForOpen();
    if (m_IsOpen)
    {
        Close();
    }
    std::remove(m_Name.c_str());
}

void FilePOSIX::CheckFile(const std::string hint) const
{
    if (m_FileDescriptor == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "CheckFile",
                                              hint + SysErrMsg());
    }
}

std::string FilePOSIX::SysErrMsg() const
{
    return std::string(": errno = " + std::to_string(m_Errno) + ": " + strerror(m_Errno));
}

void FilePOSIX::SeekToEnd()
{
    WaitForOpen();
    errno = 0;
    const int status = lseek(m_FileDescriptor, 0, SEEK_END);
    m_Errno = 0;
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "SeekToEnd",
                                              "couldn't seek to the end of file " + m_Name + " " +
                                                  SysErrMsg());
    }
}

void FilePOSIX::SeekToBegin()
{
    WaitForOpen();
    errno = 0;
    const int status = lseek(m_FileDescriptor, 0, SEEK_SET);
    m_Errno = errno;
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "transport::file::FilePOSIX", "SeekToBegin",
            "couldn't seek to the begin of file " + m_Name + " " + SysErrMsg());
    }
}

void FilePOSIX::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        WaitForOpen();
        errno = 0;
        const int status = lseek(m_FileDescriptor, start, SEEK_SET);
        m_Errno = errno;
        if (status == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "Seek",
                                                  "couldn't seek to offset " +
                                                      std::to_string(start) + " of file " + m_Name +
                                                      " " + SysErrMsg());
        }
    }
    else
    {
        SeekToEnd();
    }
}

void FilePOSIX::Truncate(const size_t length)
{
    WaitForOpen();
    errno = 0;
    const int status = ftruncate(m_FileDescriptor, static_cast<off_t>(length));
    m_Errno = errno;
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FilePOSIX", "Truncate",
                                              "couldn't truncate to " + std::to_string(length) +
                                                  " bytes of file " + m_Name + " " + SysErrMsg());
    }
}

void FilePOSIX::MkDir(const std::string &fileName) {}

void FilePOSIX::SetParameters(const Params &params)
{
    // Parameters are set from config parameters if present
    // Otherwise, they are set from environment if present
    // Otherwise, they remain at their default value

    helper::GetParameter(params, "FailOnEOF", m_FailOnEOF);
}

} // end namespace transport
} // end namespace adios2
