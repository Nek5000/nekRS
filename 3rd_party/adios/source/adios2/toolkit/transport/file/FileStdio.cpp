/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileStdio.cpp
 *
 *  Created on: Jan 6, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "FileStdio.h"
#include "adios2/helper/adiosLog.h"

/// \cond EXCLUDE_FROM_DOXYGEN
#include <cerrno>
#include <cstring>
#include <ios> //std::ios_base::failure
/// \endcond

// removes fopen warning on Windows
#ifdef _WIN32
#pragma warning(disable : 4996) // fopen
#endif

namespace adios2
{
namespace transport
{

FileStdio::FileStdio(helper::Comm const &comm) : Transport("File", "stdio", comm) {}

FileStdio::~FileStdio()
{
    if (m_IsOpen)
    {
        std::fclose(m_File);
    }
}

void FileStdio::WaitForOpen()
{
    if (m_IsOpening)
    {
        if (m_OpenFuture.valid())
        {
            m_File = m_OpenFuture.get();
        }
        m_IsOpening = false;
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to POSIX open");
        m_IsOpen = true;
        if (m_DelayedBufferSet)
        {
            SetBuffer(m_DelayedBuffer, m_DelayedBufferSize);
        }
    }
}

void FileStdio::Open(const std::string &name, const Mode openMode, const bool async,
                     const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name) -> FILE * {
        errno = 0;
        return std::fopen(name.c_str(), "wb");
    };
    m_Name = name;
    CheckName();
    m_OpenMode = openMode;

    switch (m_OpenMode)
    {
    case Mode::Write:
        if (async)
        {
            m_IsOpening = true;
            m_OpenFuture = std::async(std::launch::async, lf_AsyncOpenWrite, name);
        }
        else
        {
            errno = 0;
            m_File = std::fopen(name.c_str(), "wb");
        }
        break;
    case Mode::Append:
        errno = 0;
        m_File = std::fopen(name.c_str(), "rwb");
        std::fseek(m_File, 0, SEEK_END);
        break;
    case Mode::Read:
        errno = 0;
        m_File = std::fopen(name.c_str(), "rb");
        break;
    default:
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Open",
                                              "unknown open mode for file " + m_Name);
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to stdio open");
        m_IsOpen = true;
    }
}

void FileStdio::OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                          const bool async, const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name) -> FILE * {
        errno = 0;
        return std::fopen(name.c_str(), "wb");
    };

    int token = 1;
    m_Name = name;
    CheckName();

    if (chainComm.Rank() > 0)
    {
        chainComm.Recv(&token, 1, chainComm.Rank() - 1, 0, "Chain token in FileStdio::OpenChain");
    }

    m_OpenMode = openMode;
    switch (m_OpenMode)
    {
    case Mode::Write:
        if (async)
        {
            m_IsOpening = true;
            m_OpenFuture = std::async(std::launch::async, lf_AsyncOpenWrite, name);
        }
        else
        {
            errno = 0;
            m_File = std::fopen(name.c_str(), "wb");
        }
        break;
    case Mode::Append:
        errno = 0;
        m_File = std::fopen(name.c_str(), "rwb");
        std::fseek(m_File, 0, SEEK_END);
        break;
    case Mode::Read:
        errno = 0;
        m_File = std::fopen(name.c_str(), "rb");
        break;
    default:
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Open",
                                              "unknown open mode for file " + m_Name);
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to stdio open");
        m_IsOpen = true;
    }

    if (chainComm.Rank() < chainComm.Size() - 1)
    {
        chainComm.Isend(&token, 1, chainComm.Rank() + 1, 0,
                        "Sending Chain token in FileStdio::OpenChain");
    }
}

void FileStdio::SetBuffer(char *buffer, size_t size)
{
    if (!m_File)
    {
        m_DelayedBufferSet = true;
        m_DelayedBuffer = buffer;
        m_DelayedBufferSize = size;
        return;
    }
    m_DelayedBufferSet = false;
    m_DelayedBuffer = nullptr;
    m_DelayedBufferSize = 0;

    int status;
    if (buffer)
    {
        status = std::setvbuf(m_File, buffer, _IOFBF, size);
    }
    else
    {
        if (size != 0)
        {
            helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileStdio",
                                                 "SetBuffer",
                                                 "buffer size must be 0 when using a NULL buffer");
        }
        status = std::setvbuf(m_File, NULL, _IONBF, 0);
    }

    if (status)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "SetBuffer",
                                              "could not set FILE* buffer in file " + m_Name +
                                                  ", in call to stdio setvbuf");
    }
}

void FileStdio::Write(const char *buffer, size_t size, size_t start)
{
    auto lf_Write = [&](const char *buffer, size_t size) {
        ProfilerStart("write");
        const auto writtenSize = std::fwrite(buffer, sizeof(char), size, m_File);
        ProfilerStop("write");

        CheckFile("couldn't write to file " + m_Name + ", in call to stdio fwrite");

        if (writtenSize != size)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Write",
                                                  "written size + " + std::to_string(writtenSize) +
                                                      " is not equal to intended size " +
                                                      std::to_string(size) + " in file " + m_Name +
                                                      ", in call to stdio fwrite");
        }
    };

    WaitForOpen();
    if (start != MaxSizeT)
    {
        const auto status = std::fseek(m_File, static_cast<long int>(start), SEEK_SET);
        if (status != 0)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Write",
                                                  "couldn't move position of " + m_Name +
                                                      " file, in call to FileStdio Write fseek");
        }

        CheckFile("couldn't move to start position " + std::to_string(start) + " in file " +
                  m_Name + ", in call to stdio fseek at write ");
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

void FileStdio::Read(char *buffer, size_t size, size_t start)
{
    auto lf_Read = [&](char *buffer, size_t size) {
        ProfilerStart("read");
        const auto readSize = std::fread(buffer, sizeof(char), size, m_File);
        ProfilerStop("read");

        CheckFile("couldn't read to file " + m_Name + ", in call to stdio fread");

        if (readSize != size)
        {
            helper::Throw<std::ios_base::failure>(
                "Toolkit", "transport::file::FileStdio", "Read",
                "read size of " + std::to_string(readSize) + " is not equal to intended size " +
                    std::to_string(size) + " in file " + m_Name + ", in call to stdio fread");
        }
    };

    WaitForOpen();
    if (start != MaxSizeT)
    {
        const auto status = std::fseek(m_File, static_cast<long int>(start), SEEK_SET);
        CheckFile("couldn't move to start position " + std::to_string(start) + " in file " +
                  m_Name + ", in call to stdio fseek for read, result=" + std::to_string(status));
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

size_t FileStdio::GetSize()
{
    WaitForOpen();
    const auto currentPosition = ftell(m_File);
    if (currentPosition == -1L)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "GetSize",
                                              "couldn't get current position of " + m_Name +
                                                  " file, in call to FileStdio GetSize ftell");
    }

    fseek(m_File, 0, SEEK_END);
    const auto size = ftell(m_File);
    if (size == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "GetSize",
                                              "couldn't get size of " + m_Name +
                                                  " file, in call to FileStdio GetSize ftell");
    }
    fseek(m_File, currentPosition, SEEK_SET);
    return static_cast<size_t>(size);
}

void FileStdio::Flush()
{
    WaitForOpen();
    ProfilerStart("write");
    const int status = std::fflush(m_File);
    ProfilerStop("write");

    if (status == EOF)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Flush",
                                              "couldn't flush file " + m_Name +
                                                  ", in call to stdio fflush");
    }
}

void FileStdio::Close()
{
    WaitForOpen();
    ProfilerStart("close");
    const int status = std::fclose(m_File);
    ProfilerStop("close");

    if (status == EOF)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Close",
                                              "couldn't close file " + m_Name +
                                                  ", in call to stdio fclose");
    }

    m_IsOpen = false;
}

void FileStdio::Delete()
{
    WaitForOpen();
    if (m_IsOpen)
    {
        Close();
    };
    std::remove(m_Name.c_str());
}

void FileStdio::CheckFile(const std::string hint) const
{
    if (!m_File)
    {
        std::string errmsg;
        if (errno)
        {
            errmsg = std::strerror(errno);
        }
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "CheckFile",
                                              "ERROR: " + hint + ":" + errmsg);
    }
    else if (std::ferror(m_File))
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "CheckFile",
                                              "ERROR: " + hint);
    }
}

void FileStdio::SeekToEnd()
{
    WaitForOpen();
    const auto status = std::fseek(m_File, 0, SEEK_END);
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "SeekToEnd",
                                              "couldn't seek to the end of file " + m_Name);
    }
}

void FileStdio::SeekToBegin()
{
    WaitForOpen();
    const auto status = std::fseek(m_File, 0, SEEK_SET);
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio",
                                              "SeekToBegin",
                                              "couldn't seek to the begin of file " + m_Name);
    }
}

void FileStdio::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        WaitForOpen();
        const auto status = std::fseek(m_File, 0, SEEK_SET);
        if (status == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Seek",
                                                  "couldn't seek to offset " +
                                                      std::to_string(start) + " of file " + m_Name);
        }
    }
    else
    {
        SeekToEnd();
    }
}

#ifdef _WIN32
void FileStdio::Truncate(const size_t length)
{
    helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Truncate",
                                          "This is not supported on Windows");
}
#else
#include <unistd.h> // ftruncate
void FileStdio::Truncate(const size_t length)
{

    WaitForOpen();
    int fd = fileno(m_File);
    const auto status = ftruncate(fd, length);
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileStdio", "Truncate",
                                              "couldn't truncate to " + std::to_string(length) +
                                                  " of file " + m_Name);
    }
}
#endif

void FileStdio::MkDir(const std::string &fileName) {}

} // end namespace transport
} // end namespace adios2
