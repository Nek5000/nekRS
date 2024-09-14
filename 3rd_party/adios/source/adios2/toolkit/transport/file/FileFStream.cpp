/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileStream.cpp
 *
 *  Created on: Oct 24, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "FileFStream.h"
#include "adios2/helper/adiosLog.h"
#include <cstdio> // remove

/// \cond EXCLUDE_FROM_DOXYGEN
#include <ios> // std::ios_base::failure
/// \endcond

#if __cplusplus >= 201703L
#include <filesystem>
#endif

namespace adios2
{
namespace transport
{

FileFStream::FileFStream(helper::Comm const &comm) : Transport("File", "fstream", comm) {}

void FileFStream::WaitForOpen()
{
    if (m_IsOpening)
    {
        if (m_OpenFuture.valid())
        {
            m_OpenFuture.get();
        }
        m_IsOpening = false;
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to POSIX open");
        m_IsOpen = true;
    }
}

void FileFStream::Open(const std::string &name, const Mode openMode, const bool async,
                       const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name) -> void {
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::out | std::fstream::binary | std::fstream::trunc);
        ProfilerStop("open");
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
            ProfilerStart("open");
            m_FileStream.open(name, std::fstream::out | std::fstream::binary | std::fstream::trunc);
            ProfilerStop("open");
        }
        break;

    case Mode::Append:
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::in | std::fstream::out | std::fstream::binary);
        m_FileStream.seekp(0, std::ios_base::end);
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::in | std::fstream::binary);
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to stream open");
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to fstream open");
        m_IsOpen = true;
    }
}

void FileFStream::OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                            const bool async, const bool directio)
{
    auto lf_AsyncOpenWrite = [&](const std::string &name) -> void {
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::out | std::fstream::binary | std::fstream::trunc);
        ProfilerStop("open");
    };

    int token = 1;
    m_Name = name;
    CheckName();

    if (chainComm.Rank() > 0)
    {
        chainComm.Recv(&token, 1, chainComm.Rank() - 1, 0, "Chain token in FileFStream::OpenChain");
    }

    m_OpenMode = openMode;
    switch (m_OpenMode)
    {
    case Mode::Write:
        if (async && chainComm.Size() == 1)
        {
            m_IsOpening = true;
            m_OpenFuture = std::async(std::launch::async, lf_AsyncOpenWrite, name);
        }
        else
        {
            ProfilerStart("open");
            if (chainComm.Rank() == 0)
            {
                m_FileStream.open(name,
                                  std::fstream::out | std::fstream::binary | std::fstream::trunc);
            }
            else
            {
                m_FileStream.open(name, std::fstream::out | std::fstream::binary);
            }
            ProfilerStop("open");
        }
        break;

    case Mode::Append:
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::in | std::fstream::out | std::fstream::binary);
        m_FileStream.seekp(0, std::ios_base::end);
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        m_FileStream.open(name, std::fstream::in | std::fstream::binary);
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to stream open");
    }

    if (!m_IsOpening)
    {
        CheckFile("couldn't open file " + m_Name +
                  ", check permissions or path existence, in call to fstream open");
        m_IsOpen = true;
    }

    if (chainComm.Rank() < chainComm.Size() - 1)
    {
        chainComm.Isend(&token, 1, chainComm.Rank() + 1, 0,
                        "Sending Chain token in FileFStream::OpenChain");
    }
}

void FileFStream::SetBuffer(char *buffer, size_t size)
{
    if (!buffer && size != 0)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileFStream", "SetBuffer",
                                             "buffer size must be 0 when using a NULL buffer");
    }
    m_FileStream.rdbuf()->pubsetbuf(buffer, size);
    CheckFile("couldn't set buffer in file " + m_Name + ", in call to fstream rdbuf()->pubsetbuf");
}

void FileFStream::Write(const char *buffer, size_t size, size_t start)
{
    auto lf_Write = [&](const char *buffer, size_t size) {
        ProfilerStart("write");
        m_FileStream.write(buffer, static_cast<std::streamsize>(size));
        ProfilerStop("write");
        CheckFile("couldn't write from file " + m_Name + ", in call to fstream write");
    };

    WaitForOpen();
    if (start != MaxSizeT)
    {
        m_FileStream.seekp(start);
        CheckFile("couldn't move to start position " + std::to_string(start) + " in file " +
                  m_Name + ", in call to fstream seekp");
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

void FileFStream::Read(char *buffer, size_t size, size_t start)
{
    auto lf_Read = [&](char *buffer, size_t size) {
        ProfilerStart("read");
        m_FileStream.read(buffer, static_cast<std::streamsize>(size));
        ProfilerStop("read");
        CheckFile("couldn't read from file " + m_Name + ", in call to fstream read");
    };

    WaitForOpen();
    if (start != MaxSizeT)
    {
        m_FileStream.seekg(start);
        CheckFile("couldn't move to start position " + std::to_string(start) + " in file " +
                  m_Name + ", in call to fstream seekg");
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

size_t FileFStream::GetSize()
{
    WaitForOpen();
    const auto currentPosition = m_FileStream.tellg();
    m_FileStream.seekg(0, std::ios_base::end);
    const std::streampos size = m_FileStream.tellg();
    if (static_cast<int>(size) == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileFStream", "GetSize",
                                              "couldn't get size of " + m_Name + " file");
    }
    m_FileStream.seekg(currentPosition);
    return static_cast<size_t>(size);
}

void FileFStream::Flush()
{
    WaitForOpen();
    ProfilerStart("write");
    m_FileStream.flush();
    ProfilerStart("write");
    CheckFile("couldn't flush to file " + m_Name + ", in call to fstream flush");
}

void FileFStream::Close()
{
    WaitForOpen();
    ProfilerStart("close");
    m_FileStream.close();
    ProfilerStop("close");

    CheckFile("couldn't close file " + m_Name + ", in call to fstream close");
    m_IsOpen = false;
}

void FileFStream::Delete()
{
    WaitForOpen();
    if (m_IsOpen)
    {
        Close();
    }
    std::remove(m_Name.c_str());
}

void FileFStream::CheckFile(const std::string hint) const
{
    if (!m_FileStream)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileFStream",
                                              "CheckFile", hint);
    }
}

void FileFStream::SeekToEnd()
{
    WaitForOpen();
    m_FileStream.seekp(0, std::ios_base::end);
    CheckFile("couldn't move to the end of file " + m_Name + ", in call to fstream seekp");
}

void FileFStream::SeekToBegin()
{
    WaitForOpen();
    m_FileStream.seekp(0, std::ios_base::beg);
    CheckFile("couldn't move to the beginning of file " + m_Name + ", in call to fstream seekp");
}

void FileFStream::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        WaitForOpen();
        m_FileStream.seekp(start, std::ios_base::beg);
        CheckFile("couldn't move to offset " + std::to_string(start) + " of file " + m_Name +
                  ", in call to fstream seekp");
    }
    else
    {
        SeekToEnd();
    }
}

void FileFStream::Truncate(const size_t length)
{
#if __cplusplus >= 201703L
    // C++17 specific stuff here
    WaitForOpen();
    std::filesystem::path p(m_Name);
    std::filesystem::resize_file(p, static_cast<std::uintmax_t>(length));
    // TODO: variable start has not been defined...
    // CheckFile("couldn't move to offset " + std::to_string(start) + " of file
    // " + m_Name + ", in call to fstream seekp");
#else
    // Trunation is not supported in a portable manner pre C++17
#endif
}

void FileFStream::MkDir(const std::string &fileName) {}

} // end namespace transport
} // end namespace adios2
