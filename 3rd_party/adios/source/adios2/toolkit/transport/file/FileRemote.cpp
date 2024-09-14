/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 */
#include "FileRemote.h"
#include "adios2/core/ADIOS.h"
#include "adios2/helper/adiosLog.h"
#include "adios2/helper/adiosString.h"
#include "adios2/helper/adiosSystem.h"

#include <cstdio>  // remove
#include <cstring> // strerror
#include <errno.h> // errno
#include <fcntl.h> // open
#include <regex>
#include <sys/stat.h>  // open, fstat
#include <sys/types.h> // open
#include <unistd.h>    // write, close, ftruncate

namespace adios2
{
namespace transport
{

FileRemote::FileRemote(helper::Comm const &comm)
: Transport("File", "Remote", comm) /*, m_Impl(&m_ImplSingleton)*/
{
}

FileRemote::~FileRemote()
{
    if (m_IsOpen)
    {
        Close();
    }
}

void FileRemote::SetParameters(const Params &params)
{
    // Parameters are set from config parameters if present
    // Otherwise, they are set from environment if present
    // Otherwise, they remain at their default value

    helper::SetParameterValue("cache", params, m_CachePath);
    if (m_CachePath.empty())
    {
        if (const char *Env = std::getenv("AWS_CACHE"))
        {
            m_CachePath = std::string(Env);
        }
    }
}

void FileRemote::WaitForOpen() {}

void FileRemote::SetUpCache()
{
    if (!m_CachePath.empty())
    {
        if (helper::EndsWith(m_FileName, "md.idx") || helper::EndsWith(m_FileName, "md.0") ||
            helper::EndsWith(m_FileName, "mmd.0"))
        {
            m_CachingThisFile = true;
        }
    }

    if (m_CachingThisFile)
    {
    }
}

void FileRemote::Open(const std::string &name, const Mode openMode, const bool async,
                      const bool directio)
{
    m_Name = name;
    size_t pos = name.find(PathSeparator);
    if (pos == std::string::npos)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "transport::file::FileRemote", "Open",
                                             "invalid 'bucket/object' name " + name);
    }
    m_OpenMode = openMode;
    switch (m_OpenMode)
    {

    case Mode::Write:
    case Mode::Append:
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileRemote", "Open",
                                              "does not support writing yet " + m_Name);
        break;

    case Mode::Read: {
        ProfilerStart("open");
        m_Remote.OpenSimpleFile("localhost", RemoteCommon::ServerPort, m_Name);
        ProfilerStop("open");
        m_Size = m_Remote.m_Size;
        break;
    }
    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to Remote open");
    }
}

void FileRemote::OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                           const bool async, const bool directio)
{
    int token = 1;
    if (chainComm.Rank() > 0)
    {
        chainComm.Recv(&token, 1, chainComm.Rank() - 1, 0, "Chain token in FileRemote::OpenChain");
    }

    Open(name, openMode, async, directio);

    if (chainComm.Rank() < chainComm.Size() - 1)
    {
        chainComm.Isend(&token, 1, chainComm.Rank() + 1, 0,
                        "Sending Chain token in FileRemote::OpenChain");
    }
}

void FileRemote::Write(const char *buffer, size_t size, size_t start)
{
    helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileRemote", "Write",
                                          "does not support writing yet " + m_Name);
}

void FileRemote::Read(char *buffer, size_t size, size_t start)
{
    WaitForOpen();

    if (start != MaxSizeT)
    {
        if (start >= m_Size)
        {
            helper::Throw<std::ios_base::failure>(
                "Toolkit", "transport::file::FileRemote", "Read",
                "couldn't move to start position " + std::to_string(start) +
                    " beyond the size of " + m_Name + " which is " + std::to_string(m_Size));
        }
        m_SeekPos = start;
        errno = 0;
        m_Errno = errno;
    }

    if (m_SeekPos + size > m_Size)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileRemote", "Read",
                                              "can't read " + std::to_string(size) +
                                                  " bytes from position " +
                                                  std::to_string(m_SeekPos) + " from " + m_Name +
                                                  " whose size is " + std::to_string(m_Size));
    }

    m_Remote.Read(start, size, buffer);
    if (m_IsCached)
    {
    }
}

size_t FileRemote::GetSize()
{
    WaitForOpen();
    switch (m_OpenMode)
    {
    case Mode::Write:
    case Mode::Append:
        return 0;
    case Mode::Read:
        return m_Size;
    default:
        return 0;
    }
}

void FileRemote::Flush() {}

void FileRemote::Close()
{
    WaitForOpen();
    ProfilerStart("close");
    errno = 0;
    m_Errno = errno;
    if (m_IsCached)
    {
    }

    m_IsOpen = false;
    ProfilerStop("close");
}

void FileRemote::Delete()
{
    WaitForOpen();
    if (m_IsOpen)
    {
        Close();
    }
    std::remove(m_Name.c_str());
}

void FileRemote::CheckFile(const std::string hint) const {}

void FileRemote::SeekToEnd() { m_SeekPos = MaxSizeT; }

void FileRemote::SeekToBegin() { m_SeekPos = 0; }

void FileRemote::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        m_SeekPos = start;
    }
    else
    {
        SeekToEnd();
    }
}

void FileRemote::Truncate(const size_t length)
{
    helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileRemote", "Truncate",
                                          "does not support truncating " + m_Name);
}

void FileRemote::MkDir(const std::string &fileName) {}

} // end namespace transport
} // end namespace adios2
