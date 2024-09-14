/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileIME.cpp file I/O using IME Native library
 *
 *  Created on: Dec 1, 2019
 *      Author: Keichi Takahashi keichi@is.naist.jp
 */
#include "FileIME.h"

#include "adios2/helper/adiosFunctions.h"

#include <iostream>

#include <fcntl.h>        // O_* flags
#include <linux/limits.h> // PATH_MAX
#include <unistd.h>       // getcwd

extern "C" {
#include <im_client_native2.h>
}

/// \cond EXCLUDE_FROM_DOXYGEN
#include <ios> //std::ios_base::failure
/// \endcond

namespace adios2
{
namespace transport
{

std::atomic_uint FileIME::client_refcount(0);

FileIME::FileIME(helper::Comm const &comm) : Transport("File", "IME", comm)
{
    /** Initialize the IME client if there are no existing FileIME instances.
     * We need to check if the IME client has been already initialized because
     * initializing the IME multiple times leaves the client in a weired
     * state. */
    if (!client_refcount)
    {
        ime_client_native2_init();
    }

    client_refcount++;
}

FileIME::~FileIME()
{
    if (m_IsOpen)
    {
        /** Trigger a flush from IME to PFS. Note that fsync needs to be
         * called before bfs_sync. */
        if (m_SyncToPFS)
        {
            ime_client_native2_fsync(m_FileDescriptor);
            ime_client_native2_bfs_sync(m_FileDescriptor, true);
        }
        ime_client_native2_close(m_FileDescriptor);
    }

    client_refcount--;

    /** Finalize the IME client if there are no existing FileIME instances. */
    if (client_refcount)
    {
        ime_client_native2_finalize();
    }
}

/** Note that async mode is unsupported in FileIME. */
void FileIME::Open(const std::string &name, const Mode openMode, const bool async,
                   const bool directio)
{
    /** DEFAULT_IME_FILE_PREFIX is "ime://" */
    m_Name = DEFAULT_IME_FILE_PREFIX;

    /** Convert relative path to absolute path. The IME API only supports
     * absolute path. */
    if (name[0] != '/')
    {
        char tmp[PATH_MAX];
        m_Name += std::string(getcwd(tmp, PATH_MAX)) + "/";
    }
    m_Name += name;

    CheckName();
    m_OpenMode = openMode;
    switch (m_OpenMode)
    {
    case Mode::Write:
        ProfilerStart("open");
        m_FileDescriptor =
            ime_client_native2_open(m_Name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
        ProfilerStop("open");
        break;

    case Mode::Append:
        ProfilerStart("open");
        m_FileDescriptor = ime_client_native2_open(m_Name.c_str(), O_RDWR | O_CREAT, 0777);
        lseek(m_FileDescriptor, 0, SEEK_END);
        ProfilerStop("open");
        break;

    case Mode::Read:
        ProfilerStart("open");
        m_FileDescriptor = ime_client_native2_open(m_Name.c_str(), O_RDONLY, 0000);
        ProfilerStop("open");
        break;

    default:
        CheckFile("unknown open mode for file " + m_Name + ", in call to IME open");
    }

    CheckFile("couldn't open file " + m_Name +
              ", check permissions or path existence, in call to IME open");

    m_IsOpen = true;
}

void FileIME::SetParameters(const Params &params)
{
    auto param = params.find("synctopfs");
    if (param != params.end())
    {
        m_SyncToPFS = helper::StringTo<bool>(param->second, " in Parameter key=SyncToPFS");
    }
}

void FileIME::Write(const char *buffer, size_t size, size_t start)
{
    auto lf_Write = [&](const char *buffer, size_t size) {
        while (size > 0)
        {
            ProfilerStart("write");
            const auto writtenSize = ime_client_native2_write(m_FileDescriptor, buffer, size);
            ProfilerStop("write");

            if (writtenSize == -1)
            {
                if (errno == EINTR)
                {
                    continue;
                }

                helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME",
                                                      "Write", "couldn't write to file " + m_Name);
            }

            buffer += writtenSize;
            size -= writtenSize;
        }
    };

    if (start != MaxSizeT)
    {
        const auto newPosition = ime_client_native2_lseek(m_FileDescriptor, start, SEEK_SET);

        if (static_cast<size_t>(newPosition) != start)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Write",
                                                  "couldn't move to start position " +
                                                      std::to_string(start) + " in file " + m_Name);
        }
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

void FileIME::Read(char *buffer, size_t size, size_t start)
{
    auto lf_Read = [&](char *buffer, size_t size) {
        while (size > 0)
        {
            ProfilerStart("read");
            const auto readSize = ime_client_native2_read(m_FileDescriptor, buffer, size);
            ProfilerStop("read");

            if (readSize == -1)
            {
                if (errno == EINTR)
                {
                    continue;
                }

                helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Read",
                                                      "ERROR: couldn't read from file " + m_Name);
            }

            buffer += readSize;
            size -= readSize;
        }
    };

    if (start != MaxSizeT)
    {
        const auto newPosition = ime_client_native2_lseek(m_FileDescriptor, start, SEEK_SET);

        if (static_cast<size_t>(newPosition) != start)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Read",
                                                  "couldn't move to start position " +
                                                      std::to_string(start) + " in file " + m_Name +
                                                      ", errno " + std::to_string(errno));
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

size_t FileIME::GetSize()
{
    struct stat fileStat;
    if (fstat(m_FileDescriptor, &fileStat) == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "GetSize",
                                              "couldn't get size of file " + m_Name);
    }
    return static_cast<size_t>(fileStat.st_size);
}

void FileIME::Flush()
{
    if (m_SyncToPFS)
    {
        ime_client_native2_fsync(m_FileDescriptor);
        ime_client_native2_bfs_sync(m_FileDescriptor, true);
    }
}

void FileIME::Close()
{
    ProfilerStart("close");
    if (m_SyncToPFS)
    {
        ime_client_native2_fsync(m_FileDescriptor);
        ime_client_native2_bfs_sync(m_FileDescriptor, true);
    }
    const int status = ime_client_native2_close(m_FileDescriptor);
    ProfilerStop("close");

    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Close",
                                              "ERROR: couldn't close file " + m_Name);
    }

    m_IsOpen = false;
}

void FileIME::Delete()
{
    if (m_IsOpen)
    {
        Close();
    }
    ime_client_native2_unlink(m_Name.c_str());
}

void FileIME::CheckFile(const std::string hint) const
{
    if (m_FileDescriptor == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "CheckFile",
                                              hint);
    }
}

void FileIME::SeekToEnd()
{
    const int status = ime_client_native2_lseek(m_FileDescriptor, 0, SEEK_END);
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "SeekToEnd",
                                              "couldn't seek to the end of file " + m_Name);
    }
}

void FileIME::SeekToBegin()
{
    const int status = ime_client_native2_lseek(m_FileDescriptor, 0, SEEK_SET);
    if (status == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "SeekToBegin",
                                              "couldn't seek to the begin of file " + m_Name);
    }
}

void FileIME::Seek(const size_t start)
{
    if (start != MaxSizeT)
    {
        const int status = ime_client_native2_lseek(m_FileDescriptor, start, SEEK_SET);
        if (status == -1)
        {
            helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Seek",
                                                  "couldn't seek to offset " +
                                                      std::to_string(start) + " of file " + m_Name);
        }
    }
    else
    {
        SeekToEnd();
    }
}

void FileIME::Truncate(const size_t length)
{
    helper::Throw<std::ios_base::failure>("Toolkit", "transport::file::FileIME", "Truncate",
                                          "IME Truncate is not implemented yet");
}

void FileIME::MkDir(const std::string &fileName) {}

} // end namespace transport
} // end namespace adios2
