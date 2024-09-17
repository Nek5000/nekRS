/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * ShmSystemV.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "ShmSystemV.h"
#include "adios2/helper/adiosLog.h"

#include <cstring> //std::memcpy

#include <sys/ipc.h>   //ftok
#include <sys/shm.h>   //shmget, shmmat
#include <sys/types.h> //key_t

namespace adios2
{
namespace transport
{

ShmSystemV::ShmSystemV(const unsigned int projectID, const size_t size, helper::Comm const &comm,
                       const bool removeAtClose)
: Transport("Shm", "SystemV", comm), m_ProjectID(projectID), m_Size(size),
  m_RemoveAtClose(removeAtClose)
{
    if (projectID == 0)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "transport::shm::ShmSystemV", "ShmSystemV",
                                             "projectID can't be zero, in shared memory segment");
    }
}

ShmSystemV::~ShmSystemV() // this might not be correct
{
    if (m_IsOpen)
    {
        shmdt(m_Buffer);

        if (m_RemoveAtClose)
        {
            shmctl(m_ShmID, IPC_RMID, NULL);
        }
    }
}

void ShmSystemV::Open(const std::string &name, const Mode openMode, const bool async,
                      const bool directio)
{
    m_Name = name;
    CheckName();
    m_OpenMode = openMode;

    // not using const
    key_t key = ftok(m_Name.c_str(), static_cast<int>(m_ProjectID));

    switch (m_OpenMode)
    {
    case (Mode::Write):
        ProfilerStart("open");
        m_ShmID = shmget(key, m_Size, IPC_CREAT | 0666);
        ProfilerStop("open");
        break;

    case (Mode::Append):
        ProfilerStart("open");
        m_ShmID = shmget(key, m_Size, 0);
        ProfilerStop("open");
        break;

    case (Mode::Read):
        ProfilerStart("open");
        m_ShmID = shmget(key, m_Size, 0);
        ProfilerStop("open");
        break;

    default:
        helper::Throw<std::invalid_argument>("Toolkit", "transport::shm::ShmSystemV", "Open",
                                             "unknown open mode for shared memory segment " +
                                                 m_Name);
    }

    CheckShmID("in call to ShmSystemV shmget at Open");

    m_Buffer = static_cast<char *>(shmat(m_ShmID, nullptr, 0));
    CheckBuffer("in call to SystemV shmat at Open");
    m_IsOpen = false;
}

void ShmSystemV::Write(const char *buffer, size_t size, size_t start)
{
    CheckSizes(size, start, "in call to Write");
    ProfilerStart("write");
    std::memcpy(&m_Buffer[start], buffer, size);
    ProfilerStop("write");
}

void ShmSystemV::Read(char *buffer, size_t size, size_t start)
{
    CheckSizes(size, start, "in call to Read");
    ProfilerStart("read");
    std::memcpy(buffer, &m_Buffer[start], size);
    ProfilerStop("read");
}

void ShmSystemV::Close()
{
    ProfilerStart("close");
    int result = shmdt(m_Buffer);
    ProfilerStop("close");
    if (result < 1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::shm::ShmSystemV", "Close",
                                              "failed to detach shared memory segment of size " +
                                                  std::to_string(m_Size) + " and name " + m_Name);
    }

    if (m_RemoveAtClose)
    {
        ProfilerStart("close");
        const int remove = shmctl(m_ShmID, IPC_RMID, NULL);
        ProfilerStop("close");
        if (remove < 1)
        {
            helper::Throw<std::ios_base::failure>(
                "Toolkit", "transport::shm::ShmSystemV", "Close",
                "failed to remove shared memory segment of size " + std::to_string(m_Size) +
                    " and name " + m_Name);
        }
    }

    m_IsOpen = false;
}

void ShmSystemV::SeekToEnd()
{
    // empty function. seek operation is meaningless for shared memory
}

void ShmSystemV::SeekToBegin()
{
    // empty function. seek operation is meaningless for shared memory
}

void ShmSystemV::Seek(const size_t start)
{
    // empty function. seek operation is meaningless for shared memory
}

// PRIVATE
void ShmSystemV::CheckShmID(const std::string hint) const
{
    if (m_ShmID < 0)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "transport::shm::ShmSystemV", "CheckShmID",
                                              "Failed shared memory segment of size " +
                                                  std::to_string(m_Size) + " and name " + m_Name +
                                                  ", " + hint);
    }
}

void ShmSystemV::CheckBuffer(const std::string hint) const
{
    if (m_Buffer == nullptr)
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "transport::shm::ShmSystemV", "CheckBuffer",
            "nullptr shared memory segment of size " + std::to_string(m_Size) + " and name " +
                m_Name + " " + hint);
    }
}

void ShmSystemV::CheckSizes(const size_t start, const size_t size, const std::string hint) const
{
    if (start + size > m_Size)
    {
        helper::Throw<std::invalid_argument>(
            "Toolkit", "transport::shm::ShmSystemV", "CheckSizes",
            "final position (start + size) = (" + std::to_string(start) + " + " +
                std::to_string(size) + " ) exceeding shared memory pre-allocated size:" +
                std::to_string(m_Size) + "," + hint);
    }
}

void ShmSystemV::MkDir(const std::string &fileName) {}

} // end namespace transport
} // end namespace adios2
