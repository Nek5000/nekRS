/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BufferSystemV.cpp
 *
 *  Created on: Jul 9, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "BufferSystemV.h"
#include "adios2/helper/adiosLog.h"

#include <assert.h>
#include <cstring> //std::memcpy
#include <ios>     //std::ios_base::failure

#include <sys/ipc.h>   //ftok
#include <sys/shm.h>   //shmget, shmmat
#include <sys/types.h> //key_t

namespace adios2
{
namespace format
{

BufferSystemV::BufferSystemV(const size_t fixedSize, const std::string &name,
                             const unsigned int projectID, const bool remove)
: Buffer("BufferSystemV", fixedSize), m_Remove(remove)
{
    assert(projectID > 0); // for the developer
    key_t key = ftok(name.c_str(), static_cast<int>(projectID));
    m_ShmID = shmget(key, static_cast<unsigned long int>(fixedSize), IPC_CREAT | 0666);
    if (m_ShmID == -1)
    {
        helper::Throw<std::ios_base::failure>("Toolkit", "format::buffer::ipc::BufferSystemV",
                                              "BufferSystemV",
                                              "could not create shared memory buffer of size " +
                                                  std::to_string(fixedSize) + " with shmget");
    }

    void *data = shmat(m_ShmID, nullptr, 0);
    int *status = reinterpret_cast<int *>(data);
    if (*status == -1)
    {
        helper::Throw<std::runtime_error>("Toolkit", "format::buffer::ipc::BufferSystemV",
                                          "BufferSystemV",
                                          "could not attach shared memory buffer "
                                          "to address with shmat");
    }
    m_Data = static_cast<char *>(data);
}

BufferSystemV::~BufferSystemV()
{
    shmdt(m_Data);

    if (m_Remove)
    {
        shmctl(m_ShmID, IPC_RMID, NULL);
    }
}

char *BufferSystemV::Data() noexcept { return m_Data; }

const char *BufferSystemV::Data() const noexcept { return m_Data; }

void BufferSystemV::Reset(const bool resetAbsolutePosition, const bool zeroInitialize)
{
    m_Position = 0;
    if (resetAbsolutePosition)
    {
        m_AbsolutePosition = 0;
    }
    if (zeroInitialize)
    {
        memset(m_Data, 0, m_FixedSize);
    }
}

} // end namespace format
} // end namespace adios2
