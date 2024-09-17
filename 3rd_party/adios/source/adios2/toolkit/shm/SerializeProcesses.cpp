/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SerializeProcesses.cpp
 *
 *  Created on: Oct 12, 2021
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "SerializeProcesses.h"
#include <chrono>
#include <thread>

namespace adios2
{
namespace shm
{

SerializeProcesses::SerializeProcesses(helper::Comm *comm)
: m_NodeComm(comm), m_Rank(comm->Rank()), m_nProc(comm->Size())
{
    if (m_nProc > 1)
    {
        char *ptr;
        if (!m_Rank)
        {
            m_Win = m_NodeComm->Win_allocate_shared(sizeof(int), 1, &ptr);
        }
        else
        {
            m_Win = m_NodeComm->Win_allocate_shared(0, 1, &ptr);
            size_t shmsize;
            int disp_unit;
            m_NodeComm->Win_shared_query(m_Win, 0, &shmsize, &disp_unit, &ptr);
        }
        m_ShmValue = reinterpret_cast<int *>(ptr);

        if (!m_Rank)
        {
            *m_ShmValue = 0;
        }
    }
    else
    {
        m_ShmValue = new int;
    }
};

SerializeProcesses::~SerializeProcesses()
{
    if (m_nProc > 1)
    {
        m_NodeComm->Win_free(m_Win);
    }
    else
    {
        delete m_ShmValue;
    }
}

void SerializeProcesses::Wait()
{
    while (*m_ShmValue != m_Rank)
    {
        std::this_thread::sleep_for(std::chrono::duration<double>(0.00001));
    }
}

bool SerializeProcesses::IsMyTurn() { return (*m_ShmValue == m_Rank); }

void SerializeProcesses::Done()
{
    if (m_Rank < m_NodeComm->Size() - 1)
    {
        ++(*m_ShmValue);
    }
    else
    {
        *m_ShmValue = 0;
    }
}

} // end namespace shm
} // end namespace adios2
