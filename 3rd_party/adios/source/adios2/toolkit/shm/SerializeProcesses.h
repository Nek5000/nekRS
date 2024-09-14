/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SerializeProcesses.h
 *
 *  Created on: Oct 12, 2021
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 *
 * Tiny shared memory segment for the purpose of serializing a bunch of
 * processes. The ranks of the communicator is used for serializing. It requires
 * a communicator that connects processes on the same node only. Use
 * adios2::helper::Comm::GroupByShm() to create one if needed.
 *
 * Rank 0 does not need to Wait() first, but must call Done() to pass control to
 * Rank. Rank 0 can call Wait() after Done() and get back control when the last
 * process is Done(). Other processes can Wait() or check on IsMyTurn()
 * regularly. Nothing prevents processes to do anything they want unless they
 * enter the blocking Wait().
 */

#ifndef ADIOS2_TOOLKIT_SHM_SERIALIZEPROCESSES_H_
#define ADIOS2_TOOLKIT_SHM_SERIALIZEPROCESSES_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/helper/adiosComm.h"

#include <atomic>
#include <chrono>
#include <thread>

namespace adios2
{
namespace shm
{

class SerializeProcesses
{

public:
    SerializeProcesses(helper::Comm *comm);

    ~SerializeProcesses();

    void Wait();     // blocking wait until it's my turn
    bool IsMyTurn(); // non-blocking check if it's my turn
    void Done();     // this process is done, next please

private:
    helper::Comm *m_NodeComm;
    const int m_Rank;
    const int m_nProc;

    helper::Comm::Win m_Win;
    int *m_ShmValue; // single integer is the whole shm
};

} // end namespace shm
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_SHM_SERIALIZEPROCESSES_H_ */
