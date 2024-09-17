/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TokenChain.h
 *
 *  Created on: Oct 12, 2021
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 *
 * A shared memory segment for the purpose of serializing a bunch of
 * processes. The ranks of the communicator is used for serializing. It requires
 * a communicator that connects processes on the same node only. Use
 * adios2::helper::Comm::GroupByShm() to create one if needed.
 *
 * Rank 0 must call SendToken() to pass a token to Rank 1.
 * Rank 0 MAY call RecvToken() to wait for the last process and get back a token
 * from it. It MUST call RecvToken() however if it ever wants to start a second
 * round so that it is synced with ending the first round.
 *
 * Other processes can blocking wait on RecvToken() or check status on
 * CheckToken() regularly. Nothing prevents processes to do anything they want
 * unless they enter the blocking wait.
 *
 * Example:

 helper:Comm comm = m_Comm.GroupByShm();
 shm::TokenChain<int> tokenChain(&comm);
 if (!rank)
 {
    int token = m_Comm.Size();
    tokenChain.SendToken(token);
    ...
    token = tokenChain.RecvToken();
 }
 else
 {
    int token = tokenChain.RecvToken();
    ...
    --token;
    tokenChain.SendToken(token);
 }

 */

#ifndef ADIOS2_TOOLKIT_SHM_TOKENCHAIN_H_
#define ADIOS2_TOOLKIT_SHM_TOKENCHAIN_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/helper/adiosComm.h"

#include <assert.h>
#include <atomic>
#include <chrono>
#include <iostream>
#include <thread>

namespace adios2
{
namespace shm
{

template <class T>
class TokenChain
{

    struct segment_
    {
        int currentRank;
        T token;
    };

public:
    TokenChain<T>(helper::Comm *comm)
    : m_NodeComm(comm), m_Rank(comm->Rank()), m_nProc(comm->Size())
    {
        if (m_nProc > 1)
        {
            char *ptr;
            if (!m_Rank)
            {
                m_Win = m_NodeComm->Win_allocate_shared(sizeof(segment_), 1, &ptr);
            }
            else
            {
                m_Win = m_NodeComm->Win_allocate_shared(0, 1, &ptr);
                size_t shmsize;
                int disp_unit;
                m_NodeComm->Win_shared_query(m_Win, 0, &shmsize, &disp_unit, &ptr);
            }
            m_Shm = reinterpret_cast<segment_ *>(ptr);

            if (!m_Rank)
            {
                m_Shm->currentRank = 0;
                m_Shm->token = T();
            }
        }
        else
        {
            m_Shm = new segment_;
            m_Shm->currentRank = 0;
            m_Shm->token = T();
        }
    }

    ~TokenChain()
    {
        if (m_nProc > 1)
        {
            m_NodeComm->Win_free(m_Win);
        }
        else
        {
            delete m_Shm;
        }
    }

    /** blocking wait until it's my turn, returns token */
    T &RecvToken()
    {
        while (m_Shm->currentRank != m_Rank)
        {
            assert(0 <= m_Shm->currentRank && m_Shm->currentRank < m_nProc);
            std::this_thread::sleep_for(std::chrono::duration<double>(0.00001));
        }
        return m_Shm->token;
    }

    /** non-blocking check if it's my turn */
    bool CheckToken()
    {
        assert(0 <= m_Shm->currentRank && m_Shm->currentRank < m_nProc);
        return (m_Shm->currentRank == m_Rank);
    }

    /** this process is done, pass token to next */
    void SendToken(T &token)
    {
        if (m_Rank != m_Shm->currentRank)
        {
            helper::Throw<std::runtime_error>("Toolkit", "shm::TokenChain", "SendToken",
                                              "function can only be "
                                              "called by the Rank who last called "
                                              "RecvToken, rank = " +
                                                  std::to_string(m_Rank));
        }
        assert(0 <= m_Shm->currentRank && m_Shm->currentRank < m_nProc);
        m_Shm->token = token;
        /* Warning: flipping the currentRank may activate other process'
         * RecvToken which in turn may call SendToken and change m_Shm
         * immediately. So this action must be the very last action here */
        if (m_Rank < m_nProc - 1)
        {
            ++(m_Shm->currentRank);
        }
        else
        {
            m_Shm->currentRank = 0;
        }
    }

private:
    helper::Comm *m_NodeComm;
    const int m_Rank;
    const int m_nProc;

    helper::Comm::Win m_Win;
    segment_ *m_Shm;
};

} // end namespace shm
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_SHM_TOKENCHAIN_H_ */
