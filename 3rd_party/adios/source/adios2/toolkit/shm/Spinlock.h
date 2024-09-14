/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Spinlock.h
 *
 *  Created on: Oct 12, 2021
 *  Moved out from adios2/toolkit/aggregator/mpi/MPIShmChain.h
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 *
 */

#ifndef ADIOS2_TOOLKIT_SHM_SPINLOCK_H_
#define ADIOS2_TOOLKIT_SHM_SPINLOCK_H_

#include <atomic>

namespace adios2
{
namespace shm
{

class Spinlock
{
    /* from
     * https://wang-yimu.com/a-tutorial-on-shared-memory-inter-process-communication
     */
public:
    Spinlock();
    virtual ~Spinlock() = default;
    void lock();
    void unlock();

private:
    inline bool try_lock();
    std::atomic_flag flag_; //{ATOMIC_FLAG_INIT};
};

} // end namespace shm
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_SHM_SPINLOCK_H_ */
