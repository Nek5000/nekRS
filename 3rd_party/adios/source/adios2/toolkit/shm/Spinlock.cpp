/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Spinlock.cpp
 *
 *  Created on: Oct 12, 2021
 *  Moved out from adios2/toolkit/aggregator/mpi/MPIShmChain.h
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 *
 */

#include "Spinlock.h"

#include <chrono>
#include <thread>

namespace adios2
{
namespace shm
{

Spinlock::Spinlock() { flag_.clear(); }
void Spinlock::lock()
{
    while (!try_lock())
    {
        std::this_thread::sleep_for(std::chrono::duration<double>(0.00001));
    }
}
void Spinlock::unlock() { flag_.clear(); }

inline bool Spinlock::try_lock() { return !flag_.test_and_set(); }

} // end namespace shm
} // end namespace adios2
