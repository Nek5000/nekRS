/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosCommDummy.h : A dummy Comm that does not actually communicate.
 */

#ifndef ADIOS2_HELPER_ADIOSCOMMDUMMY_H_
#define ADIOS2_HELPER_ADIOSCOMMDUMMY_H_

#include "adiosComm.h"

namespace adios2
{
namespace helper
{

/**
 * @brief Create a dummy communicator.
 */
Comm CommDummy();

} // end namespace helper
} // end namespace adios2

#endif // ADIOS2_HELPER_ADIOSCOMMDUMMY_H_
