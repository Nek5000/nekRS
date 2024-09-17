/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosMpiHandshake.h
 *
 *  Created on: Mar 1, 2020
 *      Author: Jason Wang
 */

#ifndef ADIOS2_HELPER_ADIOSMPIHANDSHAKE_H_
#define ADIOS2_HELPER_ADIOSMPIHANDSHAKE_H_

#include "adios2/common/ADIOSConfig.h"
#ifndef ADIOS2_HAVE_MPI
#error "Do not include adiosMpiHandshake.h without ADIOS2_HAVE_MPI."
#endif

#include <mpi.h>
#include <string>
#include <vector>

namespace adios2
{
namespace helper
{

/**
 * Start the handshake operations and wait until the rendezvous conditions
 * are reached, or timeout.
 *
 * @param filename: name of the staging stream, must be within the length of
 * maxFilenameLength
 *
 * @param mode: 'r' or 'w', read or write
 *
 * @param timeoutSeconds: timeout for the handshake, will throw exception
 * when reaching this timeout
 *
 * @param rendezvousAppCountForStream: the number of apps, including both
 * writers and readers, that will work on this stream. The function will
 * block until it receives the MPI handshake information from all these
 * apps, or until timeoutSeconds is passed.
 *
 * @param localComm: local MPI communicator for the app
 *
 * @return 3 vectors of ranks. [0] is the vector of all writer and reader ranks
 * for stream *filename*. [1] is the vector of all writer ranks for stream
 * *filename*. [2] is the vector of all reader ranks for stream *filename*.
 */
const std::vector<std::vector<int>> HandshakeRank(const std::string &filename, const char mode,
                                                  const int timeoutSeconds, MPI_Comm localComm,
                                                  int verbosity = 0);

void HandshakeComm(const std::string &filename, const char mode, const int timeoutSeconds,
                   MPI_Comm localComm, MPI_Group &streamGroup, MPI_Group &writerGroup,
                   MPI_Group &readerGroup, MPI_Comm &streamComm, MPI_Comm &writerComm,
                   MPI_Comm &readerComm, int verbosity = 0);

} // end namespace helper
} // end namespace adios2

#endif // ADIOS2_HELPER_ADIOSMPIHANDSHAKE_H_
