/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MPIChain.h
 *
 *  Created on: Feb 20, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_AGGREGATOR_MPI_MPICHAIN_H_
#define ADIOS2_TOOLKIT_AGGREGATOR_MPI_MPICHAIN_H_

#include "adios2/toolkit/aggregator/mpi/MPIAggregator.h"

namespace adios2
{
namespace aggregator
{

class MPIChain : public MPIAggregator
{

public:
    /* Chain aggregator used by BP3/BP4 */
    MPIChain();

    ~MPIChain() = default;

    void Init(const size_t numAggregators, const size_t subStreams,
              helper::Comm const &parentComm) final;

    void Close() final;

    struct ExchangeRequests
    {
        helper::Comm::Req m_SendSize;
        helper::Comm::Req m_SendData;
        helper::Comm::Req m_RecvData;
    };

    ExchangeRequests IExchange(format::Buffer &buffer, const int step);

    struct ExchangeAbsolutePositionRequests
    {
        helper::Comm::Req m_Send;
        helper::Comm::Req m_Recv;
    };

    ExchangeAbsolutePositionRequests IExchangeAbsolutePosition(format::Buffer &buffer,
                                                               const int step);

    void WaitAbsolutePosition(ExchangeAbsolutePositionRequests &requests, const int step);

    void Wait(ExchangeRequests &requests, const int step);

    void SwapBuffers(const int step) noexcept;

    void ResetBuffers() noexcept;

    format::Buffer &GetConsumerBuffer(format::Buffer &buffer);

private:
    bool m_IsInExchangeAbsolutePosition = false;
    size_t m_SizeSend = 0;
    size_t m_ExchangeAbsolutePosition = 0;

    /** current sender/receiver pair order (0: sender is original buffer and
    receiver is extra buffer, 1: sender is extra buffer and receiver is the
    original buffer)*/
    unsigned int m_CurrentBufferOrder = 0;

    void HandshakeLinks();

    /**
     * Returns a reference to the sender buffer depending on
     * m_CurrentBufferOrder
     * flag
     * @param buffer original buffer from serializer
     * @return reference to sender buffer
     */
    format::Buffer &GetSender(format::Buffer &buffer);

    /**
     * Returns a reference to the receiver buffer depending on
     * m_CurrentBufferOrder
     * flag
     * @param buffer original buffer from serializer
     * @return reference to receiver buffer
     */
    format::Buffer &GetReceiver(format::Buffer &buffer);

    /**
     * Resizes and updates m_Position in a buffer, used for receiving buffers
     * @param newSize new size for receiving buffer
     * @param buffer to be resized
     * @param hint used in exception error message
     */
    void ResizeUpdateBuffer(const size_t newSize, format::Buffer &buffer, const std::string hint);
};

} // end namespace aggregator
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_AGGREGATOR_MPI_MPICHAIN_H_ */
