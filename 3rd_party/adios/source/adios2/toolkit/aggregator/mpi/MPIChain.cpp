/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MPIChain.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include "MPIChain.h"

#include "adios2/helper/adiosLog.h"
#include "adios2/toolkit/format/buffer/heap/BufferSTL.h"

namespace adios2
{
namespace aggregator
{

MPIChain::MPIChain() : MPIAggregator() {}

void MPIChain::Init(const size_t numAggregators, const size_t subStreams,
                    helper::Comm const &parentComm)
{
    /* numAggregators ignored here as BP3/BP4 uses substreams = aggregators */
    m_NumAggregators = subStreams;
    if (subStreams > 0)
    {
        InitComm(subStreams, parentComm);
        HandshakeRank(0);
    }
    else
    {
        InitCommOnePerNode(parentComm);
    }

    HandshakeLinks();

    // add a receiving buffer except for the last rank (only sends)
    if (m_Rank < m_Size)
    {
        m_Buffers.emplace_back(new format::BufferSTL()); // just one for now
    }
}

void MPIChain::Close() { MPIAggregator::Close(); }

MPIChain::ExchangeRequests MPIChain::IExchange(format::Buffer &buffer, const int step)
{
    if (m_Size == 1)
    {
        return {};
    }

    format::Buffer &sendBuffer = GetSender(buffer);
    const int endRank = m_Size - 1 - step;
    const bool sender = (m_Rank >= 1 && m_Rank <= endRank) ? true : false;
    const bool receiver = (m_Rank < endRank) ? true : false;

    ExchangeRequests requests;

    if (sender) // sender
    {
        requests.m_SendSize =
            m_Comm.Isend(&sendBuffer.m_Position, 1, m_Rank - 1, 0,
                         ", aggregation Isend size at iteration " + std::to_string(step) + "\n");

        // only send data if buffer larger than 0
        if (sendBuffer.m_Position > 0)
        {
            requests.m_SendData =
                m_Comm.Isend(sendBuffer.Data(), sendBuffer.m_Position, m_Rank - 1, 1,
                             ", aggregation Isend data at iteration " + std::to_string(step));
        }
    }
    // receive size, resize receiving buffer and receive data
    if (receiver)
    {
        size_t bufferSize = 0;
        helper::Comm::Req receiveSizeRequest =
            m_Comm.Irecv(&bufferSize, 1, m_Rank + 1, 0,
                         ", aggregation Irecv size at iteration " + std::to_string(step) + "\n");

        receiveSizeRequest.Wait(", aggregation waiting for receiver size at iteration " +
                                std::to_string(step) + "\n");

        format::Buffer &receiveBuffer = GetReceiver(buffer);
        ResizeUpdateBuffer(bufferSize, receiveBuffer,
                           "in aggregation, when resizing receiving buffer to size " +
                               std::to_string(bufferSize));

        // only receive data if buffer is larger than 0
        if (bufferSize > 0)
        {
            requests.m_RecvData =
                m_Comm.Irecv(receiveBuffer.Data(), receiveBuffer.m_Position, m_Rank + 1, 1,
                             ", aggregation Irecv data at iteration " + std::to_string(step));
        }
    }

    return requests;
}

MPIChain::ExchangeAbsolutePositionRequests
MPIChain::IExchangeAbsolutePosition(format::Buffer &buffer, const int step)
{
    if (m_Size == 1)
    {
        return {};
    }

    if (m_IsInExchangeAbsolutePosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "aggregator::mpi::MPIChain",
                                          "IExchangeAbsolutePosition",
                                          "An existing exchange is still active");
    }

    const int destination = (step != m_Size - 1) ? step + 1 : 0;
    ExchangeAbsolutePositionRequests requests;

    if (step == 0)
    {
        m_SizeSend = (m_Rank == 0) ? buffer.m_AbsolutePosition : buffer.m_Position;
    }

    if (m_Rank == step)
    {
        m_ExchangeAbsolutePosition =
            (m_Rank == 0) ? m_SizeSend : m_SizeSend + buffer.m_AbsolutePosition;

        requests.m_Send = m_Comm.Isend(&m_ExchangeAbsolutePosition, 1, destination, 0,
                                       ", aggregation Isend absolute position at iteration " +
                                           std::to_string(step) + "\n");
    }
    else if (m_Rank == destination)
    {
        requests.m_Recv = m_Comm.Irecv(&buffer.m_AbsolutePosition, 1, step, 0,
                                       ", aggregation Irecv absolute position at iteration " +
                                           std::to_string(step) + "\n");
    }

    m_IsInExchangeAbsolutePosition = true;
    return requests;
}

void MPIChain::Wait(ExchangeRequests &requests, const int step)
{
    if (m_Size == 1)
    {
        return;
    }

    const int endRank = m_Size - 1 - step;
    const bool sender = (m_Rank >= 1 && m_Rank <= endRank) ? true : false;
    const bool receiver = (m_Rank < endRank) ? true : false;

    if (receiver)
    {
        requests.m_RecvData.Wait(", aggregation waiting for receiver request at iteration " +
                                 std::to_string(step) + "\n");
    }

    if (sender)
    {
        const std::string hint =
            ", aggregation waiting for sender request at iteration " + std::to_string(step) + "\n";
        requests.m_SendSize.Wait(hint);
        requests.m_SendData.Wait(hint);
    }
}

void MPIChain::WaitAbsolutePosition(ExchangeAbsolutePositionRequests &requests, const int step)
{
    if (m_Size == 1)
    {
        return;
    }

    if (!m_IsInExchangeAbsolutePosition)
    {
        helper::Throw<std::runtime_error>("Toolkit", "aggregator::mpi::MPIChain",
                                          "WaitAbsolutePosition",
                                          "An existing exchange is not active");
    }

    const int destination = (step != m_Size - 1) ? step + 1 : 0;

    if (m_Rank == destination)
    {
        requests.m_Recv.Wait(", aggregation Irecv Wait absolute position at iteration " +
                             std::to_string(step) + "\n");
    }

    if (m_Rank == step)
    {
        requests.m_Send.Wait(", aggregation Isend Wait absolute position at iteration " +
                             std::to_string(step) + "\n");
    }
    m_IsInExchangeAbsolutePosition = false;
}

void MPIChain::SwapBuffers(const int /*step*/) noexcept
{
    m_CurrentBufferOrder = (m_CurrentBufferOrder == 0) ? 1 : 0;
}

void MPIChain::ResetBuffers() noexcept { m_CurrentBufferOrder = 0; }

format::Buffer &MPIChain::GetConsumerBuffer(format::Buffer &buffer) { return GetSender(buffer); }

// PRIVATE
void MPIChain::HandshakeLinks()
{
    int link = -1;

    helper::Comm::Req sendRequest;
    if (m_Rank > 0) // send
    {
        sendRequest = m_Comm.Isend(&m_Rank, 1, m_Rank - 1, 0,
                                   "Isend handshake with neighbor, MPIChain aggregator, at Open");
    }

    if (m_Rank < m_Size - 1) // receive
    {
        helper::Comm::Req receiveRequest = m_Comm.Irecv(
            &link, 1, m_Rank + 1, 0, "Irecv handshake with neighbor, MPIChain aggregator, at Open");

        receiveRequest.Wait("Irecv Wait handshake with neighbor, MPIChain aggregator, at Open");
    }

    if (m_Rank > 0)
    {
        sendRequest.Wait("Isend wait handshake with neighbor, MPIChain aggregator, at Open");
    }
}

format::Buffer &MPIChain::GetSender(format::Buffer &buffer)
{
    if (m_CurrentBufferOrder == 0)
    {
        return buffer;
    }
    else
    {
        return *m_Buffers.front();
    }
}

format::Buffer &MPIChain::GetReceiver(format::Buffer &buffer)
{
    if (m_CurrentBufferOrder == 0)
    {
        return *m_Buffers.front();
    }
    else
    {
        return buffer;
    }
}

void MPIChain::ResizeUpdateBuffer(const size_t newSize, format::Buffer &buffer,
                                  const std::string hint)
{
    if (buffer.m_FixedSize > 0)
    {
        if (newSize > buffer.m_FixedSize)
        {
            helper::Throw<std::invalid_argument>(
                "Toolkit", "aggregator::mpi::MPIChain", "ResizeUpdateBuffer",
                "requesting new size: " + std::to_string(newSize) +
                    " bytes, for fixed size buffer " + std::to_string(buffer.m_FixedSize) +
                    " of type " + buffer.m_Type + ", allocate more memory");
        }
        return; // do nothing if fixed size is enough
    }

    buffer.Resize(newSize, hint);
    buffer.m_Position = newSize;
}

} // end namespace aggregator
} // end namespace adios2
