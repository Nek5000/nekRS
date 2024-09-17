/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManWriter.cpp
 *
 *  Created on: Jan 10, 2017
 *      Author: Jason Wang
 */

#include "DataManWriter.tcc"

namespace adios2
{
namespace core
{
namespace engine
{

DataManWriter::DataManWriter(IO &io, const std::string &name, const Mode openMode,
                             helper::Comm comm)
: Engine("DataManWriter", io, name, openMode, std::move(comm)), m_SentSteps(0),
  m_Serializer(m_Comm, (io.m_ArrayOrder == ArrayOrdering::RowMajor)), m_ReplyThreadActive(true),
  m_PublishThreadActive(true)
{

    m_MpiRank = m_Comm.Rank();
    m_MpiSize = m_Comm.Size();

    if (m_MpiSize > 1)
    {
        helper::Throw<std::logic_error>("Engine", "DataManWriter", "Open", m_Name, m_Comm.Rank());
    }

    helper::GetParameter(m_IO.m_Parameters, "IPAddress", m_IPAddress);
    helper::GetParameter(m_IO.m_Parameters, "Port", m_Port);
    helper::GetParameter(m_IO.m_Parameters, "Timeout", m_Timeout);
    helper::GetParameter(m_IO.m_Parameters, "Verbose", m_Verbosity);
    helper::GetParameter(m_IO.m_Parameters, "RendezvousReaderCount", m_RendezvousReaderCount);
    helper::GetParameter(m_IO.m_Parameters, "Threading", m_Threading);
    helper::GetParameter(m_IO.m_Parameters, "TransportMode", m_TransportMode);
    helper::GetParameter(m_IO.m_Parameters, "Monitor", m_MonitorActive);
    helper::GetParameter(m_IO.m_Parameters, "CombiningSteps", m_CombiningSteps);
    helper::GetParameter(m_IO.m_Parameters, "FloatAccuracy", m_FloatAccuracy);

    helper::Log("Engine", "DataManWriter", "Open", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    m_HandshakeJson["Threading"] = m_Threading;
    m_HandshakeJson["Transport"] = m_TransportMode;
    m_HandshakeJson["FloatAccuracy"] = m_FloatAccuracy;

    if (m_IPAddress.empty())
    {
        helper::Throw<std::invalid_argument>("Engine", "DataManWriter", "Open",
                                             "IP address not specified");
    }

    if (m_MonitorActive)
    {
        if (m_CombiningSteps < 20)
        {
            m_Monitor.SetAverageSteps(40);
        }
        else
        {
            m_Monitor.SetAverageSteps(m_CombiningSteps * 2);
        }
    }

    std::string replierAddress = "tcp://" + m_IPAddress + ":" + std::to_string(m_Port);
    std::string publisherAddress = "tcp://" + m_IPAddress + ":" + std::to_string(m_Port + 1);

    if (m_TransportMode == "fast")
    {
        m_Publisher.OpenPublisher(publisherAddress);
    }

    m_Replier.OpenReplier(replierAddress, m_Timeout, 64);

    if (m_RendezvousReaderCount == 0)
    {
        m_ReplyThreadActive = true;
        m_ReplyThread = std::thread(&DataManWriter::ReplyThread, this);
    }
    else
    {
        Handshake();
    }

    if (m_TransportMode == "reliable" && m_RendezvousReaderCount > 0)
    {
        m_ReplyThreadActive = true;
        m_ReplyThread = std::thread(&DataManWriter::ReplyThread, this);
    }

    if (m_TransportMode == "fast")
    {
        m_ReplyThreadActive = false;
    }

    if (m_Threading && m_TransportMode == "fast")
    {
        m_PublishThreadActive = true;
        m_PublishThread = std::thread(&DataManWriter::PublishThread, this);
    }
    m_IsOpen = true;
}

DataManWriter::~DataManWriter()
{
    if (not m_IsClosed)
    {
        DoClose();
    }
}

StepStatus DataManWriter::BeginStep(StepMode mode, const float timeout_sec)
{
    ++m_CurrentStep;
    if (m_CombinedSteps == 0)
    {
        m_Serializer.NewWriterBuffer(m_SerializerBufferSize);
    }

    if (m_MonitorActive)
    {
        m_Monitor.BeginStep(m_CurrentStep);
    }

    helper::Log("Engine", "DataManWriter", "BeginStep", std::to_string(CurrentStep()), 0,
                m_Comm.Rank(), 5, m_Verbosity, helper::LogMode::INFO);

    m_Serializer.AttachTimeStamp(std::chrono::duration_cast<std::chrono::milliseconds>(
                                     std::chrono::system_clock::now().time_since_epoch())
                                     .count());

    return StepStatus::OK;
}

size_t DataManWriter::CurrentStep() const { return m_CurrentStep; }

void DataManWriter::PerformPuts() {}

void DataManWriter::EndStep()
{
    helper::Log("Engine", "DataManWriter", "EndStep", std::to_string(CurrentStep()), 0,
                m_Comm.Rank(), 5, m_Verbosity, helper::LogMode::INFO);

    if (m_CurrentStep == 0)
    {
        m_Serializer.PutAttributes(m_IO);
    }

    ++m_CombinedSteps;

    if (m_CombinedSteps >= m_CombiningSteps)
    {
        m_CombinedSteps = 0;
        m_Serializer.AttachAttributesToLocalPack();
        const auto buffer = m_Serializer.GetLocalPack();
        if (buffer->size() > m_SerializerBufferSize)
        {
            m_SerializerBufferSize = buffer->size();
        }

        if (m_Threading || m_TransportMode == "reliable")
        {
            PushBufferQueue(buffer);
        }
        else
        {
            m_Publisher.Send(buffer);
        }
    }

    if (m_MonitorActive)
    {
        m_Monitor.EndStep(m_CurrentStep);
    }
}

void DataManWriter::Flush(const int transportIndex) {}

// PRIVATE functions below

#define declare_type(T)                                                                            \
    void DataManWriter::DoPutSync(Variable<T> &variable, const T *values)                          \
    {                                                                                              \
        helper::Log("Engine", "DataManWriter", "PutSync", variable.m_Name, 0, m_Comm.Rank(), 5,    \
                    m_Verbosity, helper::LogMode::INFO);                                           \
        PutSyncCommon(variable, values);                                                           \
    }                                                                                              \
    void DataManWriter::DoPutDeferred(Variable<T> &variable, const T *values)                      \
    {                                                                                              \
        helper::Log("Engine", "DataManWriter", "PutDeferred", variable.m_Name, 0, m_Comm.Rank(),   \
                    5, m_Verbosity, helper::LogMode::INFO);                                        \
        PutDeferredCommon(variable, values);                                                       \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void DataManWriter::DoClose(const int transportIndex)
{
    helper::Log("Engine", "DataManWriter", "Close", m_Name, 0, m_Comm.Rank(), 5, m_Verbosity,
                helper::LogMode::INFO);

    if (m_CombinedSteps < m_CombiningSteps && m_CombinedSteps > 0)
    {
        m_Serializer.AttachAttributesToLocalPack();
        const auto buffer = m_Serializer.GetLocalPack();
        if (buffer->size() > m_SerializerBufferSize)
        {
            m_SerializerBufferSize = buffer->size();
        }

        if (m_TransportMode == "reliable")
        {
            PushBufferQueue(buffer);
        }
        else if (m_TransportMode == "fast")
        {
            if (m_Threading)
            {
                PushBufferQueue(buffer);
            }
            else
            {
                m_Publisher.Send(buffer);
            }
        }
    }

    nlohmann::json endSignal;
    endSignal["FinalStep"] = static_cast<int64_t>(m_CurrentStep);
    std::string s = endSignal.dump() + '\0';
    auto cvp = std::make_shared<std::vector<char>>(s.size());
    std::memcpy(cvp->data(), s.c_str(), s.size());

    if (m_TransportMode == "reliable")
    {
        PushBufferQueue(cvp);
    }
    else if (m_TransportMode == "fast")
    {
        if (m_Threading)
        {
            while (!IsBufferQueueEmpty())
            {
            }
            for (int i = 0; i < 3; ++i)
            {
                PushBufferQueue(cvp);
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                m_Publisher.Send(cvp);
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
    }

    if (m_ReplyThreadActive)
    {
        while (m_SentSteps < static_cast<size_t>(m_CurrentStep + 2))
        {
        }
        m_ReplyThreadActive = false;
    }
    if (m_ReplyThread.joinable())
    {
        m_ReplyThread.join();
    }

    m_PublishThreadActive = false;
    if (m_PublishThread.joinable())
    {
        m_PublishThread.join();
    }

    m_IsClosed = true;
}

bool DataManWriter::IsBufferQueueEmpty()
{
    std::lock_guard<std::mutex> l(m_BufferQueueMutex);
    return m_BufferQueue.empty();
}

void DataManWriter::PushBufferQueue(std::shared_ptr<std::vector<char>> buffer)
{
    std::lock_guard<std::mutex> l(m_BufferQueueMutex);
    m_BufferQueue.push(buffer);
}

std::shared_ptr<std::vector<char>> DataManWriter::PopBufferQueue()
{
    std::lock_guard<std::mutex> l(m_BufferQueueMutex);
    if (m_BufferQueue.empty())
    {
        return nullptr;
    }
    else
    {
        auto ret = m_BufferQueue.front();
        m_BufferQueue.pop();
        return ret;
    }
}

void DataManWriter::PublishThread()
{
    while (m_PublishThreadActive)
    {
        auto buffer = PopBufferQueue();
        if (buffer != nullptr && buffer->size() > 0)
        {
            m_Publisher.Send(buffer);
        }
    }
}

void DataManWriter::Handshake()
{
    int readerCount = 0;
    while (true)
    {
        auto request = m_Replier.ReceiveRequest();
        if (request != nullptr && request->size() > 0)
        {
            std::string r(request->begin(), request->end());
            if (r == "Handshake")
            {
                m_HandshakeJson["TimeStamp"] =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now().time_since_epoch())
                        .count();
                std::string js = m_HandshakeJson.dump() + '\0';
                m_Replier.SendReply(js.data(), js.size());
            }
            else if (r == "Ready")
            {
                m_Replier.SendReply("OK", 2);
                ++readerCount;
            }

            if (readerCount >= m_RendezvousReaderCount)
            {
                break;
            }
        }
    }
}

void DataManWriter::ReplyThread()
{
    while (m_ReplyThreadActive)
    {
        auto request = m_Replier.ReceiveRequest();
        if (request != nullptr && request->size() > 0)
        {
            std::string r(request->begin(), request->end());
            if (r == "Handshake")
            {
                m_HandshakeJson["TimeStamp"] =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now().time_since_epoch())
                        .count();
                std::string js = m_HandshakeJson.dump() + '\0';
                m_Replier.SendReply(js.data(), js.size());
            }
            else if (r == "Ready")
            {
                m_Replier.SendReply("OK", 2);
            }
            else if (r == "Step")
            {
                auto buffer = PopBufferQueue();
                while (buffer == nullptr)
                {
                    buffer = PopBufferQueue();
                }
                if (buffer->size() > 0)
                {
                    m_Replier.SendReply(buffer);
                    m_SentSteps = m_SentSteps + m_CombiningSteps;
                }
            }
        }
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
