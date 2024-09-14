/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManWriter.h
 *
 *  Created on: Jan 10, 2017
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_H_
#define ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_H_

#include "DataManMonitor.h"
#include "adios2/core/Engine.h"
#include "adios2/toolkit/format/dataman/DataManSerializer.tcc"
#include "adios2/toolkit/zmq/zmqpubsub/ZmqPubSub.h"
#include "adios2/toolkit/zmq/zmqreqrep/ZmqReqRep.h"
#include <atomic>

namespace adios2
{
namespace core
{
namespace engine
{

class DataManWriter : public Engine
{

public:
    DataManWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm);
    virtual ~DataManWriter();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const;
    void PerformPuts() final;
    void EndStep() final;
    void Flush(const int transportIndex = -1) final;

private:
    std::string m_IPAddress;
    int m_Port = 50001;
    int m_RendezvousReaderCount = 1;
    int m_Timeout = 5;
    int m_Verbosity = 0;
    bool m_Threading = false;
    std::string m_TransportMode = "fast";
    bool m_MonitorActive = false;
    int m_CombiningSteps = 1;
    int m_CombinedSteps = 0;
    std::string m_FloatAccuracy;

    int m_MpiRank;
    int m_MpiSize;
    size_t m_SerializerBufferSize = 1024 * 1024;
    int64_t m_CurrentStep = -1;
    std::atomic<size_t> m_SentSteps;
    nlohmann::json m_HandshakeJson;

    format::DataManSerializer m_Serializer;

    zmq::ZmqPubSub m_Publisher;
    zmq::ZmqReqRep m_Replier;

    DataManMonitor m_Monitor;

    std::thread m_ReplyThread;
    std::thread m_PublishThread;
    std::atomic<bool> m_ReplyThreadActive;
    bool m_PublishThreadActive;

    std::queue<std::shared_ptr<std::vector<char>>> m_BufferQueue;
    std::mutex m_BufferQueueMutex;

    void PushBufferQueue(std::shared_ptr<std::vector<char>> buffer);
    std::shared_ptr<std::vector<char>> PopBufferQueue();
    bool IsBufferQueueEmpty();

    void Handshake();
    void ReplyThread();
    void PublishThread();

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &, const T *) final;                                                \
    void DoPutDeferred(Variable<T> &, const T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    template <class T>
    void PutSyncCommon(Variable<T> &variable, const T *values);

    template <class T>
    void PutDeferredCommon(Variable<T> &variable, const T *values);

    void DoClose(const int transportIndex = -1) final;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DATAMAN_DATAMANWRITER_H_ */
