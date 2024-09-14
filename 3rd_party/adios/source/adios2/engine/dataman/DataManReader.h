/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManReader.h
 *
 *  Created on: Feb 21, 2017
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_DATAMAN_DATAMANREADER_H_
#define ADIOS2_ENGINE_DATAMAN_DATAMANREADER_H_

#include <atomic>

#include "adios2/core/Engine.h"
#include "adios2/engine/dataman/DataManMonitor.h"
#include "adios2/toolkit/format/dataman/DataManSerializer.tcc"
#include "adios2/toolkit/zmq/zmqpubsub/ZmqPubSub.h"
#include "adios2/toolkit/zmq/zmqreqrep/ZmqReqRep.h"

namespace adios2
{
namespace core
{
namespace engine
{

class DataManReader : public Engine
{

public:
    DataManReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm);
    virtual ~DataManReader();
    StepStatus BeginStep(StepMode stepMode, const float timeoutSeconds) final;
    size_t CurrentStep() const final;
    void PerformGets() final;
    void EndStep() final;
    void Flush(const int transportIndex = -1) final;

private:
    std::string m_IPAddress;
    int m_Port = 50001;
    int m_Timeout = 5;
    int m_Verbosity = 0;
    bool m_Threading = true;
    uint64_t m_ReceiverBufferSize = 128 * 1024 * 1024;
    std::string m_TransportMode = "fast";
    bool m_MonitorActive = false;

    int m_MpiRank;
    int m_MpiSize;
    int64_t m_CurrentStep = -1;
    std::atomic<int64_t> m_FinalStep;
    format::DmvVecPtr m_CurrentStepMetadata;

    format::DataManSerializer m_Serializer;

    zmq::ZmqPubSub m_Subscriber;
    zmq::ZmqReqRep m_Requester;

    DataManMonitor m_Monitor;

    std::thread m_SubscriberThread;
    std::thread m_RequesterThread;

    std::atomic<bool> m_RequesterThreadActive;
    std::atomic<bool> m_SubscriberThreadActive;

    void SubscribeThread();
    void RequestThread();

    void DoClose(const int transportIndex = -1) final;

    /**
     * Called if destructor is called on an open engine.  Should warn or take
     * any non-complex measure that might help recover.
     */
    void DestructorClose(bool Verbose) noexcept final{};

#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;                                                  \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> DoAllStepsBlocksInfo(              \
        const Variable<T> &variable) const final;                                                  \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    template <typename T>
    void AccumulateMinMax(T &min, T &max, const std::vector<char> &minVec,
                          const std::vector<char> &maxVec) const;

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

    template <typename T>
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>>
    AllStepsBlocksInfoCommon(const Variable<T> &variable) const;

    template <typename T>
    std::vector<typename Variable<T>::BPInfo> BlocksInfoCommon(const Variable<T> &variable,
                                                               const size_t step) const;

    template <typename T>
    void CheckIOVariable(const std::string &name, const Dims &shape, const Dims &start,
                         const Dims &count);
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DATAMAN_DATAMANREADER_H_ */
