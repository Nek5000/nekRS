/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DataManCommon.h
 *
 *  Created on: Jun 2, 2020
 *      Author: Jason Wang
 */

#ifndef ADIOS2_ENGINE_DATAMAN_DATAMANMONITOR_H_
#define ADIOS2_ENGINE_DATAMAN_DATAMANMONITOR_H_

#include <chrono>
#include <mutex>
#include <queue>
#include <string>

namespace adios2
{
namespace core
{
namespace engine
{

class DataManMonitor
{
public:
    void BeginStep(const size_t step);
    void EndStep(const size_t step);
    void AddLatencyMilliseconds(const uint64_t remoteStamp);
    void AddBytes(const size_t bytes);
    void SetAverageSteps(const size_t steps);
    void SetCombiningSteps(const size_t steps);
    void SetWriterThreading();
    void SetReaderThreading();
    void SetClockError(const uint64_t roundLatency, const uint64_t remoteTimeBase);
    void AddCompression(const std::string &method, const std::string &accuracyUsed);
    void SetRequiredAccuracy(const std::string &accuracyRequired);
    void AddTransport(const std::string &method);
    void OutputJson(const std::string &filename);
    void OutputCsv(const std::string &filename);

private:
    bool FileExisted(const std::string &filename);

    using TimePoint = std::chrono::time_point<std::chrono::system_clock>;
    TimePoint m_InitialTimer;
    std::queue<TimePoint> m_StepTimers;
    std::queue<size_t> m_TotalBytes;
    std::deque<uint64_t> m_LatencyMilliseconds;
    size_t m_StepBytes;

    std::mutex m_PrintMutex;

    size_t m_AverageSteps = 50;
    size_t m_CombiningSteps = 1;
    int64_t m_CurrentStep = -1;
    uint64_t m_ClockError = 0;

    double m_RoundLatency = 0.0;
    double m_TotalTime = 0;
    double m_AverageTime = 0;
    double m_TotalRate = 0;
    double m_AverageRate = 0;
    double m_DropRate = 0;
    double m_StepsPerSecond = 0;
    double m_AccumulatedLatency = 0;
    std::string m_CompressionMethod;
    float m_CompressionAccuracy = 0;
    float m_RequiredAccuracy = 0.00001;
    std::string m_TransportMethod;
    bool m_ReaderThreading = false;
    bool m_WriterThreading = false;

    bool m_Verbose = true;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DATAMAN_DATAMANMONITOR_H_ */
