/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IOChrono.h
 *
 *  Created on: Mar 9, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_PROFILING_IOCHRONO_IOCHRONO_H_
#define ADIOS2_TOOLKIT_PROFILING_IOCHRONO_IOCHRONO_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <unordered_map>
#include <vector>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/profiling/iochrono/Timer.h"

namespace adios2
{
namespace profiling
{

/** Class used to track different timers */
class IOChrono
{

public:
    /**
     * Create timers for each process
     * <pre>
     * Key: process name
     * Value: see Timer class public API, use to track process time as a
     * chronometer with Resume() and Pause() functions
     * </pre>
     */
    std::unordered_map<std::string, Timer> m_Timers;

    /** Create byte tracking counter for each process*/
    std::unordered_map<std::string, size_t> m_Bytes;

    /** flag to determine if IOChrono object is being used */
    bool m_IsActive = false;

    IOChrono() = default;
    ~IOChrono() = default;

    /** Start existing process in m_Timers */
    void Start(const std::string process) noexcept;

    /**
     * Stop existing process in m_Timers
     * @throws std::invalid_argument if Start wasn't called
     * */
    void Stop(const std::string process);
};

class JSONProfiler
{
public:
    JSONProfiler(helper::Comm const &comm);
    void Gather();
    void AddTimerWatch(const std::string &, const bool trace = false);

    void Start(const std::string process) { m_Profiler.Start(process); };
    void Stop(const std::string process) { m_Profiler.Stop(process); };
    void AddBytes(const std::string process, size_t bytes)
    {
        m_Profiler.m_Bytes[process] += bytes;
    };

    std::string GetRankProfilingJSON(
        const std::vector<std::string> &transportsTypes,
        const std::vector<adios2::profiling::IOChrono *> &transportsProfilers) noexcept;

    std::vector<char> AggregateProfilingJSON(const std::string &rankLog) const;

private:
    IOChrono m_Profiler;
    int m_RankMPI = 0;
    helper::Comm const &m_Comm;
};

} // end namespace profiling
} // end namespace adios

#endif /* ADIOS2_CORE_IOCHRONO_H_ */
