/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Timer.h
 *
 *  Created on: Apr 4, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_PROFILING_IOCHRONO_TIMER_H_
#define ADIOS2_TOOLKIT_PROFILING_IOCHRONO_TIMER_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <chrono>
#include <string>
/// \endcond

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"

#include <iostream> // myTimer

namespace adios2
{
namespace profiling
{

static std::chrono::time_point<std::chrono::high_resolution_clock> m_ADIOS2ProgStart =
    std::chrono::high_resolution_clock::now();

class Timer
{

public:
    /** process name */
    const std::string m_Process;

    bool m_Trace = false;
    /** process elapsed time */
    int64_t m_ProcessTime = 0;

    /** time unit for elapsed time from ADIOSTypes.h */
    const TimeUnit m_TimeUnit;

    /** creation timedate from std::ctime */
    std::string m_LocalTimeDate;

    /**
     * Timer object constructor using std::chrono class
     * @param process name of process to be measured
     * @param timeUnit (mus, ms, s, etc.) from ADIOSTypes.h TimeUnit
     */
    Timer(const std::string process, const TimeUnit timeUnit, const bool trace = false);

    ~Timer() = default;

    /** sets timer active to start counting */
    void Resume() noexcept;

    /**
     * Pauses timer (set to inactive)
     * @throws std::invalid_argument if Resume not previously called
     */
    void Pause();

    /** Returns TimeUnit as a short std::string  */
    std::string GetShortUnits() const noexcept;

    void AddDetail()
    {
        m_nCalls++;
        if (m_Trace)
        {
            auto relative = std::chrono::duration_cast<std::chrono::microseconds>(m_InitialTime -
                                                                                  m_ADIOS2ProgStart)
                                .count();
            auto micros =
                std::chrono::duration_cast<std::chrono::microseconds>(m_ElapsedTime - m_InitialTime)
                    .count();

            if ((micros > 10000))
            {
                if (m_Details.size() > 0)
                    m_Details += ",";

                std::ostringstream ss;

                ss << "\"" << relative / 1000.0 << "+" << micros / 1000.0 << "\"";

                m_Details += ss.str();
            }
        }
    }

    void AddToJsonStr(std::string &rankLog, const bool addComma = true) const
    {
        if (0 == m_nCalls)
            return;

        if (addComma)
        {
            rankLog += ", ";
        }
        rankLog += "\"" + m_Process + "\":{\"mus\":" + std::to_string(m_ProcessTime);
        rankLog += ", \"nCalls\":" + std::to_string(m_nCalls);

        if (500 > m_nCalls)
        {
            if (m_Details.size() > 2)
            {
                rankLog += ", \"trace\":[" + m_Details + "]";
            }
        }
        rankLog += "}";
    }

    std::string m_Details;
    uint64_t m_nCalls = 0;

private:
    /** Set at Resume */
    std::chrono::time_point<std::chrono::high_resolution_clock> m_InitialTime;

    /** Set at Pause */
    std::chrono::time_point<std::chrono::high_resolution_clock> m_ElapsedTime;

    /** Checks if m_InitialTime is set, timer is running */
    bool m_InitialTimeSet = false;

    /** called by Pause to get time between Pause and Resume */
    int64_t GetElapsedTime();
};

} // end namespace profiling
} // end namespace adios

#endif /* ADIOS2_CORE_TIMER_H_ */
