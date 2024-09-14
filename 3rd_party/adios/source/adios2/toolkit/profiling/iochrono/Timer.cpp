/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Timer.cpp
 *
 *  Created on: Apr 4, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "Timer.h"

#include "adios2/helper/adiosFunctions.h" //LocalTimeDate

namespace adios2
{
namespace profiling
{

Timer::Timer(const std::string process, const TimeUnit timeUnit, const bool trace)
: m_Process(process), m_TimeUnit(timeUnit), m_LocalTimeDate(helper::LocalTimeDate())
{
    if (trace)
        m_Trace = true;
}

void Timer::Resume() noexcept
{
    m_InitialTime = std::chrono::high_resolution_clock::now();
    m_InitialTimeSet = true;
}

void Timer::Pause()
{
    m_ElapsedTime = std::chrono::high_resolution_clock::now();
    m_ProcessTime += GetElapsedTime();

    AddDetail();
}

std::string Timer::GetShortUnits() const noexcept
{
    std::string units;
    switch (m_TimeUnit)
    {
    case TimeUnit::Microseconds:
        units = "mus";
        break;
    case TimeUnit::Milliseconds:
        units = "ms";
        break;
    case TimeUnit::Seconds:
        units = "s";
        break;
    case TimeUnit::Minutes:
        units = "m";
        break;
    case TimeUnit::Hours:
        units = "h";
        break;
    }
    return units;
}

// PRIVATE
int64_t Timer::GetElapsedTime()
{
    if (!m_InitialTimeSet)
    {
        helper::Throw<std::invalid_argument>("Toolkit", "profiling::iochrono::Timer",
                                             "GetElapsedTime",
                                             "Resume() in process " + m_Process + " not called");
    }

    int64_t time = -1;

    switch (m_TimeUnit)
    {

    case TimeUnit::Microseconds:
        time = std::chrono::duration_cast<std::chrono::microseconds>(m_ElapsedTime - m_InitialTime)
                   .count();
        break;

    case TimeUnit::Milliseconds:
        time = std::chrono::duration_cast<std::chrono::milliseconds>(m_ElapsedTime - m_InitialTime)
                   .count();
        break;

    case TimeUnit::Seconds:
        time =
            std::chrono::duration_cast<std::chrono::seconds>(m_ElapsedTime - m_InitialTime).count();
        break;

    case TimeUnit::Minutes:
        time =
            std::chrono::duration_cast<std::chrono::minutes>(m_ElapsedTime - m_InitialTime).count();
        break;

    case TimeUnit::Hours:
        time =
            std::chrono::duration_cast<std::chrono::hours>(m_ElapsedTime - m_InitialTime).count();
        break;
    }

    return time;
}

} // end namespace profiling
} // end namespace adios2
