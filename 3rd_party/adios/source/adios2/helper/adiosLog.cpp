/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosLog.cpp
 *
 *  Created on: Nov 15, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "adiosLog.h"
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace adios2
{
namespace helper
{

std::string timeColor = "\033[1;36m";
std::string outputColor = "\033[1;32m";
std::string warningColor = "\033[1;33m";
std::string errorColor = "\033[1;31m";
std::string exceptionColor = "\033[1;34m";
std::string defaultColor = "\033[0m";

std::unordered_set<std::string> messages;

std::string MakeMessage(const std::string &component, const std::string &source,
                        const std::string &activity, const std::string &message, const int commRank,
                        const LogMode mode)
{
    std::stringstream m;

    auto timeNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    char buf[30];
    struct tm now_tm;
#ifdef _WIN32
    localtime_s(&now_tm, &timeNow);
#else
    localtime_r(&timeNow, &now_tm);
#endif
    strftime(buf, sizeof(buf), "%a %b %d %H:%M:%S %Y", &now_tm);

    m << timeColor << "[" << buf << "]";

    if (mode == INFO)
    {
        m << outputColor << " [ADIOS2 INFO]" << defaultColor;
    }
    else if (mode == WARNING)
    {
        m << warningColor << " [ADIOS2 WARNING]" << defaultColor;
    }
    else if (mode == FATALERROR)
    {
        m << errorColor << " [ADIOS2 ERROR]" << defaultColor;
    }
    else if (mode == EXCEPTION)
    {
        m << exceptionColor << " [ADIOS2 EXCEPTION]" << defaultColor;
    }

    if (commRank >= 0)
    {
        m << " [Rank " << commRank << "]";
    }

    m << " <" << component << "> <" << source << "> <" << activity << "> : " << message
      << defaultColor << std::endl;

    return m.str();
}

void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const LogMode mode)
{
    Log(component, source, activity, message, -1, -1, 0, 0, mode);
}

void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const int priority, const int verbosity, const LogMode mode)
{
    Log(component, source, activity, message, -1, -1, priority, verbosity, mode);
}

void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const int logRank, const int commRank, const int priority,
         const int verbosity, const LogMode mode)
{

    // don't print if
    // 1. logRank does not meet commRank, or
    // 2. priority does not meet verbosity, or
    // 3. the error or warning has been already printed
    if ((logRank >= 0 && commRank >= 0 && logRank != commRank) || priority > verbosity ||
        (messages.find(message) != messages.end() &&
         (mode == LogMode::FATALERROR || mode == LogMode::WARNING)))
    {
        return;
    }

    messages.insert(message);

    auto m = MakeMessage(component, source, activity, message, commRank, mode);

    if (mode == INFO || mode == WARNING)
    {
        std::cout << m;
    }
    else if (mode == FATALERROR)
    {
        std::cerr << m;
    }
}

} // end namespace helper
} // end namespace adios2
