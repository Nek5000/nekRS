/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosLog.h
 *
 *  Created on: Nov 15, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_HELPER_ADIOSLOG_H_
#define ADIOS2_HELPER_ADIOSLOG_H_

#include <string>

namespace adios2
{
namespace helper
{

enum LogMode : char
{
    EXCEPTION = 'x',
    FATALERROR = 'e',
    WARNING = 'w',
    INFO = 'i'
};

std::string MakeMessage(const std::string &component, const std::string &source,
                        const std::string &activity, const std::string &message, const int commRank,
                        const LogMode mode);

/**
 * Print outputs, warnings, errors, and exceptions
 * @param component: Engine, Transport, Operator, etc.
 * @param source: class name of component
 * @param activity: function name where this is called
 * @param message: text message
 * @param mode: INFO, WARNING or ERROR
 */
void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const LogMode mode);

/**
 * Print outputs, warnings, errors, and exceptions
 * @param component: Engine, Transport, Operator, etc.
 * @param source: class name of component
 * @param activity: function name where this is called
 * @param message: text message
 * @param priority: only print if(priority<=verbosity)
 * @param verbosity: engine parameter for engine wide verbosity level
 * @param mode: INFO, WARNING or ERROR
 */
void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const int priority, const int verbosity, const LogMode mode);

/**
 * Print outputs, warnings, errors, and exceptions
 * @param component: Engine, Transport, Operator, etc.
 * @param source: class name of component
 * @param activity: function name where this is called
 * @param message: text message
 * @param logRank: only print if(logRank==commRank)
 * @param commRank: current MPI rank
 * @param priority: only print if(priority<=verbosity)
 * @param verbosity: engine parameter for engine wide verbosity level
 * @param mode: INFO, WARNING or ERROR
 */
void Log(const std::string &component, const std::string &source, const std::string &activity,
         const std::string &message, const int logRank, const int commRank, const int priority,
         const int verbosity, const LogMode mode);

template <class T>
void Throw(const std::string &component, const std::string &source, const std::string &activity,
           const std::string &message, const int commRank = -1)
{
    auto m = MakeMessage(component, source, activity, message, commRank, LogMode::EXCEPTION);
    throw(T(m));
}

template <class T>
void ThrowNested(const std::string &component, const std::string &source,
                 const std::string &activity, const std::string &message, const int commRank = -1)
{
    auto m = MakeMessage(component, source, activity, message, commRank, LogMode::EXCEPTION);
    throw_with_nested(T(m));
}

} // end namespace helper
} // end namespace adios2

#endif /* ADIOS2_HELPER_ADIOSLOG_H_ */
