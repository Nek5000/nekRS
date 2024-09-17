/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SkeletonReader.cpp
 *
 *  Created on: Jan 04, 2018
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "SkeletonReader.h"
#include "SkeletonReader.tcc"

#include "adios2/helper/adiosFunctions.h" // CSVToVector

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

SkeletonReader::SkeletonReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("SkeletonReader", io, name, mode, std::move(comm))
{
    m_ReaderRank = m_Comm.Rank();
    Init();
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << " Open(" << m_Name << ") in constructor."
                  << std::endl;
    }
    m_IsOpen = true;
}

SkeletonReader::~SkeletonReader()
{
    /* m_Skeleton deconstructor does close and finalize */
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << " deconstructor on " << m_Name << "\n";
    }
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus SkeletonReader::BeginStep(const StepMode mode, const float timeoutSeconds)
{
    // step info should be received from the writer side in BeginStep()
    // so this forced increase should not be here
    ++m_CurrentStep;

    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << "   BeginStep() new step "
                  << m_CurrentStep << "\n";
    }

    // If we reach the end of stream (writer is gone or explicitly tells the
    // reader)
    // we return EndOfStream to the reader application
    if (m_CurrentStep == 2)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank
                  << "   forcefully returns End of Stream at this step\n";

        return StepStatus::EndOfStream;
    }

    // We should block until a new step arrives or reach the timeout

    // m_IO Variables and Attributes should be defined at this point
    // so that the application can inquire them and start getting data

    return StepStatus::OK;
}

void SkeletonReader::PerformGets()
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << "     PerformGets()\n";
    }
    m_NeedPerformGets = false;
}

size_t SkeletonReader::CurrentStep() const { return m_CurrentStep; }

void SkeletonReader::EndStep()
{
    // EndStep should call PerformGets() if there are unserved GetDeferred()
    // requests
    if (m_NeedPerformGets)
    {
        PerformGets();
    }

    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << "   EndStep()\n";
    }
}

// PRIVATE

#define declare_type(T)                                                                            \
    void SkeletonReader::DoGetSync(Variable<T> &variable, T *data)                                 \
    {                                                                                              \
        GetSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void SkeletonReader::DoGetDeferred(Variable<T> &variable, T *data)                             \
    {                                                                                              \
        GetDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SkeletonReader::Init()
{
    InitParameters();
    InitTransports();
}

void SkeletonReader::InitParameters()
{
    for (const auto &pair : m_IO.m_Parameters)
    {
        std::string key(pair.first);
        std::transform(key.begin(), key.end(), key.begin(), ::tolower);

        std::string value(pair.second);
        std::transform(value.begin(), value.end(), value.begin(), ::tolower);

        if (key == "verbose")
        {
            m_Verbosity = std::stoi(value);
            if (m_Verbosity < 0 || m_Verbosity > 5)
                helper::Throw<std::invalid_argument>("Engine", "SkeletonReader", "InitParameters",
                                                     "Method verbose argument must be an "
                                                     "integer in the range [0,5], in call to "
                                                     "Open or Engine constructor");
        }
    }
}

void SkeletonReader::InitTransports()
{
    // Nothing to process from m_IO.m_TransportsParameters
}

void SkeletonReader::DoClose(const int transportIndex)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Reader " << m_ReaderRank << " Close(" << m_Name << ")\n";
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
