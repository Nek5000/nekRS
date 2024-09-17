/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SkeletonWriter.cpp
 * Skeleton engine from which any engine can be built.
 *
 *  Created on: Jan 04, 2018
 *      Author: Norbert Podhorszki pnorbert@ornl.gov
 */

#include "SkeletonWriter.h"
#include "SkeletonWriter.tcc"

#include "adios2/helper/adiosFunctions.h"

#include <iostream>

namespace adios2
{
namespace core
{
namespace engine
{

SkeletonWriter::SkeletonWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("SkeletonWriter", io, name, mode, std::move(comm))
{
    m_WriterRank = m_Comm.Rank();
    Init();
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << " Open(" << m_Name << ")." << std::endl;
    }
    m_IsOpen = true;
}

SkeletonWriter::~SkeletonWriter()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus SkeletonWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    m_CurrentStep++; // 0 is the first step
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "   BeginStep() new step "
                  << m_CurrentStep << "\n";
    }
    return StepStatus::OK;
}

size_t SkeletonWriter::CurrentStep() const
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "   CurrentStep() returns "
                  << m_CurrentStep << "\n";
    }
    return m_CurrentStep;
}

/* PutDeferred = PutSync, so nothing to be done in PerformPuts */
void SkeletonWriter::PerformPuts()
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "     PerformPuts()\n";
    }
    m_NeedPerformPuts = false;
}

void SkeletonWriter::EndStep()
{
    if (m_NeedPerformPuts)
    {
        PerformPuts();
    }
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "   EndStep()\n";
    }
}
void SkeletonWriter::Flush(const int transportIndex)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << "   Flush()\n";
    }
}

// PRIVATE

#define declare_type(T)                                                                            \
    void SkeletonWriter::DoPutSync(Variable<T> &variable, const T *data)                           \
    {                                                                                              \
        PutSyncCommon(variable, variable.SetBlockInfo(data, CurrentStep()));                       \
        variable.m_BlocksInfo.clear();                                                             \
    }                                                                                              \
    void SkeletonWriter::DoPutDeferred(Variable<T> &variable, const T *data)                       \
    {                                                                                              \
        PutDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SkeletonWriter::Init()
{
    InitParameters();
    InitTransports();
}

void SkeletonWriter::InitParameters()
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
                helper::Throw<std::invalid_argument>("Engine", "SkeletonWriter", "InitParameters",
                                                     "Method verbose argument must be an "
                                                     "integer in the range [0,5], in call to "
                                                     "Open or Engine constructor");
        }
    }
}

void SkeletonWriter::InitTransports()
{
    // Nothing to process from m_IO.m_TransportsParameters
}

void SkeletonWriter::DoClose(const int transportIndex)
{
    if (m_Verbosity == 5)
    {
        std::cout << "Skeleton Writer " << m_WriterRank << " Close(" << m_Name << ")\n";
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
