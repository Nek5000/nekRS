/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MhsWriter.cpp
 *
 *  Created on: Apr 6, 2019
 *      Author: Jason Wang w4g@ornl.gov
 */

#include "MhsWriter.tcc"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/operator/compress/CompressSirius.h"

namespace adios2
{
namespace core
{
namespace engine
{

MhsWriter::MhsWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("MhsWriter", io, name, mode, std::move(comm))
{
    helper::GetParameter(io.m_Parameters, "Tiers", m_Tiers);
    for (const auto &transportParams : io.m_TransportsParameters)
    {
        auto itVar = transportParams.find("variable");
        if (itVar == transportParams.end())
        {
            continue;
        }
        auto itTransport = transportParams.find("transport");
        if (itTransport == transportParams.end())
        {
            continue;
        }

        if (itTransport->second == "sirius")
        {
            m_TransportMap.emplace(itVar->second,
                                   std::make_shared<compress::CompressSirius>(io.m_Parameters));
        }
        else
        {
            helper::Throw<std::invalid_argument>("Engine", "MhsWriter", "MhsWriter",
                                                 "invalid operator");
        }
    }
    for (int i = 0; i < m_Tiers; ++i)
    {
        m_SubIOs.emplace_back(&io.m_ADIOS.DeclareIO("SubIO" + std::to_string(i)));
        m_SubEngines.emplace_back(
            &m_SubIOs.back()->Open(m_Name + ".tier" + std::to_string(i), adios2::Mode::Write));
    }
    m_IsOpen = true;
}

MhsWriter::~MhsWriter()
{
    for (int i = 0; i < m_Tiers; ++i)
    {
        m_IO.m_ADIOS.RemoveIO("SubIO" + std::to_string(i));
    }
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus MhsWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    for (auto &e : m_SubEngines)
    {
        e->BeginStep(mode, timeoutSeconds);
    }
    return StepStatus::OK;
}

size_t MhsWriter::CurrentStep() const { return m_SubEngines[0]->CurrentStep(); }

void MhsWriter::PerformPuts()
{
    for (auto &e : m_SubEngines)
    {
        e->PerformPuts();
    }
}

void MhsWriter::EndStep()
{
    for (auto &e : m_SubEngines)
    {
        e->EndStep();
    }
}

void MhsWriter::Flush(const int transportIndex)
{
    for (auto &e : m_SubEngines)
    {
        e->Flush(transportIndex);
    }
}

// PRIVATE

#define declare_type(T)                                                                            \
    void MhsWriter::DoPutSync(Variable<T> &variable, const T *data)                                \
    {                                                                                              \
        PutSyncCommon(variable, data);                                                             \
    }                                                                                              \
    void MhsWriter::DoPutDeferred(Variable<T> &variable, const T *data)                            \
    {                                                                                              \
        PutDeferredCommon(variable, data);                                                         \
    }
ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void MhsWriter::DoClose(const int transportIndex)
{
    for (auto &e : m_SubEngines)
    {
        e->Close();
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
