/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * MhsReader.cpp
 *
 *  Created on: Aug 04, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "MhsReader.tcc"
#include "adios2/helper/adiosFunctions.h"

namespace adios2
{
namespace core
{
namespace engine
{

MhsReader::MhsReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("MhsReader", io, name, mode, std::move(comm))
{
    helper::GetParameter(io.m_Parameters, "Tiers", m_Tiers);
    Params params = {{"Tiers", std::to_string(m_Tiers)}};
    m_SiriusCompressor = std::make_shared<compress::CompressSirius>(params);
    io.SetEngine("");
    m_SubIOs.emplace_back(&io);
    m_SubEngines.emplace_back(&io.Open(m_Name + ".tier0", adios2::Mode::Read));

    for (int i = 1; i < m_Tiers; ++i)
    {
        m_SubIOs.emplace_back(&io.m_ADIOS.DeclareIO("SubIO" + std::to_string(i)));
        m_SubEngines.emplace_back(
            &m_SubIOs.back()->Open(m_Name + ".tier" + std::to_string(i), adios2::Mode::Read));
    }
    m_IsOpen = true;
}

MhsReader::~MhsReader()
{
    for (int i = 1; i < m_Tiers; ++i)
    {
        m_IO.m_ADIOS.RemoveIO("SubIO" + std::to_string(i));
    }
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus MhsReader::BeginStep(const StepMode mode, const float timeoutSeconds)
{
    bool endOfStream = false;
    for (auto &e : m_SubEngines)
    {
        auto status = e->BeginStep(mode, timeoutSeconds);
        if (status == StepStatus::EndOfStream)
        {
            endOfStream = true;
        }
    }
    if (endOfStream)
    {
        return StepStatus::EndOfStream;
    }
    return StepStatus::OK;
}

size_t MhsReader::CurrentStep() const { return m_SubEngines[0]->CurrentStep(); }

void MhsReader::PerformGets()
{
    for (auto &e : m_SubEngines)
    {
        e->PerformGets();
    }
}

void MhsReader::EndStep()
{
    for (auto &e : m_SubEngines)
    {
        e->EndStep();
    }
}

// PRIVATE

#define declare_type(T)                                                                            \
    void MhsReader::DoGetSync(Variable<T> &variable, T *data) { GetSyncCommon(variable, data); }   \
    void MhsReader::DoGetDeferred(Variable<T> &variable, T *data)                                  \
    {                                                                                              \
        GetDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void MhsReader::DoClose(const int transportIndex)
{
    for (auto &e : m_SubEngines)
    {
        e->Close();
    }
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
