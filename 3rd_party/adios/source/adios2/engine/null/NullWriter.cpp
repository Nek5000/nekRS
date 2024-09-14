#include "NullWriter.h"
#include "NullWriter.tcc"

#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace core
{
namespace engine
{

struct NullWriter::NullWriterImpl
{
    size_t CurrentStep = 0;
    bool IsInStep = false;
    bool IsOpen = true;
};

NullWriter::NullWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("NullWriter", io, name, mode, std::move(comm)), Impl(new NullWriter::NullWriterImpl)
{
    m_IsOpen = true;
}

NullWriter::~NullWriter()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus NullWriter::BeginStep(StepMode mode, const float timeoutSeconds)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "BeginStep",
                                          "NullWriter::BeginStep: Engine already closed");
    }

    if (Impl->IsInStep)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "BeginStep",
                                          "NullWriter::BeginStep: Step already active");
    }

    Impl->IsInStep = true;
    ++Impl->CurrentStep;
    return StepStatus::OK;
}

size_t NullWriter::CurrentStep() const
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "CurrentStep",
                                          "NullWriter::CurrentStep: Engine already closed");
    }

    return Impl->CurrentStep;
}

void NullWriter::EndStep()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "EndStep",
                                          "NullWriter::EndStep: Engine already closed");
    }

    if (!Impl->IsInStep)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "EndStep",
                                          "NullWriter::EndStep: No active step");
    }

    Impl->IsInStep = false;
}

void NullWriter::PerformPuts()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "PerformPuts",
                                          "NullWriter::PerformPuts: Engine already closed");
    }

    return;
}

void NullWriter::Flush(const int)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "Flush",
                                          "NullWriter::Flush: Engine already closed");
    }

    return;
}

void NullWriter::DoClose(const int)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullWriter", "DoClose", "already closed");
    }

    Impl->IsOpen = false;
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
