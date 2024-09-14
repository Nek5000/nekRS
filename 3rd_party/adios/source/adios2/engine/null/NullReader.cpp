#include "NullReader.h"
#include "NullReader.tcc"

#include "adios2/helper/adiosLog.h"

namespace adios2
{
namespace core
{
namespace engine
{

struct NullReader::NullReaderImpl
{
    int64_t CurrentStep = -1;
    bool IsInStep = false;
    bool IsOpen = true;
};

NullReader::NullReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("NullReader", io, name, mode, std::move(comm)), Impl(new NullReader::NullReaderImpl)
{
    m_IsOpen = true;
}

NullReader::~NullReader()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    m_IsOpen = false;
}

StepStatus NullReader::BeginStep(StepMode mode, const float timeoutSeconds)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "BeginStep",
                                          "NullReader::BeginStep: Engine already closed");
    }

    if (Impl->IsInStep)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "BeginStep",
                                          "NullReader::BeginStep: Step already active");
    }

    Impl->IsInStep = true;
    ++Impl->CurrentStep;
    return StepStatus::EndOfStream;
}

size_t NullReader::CurrentStep() const
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "CurrentStep",
                                          "NullReader::CurrentStep: Engine already closed");
    }

    return static_cast<size_t>(Impl->CurrentStep);
}

void NullReader::EndStep()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "EndStep",
                                          "NullReader::EndStep: Engine already closed");
    }

    if (!Impl->IsInStep)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "EndStep",
                                          "NullReader::EndStep: No active step");
    }

    Impl->IsInStep = false;
}

void NullReader::PerformGets()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "PerformGets",
                                          "NullReader::PerformPuts: Engine already closed");
    }

    return;
}

void NullReader::DestructorClose(bool Verbose) noexcept { m_IsOpen = false; }

void NullReader::DoClose(const int)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Engine", "NullReader", "DoClose", "already closed");
    }

    Impl->IsOpen = false;
}

#define declare_type(T)                                                                            \
    void NullReader::DoGetSync(Variable<T> &variable, T *data) { GetSyncCommon(variable, data); }  \
    void NullReader::DoGetDeferred(Variable<T> &variable, T *data)                                 \
    {                                                                                              \
        GetDeferredCommon(variable, data);                                                         \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

} // end namespace engine
} // end namespace core
} // end namespace adios2
