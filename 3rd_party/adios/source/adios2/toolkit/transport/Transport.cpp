/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Transport.cpp
 *
 *  Created on: Dec 5, 2016
 *      Author: wfg
 */

#include "Transport.h"
#include "adios2/core/CoreTypes.h"
#include <algorithm> // max

#include "adios2/helper/adiosFunctions.h" //CreateDirectory

namespace adios2
{

Transport::Transport(const std::string type, const std::string library, helper::Comm const &comm)
: m_Type(type), m_Library(library), m_Comm(comm)
{
}

void Transport::WriteV(const core::iovec *iov, const int iovcnt, size_t start)
{
    if (iovcnt > 0)
    {
        Write(static_cast<const char *>(iov[0].iov_base), iov[0].iov_len, start);
        for (int c = 1; c < iovcnt; ++c)
        {
            Write(static_cast<const char *>(iov[c].iov_base), iov[c].iov_len);
        }
    }
    else if (start != MaxSizeT)
    {
        Seek(start);
    }
}

void Transport::InitProfiler(const Mode openMode, const TimeUnit timeUnit)
{
    m_Profiler.m_IsActive = true;

    m_Profiler.m_Timers.emplace(
        std::make_pair("open", profiling::Timer("open", TimeUnit::Microseconds)));

    if (openMode == Mode::Write)
    {
        m_Profiler.m_Timers.emplace("write", profiling::Timer("write", timeUnit));

        m_Profiler.m_Bytes.emplace("write", 0);
    }
    else if (openMode == Mode::Append)
    {
        /*
        m_Profiler.Timers.emplace(
            "append", profiling::Timer("append", timeUnit));
        m_Profiler.Bytes.emplace("append", 0);
        */
        m_Profiler.m_Timers.emplace("write", profiling::Timer("write", timeUnit));

        m_Profiler.m_Bytes.emplace("write", 0);

        m_Profiler.m_Timers.emplace("read", profiling::Timer("read", timeUnit));

        m_Profiler.m_Bytes.emplace("read", 0);
    }
    else if (openMode == Mode::Read)
    {
        m_Profiler.m_Timers.emplace("read", profiling::Timer("read", timeUnit));
        m_Profiler.m_Bytes.emplace("read", 0);
    }

    m_Profiler.m_Timers.emplace("close", profiling::Timer("close", TimeUnit::Microseconds));
}

void Transport::OpenChain(const std::string &name, const Mode openMode,
                          const helper::Comm &chainComm, const bool async, const bool directio)
{
    std::invalid_argument("ERROR: " + m_Name + " transport type " + m_Type + " using library " +
                          m_Library + " doesn't implement the OpenChain function\n");
}

void Transport::SetParameters(const Params &parameters) {}

void Transport::SetBuffer(char * /*buffer*/, size_t /*size*/)
{
    std::invalid_argument("ERROR: " + m_Name + " transport type " + m_Type + " using library " +
                          m_Library + " doesn't implement the SetBuffer function\n");
}

void Transport::Flush()
{
    std::invalid_argument("ERROR: " + m_Name + " transport type " + m_Type + " using library " +
                          m_Library + " doesn't implement the Flush function\n");
}

size_t Transport::GetSize() { return 0; }

void Transport::ProfilerWriteBytes(size_t bytes) noexcept
{
    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Bytes.at("write") += bytes;
    }
}

void Transport::ProfilerStart(const std::string process) noexcept
{
    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Timers.at(process).Resume();
    }
}

void Transport::ProfilerStop(const std::string process) noexcept
{
    if (m_Profiler.m_IsActive)
    {
        m_Profiler.m_Timers.at(process).Pause();
    }
}

void Transport::CheckName() const
{
    if (m_Name.empty())
    {
        helper::Throw<std::invalid_argument>("Toolkit", "transport::Transport", "CheckName",
                                             "name can't be empty for " + m_Library +
                                                 " transport ");
    }
}

} // end namespace adios2
