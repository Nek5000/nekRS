/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  NullTransport.cpp
 *
 *  Created on: Apr 17 2019
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 */

#include "NullTransport.h"
#include "adios2/helper/adiosLog.h"

#include <cstring> // std::memset

namespace adios2
{
namespace transport
{

struct NullTransport::NullTransportImpl
{
    bool IsOpen = false;
    size_t CurPos = 0;
    size_t Capacity = 0;
};

NullTransport::NullTransport(helper::Comm const &comm)
: Transport("NULL", "NULL", comm), Impl(new NullTransportImpl)
{
}

NullTransport::~NullTransport() = default;

void NullTransport::Open(const std::string &name, const Mode openMode, const bool async,
                         const bool directio)
{
    if (Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Open",
                                          "transport is already open");
    }

    ProfilerStart("open");
    Impl->IsOpen = true;
    Impl->CurPos = 0;
    Impl->Capacity = 0;
    ProfilerStop("open");
}

void NullTransport::SetBuffer(char *buffer, size_t size) { return; }

void NullTransport::Write(const char *buffer, size_t size, size_t start)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Write",
                                          "transport is not open yet");
    }

    ProfilerStart("write");
    Impl->CurPos = start + size;
    if (Impl->CurPos > Impl->Capacity)
    {
        Impl->Capacity = Impl->CurPos;
    }
    ProfilerStop("write");
}

void NullTransport::Read(char *buffer, size_t size, size_t start)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Read",
                                          "transport is not open yet");
    }

    ProfilerStart("read");
    if (start + size > Impl->Capacity)
    {
        helper::Throw<std::out_of_range>("Toolkit", "transport::NullTransport", "Read",
                                         "size+start exceeds capacity");
    }
    std::memset(buffer, 0, size);
    Impl->CurPos = start + size;
    ProfilerStop("read");
}

size_t NullTransport::GetSize() { return Impl->Capacity; }

void NullTransport::Flush()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Flush",
                                          "transport is not open yet");
    }
}

void NullTransport::Close()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Close",
                                          "transport is not open yet");
    }

    Impl->CurPos = 0;
    Impl->Capacity = 0;
    Impl->IsOpen = false;
}

void NullTransport::Delete() { Close(); }

void NullTransport::SeekToEnd()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "SeekToEnd",
                                          "transport is not open yet");
    }
    Impl->CurPos = Impl->Capacity - 1;
}

void NullTransport::SeekToBegin()
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "SeekToBegin",
                                          "transport is not open yet");
    }
    Impl->CurPos = 0;
}

void NullTransport::Seek(const size_t start)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Seek",
                                          "transport is not open yet");
    }
    Impl->CurPos = start;
}

void NullTransport::Truncate(const size_t length)
{
    if (!Impl->IsOpen)
    {
        helper::Throw<std::runtime_error>("Toolkit", "transport::NullTransport", "Truncate",
                                          "transport is not open yet");
    }
    Impl->Capacity = length;
}

void NullTransport::MkDir(const std::string &fileName) { return; }

void NullTransport::CheckName() const { return; }

} // end namespace transport
} // end namespace adios2
