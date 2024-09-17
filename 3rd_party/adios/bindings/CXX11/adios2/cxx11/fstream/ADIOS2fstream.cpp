/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2fstream.cpp
 *
 *  Created on: Mar 5, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "ADIOS2fstream.h"
#include "ADIOS2fstream.tcc"

namespace adios2
{

fstream::fstream(const std::string &name, const openmode mode, const std::string engineType)
: m_Stream(std::make_shared<core::Stream>(name, ToMode(mode), engineType, "C++"))
{
}

fstream::fstream(const std::string &name, const openmode mode, const std::string &configFile,
                 const std::string ioInConfigFile)
: m_Stream(std::make_shared<core::Stream>(name, ToMode(mode), configFile, ioInConfigFile, "C++"))
{
}

void fstream::open(const std::string &name, const openmode mode, const std::string engineType)
{
    CheckOpen(name);
    m_Stream = std::make_shared<core::Stream>(name, ToMode(mode), engineType, "C++");
}

void fstream::open(const std::string &name, const openmode mode, const std::string configFile,
                   const std::string ioInConfigFile)
{
    CheckOpen(name);
    m_Stream =
        std::make_shared<core::Stream>(name, ToMode(mode), configFile, ioInConfigFile, "C++");
}

fstream::operator bool() const noexcept
{
    if (!m_Stream)
    {
        return false;
    }

    return true;
}

void fstream::end_step() { m_Stream->EndStep(); }

void fstream::close()
{
    m_Stream->Close();
    m_Stream.reset();
}

void fstream::set_parameter(const std::string key, const std::string value) noexcept
{
    m_Stream->m_IO->SetParameter(key, value);
}

void fstream::set_parameters(const adios2::Params &parameters) noexcept
{
    m_Stream->m_IO->SetParameters(parameters);
}

bool getstep(adios2::fstream &stream, adios2::fstep &step)
{
    step = stream;
    return step.m_Stream->GetStep();
}

size_t fstream::current_step() const noexcept { return m_Stream->CurrentStep(); }

adios2::Mode fstream::ToMode(const openmode mode) const noexcept
{
    adios2::Mode modeCpp = adios2::Mode::Undefined;
    switch (mode)
    {
    case (openmode::out):
        modeCpp = adios2::Mode::Write;
        break;
    case (openmode::in):
        modeCpp = adios2::Mode::Read;
        break;
    case (openmode::in_random_access):
        modeCpp = adios2::Mode::ReadRandomAccess;
        break;
    case (openmode::app):
        modeCpp = adios2::Mode::Append;
        break;
    }
    return modeCpp;
}

size_t fstream::steps() const { return m_Stream->Steps(); }

// PRIVATE
void fstream::CheckOpen(const std::string &name) const
{
    if (m_Stream)
    {
        throw std::invalid_argument("ERROR: adios2::fstream with name " + name +
                                    " is already opened, in call to open");
    }
}

#define declare_template_instantiation(T)                                                          \
    template void fstream::write_attribute<T>(const std::string &, const T &, const std::string &, \
                                              const std::string, const bool);                      \
                                                                                                   \
    template void fstream::write_attribute<T>(const std::string &, const T *, const size_t,        \
                                              const std::string &, const std::string, const bool); \
                                                                                                   \
    template std::vector<T> fstream::read_attribute<T>(const std::string &, const std::string &,   \
                                                       const std::string);

ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

#define declare_template_instantiation(T)                                                          \
    template void fstream::write<T>(const std::string &, const T *, const Dims &, const Dims &,    \
                                    const Dims &, const bool);                                     \
                                                                                                   \
    template void fstream::write<T>(const std::string &, const T *, const Dims &, const Dims &,    \
                                    const Dims &, const vParams &, const bool);                    \
                                                                                                   \
    template void fstream::write<T>(const std::string &, const T &, const bool, const bool);       \
                                                                                                   \
    template std::vector<T> fstream::read<T>(const std::string &, const size_t);                   \
                                                                                                   \
    template std::vector<T> fstream::read<T>(const std::string &, const size_t, const size_t,      \
                                             const size_t);                                        \
                                                                                                   \
    template std::vector<T> fstream::read<T>(const std::string &, const Dims &, const Dims &,      \
                                             const size_t);                                        \
                                                                                                   \
    template std::vector<T> fstream::read<T>(const std::string &, const Dims &, const Dims &,      \
                                             const size_t, const size_t, const size_t);            \
                                                                                                   \
    template void fstream::read<T>(const std::string &, T *, const size_t);                        \
                                                                                                   \
    template void fstream::read(const std::string &, T *, const size_t, const size_t,              \
                                const size_t);                                                     \
                                                                                                   \
    template void fstream::read<T>(const std::string &name, T &, const size_t);                    \
                                                                                                   \
    template void fstream::read<T>(const std::string &name, T &, const size_t, const size_t);      \
                                                                                                   \
    template void fstream::read<T>(const std::string &, T *, const Dims &, const Dims &,           \
                                   const size_t);                                                  \
                                                                                                   \
    template void fstream::read<T>(const std::string &, T *, const Dims &, const Dims &,           \
                                   const size_t, const size_t, const size_t);

ADIOS2_FOREACH_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace adios2
