/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11IO.cpp
 *
 *  Created on: Mar 14, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#include "py11IO.h"

#include "adios2/common/ADIOSMacros.h"
#include "adios2/helper/adiosFunctions.h" //GetDataType<T>

#include "py11types.h"

namespace adios2
{
namespace py11
{

IO::IO(core::IO *io) : m_IO(io) {}

IO::operator bool() const noexcept { return (m_IO == nullptr) ? false : true; }

bool IO::InConfigFile() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::InConfigFile");
    return m_IO->InConfigFile();
}

void IO::SetEngine(const std::string type)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetEngine");
    m_IO->SetEngine(type);
}

void IO::SetParameter(const std::string key, const std::string value)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetParameter");
    m_IO->SetParameter(key, value);
}

void IO::SetParameters(const Params &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetParameters");
    m_IO->SetParameters(parameters);
}

Params IO::Parameters() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::Parameters");
    return m_IO->GetParameters();
}

size_t IO::AddTransport(const std::string type, const Params &parameters)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AddTransport");
    return m_IO->AddTransport(type, parameters);
}

void IO::SetTransportParameter(const size_t transportIndex, const std::string key,
                               const std::string value)
{
    helper::CheckForNullptr(m_IO, "in call to IO::SetTransportParameter");
    m_IO->SetTransportParameter(transportIndex, key, value);
}

Variable IO::DefineVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::DefineVariable");
    return Variable(&m_IO->DefineVariable<std::string>(name));
}

Variable IO::DefineVariable(const std::string &name, const pybind11::array &array,
                            const Dims &shape, const Dims &start, const Dims &count,
                            const bool isConstantDims)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::DefineVariable");
    core::VariableBase *variable = nullptr;

    if (false)
    {
    }
#define declare_type(T)                                                                            \
    else if (pybind11::isinstance<pybind11::array_t<T, pybind11::array::c_style>>(array))          \
    {                                                                                              \
        variable = &m_IO->DefineVariable<T>(name, shape, start, count, isConstantDims);            \
    }
    ADIOS2_FOREACH_NUMPY_TYPE_1ARG(declare_type)
#undef declare_type
    else
    {
        throw std::invalid_argument("ERROR: variable " + name +
                                    " can't be defined, either type is not "
                                    "supported or is not memory "
                                    "contiguous, in call to DefineVariable\n");
    }

    return Variable(variable);
}

Variable IO::DefineVariable(const std::string &name, const pybind11::object &value,
                            const Dims &shape, const Dims &start, const Dims &count,
                            const bool isConstantDims)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::DefineVariable");
    core::VariableBase *variable = nullptr;
    const auto t = value.get_type();
    const auto ts = pybind11::str(t);
    const auto tss = pybind11::cast<std::string>(ts);
    if (pybind11::isinstance<pybind11::str>(value))
    {
        variable = &m_IO->DefineVariable<std::string>(name);
    }
    else if (pybind11::isinstance<pybind11::int_>(value))
    {
        variable = &m_IO->DefineVariable<int64_t>(name, shape, start, count, isConstantDims);
    }
    else if (pybind11::isinstance<pybind11::float_>(value))
    {
        variable = &m_IO->DefineVariable<double>(name, shape, start, count, isConstantDims);
    }
    else if (tss == "<class 'complex'>")
    {
        variable =
            &m_IO->DefineVariable<std::complex<double>>(name, shape, start, count, isConstantDims);
    }
    else if (tss == "<class 'numpy.complex64'>")
    {
        variable =
            &m_IO->DefineVariable<std::complex<float>>(name, shape, start, count, isConstantDims);
    }
    else if (tss == "<class 'numpy.complex128'>")
    {
        variable =
            &m_IO->DefineVariable<std::complex<double>>(name, shape, start, count, isConstantDims);
    }
    else
    {
        throw std::invalid_argument("ERROR: variable " + name +
                                    " can't be defined with an object with type " + tss +
                                    ", in call to DefineVariable\n");
    }
    return Variable(variable);
}

Variable IO::InquireVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::InquireVariable");

    const DataType type(m_IO->InquireVariableType(name));
    core::VariableBase *variable = nullptr;

    if (type == DataType::None)
    {
    }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        variable = m_IO->InquireVariable<T>(name);                                                 \
    }
    ADIOS2_FOREACH_PYTHON_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    return Variable(variable);
}

Attribute IO::DefineAttribute(const std::string &name, const pybind11::array &array,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    core::AttributeBase *attribute = nullptr;

    if (false)
    {
    }
#define declare_type(T)                                                                            \
    else if (pybind11::isinstance<pybind11::array_t<T, pybind11::array::c_style>>(array))          \
    {                                                                                              \
        const T *data = reinterpret_cast<const T *>(array.data());                                 \
        if (array.ndim())                                                                          \
        {                                                                                          \
            const size_t size = static_cast<size_t>(array.size());                                 \
            attribute = &m_IO->DefineAttribute<T>(name, data, size, variableName, separator);      \
        }                                                                                          \
        else                                                                                       \
        {                                                                                          \
            attribute = &m_IO->DefineAttribute<T>(name, data[0], variableName, separator);         \
        }                                                                                          \
    }
    ADIOS2_FOREACH_NUMPY_ATTRIBUTE_TYPE_1ARG(declare_type)
#undef declare_type
    else
    {
        throw std::invalid_argument("ERROR: attribute " + name +
                                    " can't be defined, either type is not "
                                    "supported or is not memory "
                                    "contiguous, in call to DefineAttribute\n");
    }

    return Attribute(attribute);
}

Attribute IO::DefineAttribute(const std::string &name, const std::string &stringValue,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    return Attribute(
        &m_IO->DefineAttribute<std::string>(name, stringValue, variableName, separator));
}

Attribute IO::DefineAttribute(const std::string &name, const std::vector<std::string> &strings,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    return Attribute(&m_IO->DefineAttribute<std::string>(name, strings.data(), strings.size(),
                                                         variableName, separator));
}

Attribute IO::DefineAttribute(const std::string &name, const std::vector<int> &ints,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    return Attribute(
        &m_IO->DefineAttribute<int>(name, ints.data(), ints.size(), variableName, separator));
}

Attribute IO::DefineAttribute(const std::string &name, const std::vector<double> &doubles,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    return Attribute(&m_IO->DefineAttribute<double>(name, doubles.data(), doubles.size(),
                                                    variableName, separator));
}

Attribute IO::DefineAttribute(const std::string &name,
                              const std::vector<std::complex<double>> &complexdoubles,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");
    return Attribute(&m_IO->DefineAttribute<std::complex<double>>(
        name, complexdoubles.data(), complexdoubles.size(), variableName, separator));
}

Attribute IO::DefineAttribute(const std::string &name, const pybind11::object &value,
                              const std::string &variableName, const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::DefineAttribute");

    core::AttributeBase *attribute = nullptr;
    const auto t = value.get_type();
    const auto ts = pybind11::str(t);
    const auto tss = pybind11::cast<std::string>(ts);
    if (pybind11::isinstance<pybind11::int_>(value))
    {
        auto v = pybind11::cast<const int64_t>(value);
        attribute = &m_IO->DefineAttribute(name, v, variableName, separator);
    }
    else if (pybind11::isinstance<pybind11::float_>(value))
    {
        auto v = pybind11::cast<const double>(value);
        attribute = &m_IO->DefineAttribute(name, v, variableName, separator);
    }
    else if (tss == "<class 'complex'>")
    {
        auto v = pybind11::cast<const std::complex<double>>(value);
        attribute = &m_IO->DefineAttribute(name, v, variableName, separator);
    }
    else if (tss == "<class 'numpy.complex64'>")
    {
        auto v = pybind11::cast<const std::complex<float>>(value);
        attribute = &m_IO->DefineAttribute(name, v, variableName, separator);
    }
    else if (tss == "<class 'numpy.complex128'>")
    {
        auto v = pybind11::cast<const std::complex<double>>(value);
        attribute = &m_IO->DefineAttribute(name, v, variableName, separator);
    }
    else
    {
        throw std::invalid_argument("ERROR: attribute " + name +
                                    " can't be defined with an object with type " + tss +
                                    ", in call to DefineAttribute\n");
    }
    return Attribute(attribute);
}

Attribute IO::InquireAttribute(const std::string &name, const std::string &variableName,
                               const std::string separator)
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + ", in call to IO::InquireAttribute");

    core::AttributeBase *attribute = nullptr;
    const DataType type(m_IO->InquireAttributeType(name, variableName, separator));

    if (type == DataType::None)
    {
    }
#define declare_template_instantiation(T)                                                          \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        attribute = m_IO->InquireAttribute<T>(name, variableName, separator);                      \
    }
    ADIOS2_FOREACH_ATTRIBUTE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

    return Attribute(attribute);
}

bool IO::RemoveVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::RemoveVariable");
    return m_IO->RemoveVariable(name);
}

void IO::RemoveAllVariables()
{
    helper::CheckForNullptr(m_IO, ", in call to IO::RemoveAllVariables");
    m_IO->RemoveAllVariables();
}

bool IO::RemoveAttribute(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable " + name + ", in call to IO::RemoveAttribute");
    return m_IO->RemoveAttribute(name);
}

void IO::RemoveAllAttributes()
{
    helper::CheckForNullptr(m_IO, ", in call to IO::RemoveAllAttributes");
    m_IO->RemoveAllAttributes();
}

Engine IO::Open(const std::string &name, const int mode)
{
    helper::CheckForNullptr(m_IO, "for engine " + name + ", in call to IO::Open");
    return Engine(&m_IO->Open(name, static_cast<adios2::Mode>(mode)));
}

void IO::FlushAll()
{
    helper::CheckForNullptr(m_IO, "in call to IO::FlushAll");
    m_IO->FlushAll();
}

std::map<std::string, Params> IO::AvailableVariables()
{
    helper::CheckForNullptr(m_IO, "in call to IO::AvailableVariables");
    return m_IO->GetAvailableVariables();
}

std::map<std::string, Params> IO::AvailableAttributes(const std::string &varname,
                                                      const std::string &separator)
{
    helper::CheckForNullptr(m_IO, "in call to IO::AvailableAttributes");
    return m_IO->GetAvailableAttributes(varname, separator);
}

std::string IO::VariableType(const std::string &name) const
{
    helper::CheckForNullptr(m_IO, "for variable " + name + " in call to IO::VariableType");
    return ToString(m_IO->InquireVariableType(name));
}

std::string IO::AttributeType(const std::string &name) const
{
    helper::CheckForNullptr(m_IO, "for attribute " + name + " in call to IO::AttributeType");
    return ToString(m_IO->InquireAttributeType(name));
}

std::string IO::EngineType() const
{
    helper::CheckForNullptr(m_IO, "in call to IO::EngineType");
    return m_IO->m_EngineType;
}

} // end namespace py11
} // end namespace adios2
