/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11IO.h
 *
 *  Created on: Mar 14, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_IO_H_
#define ADIOS2_BINDINGS_PYTHON_IO_H_

#include <pybind11/numpy.h>

#include <complex>
#include <string>

#include "py11Attribute.h"
#include "py11Engine.h"
#include "py11Variable.h"
#include "py11types.h"

namespace adios2
{
namespace py11
{

class IO
{
    friend class ADIOS;

public:
    IO() = default;

    ~IO() = default;

    explicit operator bool() const noexcept;

    bool InConfigFile() const;

    void SetEngine(const std::string type);

    void SetParameter(const std::string key, const std::string value);

    void SetParameters(const Params &parameters);

    Params Parameters() const;

    size_t AddTransport(const std::string type, const Params &parameters);

    void SetTransportParameter(const size_t transportIndex, const std::string key,
                               const std::string value);

    Variable DefineVariable(const std::string &name);

    Variable DefineVariable(const std::string &name, const pybind11::array &array,
                            const Dims &shape, const Dims &start, const Dims &count,
                            const bool isConstantDims);

    Variable DefineVariable(const std::string &name, const pybind11::object &value,
                            const Dims &shape, const Dims &start, const Dims &count,
                            const bool isConstantDims);

    Variable InquireVariable(const std::string &name);

    Attribute DefineAttribute(const std::string &name, const pybind11::array &array,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name, const std::string &stringValue,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name, const std::vector<std::string> &strings,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name, const std::vector<int> &ints,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name, const std::vector<double> &doubles,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name,
                              const std::vector<std::complex<double>> &complexdoubles,
                              const std::string &variableName = "",
                              const std::string separator = "/");

    Attribute DefineAttribute(const std::string &name, const pybind11::object &value,
                              const std::string &variableName, const std::string separator);

    Attribute InquireAttribute(const std::string &name, const std::string &variableName = "",
                               const std::string separator = "/");

    bool RemoveVariable(const std::string &name);

    void RemoveAllVariables();

    bool RemoveAttribute(const std::string &name);

    void RemoveAllAttributes();

    Engine Open(const std::string &name, const int openMode);

#if ADIOS2_USE_MPI
    Engine Open(const std::string &name, const int openMode, MPI4PY_Comm comm);
#endif

    void FlushAll();

    std::map<std::string, Params> AvailableVariables();

    std::map<std::string, Params> AvailableAttributes(const std::string &varname = "",
                                                      const std::string &separator = "/");

    std::string VariableType(const std::string &name) const;

    std::string AttributeType(const std::string &name) const;

    std::string EngineType() const;

private:
    IO(core::IO *io);
    core::IO *m_IO = nullptr;
};

} // end namespace py11
} // end namespace adios2

#endif
