/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Variable.h :
 *
 *  Created on: Dec 11, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_ATTRIBUTE_H_
#define ADIOS2_BINDINGS_PYTHON_ATTRIBUTE_H_

#include <pybind11/numpy.h>

#include "adios2/core/AttributeBase.h"

namespace adios2
{
namespace py11
{

class IO;

class Attribute
{
    friend class IO;

public:
    Attribute() = default;

    ~Attribute() = default;

    explicit operator bool() const noexcept;

    std::string Name() const;

    std::string Type() const;

    bool SingleValue() const;

    pybind11::array Data();

    std::vector<std::string> DataString();

private:
    Attribute(core::AttributeBase *attribute);
    core::AttributeBase *m_Attribute = nullptr;
};

} // end namespace py11
} // end namespace adios2

#endif /* BINDINGS_PYTHON_PYATTRIBUTE_H_ */
