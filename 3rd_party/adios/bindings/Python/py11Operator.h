/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Operator.h :
 *
 *  Created on: Dec 12, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_OPERATOR_H_
#define ADIOS2_BINDINGS_PYTHON_OPERATOR_H_

#include <string>

#include "adios2/core/Operator.h"

namespace adios2
{
namespace py11
{

class ADIOS;
class IO;
class Variable;

class Operator
{
    friend class ADIOS;
    friend class IO;
    friend class Variable;

public:
    Operator() = default;

    ~Operator() = default;

    explicit operator bool() const noexcept;

    std::string Type() const noexcept;

    void SetParameter(const std::string key, const std::string value);

    Params &Parameters() const;

private:
    Params *m_Parameters;
    std::string m_Type;
    Operator(const std::string &type, Params *params);
};

} // end namespace py11
} // end namespace adios2

#endif /* BINDINGS_PYTHON_PY11OPERATOR_H_ */
