/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Operator.h :
 *
 *  Created on: Jun 7, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_CXX11_CXX11_OPERATOR_H_
#define ADIOS2_BINDINGS_CXX11_CXX11_OPERATOR_H_

#include "adios2/common/ADIOSMacros.h"
#include "adios2/common/ADIOSTypes.h"

namespace adios2
{

/// \cond EXCLUDE_FROM_DOXYGEN
// forward declare
class ADIOS; // friend class
class IO;    // friend class

template <class T>
class Variable; // friend class

namespace core
{
class Operator; // private implementation
}
/// \endcond

class Operator
{
public:
    /**
     * Empty (default) constructor, use it as a placeholder for future
     * operators from ADIOS::DefineOperator functions.
     * Can be used with STL containers.
     */
    Operator() = default;

    ~Operator() = default;

    /** true: valid object, false: invalid object */
    explicit operator bool() const noexcept;

    /**
     * Inspect current Operator type
     * @return type as string, if invalid returns an empty std::string
     */
    std::string Type() const noexcept;

    /**
     * Set a key/value parameters associated with this operator (global
     * parameter from the object it's applied to: Variable, IO).
     * If key exists, it replace the current value.
     * @param key parameter key
     * @param value parameter value
     */
    void SetParameter(const std::string key, const std::string value);

    /**
     * Inspect current operator parameters
     * @return map of key/value parameters
     */
    Params &Parameters() const;

private:
    friend class ADIOS;
    friend class IO;

    friend class VariableNT;
    template <class T>
    friend class Variable;

    Params *m_Parameters;
    std::string m_Type;
    Operator(const std::string &type, Params *params);
};

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_CXX11_OPERATOR_H_ */
