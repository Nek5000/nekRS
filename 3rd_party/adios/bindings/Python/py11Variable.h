/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Variable.h :
 *
 *  Created on: Sep 7, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_VARIABLE_H_
#define ADIOS2_BINDINGS_PYTHON_VARIABLE_H_

#include <pybind11/numpy.h>

#include "py11Operator.h"

#include "adios2/core/VariableBase.h"

namespace adios2
{
namespace py11
{

class IO;
class Engine;

class Variable
{
    friend class IO;
    friend class Engine;

public:
    Variable() = default;

    ~Variable() = default;

    explicit operator bool() const noexcept;

    void SetShape(const Dims &shape);

    void SetBlockSelection(const size_t blockID);

    void SetSelection(const Box<Dims> &selection);

    void SetStepSelection(const Box<size_t> &stepSelection);

    size_t SelectionSize() const;

    std::string Name() const;

    std::string Type() const;

    /**
     * Inspects size of the current element type, sizeof(T)
     * @return sizeof(T) for current system
     */
    size_t Sizeof() const;

    /**
     * Inspects shape id for current variable
     * @return from enum adios2::ShapeID
     */
    adios2::ShapeID ShapeID() const;

    /**
     * Inspects current shape
     * @return shape vector
     */
    Dims Shape(const size_t step = adios2::EngineCurrentStep) const;

    /**
     * Inspects current start point
     * @return start vector
     */
    Dims Start() const;

    /**
     * Inspects current count from start
     * @return count vector
     */
    Dims Count() const;

    /**
     * For read mode, inspect the number of available steps
     * @return available steps
     */
    size_t Steps() const;

    /**
     * For read mode, inspect the start step for available steps
     * @return available start step
     */
    size_t StepsStart() const;

    size_t BlockID() const;

    size_t SingleValue() const;

    /**
     * EXPERIMENTAL: Adds operation and parameters to current Variable object
     * @param op operator to be added
     * @param parameters key/value settings particular to the Variable, not to
     * be confused by op own parameters
     * @return operation index handler in Operations()
     */
    size_t AddOperation(const Operator op, const Params &parameters = Params());

    /**
     * EXPERIMENTAL: Adds operation and parameters to current Variable object
     * @param op operator to be added
     * @param parameters key/value settings particular to the Variable, not to
     * be confused by op own parameters
     * @return operation index handler in Operations()
     */
    size_t AddOperation(const std::string &op, const Params &parameters = Params());

    /**
     * EXPERIMENTAL: inspects current operators added with AddOperator
     * @return vector of Variable<T>::OperatorInfo
     */
    std::vector<Operator> Operations() const;

    /**
     * Removes all current Operations associated with AddOperation.
     * Provides the posibility to apply or not operators on a step basis.
     */
    void RemoveOperations();

    /** Contains sub-block information for a particular Variable<T> */
    struct Info
    {
        Dims Start;
        Dims Count;
        pybind11::array Min = pybind11::array();
        pybind11::array Max = pybind11::array();
        pybind11::array Value = pybind11::array();
        bool IsValue;
    };

private:
    Variable(core::VariableBase *variable);

    core::VariableBase *m_VariableBase = nullptr;
};

} // end namespace py11
} // end namespace adios2

#endif /* ADIOS2_BINDINGS_PYTHON_PY11VARIABLE_H_ */
