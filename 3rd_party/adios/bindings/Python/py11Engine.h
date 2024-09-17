/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * EnginePy.h
 *
 *  Created on: Mar 15, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_PYTHON_ENGINE_H_
#define ADIOS2_BINDINGS_PYTHON_ENGINE_H_

#include <pybind11/numpy.h>

#include <string>

#include "adios2/core/Engine.h"

#include "py11Variable.h"

namespace adios2
{
namespace py11
{

// forward declare
class IO; // friend

class Engine
{
    friend class IO;
    friend class Query;

public:
    struct Info
    {
        Dims Start;
        Dims Count;
    };

    Engine() = default;

    ~Engine() = default;

    explicit operator bool() const noexcept;

    StepStatus BeginStep(const StepMode mode, const float timeoutSeconds = -1.f);
    StepStatus BeginStep();

    void Put(Variable variable, const pybind11::array &array, const Mode launch = Mode::Deferred);
    void Put(Variable variable, const std::vector<int64_t> &ints,
             const Mode launch = Mode::Deferred);
    void Put(Variable variable, const std::vector<double> &doubles,
             const Mode launch = Mode::Deferred);
    void Put(Variable variable, const std::vector<std::complex<double>> &complexes,
             const Mode launch = Mode::Deferred);
    void Put(Variable variable, const std::string &string);
    void PerformPuts();
    void PerformDataWrite();

    void Get(Variable variable, pybind11::array &array, const Mode launch = Mode::Deferred);
    std::string Get(Variable variable, const Mode launch = Mode::Deferred);

    void PerformGets();

    void EndStep();

    /**
     * Returns current status information for each engine.
     * @return if between BeginStep/EndStep() pair
     */
    bool BetweenStepPairs() const;

    void Flush(const int transportIndex = -1);

    void Close(const int transportIndex = -1);

    size_t CurrentStep() const;

    std::string Name() const;
    std::string Type() const;
    size_t Steps() const;
    void LockWriterDefinitions() const;
    void LockReaderSelections() const;

    std::vector<std::map<std::string, std::string>> BlocksInfo(std::string &string,
                                                               const size_t step) const;

private:
    Engine(core::Engine *engine);
    core::Engine *m_Engine = nullptr;
};

} // end namespace py11
} // end namespace adios2

#endif
