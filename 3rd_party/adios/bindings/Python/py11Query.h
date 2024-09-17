/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Query.h :
 *
 *  Created on: August 5, 2021
 *      Author: Junmin Gu (jgu@lbl.gov)
 */

#ifndef ADIOS2_BINDINGS_PYTHON_QUERY_H_
#define ADIOS2_BINDINGS_PYTHON_QUERY_H_

#include <pybind11/numpy.h>

// #include "adios2/toolkit/query/Query.h"
#include "adios2/toolkit/query/Worker.h"
#include "py11Engine.h"

namespace adios2
{
namespace py11
{

class Engine;

class Query
{

public:
    Query(std::string queryFile, Engine reader);
    ~Query() = default;

    explicit operator bool() const noexcept;

    std::vector<Box<Dims>> GetResult();
    std::vector<size_t> GetBlockIDs();

private:
    Query(adios2::query::Worker *qw);
    std::shared_ptr<adios2::query::Worker> m_QueryWorker;
};

} // end namespace py11
} // end namespace adios2

#endif /* BINDINGS_PYTHON_PYQUERY_H_ */
