/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Query.h : provides type utilities for ADIOS2 C++11 bindings
 *
 *  Created on: Aug 20, 2019
 *      Author: Junmin Gu <jgu@lbl.gov>
 */

#ifndef ADIOS2_BINDINGS_CXX11_QUERY_H_
#define ADIOS2_BINDINGS_CXX11_QUERY_H_

#include <memory> // otherwise MSVC complains about std::shared_ptr

#include "Engine.h"

namespace adios2
{

namespace query
{
class Worker;
} //

class QueryWorker
{
public:
    // configFile has query, can be either xml or json
    QueryWorker(const std::string &configFile, adios2::Engine &engine);

    void GetResultCoverage(std::vector<size_t> &touched_block_ids);
    // touched_blocks is a list of regions specified by (start, count),
    // that contains data that satisfies the query file
    void GetResultCoverage(std::vector<adios2::Box<adios2::Dims>> &touched_blocks);

    // supply output bound for the results
    void GetResultCoverage(const adios2::Box<adios2::Dims> &,
                           std::vector<adios2::Box<adios2::Dims>> &touched_blocks);

private:
    std::shared_ptr<adios2::query::Worker> m_Worker;
}; // class QueryWorker

} // end namespace adios2
#endif /* ADIOS2_BINDINGS_CXX11_QUERY_H_ */
