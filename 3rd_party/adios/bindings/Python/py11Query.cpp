/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * py11Query.h :
 *
 *  Created on: August 5, 2021
 *      Author: Junmin Gu (jgu@lbl.gov)
 */

#include "py11Query.h"

namespace adios2
{
namespace py11
{

Query::Query(adios2::query::Worker *qw) : m_QueryWorker(qw) {}

Query::Query(std::string queryFile, Engine reader)
{
    adios2::query::Worker *m = adios2::query::GetWorker(queryFile, reader.m_Engine);
    if (m == nullptr)
        throw std::invalid_argument("ERROR: unable to construct query. ");
    m_QueryWorker = std::make_shared<adios2::query::Worker>(std::move(*m));
    delete m;
}

Query::operator bool() const noexcept { return (m_QueryWorker == nullptr) ? false : true; }

std::vector<Box<Dims>> Query::GetResult()
{
    std::vector<Box<Dims>> touched_blocks;
    adios2::Box<adios2::Dims> empty;
    m_QueryWorker->GetResultCoverage(empty, touched_blocks);
    return touched_blocks;
}

std::vector<size_t> Query::GetBlockIDs()
{
    std::vector<size_t> touched_block_ids;
    m_QueryWorker->GetResultCoverage(touched_block_ids);
    return touched_block_ids;
}

} // py11
} // adios2
