#include "Query.h"
#include "adios2/toolkit/query/Worker.h"

#include <utility>

namespace adios2
{
QueryWorker::QueryWorker(const std::string &configFile, adios2::Engine &reader)
{
    adios2::query::Worker *m = adios2::query::GetWorker(configFile, reader.m_Engine);
    if (m == nullptr)
        throw std::invalid_argument("ERROR: unable to construct query. ");
    m_Worker = std::make_shared<adios2::query::Worker>(std::move(*m));
    delete m;
}

void QueryWorker::GetResultCoverage(std::vector<size_t> &touched_blockIDs)
{
    m_Worker->GetResultCoverage(touched_blockIDs);
}

void QueryWorker::GetResultCoverage(std::vector<adios2::Box<adios2::Dims>> &touched_blocks)
{
    adios2::Box<adios2::Dims> empty;
    GetResultCoverage(empty, touched_blocks);
}

void QueryWorker::GetResultCoverage(const adios2::Box<adios2::Dims> &outputSelection,
                                    std::vector<adios2::Box<adios2::Dims>> &touched_blocks)
{
    if (m_Worker)
        return m_Worker->GetResultCoverage(outputSelection, touched_blocks);
}

} // namespace
