#ifndef ADIOS2_QUERY_INDEX_H
#define ADIOS2_QUERY_INDEX_H

#include "Query.h"

namespace adios2
{
namespace query
{
struct IndexInfo
{
    std::string m_IdxType; // minmax (default)
    adios2::Params m_SetupParameters;
};

template <class T>
class AbstractQueryIndex
{
public:
    void Generate(std::string &fromBPFile, const adios2::Params &inputs) = 0;
    void Evaluate(const QueryVar &query, std::vector<Box<Dims>> &touchedBlocks) = 0;
};

/*
template <class T>
  class BlockIndex : public AbstractQueryIndex<T>
{
public:
  void Generate(std::string& fromBPFile, const adios2::Params& inputs);
  void Evaluate(const Query<T>& query, std::vector<Box<Dims>>& touchedBlocks);
};
*/
}; // name space query
}; // name space adios2

#endif
