#include "Query.h"
#include "BlockIndex.h"
#include "adios2/helper/adiosFunctions.h"

#include "Query.tcc"

namespace adios2
{
namespace query
{

adios2::query::Relation strToRelation(std::string relationStr) noexcept
{
    if ((relationStr.compare("or") == 0) || (relationStr.compare("OR") == 0))
        return adios2::query::Relation::OR;

    return adios2::query::Relation::AND; // default
}

adios2::query::Op strToQueryOp(std::string opStr) noexcept
{
    if ((opStr.compare("lt") == 0) || (opStr.compare("LT") == 0))
        return adios2::query::Op::LT;
    if ((opStr.compare("gt") == 0) || (opStr.compare("GT") == 0))
        return adios2::query::Op::GT;
    if ((opStr.compare("ge") == 0) || (opStr.compare("GE") == 0))
        return adios2::query::Op::GE;
    if ((opStr.compare("le") == 0) || (opStr.compare("LE") == 0))
        return adios2::query::Op::LE;
    if ((opStr.compare("eq") == 0) || (opStr.compare("EQ") == 0))
        return adios2::query::Op::EQ;
    if ((opStr.compare("ne") == 0) || (opStr.compare("NE") == 0))
        return adios2::query::Op::NE;

    return adios2::query::Op::EQ; // default
}

adios2::Dims split(const std::string &s, char delim)
{
    adios2::Dims dim;

    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim))
    {
        std::stringstream curr(item);
        size_t val;
        curr >> val;
        dim.push_back(val);
    }

    return dim;
}

BlockHit::BlockHit(size_t id) : m_ID(id) {}

BlockHit::BlockHit(size_t id, Box<Dims> &box) : m_ID(id) { m_Regions.push_back(box); }

BlockHit::BlockHit(const BlockHit &cpy)
{
    m_ID = cpy.m_ID;
    m_Regions = cpy.m_Regions;
}

//
// return false if no intersection to apply
//  e.g. different blockID, no overlapped subblocks
// return true  if has intersection
//
// note that BlockHit comparison should match
// e.g. both are local/global arrays blocks
//
bool BlockHit::applyIntersection(const BlockHit &tmp)
{
    if (m_ID != tmp.m_ID)
        return false;

    // local array, no subblock info
    if (isLocalArrayBlock() || tmp.isLocalArrayBlock())
        return true;

    // check subblocks:
    bool overlapped = false;
    for (auto b : tmp.m_Regions)
    {
        for (auto it = m_Regions.begin(); it != m_Regions.end(); it++)
        {
            adios2::Box<Dims> curr = QueryBase::GetIntersection(*it, b);
            if (curr.first.size() != 0) // has intersection
            {
                overlapped = true;
                *it = curr;
            }
        }
    }

    return overlapped;
}

// return true  if can extended with tmp
//   e.g. at least same blockID with tmp
bool BlockHit::applyExtension(const BlockHit &tmp)
{
    if (m_ID != tmp.m_ID)
        return false;

    // check subblocks:
    for (auto b : tmp.m_Regions)
    {
        bool duplicated = false;
        for (auto box : m_Regions)
        {
            if (adios2::helper::IdenticalBoxes(box, b))
            {
                duplicated = true;
                continue;
            }
        }
        if (!duplicated)
            m_Regions.push_back(b);
    }

    return true;
}

void QueryBase::ApplyOutputRegion(std::vector<Box<Dims>> &touchedBlocks,
                                  const adios2::Box<Dims> &referenceRegion)
{
    if (m_OutputRegion.first.size() == 0)
        return;

    adios2::Dims diff;
    diff.resize(m_OutputRegion.first.size());
    bool isDifferent = false;
    for (size_t k = 0; k < m_OutputRegion.first.size(); k++)
    {
        diff[k] = m_OutputRegion.first[k] - referenceRegion.first[k];
        if (diff[k] != 0)
            isDifferent = true;
    }

    if (!isDifferent)
        return;

    // blocks are usually part of the reference region
    for (auto it = touchedBlocks.begin(); it != touchedBlocks.end(); it++)
    {
        for (size_t k = 0; k < m_OutputRegion.first.size(); k++)
            it->first[k] += diff[k];
    }
}
bool QueryComposite::AddNode(QueryBase *var)
{
    if (nullptr == var)
        return false;

    if (adios2::query::Relation::NOT == m_Relation)
    {
        // if (m_Nodes.size() > 0) return false;
        // don't want to support NOT for composite queries
        // return false;
        helper::Throw<std::ios_base::failure>("Toolkit", "query::QueryComposite", "AddNode",
                                              "Currently NOT is not suppprted for composite query");
    }
    m_Nodes.push_back(var);
    return true;
}

void QueryComposite::BlockIndexEvaluate(adios2::core::IO &io, adios2::core::Engine &reader,
                                        // std::vector<Box<Dims>> &touchedBlocks)
                                        std::vector<BlockHit> &touchedBlocks)
{
    auto lf_ApplyAND = [&](std::vector<BlockHit> &touched,
                           const std::vector<BlockHit> &curr) -> void {
        if (curr.size() == 0)
        {
            touched.clear();
            return;
        }

        for (auto i = touched.size(); i >= 1; i--)
        {
            bool intersects = false;
            for (auto b : curr)
            {
                if (touched[i].applyIntersection(b))
                {
                    intersects = true;
                    break;
                }
            }
            if (!intersects)
                touched.erase(touched.begin() + i - 1);
        }
    }; // lf_ApplyAND

    auto lf_ApplyOR = [&](std::vector<BlockHit> &touched,
                          const std::vector<BlockHit> &curr) -> void {
        if (curr.size() == 0)
            return;

        for (auto b : curr)
        {
            bool duplicated = false;
            for (auto box : touched)
            {
                if (box.applyExtension(b))
                {
                    duplicated = true;
                    continue;
                }
            }
            if (!duplicated)
                touched.push_back(b);
        }
    }; // lf_ApplyOR

    if (m_Nodes.size() == 0)
        return;

    int counter = 0;
    for (auto node : m_Nodes)
    {
        counter++;
        std::vector<BlockHit> currBlocks;
        node->BlockIndexEvaluate(io, reader, currBlocks);
        if (counter == 1)
        {
            touchedBlocks = currBlocks;
            continue;
        }

        if (currBlocks.size() == 0)
        {
            if (adios2::query::Relation::AND == m_Relation)
            {
                touchedBlocks.clear();
                break;
            }
            else
                continue;
        }

        if (adios2::query::Relation::AND == m_Relation)
            lf_ApplyAND(touchedBlocks, currBlocks);
        else if (adios2::query::Relation::OR == m_Relation)
            lf_ApplyOR(touchedBlocks, currBlocks);
    }
}

bool QueryVar::IsSelectionValid(adios2::Dims &shape) const
{
    if (0 == m_Selection.first.size())
        return true;

    if (m_Selection.first.size() != shape.size())
    {
        helper::Log("Query", "QueryVar", "IsSelectionValid",
                    "Query selection dimension is different from shape dimension",
                    helper::FATALERROR);
        return false; // different dimension
    }

    return true;
}

void QueryVar::LoadSelection(const std::string &startStr, const std::string &countStr)
{
    adios2::Dims start = split(startStr, ',');
    adios2::Dims count = split(countStr, ',');

    if (start.size() != count.size())
    {
        helper::Throw<std::ios_base::failure>(
            "Toolkit", "query::QueryVar", "LoadSelection",
            "dim of startcount does match in the bounding box definition");
    }

    // simpleQ.setSelection(box.first, box.second);
    adios2::Dims shape = this->m_Selection.second; // set at the creation for default
    this->SetSelection(start, count);
    if (!this->IsSelectionValid(shape))
        helper::Throw<std::ios_base::failure>("Toolkit", "query::QueryVar", "LoadSelection",
                                              "invalid selections for selection of var: " +
                                                  this->GetVarName());
}

bool QueryVar::TouchSelection(adios2::Dims &start, adios2::Dims &count) const
{
    if (0 == m_Selection.first.size())
        return true;

    const size_t dimensionsSize = start.size();

    for (size_t i = 0; i < dimensionsSize; i++)
    {
        size_t end = start[i] + count[i];
        size_t selEnd = m_Selection.first[i] + m_Selection.second[i];

        if (end <= m_Selection.first[i])
            return false;
        if (selEnd <= start[i])
            return false;
    }
    return true;
}

void QueryVar::BlockIndexEvaluate(adios2::core::IO &io, adios2::core::Engine &reader,
                                  std::vector<BlockHit> &touchedBlocks)
{
    const DataType varType = io.InquireVariableType(m_VarName);

    // var already exists when loading query. skipping validity checking
#define declare_type(T)                                                                            \
    if (varType == adios2::helper::GetDataType<T>())                                               \
    {                                                                                              \
        core::Variable<T> *var = io.InquireVariable<T>(m_VarName);                                 \
        BlockIndex<T> idx(var, io, reader);                                                        \
        idx.Evaluate(*this, touchedBlocks);                                                        \
    }
    // ADIOS2_FOREACH_ATTRIBUTE_TYPE_1ARG(declare_type) //skip complex types
    ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

    if (touchedBlocks.size() > 0)
    {
        LimitToSelection(touchedBlocks);

        for (auto blk : touchedBlocks)
        {
            if (!blk.isLocalArrayBlock())
            {
                ApplyOutputRegion(blk.m_Regions, m_Selection);
            }
        }
    }
}
} // namespace query
} // namespace adios2
