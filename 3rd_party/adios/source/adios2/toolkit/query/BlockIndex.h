#ifndef ADIOS2_BLOCK_INDEX_H
#define ADIOS2_BLOCK_INDEX_H

#include "Index.h"
#include "Query.h"

namespace adios2
{
namespace query
{

template <class T>
class BlockIndex
{
public:
    BlockIndex<T>(adios2::core::Variable<T> *var, adios2::core::IO &io,
                  adios2::core::Engine &reader)
    : m_VarPtr(var), m_IdxIO(io), m_IdxReader(reader)
    {
    }

    void Generate(std::string &fromBPFile, const adios2::Params &inputs) {}

    void Evaluate(const QueryVar &query, std::vector<BlockHit> &resultBlockIDs)
    {
        if (nullptr == m_VarPtr)
        {
            throw std::runtime_error("Unable to evaluate query! Invalid Variable detected");
        }

        size_t currStep = m_IdxReader.CurrentStep();
        adios2::Dims currShape = m_VarPtr->Shape();
        if (!query.IsSelectionValid(currShape))
            return;

        auto MinBlocksInfo = m_IdxReader.MinBlocksInfo(*m_VarPtr, currStep);

        if (MinBlocksInfo != nullptr)
        {
            RunStatMinBlocksInfo(query, MinBlocksInfo, resultBlockIDs);
            delete MinBlocksInfo;
        }
        else
        {
            RunStatBlocksInfo(query, currStep, resultBlockIDs);
        }
    }

    void RunStatMinBlocksInfo(const QueryVar &query, const adios2::MinVarInfo *MinBlocksInfo,
                              std::vector<BlockHit> &hitBlocks)
    {
        for (auto &blockInfo : MinBlocksInfo->BlocksInfo)
        {
            T bmin = *(T *)&blockInfo.MinMax.MinUnion;
            T bmax = *(T *)&blockInfo.MinMax.MaxUnion;
            bool isHit = query.m_RangeTree.CheckInterval(bmin, bmax);

            if (!isHit)
                continue;

            if (m_VarPtr->m_ShapeID != adios2::ShapeID::LocalArray)
            {
                Dims ss(MinBlocksInfo->Dims);
                Dims cc(MinBlocksInfo->Dims);
                for (std::vector<int>::size_type i = 0; i < ss.size(); i++)
                {
                    ss[i] = blockInfo.Start[i];
                    cc[i] = blockInfo.Count[i];
                }
                if (!query.TouchSelection(ss, cc))
                    continue;

                if (isHit)
                {
                    adios2::Box<adios2::Dims> box = {ss, cc};
                    hitBlocks.push_back(BlockHit(blockInfo.BlockID, box));
                }
            }
            else
            { // local array
                hitBlocks.push_back(BlockHit(blockInfo.BlockID));
            }
        }
    }

    void RunStatBlocksInfo(const QueryVar &query, const size_t currStep,
                           std::vector<BlockHit> &hitBlocks)
    {
        std::vector<typename adios2::core::Variable<T>::BPInfo> varBlocksInfo =
            m_IdxReader.BlocksInfo(*m_VarPtr, currStep);

        for (auto &blockInfo : varBlocksInfo)
        {
            bool isHit = query.m_RangeTree.CheckInterval(blockInfo.Min, blockInfo.Max);
            if (!isHit)
                continue;

            if (m_VarPtr->m_ShapeID == adios2::ShapeID::LocalArray)
            {
                if (isHit)
                    hitBlocks.push_back(BlockHit(blockInfo.BlockID));
            }
            else
            {
                // global array
                if (!query.TouchSelection(blockInfo.Start, blockInfo.Count))
                    continue;

                BlockHit tmp(blockInfo.BlockID);
                if (blockInfo.MinMaxs.size() > 0)
                {
                    // Consolidate to whole block If all subblocks are hits, then return the whole
                    // block
                    bool allCovered = true;

                    adios2::helper::CalculateSubblockInfo(blockInfo.Count, blockInfo.SubBlockInfo);
                    unsigned int numSubBlocks =
                        static_cast<unsigned int>(blockInfo.MinMaxs.size() / 2);
                    for (unsigned int i = 0; i < numSubBlocks; i++)
                    {
                        bool isSubblockHit = query.m_RangeTree.CheckInterval(
                            blockInfo.MinMaxs[2 * i], blockInfo.MinMaxs[2 * i + 1]);
                        if (isSubblockHit)
                        {
                            adios2::Box<adios2::Dims> currSubBlock = adios2::helper::GetSubBlock(
                                blockInfo.Count, blockInfo.SubBlockInfo, i);
                            for (size_t d = 0; d < blockInfo.Count.size(); ++d)
                                currSubBlock.first[d] += blockInfo.Start[d];

                            if (!query.TouchSelection(currSubBlock.first, currSubBlock.second))
                                continue;
                            tmp.m_Regions.push_back(currSubBlock);
                        }
                        else
                        {
                            allCovered = false;
                        }
                    } // for num subblocks

                    if (!allCovered)
                    {
                        hitBlocks.push_back(tmp);
                        continue;
                    }
                }

                // no subblock info or (allCovered = true)
                adios2::Box<adios2::Dims> box = {blockInfo.Start, blockInfo.Count};
                hitBlocks.push_back(BlockHit(blockInfo.BlockID, box));
            }
        }
    }

    // can not be unique_ptr as it changes with bp5 through steps
    // as BP5Deserializer::SetupForStep calls io.RemoveVariables()
    // must use ptr as bp5 associates ptrs with blockinfo, see MinBlocksInfo() in bp5
    adios2::core::Variable<T> *m_VarPtr;

private:
    //
    // blockid <=> vector of subcontents
    //
    // std::string& m_DataFileName;
    adios2::core::IO &m_IdxIO;
    adios2::core::Engine &m_IdxReader;

}; // class blockIndex

}; // end namespace query
}; // end namespace adios2

#endif
