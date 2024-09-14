/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * SstWriter.cpp
 *
 *  Created on: Aug 17, 2017
 *      Author: Greg Eisenhauer
 */

#include "adios2/helper/adiosComm.h"
#include <memory>

#include "SstParamParser.h"
#include "SstWriter.h"
#include "SstWriter.tcc"
#include "adios2/toolkit/format/buffer/malloc/MallocV.h"
#include <adios2-perfstubs-interface.h>

namespace adios2
{
namespace core
{
namespace engine
{

SstWriter::SstWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm)
: Engine("SstWriter", io, name, mode, std::move(comm))
{
    auto AssembleMetadata = [](void *writer, int CohortSize, struct _SstData * /*PerRankMetadata*/,
                               struct _SstData * /*PerRankAttributeData*/) {
        for (int i = 0; i < CohortSize; i++)
        {
            //            std::cout << "Rank " << i << " Metadata size is "
            //                      << PerRankMetadata[i].DataSize << ", base
            //                      address is "
            //                      << static_cast<void
            //                      *>(PerRankMetadata[i].block)
            //                      << std::endl;
        }
        /*
         *  This lambda function is provided to SST, and subsequently called
         *  from inside SstProvideTimestep on rank 0, but the
         *  PerRankMetadata and PerRankAttributeData entries correspond with
         *  the individual Metadata and AttributeData blocks that were
         *  passed to SstProvideTimestep by each writer rank.  That is, SST
         *  has brought those independent blocks to rank 0, and the purpose
         *  of this function is to assemble them into a single block of
         *  aggregate metadata.
         *
         *  This function should:
         *   - Take the PerRankMetadata blocks (which are just the Metadata
         * values that were passed to SstProvideTimestep) and assemble them into
         * a single serialized block
         *   - PerRankMetadata[0] should be modified to hold the address and
         * size of this single serialized block
         *   - If this single serialized block creates objects/data that must be
         * free'd later, return some value (cast to void*) that will enable you
         * to free it when it is passed to FreeAssembleMetadata.  (You'll also
         * get back PerRankMetadata[0], so if that's sufficient, you can return
         * NULL here.)
         *
         *  This function should *not* free any of the blocks pointer to by
         * PerRankMetadata.  It need not zero out or otherwise modify the
         * PerRankMetadata entries other than entry [0];

         *
         *  Note that this code works as is with current BP3 marshalling
         *  because it's already the case that only PerRanksMetadata[0] is
         *  valid (the rest are Null/0), and only 0 is used on the reading
         *  side.  This is the case because
         *  m_BP3Serializer->AggregateCollectiveMetadata() serializes the
         *  individual contributions of each rank, gathers them to rank 0,
         *  assembles them into a single block and then passes that to
         *  SstProvideTimestep().  However, this approach costs us extra MPI
         *  operations because AggregateCollectiveMetadata is doing them,
         *  and then we have more to do inside SstProvideTimestep().  To
         *  avoid these multiple MPI ops, we need to *not* use
         *  AggregateCollectiveMetadata.  Instead, each rank should
         *  serialize its individual contribution, and then pass that to
         *  ProvideTimestep.  Then this function should be changed to do the
         *  "assemble them into a single block" part.  Assuming that the
         *  result is the same as what AggregateCollectiveMetadata would
         *  have done, then the reader side should be able to remain
         *  unchanged.
         *
         *  Note also that this code is capable of handling marshaling
         *  systems that provide the attributes separately from the timestep
         *  metadata.  This is necessary in order to maintain attribute
         *  semantics when timesteps might be dropped on the writer side
         *  because of queue limits, or to be provided to late arriving
         *  readers, etc.  BP marshaling does not yet provide the ability to
         *  marshal attributes separately (AFAIK), so the AttributeData
         *  params can be ignored.
         */
        return (void *)malloc(1); /* return value will be passed as ClientData
                                     to registered Free routine */
    };

    auto FreeAssembledMetadata = [](void * /*writer*/, struct _SstData * /*PerRankMetadata*/,
                                    struct _SstData * /*PerRankAttributeData*/, void *ClientData) {
        //        std::cout << "Free called with client data " << ClientData
        //        << std::endl;
        free(ClientData);
        return;
    };

    Init();

    m_Output = SstWriterOpen(name.c_str(), &Params, &m_Comm);

    if (Params.MarshalMethod == SstMarshalBP)
    {
        SstWriterInitMetadataCallback(m_Output, this, AssembleMetadata, FreeAssembledMetadata);
    }
    m_IsOpen = true;
}

SstWriter::~SstWriter()
{
    if (m_IsOpen)
    {
        DestructorClose(m_FailVerbose);
    }
    SstStreamDestroy(m_Output);
}

StepStatus SstWriter::BeginStep(StepMode mode, const float timeout_sec)
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    m_WriterStep++;
    if (m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "SstWriter", "BeginStep",
                                        "BeginStep() is called a second time "
                                        "without an intervening EndStep()");
    }

    m_BetweenStepPairs = true;
    if (Params.MarshalMethod == SstMarshalFFS)
    {
        return (StepStatus)SstFFSWriterBeginStep(m_Output, (int)mode, timeout_sec);
    }
    else if (Params.MarshalMethod == SstMarshalBP)
    {
        // initialize BP serializer, deleted in
        // SstWriter::EndStep()::lf_FreeBlocks()
        m_BP3Serializer = std::unique_ptr<format::BP3Serializer>(new format::BP3Serializer(m_Comm));
        m_BP3Serializer->Init(m_IO.m_Parameters, "in call to BP3::Open for writing", "sst");
        m_BP3Serializer->ResizeBuffer(m_BP3Serializer->m_Parameters.InitialBufferSize,
                                      "in call to BP3::Open for writing by SST engine");
        m_BP3Serializer->m_MetadataSet.TimeStep = 1;
        m_BP3Serializer->m_MetadataSet.CurrentStep = m_WriterStep;
    }
    else if (Params.MarshalMethod == SstMarshalBP5)
    {
        if (!m_BP5Serializer)
        {
            m_BP5Serializer = std::unique_ptr<format::BP5Serializer>(new format::BP5Serializer());
            m_BP5Serializer->m_StatsLevel = Params.StatsLevel;
        }
        m_BP5Serializer->InitStep(new format::MallocV("SstWriter", true));
        m_BP5Serializer->m_Engine = this;
        //        m_BP5Serializer->Init(m_IO.m_Parameters,
        //                              "in call to BP5::Open for writing",
        //                              "sst");
        //        m_BP5Serializer->m_MetadataSet.TimeStep = 1;
        //        m_BP5Serializer->m_MetadataSet.CurrentStep = m_WriterStep;
    }
    else
    {
        // unknown marshaling method, shouldn't happen
    }
    return StepStatus::OK;
}

size_t SstWriter::CurrentStep() const { return m_WriterStep; }

void SstWriter::MarshalAttributes()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    const auto &attributes = m_IO.GetAttributes();

    if ((m_WriterStep == 0) && Params.UseOneTimeAttributes)
    {
        for (const auto &attributePair : attributes)
        {
            m_BP5Serializer->OnetimeMarshalAttribute(*(attributePair.second));
        }
    }

    // if there are no new attributes, nothing to do
    if (!m_MarshalAttributesNecessary)
        return;

    for (const auto &attributePair : attributes)
    {
        const std::string name(attributePair.first);
        const DataType type(attributePair.second->m_Type);

        if (type == DataType::None)
        {
        }
        else if (type == helper::GetDataType<std::string>())
        {
            core::Attribute<std::string> &attribute = *m_IO.InquireAttribute<std::string>(name);
            int element_count = -1;
            const char *data_addr = attribute.m_DataSingleValue.c_str();
            if (!attribute.m_IsSingleValue)
            {
                //
            }

            if (Params.MarshalMethod == SstMarshalFFS)
                SstFFSMarshalAttribute(m_Output, name.c_str(), (int)type, sizeof(char *),
                                       element_count, data_addr);
            else if (Params.MarshalMethod == SstMarshalBP5)
                m_BP5Serializer->MarshalAttribute(name.c_str(), type, sizeof(char *), element_count,
                                                  data_addr);
        }
#define declare_type(T)                                                                            \
    else if (type == helper::GetDataType<T>())                                                     \
    {                                                                                              \
        core::Attribute<T> &attribute = *m_IO.InquireAttribute<T>(name);                           \
        int element_count = -1;                                                                    \
        void *data_addr = &attribute.m_DataSingleValue;                                            \
        if (!attribute.m_IsSingleValue)                                                            \
        {                                                                                          \
            element_count = attribute.m_Elements;                                                  \
            data_addr = attribute.m_DataArray.data();                                              \
        }                                                                                          \
        if (Params.MarshalMethod == SstMarshalFFS)                                                 \
            SstFFSMarshalAttribute(m_Output, attribute.m_Name.c_str(), (int)type, sizeof(T),       \
                                   element_count, data_addr);                                      \
        else if (Params.MarshalMethod == SstMarshalBP5)                                            \
            m_BP5Serializer->MarshalAttribute(attribute.m_Name.c_str(), type, sizeof(T),           \
                                              element_count, data_addr);                           \
    }

        ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type
    }
    m_MarshalAttributesNecessary = false;
}

void SstWriter::NotifyEngineAttribute(std::string name, DataType type) noexcept
{
    helper::Throw<std::invalid_argument>("SstWriter", "Engine", "ThrowUp",
                                         "Engine does not support NotifyEngineAttribute");
}

void SstWriter::NotifyEngineAttribute(std::string name, AttributeBase *Attr, void *data) noexcept
{
    if (!Params.UseOneTimeAttributes)
    {
        m_MarshalAttributesNecessary = true;
        return;
    }

    m_BP5Serializer->OnetimeMarshalAttribute(*Attr);
    m_MarshalAttributesNecessary = false;
}

void SstWriter::EndStep()
{
    PERFSTUBS_SCOPED_TIMER_FUNC();
    if (!m_BetweenStepPairs)
    {
        helper::Throw<std::logic_error>("Engine", "SstWriter", "EndStep",
                                        "EndStep() is called without a successful BeginStep()");
    }
    m_BetweenStepPairs = false;
    if (m_WriterDefinitionsLocked && !m_DefinitionsNotified)
    {
        SstWriterDefinitionLock(m_Output, m_WriterStep);
        m_DefinitionsNotified = true;
    }
    if (Params.MarshalMethod == SstMarshalFFS)
    {
        PERFSTUBS_SCOPED_TIMER("Marshaling Overhead");
        PERFSTUBS_TIMER_START(timer, "SstMarshalFFS");
        MarshalAttributes();
        PERFSTUBS_TIMER_STOP(timer);
        SstFFSWriterEndStep(m_Output, m_WriterStep);
    }
    else if (Params.MarshalMethod == SstMarshalBP5)
    {
        MarshalAttributes();
        format::BP5Serializer::TimestepInfo *TSInfo =
            new format::BP5Serializer::TimestepInfo(m_BP5Serializer->CloseTimestep(m_WriterStep));
        auto lf_FreeBlocks = [](void *vBlock) {
            BP5DataBlock *BlockToFree = reinterpret_cast<BP5DataBlock *>(vBlock);
            //  Free data and metadata blocks here.  BlockToFree is the newblock
            //  value in the enclosing function.
            free(BlockToFree->MetaMetaBlocks);
            delete BlockToFree->TSInfo->DataBuffer;
            delete BlockToFree->TSInfo;
            delete BlockToFree;
        };

        BP5DataBlock *newblock = new BP5DataBlock;
        SstMetaMetaList MetaMetaBlocks = (SstMetaMetaList)malloc(
            (TSInfo->NewMetaMetaBlocks.size() + 1) * sizeof(MetaMetaBlocks[0]));
        int i = 0;
        for (const auto &MM : TSInfo->NewMetaMetaBlocks)
        {
            MetaMetaBlocks[i].BlockData = MM.MetaMetaInfo;
            MetaMetaBlocks[i].BlockSize = MM.MetaMetaInfoLen;
            MetaMetaBlocks[i].ID = MM.MetaMetaID;
            MetaMetaBlocks[i].IDSize = MM.MetaMetaIDLen;
            i++;
        }
        MetaMetaBlocks[TSInfo->NewMetaMetaBlocks.size()] = {NULL, 0, NULL, 0};
        newblock->MetaMetaBlocks = MetaMetaBlocks;
        newblock->metadata.DataSize = TSInfo->MetaEncodeBuffer->m_FixedSize;
        newblock->metadata.block = TSInfo->MetaEncodeBuffer->Data();
        std::vector<core::iovec> iovec = TSInfo->DataBuffer->DataVec();
        if (!iovec.empty())
        {
            newblock->data.DataSize = iovec[0].iov_len;
            newblock->data.block = (char *)iovec[0].iov_base;
        }
        else
        {
            newblock->data.DataSize = 0;
            newblock->data.block = nullptr;
        }
        newblock->TSInfo = TSInfo;
        if (TSInfo->AttributeEncodeBuffer)
        {
            newblock->attribute_data.DataSize = TSInfo->AttributeEncodeBuffer->m_FixedSize;
            newblock->attribute_data.block = TSInfo->AttributeEncodeBuffer->Data();
        }
        else
        {
            newblock->attribute_data.DataSize = 0;
            newblock->attribute_data.block = NULL;
        }
        SstProvideTimestepMM(m_Output, &newblock->metadata, &newblock->data, m_WriterStep,
                             lf_FreeBlocks, newblock, &newblock->attribute_data, NULL, newblock,
                             MetaMetaBlocks);
    }
    else if (Params.MarshalMethod == SstMarshalBP)
    {
        // This should finalize BP marshaling at the writer side.  All
        // marshaling methods should result in two blocks, one a block of
        // local metadata, and the second a block of data.  The metadata
        // will be aggregated across all writer ranks (by SST, inside
        // SstProvideTimestep) and made available to the readers (as a set
        // of N metadata blocks).  The Data block will be held on the writer
        // side (inside SstProvideTimestep), waiting on read requests from
        // individual reader ranks.  Any metadata or data blocks created
        // here should not be deallocated when SstProvideTimestep returns!
        // They should not be deallocated until SST is done with them
        // (explicit deallocation callback).
        PERFSTUBS_TIMER_START(timer, "Marshaling overhead");
        auto lf_FreeBlocks = [](void *vBlock) {
            BP3DataBlock *BlockToFree = reinterpret_cast<BP3DataBlock *>(vBlock);
            //  Free data and metadata blocks here.  BlockToFree is the newblock
            //  value in the enclosing function.
            delete BlockToFree->serializer;
            delete BlockToFree;
        };

        m_BP3Serializer->CloseStream(m_IO, true);
        m_BP3Serializer->AggregateCollectiveMetadata(m_Comm, m_BP3Serializer->m_Metadata, true);
        BP3DataBlock *newblock = new BP3DataBlock;
        newblock->metadata.DataSize = m_BP3Serializer->m_Metadata.m_Position;
        newblock->metadata.block = m_BP3Serializer->m_Metadata.m_Buffer.data();
        newblock->data.DataSize = m_BP3Serializer->m_Data.m_Position;
        newblock->data.block = m_BP3Serializer->m_Data.m_Buffer.data();
        newblock->serializer = m_BP3Serializer.release();
        PERFSTUBS_TIMER_STOP(timer);
        SstProvideTimestep(m_Output, &newblock->metadata, &newblock->data, m_WriterStep,
                           lf_FreeBlocks, newblock, NULL, NULL, NULL);
    }
    else
    {
        // unknown marshaling method, shouldn't happen
    }
}

void SstWriter::PerformPuts() {}

void SstWriter::Flush(const int transportIndex) {}

// PRIVATE functions below
void SstWriter::Init()
{
    SstParamParser Parser;

    Parser.ParseParams(m_IO, Params, m_UserOptions);

    if (Params.verbose < 0 || Params.verbose > 5)
    {
        helper::Throw<std::invalid_argument>("Engine", "SstWriter", "Init",
                                             "ERROR: Method verbose argument must be an "
                                             "integer in the range [0,5], in call to "
                                             "Open or Engine constructor\n");
    }
}

#define declare_type(T)                                                                            \
    void SstWriter::DoPutSync(Variable<T> &variable, const T *values)                              \
    {                                                                                              \
        PutSyncCommon(variable, values);                                                           \
    }                                                                                              \
    void SstWriter::DoPutDeferred(Variable<T> &variable, const T *values)                          \
    {                                                                                              \
        PutSyncCommon(variable, values);                                                           \
    }

ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

void SstWriter::PutStructCommon(VariableBase &variable, const void *data)
{
    size_t *Shape = NULL;
    size_t *Start = NULL;
    size_t *Count = NULL;
    size_t DimCount = 0;

    if (m_BetweenStepPairs == false)
    {
        helper::Throw<std::logic_error>("Engine", "SstWriter", "PutSyncCommon",
                                        "When using the SST engine in ADIOS2, "
                                        "Put() calls must appear between "
                                        "BeginStep/EndStep pairs");
    }

    if (Params.MarshalMethod != SstMarshalBP5)
    {
        helper::Throw<std::logic_error>(
            "Engine", "SstWriter", "PutStructCommon",
            "Support for struct types only exists when using BP5 marshalling");
    }

    if (variable.m_ShapeID == ShapeID::GlobalArray)
    {
        DimCount = variable.m_Shape.size();
        Shape = variable.m_Shape.data();
        Start = variable.m_Start.data();
        Count = variable.m_Count.data();
    }
    else if (variable.m_ShapeID == ShapeID::LocalArray)
    {
        DimCount = variable.m_Count.size();
        Count = variable.m_Count.data();
    }
    m_BP5Serializer->Marshal((void *)&variable, variable.m_Name.c_str(), variable.m_Type,
                             variable.m_ElementSize, DimCount, Shape, Count, Start, data, true,
                             nullptr);
}

void SstWriter::DoPutStructSync(VariableStruct &variable, const void *data)
{
    PutStructCommon(variable, data);
}

void SstWriter::DoPutStructDeferred(VariableStruct &variable, const void *data)
{
    PutStructCommon(variable, data);
}

void PutStruct(VariableStruct &, const void *, bool);

void SstWriter::DoClose(const int transportIndex) { SstWriterClose(m_Output); }

/**
 * Called if destructor is called on an open engine.  Should warn or take any
 * non-complex measure that might help recover.
 */
void SstWriter::DestructorClose(bool Verbose) noexcept
{
    if (Verbose)
    {
        std::cerr << "SST Writer \"" << m_Name << "\" Destroyed without a prior Close()."
                  << std::endl;
        std::cerr << "This may result in loss of data and/or disconnect "
                     "warnings for a connected SST Reader."
                  << std::endl;
    }
    m_IsOpen = false;
    // should at least call control plane to remove contact file
}

} // end namespace engine
} // end namespace core
} // end namespace adios2
