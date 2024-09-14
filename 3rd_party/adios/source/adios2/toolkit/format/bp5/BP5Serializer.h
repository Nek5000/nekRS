/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Serializer.h
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP5_BP5SERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_BP5_BP5SERIALIZER_H_

#include "BP5Base.h"
#include "adios2/core/Attribute.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/IO.h"
#include "adios2/toolkit/format/buffer/BufferV.h"
#include "adios2/toolkit/format/buffer/heap/BufferSTL.h"
#include "atl.h"
#include "ffs.h"
#include "fm.h"
#ifdef _WIN32
#pragma warning(disable : 4250)
#endif

#include <unordered_map>

namespace adios2
{
namespace format
{

class BP5Serializer : virtual public BP5Base
{

public:
    BP5Serializer();
    ~BP5Serializer();

    struct TimestepInfo
    {
        std::vector<MetaMetaInfoBlock> NewMetaMetaBlocks;
        std::shared_ptr<Buffer> MetaEncodeBuffer;
        std::shared_ptr<Buffer> AttributeEncodeBuffer;
        BufferV *DataBuffer;
    };

    typedef struct _MetadataInfo
    {
        std::vector<MetaMetaInfoBlock> NewMetaMetaBlocks;
        std::vector<size_t> MetaEncodeBufferSizes;
        std::vector<char *> MetaEncodeBuffers;

        std::vector<size_t> AttributeEncodeBufferSizes;
        std::vector<char *> AttributeEncodeBuffers;
        Buffer BackingBuffer;
    } AggregatedMetadataInfo;

    void Marshal(void *Variable, const char *Name, const DataType Type, size_t ElemSize,
                 size_t DimCount, const size_t *Shape, const size_t *Count, const size_t *Offsets,
                 const void *Data, bool Sync, BufferV::BufferPos *span);
    /*
     * BP5 has two attribute marshalling methods.  The first,
     * MarshallAttribute(), creates new MetaMetadata whenever a new
     * attribute gets marshalled, and produces AttributeData that
     * contains all extant attributes (so that only the most recent
     * need be installed to the Deserializer).  The second,
     * OnetimeMarshalAttribute(), produces MetaMetadata only on the
     * first timestep, and produces AttributeData that contains only
     * the attributes that were created or modified on that step (that
     * is, each timesteps attribute data only contains the delta from
     * the prior step).  Therefore for timestep X to have the right
     * attributes *all* attribute data created prior to timestep X
     * needs to be provided to the Deserializer (not just the most
     * recent, as with the other approach).  The first approach is
     * generally more space efficient and convenient for engines,
     * provided that attributes are few and new attributes are rare.
     * However, it performs very poorly if, for example, new
     * attributes are produced on every timestep.  In the latter case,
     * OnetimeMarshalAttribute() is by far better.
     */
    void MarshalAttribute(const char *Name, const DataType Type, size_t ElemSize, size_t ElemCount,
                          const void *Data);
    void OnetimeMarshalAttribute(const core::AttributeBase &baseAttr);
    void OnetimeMarshalAttribute(const char *Name, const DataType Type, size_t ElemCount,
                                 const void *Data);

    /*
     *  InitStep must be called with an appropriate BufferV subtype before a
     * step can begin
     */
    void InitStep(BufferV *DataBuffer);

    /*
     * ReinitStepData can be called to "flush" out already written
     * data it returns a BufferV representing already-written data and
     * provides the serializer with a new, empty BufferV This call
     * does *not* reset the data offsets generated with Marshal, so
     * those offsets are relative to the entire sequence of data
     * produced by a writer rank.
     */
    BufferV *ReinitStepData(BufferV *DataBuffer, bool forceCopyDeferred = false);

    TimestepInfo CloseTimestep(int timestep, bool forceCopyDeferred = false);
    void PerformPuts(bool forceCopyDeferred = false);

    /*
     * internal use  This calculates statistics on data that isn't available
     * until it's ready to be written to disk
     */
    void ProcessDeferredMinMax();

    core::Engine *m_Engine = NULL;

    std::vector<char>
    CopyMetadataToContiguous(const std::vector<BP5Base::MetaMetaInfoBlock> NewMetaMetaBlocks,
                             const std::vector<core::iovec> &MetaEncodeBuffers,
                             const std::vector<core::iovec> &AttributeEncodeBuffers,
                             const std::vector<uint64_t> &DataSizes,
                             const std::vector<uint64_t> &WriterDataPositions) const;

    std::vector<core::iovec>
    BreakoutContiguousMetadata(std::vector<char> &Aggregate, const std::vector<size_t> Counts,
                               std::vector<MetaMetaInfoBlock> &UniqueMetaMetaBlocks,
                               std::vector<core::iovec> &AttributeBlocks,
                               std::vector<uint64_t> &DataSizes,
                               std::vector<uint64_t> &WriterDataPositions) const;

    void *GetPtr(int bufferIdx, size_t posInBuffer);
    size_t CalcSize(const size_t Count, const size_t *Vals);

    size_t DebugGetDataBufferSize() const;

    MinVarInfo *MinBlocksInfo(const core::VariableBase &Var);

    int m_StatsLevel = 1;

    /* Variables to help appending to existing file */
    size_t m_PreMetaMetadataFileLength = 0;

    size_t m_BufferAlign = 1; // align buffers in memory
    // align buffers to integer multiples of block size
    size_t m_BufferBlockSize = sizeof(max_align_t);

private:
    void Init();
    typedef struct _BP5WriterRec
    {
        void *Key;
        int FieldID;
        ShapeID Shape;
        size_t DataOffset;
        size_t MetaOffset;
        char *OperatorType = NULL;
        int DimCount;
        int Type;
        size_t MinMaxOffset;
    } *BP5WriterRec;

    struct FFSWriterMarshalBase
    {
        int RecCount = 0;
        FMContext LocalFMContext = {0};
        int MetaFieldCount = 0;
        FMFieldList MetaFields = NULL;
        FMFormat MetaFormat;
        int AttributeFieldCount = 0;
        FMFieldList AttributeFields = NULL;
        FMFormat AttributeFormat = NULL;
        void *AttributeData = NULL;
        int AttributeSize = 0;
        std::unordered_map<std::string, _BP5WriterRec> RecNameMap;
    };

    FMFormat GenericAttributeFormat = NULL;
    std::vector<FMFormat> NewStructFormats;

    struct DeferredExtern
    {
        size_t MetaOffset;
        size_t BlockID;
        const void *Data;
        size_t DataSize;
        size_t AlignReq;
    };
    std::vector<DeferredExtern> DeferredExterns;

    struct DeferredSpanMinMax
    {
        const BufferV::BufferPos Data;
        const size_t ElemCount;
        const DataType Type;
        const MemorySpace MemSpace;
        const size_t MetaOffset;
        const size_t MinMaxOffset;
        const size_t BlockNum;
    };
    std::vector<DeferredSpanMinMax> DefSpanMinMax;

    BP5AttrStruct *PendingAttrs = nullptr;

    FFSWriterMarshalBase Info;
    void *MetadataBuf = NULL;
    bool NewAttribute = false;

    size_t MetadataSize = 0;
    BufferV *CurDataBuffer = NULL;

    std::vector<MetaMetaInfoBlock> PreviousMetaMetaInfoBlocks;

    size_t m_PriorDataBufferSizeTotal = 0;

    BP5WriterRec LookupWriterRec(void *Key) const;
    BP5WriterRec CreateWriterRec(void *Variable, const char *Name, DataType Type, size_t ElemSize,
                                 size_t DimCount);
    void ValidateWriterRec(BP5WriterRec Rec, void *Variable);
    void CollectFinalShapeValues();
    void RecalcMarshalStorageSize();
    void RecalcAttributeStorageSize();
    void AddSimpleField(FMFieldList *FieldP, int *CountP, const char *Name, const char *Type,
                        int ElementSize);
    void AddField(FMFieldList *FieldP, int *CountP, const char *Name, const DataType Type,
                  int ElementSize);
    void AddFixedArrayField(FMFieldList *FieldP, int *CountP, const char *Name, const DataType Type,
                            int ElementSize, int DimCount);
    void AddVarArrayField(FMFieldList *FieldP, int *CountP, const char *Name, const DataType Type,
                          int ElementSize, char *SizeField);
    void AddDoubleArrayField(FMFieldList *FieldP, int *CountP, const char *Name,
                             const DataType Type, int ElementSize, char *SizeField);
    char *BuildVarName(const char *base_name, const ShapeID Shape, const int type,
                       const int element_size);
    void BreakdownVarName(const char *Name, char **base_name_p, int *type_p, int *element_size_p);
    char *BuildArrayDimsName(const char *base_name, const int type, const int element_size);
    char *BuildArrayDBCountName(const char *base_name, const int type, const int element_size);
    char *BuildArrayBlockCountName(const char *base_name, const int type, const int element_size);
    char *TranslateADIOS2Type2FFS(const DataType Type);
    size_t *CopyDims(const size_t Count, const size_t *Vals);
    size_t *AppendDims(size_t *OldDims, const size_t OldCount, const size_t Count,
                       const size_t *Vals);

    void DumpDeferredBlocks(bool forceCopyDeferred = false);
    void VariableStatsEnabled(void *Variable);

    typedef struct _ArrayRec
    {
        size_t ElemCount;
        void *Array;
    } ArrayRec;

private:
    const void *SearchDeferredBlocks(size_t MetaOffset, size_t blocknum);
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_UTILITIES_FORMAT_B5_BP5Serializer_H_ */
