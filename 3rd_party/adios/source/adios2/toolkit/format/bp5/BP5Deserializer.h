/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP5Deserializer.h
 *
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP5_BP5DESERIALIZER_H_
#define ADIOS2_TOOLKIT_FORMAT_BP5_BP5DESERIALIZER_H_

#include "adios2/core/Attribute.h"
#include "adios2/core/IO.h"
#include "adios2/core/Variable.h"

#include "BP5Base.h"
#include "atl.h"
#include "ffs.h"
#include "fm.h"

#include <mutex>

#ifdef _WIN32
#pragma warning(disable : 4250)
#endif

namespace adios2
{
namespace format
{

using namespace core;

class BP5Deserializer : virtual public BP5Base
{

public:
    BP5Deserializer(bool WriterIsRowMajor, bool ReaderIsRowMajor, bool RandomAccessMode = false);
    BP5Deserializer(bool WriterIsRowMajor, bool ReaderIsRowMajor, bool RandomAccessMode,
                    bool FlattenSteps);

    ~BP5Deserializer();

    struct ReadRequest
    {
        size_t Timestep;
        size_t WriterRank;
        size_t StartOffset;
        size_t ReadLength;
        char *DestinationAddr;
        bool DirectToAppMemory;
        size_t ReqIndex;
        size_t OffsetInBlock;
        size_t BlockID;
    };
    void InstallMetaMetaData(MetaMetaInfoBlock &MMList);
    void InstallMetaData(void *MetadataBlock, size_t BlockLen, size_t WriterRank,
                         size_t Step = SIZE_MAX);
    void InstallAttributeData(void *AttributeBlock, size_t BlockLen, size_t Step = SIZE_MAX);
    void InstallAttributesV1(FFSTypeHandle FFSformat, void *BaseData, size_t Step);
    void InstallAttributesV2(FFSTypeHandle FFSformat, void *BaseData, size_t Step);

    void SetupForStep(size_t Step, size_t WriterCount);
    // return from QueueGet is true if a sync is needed to fill the data
    bool QueueGet(core::VariableBase &variable, void *DestData);
    bool QueueGetSingle(core::VariableBase &variable, void *DestData, size_t AbsStep,
                        size_t RelStep);

    /* generate read requests. return vector of requests AND the size of
     * the largest allocation block necessary for reading.
     * input flag: true allocates a temporary buffer for each read request
     * unless the request can go directly to user memory.
     * False will not allocate a temporary buffer
     * (RR.DestinationAddress==nullptr) but may also assign the user memory
     * pointer for direct read in
     */
    std::vector<ReadRequest> GenerateReadRequests(const bool doAllocTempBuffers,
                                                  size_t *maxReadSize);
    void FinalizeGet(const ReadRequest &, const bool freeAddr);
    void FinalizeGets(std::vector<ReadRequest> &);

    MinVarInfo *AllRelativeStepsMinBlocksInfo(const VariableBase &var);
    MinVarInfo *AllStepsMinBlocksInfo(const VariableBase &var);
    MinVarInfo *MinBlocksInfo(const VariableBase &Var, const size_t Step);
    bool VarShape(const VariableBase &, const size_t Step, Dims &Shape) const;
    bool VariableMinMax(const VariableBase &var, const size_t Step, MinMaxStruct &MinMax);
    void GetAbsoluteSteps(const VariableBase &variable, std::vector<size_t> &keys) const;

    const bool m_WriterIsRowMajor;
    const bool m_ReaderIsRowMajor;
    core::Engine *m_Engine = NULL;

    enum RequestTypeEnum
    {
        Global = 0,
        Local = 1
    };

    struct BP5ArrayRequest
    {
        void *VarRec = NULL;
        char *VarName;
        enum RequestTypeEnum RequestType;
        size_t Step;    // local operations use absolute steps
        size_t RelStep; // preserve Relative Step for remote
        size_t BlockID;
        Dims Start;
        Dims Count;
        MemorySpace MemSpace;
        void *Data;
    };
    std::vector<BP5ArrayRequest> PendingGetRequests;

private:
    size_t m_VarCount = 0;
    struct BP5VarRec
    {
        size_t VarNum;
        void *Variable = NULL;
        char *VarName = NULL;
        size_t DimCount = 0;
        size_t JoinedDimen = SIZE_MAX;
        size_t *LastJoinedOffset = NULL;
        size_t *LastJoinedShape = NULL;
        bool Derived = false;
        char *ExprStr = NULL;
        ShapeID OrigShapeID;
        core::StructDefinition *Def = nullptr;
        core::StructDefinition *ReaderDef = nullptr;
        char *Operator = NULL;
        DataType Type;
        int ElementSize = 0;
        size_t MinMaxOffset = SIZE_MAX;
        size_t *GlobalDims = NULL;
        size_t LastTSAdded = SIZE_MAX;
        size_t FirstTSSeen = SIZE_MAX;
        size_t LastStepAdded = SIZE_MAX;
        std::vector<size_t> AbsStepFromRel; // per relative step vector
        std::vector<size_t> PerWriterMetaFieldOffset;
        std::vector<size_t> PerWriterBlockStart;
    };

    struct ControlStruct
    {
        int FieldOffset;
        BP5VarRec *VarRec;
        ShapeID OrigShapeID;
        DataType Type;
        int ElementSize;
    };

    struct ControlInfo
    {
        FMFormat Format;
        int ControlCount;
        struct ControlInfo *Next;
        std::vector<size_t> *MetaFieldOffset;
        std::vector<size_t> *CIVarIndex;
        struct ControlStruct Controls[1];
    };

    enum WriterDataStatusEnum
    {
        Empty = 0,
        Needed = 1,
        Requested = 2,
        Full = 3
    };

    struct FFSReaderPerWriterRec
    {
        enum WriterDataStatusEnum Status = Empty;
    };

    FFSContext ReaderFFSContext;

    const bool m_RandomAccessMode;
    const bool m_FlattenSteps;

    std::vector<size_t> m_WriterCohortSize; // per step, in random mode
    size_t m_CurrentWriterCohortSize;       // valid in streaming mode
    // return the number of writers
    // m_CurrentWriterCohortSize in streaming mode
    // m_WriterCohortSize[Step] in random access mode
    size_t WriterCohortSize(size_t Step) const;

    size_t m_LastAttrStep = MaxSizeT; // invalid timestep for start

    std::unordered_map<std::string, BP5VarRec *> VarByName;
    std::unordered_map<void *, BP5VarRec *> VarByKey;

    std::vector<void *> *m_MetadataBaseAddrs =
        nullptr; // may be a pointer into MetadataBaseArray or m_FreeableMBA
    std::vector<void *> *m_FreeableMBA = nullptr;

    // for random access mode, for each timestep, for each writerrank, what
    // metameta info applies to the metadata
    std::vector<std::vector<ControlInfo *>> m_ControlArray;
    // for random access mode, for each timestep, for each writerrank, base
    // address of the metadata
    std::vector<std::vector<void *> *> MetadataBaseArray;
    // for random access mode, for each timestep, for each writerrank, base
    // address of the joined dim arrays, for streaming use 0 index
    std::vector<std::vector<size_t *>> JoinedDimArray;
    size_t JDAIdx = 0;

    ControlInfo *ControlBlocks = nullptr;
    ControlInfo *GetPriorControl(FMFormat Format);
    ControlInfo *BuildControl(FMFormat Format);
    bool NameIndicatesArray(const char *Name);
    bool NameIndicatesAttrArray(const char *Name);
    DataType TranslateFFSType2ADIOS(const char *Type, int size);
    BP5VarRec *LookupVarByKey(void *Key) const;
    BP5VarRec *LookupVarByName(const char *Name);
    BP5VarRec *CreateVarRec(const char *ArrayName);
    void ReverseDimensions(size_t *Dimensions, size_t count, size_t times);
    const char *BreakdownVarName(const char *Name, DataType *type_p, int *element_size_p);
    void BreakdownFieldType(const char *FieldType, bool &Operator, bool &MinMax);
    void BreakdownArrayName(const char *Name, char **base_name_p, DataType *type_p,
                            int *element_size_p, FMFormat *Format);
    void BreakdownV1ArrayName(const char *Name, char **base_name_p, DataType *type_p,
                              int *element_size_p, bool &Operator, bool &MinMax);
    void *VarSetup(core::Engine *engine, const char *variableName, const DataType type, void *data);
    void *ArrayVarSetup(core::Engine *engine, const char *variableName, const DataType type,
                        int DimCount, size_t *Shape, size_t *Start, size_t *Count,
                        core::StructDefinition *Def, core::StructDefinition *ReaderDef);
    void MapGlobalToLocalIndex(size_t Dims, const size_t *GlobalIndex, const size_t *LocalOffsets,
                               size_t *LocalIndex);
    size_t RelativeToAbsoluteStep(const BP5VarRec *VarRec, size_t RelStep);
    int FindOffset(size_t Dims, const size_t *Size, const size_t *Index);
    bool GetSingleValueFromMetadata(core::VariableBase &variable, BP5VarRec *VarRec, void *DestData,
                                    size_t Step, size_t WriterRank);
    void StructQueueReadChecks(core::VariableStruct *variable, BP5VarRec *VarRec);

    void *GetMetadataBase(BP5VarRec *VarRec, size_t Step, size_t WriterRank) const;
    bool IsContiguousTransfer(BP5ArrayRequest *Req, size_t *offsets, size_t *count);

    size_t CurTimestep = 0;

    /* We assume operators are not thread-safe, call Decompress() one at a time
     */
    std::mutex mutexDecompress;
};

} // end namespace format
} // end namespace adios2

#endif /* ADIOS2_UTILITIES_FORMAT_BP5_BP5Serializer_H_ */
