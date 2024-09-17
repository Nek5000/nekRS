/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DaosWriter.h
 *
 */

#ifndef ADIOS2_ENGINE_DAOS_DAOSWRITER_H_
#define ADIOS2_ENGINE_DAOS_DAOSWRITER_H_
#define DSS_PSETID "daos_server"

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/engine/daos/DaosEngine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosMemory.h" // PaddingToAlignOffset
#include "adios2/toolkit/aggregator/mpi/MPIChain.h"
#include "adios2/toolkit/aggregator/mpi/MPIShmChain.h"
#include "adios2/toolkit/burstbuffer/FileDrainerSingleThread.h"
#include "adios2/toolkit/format/bp5/BP5Serializer.h"
#include "adios2/toolkit/format/buffer/BufferV.h"
#include "adios2/toolkit/shm/Spinlock.h"
#include "adios2/toolkit/shm/TokenChain.h"
#include "adios2/toolkit/transportman/TransportMan.h"
#include <daos.h>
#include <daos_obj.h>
#include <mpi.h>

#define FAIL(fmt, ...)                                                                             \
    do                                                                                             \
    {                                                                                              \
        fprintf(stderr, "Process %d(%s): " fmt " aborting\n", m_Comm.Rank(), node, ##__VA_ARGS__); \
        MPI_Abort(MPI_COMM_WORLD, 1);                                                              \
    } while (0)
#define ASSERT(cond, ...)                                                                          \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
            FAIL(__VA_ARGS__);                                                                     \
    } while (0)

namespace adios2
{
namespace core
{
namespace engine
{

class DaosWriter : public DaosEngine, public core::Engine
{

public:
    /**
     * Constructor for file Writer in Daos format
     * @param name unique name given to the engine
     * @param openMode w (supported), r, a from OpenMode in ADIOSTypes.h
     * @param comm multi-process communicator
     */
    DaosWriter(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    ~DaosWriter();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const final;
    void PerformPuts() final;
    void PerformDataWrite() final;
    void EndStep() final;
    void Flush(const int transportIndex = -1) final;

    size_t DebugGetDataBufferSize() const final;

private:
    /** Single object controlling BP buffering */
    format::BP5Serializer m_BP5Serializer;

    /** Manage BP data files Transports from IO AddTransport */
    transportman::TransportMan m_FileDataManager;

    /** Manages the optional collective metadata files */
    transportman::TransportMan m_FileMetadataManager;

    /* transport manager for managing the metadata index file */
    transportman::TransportMan m_FileMetadataIndexManager;

    transportman::TransportMan m_FileMetaMetadataManager;

    /* DAOS declarations */

    uuid_t pool_uuid, cont_uuid;
    char *pool_label = "pool_ranjansv";
    char *cont_label = "adios-daos-engine-cont";

    /* Declare variables for pool and container handles */
    daos_handle_t poh, coh;

    enum DAOS_handleType
    {
        HANDLE_POOL,
        HANDLE_CO,
    };

    /* Declare variables for the KV object */
    daos_handle_t oh;
    daos_obj_id_t oid;

    char node[128] = "unknown";

    int64_t m_WriterStep = 0;
    /*
     *  Burst buffer variables
     */
    /** true if burst buffer is used to write */
    bool m_WriteToBB = false;
    /** true if burst buffer is drained to disk  */
    bool m_DrainBB = true;
    /** File drainer thread if burst buffer is used */
    burstbuffer::FileDrainerSingleThread m_FileDrainer;
    /** m_Name modified with burst buffer path if BB is used,
     * == m_Name otherwise.
     * m_Name is a constant of Engine and is the user provided target path
     */
    std::string m_BBName;
    /* Name of subfiles to directly write to (for all transports)
     * This is either original target or burst buffer if used */
    std::vector<std::string> m_SubStreamNames;
    /* Name of subfiles on target if burst buffer is used (for all transports)
     */
    std::vector<std::string> m_DrainSubStreamNames;
    std::vector<std::string> m_MetadataFileNames;
    std::vector<std::string> m_DrainMetadataFileNames;
    std::vector<std::string> m_MetaMetadataFileNames;
    std::vector<std::string> m_MetadataIndexFileNames;
    std::vector<std::string> m_DrainMetadataIndexFileNames;
    std::vector<std::string> m_ActiveFlagFileNames;

    bool m_BetweenStepPairs = false;

    void Init() final;

    /** Parses parameters from IO SetParameters */
    void InitParameters() final;
    /** Set up the aggregator */
    void InitAggregator();
    /** Complete opening/createing metadata and data files */
    void InitTransports() final;
    /** DAOS pool connection and container opening */
    void InitDAOS();
    /** Allocates memory and starts a PG group */
    void InitBPBuffer();
    void NotifyEngineAttribute(std::string name, DataType type) noexcept;
    /** Notify the engine when a new attribute is defined or modified. Called
     * from IO.tcc
     */
    void NotifyEngineAttribute(std::string name, AttributeBase *Attr, void *data) noexcept;

    void EnterComputationBlock() noexcept;
    /** Inform about computation block through User->ADIOS->IO */
    void ExitComputationBlock() noexcept;

#define declare_type(T)                                                                            \
    void DoPut(Variable<T> &variable, typename Variable<T>::Span &span, const bool initialize,     \
               const T &value) final;

    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

    template <class T>
    void PutCommonSpan(Variable<T> &variable, typename Variable<T>::Span &span,
                       const bool initialize, const T &value);

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &, const T *) final;                                                \
    void DoPutDeferred(Variable<T> &, const T *) final;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void PutCommon(VariableBase &variable, const void *data, bool sync);

#define declare_type(T, L)                                                                         \
    T *DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,                         \
                        const size_t bufferID = 0) noexcept final;

    ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

    void DoPutStructSync(VariableStruct &, const void *) final;
    void DoPutStructDeferred(VariableStruct &, const void *) final;

    void PutStruct(VariableStruct &, const void *, bool);

    void FlushData(const bool isFinal = false);

    void DoClose(const int transportIndex = -1) final;

    /** Write a profiling.json file from m_BP1Writer and m_TransportsManager
     * profilers*/
    void WriteProfilingJSONFile();

    void WriteMetaMetadata(const std::vector<format::BP5Base::MetaMetaInfoBlock> MetaMetaBlocks);

    void WriteMetadataFileIndex(uint64_t MetaDataPos, uint64_t MetaDataSize);

    uint64_t WriteMetadata(const std::vector<core::iovec> &MetaDataBlocks,
                           const std::vector<core::iovec> &AttributeBlocks);

    /** Write Data to disk, in an aggregator chain */
    void WriteData(format::BufferV *Data);
    void WriteData_EveryoneWrites(format::BufferV *Data, bool SerializedWriters);
    void WriteData_EveryoneWrites_Async(format::BufferV *Data, bool SerializedWriters);
    void WriteData_TwoLevelShm(format::BufferV *Data);
    void WriteData_TwoLevelShm_Async(format::BufferV *Data);

    void UpdateActiveFlag(const bool active);

    void WriteCollectiveMetadataFile(const bool isFinal = false);

    void MarshalAttributes();

    /* Two-level-shm aggregator functions */
    void WriteMyOwnData(format::BufferV *Data);
    void SendDataToAggregator(format::BufferV *Data);
    void WriteOthersData(const size_t TotalSize);

    template <class T>
    void PerformPutCommon(Variable<T> &variable);

    void FlushProfiler();

    /** manages all communication tasks in aggregation */
    aggregator::MPIAggregator *m_Aggregator; // points to one of these below
    aggregator::MPIShmChain m_AggregatorTwoLevelShm;
    aggregator::MPIChain m_AggregatorEveroneWrites;
    bool m_IAmDraining = false;
    bool m_IAmWritingData = false;
    helper::Comm *DataWritingComm; // processes that write the same data file
    // aggregators only (valid if m_Aggregator->m_Comm.Rank() == 0)
    helper::Comm m_CommAggregators;
    adios2::profiling::JSONProfiler m_Profiler;

protected:
    virtual void DestructorClose(bool Verbose) noexcept;

private:
    // updated during WriteMetaData
    uint64_t m_MetaDataPos = 0;

    /** On every process, at the end of writing, this holds the offset
     *  where they started writing (needed for global metadata)
     */
    uint64_t m_StartDataPos = 0;
    /** On aggregators, at the end of writing, this holds the starting offset
     *  to the next step's writing; otherwise used as temporary offset variable
     *  during writing on every process and points to the end of the process'
     *  data block in the file (not used for anything)
     */
    uint64_t m_DataPos = 0;

    /*
     *  Total data written this timestep
     */
    uint64_t m_ThisTimestepDataSize = 0;

    /** rank 0 collects m_StartDataPos in this vector for writing it
     *  to the index file
     */
    std::vector<uint64_t> m_WriterDataPos;

    bool m_MarshalAttributesNecessary = true;

    std::vector<std::vector<size_t>> FlushPosSizeInfo;

    void MakeHeader(std::vector<char> &buffer, size_t &position, const std::string fileType,
                    const bool isActive);

    std::vector<uint64_t> m_WriterSubfileMap; // rank => subfile index

    // Append helper data
    std::vector<size_t> m_AppendDataPos;  // each subfile append pos
    size_t m_AppendMetadataPos;           // metadata file append pos
    size_t m_AppendMetaMetadataPos;       // meta-metadata file append pos
    size_t m_AppendMetadataIndexPos;      // index file append pos
    uint32_t m_AppendWriterCount;         // last active number of writers
    unsigned int m_AppendAggregatorCount; // last active number of aggr
    unsigned int m_AppendSubfileCount;    // last active number of subfiles
    /* Process existing index, fill in append variables,
     * and return the actual step we land after appending.
     * Uses parameter AppendAfterStep
     * It resets m_Aggregator->m_NumAggregators so init aggregators later
     */
    uint64_t CountStepsInMetadataIndex(format::BufferSTL &bufferSTL);

    /* Async write's future */
    std::future<int> m_WriteFuture;
    // variables to delay writing to index file
    uint64_t m_LatestMetaDataPos;
    uint64_t m_LatestMetaDataSize;
    Seconds m_LastTimeBetweenSteps = Seconds(0.0);
    Seconds m_TotalTimeBetweenSteps = Seconds(0.0);
    Seconds m_AvgTimeBetweenSteps = Seconds(0.0);
    Seconds m_ExpectedTimeBetweenSteps = Seconds(0.0);
    TimePoint m_EndStepEnd;
    TimePoint m_EngineStart;
    TimePoint m_BeginStepStart;
    bool m_flagRush;                   // main thread flips this in Close, async thread watches it
    bool m_InComputationBlock = false; // main thread flips this in Clos
    TimePoint m_ComputationBlockStart;
    /* block counter and length in seconds */
    size_t m_ComputationBlockID = 0;

    struct ComputationBlockInfo
    {
        size_t blockID;
        double length; // seconds
        ComputationBlockInfo(const size_t id, const double len) : blockID(id), length(len){};
    };

    std::vector<ComputationBlockInfo> m_ComputationBlockTimes;
    /* sum of computationBlockTimes at start of async IO; */
    double m_ComputationBlocksLength = 0.0;

    /* struct of data passed from main thread to async write thread at launch */
    struct AsyncWriteInfo
    {
        adios2::aggregator::MPIAggregator *aggregator;
        int rank_global;
        helper::Comm comm_chain;
        int rank_chain;
        int nproc_chain;
        TimePoint tstart;
        adios2::shm::TokenChain<uint64_t> *tokenChain;
        transportman::TransportMan *tm;
        adios2::format::BufferV *Data;
        uint64_t startPos;
        uint64_t totalSize;
        double deadline;          // wall-clock time available in seconds
        bool *flagRush;           // flipped from false to true by main thread
        bool *inComputationBlock; // flipped back and forth by main thread
        // comm-free time within deadline in seconds
        double computationBlocksLength;
        std::vector<ComputationBlockInfo> expectedComputationBlocks; // a copy
        std::vector<ComputationBlockInfo> *currentComputationBlocks; // extended by main thread
        size_t *currentComputationBlockID;                           // increased by main thread
        shm::Spinlock *lock; // race condition over currentComp* variables
    };

    AsyncWriteInfo *m_AsyncWriteInfo;
    /* lock to handle race condition over the following currentComp* variables
         m_InComputationBlock / AsyncWriteInfo::inComputationBlock
         m_ComputationBlockID / AsyncWriteInfo::currentComputationBlockID
         m_flagRush / AsyncWriteInfo::flagRush
       Currently not used
         m_ComputationBlockTimes / AsyncWriteInfo::currentComputationBlocks
       Note: The rush flag does not need protection but CI TSAN sanitizer
       screams data race if not protected.
    */
    shm::Spinlock m_AsyncWriteLock;

    /* Static functions that will run in another thread */
    static int AsyncWriteThread_EveryoneWrites(AsyncWriteInfo *info);
    static int AsyncWriteThread_TwoLevelShm(AsyncWriteInfo *info);
    static void AsyncWriteThread_TwoLevelShm_Aggregator(AsyncWriteInfo *info);
    static void AsyncWriteThread_TwoLevelShm_SendDataToAggregator(aggregator::MPIShmChain *a,
                                                                  format::BufferV *Data);

    /* write own data used by both
       EveryoneWrites and TwoLevelShm  async threads  */
    static void AsyncWriteOwnData(AsyncWriteInfo *info, std::vector<core::iovec> &DataVec,
                                  const size_t totalsize, const bool seekOnFirstWrite);
    enum class ComputationStatus
    {
        InComp,
        NotInComp_ExpectMore,
        NoMoreComp
    };
    static ComputationStatus IsInComputationBlock(AsyncWriteInfo *info, size_t &compBlockIdx);

    void AsyncWriteDataCleanup();
    void AsyncWriteDataCleanup_EveryoneWrites();
    void AsyncWriteDataCleanup_TwoLevelShm();

    void daos_handle_share(daos_handle_t *, int);
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DAOS_DAOSWRITER_H_ */
