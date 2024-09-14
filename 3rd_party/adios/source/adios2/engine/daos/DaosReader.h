/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DaosReader.h
 *
 */

#ifndef ADIOS2_ENGINE_DAOS_DAOSREADER_H_
#define ADIOS2_ENGINE_DAOS_DAOSREADER_H_
#define DSS_PSETID "daos_server"

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/engine/daos/DaosEngine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/helper/adiosRangeFilter.h"
#include "adios2/toolkit/format/bp5/BP5Deserializer.h"
#include "adios2/toolkit/transportman/TransportMan.h"

#include <chrono>
#include <daos.h>
#include <daos_obj.h>
#include <map>
#include <mpi.h>
#include <vector>

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

class DaosReader : public DaosEngine, public Engine
{

public:
    /**
     * Unique constructor
     * @param io
     * @param name
     * @param openMode only read
     * @param comm
     */
    DaosReader(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    ~DaosReader();

    StepStatus BeginStep(StepMode mode = StepMode::Read, const float timeoutSeconds = -1.0) final;

    size_t CurrentStep() const final;

    void EndStep() final;

    void PerformGets() final;

    MinVarInfo *MinBlocksInfo(const VariableBase &, const size_t Step) const;
    bool VarShape(const VariableBase &Var, const size_t Step, Dims &Shape) const;
    bool VariableMinMax(const VariableBase &, const size_t Step, MinMaxStruct &MinMax);

private:
    format::BP5Deserializer *m_BP5Deserializer = nullptr;
    /* transport manager for metadata file */
    transportman::TransportMan m_MDFileManager;
    /* How many bytes of metadata have we already read in? */
    size_t m_MDFileAlreadyReadSize = 0;
    /* How many bytes of metadata have we already processed?
     * It is <= m_MDFileAlreadyReadSize, at = we need to read more */
    size_t m_MDFileProcessedSize = 0;
    /* The file position of the first byte that is currently
     * residing in memory. Needed for skewing positions when
     * processing metadata index.
     */
    size_t m_MDFileAbsolutePos = 0;
    /* m_MDFileAbsolutePos <= m_MDFileProcessedSize <= m_MDFileAlreadyReadSize
     */

    /* transport manager for managing data file(s) */
    transportman::TransportMan m_DataFileManager;

    /* transport manager for managing the metadata index file */
    transportman::TransportMan m_MDIndexFileManager;
    /* transport manager for managing the metadata index file */
    transportman::TransportMan m_FileMetaMetadataManager;
    /* How many bytes of metadata index have we already read in? */
    size_t m_MDIndexFileAlreadyReadSize = 0;

    /* How many bytes of meta-metadata have we already read in? */
    size_t m_MetaMetaDataFileAlreadyReadSize = 0;
    /* How many bytes of meta-metadata have we already processed? */
    size_t m_MetaMetaDataFileAlreadyProcessedSize = 0;

    /* transport manager for managing the active flag file */
    transportman::TransportMan m_ActiveFlagFileManager;
    bool m_WriterIsActive = true;

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

    /** used for per-step reads, TODO: to be moved to BP5Deserializer */
    size_t m_CurrentStep = 0;
    size_t m_StepsCount = 0;
    size_t m_AbsStepsInFile = 0;    // all steps parsed including unselected
    uint64_t m_LastMapStep = 0;     // remember last step that had writer map
    uint64_t m_LastWriterCount = 0; // remember writer count in that step
    bool m_FirstStep = true;

    /** used to filter steps */
    helper::RangeFilter m_SelectedSteps;

    // offset/size pairs to read sections of metadata from file in InitBuffer
    std::vector<std::pair<uint64_t, uint64_t>> m_FilteredMetadataInfo;

    Minifooter m_Minifooter;

    void Init();
    void InitParameters();
    void InitTransports();

    /** DAOS pool connection and container opening */
    void InitDAOS();

    /* Sleep up to pollSeconds time if we have not reached timeoutInstant.
     * Return true if slept
     * return false if sleep was not needed because it was overtime
     */
    bool SleepOrQuit(const TimePoint &timeoutInstant, const Seconds &pollSeconds);
    /** Open one category of files within timeout.
     * @return: 0 = OK, 1 = timeout, 2 = error
     * lasterrmsg contains the error message in case of error
     */
    size_t OpenWithTimeout(transportman::TransportMan &tm,
                           const std::vector<std::string> &fileNames,
                           const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                           std::string &lasterrmsg /*INOUT*/);

    /** Open files within timeout.
     * @return True if files are opened, False in case of timeout
     */
    void OpenFiles(TimePoint &timeoutInstant, const Seconds &pollSeconds,
                   const Seconds &timeoutSeconds);

    /** Read in metadata if exist (throwing away old).
     *  It reads and parses metadata-index, and reads metadata into memory.
     *  In streaming mode, only a limited size of metadata is read in.
     *  Changes in m_StepsCount before and after calling can be used to
     *  track if new steps (after filtering with SelectSteps) are read in
     *  and are ready to be processed.
     */
    void UpdateBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                      const Seconds &timeoutSeconds);

    bool ReadActiveFlag(std::vector<char> &buffer);

    /* Parse metadata.
     *
     * Return the size of metadataindex where parsing stopped. In streaming mode
     * parsing is limited to read only a certain size of metadata at once.
     *
     * As a side effect, the following variables are filled out:
     *   m_MetadataIndexTable
     *   m_WriterMapIndex
     *   m_FilteredMetadataInfo
     */
    size_t ParseMetadataIndex(format::BufferSTL &bufferSTL, const size_t absoluteStartPos,
                              const bool hasHeader);

    /** Process the new metadata coming in (in UpdateBuffer)
     *  @param newIdxSize: the size of the new content from Index Table
     */
    void ProcessMetadataForNewSteps(const size_t newIdxSize);

    /** Check the active status of the writer.
     *  @return true if writer is still active.
     *  It sets m_WriterIsActive.
     */
    bool CheckWriterActive();

    /** Check for a step that is already in memory but haven't
     * been processed yet.
     *  @return true: if new step has been found and processed, false otherwise
     *  Used by CheckForNewSteps() to get the next step from memory if there is
     * one.
     */
    bool ProcessNextStepInMemory();

    /** Check for new steps withing timeout and only if writer is active.
     *  @return the status flag
     *  Used by BeginStep() to get new steps from file when it reaches the
     *  end of steps in memory.
     */
    StepStatus CheckForNewSteps(Seconds timeoutSeconds);

    /** Notify the engine when InquireVariable is called when the IO is empty.
     * Called from IO.tcc
     */
    void NotifyEngineNoVarsQuery();

#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex = -1) final;

    void GetSyncCommon(VariableBase &variable, void *data);

    void GetDeferredCommon(VariableBase &variable, void *data);

    void DoGetStructSync(VariableStruct &, void *);
    void DoGetStructDeferred(VariableStruct &, void *);

    template <class T>
    void ReadVariableBlocks(Variable<T> &variable);

    size_t DoSteps() const final;

    void DoGetAbsoluteSteps(const VariableBase &variable, std::vector<size_t> &keys) const final;

    uint32_t m_WriterColumnMajor = 0;
    bool m_ReaderIsRowMajor = true;
    bool m_WriterIsRowMajor = true;

    format::BufferSTL m_MetadataIndex;
    format::BufferSTL m_MetaMetadata;
    format::BufferSTL m_Metadata;

    void InstallMetaMetaData(format::BufferSTL MetaMetadata);
    void InstallMetadataForTimestep(size_t Step);
    std::pair<double, double> ReadData(adios2::transportman::TransportMan &FileManager,
                                       const size_t maxOpenFiles, const size_t WriterRank,
                                       const size_t Timestep, const size_t StartOffset,
                                       const size_t Length, char *Destination);

    struct WriterMapStruct
    {
        uint32_t WriterCount = 0;
        uint32_t AggregatorCount = 0;
        uint32_t SubfileCount = 0;
        std::vector<uint64_t> RankToSubfile; // size WriterCount
    };

    // step -> writermap but not for all steps
    std::map<uint64_t, WriterMapStruct> m_WriterMap;
    // step -> writermap index (for all steps)
    std::vector<uint64_t> m_WriterMapIndex;

    void DestructorClose(bool Verbose) noexcept;

    /* Communicator connecting ranks on each Compute Node.
       Only used to calculate the number of threads available for reading */
    helper::Comm m_NodeComm;
    helper::Comm singleComm;
    unsigned int m_Threads;
    std::vector<transportman::TransportMan> fileManagers; // manager per thread

    void daos_handle_share(daos_handle_t *, int);
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_DAOS_DAOSREADER_H_ */
