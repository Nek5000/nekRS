/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Reader.h
 *
 *  Created on: Aug 1, 2018
 *      Author: Lipeng Wan wanl@ornl.gov
 */

#ifndef ADIOS2_ENGINE_BP4_BP4READER_H_
#define ADIOS2_ENGINE_BP4_BP4READER_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/CoreTypes.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/format/bp/bp4/BP4Deserializer.h"
#include "adios2/toolkit/transportman/TransportMan.h"

namespace adios2
{
namespace core
{
namespace engine
{

class BP4Reader : public Engine
{

public:
    /**
     * Unique constructor
     * @param io
     * @param name
     * @param openMode only read
     * @param comm
     */
    BP4Reader(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    virtual ~BP4Reader();

    StepStatus BeginStep(StepMode mode = StepMode::Read, const float timeoutSeconds = -1.0) final;

    size_t CurrentStep() const final;

    void EndStep() final;

    void PerformGets() final;

protected:
    void DestructorClose(bool Verbose) noexcept {};

private:
    format::BP4Deserializer m_BP4Deserializer;
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
    /* How many bytes of metadata index have we already read in? */
    size_t m_MDIndexFileAlreadyReadSize = 0;

    /* transport manager for managing the active flag file */
    transportman::TransportMan m_ActiveFlagFileManager;
    bool m_WriterIsActive = true;

    /** used for per-step reads, TODO: to be moved to BP4Deserializer */
    size_t m_CurrentStep = 0;
    bool m_FirstStep = true;
    bool m_IdxHeaderParsed = false; // true after first index parsing

    int m_Verbosity = 0;

    void Init();
    void InitTransports();

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
    void InitBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds,
                    const Seconds &timeoutSeconds);

    /** Read in more metadata if exist (throwing away old).
     *  For streaming only.
     *  @return size of new content from Index Table
     */
    size_t UpdateBuffer(const TimePoint &timeoutInstant, const Seconds &pollSeconds);

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

#define declare_type(T)                                                                            \
    void DoGetSync(Variable<T> &, T *) final;                                                      \
    void DoGetDeferred(Variable<T> &, T *) final;
    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    void DoClose(const int transportIndex = -1) final;

    template <class T>
    void GetSyncCommon(Variable<T> &variable, T *data);

    template <class T>
    void GetDeferredCommon(Variable<T> &variable, T *data);

    template <class T>
    void ReadVariableBlocks(Variable<T> &variable);

#define declare_type(T)                                                                            \
    std::map<size_t, std::vector<typename Variable<T>::BPInfo>> DoAllStepsBlocksInfo(              \
        const Variable<T> &variable) const final;                                                  \
                                                                                                   \
    std::vector<std::vector<typename Variable<T>::BPInfo>> DoAllRelativeStepsBlocksInfo(           \
        const Variable<T> &) const final;                                                          \
                                                                                                   \
    std::vector<typename Variable<T>::BPInfo> DoBlocksInfo(const Variable<T> &variable,            \
                                                           const size_t step) const final;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    size_t DoSteps() const final;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_BP4_BP4READER_H_ */
