/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP4Writer.h
 *
 *  Created on: Aug 1, 2018
 *      Author: Lipeng Wan wanl@ornl.gov
 */

#ifndef ADIOS2_ENGINE_BP4_BP4WRITER_H_
#define ADIOS2_ENGINE_BP4_BP4WRITER_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/core/Engine.h"
#include "adios2/helper/adiosComm.h"
#include "adios2/toolkit/burstbuffer/FileDrainerSingleThread.h"
#include "adios2/toolkit/format/bp/bp4/BP4Serializer.h"
#include "adios2/toolkit/transportman/TransportMan.h"

namespace adios2
{
namespace core
{
namespace engine
{

class BP4Writer : public core::Engine
{

public:
    /**
     * Constructor for file Writer in BP4 format
     * @param name unique name given to the engine
     * @param openMode w (supported), r, a from OpenMode in ADIOSTypes.h
     * @param comm multi-process communicator
     */
    BP4Writer(IO &io, const std::string &name, const Mode mode, helper::Comm comm);

    ~BP4Writer();

    StepStatus BeginStep(StepMode mode, const float timeoutSeconds = -1.0) final;
    size_t CurrentStep() const final;
    void PerformPuts() final;
    void EndStep() final;
    void Flush(const int transportIndex = -1) final;

    size_t DebugGetDataBufferSize() const final;

protected:
    void DestructorClose(bool Verbose) noexcept;

private:
    /** Single object controlling BP buffering */
    format::BP4Serializer m_BP4Serializer;

    /** Manage BP data files Transports from IO AddTransport */
    transportman::TransportMan m_FileDataManager;

    /** Manages the optional collective metadata files */
    transportman::TransportMan m_FileMetadataManager;

    /* transport manager for managing the metadata index file */
    transportman::TransportMan m_FileMetadataIndexManager;

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
    std::vector<std::string> m_MetadataIndexFileNames;
    std::vector<std::string> m_DrainMetadataIndexFileNames;
    std::vector<std::string> m_ActiveFlagFileNames;

    int m_Verbosity = 0;
    // true if BeginStep was ever called
    bool m_DidBeginStep = false;

    void Init() final;

    /** Parses parameters from IO SetParameters */
    void InitParameters() final;
    /** Parses transports and parameters from IO AddTransport */
    void InitTransports() final;
    /** Allocates memory and starts a PG group */
    void InitBPBuffer();

#define declare_type(T)                                                                            \
    void DoPut(Variable<T> &variable, typename Variable<T>::Span &span, const bool initialize,     \
               const T &value) final;

    ADIOS2_FOREACH_PRIMITIVE_STDTYPE_1ARG(declare_type)
#undef declare_type

#define declare_type(T)                                                                            \
    void DoPutSync(Variable<T> &, const T *) final;                                                \
    void DoPutDeferred(Variable<T> &, const T *) final;

    ADIOS2_FOREACH_STDTYPE_1ARG(declare_type)
#undef declare_type

    template <class T>
    void PutCommon(Variable<T> &variable, typename Variable<T>::Span &span, const size_t bufferID,
                   const T &value);

    template <class T>
    void PutSyncCommon(Variable<T> &variable, const typename Variable<T>::BPInfo &blockInfo,
                       const bool resize = true);

    template <class T>
    void PutDeferredCommon(Variable<T> &variable, const T *data);

    void DoFlush(const bool isFinal = false, const int transportIndex = -1);

    void DoClose(const int transportIndex = -1) final;

    /** Write a profiling.json file from m_BP1Writer and m_TransportsManager
     * profilers*/
    void WriteProfilingJSONFile();

    void PopulateMetadataIndexFileContent(format::BufferSTL &buffer, const uint64_t currentStep,
                                          const uint64_t mpirank, const uint64_t pgIndexStart,
                                          const uint64_t variablesIndexStart,
                                          const uint64_t attributesIndexStart,
                                          const uint64_t currentStepEndPos,
                                          const uint64_t currentTimeStamp);

    void UpdateActiveFlag(const bool active);

    void WriteCollectiveMetadataFile(const bool isFinal = false);

    /**
     * N-to-N data buffers writes, including metadata file
     * @param transportIndex
     */
    void WriteData(const bool isFinal, const int transportIndex = -1);

    /**
     * N-to-M (aggregation) data buffers writes, including metadata file
     * @param transportIndex
     */
    void AggregateWriteData(const bool isFinal, const int transportIndex = -1);

#define declare_type(T, L)                                                                         \
    T *DoBufferData_##L(const int bufferIdx, const size_t payloadPosition,                         \
                        const size_t bufferID = 0) noexcept final;

    ADIOS2_FOREACH_PRIMITVE_STDTYPE_2ARGS(declare_type)
#undef declare_type

    template <class T>
    T *BufferDataCommon(const int bufferIdx, const size_t payloadOffset,
                        const size_t bufferID) noexcept;

    template <class T>
    void PerformPutCommon(Variable<T> &variable);

    void NotifyEngineAttribute(std::string name, DataType type) noexcept;
    virtual void NotifyEngineAttribute(std::string name, AttributeBase *attr, void *Data) noexcept;
};

} // end namespace engine
} // end namespace core
} // end namespace adios2

#endif /* ADIOS2_ENGINE_BP4_BP4WRITER_H_ */
