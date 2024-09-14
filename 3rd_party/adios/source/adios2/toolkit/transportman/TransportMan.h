/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TransportMan.h : manages a vector of transports
 *
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_
#define ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_

#include <future> //std::async, std::future
#include <memory> //std::shared_ptr
#include <string>
#include <unordered_map>
#include <vector>

#include "adios2/core/CoreTypes.h"
#include "adios2/toolkit/transport/Transport.h"

namespace adios2
{
namespace helper
{
class Comm;
}
namespace core
{
class IO;
}
namespace transportman
{

class TransportMan
{

public:
    /**
     * Contains all transports
     * <pre>
     * key : unique id
     * value : object derived from Transport base class
     * </pre>
     */
    std::unordered_map<size_t, std::shared_ptr<Transport>> m_Transports;

    /**
     * Unique base constructor
     * @param comm
     */
    TransportMan(core::IO &IO, helper::Comm &comm);

    virtual ~TransportMan() = default;

    /**
     * Function that will be called from all ranks in communicator, only rank
     * zero creates directories
     * @param fileNames extract directory if needed to be created
     * @param nodeLocal true: all ranks create a directory
     */
    void MkDirsBarrier(const std::vector<std::string> &fileNames,
                       const std::vector<Params> &parametersVector, const bool nodeLocal);

    /**
     * OpenFiles passed from fileNames
     * @param fileNames
     * @param openMode
     * @param parametersVector from IO
     * @param profile
     */
    void OpenFiles(const std::vector<std::string> &fileNames, const Mode openMode,
                   const std::vector<Params> &parametersVector, const bool profile);

    /**
     * OpenFiles passed from fileNames, in a chain to avoid DOS attacking of the
     * file system
     * @param fileNames
     * @param openMode
     * @param parametersVector from IO
     * @param profile
     * @oaram chainComm
     */
    void OpenFiles(const std::vector<std::string> &fileNames, const Mode openMode,
                   const std::vector<Params> &parametersVector, const bool profile,
                   const helper::Comm &chainComm);

    /**
     * Used for sub-files defined by index
     * @param name
     * @param id
     * @param openMode
     * @param parameters
     * @param profile
     */
    void OpenFileID(const std::string &name, const size_t id, const Mode mode,
                    const Params &parameters, const bool profile);

    /**
     * Gets each transport base name from either baseName at Open or name
     * key in
     * parameters
     * Checks if transport name rules IO AddTransport have unique names for
     * every type (for now)
     * @param baseName from Open
     * @param parameters from IO TransportsParameters (from AddTransport
     * function)
     * @return transport base names
     */
    std::vector<std::string> GetFilesBaseNames(const std::string &baseName,
                                               const std::vector<Params> &parametersVector) const;

    /**
     * m_Type from m_Transports based on derived classes of Transport
     * @return m_Type for each transport in m_Transports (e.g.
     * {FileDescriptor,
     * FilePointer} )
     */
    std::vector<std::string> GetTransportsTypes() noexcept;

    /** Returns a vector of pointer references (not owning the memory) to
     * m_Transports.m_Profiler */
    std::vector<profiling::IOChrono *> GetTransportsProfilers() noexcept;

    /**
     * Write to file transports
     * @param transportIndex
     * @param buffer
     * @param size
     */
    void WriteFiles(const char *buffer, const size_t size, const int transportIndex = -1);

    /**
     * Write data to a specific location in files
     * @param transportIndex
     * @param buffer
     * @param size
     * @param start offset in file
     */
    void WriteFileAt(const char *buffer, const size_t size, const size_t start,
                     const int transportIndex = -1);

    /**
     * Write to file transports, writev version
     * @param transportIndex
     * @param iovec array pointer
     * @param iovcnt number of entries
     */
    void WriteFiles(const core::iovec *iov, const size_t iovcnt, const int transportIndex = -1);

    /**
     * Write data to a specific location in files, writev version
     * @param transportIndex
     * @param iovec array pointer
     * @param iovcnt number of entries
     * @param start offset in file
     */
    void WriteFileAt(const core::iovec *iov, const size_t iovcnt, const size_t start,
                     const int transportIndex = -1);

    size_t GetFileSize(const size_t transportIndex = 0) const;

    /**
     * Read contents from a single file and assign it to buffer
     * @param buffer
     * @param size
     * @param start
     * @param transportIndex
     */
    void ReadFile(char *buffer, const size_t size, const size_t start = 0,
                  const size_t transportIndex = 0);

    /**
     * Flush file or files depending on transport index. Throws an exception
     * if transport is not a file when transportIndex > -1.
     * @param transportIndex -1: all transports, otherwise index in m_Transports
     */
    void FlushFiles(const int transportIndex = -1);

    /**
     * Close file or files depending on transport index. Throws an exception
     * if transport is not a file when transportIndex > -1.
     * @param transportIndex -1: all transports, otherwise index in m_Transports
     */
    void CloseFiles(const int transportIndex = -1);

    /**
     * Delete file or files depending on transport index. The files
     * must be open for this function to have an effect.
     */
    void DeleteFiles(const int transportIndex = -1);

    /** Checks if all transports are closed */
    bool AllTransportsClosed() const noexcept;

    void SeekToFileEnd(const int transportIndex = -1);

    void SeekToFileBegin(const int transportIndex = -1);

    void SeekTo(const size_t start, const int transportIndex = -1);

    void Truncate(const size_t length, const int transportIndex = -1);

    /**
     * Check if a file exists.
     * @param name
     * @param parameters
     * @param profile
     */
    bool FileExists(const std::string &name, const Params &parameters, const bool profile);

    /**
     * Set Transport Paramers
     * @param params
     * @param transportIndex
     */
    void SetParameters(const Params &params, const int transportIndex = -1);

protected:
    core::IO &m_IO;
    helper::Comm const &m_Comm;

    std::shared_ptr<Transport> OpenFileTransport(const std::string &fileName, const Mode openMode,
                                                 const Params &parameters, const bool profile,
                                                 const bool useComm, const helper::Comm &chainComm);

    void
    CheckFile(std::unordered_map<size_t, std::shared_ptr<Transport>>::const_iterator itTransport,
              const std::string hint) const;

    void WaitForAsync() const;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_ */
