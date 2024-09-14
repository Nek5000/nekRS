/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileDescriptor.h wrapper of DAOS library functions for file I/O
 *
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_FILE_DAOS_H_
#define ADIOS2_TOOLKIT_TRANSPORT_FILE_DAOS_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/toolkit/transport/Transport.h"

#include <future>
#include <memory>

namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File descriptor transport using the DAOS library */
class FileDaos : public Transport
{

public:
    FileDaos(helper::Comm const &comm);

    ~FileDaos();

    void SetParameters(const Params &parameters);

    /** directio option is ignored in this transport */
    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) final;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) final;

    void Read(char *buffer, size_t size, size_t start = MaxSizeT) final;

    size_t GetSize() final;

    /** Does nothing, each write is supposed to flush */
    void Flush() final;

    void Close() final;

    void Delete() final;

    void SeekToEnd() final;

    void SeekToBegin() final;

    void Seek(const size_t start = MaxSizeT) final;

    void Truncate(const size_t length) final;

    void MkDir(const std::string &fileName) final;

private:
    class Impl;
    std::unique_ptr<Impl> m_Impl;
    /** DAOS file handle reterned by Open */
    bool m_DAOSOpenSucceed = false;
    size_t m_GlobalOffset = 0;
    int m_Errno = 0;
    bool m_IsOpening = false;
    std::future<bool> m_OpenFuture;

    /**
     * Check if m_FileDescriptor is -1 after an operation
     * @param hint exception message
     */
    void CheckFile(const std::string hint) const;
    void WaitForOpen();
    std::string SysErrMsg() const;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TRANSPORT_FILE_DAOS_H_ */
