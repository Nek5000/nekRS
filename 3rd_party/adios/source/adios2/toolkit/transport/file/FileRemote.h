/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_FILE_Remote_H_
#define ADIOS2_TOOLKIT_TRANSPORT_FILE_Remote_H_

#include <future> //std::async, std::future

#include "adios2/common/ADIOSConfig.h"
#include "adios2/toolkit/remote/Remote.h"
#include "adios2/toolkit/transport/Transport.h"

#include <evpath.h>
namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File descriptor transport using the AWSSDK IO library */
class FileRemote : public Transport
{

public:
    FileRemote(helper::Comm const &comm);

    ~FileRemote();

    void SetParameters(const Params &parameters);

    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) final;

    void OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                   const bool async = false, const bool directio = false) final;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) final;

#ifdef REALLY_WANT_WRITEV
    /* Actual writev() function, inactive for now */
    void WriteV(const core::iovec *iov, const int iovcnt, size_t start = MaxSizeT) final;
#endif

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
    // class Impl;
    // static class Impl m_ImplSingleton;
    // Impl *m_Impl;
    // std::unique_ptr<Impl> m_Impl;
    Remote m_Remote;
    int m_Errno = 0;
    bool m_IsOpening = false;
    std::future<int> m_OpenFuture;
    size_t m_SeekPos = 0;
    size_t m_Size = 0;

    void SetUpCache();
    std::string m_FileName;
    std::string m_CachePath;        // local cache directory
    bool m_CachingThisFile = false; // save content to local cache
    bool m_IsCached = false;        // true if file is already in cache
    void CheckFile(const std::string hint) const;
    void WaitForOpen();
    std::string SysErrMsg() const;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TRANSPORT_FILE_Remote_H_ */
