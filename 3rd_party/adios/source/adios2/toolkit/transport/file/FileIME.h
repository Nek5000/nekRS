/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileIME.h file I/O using IME Native library
 *
 *  Created on: Dec 1, 2019
 *      Author: Keichi Takahashi keichi@is.naist.jp
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_FILE_IME_H_
#define ADIOS2_TOOLKIT_TRANSPORT_FILE_IME_H_

#include <atomic>

#include "adios2/common/ADIOSConfig.h"
#include "adios2/toolkit/transport/Transport.h"

namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File descriptor transport using the IME Native library */
class FileIME : public Transport
{

public:
    FileIME(helper::Comm const &comm);

    ~FileIME();

    /** Async option is ignored in FileIME transport */
    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) final;

    void SetParameters(const Params &parameters) final;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) final;

    void Read(char *buffer, size_t size, size_t start = MaxSizeT) final;

    size_t GetSize() final;

    /** Triggers a flush from IME to PFS if SyncToPFS=ON */
    void Flush() final;

    void Close() final;

    void Delete() final;

    void SeekToEnd() final;

    void SeekToBegin() final;

    void Seek(const size_t start = MaxSizeT) final;

    void Truncate(const size_t length) final;

    void MkDir(const std::string &fileName) final;

private:
    /** POSIX file handle returned by Open */
    int m_FileDescriptor = -1;

    /**
     * Check if m_FileDescriptor is -1 after an operation
     * @param hint exception message
     */
    void CheckFile(const std::string hint) const;

    /** Number of existing FileIME instances. We use atomic type for thread
     * safety but it is unclear if the IME client itself is thread safe. */
    static std::atomic_uint client_refcount;

    bool m_SyncToPFS = true;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TRANSPORT_FILE_FILEDESCRIPTOR_H_ */
