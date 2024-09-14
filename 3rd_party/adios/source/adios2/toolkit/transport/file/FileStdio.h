/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * FileStdio.h wrapper of C/C++ stdio.h for file I/O
 *
 *  Created on: Jan 6, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_FILE_FILEPOINTER_H_
#define ADIOS2_TOOLKIT_TRANSPORT_FILE_FILEPOINTER_H_

#include <cstdio> // FILE*
#include <future> //std::async, std::future

#include "adios2/toolkit/transport/Transport.h"

namespace adios2
{
namespace helper
{
class Comm;
}
namespace transport
{

/** File transport using C stdio FILE* */
class FileStdio : public Transport
{

public:
    FileStdio(helper::Comm const &comm);

    ~FileStdio();

    /** directio option is ignored in this transport */
    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) final;

    void OpenChain(const std::string &name, Mode openMode, const helper::Comm &chainComm,
                   const bool async = false, const bool directio = false) final;

    void SetBuffer(char *buffer, size_t size) final;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) final;

    void Read(char *buffer, size_t size, size_t start = MaxSizeT) final;

    size_t GetSize() final;

    void Flush() final;

    void Close() final;

    void Delete() final;

    void SeekToEnd() final;

    void SeekToBegin() final;

    void Seek(const size_t start) final;

    void Truncate(const size_t length) final;

    void MkDir(const std::string &fileName) final;

private:
    /** C File pointer */
    FILE *m_File = nullptr;
    bool m_IsOpening = false;
    std::future<FILE *> m_OpenFuture;

    /** Buffer settings need to be delayed until after opening */
    bool m_DelayedBufferSet = false;
    char *m_DelayedBuffer = nullptr;
    size_t m_DelayedBufferSize = 0;

    /**
     * Check for std::ferror and throw an exception if true
     * @param hint exception message
     */
    void CheckFile(const std::string hint) const;
    void WaitForOpen();
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_TRANSPORT_FILE_FILEPOINTER_H_ */
