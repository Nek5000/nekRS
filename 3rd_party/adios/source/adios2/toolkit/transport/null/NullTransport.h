/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 *  NullTransport.h
 *
 *  Created on: Apr 17 2019
 *      Author: Chuck Atkins <chuck.atkins@kitware.com>
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_NULL_NULLTRANSPORT_H_
#define ADIOS2_TOOLKIT_TRANSPORT_NULL_NULLTRANSPORT_H_

#include <memory>

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

/** File stream transport using C++ fstream */
class NullTransport : public Transport
{

public:
    NullTransport(helper::Comm const &comm);

    virtual ~NullTransport();

    void Open(const std::string &name, const Mode openMode, const bool async = false,
              const bool directio = false) override;

    void SetBuffer(char *buffer, size_t size) override;

    void Write(const char *buffer, size_t size, size_t start = MaxSizeT) override;

    void Read(char *buffer, size_t size, size_t start = MaxSizeT) override;

    size_t GetSize() override;

    void Flush() override;

    void Close() override;

    void Delete() final;

    void SeekToEnd() override;

    void SeekToBegin() override;

    void Seek(const size_t start = MaxSizeT) override;

    void Truncate(const size_t length) override;

protected:
    struct NullTransportImpl;
    std::unique_ptr<NullTransportImpl> Impl;

    void MkDir(const std::string &fileName) override;

    void CheckName() const override;
};

} // end namespace transport
} // end namespace adios2

#endif /* ADIOS2_TOOLKIT_TRANSPORT_NULL_NULLTRANSPORT_H_ */
