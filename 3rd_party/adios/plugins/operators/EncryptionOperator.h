/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * EncryptionOperator.h
 *
 *  Created on: Dec 7, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#ifndef ENCRYPTIONOPERATOR_H_
#define ENCRYPTIONOPERATOR_H_

#include <memory>

#include "adios2/common/ADIOSTypes.h"
#include "adios2/operator/plugin/PluginOperatorInterface.h"

namespace adios2
{
namespace plugin
{

/** EncryptionOperator that uses libsodium. Secret-key encryption is used.
 * The user must provide the 'SecretKeyFile' param with a path to the secret
 * key. If this file exists, it will use the key in the file (which should
 * have been generated using libsodium's crypto_secretbox_keygen().)
 * If the file doesn't exist, the operator will generate it using this call
 * and write it out to the specified file.
 */
class EncryptionOperator : public PluginOperatorInterface
{
public:
    EncryptionOperator(const Params &parameters);
    virtual ~EncryptionOperator();

    size_t Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                   const DataType type, char *bufferOut) override;

    size_t InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut) override;

    bool IsDataTypeValid(const DataType type) const override;

    size_t GetEstimatedSize(const size_t ElemCount, const size_t ElemSize, const size_t ndims,
                            const size_t *dims) const override;

private:
    struct EncryptImpl;
    std::unique_ptr<EncryptImpl> Impl;
};

} // end namespace plugin
} // end namespace adios2

extern "C" {

adios2::plugin::EncryptionOperator *OperatorCreate(const adios2::Params &parameters);
void OperatorDestroy(adios2::plugin::EncryptionOperator *obj);
}

#endif /* ENCYRPTIONOPERATOR_H_ */
