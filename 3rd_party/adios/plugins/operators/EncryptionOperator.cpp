/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * EncryptionOperator.cpp
 *
 *  Created on: Dec 7, 2021
 *      Author: Caitlin Ross <caitlin.ross@kitware.com>
 */

#include "EncryptionOperator.h"

#include <fstream>
#include <string>

#include <sodium.h>

namespace adios2
{
namespace plugin
{

struct EncryptionOperator::EncryptImpl
{
    std::string KeyFilename;
    unsigned char Key[crypto_secretbox_KEYBYTES] = {0};
    bool KeyValid = false;

    ~EncryptImpl()
    {
        // unlocks the memory location as well as zeroing out the data
        sodium_munlock(Key, crypto_secretbox_KEYBYTES);
    }

    void GenerateOrReadKey()
    {
        // if the key file already exists, we'll use that key,
        // otherwise we'll generate one and write it out
        std::fstream keyFile(KeyFilename.c_str());
        if (keyFile)
        {
            keyFile.read(reinterpret_cast<char *>(&Key), crypto_secretbox_KEYBYTES);
            keyFile.close();
        }
        else
        {
            keyFile.open(KeyFilename.c_str(), std::fstream::out);
            if (!keyFile)
            {
                throw std::runtime_error("couldn't open file to write key");
            }
            crypto_secretbox_keygen(Key);
            keyFile.write(reinterpret_cast<char *>(&Key), crypto_secretbox_KEYBYTES);
            keyFile.close();
        }

        // lock the key to avoid swapping to disk
        if (sodium_mlock(Key, crypto_secretbox_KEYBYTES) == -1)
        {
            throw std::runtime_error("Unable to lock memory location of secret key,"
                                     " due to system limit on amount of memory that can be locked "
                                     "by a process.");
        }
        KeyValid = true;
    }
};

EncryptionOperator::EncryptionOperator(const Params &parameters)
: PluginOperatorInterface(parameters), Impl(new EncryptImpl)
{
    if (sodium_init() < 0)
    {
        throw std::runtime_error("libsodium could not be initialized");
    }

    // in the case "secretkeyfile" is found, so we know the operator should
    // calling Operate(). If "secretkeyfile" is not found, then the operator
    // should be calling InverseOperate(), due to ADIOS calling InverseOperate()
    // not allowing Parameters to be passed.
    auto skFileIt = m_Parameters.find("secretkeyfile");
    if (skFileIt != m_Parameters.end())
    {
        Impl->KeyFilename = skFileIt->second;
        Impl->GenerateOrReadKey();
    }
}

EncryptionOperator::~EncryptionOperator() {}

#if defined(__clang__)
#if __has_feature(memory_sanitizer)
// Memory Sanitizer has an issue with some libsodium calls.
__attribute__((no_sanitize("memory")))
#endif
#endif
size_t
EncryptionOperator::Operate(const char *dataIn, const Dims &blockStart, const Dims &blockCount,
                            const DataType type, char *bufferOut)
{
    if (!Impl->KeyValid)
    {
        throw std::runtime_error("EncryptionOperator::Operate was called, but"
                                 " a valid secret key has not been generated. "
                                 "Did you add the SecretKeyFile"
                                 " param when setting up the operator?");
    }

    // offset for writing into bufferOut
    size_t offset = 0;

    // write any parameters we need to save for the InverseOperate() call
    // In this case, we just write out the size of the data
    size_t sizeIn = helper::GetTotalSize(blockCount, helper::GetDataTypeSize(type));
    PutParameter(bufferOut, offset, sizeIn);

    // create the nonce directly in the output buffer, since we'll need it for
    // decryption
    unsigned char *nonce = reinterpret_cast<unsigned char *>(bufferOut + offset);
    randombytes_buf(nonce, crypto_secretbox_NONCEBYTES);
    offset += crypto_secretbox_NONCEBYTES;

    // encrypt data directly into the output buffer
    size_t cipherTextSize = sizeIn + crypto_secretbox_MACBYTES;
    unsigned char *cipherText = reinterpret_cast<unsigned char *>(bufferOut + offset);
    crypto_secretbox_easy(cipherText, reinterpret_cast<const unsigned char *>(dataIn), sizeIn,
                          nonce, Impl->Key);
    offset += cipherTextSize;

    // need to return the size of data in the buffer
    return offset;
}

#if defined(__clang__)
#if __has_feature(memory_sanitizer)
// Memory Sanitizer has an issue with some libsodium calls.
__attribute__((no_sanitize("memory")))
#endif
#endif
size_t
EncryptionOperator::InverseOperate(const char *bufferIn, const size_t sizeIn, char *dataOut)
{
    size_t offset = 0;

    // need to grab any parameter(s) we saved in Operate()
    const size_t dataBytes = GetParameter<size_t>(bufferIn, offset);

    // grab the nonce ptr
    const unsigned char *nonce = reinterpret_cast<const unsigned char *>(bufferIn + offset);
    offset += crypto_secretbox_NONCEBYTES;

    // grab the cipher text ptr
    size_t cipherTextSize = dataBytes + crypto_secretbox_MACBYTES;
    const unsigned char *cipherText = reinterpret_cast<const unsigned char *>(bufferIn + offset);
    offset += cipherTextSize;

    // decrypt directly into dataOut buffer
    if (crypto_secretbox_open_easy(reinterpret_cast<unsigned char *>(dataOut), cipherText,
                                   cipherTextSize, nonce, Impl->Key) != 0)
    {
        throw std::runtime_error("message forged!");
    }

    // return the size of the data
    return dataBytes;
}

bool EncryptionOperator::IsDataTypeValid(const DataType type) const { return true; }

size_t EncryptionOperator::GetEstimatedSize(const size_t ElemCount, const size_t ElemSize,
                                            const size_t ndims, const size_t *dims) const
{
    size_t sizeIn = ElemCount * ElemSize;
    return (sizeof(size_t)                // Data size
            + crypto_secretbox_NONCEBYTES // Nonce
            + sizeIn                      // Data
            + crypto_secretbox_MACBYTES   // MAC
    );
}
} // end namespace plugin
} // end namespace adios2

extern "C" {

adios2::plugin::EncryptionOperator *OperatorCreate(const adios2::Params &parameters)
{
    return new adios2::plugin::EncryptionOperator(parameters);
}

void OperatorDestroy(adios2::plugin::EncryptionOperator *obj) { delete obj; }
}
