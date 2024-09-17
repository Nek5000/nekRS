/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * OperatorFactory.cpp :
 *
 *  Created on: Sep 29, 2021
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#include "OperatorFactory.h"
#include "adios2/helper/adiosFunctions.h"
#include "adios2/operator/compress/CompressNull.h"
#include "adios2/operator/plugin/PluginOperator.h"
#include <numeric>

#ifdef ADIOS2_HAVE_BLOSC2
#include "adios2/operator/compress/CompressBlosc.h"
#endif

#ifdef ADIOS2_HAVE_BZIP2
#include "adios2/operator/compress/CompressBZIP2.h"
#endif

#ifdef ADIOS2_HAVE_LIBPRESSIO
#include "adios2/operator/compress/CompressLibPressio.h"
#endif

#ifdef ADIOS2_HAVE_MGARD
#include "adios2/operator/compress/CompressMGARD.h"
#include "adios2/operator/compress/CompressMGARDPlus.h"
#endif

#ifdef ADIOS2_HAVE_MGARD_MDR
#include "adios2/operator/refactor/RefactorMDR.h"
#endif

#ifdef ADIOS2_HAVE_PNG
#include "adios2/operator/compress/CompressPNG.h"
#endif

#ifdef ADIOS2_HAVE_MHS
#include "adios2/operator/compress/CompressSirius.h"
#endif

#ifdef ADIOS2_HAVE_SZ
#include "adios2/operator/compress/CompressSZ.h"
#endif

#ifdef ADIOS2_HAVE_ZFP
#include "adios2/operator/compress/CompressZFP.h"
#endif

namespace adios2
{
namespace core
{

std::string OperatorTypeToString(const Operator::OperatorType type)
{
    switch (type)
    {
    case Operator::COMPRESS_BLOSC:
        return "blosc";
    case Operator::COMPRESS_BZIP2:
        return "bzip2";
    case Operator::COMPRESS_LIBPRESSIO:
        return "libpressio";
    case Operator::COMPRESS_MGARD:
        return "mgard";
    case Operator::COMPRESS_MGARDPLUS:
        return "mgardplus";
    case Operator::COMPRESS_PNG:
        return "png";
    case Operator::COMPRESS_SIRIUS:
        return "sirius";
    case Operator::COMPRESS_SZ:
        return "sz";
    case Operator::COMPRESS_ZFP:
        return "zfp";
    case Operator::REFACTOR_MDR:
        return "mdr";
    case Operator::PLUGIN_INTERFACE:
        return "plugin";
    default:
        return "null";
    }
}

std::shared_ptr<Operator> MakeOperator(const std::string &type, const Params &parameters)
{
    std::shared_ptr<Operator> ret = nullptr;

    const std::string typeLowerCase = helper::LowerCase(type);

    if (typeLowerCase == "blosc")
    {
#ifdef ADIOS2_HAVE_BLOSC2
        ret = std::make_shared<compress::CompressBlosc>(parameters);
#endif
    }
    else if (typeLowerCase == "bzip2")
    {
#ifdef ADIOS2_HAVE_BZIP2
        ret = std::make_shared<compress::CompressBZIP2>(parameters);
#endif
    }
    else if (typeLowerCase == "libpressio")
    {
#ifdef ADIOS2_HAVE_LIBPRESSIO
        ret = std::make_shared<compress::CompressLibPressio>(parameters);
#endif
    }
    else if (typeLowerCase == "mgard")
    {
#ifdef ADIOS2_HAVE_MGARD
        ret = std::make_shared<compress::CompressMGARD>(parameters);
#endif
    }
    else if (typeLowerCase == "mgardplus")
    {
#ifdef ADIOS2_HAVE_MGARD
        ret = std::make_shared<compress::CompressMGARDPlus>(parameters);
#endif
    }
    else if (typeLowerCase == "png")
    {
#ifdef ADIOS2_HAVE_PNG
        ret = std::make_shared<compress::CompressPNG>(parameters);
#endif
    }
    else if (typeLowerCase == "sirius")
    {
#ifdef ADIOS2_HAVE_MHS
        ret = std::make_shared<compress::CompressSirius>(parameters);
#endif
    }
    else if (typeLowerCase == "sz")
    {
#ifdef ADIOS2_HAVE_SZ
        ret = std::make_shared<compress::CompressSZ>(parameters);
#endif
    }
    else if (typeLowerCase == "zfp")
    {
#ifdef ADIOS2_HAVE_ZFP
        ret = std::make_shared<compress::CompressZFP>(parameters);
#endif
    }
    else if (typeLowerCase == "mdr")
    {
#ifdef ADIOS2_HAVE_MGARD_MDR
        ret = std::make_shared<refactor::RefactorMDR>(parameters);
#endif
    }
    else if (typeLowerCase == "plugin")
    {
        ret = std::make_shared<plugin::PluginOperator>(parameters);
    }
    else if (typeLowerCase == "null")
    {
        ret = std::make_shared<compress::CompressNull>(parameters);
    }
    else
    {
        helper::Throw<std::invalid_argument>("Operator", "OperatorFactory", "MakeOperator",
                                             "ADIOS2 does not support " + typeLowerCase +
                                                 " operation");
    }

    if (ret == nullptr)
    {
        helper::Throw<std::invalid_argument>("Operator", "OperatorFactory", "MakeOperator",
                                             "ADIOS2 didn't compile with " + typeLowerCase +
                                                 " library, operator not added");
    }

    return ret;
}

size_t Decompress(const char *bufferIn, const size_t sizeIn, char *dataOut, MemorySpace memSpace,
                  std::shared_ptr<Operator> op)
{
    Operator::OperatorType compressorType;
    std::memcpy(&compressorType, bufferIn, 1);
    if (op == nullptr || op->m_TypeEnum != compressorType)
    {
        op = MakeOperator(OperatorTypeToString(compressorType), {});
    }
    size_t sizeOut = op->InverseOperate(bufferIn, sizeIn, dataOut);
    if (sizeOut == 0) // the inverse operator was not applied
    {
        size_t headerSize = op->GetHeaderSize();
        sizeOut = sizeIn - headerSize;
        helper::CopyContiguousMemory(bufferIn + headerSize, sizeOut, dataOut,
                                     /*endianReverse*/ false, memSpace);
    }
    return sizeOut;
}

} // end namespace core
} // end namespace adios2
