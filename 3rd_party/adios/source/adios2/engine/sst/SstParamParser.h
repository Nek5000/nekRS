#ifndef ADIOS2_ENGINE_SST_SSTPARAMPARSER_H_
#define ADIOS2_ENGINE_SST_SSTPARAMPARSER_H_

#include "adios2/common/ADIOSConfig.h"
#include "adios2/common/ADIOSTypes.h"
#include "adios2/core/IO.h"
#include <iostream>
#include <unistd.h>

#include "adios2/toolkit/sst/sst.h"

using namespace adios2::core;

class SstParamParser
{
public:
    void ParseParams(adios2::core::IO &io, struct _SstParams &Params,
                     const adios2::UserOptions &userOptions);
};

#endif
