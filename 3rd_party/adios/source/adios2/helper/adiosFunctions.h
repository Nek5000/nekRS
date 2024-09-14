/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adiosFunctions.h helper functions used by ADIOS class. This is the only
 * required "public" header.
 *
 *  Created on: Oct 10, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_HELPER_ADIOSFUNCTIONS_H_
#define ADIOS2_HELPER_ADIOSFUNCTIONS_H_

#include "adios2/helper/adiosGPUFunctions.h" //GPU functions
#include "adios2/helper/adiosLog.h"          //logging functions
#include "adios2/helper/adiosMath.h"         //math functions (cmath, algorithm)
#include "adios2/helper/adiosMemory.h"       //memcpy, std::copy, insert, resize
#include "adios2/helper/adiosNetwork.h"      //network and staging functions
#include "adios2/helper/adiosString.h"       //std::string manipulation
#include "adios2/helper/adiosSystem.h"       //OS functionality, POSIX, filesystem
#include "adios2/helper/adiosType.h"         //Type casting, conversion, checks, etc.
#include "adios2/helper/adiosXML.h"          //XML parsing
#include "adios2/helper/adiosYAML.h"         //YAML parsing

#endif /* ADIOS2_HELPER_ADIOSFUNCTIONS_H_ */
