/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Types.h : provides type utilities for ADIOS2 C++11 bindings
 *
 *  Created on: Feb 11, 2019
 *      Author: Kai Germaschewski <kai.germaschewski@unh.edu>
 *              William F Godoy <godoywf@ornl.gov>
 */

#ifndef ADIOS2_BINDINGS_CXX11_TYPES_H_
#define ADIOS2_BINDINGS_CXX11_TYPES_H_

#include "Variable.h"

#include <vector>

namespace adios2
{

template <class T>
std::string GetType() noexcept;

} // end namespace adios2
#endif /* ADIOS2_BINDINGS_CXX11_TYPES_H_ */
