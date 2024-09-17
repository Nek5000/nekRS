/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DynamicBinder.h
 *
 * Created on: Jul 21, 2017
 *     Author: Chuck Atkins
 */

#ifndef ADIOS2_HELPER_DYNAMICBINDER_H_
#define ADIOS2_HELPER_DYNAMICBINDER_H_

#include <memory>
#include <string>
#include <type_traits>

namespace adios2
{
namespace helper
{

class DynamicBinder
{
public:
    using VoidSymbolPointer = std::add_pointer<void()>::type;

public:
    DynamicBinder(std::string libName);
    DynamicBinder(std::string libName, std::string libPath);
    ~DynamicBinder();

    VoidSymbolPointer GetSymbol(std::string symbolName);

private:
    struct Impl;
    std::unique_ptr<Impl> m_Impl;
};

} // end namespace helper
} // end namespace adios2

#endif // ADIOS2_HELPER_DYNAMICBINDER_H_
