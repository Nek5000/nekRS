/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * DynamicBinder.cpp
 *
 * Created on: Jul 21, 2017
 *     Author: Chuck Atkins
 */

#include "adiosDynamicBinder.h"
#include "adiosLog.h"

#include <algorithm> // for copy
#include <iostream>  // for operator<<, stringstream, bas...
#include <iterator>  // for ostream_iterator
#include <sstream>   // for stringstream
#include <stdexcept> // for runtime_error
#include <vector>    // for vector

#include <adios2sys/DynamicLoader.hxx>
#include <adios2sys/SystemTools.hxx>

namespace adios2
{
namespace helper
{

struct DynamicBinder::Impl
{
    adios2sys::DynamicLoader::LibraryHandle m_LibraryHandle;
};

DynamicBinder::DynamicBinder(std::string libName) : DynamicBinder(libName, "") {}

DynamicBinder::DynamicBinder(std::string libName, std::string libPath) : m_Impl(new Impl)
{
    std::vector<std::string> libPrefixes;
    libPrefixes.emplace_back("");
    libPrefixes.emplace_back("lib");
#ifdef __CYGWIN__
    libPrefixes.emplace_back("cyg");
#endif

    std::vector<std::string> libSuffixes;
    libSuffixes.emplace_back("");
#ifdef __APPLE__
    libSuffixes.emplace_back(".dylib");
    libSuffixes.emplace_back(".so");
#endif
#ifdef __hpux
    libSuffixes.emplace_back(".sl");
#endif
#ifdef __unix__
    libSuffixes.emplace_back(".so");
#endif
#ifdef _WIN32
    libSuffixes.emplace_back(".dll");
#endif

    std::vector<std::string> searchedLibs;
    std::string fileName;

    // Test the various combinations of library names
    for (const std::string &prefix : libPrefixes)
    {
        for (const std::string &suffix : libSuffixes)
        {
            if (!libPath.empty())
            {
                fileName = libPath + "/" + prefix + libName + suffix;
                // Slashes in fileName is correct for unix-like systems
                // ConvertToOutputPath() will change slashes if we're running on
                // a Windows system
                fileName = adios2sys::SystemTools::ConvertToOutputPath(fileName);
            }
            else
            {
                fileName = prefix + libName + suffix;
            }
            m_Impl->m_LibraryHandle = adios2sys::DynamicLoader::OpenLibrary(fileName);
            searchedLibs.push_back(fileName);
            if (m_Impl->m_LibraryHandle)
            {
                break;
            }
        }
        if (m_Impl->m_LibraryHandle)
        {
            break;
        }
    }
    if (!m_Impl->m_LibraryHandle)
    {
        std::stringstream errString;
        errString << "Unable to locate the " << libName << " library; searched for ";
        std::copy(searchedLibs.begin(), searchedLibs.end(),
                  std::ostream_iterator<std::string>(errString, " "));

        helper::Throw<std::runtime_error>("Helper", "adiosDynamicBinder", "DynamicBinder",
                                          errString.str());
    }
}

DynamicBinder::~DynamicBinder() { adios2sys::DynamicLoader::CloseLibrary(m_Impl->m_LibraryHandle); }

DynamicBinder::VoidSymbolPointer DynamicBinder::GetSymbol(std::string symbolName)
{
    return adios2sys::DynamicLoader::GetSymbolAddress(m_Impl->m_LibraryHandle, symbolName);
}

} // end namespace helper
} // end namespace adios2
