/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * IO.tcc :
 *
 *  Created on: Jun 4, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */
#ifndef ADIOS2_BINDINGS_CXX11_CXX11_IO_TCC_
#define ADIOS2_BINDINGS_CXX11_CXX11_IO_TCC_

#include "IO.h"

#include "adios2/core/IO.h"

namespace adios2
{

template <class T>
Variable<T> IO::DefineVariable(const std::string &name, const Dims &shape, const Dims &start,
                               const Dims &count, const bool constantDims)
{
    helper::CheckForNullptr(m_IO, "for variable name " + name + ", in call to IO::DefineVariable");
    return Variable<T>(&m_IO->DefineVariable<typename TypeInfo<T>::IOType>(name, shape, start,
                                                                           count, constantDims));
}

template <class T>
Variable<T> IO::InquireVariable(const std::string &name)
{
    helper::CheckForNullptr(m_IO, "for variable name " + name + ", in call to IO::InquireVariable");
    return Variable<T>(m_IO->InquireVariable<typename TypeInfo<T>::IOType>(name));
}

template <class T>
Attribute<T> IO::DefineAttribute(const std::string &name, const T *data, const size_t size,
                                 const std::string &variableName, const std::string separator,
                                 const bool allowModification)
{
    using IOType = typename TypeInfo<T>::IOType;
    helper::CheckForNullptr(m_IO, "for attribute name " + name + " and variable name " +
                                      variableName + ", in call to IO::DefineAttribute");
    return Attribute<T>(&m_IO->DefineAttribute(name, reinterpret_cast<const IOType *>(data), size,
                                               variableName, separator, allowModification));
}

template <class T>
Attribute<T> IO::DefineAttribute(const std::string &name, const T &value,
                                 const std::string &variableName, const std::string separator,
                                 const bool allowModification)
{
    using IOType = typename TypeInfo<T>::IOType;
    helper::CheckForNullptr(m_IO,
                            "for attribute name " + name + ", in call to IO::DefineAttribute");
    return Attribute<T>(&m_IO->DefineAttribute(name, reinterpret_cast<const IOType &>(value),
                                               variableName, separator, allowModification));
}

template <class T>
Attribute<T> IO::InquireAttribute(const std::string &name, const std::string &variableName,
                                  const std::string separator)
{
    using IOType = typename TypeInfo<T>::IOType;
    helper::CheckForNullptr(m_IO,
                            "for attribute name " + name + ", in call to IO::InquireAttribute");
    return Attribute<T>(m_IO->InquireAttribute<IOType>(name, variableName, separator));
}

} // end namespace adios2

#endif // ADIOS2_BINDINGS_CXX11_CXX11_IO_TCC_
