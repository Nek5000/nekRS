/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * adios2fstream.tcc
 *
 *  Created on: Mar 5, 2018
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_BINDINGS_CXX11_CXX11_FSTREAM_ADIOS2FSTREAM_TCC_
#define ADIOS2_BINDINGS_CXX11_CXX11_FSTREAM_ADIOS2FSTREAM_TCC_

#include "ADIOS2fstream.h"

#include "adios2/core/Stream.h"

namespace adios2
{

template <class T>
void fstream::write_attribute(const std::string &name, const T &value,
                              const std::string &variableName, const std::string separator,
                              const bool endStep)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->WriteAttribute(name, reinterpret_cast<const IOType &>(value), variableName, separator,
                             endStep);
}

template <class T>
void fstream::write_attribute(const std::string &name, const T *data, const size_t size,
                              const std::string &variableName, const std::string separator,
                              const bool endStep)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->WriteAttribute(name, reinterpret_cast<const IOType *>(data), size, variableName,
                             separator, endStep);
}

template <class T>
void fstream::write(const std::string &name, const T *data, const Dims &shape, const Dims &start,
                    const Dims &count, const bool endStep)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Write(name, reinterpret_cast<const IOType *>(data), shape, start, count, vParams(),
                    endStep);
}

template <class T>
void fstream::write(const std::string &name, const T *data, const adios2::Dims &shape,
                    const adios2::Dims &start, const adios2::Dims &count,
                    const adios2::vParams &operations, const bool endStep)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Write(name, reinterpret_cast<const IOType *>(data), shape, start, count, operations,
                    endStep);
}

template <class T>
void fstream::write(const std::string &name, const T &value, const bool isLocalValue,
                    const bool endStep)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Write(name, reinterpret_cast<const IOType &>(value), isLocalValue, endStep);
}

template <class T>
void fstream::read(const std::string &name, T *data, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, reinterpret_cast<IOType *>(data), blockID);
}

template <class T>
void fstream::read(const std::string &name, T *data, const size_t stepStart, const size_t stepCount,
                   const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, reinterpret_cast<IOType *>(data), Box<size_t>{stepStart, stepCount},
                   blockID);
}

template <class T>
void fstream::read(const std::string &name, T &value, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, &reinterpret_cast<IOType &>(value), blockID);
}

template <class T>
void fstream::read(const std::string &name, T *data, const Dims &selectionStart,
                   const Dims &selectionCount, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, reinterpret_cast<IOType *>(data),
                   Box<Dims>(selectionStart, selectionCount), blockID);
}

template <class T>
void fstream::read(const std::string &name, T *data, const Dims &start, const Dims &count,
                   const size_t stepStart, const size_t stepCount, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, reinterpret_cast<IOType *>(data), Box<Dims>(start, count),
                   Box<size_t>(stepStart, stepCount), blockID);
}

template <class T>
void fstream::read(const std::string &name, T &value, const size_t stepStart, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    m_Stream->Read(name, &reinterpret_cast<IOType &>(value), Box<size_t>(stepStart, 1), blockID);
}

template <class T>
std::vector<T> fstream::read(const std::string &name, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    auto vec = m_Stream->Read<IOType>(name, blockID);
    return reinterpret_cast<std::vector<T> &>(vec);
}

template <class T>
std::vector<T> fstream::read(const std::string &name, const size_t stepsStart,
                             const size_t stepsCount, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    auto vec = m_Stream->Read<IOType>(name, Box<size_t>(stepsStart, stepsCount), blockID);
    return reinterpret_cast<std::vector<T> &>(vec);
}

template <class T>
std::vector<T> fstream::read(const std::string &name, const Dims &start, const Dims &count,
                             const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    auto vec = m_Stream->Read<IOType>(name, Box<Dims>(start, count), blockID);
    return reinterpret_cast<std::vector<T> &>(vec);
}

template <class T>
std::vector<T> fstream::read(const std::string &name, const Dims &start, const Dims &count,
                             const size_t stepStart, const size_t stepCount, const size_t blockID)
{
    using IOType = typename TypeInfo<T>::IOType;
    auto vec = m_Stream->Read<IOType>(name, Box<Dims>(start, count),
                                      Box<size_t>(stepStart, stepCount), blockID);
    return reinterpret_cast<std::vector<T> &>(vec);
}

template <class T>
std::vector<T> fstream::read_attribute(const std::string &name, const std::string &variableName,
                                       const std::string separator)
{
    using IOType = typename TypeInfo<T>::IOType;
    std::vector<T> data;
    core::Attribute<IOType> *attribute = m_Stream->m_IO->InquireAttribute<IOType>(name);

    if (attribute == nullptr)
    {
        return data;
    }

    data.resize(attribute->m_Elements);
    m_Stream->ReadAttribute(name, reinterpret_cast<IOType *>(data.data()), variableName, separator);
    return data;
}

} // end namespace adios2

#endif /* ADIOS2_BINDINGS_CXX11_CXX11_FSTREAM_ADIOS2FSTREAM_TCC_ */
