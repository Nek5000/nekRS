/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Span.h
 *
 *  Created on: Apr 17, 2022
 *      Author: Jason Wang jason.ruonan.wang@gmail.com
 */

#ifndef ADIOS2_CORE_SPAN_H_
#define ADIOS2_CORE_SPAN_H_

#include "adios2/core/VariableBase.h"

namespace adios2
{
namespace core
{

template <class T>
class Span
{
public:
    std::pair<size_t, size_t> m_MinMaxMetadataPositions;

    // internal position variables from which the engine
    // can return a valid pointer any time
    // BP5 needs two levels of reference, BP3/4 uses only one
    size_t m_PayloadPosition = 0;
    int m_BufferIdx = -1;

    T m_Value = T{};

    Span(Engine &engine, const size_t size);
    ~Span() = default;

    size_t Size() const noexcept;
    T *Data() const noexcept;

    T &At(const size_t position);

    T &operator[](const size_t position);

private:
    Engine &m_Engine;
    size_t m_Size = 0;
};

} // end namespace core
} // end namespace adios2

#endif // ADIOS2_CORE_SPAN_H_
